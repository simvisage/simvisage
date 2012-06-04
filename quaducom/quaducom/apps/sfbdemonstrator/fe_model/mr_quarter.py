
#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Dec 9, 2009 by: rch

'''
MUSHROOF - Demostrator SFB532

@todos 
1) geo_transform for the bottom bar [Done]

'''
from etsproxy.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool, HasTraits

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel


from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic


from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U
from ibvpy.mesh.fe_grid import \
    FEGrid
from mathkit.mfn import MFnLineArray

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, shape, \
    cos, sin, arctan, where, abs, all, any, diag

from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos, pi as Pi

from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

# Interpolation
from scipy.interpolate import Rbf

from geo_column import GEOColumn

from simiter.sim_pstudy import ISimModel, SimOut, SimPStudy

from hp_shell import HPShell

from mush_roof_model import MushRoofModel


class MRquarter(MushRoofModel):

    implements(ISimModel)
    mushroof_part = 'quarter'

    n_elems_xy_quarter = Int(4, ps_levels = [3, 15, 5])
    n_elems_z = Int(2, ps_levels = [1, 4, 2])

    #----------------------------------------------------
    # elements
    #----------------------------------------------------
    vtk_r = Float(1.0)
    #default roof
    fe_roof = Instance(FETSEval,
                        ps_levels = ['fe_linear',
                                     'fe2d5_quad_serendipity',
                                     'fe_quad_serendipity',
                                     'fe_quad_lagrange' ] ,
                                     depends_on = '+initial_strain_roof, +initial_strain_col, +vtk_r')
    def _fe_roof_default(self):
        fets = self.fe_quad_serendipity
        fets.vtk_r *= self.vtk_r
        return fets

    #----------------------------------------------------
    # grid and geometric transformation
    #----------------------------------------------------
    fe_grid_roof = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_grid_roof(self):
        return FEGrid(coord_min = (0.0, 0.0, 0.0),
                       coord_max = (1.0, 1.0, 1.0),
                       geo_transform = self.hp_shell,
                       shape = (self.n_elems_xy, self.n_elems_xy, self.n_elems_z),
                       fets_eval = self.fe_quad_serendipity)

    fe_quad_serendipity = Property(Instance(FETSEval, transient = True), depends_on = '+input')
    def _get_fe_quad_serendipity(self):
        return FETS3D8H20U(mats_eval = self.mats_roof)

    shrink_factor = Float(1.0)
    # shell
    #
    hp_shell = Property(Instance(HPShell) , depends_on = '+ps_levels, +input')
    @cached_property
    def _get_hp_shell(self):
        return HPShell(length_xy_quarter = self.length_xy_quarter / self.shrink_factor ,
                        length_z = self.length_z / self.shrink_factor,
                        n_elems_xy_quarter = self.n_elems_xy_quarter,
                        n_elems_z = self.n_elems_z,
                        scalefactor_delta_h = self.scalefactor_delta_h,
                        mushroof_part = 'quarter',
                        shift_elems = False,
                        X0 = self.X0)

    #----------------------------------------------------
    # ps_study
    #----------------------------------------------------
    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()
        u_center_top_z = U[ self.center_top_dof ][0, 0, 2]
        max_princ_stress = max(self.max_princ_stress._get_field_data().flatten())
        return array([ u_center_top_z, max_princ_stress ],
                        dtype = 'float_')

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name = 'u_z_free_corner', unit = 'm'),
                 SimOut(name = 'maximum principle stress', unit = 'MPa'), ]

    #----------------------------------------------------
    # response tracer
    #----------------------------------------------------

    rtrace_list = List
    def _rtrace_list_default(self):
        return [  self.max_princ_stress, self.sig_app, self.u]

    #----------------------------------------------------
    # time loop
    #----------------------------------------------------

    tloop = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_tloop(self):
        domain = self.fe_grid_roof

        self.center_top_dof = domain[-1, -1, -1, -1, -1, -1].dofs

        #----------------------------------------------------
        # loading and boundaries
        #----------------------------------------------------


        #--- LC1: dead load
        # g = 22.4 kN/m^3 
        # orientation: global z-direction; 
        material_density_roof = -22.4e-3    # [MN/m^3]

        #--- LC2 additional dead load 
        # gA = 0,20 kN/m^2 
        # orientation: global z-direction (following the curved structure); 
        additional_dead_load = -0.20e-3 # [MN/m^2]

        #--- LC3 snow
        # s = 0,79 kN/m^2 
        # orientation: global z-direction (projection); 
        surface_load_s = -0.85e-3 # [MN/m^2]

        #--- LC4 wind (pressure) 
        # w = 0,13 kN/m^2 
        # orientation: local t-direction (surface normal); 
        surface_load_w = -0.13e-3 # [MN/m^2]

        # NOTE: additional line-loads at the edge of the roof need to be considered!  

        upper_surface = domain[:, :, -1, :, :, -1]
        whole_domain = domain[:, :, :, :, :, :]
        force_bc = [
                    # own weight
                     BCSlice(var = 'f', value = material_density_roof, dims = [2],
                          integ_domain = 'global',
                              slice = whole_domain),

                     # LC2: additional dead-load
#                     BCSlice( var = 'f', value = additional_dead_load, dims = [2],
#                              integ_domain = 'global',
#                              slice = upper_surface ),
#
#                     # LC3: snow load         
#                     BCSlice( var = 'f', value = surface_load_s, dims = [2],
#                              integ_domain = 'global',
#                              slice = upper_surface ),
#
#                     # LC3: wind
#                     BCSlice( var = 'f', value = surface_load_w, dims = [2],
#                              integ_domain = 'global',
#                              slice = upper_surface )
                   ]

        bc_symplane_yz = BCSlice(var = 'u', value = 0.  , dims = [0], slice = domain[0, :, :, 0, :, :])
        bc_symplane_xz = BCSlice(var = 'u', value = 0.  , dims = [1], slice = domain[:, 0, :, :, 0, :])


        bc_support_000 = BCSlice(var = 'u', value = 0.  , dims = [2], slice = domain[0, 0, 0, :, : , 0])

#        bc_column = [
#                     BCSlice( var = 'u'  , dims = [0, 1, 2],
#                              slice = domain[self.n_elems_xy_quarter - 1,
#                                               self.n_elems_xy_quarter - 1,
#                                               0,
#                                               0, -1, 0 ],
#                              value = 0. ),
#                    BCSlice( var = 'u'  , dims = [0, 1, 2],
#                            slice = domain[self.n_elems_xy_quarter - 1,
#                                           self.n_elems_xy_quarter - 1 ,
#                                           0,
#                                           - 1, 0, 0],
#                            value = 0. )]

        #bc_corner_load   = BCSlice( var = 'f', value = -nodal_load, dims = [2], slice = domain[-1,-1,-1,-1,-1,-1] )
        #bc_topface_load  = BCSlice( var = 'f', value = -nodal_load, dims = [2], slice = domain[:,:,-1,:,:,-1] )

        w_z = domain[-1, -1, -1, -1, -1, -1].dofs[0, 0, 2]

        self.f_w_diagram = RTraceGraph(name = 'load - corner deflection',
                                           var_x = 'U_k', idx_x = w_z,
                                           var_y = 'time', idx_y = 0,
                                           record_on = 'update')

        rtrace_list = [ self.f_w_diagram ] + self.rtrace_list

        ts = TS(sdomain = [domain],
                 dof_resultants = True,
                 bcond_list = [ bc_symplane_yz,
                                 bc_symplane_xz,
                               bc_support_000] + force_bc,
                 rtrace_list = rtrace_list
               )

        # Add the time-loop control
        tloop = TLoop(tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline)

        return tloop

if __name__ == '__main__':

    sim_model = MRquarter(n_elems_xy_quarter = 10,
                          n_elems_z = 1,
                          shrink_factor = 1.0
                          )
    #sim_model.initial_strain_roof = True
#    interior_elems = sim_model.fe_grid_column[ 1:-1, 1:-1, :, :, :, : ].elems
#    sim_model.fe_grid_column.inactive_elems = list( interior_elems )

    do = 'ui'

    if do == 'eval':
        print 'eval', sim_model.peval()

    if do == 'ui':

        sim_model.peval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource = sim_model)
        app.main()

    elif do == 'ps':

        sim_ps = SimPStudy(sim_model = sim_model)
        sim_ps.configure_traits()

    elif do == 'pickle':

        import pickle
        filename = '/tmp/sim.pickle'
        file = open(filename, 'w')
        pickle.dump(sim_model, file)
        file.close()
        file = open(filename, 'r')
        sm = pickle.load(file)
        file.close()

