

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
from enthought.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel

from ibvpy.rtrace.rt_domain_list_field import \
    RTraceDomainListField

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic
from ibvpy.mats.mats3D.mats3D_sdamage.mats3D_sdamage import \
    MATS3DScalarDamage
from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import \
    MATS3DMicroplaneDamage
from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
    MATS2D5MicroplaneDamage, PhiFnGeneral, PhiFnStrainHardening

from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U
from ibvpy.fets.fets2D5.fets2D58h import \
    FETS2D58H
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U

from ibvpy.mesh.fe_grid import \
    FEGrid

from mathkit.mfn import MFnLineArray
from numpy import array, tensordot, dot, zeros, c_, ix_, max, diag
from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos
from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

# Interpolation
from scipy.interpolate import Rbf

from simiter.sim_pstudy import ISimModel, SimOut, SimPStudy

from hp_shell import HPShell

def temperature_strain(X_pnt, x_pnt):

    alpha = 1.3e-5;
    t_up = +10.5
    t_lo = +9.5
    delta_t = t_lo + (t_up - t_lo) * x_pnt[2]
    epsilon_0 = alpha * delta_t

    return diag([ epsilon_0 for i in range(3) ])

class SFBMushRoofModel(IBVModel):
    '''SFB - Demonstrator model specification.
    '''
    implements(ISimModel)

    # dimensions of one quarter of the shell structure [m]
    #
    length_xy = Float(4.)
    length_z = Float(1.062)
    t_shell = Float(0.06)

    # choose model and discretization:
    #
    mushroof_part = 'quarter'
    n_elems_xy_quarter = Int(5, ps_levels = (17, 29, 1))
    n_elems_z = Int(3, input = True)
    n_elems_col = Int(1, input = True)# number of dofs used for column refinement   
    n_elems_edge = Int(1, input = True)
    t_edge = 0.03
    # grid parameters:
    #
    shift_elems_column = Bool(False, input = True)
    const_edge_elem = Bool (False, input = True)
    width_column = Float(0.45, input = True) # [m]
    delta_h_scalefactor = Float(1.30, ps_levels = (1.0, 1.6 , 1))

    # Material properties: Youngs-modulus, poision ratio 
    #
    E = Float(28700) # [MPa]
    nu = Float(0.2) # [-]

    def _elem_to_dofs(self, elems):
        if self.fets == self.fe_linear:
            return int(2 * elems + 1)
        elif self.fets == self.fe_quad_serendipity \
            or self.fets == self.fe2d5_quad_serendipity \
            or self.fets == self.fe_quad_lagrange:
            return int (elems * 2 + 1)
        else:
            raise ValueError

    # variable type of the finite element
    fets = Instance(FETSEval,
                     ps_levels = [ 'fe_linear',
                                  'fe2d5_quad_serendipity',
                                  'fe_quad_serendipity',
                                  'fe_quad_lagrange' ])
    def _fets_default(self):
        return self.fe_quad_serendipity
        #return self.fe_linear

    mats = Instance(MATS3DElastic)
    def _mats_default(self):
        return MATS3DElastic(E = self.E, nu = self.nu, initial_strain = temperature_strain)

    fe_linear = Instance(FETSEval, transient = True)
    def _fe_linear_default(self):
        return FETS3D8H(mats_eval = self.mats)

    fe_quad_serendipity = Instance(FETSEval, transient = True)
    def _fe_quad_serendipity_default(self):
        return FETS3D8H20U(mats_eval = self.mats)

    fe2d5_quad_serendipity = Instance(FETSEval, transient = True)
    def _fe2d5_quad_serendipity_default(self):
        return FETS2D58H20U(mats_eval = self.mats)

    fe_quad_lagrange = Instance(FETSEval, transient = True)
    def _fe_quad_lagrange_default(self):
        return FETS3D8H27U(mats_eval = self.mats)



    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name = 'u_z_free_corner', unit = 'm'),
                 SimOut(name = 'maximum principle stress', unit = 'MPa'), ]

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

    tline = Instance(TLine)
    def _tline_default(self):
        return TLine(min = 0.0, step = 1.0, max = 1.0)

    max_princ_stress = Instance(RTraceDomainListField)
    def _max_princ_stress_default(self):
        return RTraceDomainListField(name = 'max principle stress' , idx = 0,
                                      var = 'max_principle_sig', warp = True,
                                      record_on = 'update',)

    u = Property(Instance(RTraceDomainListField), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_u(self):
        return RTraceDomainListField(name = 'displacement' ,
                                      var = 'u', warp = True,
                                      record_on = 'update',)
    rtrace_list = List
    def _rtrace_list_default(self):
        return [  self.max_princ_stress, self.u ]
#                                 RTraceDomainListField( name = 'Displacement' ,
#                                                var = 'u', idx = 0, warp = True ),
#                                 RTraceDomainListField( name = 'Stress' ,
#                                                var = 'sig_app', idx = 0, warp = True,
#                                                record_on = 'update', ),

    hp_shell = Property(Instance(HPShell), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_hp_shell(self):
        return HPShell(mushroof_part = self.mushroof_part,
                        shift_elems_column = self.shift_elems_column,
                        n_elems_xy = self.n_elems_xy_quarter,
                        n_elems_col = self.n_elems_col,
                        n_elems_z = self.n_elems_z,
                        n_elems_edge = self.n_elems_edge,
                        const_edge_elem = self.const_edge_elem,
                        delta_h_scalefactor = self.delta_h_scalefactor)

    fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_grid(self):
        return FEGrid(coord_min = (0.0, 0.0, 0.0),
                       coord_max = (1.0, 1.0, 1.0),
                       geo_transform = self.hp_shell,
                       shape = (self.n_elems_xy_quarter, self.n_elems_xy_quarter, self.n_elems_z),
                       fets_eval = self.fets)

    # time loop
    tloop = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_tloop(self):

        self.fets.vtk_r *= 0.95

        domain = self.fe_grid

        self.center_top_dof = domain[-1, -1, -1, -1, -1, -1].dofs

        #----------------------------------------------------
        # loading cases (LC):
        #----------------------------------------------------

        #--- LC1: dead load
        # g = 22.4 kN/m^3 
        # orientation: global z-direction; 
        material_density = -0.0224 # [MN/m^3]

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
                     # LC1: dead load
                     BCSlice(var = 'f', value = material_density, dims = [2],
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
#                     # LC3: snow load         
#                     BCSlice( var = 'f', value = surface_load_w, dims = [2],
#                              integ_domain = 'global',
#                              slice = upper_surface )
                   ]

        bc_symplane_yz = BCSlice(var = 'u', value = 0.  , dims = [0], slice = domain[0, :, :, 0, :, :])
        bc_symplane_xz = BCSlice(var = 'u', value = 0.  , dims = [1], slice = domain[:, 0, :, :, 0, :])
        bc_support_000 = BCSlice(var = 'u', value = 0.  , dims = [2], slice = domain[0, 0, 0, :, :, 0])
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
                                 bc_support_000,
                               ] + force_bc,
                 rtrace_list = rtrace_list
               )
        # Add the time-loop control
        tloop = TLoop(tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline)

        return tloop

if __name__ == '__main__':

    sim_model = SFBMushRoofModel(n_elems_xy_quarter = 8,
                                  n_elems_z = 3,
                                  const_edge_elem = False,
                                  shift_elems_column = True,
                                  t_edge = 0.03,
                                  width_column = 0.45,
                                  n_elems_col = 1,
                                  n_elems_edge = 1,
                                  delta_h_scalefactor = 1.0)

#    interior_elems = sim_model.fe_grid[ 1:-1, 1:-1, 1:-1, :, :, : ].elems
#    sim_model.fe_grid.inactive_elems = list( interior_elems )

    do = 'ui'

    if do == 'eval':
        print 'eval', sim_model.peval()

    if do == 'ui':

        sim_model.peval()
        print sim_model.n_elems_z

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
