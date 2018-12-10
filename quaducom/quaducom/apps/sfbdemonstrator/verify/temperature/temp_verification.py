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
    Int, List

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel

from ibvpy.rtrace.rt_domain_list_field import \
    RTraceDomainListField

#from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
#    MATS2DElastic
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

from ibvpy.fets.fets2D import \
    FETS2D4Q, FETS2D4Q8U, FETS2D4Q16U

from ibvpy.mats.mats2D import MATS2DElastic

from ibvpy.mesh.fe_grid import \
    FEGrid

from mathkit.mfn import MFnLineArray
from numpy import array, tensordot, dot, zeros, c_, ix_, max
from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos
#from rsurface_reader import \
#    read_rsurface, normalize_rsurfaces

# Interpolation
from scipy.interpolate import Rbf

from simiter.sim_pstudy import ISimModel, SimOut, SimPStudy

from ibvpy.examples.ts09_initial_strain_for_123D import TemperatureLinFn

from numpy import shape



class SFBMushRoofModel(IBVModel):
    '''SFB - Demontrator model specification.
    '''
    implements(ISimModel)

    # number of elements in all dims

    n_elems_xy = Int(10, ps_levels = (20, 80, 3))


    n_dofs_xy = Property(Int, depends_on = '+ps_levels')
    def _get_n_dofs_xy(self):
        if self.fets == self.fe_2D_linear:
            return self.n_elems_xy + 1
        elif self.fets == self.fe_2D_quadratic:
            return int(self.n_elems_xy * 2)
        else:
            raise ValueError

    rtrace_list = Property(List, depends_on = '+ps_levels')
    @cached_property
    def _get_rtrace_list(self):
        return [  self.max_princ_stress, self.sig_app, self.u ]

#    sig_trace = RTraceDomainListField( name = 'Stress' ,
#                               var = 'sig_app', warp = False,
#                               record_on = 'update' )
#    eps_trace = RTraceDomainListField( name = 'Epsilon' ,
#                                       var = 'eps_app', warp = True,
#                                       record_on = 'update' )
#    eps0_trace = RTraceDomainListField( name = 'Epsilon 0' ,
#                                       var = 'eps0_app', warp = True,
#                                       record_on = 'update' )
#    eps1t_trace = RTraceDomainListField( name = 'Epsilon 1-t' ,
#                                       var = 'eps1t_app', warp = True,
#                                       record_on = 'update' )
#    u_trace = RTraceDomainListField( name = 'Displacement' ,
#                                       var = 'u', warp = True,
#                                       record_on = 'update' )

    # dimensions of the shell structure
    length_xy = Float(1.) # [m]

    E = Float(30000) # [MPa]
    nu = Float(0.2) # [-]
    alpha = Float(1e-3)
    # variable type of the finite element
    fets = Instance(FETSEval,
                     ps_levels = [ 'fe_2D_linear', 'fe_2D_quadratic'])

    def _fets_default(self):
        return self.fe_2D_quadratic

    mats = Instance(MATS2DElastic)
    def _mats_default(self):
        return MATS2DElastic(E = self.E, nu = self.nu,
                              initial_strain = TemperatureLinFn(length = self.length_xy
                                                                , n_dims = 2
                                                                , T_right = 50
                                                                , T_left = 50
                                                                , offset = 0.5
                                                                , alpha = self.alpha))

    fe_2D_linear = Instance(FETSEval, transient = True)
    def _fe_2D_linear_default(self):
        return FETS2D4Q(mats_eval = self.mats)

    fe_2D_quadratic = Instance(FETSEval, transient = True)
    def _fe_2D_quadratic_default(self):
        return FETS2D4Q8U(mats_eval = self.mats)


    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name = 'u_z_free_corner', unit = 'm'),
                 SimOut(name = 'max principle stress', unit = 'MPa'),
                 SimOut(name = 'max sig_yy', unit = 'MPa'),
                 SimOut(name = 'max sig_xx', unit = 'MPa'), ]

    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()
        u_corner = U[ self.center_top_dof][-1, -1, 0]
        max_princ_stress = max(self.max_princ_stress._get_field_data().flatten())

        max_sig_yy = max(self.sig_app._get_field_data()[:, 4])
        max_sig_xx = max(self.sig_app._get_field_data()[:, 0])

        return  array([ u_corner, max_princ_stress, max_sig_yy, max_sig_xx ], dtype = 'float_')


    tline = Instance(TLine)
    def _tline_default(self):
        return TLine(min = 0.0, step = 1.0, max = 1.0)

    max_princ_stress = Property(Instance(RTraceDomainListField), depends_on = '+ps_levels')
    @cached_property
    def _get_max_princ_stress(self):
        return RTraceDomainListField(name = 'max principle stress' , idx = 0,
   #                                   position = 'int_pnts',
                                      var = 'max_principle_sig',
                                      record_on = 'update',)

    sig_app = Property(Instance(RTraceDomainListField), depends_on = '+ps_levels')
    @cached_property
    def _get_sig_app(self):
        return RTraceDomainListField(name = 'sig_app' ,
    #                                  position = 'int_pnts',
                                      var = 'sig_app',
                                      record_on = 'update',)

    u = Property(Instance(RTraceDomainListField), depends_on = '+ps_levels')
    @cached_property
    def _get_u(self):
        return RTraceDomainListField(name = 'displacement' ,
                                      var = 'u', warp = True,
                                      record_on = 'update',)

    #[ self.sig_trace, self.eps_trace, self.eps0_trace, self.eps1t_trace, self.u_trace]#, self.f_w_diagram ]

    fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels')
    def _get_fe_grid(self):
        return FEGrid(coord_min = (0.0, 0.0),
                       coord_max = (1.0, 1.0),
                       shape = (self.n_elems_xy, self.n_elems_xy),
                       fets_eval = self.fets)

    # time loop
    tloop = Property(depends_on = '+ps_levels')
    @cached_property
    def _get_tloop(self):

        self.fets.vtk_r *= 0.95

        domain = self.fe_grid

        self.center_top_dof = domain[-1, -1, -1, -1].dofs

        # NOTE: additional line-loads at the edge of the roof need to be considered!  

#        upper_surface = domain[:, :, -1, :, :, -1]
#        whole_domain = domain[:, :, :, :, :, :]

        bcond_list = [BCSlice(var = 'u', dims = [0, 1], slice = domain[0, 0, 0, 0], value = 0),
                      BCSlice(var = 'u', dims = [0], slice = domain[0, -1, 0, -1], value = 0),
                      ]


#       w_z = domain[-1, -1, -1, -1].dofs[0]

#       self.f_w_diagram = RTraceGraph( name = 'load - corner deflection',
#                                           var_x = 'U_k', idx_x = w_z,
#                                           var_y = 'time', idx_y = 0,
#                                           record_on = 'update' )
#        rtrace_list = self.rtrace_list#[ self.f_w_diagram ] + self.rtrace_list


        ts = TS(sdomain = [domain],
                 dof_resultants = True,
                 bcond_list = bcond_list,
                 rtrace_list = self.rtrace_list
               )
        # Add the time-loop control
        tloop = TLoop(tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline)
        return tloop

if __name__ == '__main__':

    sim_model = SFBMushRoofModel(n_elems_xy = 10, n_elems_z = 1)
#    interior_elems = sim_model.fe_grid[ 1:-1, 1:-1, 1:-1, :, :, : ].elems
#    sim_model.fe_grid.inactive_elems = list( interior_elems )

    do = 'ps'

    if do == 'eval':
        sim_model.n_elems_xy = 14
        print('eval', sim_model.peval())

    if do == 'ui':

        sim_model.n_elems_xy = 10
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
