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
    Int, List, Bool, HasTraits, Enum, Str

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
    MATS2D5MicroplaneDamage

from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
    IPhiFn, PhiFnGeneralExtended, \
    PhiFnGeneral, PhiFnStrainHardening, PhiFnStrainHardeningLinear, \
    PhiFnGeneralExtendedExp

from ibvpy.api import \
    TStepper as TS, TLoop, TLine, \
    IBVModel, DOTSEval, \
    RTraceGraph, RTraceDomainListField, \
    BCDof, BCDofGroup, BCSlice

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
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U
from ibvpy.mesh.fe_grid import \
    FEGrid
from mathkit.mfn import MFnLineArray

import numpy as np

# Interpolation
from simiter.sim_pstudy import ISimModel, SimOut, SimPStudy

from hp_shell import HPShell

from mush_roof_model import MushRoofModel

from matresdev.db.exdb.ex_run_view import \
    ExRunView

from matresdev.db.matdb.trc.ccs_unit_cell import \
    CCSUnitCell, DamageFunctionEntry

from matresdev.db.simdb import \
    SimDB

from matresdev.db.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

simdb = SimDB()

from pickle import dump, load



class MRquarter(MushRoofModel):

    implements(ISimModel)
    mushroof_part = 'quarter'

    n_elems_xy_quarter = Int(6, ps_levels=[3, 15, 5])
    n_elems_z = Int(2, ps_levels=[1, 4, 2])

    #----------------------------------------------------
    # elements
    #----------------------------------------------------
    vtk_r = Float(0.9)

    # default roof
    fe_roof = Instance(FETSEval,
                        ps_levels=['fe2d5_quad_serendipity',
                                     'fe_quad_serendipity',
                                     'fe_linear',
                                     'fe_quad_lagrange' ] ,
                                     depends_on='+initial_strain_roof, +initial_strain_col, +vtk_r')
    def _fe_roof_default(self):
        fets = self.fe2d5_quad_serendipity
#        fets = self.fe_quad_serendipity
        fets.vtk_r *= self.vtk_r
        return fets

    #----------------------------------------------------
    # grid and geometric transformation
    #----------------------------------------------------
    fe_grid_roof = Property(Instance(FEGrid), depends_on='+ps_levels, +input')
    @cached_property
    def _get_fe_grid_roof(self):
        return FEGrid(coord_min=(0.0, 0.0, 0.0),
                       coord_max=(1.0, 1.0, 1.0),
                       geo_transform=self.hp_shell,
                       shape=(self.n_elems_xy, self.n_elems_xy, self.n_elems_z),
                       fets_eval=self.fe2d5_quad_serendipity)
#                       fets_eval = self.fe_roof)
#                       fets_eval = self.fe_quad_serendipity)

    mats_roof = Property(Instance(MATS2D5MicroplaneDamage), depends_on='+input')
#    mats_roof = Property( Instance( MATS3DElastic), depends_on = '+input' )
    @cached_property
    def _get_mats_roof(self):
        # return MATS3DElastic(E=self.E_roof, nu=self.nu)
        return MATS2D5MicroplaneDamage(
                                E=29100.0,
                                nu=0.2,
                                n_mp=30,
                                symmetrization='sum-type',
                                model_version='compliance',
                                phi_fn=self.phi_fn)

    fe2d5_quad_serendipity = Property(Instance(FETSEval, transient=True), depends_on='+input')
    def _get_fe2d5_quad_serendipity(self):
        return FETS2D58H20U(mats_eval=self.mats_roof)

    fe_quad_serendipity = Property(Instance(FETSEval, transient=True), depends_on='+input')
    def _get_fe_quad_serendipity(self):
        return FETS3D8H20U(mats_eval=self.mats_roof)

    shrink_factor = Float(1.0)
    # shell
    #
    hp_shell = Property(Instance(HPShell) , depends_on='+ps_levels, +input')
    @cached_property
    def _get_hp_shell(self):
        return HPShell(length_xy_quarter=self.length_xy_quarter / self.shrink_factor ,
                        length_z=self.length_z / self.shrink_factor,
                        n_elems_xy_quarter=self.n_elems_xy_quarter,
                        n_elems_z=self.n_elems_z,
                        scalefactor_delta_h=self.scalefactor_delta_h,
                        mushroof_part='quarter',
                        shift_elems=False,
                        X0=self.X0)

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
        print 'u_center_top_z', u_center_top_z
        return np.array([ u_center_top_z], dtype='float_')
#        max_princ_stress = max(self.max_princ_stress._get_field_data().flatten())
#        return np.array([ u_center_top_z, max_princ_stress ],
#                        dtype = 'float_')

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name='u_z_free_corner', unit='m'),
                 SimOut(name='maximum principal stress', unit='MPa'), ]

    #----------------------------------------------------
    # response tracer
    #----------------------------------------------------

    rtrace_list = List
    def _rtrace_list_default(self):
        return [  self.eps_app, self.sig_app, self.max_omega_i, self.phi_pdc ]

    max_lambda = Float(1.0, input=True)
    '''Maximum lambda factor to impose on the structure.
    The final loading loading level is calculated
    as the reference boundary conditions multiplied by the lambda_factor
    and uls_factor. Thus, lambda factor defines the load level as a multiple
    of the load level predicted by the linear analysis.
    '''

    f_w_diagram = Property(Instance(RTraceGraph), depends_on='+ps_levels, +input')
    @cached_property
    def _get_f_w_diagram(self):
        domain = self.fe_grid_roof
        w_z = domain[-1, -1, -1, -1, -1, -1].dofs[0, 0, 2]
        return RTraceGraph(name='load - corner deflection',
                           var_x='U_k', idx_x=w_z,
                           transform_x='-x',
                           var_y='time', idx_y=0,
                           # transform_y='y * %g' % self.lambda_factor,
                           record_on='update')

    n_steps = Int(15.0, auto_set=False, enter_set=False, input=True)
    time_fn_load = Instance(MFnLineArray, input=True)
    def _time_fn_load_default(self):
        return MFnLineArray(xdata=[0.0, 1.0, 3.0, 6.0, 15.0],
                            ydata=[0.0, 1.0, 1.0 / 0.22, 1.78 / 0.22, 9.45])

    boundary_x1 = Property(depends_on='+input')
    @cached_property
    def _get_boundary_x1(self):
        return self.fe_grid_roof.domain[-1, :, -1, -1, :, -1]

    #----------------------------------------------------
    # time loop
    #----------------------------------------------------

    tloop = Property(depends_on='+ps_levels, +input')
    @cached_property
    def _get_tloop(self):
        domain = self.fe_grid_roof

        #----------------------------------------------------
        # loading and boundaries
        #----------------------------------------------------

        #--- LC1: dead load
        # g = 22.4 kN/m^3
        # orientation: global z-direction;
        material_density_roof = -22.43e-3  # [MN/m^3]

        #--- LC2 additional dead load
        # gA = 0,20 kN/m^2
        # orientation: global z-direction (following the curved structure);
        additional_dead_load = -0.20e-3  # [MN/m^2]

        #--- LC2 additional boundary load
        # gA = 0,35 kN/m^2
        # orientation: global z-direction (following the curved structure);
        boundary_dead_load = -0.35e-3  # [MN/m]

        #--- LC3 snow
        # s = 0,79 kN/m^2
        # orientation: global z-direction (projection);
        surface_load_s = -0.85e-3  # [MN/m^2]

        #--- LC4 wind (pressure)
        # w = 0,13 kN/m^2
        # orientation: local t-direction (surface normal);
        surface_load_w = -0.13e-3  # [MN/m^2]

        # NOTE: additional line-loads at the edge of the roof need to be considered!

        upper_surface = domain[:, :, -1, :, :, -1]
        whole_domain = domain[:, :, :, :, :, :]
        boundary_x1 = domain[-1, :, -1, -1, :, -1]
        boundary_y1 = domain[:, -1, -1, :, -1, -1]

        time_fn_load = self.time_fn_load

        time_fn_permanent_load = MFnLineArray(xdata=[0.0, 1.0], ydata=[0.0, 1.0])
        time_fn_snow_load = MFnLineArray(xdata=[0.0, 1.0], ydata=[0.0, 0.0])

        force_bc = [
                    # own weight
                     BCSlice(name='self weight', var='f', value=material_density_roof, dims=[2],
                             integ_domain='global',
                             time_function=time_fn_load.get_value,
                             slice=whole_domain),

                     # LC2: additional dead-load
                     BCSlice(name='additional load', var='f', value=additional_dead_load, dims=[2],
                              integ_domain='global',
                              time_function=time_fn_load.get_value,
                              slice=upper_surface),

                     # LC2: additional boundary-load
                     BCSlice(name='additional boundary load 1', var='f', value=boundary_dead_load, dims=[2],
                              integ_domain='global',
                              time_function=time_fn_load.get_value,
                              slice=boundary_x1),

                     # LC2: additional boundary-load
                     BCSlice(name='additional boundary load 2', var='f', value=boundary_dead_load, dims=[2],
                              integ_domain='global',
                              time_function=time_fn_load.get_value,
                              slice=boundary_y1),
                     # LC3: snow load
                     BCSlice(name='snow load', var='f', value=surface_load_s, dims=[2],
                              integ_domain='global',
                              time_function=time_fn_snow_load.get_value,
                              slice=upper_surface),

#                     # LC3: wind
#                     BCSlice( var = 'f', value = surface_load_w, dims = [2],
#                              integ_domain = 'global',
#                              slice = upper_surface )
                   ]

        bc_symplane_yz = BCSlice(var='u', value=0.  , dims=[0], slice=domain[0, :, :, 0, :, :])
        bc_symplane_xz = BCSlice(var='u', value=0.  , dims=[1], slice=domain[:, 0, :, :, 0, :])
        bc_support_000 = BCSlice(var='u', value=0.  , dims=[2], slice=domain[0, 0, 0, :, : , 0])

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

        # bc_corner_load   = BCSlice( var = 'f', value = -nodal_load, dims = [2], slice = domain[-1,-1,-1,-1,-1,-1] )
        # bc_topface_load  = BCSlice( var = 'f', value = -nodal_load, dims = [2], slice = domain[:,:,-1,:,:,-1] )

#        support_z_dofs = domain[0, 0, 0, :, : , 0].dofs[:, :, 2]
#        support_f_w = RTraceGraph(name='force - corner deflection',
#                           var_x='time', idx_x=0,
#                           transform_x='x * %g' % lambda_failure,
#                           var_y='F_int', idx_y_arr=np.unique(support_z_dofs.flatten()),
#                           transform_y='-y',
#                           record_on='update')

        rtrace_list = [ self.f_w_diagram ] + self.rtrace_list

        ts = TS(sdomain=[domain],
                 dof_resultants=True,
                 bcond_list=[ bc_symplane_yz,
                                bc_symplane_xz,
                                bc_support_000] + force_bc,
                 rtrace_list=rtrace_list
               )

        step = self.n_steps
        # Add the time-loop control
        tloop = TLoop(tstepper=ts, RESETMAX=0, KMAX=70,
                       tolerance=0.5e-3,
                       tline=TLine(min=0.0, step=step, max=self.max_lambda))
        return tloop

class MRquarterDB(MRquarter):
    '''get 'phi_fn' as calibrated by the fitter and stored in the DB
    '''

    # vary the failure strain in PhiFnGeneralExtended:
    factor_eps_fail = Float(1.4, input=True,
                             ps_levels=(1.0, 1.2, 3))

    #-----------------
    # composite cross section unit cell:
    #-----------------
    #
    ccs_unit_cell_key = Enum('FIL-10-09_2D-05-11_0.00462_all0',
                              CCSUnitCell.db.keys(),
                              simdb=True, input=True,
                              auto_set=False, enter_set=True)

    ccs_unit_cell_ref = Property(Instance(SimDBClass),
                                  depends_on='ccs_unit_cell_key')
    @cached_property
    def _get_ccs_unit_cell_ref(self):
        return CCSUnitCell.db[ self.ccs_unit_cell_key ]

    #-----------------
    # damage function:
    #-----------------
    #
    material_model = Str(input=True)
    def _material_model_default(self):
        # return the material model key of the first DamageFunctionEntry
        # This is necessary to avoid an ValueError at setup
        return self.ccs_unit_cell_ref.damage_function_list[0].material_model

    calibration_test = Str(input=True)
    def _calibration_test_default(self):
        # return the material model key of the first DamageFunctionEntry
        # This is necessary to avoid an ValueError at setup
        return self.ccs_unit_cell_ref.damage_function_list[0].calibration_test

    damage_function = Property(Instance(MFnLineArray),
                                depends_on='+input')
    @cached_property
    def _get_damage_function(self):
        print 'getting damage function'
        return self.ccs_unit_cell_ref.get_param(self.material_model, self.calibration_test)

    #-----------------
    # phi function extended:
    #-----------------
    #
    phi_fn = Property(Instance(PhiFnGeneralExtended),
                       depends_on='+input,+ps_levels')
    @cached_property
    def _get_phi_fn(self):
        return PhiFnGeneralExtendedExp(mfn=self.damage_function, Dfp=0.01, Efp_frac=0.007)
#        return PhiFnGeneralExtended( mfn = self.damage_function,
#                                     factor_eps_fail = self.factor_eps_fail )

    #----------------------------------------------------------------------------------
    # mats_eval
    #----------------------------------------------------------------------------------

    # age of the plate at the time of testing
    # NOTE: that the same phi-function is used independent of age. This assumes a
    # an afine/proportional damage evolution for different ages.
    #
    age = Int(28,  # input = True
                )

    # composite E-modulus
    #
    E_c = Property(Float, depends_on='+input')
    @cached_property
    def _get_E_c(self):
        return self.ccs_unit_cell_ref.get_E_c_time(self.age)

    # Poisson's ratio
    #
    nu = Property(Float, depends_on='+input')
    @cached_property
    def _get_nu(self):
        return self.ccs_unit_cell_ref.nu



if __name__ == '__main__':

    do = 'load_factor_plot'
    # do = 'ui'
    do = 'mxn_cut'

    if do == 'mxn_cut':
        sim_model = MRquarterDB(ccs_unit_cell_key='FIL-10-09_2D-05-11_0.00462_all0',
                                calibration_test='TT-12c-6cm-0-TU-SH2-V1_age26_Ec29100_nu0.2_nsteps100_maxeps0.007_smoothed',
                                age=26,
                                max_lambda=1.0,
                                n_steps=1,
                                n_elems_xy_quarter=2,
                                n_elems_z=1,
                                )

        u = sim_model.tloop.eval()
        F_int = sim_model.tloop.tstepper.F_int

        dof_X = sim_model.fe_grid_roof[ 0, :, :, 0, :, : ].dof_X
        dof_X_1 = dof_X[:, :, 1].flatten()
        dof_X_2 = dof_X[:, :, 2].flatten()
        dofs = sim_model.fe_grid_roof[ 0, :, :, 0, :, : ].dofs
        dofs_0 = dofs[:, :, 0].flatten()
        dof_X_1_idx_sorted = np.argsort(dof_X_1)
        dof_X_1_sorted = dof_X_1[dof_X_1_idx_sorted]
        dof_X_2_sorted = dof_X_2[dof_X_1_idx_sorted]
        dofs_0_sorted = dofs_0[dof_X_1_idx_sorted]

        print 'dof_X_1'
        print dof_X_1

        print 'dof_X_2'
        print dof_X_2

        print 'dofs_0'
        print dofs_0

        print 'dof_X_1_sorted'
        print dof_X_1_sorted

        print 'dof_X_2_sorted'
        print dof_X_2_sorted

        print 'dofs_0_sorted'
        print dofs_0_sorted

        dofs_0_idx_unique = np.argsort(dofs_0_sorted)
        dofs_0_unique = dofs_0_sorted[dofs_0_idx_unique]
        print 'dofs_0_idx_unique'
        print dofs_0_idx_unique
        print 'dofs_0_unique'
        print dofs_0_unique

        print 'F_int_sorted'
        F_int_unique = F_int[np.unique(dofs_0_unique)]
        print F_int_unique


    else:
        sim_model = MRquarterDB(ccs_unit_cell_key='FIL-10-09_2D-05-11_0.00462_all0',
                                calibration_test='TT-12c-6cm-0-TU-SH2-V1_age26_Ec29100_nu0.2_nsteps100_maxeps0.007_smoothed',
                                age=26,
                                max_lambda=15.0,
                                n_steps=15,
                                n_elems_xy_quarter=6,
                                n_elems_z=2,
                                )

    # sim_model.initial_strain_roof = True
#    interior_elems = sim_model.fe_grid_column[ 1:-1, 1:-1, :, :, :, : ].elems
#    sim_model.fe_grid_column.inactive_elems = list( interior_elems )

    if do == 'eval':
#        sim_model.tloop.eval()

        print 'eval', sim_model.peval()

    if do == 'ui':
        sim_model.tloop.eval()
        # sim_model.peval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource=sim_model)
        app.main()

    elif do == 'ps':
        sim_ps = SimPStudy(sim_model=sim_model)
        sim_ps.configure_traits()

    elif do == 'load_factor_plot':
        import pylab as p

        sim_model.time_fn_load.plot(p)
        p.show()

        sim_model.tloop.eval()

        f_w = sim_model.f_w_diagram
        f_w.redraw()

        # f_w.trace.plot(p, color='black')
        w = f_w.trace.xdata
        time_ = f_w.trace.ydata
        lambda_ = sim_model.time_fn_load.get_values(time_)
        print 'lambda_', lambda_

        eta_nmd = 0.22
        chi = 0.81
        gamma = 1.5

        eta = lambda_ * eta_nmd
        eta_el = gamma / chi

        max_w = w[-1]
        max_lambda = np.max(lambda_)
        max_eta = np.max(eta)

        fig, ax1 = p.subplots()

        ax1.set_ylim(0, max_lambda * 1.1)
        ax1.plot(w, lambda_, color='black')

        ax2 = ax1.twinx()
        ax2.set_ylim(0, max_eta * 1.1)

        ax2.plot([0, max_w], [1.0, 1.0], color='black', linestyle='dashed')
        ax2.plot([0, max_w], [eta_el, eta_el], color='black', linestyle='dashed')
        ax2.plot([0, max_w], [ eta_nmd, eta_nmd], color='black', linestyle='dashed')
        p.show()

        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource=sim_model)
        app.main()

    elif do == 'pickle':

        import pickle
        filename = '/tmp/sim.pickle'
        fl = open(filename, 'w')
        pickle.dump(sim_model, fl)
        fl.close()
        fl = open(filename, 'r')
        sm = pickle.load(fl)
        fl.close()

