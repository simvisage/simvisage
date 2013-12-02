
from etsproxy.traits.api import \
    Array, Bool, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, \
    Callable, List, TraitDict, Any, Range, \
    Delegate, Event, on_trait_change, Button, \
    Interface, implements, Property, cached_property

from ibvpy.api import \
    TStepper as TS, TLoop, TLine, \
    IBVModel, DOTSEval, \
    RTraceGraph, RTraceDomainListField, \
    BCDof, BCDofGroup, BCSlice, FERefinementGrid, FEDomain

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
from ibvpy.fets.fets2D5.fets2D58h import \
    FETS2D58H

from ibvpy.mesh.fe_grid import \
    FEGrid

from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
    MATS2D5MicroplaneDamage

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
    IPhiFn, PhiFnGeneral, PhiFnGeneralExtended, PhiFnGeneralExtendedExp

from mathkit.geo.geo_ndgrid import \
    GeoNDGrid

from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import \
    MFnNDGrid, GridPoint

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

import numpy as np

from numpy import \
    sin, cos, c_, arange, hstack, array, max, frompyfunc, \
    linspace, argsort, zeros_like, append, size

from time import time
from os.path import join

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from matresdev.simiter.sim_pstudy import \
    SimPStudy, SimOut, ISimModel

from matresdev.db.matdb.trc.ccs_unit_cell import \
    CCSUnitCell, DamageFunctionEntry

from matresdev.db.simdb import \
    SimDB

from matresdev.db.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

from matresdev.db.exdb.ex_run import ExRun

import pylab as p

simdb = SimDB()

from pickle import dump, load

from quaducom.devproc.show_results import format_plot


class SimBT4PT(IBVModel):
    '''Simulation: four point bending test.
    '''

    input_change = Event
    @on_trait_change('+input,ccs_unit_cell.input_change')
    def _set_input_change(self):
        self.input_change = True

    implements(ISimModel)

    #-----------------
    # discretization:
    #-----------------

    # specify weather the elastomer is to be modeled or if the load is
    # introduced as line load
    #
    elstmr_flag = Bool(True)

    # discretization in x-direction (longitudinal):
    #
    outer_zone_shape_x = Int(6, input=True,
                             ps_levels=(4, 12, 3))

    # discretization in x-direction (longitudinal):
    #
    load_zone_shape_x = Int(2, input=True,
                            ps_levels=(1, 4, 1))

    # middle part discretization in x-direction (longitudinal):
    #
    mid_zone_shape_x = Int(3, input=True,
                           ps_levels=(1, 4, 1))

    # discretization in y-direction (width):
    #
    shape_y = Int(2, input=True,
                     ps_levels=(1, 4, 2))

    # discretization in z-direction:
    #
    shape_z = Int(2, input=True,
                     ps_levels=(1, 3, 3))

    #-----------------
    # geometry:
    #-----------------
    #
    # edge length of the bending specimen (beam) (entire length without symmetry)
    #
    length = Float(1.50, input=True)

    elstmr_length = Float(0.05, input=True)
    mid_zone_length = Float(0.50, input=True)
    elstmr_thickness = Float(0.005, input=True)
    width = Float(0.20, input=True)
    thickness = Float(0.06, input=True)

    #-----------------
    # derived geometric parameters
    #-----------------
    #
    # half the length of the elastomer (load introduction 
    # with included symmetry)
    #
    sym_specmn_length = Property
    def _get_sym_specmn_length(self):
        return self.length / 2.

    sym_mid_zone_specmn_length = Property
    def _get_sym_mid_zone_specmn_length(self):
        return self.mid_zone_length / 2.

    # half the length of the elastomer (load introduction
    # with included symmetry
    #
    sym_elstmr_length = Property
    def _get_sym_elstmr_length(self):
        return self.elstmr_length / 2.

    # half the specimen width
    #
    sym_width = Property
    def _get_sym_width(self):
        return self.width / 2.

    #----------------------------------------------------------------------------------
    # mats_eval
    #----------------------------------------------------------------------------------

    # age of the plate at the time of testing
    # NOTE: that the same phi-function is used independent of age. This assumes a 
    # an afine/proportional damage evolution for different ages. 
    #
    age = Int(28, input=True)

    # time stepping params 
    #  
    tstep = Float(0.05, auto_set=False, enter_set=True, input=True)
    tmax = Float(1.0, auto_set=False, enter_set=True, input=True)
    tolerance = Float(0.001, auto_set=False, enter_set=True, input=True)

    # specify type of 'linalg.norm'
    # default value 'None' sets norm to 2-norm,
    # i.e "norm = sqrt(sum(x_i**2))
    #
    # set 'ord=np.inf' to switch norm to
    # "norm = max(abs(x_i))"
    # 
    ord = Enum(np.inf, None)

    n_mp = Int(30, input=True)

    # @todo: for mats_eval the information of the unit cell should be used
    # in order to use the same number of microplanes and model version etc...
    #
    specmn_mats = Property(Instance(MATS2D5MicroplaneDamage),
                          depends_on='input_change')
    @cached_property
    def _get_specmn_mats(self):
        return MATS2D5MicroplaneDamage(
                                E=self.E_c,
#                                 E=self.E_m,
                                nu=self.nu,
                                # corresponding to settings in "MatsCalib"
                                n_mp=self.n_mp,
                                symmetrization='sum-type',
                                model_version='compliance',
                                phi_fn=self.phi_fn)

    if elstmr_flag:
        elstmr_mats = Property(Instance(MATS3DElastic),
                               depends_on='input_change')
        @cached_property
        def _get_elstmr_mats(self):
            # specify a small elastomer stiffness (approximation)
            E_elast = self.E_c / 10.
            print 'effective elastomer E_modulus', E_elast
            return MATS3DElastic(E=E_elast,
                                 nu=0.4)

    #-----------------
    # fets:
    #-----------------

    # specify element shrink factor in plot of fe-model
    #
    vtk_r = Float(0.95)

    # use quadratic serendipity elements
    #
    specmn_fets = Property(Instance(FETSEval),
                           depends_on='input_change')
    @cached_property
    def _get_specmn_fets(self):
        fets = FETS2D58H20U(mats_eval=self.specmn_mats)
        fets.vtk_r *= self.vtk_r
        return fets

    # use quadratic serendipity elements
    #
    elstmr_fets = Property(Instance(FETSEval),
                           depends_on='input_change')
    @cached_property
    def _get_elstmr_fets(self):
        fets = FETS2D58H20U(mats_eval=self.elstmr_mats)
        fets.vtk_r *= self.vtk_r
        return fets

    fe_domain = Property(depends_on='+ps_levels, +input')
    @cached_property
    def _get_fe_domain(self):
        return FEDomain()

    #===========================================================================
    # fe level 
    #===========================================================================

    outer_zone_specmn_fe_level = Property(depends_on='+ps_levels, +input')
    def _get_outer_zone_specmn_fe_level(self):
        return  FERefinementGrid(name='outer zone specimen patch',
                                 fets_eval=self.specmn_fets,
                                 domain=self.fe_domain)

    load_zone_specmn_fe_level = Property(depends_on='+ps_levels, +input')
    def _get_load_zone_specmn_fe_level(self):
        return  FERefinementGrid(name='load zone specimen patch',
                                 fets_eval=self.specmn_fets,
                                 domain=self.fe_domain)

    mid_zone_specmn_fe_level = Property(depends_on='+ps_levels, +input')
    def _get_mid_zone_specmn_fe_level(self):
        return  FERefinementGrid(name='mid zone specimen patch',
                                 fets_eval=self.specmn_fets,
                                 domain=self.fe_domain)

    elstmr_fe_level = Property(depends_on='+ps_levels, +input')
    def _get_elstmr_fe_level(self):
        return  FERefinementGrid(name='elastomer patch',
                                 fets_eval=self.elstmr_fets,
                                 domain=self.fe_domain)

    #===========================================================================
    # Grid definition 
    #===========================================================================

    mid_zone_specmn_fe_grid = Property(Instance(FEGrid), depends_on='+ps_levels, +input')
    @cached_property
    def _get_mid_zone_specmn_fe_grid(self):
        # only a quarter of the beam is simulated due to symmetry:
        fe_grid = FEGrid(coord_min=(0.,
                                      0.,
                                      0.),
                         coord_max=(self.sym_mid_zone_specmn_length - self.sym_elstmr_length,
                                      self.sym_width,
                                      self.thickness),
                         shape=(self.mid_zone_shape_x, self.shape_y, self.shape_z),
                         level=self.mid_zone_specmn_fe_level,
                         fets_eval=self.specmn_fets)
        return fe_grid

    load_zone_specmn_fe_grid = Property(Instance(FEGrid), depends_on='+ps_levels, +input')
    @cached_property
    def _get_load_zone_specmn_fe_grid(self):
        # only a quarter of the beam is simulated due to symmetry:
        fe_grid = FEGrid(coord_min=(self.sym_mid_zone_specmn_length - self.sym_elstmr_length,
                                      0.,
                                      0.),
                         coord_max=(self.sym_mid_zone_specmn_length + self.sym_elstmr_length,
                                      self.sym_width,
                                      self.thickness),
                         shape=(self.load_zone_shape_x, self.shape_y, self.shape_z),
                         level=self.load_zone_specmn_fe_level,
                         fets_eval=self.specmn_fets)
        return fe_grid

#    if elstmr_flag:
    elstmr_fe_grid = Property(Instance(FEGrid), depends_on='+ps_levels, +input')
    @cached_property
    def _get_elstmr_fe_grid(self):
        fe_grid = FEGrid(coord_min=(self.sym_mid_zone_specmn_length - self.sym_elstmr_length,
                                      0.,
                                      self.thickness),
                         coord_max=(self.sym_mid_zone_specmn_length + self.sym_elstmr_length,
                                      self.sym_width,
                                      self.thickness + self.elstmr_thickness),
                         level=self.elstmr_fe_level,
                         shape=(self.load_zone_shape_x, self.shape_y, 1),
                         fets_eval=self.elstmr_fets)
        return fe_grid

    outer_zone_specmn_fe_grid = Property(Instance(FEGrid), depends_on='+ps_levels, +input')
    @cached_property
    def _get_outer_zone_specmn_fe_grid(self):
        # only a quarter of the plate is simulated due to symmetry:
        fe_grid = FEGrid(coord_min=(self.sym_mid_zone_specmn_length + self.sym_elstmr_length,
                                      0.,
                                      0.),
                         coord_max=(self.sym_specmn_length,
                                      self.sym_width,
                                      self.thickness),
                         shape=(self.outer_zone_shape_x, self.shape_y, self.shape_z),
                         level=self.outer_zone_specmn_fe_level,
                         fets_eval=self.specmn_fets)
        return fe_grid

    #===========================================================================
    # Boundary conditions
    #===========================================================================
    w_max = Float(-0.030, input=True) # [m]

    w_max = Float(-0.030, input=True)  # [m]

    bc_list = Property(depends_on='+ps_levels, +input')
    @cached_property
    def _get_bc_list(self):
        mid_zone_specimen = self.mid_zone_specmn_fe_grid
        load_zone_specimen = self.load_zone_specmn_fe_grid
        outer_zone_specimen = self.outer_zone_specmn_fe_grid

        if self.elstmr_flag:
            elastomer = self.elstmr_fe_grid

        #--------------------------------------------------------------
        # boundary conditions for the symmetry 
        #--------------------------------------------------------------
        # symmetry in the xz-plane
        # (Note: the x-axis corresponds to the axis of symmetry along the longitudinal axis of the beam)
        #
        bc_outer_zone_symplane_xz = BCSlice(var='u', value=0., dims=[1],
                                            slice=outer_zone_specimen[:, 0, :, :, 0, :])
        bc_load_zone_symplane_xz = BCSlice(var='u', value=0., dims=[1],
                                           slice=load_zone_specimen[:, 0, :, :, 0, :])
        bc_mid_zone_symplane_xz = BCSlice(var='u', value=0., dims=[1],
                                          slice=mid_zone_specimen[:, 0, :, :, 0, :])

        if self.elstmr_flag:
            bc_el_symplane_xz = BCSlice(var='u', value=0., dims=[1],
                                        slice=elastomer[:, 0, :, :, 0, :])
        # symmetry in the yz-plane
        #
        bc_mid_zone_symplane_yz = BCSlice(mid='u', value=0., dims=[0],
                                          slice=mid_zone_specimen[0, :, :, 0, :, :])

        #--------------------------------------------------------------
        # boundary conditions for the support
        #--------------------------------------------------------------
        bc_support_0y0 = BCSlice(var='u', value=0., dims=[2],
                                 slice=outer_zone_specimen[-1, :, 0, -1, :, 0])

        #--------------------------------------------------------------
        # connect all grids 
        #--------------------------------------------------------------
        link_loadzn_outerzn = BCDofGroup(var='u', value=0., dims=[0, 1, 2],
                                get_dof_method=load_zone_specimen.get_right_dofs,
                                get_link_dof_method=outer_zone_specimen.get_left_dofs,
                                link_coeffs=[1.])
        link_midzn_loadzn = BCDofGroup(var='u', value=0., dims=[0, 1, 2],
                                get_dof_method=mid_zone_specimen.get_right_dofs,
                                get_link_dof_method=load_zone_specimen.get_left_dofs,
                                link_coeffs=[1.])

        if self.elstmr_flag:
            link_elstmr_loadzn_z = BCDofGroup(var='u', value=0., dims=[2],
                                    get_dof_method=elastomer.get_back_dofs,
                                    get_link_dof_method=load_zone_specimen.get_front_dofs,
                                    link_coeffs=[1.])

            # hold elastomer in a single point in order to avoid kinematic movement yielding singular K_mtx 
            #
            bc_elstmr_fix = BCSlice(var='u', value=0., dims=[0],
                                 slice=elastomer[0, 0, 0, 0, 0, 0])

        #--------------------------------------------------------------
        # loading
        #--------------------------------------------------------------
        # w_max = center displacement:
        w_max = self.w_max

        if self.elstmr_flag:
            # apply displacement at all top nodes of the elastomer (surface load)
            #
            bc_w = BCSlice(var='u', value=w_max, dims=[2],
                           slice=elastomer[:, :, -1, :, :, -1])
        else:
            bc_w = BCSlice(var='u', value=w_max, dims=[2],
                           # slice is only exactly in the center of the loading zone for 'load_zone_shape_x' = 2
                           # center line of the load zone
                           slice=load_zone_specimen[0, :, -1, -1, :, -1])

#        f_max = 0.010 / 4. / self.sym_width 
#        bc_line_f = BCSlice(var = 'f', value = f_max, dims = [2],
#                            # slice is only valid for 'load_zone_shape_x' = 2
#                            # center line of the load zone
#                            slice = load_zone_specimen[0, :, -1, -1, :, -1])

        bc_list = [bc_outer_zone_symplane_xz,
                   bc_load_zone_symplane_xz,
                   bc_mid_zone_symplane_xz,
                   bc_mid_zone_symplane_yz,
                   #
                   link_midzn_loadzn,
                   link_loadzn_outerzn,
                   bc_support_0y0,
                   #
                   bc_w,
                   ]

        if self.elstmr_flag:
            bc_list += [ bc_el_symplane_xz, link_elstmr_loadzn_z, bc_elstmr_fix ]

        return bc_list

    #----------------------
    # tloop
    #----------------------

    tloop = Property(depends_on='input_change')
    @cached_property
    def _get_tloop(self):

        #--------------------------------------------------------------
        # ts 
        #--------------------------------------------------------------

        mid_zone_spec = self.mid_zone_specmn_fe_grid
        load_zone_spec = self.load_zone_specmn_fe_grid
        outer_zone_spec = self.outer_zone_specmn_fe_grid

        if self.elstmr_flag:
            # ELSTRMR TOP SURFACE
            # dofs at elastomer top surface (used to integrate the force)
            #
            elastomer = self.elstmr_fe_grid
            elstmr_top_dofs_z = elastomer[:, :, -1, :, :, -1].dofs[:, :, 2].flatten()
            load_dofs_z = np.unique(elstmr_top_dofs_z)
            print 'load_dofs_z', load_dofs_z
        else:
            # LINE LOAD TOP OF LOAD ZONE
            # dofs at center line of the specmn load zone (used to integrate the force)
            # note slice index in x-direction is only valid for load_zone_shape_x = 2 !
            #
            load_zone_spec_topline_dofs_z = load_zone_spec[0, :, -1, -1, :, -1].dofs[:, :, 2].flatten()
            load_dofs_z = np.unique(load_zone_spec_topline_dofs_z)
            print 'load_dofs_z', load_dofs_z

        # SUPPRT LINE
        # dofs at support line of the specmn (used to integrate the force)
        #
        outer_zone_spec_supprtline_dofs_z = outer_zone_spec[-1, :, 0, -1, :, 0].dofs[:, :, 2].flatten()
        supprt_dofs_z = np.unique(outer_zone_spec_supprtline_dofs_z)
        print 'supprt_dofs_z', supprt_dofs_z

        # CENTER DOF (used for tracing of the displacement)
        #
        center_bottom_dof = mid_zone_spec[0, 0, 0, 0, 0, 0].dofs[0, 0, 2]
        print 'center_bottom_dof', center_bottom_dof

        # THIRDPOINT DOF (used for tracing of the displacement) 
        # dofs at center middle of the laod zone at the bottom side
        #
        # NOTE: slice index in x-direction is only valid for load_zone_shape_x = 2 !
        thirdpoint_bottom_dof = load_zone_spec[0, 0, 0, -1, 0, 0].dofs[0, 0, 2]
        print 'thirdpoint_bottom_dof', thirdpoint_bottom_dof

        # force-displacement-diagram (CENTER) 
        # 
        self.f_w_diagram_center = RTraceGraph(name='displacement_elasttop (center) - force',
                                       var_x='U_k'  , idx_x=center_bottom_dof,
                                       var_y='F_int', idx_y_arr=load_dofs_z,
                                       record_on='update',
                                       transform_x='-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the beam (2 symmetry axis):
                                       #
                                       transform_y='-4000. * y')

        # force-displacement-diagram_supprt (SUPPRT)
        # 
        self.f_w_diagram_supprt = RTraceGraph(name='displacement_supprtline (center) - force',
                                       var_x='U_k'  , idx_x=center_bottom_dof,
                                       var_y='F_int', idx_y_arr=supprt_dofs_z,
                                       record_on='update',
                                       transform_x='-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the beam (2 symmetry axis):
                                       #
                                       transform_y='4000. * y')

        # force-displacement-diagram (THIRDPOINT) 
        # 
        self.f_w_diagram_thirdpoint = RTraceGraph(name='displacement_elasttop (thirdpoint) - force',
                                       var_x='U_k'  , idx_x=thirdpoint_bottom_dof,
                                       var_y='F_int', idx_y_arr=load_dofs_z,
                                       record_on='update',
                                       transform_x='-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the beam (2 symmetry axis):
                                       #
                                       transform_y='-4000. * y')

        ts = TS(
                sdomain=self.fe_domain,
                bcond_list=self.bc_list,
                rtrace_list=[
                             self.f_w_diagram_center,
                             self.f_w_diagram_thirdpoint,
                             self.f_w_diagram_supprt,
                             RTraceDomainListField(name='Displacement' ,
                                            var='u', idx=0, warp=True),
#                             RTraceDomainListField(name = 'Stress' ,
#                                            var = 'sig_app', idx = 0, warp = True,
#                                            record_on = 'update'),
#                             RTraceDomainListField(name = 'Strain' ,
#                                        var = 'eps_app', idx = 0, warp = True,
#                                        record_on = 'update'),
#                             RTraceDomainListField(name = 'Damage' ,
#                                        var = 'omega_mtx', idx = 0, warp = True,
#                                        record_on = 'update'),
                             RTraceDomainListField(name='max_omega_i', warp=True,
                                        var='max_omega_i', idx=0,
                                        record_on='update'),
#                             RTraceDomainListField(name = 'IStress' ,
#                                            position = 'int_pnts',
#                                            var = 'sig_app', idx = 0,
#                                            record_on = 'update'),
#                             RTraceDomainListField(name = 'IStrain' ,
#                                            position = 'int_pnts',
#                                            var = 'eps_app', idx = 0,
#                                            record_on = 'update'),
                              ]
                )

        # Add the time-loop control
        tloop = TLoop(tstepper=ts,
                      KMAX=50,
                      tolerance=self.tolerance,
                      RESETMAX=0,
                      tline=TLine(min=0.0, step=self.tstep, max=self.tmax),
                      ord=self.ord
                      )

        return tloop

    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()

        self.f_w_diagram_center.refresh()
        F_max = max(self.f_w_diagram_center.trace.ydata)

        u_center_top_z = U[ self.center_top_dofs ][0, 0, 2]
        return array([ u_center_top_z, F_max ],
                        dtype='float_')

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name='u_center_top_z', unit='m'),
                 SimOut(name='F_max', unit='kN') ]


class SimBT4PTDB(SimBT4PT):

    # vary the failure strain in PhiFnGeneralExtended:
    factor_eps_fail = Float(1.0, input=True,
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
                               depends_on='input_change')
    @cached_property
    def _get_damage_function(self):
        return self.ccs_unit_cell_ref.get_param(self.material_model, self.calibration_test)

    #-----------------
    # phi function extended:
    #-----------------
    #
#    phi_fn = Property(Instance(PhiFnGeneralExtended),
#                       depends_on = 'input_change,+ps_levels')
#    @cached_property
#    def _get_phi_fn(self):
#        return PhiFnGeneralExtended(mfn = self.damage_function,
#                                     factor_eps_fail = self.factor_eps_fail)

    #-----------------
    # phi function General:
    #-----------------
    #
    # definition of 'phi_fn' in order to run scripting
    # @todo: how can the implementation be changed in order to allow for parameter passing, e.g. 'factor_eps_fail', etc.
    #
    phi_fn_class = Enum(PhiFnGeneral, [PhiFnGeneral, PhiFnGeneralExtended, PhiFnGeneralExtendedExp], input=True,)
    phi_fn = Property(Instance(phi_fn_class), depends_on='phi_fn_class, input_change, +ps_levels')
    @cached_property
    def _get_phi_fn(self):
        return self.phi_fn_class(mfn=self.damage_function)

    #----------------------------------------------------------------------------------
    # material properties stored in 'ccs'
    #----------------------------------------------------------------------------------

    # concrete matrix E-modulus (taken from 'ccs_unit_cell')
    #
    E_m = Property(Float, depends_on='input_change')
    @cached_property
    def _get_E_m(self):
        E_m = self.ccs_unit_cell_ref.get_E_m_time(self.age)
        print 'E_m (from ccs)', E_m
        return E_m

    # composite E-modulus (taken from 'ccs_unit_cell')
    #
    E_c = Property(Float, depends_on='input_change')
    @cached_property
    def _get_E_c(self):
        E_c = self.ccs_unit_cell_ref.get_E_c_time(self.age)
        print 'E_c (from ccs)', E_c
        return E_c

    # Poisson's ratio 
    #
    nu = Property(Float, depends_on='input_change')
    @cached_property
    def _get_nu(self):
        nu = self.ccs_unit_cell_ref.nu
        print 'nu (from ccs)', nu
        # set nu explicitly corresponding to settings in 'mats_calib_damage_fn'
        #
        print 'nu set explicitly to 0.20'
        nu = 0.2
        return nu

if __name__ == '__main__':

    sim_model = SimBT4PTDB(ccs_unit_cell_key='FIL-10-09_2D-05-11_0.00462_all0',
                                calibration_test='TT-12c-6cm-0-TU-SH2F-V3',
                                age=23)

    do = 'ui'
#    do = 'validation'
#    do = 'show_last_results'

    if do == 'ui':
        sim_model.tloop.eval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource=sim_model)
        app.main()

    if do == 'param_study':

        # influence of the calibration test
        #
        param_list = ['TT-12c-6cm-TU-SH1F-V1', 'TT-12c-6cm-0-TU-SH2F-V2', 'TT-12c-6cm-0-TU-SH2F-V3']
        for param_key in param_list:
            run_key = 'f_w_diagram_x_olmyz-22211_Ec-28600_nu-025_tol-m_nsteps-20_' + param_key
            sim_model = SimFourPointBendingDB(ccs_unit_cell_key='FIL-10-09_2D-05-11_0.00462_all0',
                                              calibration_test=param_key,
#                                              outer_zone_shape_x = 2,
#                                              load_zone_shape_x = 2,
#                                              mid_zone_shape_x = 2,
#                                              shape_y = 1,
#                                              shape_z = 1,
                                              age=26)
            sim_model.run_study(run_key)

        # influence of the z-discretization
        #
        param_list = ['1', '2', '3', '4']
        for param_key in param_list:
            shape_z = int(param_key)
            run_key = 'f_w_diagram_x_olmyz-2221' + param_key + '_Ec-28600_nu-025_tol-m_nsteps-20_TT-12c-6cm-TU-SH2F-V3'  # discretization fineness
            sim_model = SimFourPointBendingDB(ccs_unit_cell_key='FIL-10-09_2D-05-11_0.00462_all0',
                                              calibration_test='TT-12c-6cm-0-TU-SH2F-V3',
#                                              outer_zone_shape_x = 2,
#                                              load_zone_shape_x = 2,
#                                              mid_zone_shape_x = 2,
#                                              shape_y = 1,
                                              shape_z=shape_z,
                                              age=28)
            sim_model.run_study(run_key)

    if do == 'validation':

        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        #------------------------------
        # simulation:
        #------------------------------

        sim_model.tloop.eval()

        pickle_path = 'pickle_files'
        png_path = 'png_files'
        param_key = 'TT-12c-6cm-TU-SH2F-V3'

#        # F-w-diagram_center_elasttop
#        #
#        sim_model.f_w_diagram_center_elasttop.refresh()
#        file_name = 'f_w_diagram_c_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'w')
#        dump(sim_model.f_w_diagram_center_elasttop.trace, file)
#        file.close()
#        sim_model.f_w_diagram_center_elasttop.trace.mpl_plot(p, color = 'red')

        # F-w-diagram_center_elasttop (ELAST TOP)
        #
#        sim_model.f_w_diagram_center_elasttop.refresh()
#        sim_model.f_w_diagram_center_elasttop.trace.mpl_plot(p, color = 'red', linestyle = '-', label = 'ET', linewidth = 2.)

        # F-w-diagram_center_elastbottom (ELAST BOTTOM)
        #
#        sim_model.f_w_diagram_center_elastbottom.refresh()
#        sim_model.f_w_diagram_center_elastbottom.trace.mpl_plot(p, color = 'blue', linestyle = '-', label = 'EB')

        # F-w-diagram_center_spectop (TOP SURFACE)
        #
#        sim_model.f_w_diagram_center_spectop.refresh()
#        sim_model.f_w_diagram_center_spectop.trace.mpl_plot(p, color = 'green', linestyle = '-', label = 'ST')

        # F-w-diagram_center_supprtline (SUPPRT LINE)
        #
        sim_model.f_w_diagram_center.refresh()
        sim_model.f_w_diagram_center.trace.mpl_plot(p, color='black', linestyle='-', label='SL')

        # F-w-diagram_center_topline (TOP LINE)
        #
#        sim_model.f_w_diagram_center_topline.refresh()
#        sim_model.f_w_diagram_center_topline.trace.mpl_plot(p, color = 'red', linestyle = '--')

#        # F-w-diagram_load (thirdpoint of the beam, e.g. loading point)
#        #
#        sim_model.f_w_diagram_load.refresh()
#        file_name = 'f_w_diagram_load_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'w')
#        dump(sim_model.f_w_diagram_load.trace, file)
#        file.close()
#        sim_model.f_w_diagram_load.trace.mpl_plot(p, color = 'blue')

        #------------------------------
        # experiment:
        #------------------------------
        path = join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU', 'BT-4PT-12c-6cm-SH4')
        tests = [ 'BT-4PT-12c-6cm-SH4-V1.DAT']  # , 'BT-4PT-12c-6cm-SH4-V2.DAT' ]
        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
#            ex_run.ex_type._plot_smoothed_force_deflection_center( p )
            ex_run.ex_type._plot_force_deflection_center(p)
            ex_run.ex_type._plot_force_deflection_thirdpoints(p)

#        format_plot(p, xlabel = 'displacement [mm]', ylabel = 'applied force [kN]', xlim = 50., ylim = 20.)
#        run_key = 'BT-4PT_' + param_key
#        png_file_path = join( png_path, run_key )
#        p.title( run_key )
#        p.savefig( run_key, dpi=600.)
        p.show()

    if do == 'show_last_results':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = 'pickle_files'
        param_key = 'SH2F-V2_nelems14-3-2_mtol_t0-01_w30mm_PhiExp0-20'

        #------------------------------
        # simulation:
        #------------------------------

        # f-w-diagram_center
        #
#        file_name = 'f_w_diagram_c_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'r')
#        trace = load(file)
#        p.plot(trace.xdata, trace.ydata, color = 'red')

        # f-w-diagram_center-edge
        #
#        file_name = 'f_w_diagram_ce_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'r')
#        trace = load(file)
#        p.plot(trace.xdata, trace.ydata, color = 'red')

        #------------------------------
        # experiment:
        #------------------------------
        path = join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU', 'BT-4PT-12c-6cm-SH4')
        tests = [ 'BT-4PT-12c-6cm-SH4-V1.DAT', 'BT-4PT-12c-6cm-SH4-V2.DAT' ]
        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
#            ex_run.ex_type._plot_smoothed_force_deflection_center( p )
            ex_run.ex_type._plot_ironed_orig_force_deflection_center(p)

        # plot sim curve as time new roman black and white plot 
        #
#        format_plot(p, xlim = 34, ylim = 54, xlabel = 'displacement [mm]', ylabel = 'force [kN]')

        p.show()

    if do == 'pstudy':
        sim_ps = SimPStudy(sim_model=sim_model)
        sim_ps.configure_traits()


