
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

from ibvpy.mesh.fe_grid import \
    FEGrid

from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
    MATS2D5MicroplaneDamage

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
    IPhiFn, PhiFnGeneralExtended, \
    PhiFnGeneral, PhiFnStrainHardening, PhiFnStrainHardeningLinear

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

simdb = SimDB()

from pickle import dump, load


class SimFourPointBending(IBVModel):
    '''Plate test prepared for parametric study.
    '''

    input_change = Event
    @on_trait_change('+input,ccs_unit_cell.input_change')
    def _set_input_change(self):
        self.input_change = True

    implements(ISimModel)

    #-----------------
    # discretization:
    #-----------------
    #
    # discretization in x-direction (longitudinal):
    shape_x = Int(8, input = True,
                      ps_levels = (4, 12, 3))

    # discretization in x-direction (longitudinal):
    load_zone_shape_x = Int(2, input = True,
                      ps_levels = (1, 4, 1))

    # middle part discretization in x-direction (longitudinal):
    mid_zone_shape_x = Int(4, input = True,
                      ps_levels = (1, 4, 1))

    # discretization in y-direction (width):
    shape_y = Int(2, input = True,
                      ps_levels = (1, 4, 2))

    # discretization in z-direction:
    shape_z = Int(3, input = True,
                      ps_levels = (1, 3, 3))

    #-----------------
    # geometry:
    #-----------------
    #
    # edge length of the bending specimen (beam) (entire length without symmetry)
    length = Float(1.5, input = True)
    elstmr_length = Float(0.05, input = True)
    mid_zone_length = Float(0.45, input = True)
    elstmr_thickness = Float(0.005, input = True,
                                enter_set = True, auto_set = False)
    width = Float(0.20, input = True)
    thickness = Float(0.06, input = True)

    #-----------------
    # derived geometric parameters
    #-----------------
    #
    # half the length of the elastomer (load introduction 
    # with included symmetry)
    #
    sym_specmn_length = Property
    def _get_sym_specmn_length(self):
        return (self.length - self.elstmr_length) / 2.

    sym_mid_zone_specmn_length = Property
    def _get_sym_mid_zone_specmn_length(self):
        return (self.mid_zone_length) / 2.

    # half the length of the elastomer (load introduction
    # with included symmetry
    #
    sym_elstmr_length = Property
    def _get_sym_elstmr_length(self):
        return (self.elstmr_length) / 2.

    #-----------------
    # phi function extended:
    #-----------------
    #
    phi_fn = Instance(IPhiFn, input = True)
    def _phi_fn_default(self):
        return PhiFnStrainHardening()

    #----------------------------------------------------------------------------------
    # mats_eval
    #----------------------------------------------------------------------------------

    # age of the plate at the time of testing
    # NOTE: that the same phi-function is used independent of age. This assumes a 
    # an afine/proportional damage evolution for different ages. 
    #
    age = Int(28, input = True)

    # composite E-modulus 
    #
    E_c = Float(28e5, input = True)

    # Poisson's ratio 
    #
    nu = Float(0.2, input = True)

    # @todo: for mats_eval the information of the unit cell should be used
    # in order to use the same number of microplanes and model version etc...
    #
    specmn_mats = Property(Instance(MATS2D5MicroplaneDamage),
                          depends_on = 'input_change')
    @cached_property
    def _get_specmn_mats(self):

    # used for basic testing:
#        return MATS3DElastic( 
#                                E = self.E_c,
#                                nu = self.nu )
        return MATS2D5MicroplaneDamage(
                                E = self.E_c,
                                nu = self.nu,

                                # corresponding to settings in "MatsCalib"
                                n_mp = 30,
                                symmetrization = 'sum-type',
                                model_version = 'compliance',
                                phi_fn = self.phi_fn)

    elstmr_mats = Property(Instance(MATS3DElastic),
                                   depends_on = 'input_change')
    @cached_property
    def _get_elstmr_mats(self):

        max_eps = self.elstmr_thickness
        max_f = 0.020 # MN
        max_eps = 1.0 # [-]
        area = self.elstmr_length * self.width / 2.
        sig = max_f / area
        E_elast = sig / max_eps
        print 'effective elastomer E_modulus', E_elast
        return MATS3DElastic(E = E_elast,
                              nu = 0.4)

    #-----------------
    # fets:
    #-----------------

    # use quadratic serendipity elements
    #
    specmn_fets = Property(Instance(FETSEval),
                             depends_on = 'input_change')
    @cached_property
    def _get_specmn_fets(self):
        return FETS2D58H20U(mats_eval = self.specmn_mats)

    # use quadratic serendipity elements
    #
    elstmr_fets = Property(Instance(FETSEval),
                             depends_on = 'input_change')
    @cached_property
    def _get_elstmr_fets(self):
        return FETS2D58H20U(mats_eval = self.elstmr_mats)

# used for basic testing:
#        return FETS3D8H20U( mats_eval = self.mats_eval )
#        return FETS3D8H( mats_eval = self.mats_eval )

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
                        dtype = 'float_')

    fe_domain = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_domain(self):
        return FEDomain()

    specmn_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_specmn_fe_level(self):
        return  FERefinementGrid(name = 'specimen patch',
                                 fets_eval = self.specmn_fets,
                                 domain = self.fe_domain)

    load_zone_specmn_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_load_zone_specmn_fe_level(self):
        return  FERefinementGrid(name = 'load_zonedle specimen patch',
                                 fets_eval = self.specmn_fets,
                                 domain = self.fe_domain)

    mid_zone_specmn_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_mid_zone_specmn_fe_level(self):
        return  FERefinementGrid(name = 'load_zonedle specimen patch',
                                 fets_eval = self.specmn_fets,
                                 domain = self.fe_domain)

    elstmr_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_elstmr_fe_level(self):
        return  FERefinementGrid(name = 'elastomer patch',
                                 fets_eval = self.elstmr_fets,
                                 domain = self.fe_domain)

    #===========================================================================
    # Grid definition 
    #===========================================================================
    
    mid_zone_specmn_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_mid_zone_specmn_fe_grid(self):
        # only a quarter of the plate is simulated due to symmetry:
        fe_grid = FEGrid(coord_max = (self.sym_mid_zone_specmn_length,
                                      self.width / 2, self.thickness),
                         shape = (self.mid_zone_shape_x, self.shape_y, self.shape_z),
                         level = self.specmn_fe_level,
                         fets_eval = self.specmn_fets)

        return fe_grid

    load_zone_specmn_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_load_zone_specmn_fe_grid(self):
        # only a quarter of the plate is simulated due to symmetry:
        fe_grid = FEGrid(coord_min = (self.sym_mid_zone_specmn_length,
                                      0, 0),
                         coord_max = (self.sym_mid_zone_specmn_length + self.sym_elstmr_length,
                                      self.width / 2, self.thickness),
                         shape = (self.load_zone_shape_x, self.shape_y, self.shape_z),
                         level = self.specmn_fe_level,
                         fets_eval = self.specmn_fets)

        return fe_grid

    elstmr_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_elstmr_fe_grid(self):
        x_max = self.sym_mid_zone_specmn_length + self.sym_elstmr_length
        y_max = self.width / 2.
        z_max = self.thickness + self.elstmr_thickness
        fe_grid = FEGrid(coord_min = (self.sym_mid_zone_specmn_length,
                                      0, self.thickness),
                         coord_max = (x_max, y_max, z_max),
                         level = self.elstmr_fe_level,
                         shape = (self.load_zone_shape_x, self.shape_y, 1),
                         fets_eval = self.elstmr_fets)
        return fe_grid

    specmn_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_specmn_fe_grid(self):
        # only a quarter of the plate is simulated due to symmetry:
        fe_grid = FEGrid(coord_min = (self.sym_mid_zone_specmn_length + self.sym_elstmr_length,
                                      0, 0),
                         coord_max = (self.sym_specmn_length,
                                      self.width / 2,
                                      self.thickness),
                         shape = (self.shape_x, self.shape_y, self.shape_z),
                         level = self.specmn_fe_level,
                         fets_eval = self.specmn_fets)

        return fe_grid

    #===========================================================================
    # Boundary conditions
    #===========================================================================
    bc_list = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_list(self):
        specimen = self.specmn_fe_grid
        load_zone_specimen = self.load_zone_specmn_fe_grid
        elastomer = self.elstmr_fe_grid
        mid_zone_specimen = self.mid_zone_specmn_fe_grid

        #--------------------------------------------------------------
        # boundary conditions for the symmetry and the single support
        #--------------------------------------------------------------
        # the x-axis correspons to the axis of symmetry along the longitudinal axis of the beam:
        bc_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1],
                                 slice = specimen[:, 0, :, :, 0, :])
        bc_load_zone_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1],
                                 slice = load_zone_specimen[:, 0, :, :, 0, :])
#        bc_load_zone_symplane_yz = BCSlice(var = 'u', value = 0., dims = [0],
#                                 slice = load_zone_specimen[0, :, :, 0, :, :])
        bc_mid_zone_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1],
                                 slice = mid_zone_specimen[:, 0, :, :, 0, :])
        bc_mid_zone_symplane_yz = BCSlice(mid = 'u', value = 0., dims = [0],
                                 slice = mid_zone_specimen[0, :, :, 0, :, :])

        bc_support_0y0 = BCSlice(var = 'u', value = 0., dims = [2],
                                 slice = specimen[-1, :, 0, -1, :, 0])

        #--------------------------------------------------------------
        # connect elastomer
        #--------------------------------------------------------------
        link_loadzn_outerzn = BCDofGroup(var = 'u', value = 0., dims = [0, 1, 2],
                                get_dof_method = load_zone_specimen.get_right_dofs,
                                get_link_dof_method = specimen.get_left_dofs,
                                link_coeffs = [1.])

        bc_el_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1],
                                 slice = elastomer[:, 0, :, :, 0, :])
#        bc_el_symplane_yz = BCSlice(var = 'u', value = 0., dims = [0],
#                                 slice = elastomer[0, :, :, 0, :, :])
        link_elstmr_loadzn_x = BCDofGroup(var = 'u', value = 0., dims = [0],
                                get_dof_method = elastomer.get_back_dofs,
                                get_link_dof_method = load_zone_specimen.get_front_dofs,
                                link_coeffs = [1.])
        link_elstmr_loadzn_z = BCDofGroup(var = 'u', value = 0., dims = [2],
                                get_dof_method = elastomer.get_back_dofs,
                                get_link_dof_method = load_zone_specimen.get_front_dofs,
                                link_coeffs = [1.])

        link_midzn_loadzn = BCDofGroup(var = 'u', value = 0., dims = [0, 1, 2],
                                get_dof_method = mid_zone_specimen.get_right_dofs,
                                get_link_dof_method = load_zone_specimen.get_left_dofs,
                                link_coeffs = [1.])

        #--------------------------------------------------------------
        # loading
        #--------------------------------------------------------------
        # w_max = center displacement:
        w_max = -0.030 # [m]
        f_max = -0.010 / 0.10 # [MN/m]

        #--------------------------------------------------------------
        # nodal displacement applied at top at the center of the plate
        #--------------------------------------------------------------
        # NOTE: the entire symmetry axis (yz)-plane is moved downwards 
        # in order to avoid large indentations at the top nodes
        #
#        bc_center_w = BCSlice( var = 'u', value = w_max, dims = [2], slice = domain[-1, :, :, -1, :, :] )
#        bc_center_f = BCSlice( var = 'f', value = 1.0, dims = [2], slice = domain[-1, :, :, -1, :, :] )

        # apply displacement at all center node (line load)
        #
        bc_center_w = BCSlice(var = 'u', value = w_max, dims = [2],
                              slice = elastomer[:, :, -1, :, :, -1])

        return [bc_symplane_xz,
                link_midzn_loadzn,
                bc_load_zone_symplane_xz,
#                bc_load_zone_symplane_yz,
                bc_support_0y0,
                bc_center_w,
                link_elstmr_loadzn,
                link_loadzn_outerzn,
                bc_mid_zone_symplane_xz,
                bc_mid_zone_symplane_yz,
                bc_el_symplane_xz,
                bc_el_symplane_yz
#                              bc_center_f,
                ]

    tloop = Property(depends_on = 'input_change')
    @cached_property
    def _get_tloop(self):

        #--------------------------------------------------------------
        # ts 
        #--------------------------------------------------------------

        specimen = self.specmn_fe_grid
        elastomer = self.elstmr_fe_grid

        # center_top_line_dofs
        #
        ctl_dofs = elastomer[:, :, -1, :, :, -1].dofs[:, :, 2].flatten()

#        # right_bottom_line_dofs
#        #
#        ctl_dofs = domain[0, :, 0, 0, :, 0].dofs[:, :, 2].flatten()
#        print 'ctl_dofs', ctl_dofs

        self.ctl_dofs = np.unique(ctl_dofs)

        # center_dof (used for tracing of the displacement)
        #
        mid_zone_spec = self.mid_zone_specmn_fe_grid
        
        underneath_load_dof = mid_zone_spec[-1, 0, 0, -1, 0, 0].dofs[0, 0, 2] 
        center_bottom_dof = mid_zone_spec[0, 0, 0, 0, 0, 0].dofs[0, 0, 2]

        # force-displacement-diagram 
        # 
        self.f_w_diagram_center = RTraceGraph(name = 'displacement (center) - reaction 2',
                                       var_x = 'U_k'  , idx_x = center_bottom_dof,

                                       # nodal displacement at center node
                                       #
#                                       var_y = 'F_int', idx_y = center_dof,

#                                       # line load
#                                       #
                                       var_y = 'F_int', idx_y_arr = self.ctl_dofs,

                                       record_on = 'update',
                                       transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the beam (2 symmetry axis):
                                       #
                                       transform_y = '-4000. * y')

        # force-displacement-diagram 
        # 
        self.f_w_diagram_load = RTraceGraph(name = 'displacement (center) - reaction 2',
                                       var_x = 'U_k'  , idx_x = underneath_load_dof,

                                       # nodal displacement at center node
                                       #
#                                       var_y = 'F_int', idx_y = center_dof,

#                                       # line load
#                                       #
                                       var_y = 'F_int', idx_y_arr = self.ctl_dofs,

                                       record_on = 'update',
                                       transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the beam (2 symmetry axis):
                                       #
                                       transform_y = '-4000. * y')

        ts = TS(
                sdomain = self.fe_domain,
                bcond_list = self.bc_list,
                rtrace_list = [
                             self.f_w_diagram_center,
                             self.f_w_diagram_load,
                             RTraceDomainListField(name = 'Displacement' ,
                                            var = 'u', idx = 0, warp = True),
                             RTraceDomainListField(name = 'Stress' ,
                                            var = 'sig_app', idx = 0, warp = True,
                                            record_on = 'update'),
                             RTraceDomainListField(name = 'Strain' ,
                                        var = 'eps_app', idx = 0, warp = True,
                                        record_on = 'update'),
                             RTraceDomainListField(name = 'Damage' ,
                                        var = 'omega_mtx', idx = 0, warp = True,
                                        record_on = 'update'),
                             RTraceDomainListField(name = 'IStress' ,
                                            position = 'int_pnts',
                                            var = 'sig_app', idx = 0,
                                            record_on = 'update'),
                             RTraceDomainListField(name = 'IStrain' ,
                                            position = 'int_pnts',
                                            var = 'eps_app', idx = 0,
                                            record_on = 'update'),
                              ]
                )

        # Add the time-loop control
        tloop = TLoop(tstepper = ts,

                       # allow only a low tolerance 
                       #
                       KMAX = 50,
                       tolerance = 5e-4,

#                       # allow a high tolerance 
#                       #
#                       KMAX = 50,
#                       tolerance = 0.001,

                       RESETMAX = 0,
                       debug = False,
#                       tline = TLine(min = 0.0, step = 0.05, max = 0.05)
                       tline = TLine(min = 0.0, step = 0.1, max = 1.0)
                       )

        return tloop

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name = 'u_center_top_z', unit = 'm'),
                 SimOut(name = 'F_max', unit = 'kN') ]


class SimFourPointBendingDB(SimFourPointBending):

    # vary the failure strain in PhiFnGeneralExtended:
    factor_eps_fail = Float(1.0, input = True,
                             ps_levels = (1.0, 1.2, 3))

    #-----------------
    # composite cross section unit cell:
    #-----------------
    #
    ccs_unit_cell_key = Enum(CCSUnitCell.db.keys(),
                              simdb = True, input = True,
                              auto_set = False, enter_set = True)

    ccs_unit_cell_ref = Property(Instance(SimDBClass),
                                  depends_on = 'ccs_unit_cell_key')
    @cached_property
    def _get_ccs_unit_cell_ref(self):
        return CCSUnitCell.db[ self.ccs_unit_cell_key ]

    #-----------------
    # damage function:
    #-----------------
    #
    material_model = Str(input = True)
    def _material_model_default(self):
        # return the material model key of the first DamageFunctionEntry
        # This is necessary to avoid an ValueError at setup  
        return self.ccs_unit_cell_ref.damage_function_list[0].material_model

    calibration_test = Str(input = True
                            )
    def _calibration_test_default(self):
        # return the material model key of the first DamageFunctionEntry
        # This is necessary to avoid an ValueError at setup  
        return self.ccs_unit_cell_ref.damage_function_list[0].calibration_test

    damage_function = Property(Instance(MFnLineArray),
                                depends_on = 'input_change')
    @cached_property
    def _get_damage_function(self):
        return self.ccs_unit_cell_ref.get_param(self.material_model, self.calibration_test)

    #-----------------
    # phi function extended:
    #-----------------
    #
    phi_fn = Property(Instance(PhiFnGeneralExtended),
                       depends_on = 'input_change,+ps_levels')
    @cached_property
    def _get_phi_fn(self):
        return PhiFnGeneralExtended(mfn = self.damage_function,
                                     factor_eps_fail = self.factor_eps_fail)

    #----------------------------------------------------------------------------------
    # mats_eval
    #----------------------------------------------------------------------------------

    # age of the plate at the time of testing
    # NOTE: that the same phi-function is used independent of age. This assumes a 
    # an afine/proportional damage evolution for different ages. 
    #
    age = Int(28, #input = True
                )

    # composite E-modulus 
    #
    E_c = Property(Float, depends_on = 'input_change')
    @cached_property
    def _get_E_c(self):
        return self.ccs_unit_cell_ref.get_E_c_time(self.age)

    # Poisson's ratio 
    #
    nu = Property(Float, depends_on = 'input_change')
    @cached_property
    def _get_nu(self):
        return self.ccs_unit_cell_ref.nu


if __name__ == '__main__':

    sim_model = SimFourPointBendingDB(ccs_unit_cell_key = 'FIL-10-09_2D-05-11_0.00462_all0',
                                calibration_test = 'TT-12c-6cm-TU-SH1F-V1',
#                                age = 27 )
                                age = 26)

    do = 'ui'
#    do = 'pstudy'
#    do = 'validation'
#    do = 'show_last_results'

    if do == 'ui':
        sim_model.tloop.eval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource = sim_model)
        app.main()

    if do == 'pstudy':
        sim_ps = SimPStudy(sim_model = sim_model)
        sim_ps.configure_traits()

    if do == 'validation':

        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        # BT-12c-6cm-TU-0
        path = join(simdb.exdata_dir, 'bending_tests', 'ZiE_2011-06-08_BT-12c-6cm-0-TU')
        tests = [ 'BT-12c-6cm-0-Tu-V4.raw' ]

        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
            ex_run.ex_type._plot_smoothed_force_deflection_center(p)

        sim_model.tloop.eval()

        sim_model.f_w_diagram_center.refresh()
        pickle_path = 'pickle_files'
        file_name = 'f_w_diagram_bending_0-V1.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'w')
        dump(sim_model.f_w_diagram_center.trace, file)
        file.close()

        sim_model.f_w_diagram_center.trace.mpl_plot(p, color = 'red')

        p.show()

    if do == 'show_last_results':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = 'pickle_files'

        ### shape ###:

        # plate elem(12,12,4); w=0.08; step = 0.05; KMAX = 20, tol = 0.001
        file_name = 'f_w_diagram_c_average_step0-05_KMAX20_tol0-001_nelems12-12-4.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color = 'blue')

#        # plate elem(8,8,4); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_nelems8-8-4.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )

        # plate elem(10,10,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_nelems10-10-2.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color = 'black')

        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color = 'red')


        ### eps_factor ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001, factor_eps_fail = 1.4
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_epsfactor_1-4.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )
#        
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001, factor_eps_fail = 1.2
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_epsfactor_1-2.pickle'
#        pickle_file_path = join( pickle_path, file_name )
        file = open(pickle_file_path, 'r')
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )
#        
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )


        ### average - V1 (c) ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_V1_step0-1_KMAX20_tol0-001.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )


        ### c ce e - average ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_ce_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )
#        
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_e_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )


        ### c ce e - V1 ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_ce_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )
#        
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_e_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )


        ### step ###:

#        # plate elem(8,8,2); w=0.08; step = 0.02; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_V1_step0-02_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )


        ### tolerance ###:

#        # plate elem(8,8,2); w=0.08; step = 0.02; KMAX = 50, tol = 0.0005
#        file_name = 'f_w_diagram_V1_step0-02_KMAX50_tol0-0005.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )
#        
#        # plate elem(8,8,2); w=0.08; step = 0.02; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_V1_step0-02_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'yellow' )


        path = join(simdb.exdata_dir, 'plate_tests', 'PT-10a')
        tests = [ 'PT10-10a.DAT', 'PT11-10a.DAT' , 'PT12-10a.DAT' ]
#        tests = [ 'PT10-10a.DAT' ]

        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
            ex_run.ex_type._plot_smoothed_force_deflection_center(p)
#            ex_run.ex_type._plot_force_edge_deflection( p )
#            ex_run.ex_type._plot_force_center_edge_deflection( p )

        p.show()
