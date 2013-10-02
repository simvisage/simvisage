
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
    IPhiFn, PhiFnGeneralExtended, PhiFnGeneralExtendedExp, \
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

class MATS3DElastomerZ(MATS3DElastic):
    '''
    3D material model for given elastomer with predefined stress-strain law.
    '''
    # Piece wise linear stress strain curve
    #
    _stress_strain_curve = Instance(MFnLineArray)
    def __stress_strain_curve_default(self):
        return MFnLineArray(ydata = [ 0., self.E ],
                            xdata = [ 0., 1.])

    @on_trait_change('E')
    def reset_stress_strain_curve(self):
        self._stress_strain_curve = MFnLineArray(ydata = [ 0., self.E ],
                                                 xdata = [ 0., 1.])

    stress_strain_curve = Property
    def _get_stress_strain_curve(self):
        return self._stress_strain_curve

    def _set_stress_strain_curve(self, curve):
        self._stress_strain_curve = curve

    #-----------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------

    def get_corr_pred(self, sctx, eps_app_eng, d_eps, tn, tn1):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        eps_z = eps_app_eng[2]
        E = array([[self.stress_strain_curve.get_diff(eps_z)]])
        D_el = self.get_D_el(E, tn1)
        sigma = np.dot(D_el, eps_app_eng)
        return  sigma, self.D_el


class SimBT3PT(IBVModel):
    '''Simulation: Bending Test Three Point 
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
    mid_shape_x = Int(2, input = True,
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
    length = Float(1.15, input = True)
    elstmr_length = Float(0.05, input = True)
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

    # half the length of the elastomer (load introduction
    # with included symmetry
    #
    sym_elstmr_length = Property
    def _get_sym_elstmr_length(self):
        return (self.elstmr_length) / 2.

    # half the specimen width
    #
    sym_width = Property
    def _get_sym_width(self):
        return self.width / 2.

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
    E_c = Float(28600., input = True)

    # Poisson's ratio 
    #
    nu = Float(0.25, input = True)

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

    mid_specmn_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_mid_specmn_fe_level(self):
        return  FERefinementGrid(name = 'middle specimen patch',
                                 fets_eval = self.specmn_fets,
                                 domain = self.fe_domain)

    elstmr_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_elstmr_fe_level(self):
        return  FERefinementGrid(name = 'elastomer patch',
                                 fets_eval = self.elstmr_fets,
                                 domain = self.fe_domain)

    mid_specmn_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_mid_specmn_fe_grid(self):
        # only a quarter of the beam is simulated due to symmetry:
        fe_grid = FEGrid(coord_max = (self.sym_elstmr_length,
                                      self.width / 2, self.thickness),
                         shape = (self.mid_shape_x, self.shape_y, self.shape_z),
                         level = self.specmn_fe_level,
                         fets_eval = self.specmn_fets)
        return fe_grid
    
    specmn_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_specmn_fe_grid(self):
        # only a quarter of the beam is simulated due to symmetry:
        fe_grid = FEGrid(coord_min = (self.sym_elstmr_length, 
                                      0., 
                                      0.),
                         coord_max = (self.sym_specmn_length,
                                      self.sym_width,
                                      self.thickness),
                         shape = (self.shape_x, self.shape_y, self.shape_z),
                         level = self.specmn_fe_level,
                         fets_eval = self.specmn_fets)
        return fe_grid

    elstmr_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_elstmr_fe_grid(self):
        x_max = self.sym_elstmr_length
        y_max = self.width / 2.
        z_max = self.thickness + self.elstmr_thickness
        fe_grid = FEGrid(coord_min = (0, 0, self.thickness),
                         coord_max = (x_max, y_max, z_max),
                         level = self.elstmr_fe_level,
                         shape = (self.mid_shape_x, self.shape_y, 1),
                         fets_eval = self.elstmr_fets)
        return fe_grid

    #===========================================================================
    # Boundary conditions
    #===========================================================================
    bc_list = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_list(self):
        specimen = self.specmn_fe_grid
        mid_specimen = self.mid_specmn_fe_grid
        elastomer = self.elstmr_fe_grid

        #--------------------------------------------------------------
        # boundary conditions for the symmetry
        #--------------------------------------------------------------
        # the x-axis corresponds to the axis of symmetry along the longitudinal axis of the beam:
        bc_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1],
                                 slice = specimen[:, 0, :, :, 0, :])
        
        bc_mid_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1],
                                 slice = mid_specimen[:, 0, :, :, 0, :])
        
        bc_mid_symplane_yz = BCSlice(var = 'u', value = 0., dims = [0],
                                 slice = mid_specimen[0, :, :, 0, :, :])
        
        bc_el_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1],
                                 slice = elastomer[:, 0, :, :, 0, :])
        
        bc_el_symplane_yz = BCSlice(var = 'u', value = 0., dims = [0],
                                 slice = elastomer[0, :, :, 0, :, :])
        #--------------------------------------------------------------
        # boundary conditions for the support
        #--------------------------------------------------------------
        bc_support_0y0 = BCSlice(var = 'u', value = 0., dims = [2],
                                 slice = specimen[-1, :, 0, -1, :, 0])

        #--------------------------------------------------------------
        # connect elastomer
        #--------------------------------------------------------------
        link_sp_sp = BCDofGroup(var = 'u', value = 0., dims = [0, 1, 2],
                                get_dof_method = mid_specimen.get_right_dofs,
                                get_link_dof_method = specimen.get_left_dofs,
                                link_coeffs = [1.])

        link_el_sp = BCDofGroup(var = 'u', value = 0., dims = [2],
                                get_dof_method = elastomer.get_back_dofs,
                                get_link_dof_method = mid_specimen.get_front_dofs,
                                link_coeffs = [1.])

        #--------------------------------------------------------------
        # loading
        #--------------------------------------------------------------
        # w_max = center displacement:
        w_max = -0.010 # [m]
#        f_max = -0.010 / 0.10 # [MN/m]

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
                bc_mid_symplane_xz,
                bc_mid_symplane_yz,
                bc_support_0y0,
                bc_center_w,
                link_el_sp,
                link_sp_sp,
                bc_el_symplane_xz,
                bc_el_symplane_yz
                ]

    tstep = Float(0.05, auto_set = False, enter_set = True, input = True)
    tmax = Float(1.0, auto_set = False, enter_set = True, input = True)

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
        center_dof = self.ctl_dofs[0]

        # force-displacement-diagram (ELASTOMER TOP)
        # 
        self.f_w_diagram_center = RTraceGraph(name = 'displacement (center) - reaction 2',
                                       var_x = 'U_k'  , idx_x = center_dof,

                                       # nodal displacement at center node
                                       #
#                                       var_y = 'F_int', idx_y = center_dof,

#                                       # surf load
#                                       #
                                       var_y = 'F_int', idx_y_arr = self.ctl_dofs,

                                       record_on = 'update',
                                       transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the beam (2 symmetry axis):
                                       #
                                       transform_y = '-4000. * y')

        # force-displacement-diagram (SUPPORT)
        # 
        # SUPPRT LINE
        # dofs at support line of the specmn (used to integrate the force)
        #
        spec_supprtline_dofs_z = specimen[-1, :, 0, -1, :, 0].dofs[:, :, 2].flatten()
        spec_supprtline_dofs_z = np.unique( spec_supprtline_dofs_z )
        print 'spec_supprtline_dofs_z (unique)', spec_supprtline_dofs_z

        self.f_w_diagram_supprt = RTraceGraph(name = 'displacement (center) - reaction 2',
                                       var_x = 'U_k'  , idx_x = center_dof,

                                       # support line
                                       #
                                       var_y = 'F_int', idx_y_arr = spec_supprtline_dofs_z,

                                       record_on = 'update',
                                       transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the beam (2 symmetry axis):
                                       #
                                       transform_y = '4000. * y')

        ts = TS(
                sdomain = self.fe_domain,
                bcond_list = self.bc_list,
                rtrace_list = [
                             self.f_w_diagram_center,
                             self.f_w_diagram_supprt,
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
                             RTraceDomainListField(name = 'max_omega_i', warp = True,
                                        var = 'max_omega_i', idx = 0,
                                        record_on = 'update'),
                             RTraceDomainListField(name = 'IStress' ,
                                            position = 'int_pnts',
                                            var = 'sig_app', idx = 0,
                                            record_on = 'update'),
                             RTraceDomainListField(name = 'IStrain' ,
                                            position = 'int_pnts',
                                            var = 'eps_app', idx = 0,
                                            record_on = 'update'),
                             RTraceDomainListField(name = 'max_omega_i' ,
                                        var = 'max_omega_i', idx = 0,
                                        record_on = 'update'),                              ]
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
#                       tline = TLine( min = 0.0, step = 0.05, max = 0.05 )
                       tline = TLine(min = 0.0, step = self.tstep, max = self.tmax)
                       )

        return tloop

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name = 'u_center_top_z', unit = 'm'),
                 SimOut(name = 'F_max', unit = 'kN') ]


class SimBT3PTDB(SimBT3PT):
    '''Simulation: Bending Test Three Point Data Base
    '''
    # vary the failure strain in PhiFnGeneralExtended:
    factor_eps_fail = Float(1.0, input = True,
                             ps_levels = (1.0, 1.2, 3))

    Efp_frac = Float(0.2, input = True)

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
        print 'XXX', self.ccs_unit_cell_ref.damage_function_list[0].calibration_test
        print 'XXX', self.ccs_unit_cell_ref.damage_function_list[0]
        print 'XXX', len( self.ccs_unit_cell_ref.damage_function_list )
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
#        return PhiFnGeneralExtended(mfn = self.damage_function,
#                                     factor_eps_fail = self.factor_eps_fail)
        return PhiFnGeneralExtendedExp(mfn = self.damage_function,
                                       Efp_frac = self.Efp_frac )

    #----------------------------------------------------------------------------------
    # mats_eval
    #----------------------------------------------------------------------------------
    # age of the beam at the time of testing
    # NOTE: that the same phi-function is used independent of age. This assumes a 
    # an afine/proportional damage evolution for different ages. 
    #
    age = Int(28, input = True )

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

    sim_model = SimBT3PTDB(
                           
#                           ccs_unit_cell_key = 'FIL-10-09_2D-05-11_0.00462_all0',
##                           calibration_test = 'TT-12c-6cm-0-TU-SH2F-V2',
#                           calibration_test = 'TT-12c-6cm-0-TU-SH2F-V3',
#                           thickness = 0.06,
#                           length = 1.15,
                           
                           ccs_unit_cell_key = 'barrelshell_2D-05-11_0.00286_all0',
                           calibration_test = 'TTb-6c-2cm-0-TU-V1_bs5_a23d-nu02',
                           thickness = 0.02,
                           length = 0.46,
                           width = 0.10,
                           
                           # calibration performed with:
                           # 'n_mp' = 30; 
                           # 'E_c' = 28600; 
                           # 'nu' = 0.25; 
                           # 'compliance_version'

                           tstep = 0.05,
                           tmax = 1.0,

#                           calibration_test = 'TT-12c-6cm-TU-SH1F-V1',
                           age = 28 )
#                                age = 9)

#    sim_model = SimTreePointBendingDB( ccs_unit_cell_key = 'FIL-10-09_2D-02-06a_0.00273_90_0',
##                                calibration_test = 'TT11-10a-average',
#                                calibration_test = 'TT11-10a-V2',
#                                age = 28 )



#    do = 'ui'
    do = 'validation'
#    do = 'pstudy'
#    do = 'show_last_results'

    if do == 'ui':
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource = sim_model)
#        sim_model.tloop.eval()
        app.main()

    if do == 'pstudy':
        sim_ps = SimPStudy(sim_model = sim_model)
        sim_ps.configure_traits()

    if do == 'validation':

        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        # calculate F-w-curve and store as pickle file
        #
        sim_model.tloop.eval()

        pickle_path = 'pickle_files'

        ccs_unit_cell_key = sim_model.ccs_unit_cell_key
        calibration_test = sim_model.calibration_test
        L = sim_model.length
        w = sim_model.width
        t = sim_model.thickness
        shape_z = sim_model.shape_z

        param_key = ccs_unit_cell_key + '_' + calibration_test +'Lwt-%g%g%g_shapez-%g' %(L, w, t, shape_z)
        
        # F-w (simulation)
        # ELAST SURF
        #
        file_name = 'f_w_diagram_c_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'w')
        sim_model.f_w_diagram_center.refresh()
        dump(sim_model.f_w_diagram_center.trace, file)
        file.close()
        sim_model.f_w_diagram_center.trace.mpl_plot(p, color = 'red')
        
        # SUPPRT LINE
        #
        file_name = 'f_w_diagram_supprt_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'w')
        sim_model.f_w_diagram_supprt.refresh()
        dump(sim_model.f_w_diagram_supprt.trace, file)
        file.close()
        sim_model.f_w_diagram_supprt.trace.mpl_plot(p, color = 'blue')

        # F-w (experiment)
        #
        path = join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE')
        tests = [ 
                  'BT-3PT-12c-6cm-0-Tu-V1.raw',
                  'BT-3PT-12c-6cm-0-Tu-V2.raw', 
                  'BT-3PT-12c-6cm-0-Tu-V3.raw',
                  'BT-3PT-12c-6cm-0-Tu-V4.raw' ]

        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
            ex_run.ex_type._plot_force_machine_displacement_wo_elast_interpolated(p)
#            ex_run.ex_type._plot_force_machine_displacement_wo_elast(p)
#            ex_run.ex_type._plot_force_machine_displacement(p)

        p.show()

    if do == 'show_last_results':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = 'pickle_files'

        param_key = 'PhiFn-SH2F-V2'

        # f-w-diagram_center
        #
        sim_model.f_w_diagram_center.refresh()
        file_name = 'f_w_diagram_c_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color = 'blue')

        # F-w-curves (experiment)
        #
        path = join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE')
        tests = [ 
                  'BT-3PT-12c-6cm-0-Tu-V1.raw',
#                  'BT-3PT-12c-6cm-0-Tu-V2.raw', 
#                  'BT-3PT-12c-6cm-0-Tu-V3.raw',
                  'BT-3PT-12c-6cm-0-Tu-V4.raw' ]

        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
            ex_run.ex_type._plot_force_machine_displacement_wo_elast(p)
#            ex_run.ex_type._plot_force_machine_displacement_wo_elast_interpolated(p)
#            ex_run.ex_type._plot_force_machine_displacement(p)

        p.show()

        ### shape ###:
        ### eps_factor ###:
        ### average - V1 (c) ###:
        ### c ce e - average ###:
        ### c ce e - V1 ###:
        ### step ###:
        ### tolerance ###:
