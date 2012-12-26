
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
    PhiFnGeneral, PhiFnStrainHardening, PhiFnStrainHardeningLinear, \
    PhiFnGeneralExtendedExp

from mathkit.geo.geo_ndgrid import \
    GeoNDGrid

from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import \
    MFnNDGrid, GridPoint

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

from numpy import \
    sin, cos, c_, arange, hstack, array, max, frompyfunc, linspace

import numpy as np

from time import time
from os.path import join

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from simiter.sim_pstudy import \
    SimPStudy, SimOut, ISimModel

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


class SimST(IBVModel):
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
    # discretization in x,y-direction:
    shape_xy = Int( 5, input = True,
                      ps_levels = ( 8, 12, 3 ) )

    # discretization in z-direction:
    shape_z = Int( 2, input = True,
                      ps_levels = ( 2, 3, 3 ) )

    #-----------------
    # geometry:
    #-----------------
    #
    # edge length of the quadratic plate (entire length without symmetry)
    # length = 1.25m - 2*0.05m supports
    length = Float(1.15, input = True)

    thickness = Float(0.06, input = True)

    elem_length_xy = Property(Float, depends_on = 'shape_xy, length')
    @cached_property
    def _get_elem_length_xy(self):
        return self.length / self.shape_xy

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
    age = Int(28, auto_set = False, enter_set = True, input = True)

    # composite E-modulus 
    # 28 GPa = 28000 [MPa] 
    E_c = Float(28e3, auto_set = False, enter_set = True, input = True)

    # composite E-modulus 
    # 210 GPa = 210000 [MPa] 
    E_s = Float(210e3, auto_set = False, enter_set = True, input = True)

    # Poisson's ratio 
    #
    nu = Float(0.2, auto_set = False, enter_set = True, input = True)

    tstep = Float(0.05, auto_set = False, enter_set = True, input = True)


    # @todo: for mats_eval the information of the unit cell should be used
    # in order to use the same number of microplanes and model version etc...
    #
    mats_eval = Property(Instance(MATS2D5MicroplaneDamage),
                          depends_on = 'input_change')
    @cached_property
    def _get_mats_eval(self):
        mats_eval = MATS2D5MicroplaneDamage(
                                E = self.E_c,
                                nu = self.nu,
                                n_mp = 30,
                                symmetrization = 'sum-type',
                                model_version = 'compliance',
                                phi_fn = self.phi_fn)

        return mats_eval

    elstmr_mats = Property(Instance(MATS3DElastic),
                                   depends_on = 'input_change')
    @cached_property
    def _get_elstmr_mats(self):

        max_eps = self.thickness * 1.1
        max_f = 0.020 # MN
        max_eps = 1.0 # [-]
        area = self.elem_length_xy ** 2
        sig = max_f / area
        E_elast = sig / max_eps
        print 'effective elastomer E_modulus', E_elast
        return MATS3DElastic(E = E_elast,
                              nu = 0.2)

    supprt_mats = Property(Instance(MATS3DElastic),
                                   depends_on = 'input_change')
    @cached_property
    def _get_supprt_mats(self):
        return MATS3DElastic(E = self.E_s,
                              nu = 0.2)


    #-----------------
    # fets:
    #-----------------

    # use quadratic serendipity elements
    #
    specmn_fets = Property(Instance(FETSEval),
                     depends_on = 'input_change')
    @cached_property
    def _get_specmn_fets(self):
        return FETS2D58H20U(mats_eval = self.mats_eval)

    # use quadratic serendipity elements
    #
    elstmr_fets = Property(Instance(FETSEval),
                             depends_on = 'input_change')
    @cached_property
    def _get_elstmr_fets(self):
        return FETS2D58H20U(mats_eval = self.elstmr_mats)

    # use quadratic serendipity elements
    #
    supprt_fets = Property(Instance(FETSEval),
                             depends_on = 'input_change')
    @cached_property
    def _get_supprt_fets(self):
        return FETS2D58H20U(mats_eval = self.supprt_mats)

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

    specmn_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_specmn_fe_grid(self):
        # only a quarter of the plate is simulated due to symmetry:
        # only a quarter of the plate is simulated due to symmetry:
        fe_grid = FEGrid(coord_max = (self.length / 2, self.length / 2, self.thickness),
                         shape = (self.shape_xy, self.shape_xy, self.shape_z),
                         level = self.specmn_fe_level,
                         fets_eval = self.specmn_fets)

        return fe_grid

    elstmr_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_elstmr_fe_level(self):
        return  FERefinementGrid(name = 'elastomer patch',
                                 fets_eval = self.elstmr_fets,
                                 domain = self.fe_domain)

    elstmr_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_elstmr_fe_grid(self):
        # coordinates adapt to discretization
        elstmr_min = self.length / 2.0 - self.elem_length_xy
        elstmr_max = self.length / 2.0
        z_max = self.thickness * 1.1
        return FEGrid(coord_min = (elstmr_min, elstmr_min, self.thickness),
                      coord_max = (elstmr_max, elstmr_max, z_max),
                      level = self.elstmr_fe_level,
                      shape = (1, 1, 1),
                      fets_eval = self.elstmr_fets)

    supprt_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_supprt_fe_level(self):
        return  FERefinementGrid(name = 'elastomer patch',
                                 fets_eval = self.supprt_fets,
                                 domain = self.fe_domain)

    supprt_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_supprt_fe_grid(self):
        # coordinates adapt to discretization
        supprt_min = 0
        supprt_max = self.elem_length_xy
        z_max = self.thickness * 0.1

        def tappered_support(pts):
            x_, y_, z_ = pts.T
            dx = np.max(x_) - np.min(y_)
            dy = np.max(y_) - np.min(y_)
            dz = np.max(z_) - np.min(z_)

            dz2 = dz / 2 / 2 / 2
            dz4 = dz / 6 / 2 / 2

            dz_red = (dz2 * x_ / dx + dz2 * y_ / dy - dz4 * x_ / dx * y_ / dy)

            pts = c_[x_, y_, z_ + dz_red]
            return pts

        return FEGrid(coord_min = (supprt_min, supprt_min, -z_max),
                      coord_max = (supprt_max, supprt_max, 0),
                      level = self.supprt_fe_level,
                      shape = (1, 1, 1),
                      geo_transform = tappered_support,
                      fets_eval = self.supprt_fets)

    tloop = Property(depends_on = 'input_change')
    @cached_property
    def _get_tloop(self):

        #fets_eval.mats_eval.nu = self.nu
        elstmr = self.elstmr_fe_grid
        specmn = self.specmn_fe_grid
        supprt = self.supprt_fe_grid

        self.fe_domain.n_dofs

        self.center_top_dofs = specmn[-1, -1, -1, -1, -1, -1].dofs
        self.edge_top_dofs = specmn[ 0, -1, -1, 0, -1, -1].dofs

        # if shape_xy is an even number:
        if self.shape_xy % 2 == 0:
            x_idx_center_edge = (self.shape_xy / 2) - 1
            self.center_edge_top_dofs = specmn[x_idx_center_edge, -1, -1, -1, -1, -1].dofs
        # if shape_xy is an odd number:
        else:
            # get the midside node of the center-edge-element
            x_idx_center_edge = (self.shape_xy - 1) / 2
            # valid only for quadratic elements

            self.center_edge_top_dofs = specmn[x_idx_center_edge, -1, -1, 1, -1, -1].dofs

        #--------------------------------------------------------------
        # boundary conditions for the symmetry and the single support
        #--------------------------------------------------------------
        bc_symplane_yz = BCSlice(var = 'u', value = 0., dims = [0], slice = specmn[-1, :, :, -1, :, :])
        bc_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1], slice = specmn[:, -1, :, :, -1, :])

        bc_supprt_symplane_yz = BCSlice(var = 'u', value = 0., dims = [0], slice = supprt[0, :, :, 0, :, :])
        bc_supprt_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1], slice = supprt[:, 0, :, :, 0, :])
        link_sp_sp_z = BCSlice(var = 'u', value = 0., dims = [2],
                             slice = specmn[0, 0, 0, :, :, 0],
                             link_slice = supprt[0, 0, -1, :, :, -1],
                             link_dims = [2],
                             link_coeffs = [1.])

        #--------------------------------------------------------------
        # boundary conditions for the symmetry and the single support
        #--------------------------------------------------------------
        support_000 = BCSlice(var = 'u', value = 0., dims = [2], slice = supprt[0, 0, 0, 0, 0, 0])

        #--------------------------------------------------------------
        # loading
        #--------------------------------------------------------------
        # w_max = center displacement:
        w_max = -0.028 # [m]

        #=======================================================================
        # elastomer conditions
        #=======================================================================
        bc_el_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1],
                                 slice = elstmr[:, -1, :, :, -1, :])
        bc_el_symplane_yz = BCSlice(var = 'u', value = 0., dims = [0],
                                 slice = elstmr[-1, :, :, -1, :, :])
#        link_el_sp = BCDofGroup(var = 'u', value = 0., dims = [2],
#                                get_dof_method = elstmr.get_bottom_dofs,
#                                get_link_dof_method = specmn.get_top_dofs,
#                                link_coeffs = [1.])
        link_el_sp = BCSlice(var = 'u', value = 0., dims = [2],
                             slice = elstmr[0, 0, 0, :, :, 0],
                             link_slice = specmn[-1, -1, -1, :, :, -1],
                             link_dims = [2],
                             link_coeffs = [1.])
        time_function = MFnLineArray(xdata = [0.0, 0.2, 0.4, 1.0], ydata = [0.0, 0.2, 0.75, 1.0])
        bc_el_w = BCSlice(var = 'u', value = w_max, dims = [2], time_function = time_function.get_value,
                                 slice = elstmr[:, :, -1, :, :, -1])
        #--------------------------------------------------------------
        # var 1: nodal displacement applied at top at the center of the plate
        #--------------------------------------------------------------
        bc_center_w = BCSlice(var = 'u', value = w_max, dims = [2], slice = specmn[-1, -1, :, -1, -1, :])

        #--------------------------------------------------------------
        # var 2: apply nodal displacements at top nodes of center element 
        #--------------------------------------------------------------
        bc_center_w_elem = BCSlice(var = 'u', value = w_max, dims = [2], slice = specmn[-1, -1, -1, :, :, -1])

        #--------------------------------------------------------------
        # var 3: displacement applied along the edge of the 2nd element row (from the middle) for 
        # shape_xy = 10 and an edge length of L_quarter = 0.625m the distance from the center equals 0.125m
        #--------------------------------------------------------------
        bc_center_w_xline = BCSlice(var = 'u', value = w_max, dims = [2], slice = specmn[-2:, -2, -1, :, 0, -1])
        bc_center_w_yline = BCSlice(var = 'u', value = w_max, dims = [2], slice = specmn[-2, -2:, -1, 0, :, -1])


        #--------------------------------------------------------------
        # ts 
        #--------------------------------------------------------------

        # force-displacement-diagram 
        # 
        center_dof = self.center_top_dofs[0, 0, 2]
        # center_top_line_dofs
        #
        ctl_dofs = elstmr[:, :, -1, :, :, -1].dofs[:, :, 2].flatten()

        # force-displacement-diagram 
        # 
        self.f_w_diagram_center = RTraceGraph(name = 'displacement (center) - reaction 2',
                                       var_x = 'U_k'  , idx_x = center_dof,
                                       # elastomer load
                                       var_y = 'F_int', idx_y_arr = ctl_dofs,
                                       record_on = 'update',
                                       transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the beam (2 symmetry axis):
                                       #
                                       transform_y = '-4 * 1000. * y')

        center_edge_dof = self.center_edge_top_dofs[0, 0, 2]
        self.f_w_diagram_center_edge = RTraceGraph(name = 'displacement (center/edge) - reaction 2',
                                                   var_x = 'U_k'  , idx_x = center_edge_dof,
                                                   var_y = 'F_int', idx_y_arr = ctl_dofs,
                                                   record_on = 'update',
                                                   transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                                   transform_y = '-4 * 1000 * y')

        edge_dof = self.edge_top_dofs[0, 0, 2]
        self.f_w_diagram_edge = RTraceGraph(name = 'displacement (edge) - reaction 2',
                                            var_x = 'U_k'  , idx_x = edge_dof,
                                            var_y = 'F_int', idx_y_arr = ctl_dofs,
                                            record_on = 'update',
                                            transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                            transform_y = '-4 * 1000 * y')

        ts = TS(
                sdomain = self.fe_domain,
                bcond_list = [bc_symplane_yz, bc_symplane_xz,
                              bc_el_symplane_yz, bc_el_symplane_xz, link_el_sp,
                              link_sp_sp_z, bc_supprt_symplane_yz, bc_supprt_symplane_xz,
                              support_000,
                               # var 1:
                               bc_el_w,
#                               bc_center_w,
#                               # var 2:
#                               bc_center_w_elem,
#                               # var 3:
#                               bc_center_w_xline, bc_center_w_yline
                              ],
                rtrace_list = [
                             self.f_w_diagram_center,
                             self.f_w_diagram_center_edge,
                             self.f_w_diagram_edge,
                             RTraceDomainListField(name = 'Displacement' ,
                                            var = 'u', idx = 0, warp = True),

#                             RTraceDomainListField(name = 'Stress' ,
#                                            var = 'sig_app', idx = 0, warp = True, 
#                                            record_on = 'update'),
#                             RTraceDomainListField(name = 'Strain' ,
#                                        var = 'eps_app', idx = 0, warp = True, 
#                                        record_on = 'update'),
                             RTraceDomainListField(name = 'Damage' ,
                                        var = 'omega_mtx', idx = 0, warp = True,
                                        record_on = 'update'),

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

        print 'tstep', self.tstep
        # Add the time-loop control
        tloop = TLoop(tstepper = ts,

#                       # allow only a low tolerance 
#                       #
#                       KMAX = 50, 
#                       tolerance = 5e-4, 

                       # allow a high tolerance 
                       #
                       KMAX = 100,
                       tolerance = 0.001,

#                       # allow a very high tolerance 
#                       #
#                       KMAX = 50,
#                       tolerance = 0.01,

                       RESETMAX = 0,
                       debug = False,
                       tline = TLine(min = 0.0, step = self.tstep, max = 1.0))

        return tloop

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name = 'u_center_top_z', unit = 'm'),
                 SimOut(name = 'F_max', unit = 'kN') ]


class SimSTDB(SimST):

    # vary the failure strain in PhiFnGeneralExtended:
    factor_eps_fail = Float(1.4, input = True,
                             ps_levels = (1.0, 1.2, 3))

    #-----------------
    # composite cross section unit cell:
    #-----------------
    #
    ccs_unit_cell_key = Enum('FIL-10-09_2D-05-11_0.00462_all0',
                              CCSUnitCell.db.keys(),
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

    calibration_test = Str(input = True)
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
    phi_fn = Property(Instance(PhiFnGeneralExtendedExp),
                       depends_on = 'input_change,+ps_levels')
    @cached_property
    def _get_phi_fn(self):
        return PhiFnGeneralExtendedExp(mfn = self.damage_function,
                                        Efp_frac = 0.2)
#        return PhiFnStrainHardening()
#        return PhiFnGeneralExtended( mfn = self.damage_function,
#                                     factor_e    ps_fail = self.factor_eps_fail )

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

#                                ccs_unit_cell_key = 'FIL-10-09_2D-02-06a_0.00273_90_0',
#                                calibration_test = 'TT11-10a-average',
#                                calibration_test = 'TT11-10a-V2',
#                                age = 28 )

    # s_tex,z = 4.62 mm corresponds to: n_tex = 12; h = 60 mm
    #
    sim_model = SimSTDB(ccs_unit_cell_key = 'FIL-10-09_2D-05-11_0.00462_all0',
                         #calibration_test = 'TT-12c-6cm-TU-SH1F-V1',
                         calibration_test = 'TT-12c-6cm-0-TU-SH2F-V3',
                         age = 27,
                         shape_xy = 10, #10,
                         shape_z = 4,
                         tstep = 0.1
                         )

    print 'sim_model.ccs_unit_cell_key', sim_model.ccs_unit_cell_key
    print 'sim_model.damage_function', sim_model.damage_function

    # plot the used 'phi_fn'
    #
    import pylab as p
    mfn = sim_model.damage_function
    eps_last = sim_model.damage_function.xdata[-1]
    print 'eps_last', eps_last
#    mfn.mpl_plot(p)
#    sim_model.phi_fn.mfn.mpl_plot(p)
    phi_fn_vect = np.vectorize( sim_model.phi_fn.get_value )
    xdata = np.linspace(0., eps_last*1.1, 400)
    ydata = phi_fn_vect( xdata )
#    p.plot(xdata, ydata)
    #p.show()

    do = 'ui'
#    do = 'pstudy'
#    do = 'validation'
#    do = 'show_last_results'

    if do == 'ui':
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource = sim_model)
        sim_model.tloop.eval()
        app.main()

    if do == 'pstudy':
        sim_ps = SimPStudy(sim_model = sim_model)
        sim_ps.configure_traits()

    if do == 'validation':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        # PT-12c-6cm-TU
        path = join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU')
        tests = [ 'ST-12c-6cm-u-TU.DAT' ]

        # PT-9a
#        path = join( simdb.exdata_dir, 'plate_tests', 'PT-9a' )
#        tests = [ 'PT01-9a.DAT', 'PT02-9a.DAT' , 'PT09-9a.DAT' ]

        # PT-10a
#        path = join( simdb.exdata_dir, 'plate_tests', 'PT-10a' )
#        tests = [ 'PT10-10a.DAT', 'PT11-10a.DAT' , 'PT12-10a.DAT' ]

        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
            ex_run.ex_type._plot_smoothed_force_deflection_center(p)

        sim_model.tloop.eval()

        sim_model.f_w_diagram_center.refresh()
        pickle_path = 'pickle_files'
        file_name = 'f_w_diagram_SH1F-V1_nelems10-10-4.pickle'
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

        file_name = 'f_w_diagram_SH1F-V1_nelems10-10-4.pickle'
#        file_name = 'f_w_diagram_SH1F-V1_step0-05_KMAX100_tol0-01_nelems8-8-3.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color = 'red')

        ### shape ###:

        # plate elem(12,12,4); w=0.08; step = 0.05; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-05_KMAX20_tol0-001_nelems12-12-4.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'blue' )

#        # plate elem(8,8,4); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_nelems8-8-4.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )

#        # plate elem(10,10,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_nelems10-10-2.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )


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
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
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


#        path = join( simdb.exdata_dir, 'plate_tests', 'PT-10a' )
#        tests = [ 'PT10-10a.DAT', 'PT11-10a.DAT' , 'PT12-10a.DAT' ]
#        tests = [ 'PT10-10a.DAT' ]

        path = join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU')
        tests = [ 'ST-12c-6cm-u-TU.DAT' ]

        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
            ex_run.ex_type._plot_smoothed_force_deflection_center(p)
#            ex_run.ex_type._plot_force_edge_deflection( p )
#            ex_run.ex_type._plot_force_center_edge_deflection( p )

        p.show()
