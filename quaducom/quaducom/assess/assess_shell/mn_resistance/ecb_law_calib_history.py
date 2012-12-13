'''
Created on Jun 23, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    Float, Instance, Array, Property, cached_property, \
    HasStrictTraits, DelegatesTo, Int, Event, Callable, \
    Button, on_trait_change, List, Constant

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, ModelView, VGroup, HGroup, RangeEditor, InstanceEditor, TabularEditor

from etsproxy.traits.ui.tabular_adapter \
    import TabularAdapter

import numpy as np
import pylab as p

from util.traits.editors.mpl_figure_editor import  \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from scipy.optimize import fsolve
from scipy.optimize import broyden1

from constitutive_law import \
    ConstitutiveLawModelView

from ecb_cross_section_state import \
    ECBCrossSectionState

from ecb_cross_section_geo import \
    ECBCrossSectionGeo

from ecb_law import ECBLBase

from matresdev.db.simdb import SimDB
simdb = SimDB()

class ECBLCalibHist(HasStrictTraits):

    # History of measured bending moments and compressive strains
    # 
    eps_c_arr = Array(float)
    M_arr = Array(float)

    @on_trait_change('eps_c_arr, M_arr')
    def _reset_ecb_law(self):
        self.cs_geo.ecb_law.n_eps = len(self.M_arr) + 1

    n_states = Property(depends_on = 'eps_c_arr, M_arr')
    @cached_property
    def _get_n_states(self):
        # assert that both arrays have the same length
        return self.eps_c_arr.shape[0]

    #===========================================================================
    # Cross Section Specification (Geometry and Layout)
    #===========================================================================

    cs_geo = Instance(ECBCrossSectionGeo)
    def _cs_geo_default(self):
        cs_geo = ECBCrossSectionGeo(notify_change = self.set_modified,
                                   ecb_law_type = 'piecewise_linear')
        cs_geo.ecb_law.n_eps = len(self.M_arr) + 1

    def _cs_geo_changed(self):
        self.cs_geo.ecb_law_type = 'piecewise_linear'
        self.cs_geo.ecb_law.n_eps = len(self.M_arr) + 1

    cs_states = Property(List(Instance(ECBCrossSectionState)), depends_on = 'cs_geo')
    @cached_property
    def _get_cs_states(self):
        return [ECBCrossSectionState(cs_geo = self.cs_geo,
                                     notify_change = self.set_modified)
                for k in range(self.n_states)]

    notify_change = Callable(None)

    modified = Event
    @on_trait_change('eps_c_arr, M_arr')
    def set_modified(self):
        self.modified = True
        if self.notify_change != None:
            self.notify_change()

    u0 = Property(Array(float), depends_on = 'cs_geo.modified')
    @cached_property
    def _get_u0(self):
        ecb_law = self.cs_geo.ecb_law
        eps_tex_u = ecb_law.eps_tex_u
        u0 = np.zeros((2 * self.n_states), dtype = 'float')
        d_eps = eps_tex_u / self.n_states
        u0[:self.n_states] = np.repeat(d_eps, n_states)
        u0[self.n_states:] = ecb_law.sig_level_arr[1:]
        #u0[:] = [1.233333e-02, 8.0e+4]
        return u0 * 0.1

    # iteration counter
    #
    n = Int(0)
    def get_lack_of_fit(self, u):
        '''Return the difference between 'N_external' and 'N_internal' as well as 'M_external' and 'M_internal'
        N_c (=compressive force of the compressive zone of the concrete)
        N_t (=total tensile force of the reinforcement layers)
        '''

        print '--------------------iteration', self.n, '------------------------'

        print 'u', u

        sig_arr = np.hstack([[0.0], u[self.n_states:]])
        d_eps_arr = np.hstack([[0.0], u[:self.n_states]])

        eps_arr = np.cumsum(d_eps_arr)

        self.n += 1
        # set iteration counter
        #
        for k, (state, eps_c) in enumerate(zip(self.cs_states, self.eps_c_arr)):
            state.eps_up = eps_c
            eps_lo = state.convert_eps_tex_u_2_lo(eps_arr[k + 1])
            state.eps_lo = eps_lo

        #self.cs_geo.ecb_law.set_cparams(eps_tex_u, c_params)
        self.cs_geo.ecb_law.set_sig_eps_arr(eps_arr, sig_arr)

        d_MN_list = []
        for k, (state, M) in enumerate(zip(self.cs_states, self.M_arr)):

            N_internal = state.N
            M_internal = state.M

            d_N = N_internal - 0.0
            d_M = M_internal - M

            d_MN_list += [d_M, d_N]

        dR = np.array(d_MN_list, dtype = float)

        print 'R', dR
        return dR

    # solution vector returned by 'fit_response'
    #
    u_sol = Property(Array(Float), depends_on = 'modified')
    @cached_property
    def _get_u_sol(self):
        '''iterate 'eps_t' such that the lack of fit between the calculated
        normal forces in the tensile reinforcement and the compressive zone (concrete)
        is smaller then 'xtol' defined in function 'brentq'.
        NOTE: the method 'get_lack_of_fit' returns the relative error.
        '''

        # use scipy-functionality to get the iterated value of 'eps_t'
        # NOTE: get_lack_of_fit must have a sign change as a requirement
        # for the function call 'brentq' to work property. 

        # The method brentq has optional arguments such as
        #   'xtol'    - absolut error (default value = 1.0e-12)
        #   'rtol'    - relative error (not supported at the time)
        #   'maxiter' - maximum numbers of iterations used
        #
        return fsolve(self.get_lack_of_fit, self.u0, xtol = 1.0e-5)

    #===========================================================================
    # Calibrated ecb_law_mfn
    #===========================================================================

    calibrated_ecb_law = Property(depends_on = 'modified')
    @cached_property
    def _get_calibrated_ecb_law(self):
        print 'NEW CALIBRATION'
        u = self.u_sol
        eps_tex_u = self.cs_states[-1].convert_eps_lo_2_tex_u(u[self.n_states - 1])
        c_params = np.hstack([[0.0], u[self.n_states:]])
        self.cs_geo.ecb_law.set_cparams(eps_tex_u, c_params)
        return self.cs_geo.ecb_law

    view = View(buttons = ['OK', 'Cancel']
                )

#----------------------------------------------------------------------------------
# Tabular Adapter Definition 
#----------------------------------------------------------------------------------
class CSStatesTabularAdapter (TabularAdapter):

    columns = [ ('eps_up', 'eps_up'),
                ('eps_lo', 'eps_lo'),
                ]

    font = 'Courier 10'
    variable_alignment = Constant('right')

#----------------------------------------------------------------------------------
# Tabular Editor Construction 
#----------------------------------------------------------------------------------
cs_states_editor = TabularEditor(
    selected = 'current_cs_state',
    adapter = CSStatesTabularAdapter(),
    operations = [ 'move' ],
    auto_update = True
 )

class ECBLCalibHistModelView(ModelView):
    '''Model in a viewable window.
    '''
    model = Instance(ECBLCalibHist)
    def _model_default(self):
        return ECBLCalibHist()

    cs_states = Property(List(Instance(ECBCrossSectionState)), depends_on = 'model')
    @cached_property
    def _get_cs_states(self):
        return self.model.cs_states

    current_cs_state = Instance(ECBCrossSectionState)
    def _current_cs_default(self):
        self.model.cs_states[0]

    def _current_cs_state_changed(self):
        self._replot_fired()

    eps_range = Property
    def _get_eps_range(self):
        eps_arr = [[cs_state.eps_up, cs_state.eps_lo]
                     for cs_state in self.cs_states ]
        eps_range = np.asarray(eps_arr)
        print 'eps_range', eps_range
        return (-np.max(eps_range[:, 1]), -np.min(eps_range[:, 0]))

    f_range = Property
    def _get_f_range(self):
        f_arr = [[np.max(cs_state.f_ti_arr), np.min(cs_state.f_cj_arr)]
                     for cs_state in self.cs_states ]
        f_range = np.asarray(f_arr)
        print 'eps_range', f_range
        return (-np.max(f_range[:, 0]), -np.min(f_range[:, 1]))

    data_changed = Event

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor = 'white')
        return figure

    replot = Button()
    def _replot_fired(self):
        if self.current_cs_state:
            cs = self.current_cs_state
        else:
            cs = self.cs_states[0]
        cs.plot(self.figure,
                eps_range = self.eps_range,
                f_range = self.f_range)
        self.data_changed = True

    clear = Button()
    def _clear_fired(self):
        self.figure.clear()
        self.data_changed = True

    view = View(HSplit(VGroup(
                       Item('cs_states', editor = cs_states_editor, label = 'Cross section', show_label = False),
                       Item('model@', show_label = False),
                       ),
                       Group(HGroup(
                             Item('replot', show_label = False),
                             Item('clear', show_label = False),
                      ),
                      Item('figure', editor = MPLFigureEditor(),
                           resizable = True, show_label = False),
                      id = 'simexdb.plot_sheet',
                      label = 'plot sheet',
                      dock = 'tab',
                      ),
                       ),
                width = 0.8,
                height = 0.7,
                buttons = ['OK', 'Cancel'],
                resizable = True)


if __name__ == '__main__':

    #------------------------------------------------
    # 1) CALIBRATION:
    # get 'eps_t' and the parameter of the effective 
    # crack bridge function 'var_a' for a given 'eps_c_u'
    #------------------------------------------------
    #

    import os.path

    from matresdev.db.simdb import \
        SimDB

    simdb = SimDB()

    from mathkit.mfn import MFnLineArray
    from matresdev.db.exdb.ex_run_view import ExRunView
    from matresdev.db.exdb.ex_run import ExRun

    ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2012-10-11_BT-3PT-6c-2cm-0-TU_bs',
                            'BT-6c-2cm-0-TU-V1.csv')

    exrun = ExRun(ex_path)

    print 'F_max', np.max(exrun.ex_type.F)

    eps_c = exrun.ex_type.eps_c_smoothed
    M = exrun.ex_type.M_smoothed

    eps_c_M_fn = MFnLineArray(xdata = M, ydata = eps_c)
    eps_c_M_vct_fn = np.vectorize(eps_c_M_fn.get_value)
    M_max = M[-1]
    eps_c_max = np.min(eps_c)
    n_states = 4

    #M_arr = np.linspace(0, M_max, n_states + 1)[1:]
    cf = np.linspace(0.8, 1.0, n_states)
    print 'cf', cf
    M_arr = cf * M_max
    eps_c_arr = eps_c_M_vct_fn(M_arr)

    print 'eps_c_max', eps_c_max
    print 'M_arr', M_arr
    print 'eps_c_arr', eps_c_arr

    cs_geo = ECBCrossSectionGeo(
                   # mean concrete strength after 9 days
                   # 7d: f_ck,cube = 62 MPa; f_ck,cyl = 62/1.2=52
                   # 9d: f_ck,cube = 66.8 MPa; f_ck,cyl = 55,7
                   thickness = 0.02,

                   width = 0.10,

                   n_layers = 6,

                   n_rovings = 11,

                   cc_law_params = dict(bilinear = dict(eps_c_u = -eps_c_max,
                                                        f_ck = 75.7,
                                                        E_c = 29e+3)
                                        ),

                   cc_law_type = 'bilinear',
                   ecb_law_type = 'piecewise_linear'
                   )

    ec = ECBLCalibHist(cs_geo = cs_geo,
                       M_arr = M_arr,
                       eps_c_arr = eps_c_arr,
                   )


    print ec.cs_geo.ecb_law.n_eps
    print ec.n_states
    u0 = ec.u0
    print 'u0', u0

    u_sol = ec.u_sol
    print 'u_sol', u_sol
    print 'E_eff', u_sol[ ec.n_states] / u_sol[0]

    ecv = ECBLCalibHistModelView(model = ec)

    ecv.configure_traits(kind = 'live')
