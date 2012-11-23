'''
Created on Jun 23, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    Float, Instance, Array, Property, cached_property, \
    HasStrictTraits, DelegatesTo, Int, Event, Callable, \
    Button, on_trait_change

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, ModelView, VGroup, HGroup, RangeEditor, InstanceEditor

import numpy as np
import pylab as p

from util.traits.editors.mpl_figure_editor import  \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from scipy.optimize import fsolve

from constitutive_law import \
    ConstitutiveLawModelView

from ecb_cross_section_state import \
    ECBCrossSectionState

from ecb_cross_section_geo import \
    ECBCrossSectionGeo

from ecb_law import ECBLBase

from matresdev.db.simdb import SimDB
simdb = SimDB()

class ECBLCalibState(HasStrictTraits):

    # rupture moment and normal force measured in the calibration experiment
    # (three point bending test)
    #
    Mu = Float(3.5, enter_set = True, auto_set = False, input = True) # [kNm]
    Nu = Float(0.0, enter_set = True, auto_set = False, input = True) # [kN]

    #===========================================================================
    # Cross Section Specification (Geometry and Layout)
    #===========================================================================

    cs_geo = Instance(ECBCrossSectionGeo)
    def _cs_geo_default(self):
        return ECBCrossSectionGeo(notify_change = self.set_modified)

    cs_state = Property(Instance(ECBCrossSectionState), depends_on = 'cs_geo')
    @cached_property
    def _get_cs_state(self):
        return ECBCrossSectionState(cs_geo = self.cs_geo,
                                    notify_change = self.set_modified)

    notify_change = Callable(None)

    modified = Event
    @on_trait_change('Mu, Nu')
    def set_modified(self):
        self.modified = True
        if self.notify_change != None:
            self.notify_change()

    u0 = Property(Array(float), depends_on = 'cs_geo.modified, cs_state.modified')
    @cached_property
    def _get_u0(self):
        u0 = self.cs_state.ecb_law.u0
        eps_up = self.cs_geo.cc_law.eps_c_u
        self.cs_state.set(eps_up = -eps_up)
        eps_lo = self.cs_state.convert_eps_tex_u_2_lo(u0[0])
        print 'eps_up', self.cs_state.eps_up
        print 'eps_lo', eps_lo
        return np.array([eps_lo, u0[1] ], dtype = 'float')

    # iteration counter
    #
    n = Int(0)
    def get_lack_of_fit(self, u):
        '''Return the difference between 'N_external' and 'N_internal' as well as 'M_external' and 'M_internal'
        N_c (=compressive force of the compressive zone of the concrete)
        N_t (=total tensile force of the reinforcement layers)
        '''

        print '--------------------iteration', self.n, '------------------------'
        self.n += 1
        # set iteration counter
        #
        eps_up = -self.cs_geo.cc_law.eps_c_u
        eps_lo = u[0]
        self.cs_state.set(eps_lo = u[0], eps_up = -self.cs_geo.cc_law.eps_c_u)

        eps_tex_u = self.cs_state.convert_eps_lo_2_tex_u(u[0])

        print 'eps_up', eps_up
        print 'eps_lo', eps_lo

        print 'u', eps_tex_u, u[1] / eps_tex_u
        self.cs_geo.ecb_law.set_cparams(eps_tex_u, u[1])

        N_internal = self.cs_state.N
        M_internal = self.cs_state.M

        d_N = N_internal - self.Nu
        d_M = M_internal - self.Mu

        print 'R', d_M, d_N
        return np.array([ d_M, d_N ], dtype = float)

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
        self.cs_geo.ecb_law.set_cparams(*self.u_sol)
        return self.cs_geo.ecb_law

    view = View(Item('Mu'),
                Item('Nu'),
                buttons = ['OK', 'Cancel']
                )

class ECBLCalibStateModelView(ModelView):
    '''Model in a viewable window.
    '''
    model = Instance(ECBLCalibState)
    def _model_default(self):
        return ECBLCalibState()

    cs_state = Property(Instance(ECBCrossSectionState), depends_on = 'model')
    @cached_property
    def _get_cs_state(self):
        return self.model.cs_state

    data_changed = Event

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor = 'white')
        return figure

    replot = Button()
    def _replot_fired(self):
        ax = self.figure.add_subplot(1, 1, 1)
        self.model.calibrated_ecb_law.plot(ax)
        self.data_changed = True

    clear = Button()
    def _clear_fired(self):
        self.figure.clear()
        self.data_changed = True

    calibrated_ecb_law = Property(Instance(ECBLBase), depends_on = 'model')
    @cached_property
    def _get_calibrated_ecb_law(self):
        return self.model.calibrated_ecb_law

    view = View(HSplit(VGroup(
                       Item('cs', label = 'Cross section', show_label = False),
                       Item('model@', show_label = False),
                       Item('calibrated_ecb_law@', show_label = False, resizable = True),
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
                width = 0.5,
                height = 0.4,
                buttons = ['OK', 'Cancel'],
                resizable = True)

if __name__ == '__main__':

    #------------------------------------------------
    # 1) CALIBRATION:
    # get 'eps_t' and the parameter of the effective 
    # crack bridge function 'var_a' for a given 'eps_c_u'
    #------------------------------------------------
    #
    print '\n'
    print 'setup ECBLCalibState'
    print '\n'
    p.plot([0, 0], [0, 2.4e3])

    cs_geo = ECBCrossSectionGeo(
                   # mean concrete strength after 9 days
                   # 7d: f_ck,cube = 62 MPa; f_ck,cyl = 62/1.2=52
                   # 9d: f_ck,cube = 66.8 MPa; f_ck,cyl = 55,7
                   cc_law_params = dict(bilinear = dict(eps_c_u = 0.0033,
                                                        f_ck = 55.7,
                                                        E_c = 29e+3),
                                        quadratic = dict(eps_c_u = 0.0033,
                                                        f_ck = 55.7,
                                                        ),
                                        constant = dict(eps_c_u = 0.0033,
                                                        f_ck = 55.7,
                                                        )
                                        ),

                   thickness = 0.06,

                   width = 0.2,

                   n_layers = 12,

                   n_rovings = 23,
                   # measured strain at bending test rupture (0-dir)
                   #

                   cc_law_type = 'constant',
                   ecb_law_type = 'fbm'
                   )

    ec = ECBLCalibState(cs_geo = cs_geo,
                   # measured value in bending test [kNm]
                   # value per m: M = 5*3.49
                   #
                   Mu = 3.49,
                   )

    u_sol = ec.u_sol
    print 'u_sol', u_sol
    eps_tex = ec.cs_state.convert_eps_lo_2_tex_u(u_sol[0])
    print 'sig_max', u_sol[1] * eps_tex
    ecw = ECBLCalibStateModelView(model = ec)
    ecw.configure_traits()
