'''
Created on Jun 23, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    Float, Instance, Array, Property, cached_property, \
    HasStrictTraits, DelegatesTo, Int, Event, Callable, Button, \
    on_trait_change

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, ModelView, VGroup, HGroup, RangeEditor, InstanceEditor

import numpy as np
import pylab as p

from scipy.optimize import fsolve

from ecb_cross_section import \
    ECBCrossSection

from ecb_reinf_tex_uniform import \
    ECBReinfTexUniform

from ecb_matrix_cross_section import \
    ECBMatrixCrossSection

from ecb_law import ECBLBase

from util.traits.editors.mpl_figure_editor import  \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from matresdev.db.simdb import SimDB
simdb = SimDB()

class ECBLCalib(HasStrictTraits):

    # rupture moment and normal force measured in the calibration experiment
    # (three point bending test)
    #
    Mu = Float(3.5, calib_input=True) # [kNm]
    Nu = Float(0.0, calib_input=True) # [kN]

    #===========================================================================
    # Cross Section Specification (Geometry and Layout)
    #===========================================================================

    cs = Instance(ECBCrossSection)
    def _cs_default(self):
        reinf = [ECBReinfTexUniform(n_layers=3)]
        matrix = ECBMatrixCrossSection(width=0.1, n_cj=20)
        return ECBCrossSection(reinf=reinf,
                               matrix=matrix)

    ecb_law_type = DelegatesTo('cs')
    ecb_law = DelegatesTo('cs')
    cc_law_type = DelegatesTo('cs')
    cc_law = DelegatesTo('cs')
    width = DelegatesTo('cs')
    f_ck = DelegatesTo('cs')
    eps_c_u = DelegatesTo('cs')
    n_rovings = DelegatesTo('cs')
    n_layers = DelegatesTo('cs')

    notify_change = Callable(None)

    modified = Event
    @on_trait_change('cs.modified,+calib_input')
    def _set_modified(self):
        self.modified = True
        if self.notify_change != None:
            self.notify_change()

    u0 = Property(Array(float), depends_on='cs.modified')
    '''Construct the initial vector.
    '''
    @cached_property
    def _get_u0(self):
        u0 = self.ecb_law.u0
        #eps_up = u0[1]
        eps_up = -self.cs.eps_c_u
        eps_lo = self.cs.convert_eps_tex_u_2_lo(u0[0])

        print 'eps_up', eps_up
        print 'eps_lo', eps_lo

        return np.array([eps_lo, u0[1] ], dtype='float')

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
        eps_up = -self.cs.eps_c_u
        eps_lo = u[0]

        self.cs.set(eps_lo=eps_lo, eps_up=eps_up)

        eps_tex_u = self.cs.convert_eps_lo_2_tex_u(u[0])
        self.cs.ecb_law.set_cparams(eps_tex_u, u[1])

        N_internal = self.cs.N
        M_internal = self.cs.M

        d_N = N_internal - self.Nu
        d_M = M_internal - self.Mu

        return np.array([ d_N, d_M ], dtype=float)

    u_sol = Property(Array(Float), depends_on='cs.modified,+calib_input')
    '''Solution vector returned by 'fit_response'.'''
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
        return fsolve(self.get_lack_of_fit, self.u0, xtol=1.0e-5)

    calibrated_ecb_law = Property(depends_on='cs.modified,+calib_input')
    '''Calibrated ecbl_mfn
    '''
    @cached_property
    def _get_calibrated_ecb_law(self):
        print 'NEW CALIBRATION'
        self.ecb_law.set_cparams(*self.u_sol)
        return self.ecb_law

    view = View(Item('Mu'),
                Item('Nu'),
                buttons=['OK', 'Cancel']
                )

class ECBLCalibModelView(ModelView):
    '''Model in a viewable window.
    '''
    model = Instance(ECBLCalib)
    def _model_default(self):
        return ECBLCalib()

    cs_state = Property(Instance(ECBCrossSection), depends_on='model')
    @cached_property
    def _get_cs_state(self):
        return self.model.cs

    data_changed = Event

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor='white')
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

    calibrated_ecb_law = Property(Instance(ECBLBase), depends_on='model')
    @cached_property
    def _get_calibrated_ecb_law(self):
        return self.model.calibrated_ecb_law

    view = View(HSplit(VGroup(
                       Item('cs_state', label='Cross section', show_label=False),
                       Item('model@', show_label=False),
                       Item('calibrated_ecb_law@', show_label=False, resizable=True),
                       ),
                       Group(HGroup(
                             Item('replot', show_label=False),
                             Item('clear', show_label=False),
                      ),
                      Item('figure', editor=MPLFigureEditor(),
                           resizable=True, show_label=False),
                      id='simexdb.plot_sheet',
                      label='plot sheet',
                      dock='tab',
                      ),
                       ),
                width=0.5,
                height=0.4,
                buttons=['OK', 'Cancel'],
                resizable=True)

if __name__ == '__main__':

    #------------------------------------------------
    # 1) CALIBRATION:
    # get 'eps_t' and the parameter of the effective 
    # crack bridge function 'var_a' for a given 'eps_c_u'
    #------------------------------------------------
    #
    print '\n'
    print 'setup ECBLCalib'
    print '\n'
    p.plot([0, 0], [0, 2.4e3])

    ec = ECBLCalib(# mean concrete strength after 9 days
                           # 7d: f_ck,cube = 62 MPa; f_ck,cyl = 62/1.2=52
                           # 9d: f_ck,cube = 66.8 MPa; f_ck,cyl = 55,7
                           f_ck=55.7,

                           # measured strain at bending test rupture (0-dir)
                           #
                           eps_c_u=3.3 / 1000.,

                           # measured value in bending test [kNm]
                           # value per m: M = 5*3.49
                           #
                           Mu=3.49,
                       )

    ecw = ECBLCalibModelView(model=ec)
    ecw.configure_traits()

    for sig_tex_u, color in zip([1200, 1300, 1400], ['red', 'green', 'blue', 'black', 'orange', 'brown']):
    #for sig_tex_u, color in zip([1216], ['red']):

        #for ecbl_type in ['linear', 'cubic', 'fbm']:
        for ecbl_type in ['cubic']:
            print 'CALIB TYPE', ecbl_type
            ec.n = 0
            ec.cs.ecb_law_type = ecbl_type
            ec.cs.ecb_law.sig_tex_u = sig_tex_u
            ec.get_lack_of_fit(ec.u0)

            ec.calibrated_ecb_law.mfn.plot(p, color=color, linewidth=8)
            print 'E_yarn', ec.calibrated_ecb_law.mfn.get_diff(0.00001)
            print 'INTEG', ec.calibrated_ecb_law.mfn.integ_value

#            ec.ecbl_type = 'bilinear'
#            ec.ecbl.sig_tex_u = sig_tex_u
#            for eps_el_fraction in np.linspace(0.25, 0.99999, 4):
#                ec.n = 0
#                ec.ecbl.eps_el_fraction = eps_el_fraction
#                ec.ecbl_mfn.plot(p, color = color)
#                print 'E_yarn', ec.ecbl_mfn.get_diff(0.00001)
#                print 'INTEG', ec.ecbl_mfn.integ_value
    p.plot([0.0, 0.01], [0.0, 2400], color='black')

    p.show()
