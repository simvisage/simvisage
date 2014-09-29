'''
Created on 3.12.2010

Response function for a double sided pullout of a continuous filament
with the same friction at both sides.

@author: Q
'''

from etsproxy.traits.api import \
    Bool, Float, Str, implements, Int, Instance, Event, on_trait_change

from etsproxy.traits.ui.api import \
    View, Item, HSplit, VGroup

from math import pi, e

from numpy import \
    frompyfunc, sign, linspace, array, cos, sqrt, argmax, hstack, max
import numpy as np

from matplotlib import pyplot as plt

from spirrid.i_rf import IRF

from spirrid.rf import RF

def H(x):
    return (sign(x) + 1.0) / 2.0

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import Figure


class DoublePulloutSym(RF):

    implements(IRF)

    title = Str('symetrical yarn pullout')

    xi = Float(0.0179, auto_set=False, enter_set=True, input=True,
                distr=['weibull_min', 'uniform'])

    tau_fr = Float(2.5, auto_set=False, enter_set=True, input=True,
                distr=['uniform', 'norm'])
    # free length
    l = Float(0.0, auto_set=False, enter_set=True, input=True,
              distr=['uniform'])

    d = Float(26e-3, auto_set=False, input=True,
              enter_set=True, distr=['uniform', 'weibull_min'])

    E_mod = Float(72.0e3, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])
    # slack
    theta = Float(0.01, auto_set=False, enter_set=True, input=True,
                  distr=['uniform', 'norm'])

    phi = Float(1., auto_set=False, enter_set=True, input=True,
                  distr=['uniform', 'norm'])

    # embedded length
    L = Float(1., auto_set=False, enter_set=True, input=True,
              distr=['uniform'])

    free_fiber_end = Bool(True, input=True)

    w = Float(enter_set=True, input=True, ctrl_range=(0, 1, 10))

    weave_code = '''
        '''

    def __call__(self, w, tau_fr, l, d, E_mod, theta, xi, phi, L):
        '''Return the force for a prescribed crack opening displacement w.
        '''
        A = pi * d ** 2 / 4.
        l = l * (1 + theta)
        w = w - theta * l
        Tau = tau_fr * phi * d * pi
        P_ = 0.5 * (-l * Tau + sqrt(l ** 2 * Tau ** 2 + 4 * w * H(w) * E_mod * A * Tau))
        # one sided pullout P_ = ( -l * Tau + sqrt( l ** 2 * Tau ** 2 + 2 * w * H( w ) * E_mod * A * Tau ) )

        if self.free_fiber_end:
            # ------ FREE LENGTH -------

            # frictional force along the bond length
            P_fr = Tau * (L - l)

            # if pullout_criterion positive - embedded
            # otherwise pulled out
            #
            pull_out_criterion = P_fr - P_
            P_ = P_ * H(pull_out_criterion) + P_fr * H(-pull_out_criterion)
        else:
            # --------------------------
            # ------ clamped fiber end ---------
            v = L * (l * Tau + Tau * L) / E_mod / A
            P_ = P_ * H(Tau * L - P_) + (Tau * L + (w - v) / (l + 2 * L) * A * E_mod) * H(P_ - Tau * L)
            # ----------------------------------
        P = P_ * H(A * E_mod * xi - P_)
        return P

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor='white')
        return figure

    changed = Event
    @on_trait_change('+input')
    def _set_changed(self):
        self.changed = True

    data_changed = Event
    @on_trait_change('+input')
    def refresh(self):
        figure = self.figure
        figure.clear()
        axes = figure.gca()

        P_fn = lambda w: self.__call__(w, self.tau_fr, self.l, self.d,
                                        self.E_mod, self.theta, self.xi, self.phi, self.L)
        pyF = frompyfunc(P_fn, 1, 1)

        w_arr = linspace(0.0, 1.0, 100)
        P_arr = array(pyF(w_arr), dtype='float_')

        axes.plot(w_arr, P_arr, lw=1.0, color='blue')

        self.data_changed = True

    group_attribs = VGroup(Item('tau_fr'),
                                     Item('l'),
                                     Item('d'),
                                     Item('E_mod'),
                                     Item('theta'),
                                     Item('xi'),
                                     Item('phi'),
                                     Item('L'),
                                     Item('free_fiber_end'),
                                     ),

    traits_view = View(group_attribs,
                        scrollable=True,
                        resizable=True,
                        id='mk.figure.attribs',
                        dock='tab',
                        )

    traits_view_diag = View(HSplit(group_attribs,
                              VGroup(Item('figure',
                                            editor=MPLFigureEditor(),
                                            show_label=False,
                                            resizable=True),
                                id='mk.figure.view'
                                       ),),
                                id='mk.view',
                                buttons=['OK', 'Cancel'],
                                resizable=True,
                                width=600, height=400
                        )


if __name__ == '__main__':
    ds = DoublePulloutSym(tau_fr=0.1, l=10, D=0.026, E=72000, theta=0.0,
                           xi=0.019, phi=1.0, L=30)
    ds.configure_traits(view='traits_view_diag')
#    X = linspace( 0, .05, 100 )
#    Y = ds( X, 2.5, 0.01, 26e-3, 70.0e3, 0.01, 0.014, 1., 2. )
#    plt.plot( X, Y, linewidth = 2 )
#    plt.show()
