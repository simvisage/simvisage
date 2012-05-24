'''
Created on Aug 17, 2011

@author: rostar
'''

from etsproxy.traits.api import HasTraits, implements, Float, Property, \
    on_trait_change, cached_property
from etsproxy.traits.ui.api import View, Item, VGroup
from math import pi
from quaducom.ctt.homogenized_crack_bridges.i_homogenized_cb import ICB
import numpy as np
from quaducom.ctt import ICBM

def H(x):
    return np.sign(np.sign(x) + 1.)

class SteelBar(HasTraits):

    implements(ICB)
    implements(ICBM)

    # applied force
    P = Float(modified = True) # [N]

    # closest crack from left
    Ll = Float(modified = True) # [mm]

    # closest crack from right
    Lr = Float(modified = True) # [mm]

    dr = Float(6., auto_set = False, enter_set = True,
               desc = 'steel bar diameter in [mm]', modified = True, param = True)

    Er = Float(200e3, auto_set = False, enter_set = True,
               desc = 'steel modulus of elasticity [N/mm2]', modified = True, param = True)

    Em = Float(30e3, auto_set = False, enter_set = True,
               desc = 'matrix modulus of elasticity [N/mm2]', modified = True, param = True)

    tau = Float(10., auto_set = False, enter_set = True,
               desc = 'sheer force per unit area [N/mm2]', modified = True, param = True)

    Ac = Float(900., auto_set = False, enter_set = True,
             desc = 'composite cross section [mm2]', modified = True, param = True)

    Ar = Property(depends_on = 'dr')
    @cached_property
    def _get_Ar(self):
        return self.dr ** 2 / 4. * pi

    Am = Property(depends_on = 'Ac, dr')
    @cached_property
    def _get_Am(self):
        return self.Ac - self.Ar

    Kr = Property(depends_on = 'dr, Er')
    @cached_property
    def _get_Kr(self):
        return self.Ar * self.Er

    Km = Property(depends_on = 'Ac, dr, Em')
    @cached_property
    def _get_Km(self):
        return self.Am * self.Em

    Kc = Property(depends_on = 'Ac, dr, Er, Em')
    @cached_property
    def _get_Kc(self):
        return self.Kr + self.Km

    T = Property(depends_on = 'dr, tau')
    @cached_property
    def _get_T(self):
        return self.dr * pi * self.tau

    def get_eps_x_reinf(self, x):
        '''
        evaluation of strain profile in the vicinity of a crack bridge
        '''
        F = self.P - self.T * np.abs(x)
        eps = F / self.Er / self.Ar
        eps = eps * H(x + self.Ll) * H(self.Lr - x)

        return eps * H(eps)

    def get_force_x_reinf(self, x):
        '''
        evaluation of force profile in the vicinity of a crack bridge
        '''
        return self.get_eps_x_reinf(x) * self.Kr

    def get_sigma_x_reinf(self, x):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        return self.get_eps_x_reinf(x) * self.Er

    P_w = Property(depends_on = '+modified')
    @cached_property
    def _get_P_w(self):
        '''
        evaluation of force-crack width relation for a crack bridge
        '''

        P = self.P
        Lmin = np.minimum(self.Ll, self.Lr)
        Lmax = np.maximum(self.Ll, self.Lr)
        T = self.T

        # double sided debonding stage
        def w0(x):
            return x ** 2 / self.Er / self.Ar / T

        # force at which the double sided debonding is completed
        P0 = T * Lmin
        w00 = w0(P) * H(P0 - P)

        # one sided debonding stage
        def w1(x):
            return (((x - P0) + 2 * Lmin * T) ** 2 - (2 * Lmin * T) ** 2) / 2. / T / self.Ar / self.Er + w0(P0)

        # force at which the one sided debonding is completed
        P1 = T * Lmax
        w11 = w1(P) * H(P - P0) * H(P1 - P)
        # linear elastic stage
        def w2(x):
            return (P - P1) * (Lmin + Lmax) / self.Er / self.Ar + w1(P1)

        w22 = w2(P) * H(P - P1)

        return (w00 + w11 + w22) * H(P)

    traits_view = View(
                       VGroup(
                           Item('dr', resizable = False, springy = True),
                           Item('Er', resizable = False, springy = False),
                           Item('Em', resizable = False, springy = False),
                           Item('tau', resizable = False, springy = False),
                           Item('Ac', resizable = False, springy = False),
                           springy = True,
                           label = 'CB parameters',
                           dock = 'tab',
                           id = 'cb.steel_bar.params',
                        ),
                            id = 'cb.steel_bar',
                            dock = 'fixed',
                            scrollable = True,
                            resizable = True,
                            height = 0.8, width = 0.8
                                   )

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    sb = SteelBar(Ll = 60., Lr = 10., tau = 10.0, P = 1500)
    #sb.configure_traits()
    x = np.linspace(-70, 50, 200)
    eps = sb.get_eps_x_reinf(x)
    plt.plot(x, eps, lw = 2, color = 'black')
    plt.show()
