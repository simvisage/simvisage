#-------------------------------------------------------------------------
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
# Created on Jun 2, 2010 by: rch


from math import pi, e
import types

from matplotlib import pyplot as plt
from numpy import sign, linspace, cos, sqrt, exp
from traits.api import \
    Float, Str, implements, Range  # @UnresolvedImport

from stats.spirrid import IRF, RF


def Heaviside(x):
    return (sign(x) + 1.0) / 2.0


class ConstantFrictionFiniteFiber(RF):

    implements(IRF)

    title = Str('pull-out with constant friction')

    fu = Float(1200e6, auto_set=False, enter_set=True,
               distr=['weibull_min'])

    qf = Float(1500, auto_set=False, enter_set=True,
               distr=['uniform', 'norm'])

    L = Float(0.02, auto_set=False, enter_set=True,
              distr=['uniform'])

    A = Float(5.30929158457e-10, auto_set=False, enter_set=True,
              distr=['uniform', 'weibull_min'])

    E_mod = Float(70.0e9, auto_set=False, enter_set=True,
                  distr=['uniform'])

    phi = Range(0, pi, auto_set=False, enter_set=True,
                distr=['sin_distr', 'cos_distr'])

    z = Float(0, auto_set=False, enter_set=True,
              distr=['uniform'])

    f = Float(0.01, auto_set=False, enter_set=True,
              distr=['uniform'])

    w = Float(0.0, auto_set=False, enter_set=True,
              ctrl_range=(0, 0.01, 10))

    c_code = '''
            double w = eps;
            double Le = L / 2. - z / cos( phi );
            double w_deb = exp( f * phi ) * qf * pow(Le,2.0) / E_mod / A;
            double P_deb_full = sqrt( 2. * w / 2. * E_mod * A * qf ) * exp( f * phi );
            double P_deb;
            
            // Heaviside
            if ( Le < 0 || P_deb_full > fu * A || w > w_deb ){
                P_deb = 0;
            }else{
                P_deb =P_deb_full;
            }
            
            double P_pull_x = ( Le * qf - Le * qf / ( Le - w_deb ) * ( w - w_deb ) ) * exp( f * phi );
            double P_pull;
            
            // Heaviside 
            if ( P_pull_x < 0 || w_deb > w ){
                P_pull = 0;
            }else{
                P_pull = P_pull_x;
            }
            
            // Computation of the q( ... ) function
            q = P_deb + P_pull;
        '''

    def __callx__(self, w, fu, qf, L, A, E_mod, z, phi, f):
        '''Intial vectorized implementation - without regarding
        the lexical structure of the expression.
        '''
        Le = L / 2. - z / cos(phi)
        w_deb = e ** (f * phi) * qf * Le ** 2.0 / E_mod / A
        P_deb_full = sqrt(2. * w / 2. * E_mod * A * qf) * e ** (f * phi)
        P_deb = P_deb_full * \
            Heaviside(fu * A - P_deb_full) * \
            Heaviside(w_deb - w) * Heaviside(Le)
        P_pull_x = (
            Le * qf - Le * qf / (Le - w_deb) * (w - w_deb)) * e ** (f * phi)
        P_pull = P_pull_x * Heaviside(P_pull_x) * Heaviside(w - w_deb)
        return P_deb + P_pull

    def __call__(self, w, fu, qf, L, A, E_mod, z, phi, f):
        '''Lexically optimized expresseion - each result is 
        calculated only once. 
        Further optimization possible by printing out the shape
        and doing inplace operation where possible.
        However, this does not seem to have a significant effect.
        '''

        t4 = sqrt(w * E_mod * A * qf)
        t5 = f * phi
        t6 = exp(t5)
        t7 = t4 * t6
        t11 = Heaviside(fu * A - 0.1000000000e1 * t7)
        t12 = exp(t5)
        t14 = 0.5000000000e0 * L
        t15 = cos(phi)
        t17 = z / t15
        t18 = t14 - t17
        t19 = pow(t18, 0.20e1)
        t24 = t12 * qf * t19 / E_mod / A
        t25 = t24 - w
        t26 = Heaviside(t25)
        t28 = Heaviside(t18)
        t32 = t18 * qf
        t38 = (t32 + t32 / (t14 - t17 - t24) * t25) * t6
        t39 = Heaviside(t38)
        t40 = Heaviside(-t25)
        res = 0.1000000000e1 * t7 * t11 * t26 * t28 + t38 * t39 * t40
        return res


class ConstantFrictionAndFreeLength(RF):
    '''
    '''

    implements(IRF)

    title = Str('pull-out with constant friction and free length ')

    tau = Float(8, auto_set=False, enter_set=True,
                distr=['uniform'])

    # free length
    l = Float(1, auto_set=False, enter_set=True,
              distr=['uniform', 'norm'])

    E = Float(70e9, auto_set=False, enter_set=True,
              distr=['uniform'])

    A = Float(5.30929158457e-10, auto_set=False, enter_set=True,
              distr=['uniform', 'weibull_min'])

    # waviness in strains
    slack = Float(0.1, auto_set=False, enter_set=True,
                  distr=['uniform'])

    u = Float(auto_set=False, enter_set=True,
              ctrl_range=(0.0, 1.0, 10))

    def __call__(self, u, tau, l, E, A, slack):
        return -l * (1 + slack) * tau * Heaviside(u - l * (slack)) + \
            + sqrt((l * (1 + slack) * tau) ** 2
                   + 2 * E * A * (u - (l * slack)) * Heaviside(u - l * (slack)))


if __name__ == '__main__':
    cf = ConstantFrictionAndFreeLength(slack=0)
    X = linspace(0, 0.0012, 100)
    Y = cf(X, 8.0, 0.1, 210.e9, 0.004, 0.0)
    plt.plot(X, Y, linewidth=2)
    plt.show()
