#-------------------------------------------------------------------------------
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
# Created on Jun 14, 2010 by: rch

from etsproxy.traits.api import \
    Float, Str, implements

from etsproxy.traits.ui.ui_traits import Image

from math import pi

from numpy import \
    sign, sqrt, linspace

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf import \
    RF

from matplotlib import pyplot as plt

def H(x):
    return sign(sign(x) + 1.)

class POInfiniteFiber(RF):
    '''
    Pullout of fiber from a stiff matrix;
    stress criterion for debonding, infinite fiber
    '''

    implements(IRF)

    title = Str('pullout - infinite fiber with constant friction')

    xi = Float(0.0179, auto_set = False, enter_set = True, input = True,
                distr = ['weibull_min', 'uniform'])

    tau = Float(2.5, auto_set = False, enter_set = True, input = True,
                distr = ['uniform', 'norm'])

    l = Float(0.0, auto_set = False, enter_set = True, input = True,
              distr = ['uniform'], desc = 'free length')

    D_f = Float(26e-3, auto_set = False, input = True,
              enter_set = True, distr = ['uniform', 'weibull_min'])

    E_f = Float(72.0e3, auto_set = False, enter_set = True, input = True,
                  distr = ['uniform'])

    theta = Float(0.01, auto_set = False, enter_set = True, input = True,
                  distr = ['uniform', 'norm'], desc = 'slack')

    phi = Float(1., auto_set = False, enter_set = True, input = True,
                  distr = ['uniform', 'norm'], desc = 'bond quality')

    u = Float(auto_set = False, enter_set = True,
               ctrl_range = (0.0, 0.04, 100))

    x_label = Str('displacement [mm]')
    y_label = Str('force [N]')

    C_code = Str('')

    def __call__(self, u, tau, l, D_f, E_f, theta, xi, phi):

        A = pi * D_f ** 2 / 4.
        l = l * (1 + theta)
        u = u - theta * l
        T = tau * phi * D_f * pi
        q = (-l * T + sqrt(l ** 2 * T ** 2 + 2 * u * H(u) * E_f * A * T))

        # include breaking strain
        q = q * H(A * E_f * xi - q)
        return q

if __name__ == '__main__':
    q = POInfiniteFiber()
    q.plot(plt, linewidth = 2, color = 'navy')
    plt.show()
