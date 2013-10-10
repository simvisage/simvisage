'''
Created on Jan 15, 2013

@author: rostar
'''
#---#-------------------------------------------------------------------------------
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
    Float, Str, implements, Bool

import numpy as np

from spirrid.i_rf import \
    IRF

from spirrid.rf import \
    RF

from math import pi
from scipy.special import gammainc, gamma
from stats.pdistrib.weibull_fibers_composite_distr import fibers_MC


def H(x):
    return x >= 0.0


class CBResidualRandXi(RF):
    '''
    Crack bridged by a fiber with constant
    frictional interface to rigid; free fiber end;
    '''

    implements(IRF)

    title = Str('crack bridge with rigid matrix')

    tau = Float(2.5, auto_set=False, enter_set=True, input=True,
                distr=['uniform', 'norm'])

    r = Float(0.013, auto_set=False, enter_set=True, input=True,
              distr=['uniform', 'norm'], desc='fiber radius')

    E_f = Float(72e3, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    m = Float(5., auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    sV0 = Float(3.e-3, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    V_f = Float(0.0175, auto_set=False, enter_set=True, input=True,
              distr=['uniform'])

    w = Float(auto_set=False, enter_set=True, input=True,
               distr=['uniform'], desc='crack width',
               ctrl_range=(0.0, 1.0, 10))

    include_pullout = Bool(True)

    x_label = Str('crack opening [mm]')
    y_label = Str('composite stress [MPa]')

    C_code = Str('')

    def __call__(self, w, tau, E_f, V_f, r, m, sV0):
        #strain and debonded length of intact fibers
        T = 2. * tau / r
        #scale parameter with respect to a reference volume
        s = ((T * (m+1) * sV0**m)/(2. * E_f * pi * r ** 2))**(1./(m+1))
        ef0 = np.sqrt(w*T/E_f)
        Gxi = 1 - np.exp(-(ef0/s)**(m+1))
        mu_int = ef0 * (1-Gxi)
        I = s * gamma(1 + 1./(m+1)) * gammainc(1 + 1./(m+1), (ef0/s)**(m+1))
        mu_broken = I / (m+1)
        return (mu_int + mu_broken) * E_f * V_f * r**2

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    cb = CBResidualRandXi()
    w = np.linspace(0.0, 2., 300)
    plt.plot(w, cb(w, .1, 240e3, 0.01, 0.0035, 5.0, 0.0026))
    plt.show()