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


def H(x):
    return x >= 0.0


class CBFinResidualRandXi(RF):
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

    lm = Float(np.inf, auto_set=False, enter_set=True, input=True,
              distr=['uniform'])

    w = Float(auto_set=False, enter_set=True, input=True,
               distr=['uniform'], desc='crack width',
               ctrl_range=(0.0, 1.0, 10))

    x_label = Str('crack opening [mm]')
    y_label = Str('composite stress [MPa]')

    C_code = Str('')

    def __call__(self, w, tau, E_f, V_f, r, m, sV0, lm):
        #strain and debonded length of intact fibers
        T = 2. * tau / r
        #scale parameter with respect to a reference volume
        s0cb = ((T * (m+1) * sV0**m)/(2. * E_f * pi * r ** 2))**(1./(m+1))
        k = np.sqrt(T/E_f)
        ef0cb = k*np.sqrt(w)
        Gxicb = 1 - np.exp(-(ef0cb/s0cb)**(m+1))
        
        s0lin = ((T**2 * (m+1) * sV0**m * lm)/(4. * E_f**2 * pi * r ** 2))**(1./(m+2))
        ef0lin = (w - T*lm**2/4./E_f)/lm + T*lm/2./E_f
        a = ef0lin * E_f / T
        Gxilin = 1. - np.exp(-(ef0lin/s0lin)**(m+2) * (1.-(1.-lm/2./a)**(m+1)))
        
        ef0 = ef0cb * H(w) * H(T*lm**2/4./E_f - w) + ef0lin * H(w - T*lm**2/4./E_f)
        Gxi = Gxicb * H(w) * H(T*lm**2/4./E_f - w) + Gxilin * H(w - T*lm**2/4./E_f)
        mu_int = ef0 * (1.-Gxi)
        s00 = ((T * sV0**m)/(2. * E_f * pi * r ** 2))**(1./(m+1))
        I = s00 * gamma(1 + 1./(m+1)) * gammainc(1 + 1./(m+1), (ef0/s00)**(m+1))
        mu_broken = I / 2.
        return (mu_int + mu_broken) * E_f * V_f * r**2

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    cb = CBFinResidualRandXi()
    w = np.linspace(0.0, .8, 300)
    plt.plot(w, cb(w, .1, 240e3, 0.01, 0.0035, 5.0, 0.0026, 100.))
    plt.show()