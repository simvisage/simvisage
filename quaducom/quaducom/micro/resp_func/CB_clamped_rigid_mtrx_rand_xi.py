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

    def cdf(self, e, depsf, r, lm, m, sV0, mask):
        '''weibull_fibers_cdf_mc'''
        s = ((depsf*(m+1.)*sV0**m)/(2.*pi*r**2.))**(1./(m+1.))
        a0 = (e+1e-15)/depsf
        expfree = (e/s) ** (m + 1)
        expfixed = a0 / (lm/2.0) * (e/s) ** (m + 1) * (1.-(1.-lm/2.0/a0)**(m+1.))
        mask = a0 < lm/2.0
        exp = expfree * mask + expfixed * (mask == False)
        return 1. - np.exp(- exp)

    def __call__(self, w, tau, E_f, V_f, r, m, sV0, lm):
        #strain and debonded length of intact fibers
        T = 2. * tau / r
        k = np.sqrt(T/E_f)
        ef0cb = k*np.sqrt(w)  
        ef0lin = w/lm + T*lm/4./E_f
        depsf = T/E_f
        a0 = ef0cb/depsf
        mask = a0 < lm/2.0
        e = ef0cb * mask + ef0lin * (mask == False)
        Gxi = self.cdf(e, depsf, r, lm, m, sV0, mask)
        ef0 = ef0cb * mask + ef0lin * (mask == False)
        mu_int = e * (1.-Gxi)
        APPROXIMATION FROM PHOENIX OR CURTIN
        I = s00 * gamma(1 + 1./(m+1)) * gammainc(1 + 1./(m+1), (ef0/s00)**(m+1))
        #mu_broken = I / 2.
        return (mu_int) * E_f * V_f * r**2

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    cb = CBFinResidualRandXi()
    w = np.linspace(0.0, .8, 300)
    plt.plot(w, cb(w, .1, 240e3, 0.01, 0.0035, 5.0, 0.0026, 100.))
    plt.show()