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
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

def H(x):
    return x >= 0.0


class CBClampedRandXi(RF):
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

    pullout = Bool(False)

    def mu_broken(self, e, depsf, r, lm, m, sV0, mask):
        n = 200
        shape = [1] * len(mask.shape) + [n]
        fact = np.linspace(0.0, 1.0, n).reshape(tuple(shape))
        e = e.reshape(tuple(list(e.shape) + [1]))
        if isinstance(depsf, np.ndarray):
            depsf = depsf.reshape(tuple(list(depsf.shape) + [1]))
        if isinstance(r, np.ndarray):
            r = r.reshape(tuple(list(r.shape) + [1]))
        e_arr = e * fact
        a0 = (e_arr+1e-15)/depsf
        mask = a0 < lm/2.0
        pdf = np.gradient(self.cdf(e_arr, depsf, r, lm, m, sV0), e/float(n))
        if isinstance(pdf, list):
            pdf = pdf[-1]
        e_broken = e_arr/(m+1) * mask + e_arr / 2. * (mask == False)
        return np.trapz(np.nan_to_num(pdf) * e_broken, e_arr)

    def cdf(self, e, depsf, r, lm, m, sV0):
        '''weibull_fibers_cdf_mc'''
        s = ((depsf * (m+1) * sV0**m)/(2. * pi * r ** 2))**(1./(m+1))
        a0 = (e+1e-15)/depsf
        mask = a0 < lm/2.0
        ef0cb = e * mask
        ef0lin = e * (mask == False)
        Gxi_deb = 1 - np.exp(-(ef0cb/s)**(m+1))
        Gxi_clamp = 1 - np.exp(-(ef0lin/s)**(m+1) * (1-(1-lm/(2.*(ef0lin/depsf)))**(m+1)))
        return np.nan_to_num(Gxi_deb) + np.nan_to_num(Gxi_clamp)

    def __call__(self, w, tau, E_f, V_f, r, m, sV0, lm):
        '''free and fixed fibers combined
        the failure probability of fixed fibers
        is evaluated by integrating only
        between -lm/2 and lm/2.
        '''
        T = 2. * tau / r + 1e-10
        k = np.sqrt(T/E_f)
        ef0cb = k*np.sqrt(w)  
        ef0lin = w/lm + T*lm/4./E_f
        depsf = T/E_f
        a0 = ef0cb/depsf
        mask = a0 < lm/2.0
        e = ef0cb * mask + ef0lin * (mask == False)
        Gxi = self.cdf(e, depsf, r, lm, m, sV0)
        mu_int = e * (1-Gxi)
        if self.pullout:
            mu_broken = self.mu_broken(e, depsf, r, lm, m, sV0, mask)
            return (mu_int + mu_broken) * E_f * V_f * r**2
        else:
            return mu_int * E_f * V_f * r**2

    def free_deb(self, w, tau, E_f, V_f, r, m, sV0):
        '''free debonding only = __call__ with lm=infty'''
        T = 2. * tau / r
        #scale parameter with respect to a reference volume
        s = ((T * (m+1) * sV0**m)/(2. * E_f * pi * r ** 2))**(1./(m+1))
        ef0 = np.sqrt(w*T/E_f)
        Gxi = 1 - np.exp(-(ef0/s)**(m+1))
        mu_int = ef0 * (1-Gxi)
        I = s * gamma(1 + 1./(m+1)) * gammainc(1 + 1./(m+1), (ef0/s)**(m+1))
        mu_broken = I / (m+1)
        if self.pullout:
            return (mu_int + mu_broken) * E_f * V_f * r**2
        else:
            return mu_int * E_f * V_f * r**2

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    cb = CBClampedRandXi(pullout=True)
    w = np.linspace(0.0, 2., 300)
    sigma = []
    for wi in w:
        sigma.append(cb(wi, .1, 240e3, 0.01, 0.0035, 5.0, 0.0026, 1000.))
    plt.plot(w, np.array(sigma)/0.0035**2)
    plt.show()