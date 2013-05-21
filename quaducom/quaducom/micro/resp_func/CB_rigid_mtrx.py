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

from matplotlib import pyplot as plt
from math import pi
from scipy.special import gammainc, gamma


def H(x):
    return x >= 0.0


class CBResidual(RF):
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

    def __call__(self, w, tau, E_f, V_f, r, m, sV0, Pf):
        #strain and debonded length of intact fibers
        T = 2. * tau / r
        ef0_inf = np.sqrt(T * w / E_f)
        #scale parameter with respect to a reference volume
        s0 = ((T * (m+1) * sV0**m)/(2. * E_f * pi * r ** 2))**(1./(m+1))
        # strain at fiber breakage
        ef0_break = s0 * (-np.log(1.-Pf)) ** (1./(m+1))
        # debonded length at fiber breakage
        a_break = ef0_break * E_f / T
        #mean pullout length of broken fibers
        mu_Lpo = a_break / (m + 1)
        # strain carried by broken fibers
        ef0_residual = T / E_f * mu_Lpo

        if self.include_pullout == True:
            ef0_tot = ef0_residual * H(ef0_inf - ef0_break) + ef0_inf * H(ef0_break - ef0_inf)
        else:
            ef0_tot = ef0_inf * H(ef0_break - ef0_inf)
        return ef0_tot * E_f * V_f * r**2

