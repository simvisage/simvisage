# -------------------------------------------------------------------------------
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

from etsproxy.traits.api import \
    Float, Str, implements
import numpy as np
from spirrid.i_rf import IRF
from spirrid.rf import RF


class CBShortFiber(RF):
    '''
    Crack bridged by a short fiber with constant
    frictional interface to the matrix
    '''

    implements(IRF)

    title = Str('crack bridge - short fiber with constant friction')

    xi = Float(.1, auto_set=False, enter_set=True, input=True,
                distr=['weibull_min', 'uniform'])

    E_f = Float(200e+3 , auto_set=False, enter_set=True,
                desc='filament stiffness [N/mm2]',
                distr=['uniform', 'norm'],
                scale=210e3, shape=0)

    r = Float(0.3, auto_set=False, enter_set=True,
                desc='fiber radius [mm]',
                distr=['uniform', 'norm'],
                scale=0.5, shape=0)

    le = Float(8.5, auto_set=False, enter_set=True,
                desc='shorter embedded length [mm]',
                distr=['uniform'],
                scale=8.5, shape=0)

    L_f = Float(17.0, auto_set=False, enter_set=True,
                desc='fiber length [mm]',
                distr=['uniform', 'norm'],
                scale=30, shape=0)

    tau = Float(1.76, auto_set=False, enter_set=True,
                desc='bond shear stress [N/mm2]',
                distr=['norm', 'uniform', 'weibull_min'],
                scale=1.76, shape=0.5)

    snub = Float(0.03, auto_set=False, enter_set=True,
            desc='snubbing coefficient',
            distr=['uniform', 'norm'],
                scale=0.05, shape=0)

    spall = Float(0.03, auto_set=False, enter_set=True,
            desc='spalling coefficient',
            distr=['uniform', 'norm'],
                scale=0.05, shape=0)

    phi = Float(0.0, auto_set=False, enter_set=True,
       desc='inclination angle',
       distr=['sin2x', 'uniform'],
                scale=1.0, shape=0)

    w = Float(ctrl_range=(0, 0.01, 100), auto_set=False, enter_set=True)
    x_label = Str('crack_opening [mm]', enter_set=True, auto_set=False)
    y_label = Str('force [N]', enter_set=True, auto_set=False)

    C_code = ''

    def __call__(self, w, tau, r, E_f, le, phi, snub, xi, spall):
        T = 2. * tau / r
        # debonding stage
        ef0_deb = np.sqrt(T * w / E_f)
        ef0_deb *= np.exp(phi*snub) * np.cos(phi)**spall
        mask_deb = ef0_deb < xi
        ef0_deb = ef0_deb * mask_deb
        # crack opening at which debonding is finished
        w0 = le ** 2 * T / E_f
        mask_pullout = np.sqrt(T * w0 / E_f) < xi
        # pulling out stage - the fiber is pulled out from the
        # side with the shorter embedded length only
        ef0_pull = (le + w0 - w) * T / E_f
        ef0_pull *= np.exp(phi*snub) * np.cos(phi)**spall
        ef0_pull = ef0_pull * mask_pullout # fibers have ruptured during debonding
        
        ef0 = ef0_deb * (w < w0) + ef0_pull * (w > w0)
        # include breaking strain
        return ef0

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    frc = CBShortFiber()
    w = np.linspace(0.0, .06, 200)
    ef0 = frc(w, .3, 0.013, 70e3, 7.0, 0.5, 0.03, 20.0, 0.)
    plt.plot(w, ef0 * 70e3, linewidth=2, color='navy')
    #plt.ylim(0,0.01 * 70e3)
    plt.show()
