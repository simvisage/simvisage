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

from matplotlib import pyplot as plt
from etsproxy.traits.api import \
    Float, Str, implements

from etsproxy.traits.ui.ui_traits import Image

from etsproxy.traits.ui.menu import OKButton, CancelButton

from etsproxy.traits.ui.api import \
    View, Item

from math import e, pi
from numpy import sqrt, linspace, sign, abs, cos
from spirrid.i_rf import IRF
from spirrid.rf import RF



def H(x):
    return sign(sign(x) + 1.)

class CBShortFiber(RF):
    '''
    Crack bridged by a short fiber with constant
    frictional interface to the matrix
    '''

    implements(IRF)

    title = Str('crack bridge - short fiber with constant friction')
    image = Image('pics/cb_short_fiber.jpg')

    xi = Float(.1, auto_set=False, enter_set=True, input=True,
                distr=['weibull_min', 'uniform'])

    E_f = Float(200e+3 , auto_set=False, enter_set=True,
                desc='filament stiffness [N/mm2]',
                distr=['uniform', 'norm'],
                scale=210e3, shape=0)

    D_f = Float(0.3, auto_set=False, enter_set=True,
                desc='filament diameter [mm]',
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
                distr=['norm', 'uniform'],
                scale=1.76, shape=0.5)

    f = Float(0.03, auto_set=False, enter_set=True,
            desc='snubbing coefficient',
            distr=['uniform', 'norm'],
                scale=0.05, shape=0)

    phi = Float(0.0, auto_set=False, enter_set=True,
       desc='inclination angle',
       distr=['sin2x', 'sin_distr'],
                scale=1.0, shape=0)

    l = Float(0.0, auto_set=False, enter_set=True,
              distr=['uniform'], desc='free length')

    theta = Float(0.0, auto_set=False, enter_set=True,
                  distr=['uniform', 'norm'], desc='slack')

    w = Float(ctrl_range=(0, 0.01, 100), auto_set=False, enter_set=True)

    x_label = Str('crack_opening [mm]', enter_set=True, auto_set=False)
    y_label = Str('force [N]', enter_set=True, auto_set=False)

    C_code = ''


    def __call__(self, w, tau, L_f, D_f, E_f, le, phi, f, l, theta, xi):

        l = l * (1 + theta)
        w = w - theta * l
        T = tau * pi * D_f
        E = E_f
        A = D_f ** 2 / 4. * pi

        # debonding stage
        q_deb = -l * T + sqrt((l * T) ** 2 + E * A * T * w * H(w))

        # displacement at which debonding is finished
        w0 = le * T * (le + 2 * l) / E_f / A

        # pulling out stage - the fiber is pulled out from the
        # side with the shorter embedded length only

        q_pull = le * T * ((w0 - w) / ((le + 1e-15) - w0) + 1)
        q = q_deb * H(le * T - q_deb) + q_pull * H(q_deb - le * T)

        # include inclination influence
        q = q * H(q) * e ** (f * phi)

        # include breaking strain
        q = q * H(A * E_f * xi - q)
        return q


    traits_view2 = View(Item('E_f', label='fiber E-mod'),
                        Item('D_f', label='fiber diameter'),
                        Item('f' , label='snubbing coef.'),
                        Item('phi', label='inclination angle'),
                        Item('le', label='shorter embedded length'),
                        Item('tau', label='frictional coef.'),
                        resizable=True,
                        scrollable=True,
                        height=0.8, width=0.8,
                        buttons=[OKButton, CancelButton]
                        )

if __name__ == '__main__':
    q = CBShortFiber()
    # q.configure_traits(view='traits_view2')
    q.plot(plt, linewidth=2, color='navy')
    plt.show()
