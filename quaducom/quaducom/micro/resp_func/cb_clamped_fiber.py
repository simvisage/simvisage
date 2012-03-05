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

from enthought.traits.api import Float, Str, implements
from enthought.traits.ui.ui_traits import Image
from math import pi
from matplotlib import pyplot as plt
from numpy import sqrt, linspace, minimum, maximum
from stats.spirrid.i_rf import IRF
from stats.spirrid.rf import RF

def H(x):
    return x >= 0

class CBClampedFiber(RF):
    '''
    Crack bridged by a short fiber with constant
    frictional interface to the matrix; clamped fiber end
    '''

    implements(IRF)

    title = Str('crack bridge - clamped fiber with constant friction')
    image = Image('pics/cb_short_fiber.jpg')

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

    Ll = Float(1., auto_set = False, enter_set = True, input = True,
              distr = ['uniform'], desc = 'embedded length - left')

    Lr = Float(.5, auto_set = False, enter_set = True, input = True,
              distr = ['uniform'], desc = 'embedded length - right')

    w = Float(auto_set = False, enter_set = True, input = True,
               desc = 'crack width',
               ctrl_range = (0, 0.01, 100))

    x_label = Str('crack opening [mm]')
    y_label = Str('force [N]')

    C_code = Str('')

    # TODO: case where Lmin is zero - gives a double sided pullout
    # should be one sided though
    def __call__(self, w, tau, l, D_f, E_f, theta, xi, phi, Ll, Lr):

        A = pi * D_f ** 2 / 4.
        Lmin = minimum(Ll, Lr)
        Lmax = maximum(Ll, Lr)

        Lmin = maximum(Lmin - l / 2., 0)
        Lmax = maximum(Lmax - l / 2., 0)

        l = minimum(Lr + Ll, l)

        l = l * (1 + theta)
        w = w - theta * l

        T = tau * phi * D_f * pi

        # double sided debonding
        l0 = l / 2.
        q0 = (-l0 * T + sqrt((l0 * T) ** 2 + w * H(w) * E_f * A * T))

        # displacement at which the debonding to the closer clamp is finished
        # the closer distance is min(L1,L2)

        w0 = Lmin * T * (Lmin + 2 * l0) / E_f / A

        # debonding from one side; the other side is clamped
        # equal to L1*T + one sided pullout with embedded length Lmax - Lmin and free length 2*L1 + l

        # force at w0
        Q0 = Lmin * T
        l1 = 2 * Lmin + l
        q1 = (-(l1) * T + sqrt((l1 * T) ** 2 +
            2 * (w - w0) * H(w - w0) * E_f * A * T)) + Q0

        # displacement at debonding finished at both sides
        # equals a force of T*(larger clamp distance)


        # displacement, at which both sides are debonded
        w1 = w0 + (Lmax - Lmin) * T * ((Lmax - Lmin) + 2 * (l + 2 * Lmin)) / 2 / E_f / A
        # linear elastic response with stiffness EA/(clamp distance)
        q2 = E_f * A * (w - w1) / (Lmin + Lmax + l) + (Lmax) * T

        q0 = q0 * H(w0 - w)
        q1 = q1 * H(w - w0) * H(w1 - w)
        q2 = q2 * H(w - w1)

        q = q0 + q1 + q2

        # include breaking strain
        q = q * H(A * E_f * xi - q)
        #return q0, q1, q2 * H( A * E_f * xi - q2 ), w0 + theta * l, w1 + theta * l
        return q

class CBClampedFiberSP(CBClampedFiber):
    '''
    stress profile for a crack bridged by a short fiber with constant
    frictional interface to the matrix; clamped fiber end
    '''

    x = Float(0.0, auto_set = False, enter_set = True, input = True,
              distr = ['uniform'], desc = 'distance from crack')

    x_label = Str('position [mm]')
    y_label = Str('force [N]')

    C_code = Str('')

    def __call__(self, w, x, tau, l, D_f, E_f, theta, xi, phi, Ll, Lr):

        T = tau * phi * D_f * pi

        q = super(CBClampedFiberSP, self).__call__(w, tau, l, D_f, E_f, theta, xi, phi, Ll, Lr)
        q_x = q * H(l / 2. - abs(x)) + (q - T * (abs(x) - l / 2.)) * H(abs(x) - l / 2.)
        q_x = q_x * H(x + Ll) * H (Lr - x)
        q_x = q_x * H(q_x)

        return q_x

if __name__ == '__main__':

    def Pw():
        #from numpy import argwhere
        q = CBClampedFiber()
        q.plot(plt, linewidth = 2, color = 'navy')
        plt.show()
    #    i1 = argwhere( w > q[3] )[0]
    #    i2 = argwhere( w > q[4] )[0]
    #    a = i1 - 1
    #    b = i1
    #    c = i2 - 1
    #    d = i2
    #    plt.plot( w[:a], q[0][:a], lw = 2, ls = '-', color = 'black', label = 'double sided pullout' )
    #    plt.plot( w[b:c], q[1][b:c], lw = 2, ls = '--', color = 'black', label = 'one sided pullout' )
    #    plt.plot( w[d:], q[2][d:], lw = 2, ls = 'dotted', color = 'black', label = 'linear elastic' )
    #    plt.legend( loc = 'best' )
    #    #plt.xlabel( 'crack width', fontsize = 20 )
    #    #plt.ylabel( 'force [N]', fontsize = 20 )
    #    plt.xticks( [w[41], w[177]], ['', ''], fontsize = 20 )
    #    plt.yticks( fontsize = 14 )
        plt.show()

    def SP():

        cbcsp = CBClampedFiberSP()
        x = linspace(-60, 30, 300)
        q = cbcsp(.9, x, .01, 3.5, 26e-3, 72e3, 0.0, .0179, 1., 40., 50.)
        plt.plot(x, q, lw = 2, color = 'black', label = 'force along filament')
        plt.xticks(fontsize = 14)
        plt.yticks(fontsize = 14)
        plt.legend(loc = 'best')
        plt.show()

    Pw()
    #SP()
