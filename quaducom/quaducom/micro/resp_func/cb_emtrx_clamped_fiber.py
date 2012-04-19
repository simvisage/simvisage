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
    Float, Str, implements, cached_property, Property

from etsproxy.traits.ui.ui_traits import Image

from math import pi

import numpy as np
from numpy import \
    sign, sqrt, linspace, minimum, maximum

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf import \
    RF

from matplotlib import pyplot as plt

def H(x):
    return 0.5*(sign(x) + 1.)

class CBEMClampedFiber(RF):
    '''
    Crack bridged by a short fiber with constant
    frictional interface to the elastic matrix; clamped fiber end
    '''

    implements(IRF)

    title = Str('crack bridge - clamped fiber with constant friction')
    image = Image('pics/cb_short_fiber.jpg')


    xi = Float(0.0179, auto_set = False, enter_set = True, input = True,
                distr = ['weibull_min', 'uniform'])

    tau = Float(2.5, auto_set = False, enter_set = True, input = True,
                distr = ['uniform', 'norm'])

    l = Float(10.0, auto_set = False, enter_set = True, input = True,
              distr = ['uniform'], desc = 'free length')

    A_r = Float(.89, auto_set = False, input = True,
              enter_set = True, distr = ['uniform', 'weibull_min'],
              desc = 'CS area of a the reinforcement')

    E_r = Float(72e3, auto_set = False, enter_set = True, input = True,
                  distr = ['uniform'])

    E_m = Float(30e3, auto_set = False, enter_set = True, input = True,
                  distr = ['uniform'])

    A_m = Float(50.0, auto_set = False, enter_set = True, input = True,
              distr = ['uniform'])

    theta = Float(0.01, auto_set = False, enter_set = True, input = True,
                  distr = ['uniform', 'norm'], desc = 'slack')

    phi = Float(1., auto_set = False, enter_set = True, input = True,
                  distr = ['uniform', 'norm'], desc = 'bond quality')

    Ll = Float(1., auto_set = False, enter_set = True, input = True,
              distr = ['uniform'], desc = 'embedded length - left')

    Lr = Float(.5, auto_set = False, enter_set = True, input = True,
              distr = ['uniform'], desc = 'embedded length - right')

    Nf = Float(1000., auto_set = False, enter_set = True, input = True,
              desc = 'number of parallel fibers', distr = ['uniform'])

    w = Float(auto_set = False, enter_set = True, input = True,
               distr = ['uniform'], desc = 'crack width',
               ctrl_range = (0.0, 1.0, 10))

    Kr = Property(depends_on = 'A_r, E_r')
    @cached_property
    def _get_Kr(self):
        return self.A_r * self.E_r

    Km = Property(depends_on = 'A_m, E_m')
    @cached_property
    def _get_Km(self):
        return self.A_m * self.E_m

    Kc = Property(depends_on = 'A_r, E_r, E_m, A_m')
    @cached_property
    def _get_Kc(self):
        return self.Kr + self.Km

    x_label = Str('crack opening [mm]')
    y_label = Str('force [N]')

    C_code = Str('')

    def crackbridge(self, w, l, T, Kr, Km, Lmin, Lmax):
        c = (2 * l * (Km + Kr) + Kr * (Lmin + Lmax))
        P0 = (T * (Km + Kr)) / (2. * Km ** 2) * (sqrt(c ** 2 + 4 * w * Kr * Km ** 2 / T) - c)
        return P0

    def pullout(self, u, l, T, Kr, Km, L):
        c = l * (Km + Kr) + L * Kr
        P1 = (T * (Km + Kr)) / Km ** 2 * (sqrt(c ** 2 + 2 * u * Kr * Km ** 2 / T) - c)
        return P1

    def linel(self, u, l, T, Kr, Km, L):
        P2 = (T * L ** 2 + 2 * u * Kr) / 2. / (L + l)
        return P2

    def __call__(self, w, tau, l, A_r, E_r, E_m, A_m, theta, xi, phi, Ll, Lr, Nf):

        # cross sectional area of a single fiber

        Lmin = minimum(Ll, Lr)
        Lmax = maximum(Ll, Lr)

        Lmin = maximum(Lmin - l / 2., 0)
        Lmax = maximum(Lmax - l / 2., 0)

        l = minimum(l / 2., Lr) + minimum(l / 2., Ll)

        l = l * (1 + theta)
        w = w - theta * l
        w = H(w) * w
        o = sqrt(4. * A_r * Nf ** 2 * pi)
        T = tau * phi * o

        Km = A_m * E_m
        Kr = A_r * Nf * E_r


        # double sided debonding
        l0 = l / 2.
        q0 = self.crackbridge(w, l0, T, Kr, Km, Lmin, Lmax)


        # displacement at which the debonding to the closer clamp is finished
        # the closer distance is min(L1,L2)

        w0 = T * Lmin * ((2 * l0 + Lmin) * (Kr + Km) + Kr * Lmax) / (Km * Kr)
        # force at w0
        Q0 = Lmin * T * (Km + Kr) / (Km)

        # debonding from one side; the other side is clamped
        # equal to a one sided pullout with embedded length Lmax - Lmin and free length 2*Lmin + l


        l1 = 2 * Lmin + l
        L1 = Lmax - Lmin
        q1 = self.pullout(w - w0, l1, T, Kr, Km, L1) + Q0

        # debonding completed at both sides, response becomes linear elastic

        # displacement at which the debonding is finished at both sides
        w1 = (L1 * T * (Kr + Km) * (L1 + l1)) / Kr / Km - T * L1 ** 2 / 2. / Kr

        q2 = self.linel(w - w0, l1, T, Kr, Km, L1) + Q0

        # cut out definition ranges and add all parts 
        q0 = H(w) * (q0 + 1e-15) * H(w0 - w)
        q1 = H(w - w0) * q1 * H(w1 + w0 - w)
        q2 = H(w - w1 - w0) * q2

        q = q0 + q1 + q2

        # include breaking strain
        q = q * H(Kr * xi - q)
        mask0 = q0 * H(Kr * xi - q0) > 0
        mask1 = q1 * H(Kr * xi - q1) > 0
        mask2 = q2 > 0
        #return (mask0, q0[mask0]), (mask1, q1[mask1]), (mask2, (q2 * H(Kr * xi - q2))[mask2]), q
        return q

class CBEMClampedFiberSP(CBEMClampedFiber):
    '''
    stress profile for a crack bridged by a short fiber with constant
    frictional interface to the matrix; clamped fiber end
    '''

    x = Float(0.0, auto_set = False, enter_set = True, input = True,
              distr = ['uniform'], desc = 'distance from crack')

    x_label = Str('position [mm]')
    y_label = Str('force [N]')

    C_code = Str('')

    def __call__(self, w, x, tau, l, A_r, E_r, A_m, E_m, theta, xi, phi, Ll, Lr, Nf):

        o = sqrt(4. * A_r * Nf ** 2 * pi)
        T = tau * phi * o
        Km = A_m * E_m
        Kr = A_r * Nf * E_r

        q = super(CBEMClampedFiberSP, self).__call__(w, tau, l, A_r, E_r, A_m, E_m, theta, xi, phi, Ll, Lr, Nf)
        q_x = q * H(l / 2. - abs(x)) + (q - T * (abs(x) - l / 2.)) * H(abs(x) - l / 2.)
        q_const = q * Kr / (Km + Kr)
        q_x = maximum(q_x, q_const)
        #q_x = q_x * H(x + Ll) * H (Lr - x)
        return q_x

if __name__ == '__main__':

    t = 7.1
    Af = 5.31e-4
    Ef = 72e3
    Am = 50. / 1700
    Em = 30e3
    l = 0.
    theta = 0.01
    xi = 10.0179
    phi = 1.
    Ll = 0.
    Lr = 0.1
    Nf = 1.

    def Pw():
        w = linspace(0, 1.5, 300)
        P = CBEMClampedFiber()
        q = P(w, t, l, Af, Ef, Em, Am, theta, xi, phi, Ll, Lr, Nf)
        plt.plot(w, q, lw = 2, ls = '-', color = 'black', label = 'CB_emtrx')
        plt.xlabel('crack width $w$ [mm]')
        plt.ylabel('filament force $P_\mathrm{f}$ [N]')
        plt.legend(loc = 'best')
        plt.show()

    def SP():
        plt.figure()
        cbcsp = CBEMClampedFiberSP()
        x = linspace(-100, 40, 300)
        q = cbcsp(.1, x, t, l, Af, Ef, Em, Am, theta, xi, phi, Ll, Lr, Nf)
        plt.plot(x, q, lw = 2, color = 'black', label = 'force along filament')
        plt.xticks(fontsize = 14)
        plt.yticks(fontsize = 14)
        plt.legend(loc = 'best')
        plt.ylim(0,60)

    def SP2():
        plt.figure()
        cbcsp = CBEMClampedFiberSP()
        x = linspace(-100, 40, 7)
        print x
        q = cbcsp(.1, x, t, l, Af, Ef, Em, Am, theta, xi, phi, Ll, Lr, Nf)
        print q
        plt.plot(x, q, lw = 2, color = 'black', label = 'force along filament')
        plt.xticks(fontsize = 14)
        plt.yticks(fontsize = 14)
        plt.legend(loc = 'best')
        plt.ylim(0,60)
        plt.show()
    #Pw()
    SP()
    SP2()

