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

from enthought.traits.api import \
    Float, Str, implements, cached_property, Property

from enthought.traits.ui.ui_traits import Image

from math import pi

from numpy import \
    sign, sqrt, linspace, minimum, maximum

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf import \
    RF

from matplotlib import pyplot as plt

def H(x):
    return sign(sign(x) + 1.)

class CBEMClampedFiber(RF):
    '''
    Crack bridged by a short fiber with constant
    frictional interface to the elastic matrix; clamped fiber end
    '''

    implements(IRF)

    title = Str('crack bridge - clamped fiber with constant friction')
    image = Image('pics/cb_short_fiber.jpg')

    xi = Float(0.0179, auto_set=False, enter_set=True, input=True,
                distr=['weibull_min', 'uniform'])

    tau = Float(2.5, auto_set=False, enter_set=True, input=True,
                distr=['uniform', 'norm'])

    l = Float(0.0, auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='free length')

    A_r = Float(0.89, auto_set=False, input=True,
              enter_set=True, distr=['uniform', 'weibull_min'],
              desc='CS area of a the reinforcement')

    E_r = Float(72.0e3, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    E_m = Float(30.0e3, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    A_m = Float(50.0, auto_set=False, enter_set=True, input=True,
              distr=['uniform'])

    theta = Float(0.01, auto_set=False, enter_set=True, input=True,
                  distr=['uniform', 'norm'], desc='slack')

    phi = Float(1., auto_set=False, enter_set=True, input=True,
                  distr=['uniform', 'norm'], desc='bond quality')

    Ll = Float(1., auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='embedded length - left')

    Lr = Float(.5, auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='embedded length - right')

    Nf = Float(1., auto_set=False, enter_set=True, input=True,
              desc='number of parallel fibers', distr=['uniform'])

    w = Float(auto_set=False, enter_set=True, input=True,
               distr=['uniform'], desc='crack width',
               ctrl_range=(0.0, 1.0, 10))


    Kr = Property(depends_on='A_r, E_r')
    @cached_property
    def _get_Kr(self):
        return self.A_r * self.E_r

    Km = Property(depends_on='A_r, E_m')
    @cached_property
    def _get_Km(self):
        return self.A_m * self.E_m

    Kc = Property(depends_on='A_r, E_r, E_m')
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

    def __call__(self, w, tau, l, A_r, E_f, E_m, A_m, theta, xi, phi, Ll, Lr, Nf):

        # cross sectional area of a single fiber

        Lmin = minimum(Ll, Lr)
        Lmax = maximum(Ll, Lr)

        Lmin = maximum(Lmin - l / 2., 0)
        Lmax = maximum(Lmax - l / 2., 0)

        l = minimum(l / 2., Lr) + minimum(l / 2., Ll)

        l = l * (1 + theta)
        w = w - theta * l
        w = H(w) * w
        D = sqrt(A_r * Nf / pi) * 2
        T = tau * phi * D * pi

        Km = A_m * E_m
        Kr = A_r * E_f


        # double sided debonding
        l0 = l / 2.
        q0 = self.crackbridge(w, l0, T, Kr, Km, Lmin, Lmax)


        # displacement at which the debonding to the closer clamp is finished
        # the closer distance is min(L1,L2)

        w0 = T * Lmin * ((2 * l0 + Lmin) * (Kr + Km) + Kr * Lmax) / (Km * Kr)
        # force at w0
        Q0 = Lmin * T * (Km + Kr) / (Km)
        #print Q0
        # debonding from one side; the other side is clamped
        # equal to a one sided pullout with embedded length Lmax - Lmin and free length 2*Lmin + l


        l1 = 2 * Lmin + l
        L1 = Lmax - Lmin
        q1 = self.pullout(w - w0, l1, T, Kr, Km, L1) + Q0
        

        # debonding completed at both sides, response becomes linear elastic

        # displacement at which the debonding is finished at both sides
        w1 = (L1 * T * (Kr + Km) * (L1 + l1)) / Kr / Km - T * L1 ** 2 / 2. / Kr
        #print 'alt w0', w0, 'alt w1', w1 , 'alt w0+w1', w0 + w1
        q2 = self.linel(w - w0, l1, T, Kr, Km, L1) + Q0
        #print self.linel(0, l1, T, Kr, Km, L1)

        # cut out definition ranges and add all parts 
        q0 = H(w) * q0 * H(w0 - w)
        q1 = H(w - w0) * q1 * H(w1 + w0 - w)
        q2 = H(w - w1 - w0) * q2

        q = q0 + q1 + q2

        # include breaking strain
        q = q * H(Kr * xi - q)
    
        return q

class CBEMClampedFiberSP(CBEMClampedFiber):
    '''
    stress profile for a crack bridged by a short fiber with constant
    frictional interface to the matrix; clamped fiber end
    '''

    x = Float(0.0, auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='distance from crack')

    x_label = Str('position [mm]')
    y_label = Str('force [N]')

    C_code = Str('')

    def __call__(self, w, x, tau, l, A_r, E_f, A_m, E_m, theta, xi, phi, Ll, Lr, Nf):

        D = sqrt(A_r * Nf / pi) * 2
        T = tau * phi * D * pi
        Km = A_m * E_m
        Kr = A_r * E_f

        q = super(CBEMClampedFiberSP, self).__call__(w, tau, l, A_r, E_f, A_m, E_m, theta, xi, phi, Ll, Lr, Nf)
        q_x = q * H(l / 2. - abs(x)) + (q - T * (abs(x) - l / 2.)) * H(abs(x) - l / 2.)
        #q_x = q_x * H(x + Ll) * H (Lr - x)
        a = q * Km / (T * (Km + Kr))
        q_const = (q - T * a)
        q_x = maximum(q_x, q_const)
        return q_x

if __name__ == '__main__':

    t = 0.4
    Af = 3.84e-5
    Ef = 200e3
    Am = 8400
    Em = 25e3
    l = 10.
    theta = 0.0
    xi = 50000
    phi = 1.
    Ll = 10.
    Lr = 10.
    Nf = 2304000.

    def Pw():
        w = linspace(0, .5, 300)
        P = CBEMClampedFiber()
        q = P(w, t, 10., .89, 72e3, 30000., 50., 0.01, 999, 1., 15., 30., 10)
        plt.plot(w, q[0], label='CB')
        plt.legend()
        plt.show()

    def SP():
        cbcsp = CBEMClampedFiberSP()
        x = linspace(-60, 30, 300)
        q = cbcsp(.02, x, t, 10., .89, 72e3, 30000., 50., 0.0, 999, 1., 50., 30., 10)
        plt.plot(x, q, lw=2, color='black', label='force along filament')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(loc='best')
        plt.show()

    def SP2():
        plt.figure()
        cbcsp = CBEMClampedFiberSP()
        x = linspace(-100, 40, 7)
        print x
        q = cbcsp(.1, x, t, l, Af, Ef, Am, Em, theta, xi, phi, Ll, Lr, Nf)
        print q
        plt.plot(x, q, lw = 2, color = 'black', label = 'force along filament')
        plt.xticks(fontsize = 14)
        plt.yticks(fontsize = 14)
        plt.legend(loc = 'best')
        plt.ylim(0,60)
        
    Pw()
    #SP()
    #SP2()
    #plt.show()

