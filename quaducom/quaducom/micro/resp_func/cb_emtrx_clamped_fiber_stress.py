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
    Float, Str, implements, cached_property, Property

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
    return 0.5 * (sign(x) + 1.)

class CBEMClampedFiberStress(RF):
    '''
    Crack bridged by a fiber with constant
    frictional interface to the elastic matrix; clamped fiber end;
    Gives tension.
    '''

    implements(IRF)

    title = Str('crack bridge - clamped fiber with constant friction')

    xi = Float(0.0179, auto_set=False, enter_set=True, input=True,
                distr=['weibull_min', 'uniform'])

    tau = Float(2.5, auto_set=False, enter_set=True, input=True,
                distr=['uniform', 'norm'])

    l = Float(10.0, auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='free length')

    r = Float(0.013, auto_set=False, input=True,
              enter_set=True, desc='fiber radius in mm')

    E_r = Float(72e3, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    E_m = Float(30e3, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    V_f = Float(0.0175, auto_set=False, enter_set=True, input=True,
              distr=['uniform'])

    theta = Float(0.01, auto_set=False, enter_set=True, input=True,
                  distr=['uniform', 'norm'], desc='slack')

    phi = Float(1., auto_set=False, enter_set=True, input=True,
                  distr=['uniform', 'norm'], desc='bond quality')

    Ll = Float(1., auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='embedded length - left',
               ctrl_range=(0.0, 1.0, 10))

    Lr = Float(.5, auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='embedded length - right',
               ctrl_range=(0.0, 1.0, 10))

    w = Float(auto_set=False, enter_set=True, input=True,
               distr=['uniform'], desc='crack width',
               ctrl_range=(0.0, 1.0, 10))

    x = Float(auto_set=False, enter_set=True, input=True,
               distr=['uniform'], desc='crack width',
               ctrl_range=(0.0, 1.0, 10))

    x_label = Str('crack opening [mm]')
    y_label = Str('force [N]')

    C_code = Str('')

    def crackbridge(self, w, l, T, Kr, Km, Lmin, Lmax):
        #Phase A : Both sides debonding .
        Ec = Kr + Km
        c = Kr * T * (Lmin + Lmax) + l * T * Ec
        t1 = 0.5 * (Kr * Ec / Km ** 2)
        t2 = 4 * Km ** 2. * w * T
        P0 = t1 * ((c ** 2 + t2) ** 0.5 - c)
        return P0

    def pullout(self, w, l, T, Kr, Km, Lmin, Lmax):
        #Phase B : Debonding of shorter side is finished
        Ec = Kr + Km
        c = T * Kr * Lmax + T * Ec * (Lmin + l)
        t1 = Ec * Kr / Km ** 2.
        t2 = 2. * Km ** 2. * w * T
        t3 = (Km ** 2.) * (Lmin ** 2.) * (T ** 2.)
        P1 = t1 * (sqrt(c ** 2. + t2 + t3) - c)
        return P1

    def linel(self, w, l, T, Kr, Km, Lmax, Lmin):
        #Phase C: Both sides debonded - linear elastic behavior.
        t1 = Lmax ** 2. + Lmin ** 2.
        P2 = 0.5 * (2. * w + T * t1) * Kr / (Lmax + l + Lmin)
        return P2

    def __call__(self, w, tau, l, E_r, E_m, theta, xi, phi, Ll, Lr, V_f, r):
        #assigning short and long embedded length
        Lmin = minimum(Ll, Lr)
        Lmax = maximum(Ll, Lr)
        
        Lmin = maximum(Lmin - l / 2., 0.)
        Lmax = maximum(Lmax - l / 2., 0.)

        #maximum condition for free length
        l = minimum(l / 2., Lr) + minimum(l / 2., Ll)
        
        #defining variables
        l = l * (1 + theta)
        w = w - theta * l
        w = H(w) * w
        T = 2. * tau / (r * E_r)
        Km = (1. - V_f) * E_m
        Kr = V_f * E_r
        Ec = Km + Kr

        # double sided debonding
        q0 = self.crackbridge(w, l, T, Kr, Km, Lmin, Lmax)

        # displacement at which the debonding to the closer clamp is finished
        w0 = T * Lmin * (Lmin * Km + Kr * (Lmin + Lmax) + l * Ec) / Km
        
        # debonding of one side; the other side is clamped
        q1 = self.pullout(w , l, T, Kr, Km, Lmin, Lmax) 
        
        # displacement at which the debonding is finished at both sides
        w1 = (1. / 2.) * T / Km * (2. * Kr * Lmax ** 2. - Km * Lmin ** 2. + 2. * Lmin * Ec * Lmax + 2 * l * Ec * Lmax + Km * Lmax ** 2)
        
        # debonding completed at both sides, response becomes linear elastic
        q2 = self.linel(w , l, T, Kr, Km, Lmax, Lmin) 
        
        # cut out definition ranges 
        q0 = H(w) * (q0 + 1e-15) * H(w0 - w)
        q1 = H(w - w0) * q1 * H(w1 - w)
        q2 = H(w - w1) * q2
        
        #add all parts
        q = q0 + q1 + q2

        # include breaking strain
        q = q * H(Kr * xi - q)
   
        return q
    
class CBEMClampedFiberStressSP(CBEMClampedFiberStress):
        '''
        stress profile for a crack bridged by a fiber with constant
        frictional interface to the matrix; clamped fiber end
        '''
        x = Float(0.0, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'], desc='distance from crack')
    
        x_label = Str('position [mm]')
        y_label = Str('force [N]')
    
        C_code = Str('')
        
    

        def __call__(self, w, x, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r):
            T = 2. * tau / r        
            q = super(CBEMClampedFiberStressSP, self).__call__(w, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r)            
            #stress in the free length
            q_l = q / V_f * H(l / 2 - abs(x))
            
            #stress in the part, where fiber transmits stress to the matrix
            q_e = (q / V_f - T * (abs(x) - l / 2.)) * H(abs(x) - l / 2.)
            #q_e = q_e * H(x + Ll) * H (Lr - x)
            
            #far field stress
            E_c = E_m * (1-V_f) + E_f * V_f
            q_const = q * E_f/E_c
            
            #putting all parts together
            q_x = q_l + q_e
            q_x = maximum(q_x, q_const)
            return q_x

if __name__ == '__main__':

    r = 0.00345
    V_f = 0.0103
    t = 0.5
    Ef = 200e3
    Em = 25e3
    l = 10.
    theta = 0.0
    xi = 0.017
    phi = 1.
    Ll = 20.
    Lr = 20.
    
    def Pw():
        plt.figure()
        w = linspace(0, 1, 300)
        P = CBEMClampedFiberStress()
        q = P(w, t, l, Ef, Em, theta, xi, phi, Ll, Lr, V_f, r) 
        plt.plot(w, q , lw=2, ls='-', color='black', label='CB_emtrx_stress')
        #plt.legend(loc='best')
        plt.ylim(0,)
        plt.ylabel('stress ', fontsize=14)
        plt.xlabel('w', fontsize=14)
        plt.title('Pullout Resp Func Clamped Fiber EMTRX')

    def SP():
        plt.figure()
        cbcsp = CBEMClampedFiberStressSP()
        x = linspace(-40, 40, 300)
        w = linspace(0, 1, 300)[113]
        q = cbcsp(w, x, t, l, Ef, Em, theta, xi, phi, Ll, Lr, V_f, r)
        plt.plot(x, q, lw=2, color='black', label='stress along filament')
        plt.ylabel('stress', fontsize=14)
        plt.xlabel('position', fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.title('Stress Along Filament EMTRX at w = %.1f' %w)
        #plt.legend(loc='best')
        plt.ylim(0,)

    Pw()
    SP()
    plt.show()
