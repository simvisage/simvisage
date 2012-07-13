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
    Float, Str, implements

import numpy as np

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf import \
    RF

from matplotlib import pyplot as plt

def H(x):
    return 0.5 * (np.sign(x) + 1.)

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

    E_f = Float(72e3, auto_set=False, enter_set=True, input=True,
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
    
    def crackbridge(self, w, l, T , Kf, Km, Vf):
        #Phase A : Both sides debonding .
        Kc = Kf + Km
        c = Kc * T * l/2.
        q0 = (np.sqrt(c**2 + w*Kf*Km*Kc*T) - c)/Km
        return q0/Vf

    def pullout(self, w, l, T, Kf, Km, Vf, Lmin, Lmax):
        #Phase B : Debonding of shorter side is finished
        Kc = Kf + Km
        c = Kc*T*(Lmin + l)
        f = T**2*Lmin**2*Kc**2
        q1 = (np.sqrt(c ** 2. + f + 2*w*Kc*T*Kf*Km) - c)/Km
        return q1/Vf

    def linel(self, w, l, T, Kf, Km, Vf, Lmax, Lmin):
        #Phase C: Both sides debonded - linear elastic behavior.
        Kc = Kf + Km
        q2 = (2.*w*Kf*Km + T*Kc*(Lmin**2+Lmax**2))/(2.*Km *(Lmax + l + Lmin))
        return q2/Vf

    def __call__(self, w, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r):
        #assigning short and long embedded length
        Lmin = np.minimum(Ll, Lr)
        Lmax = np.maximum(Ll, Lr)
        
        Lmin = np.maximum(Lmin - l / 2., 1e-15)
        Lmax = np.maximum(Lmax - l / 2., 1e-15)

        #maximum condition for free length
        l = np.minimum(l / 2., Lr) + np.minimum(l / 2., Ll)
        
        #defining variables
        l = l * (1 + theta)
        w = w - theta * l
        w = H(w) * w
        T = 2. * tau * V_f / r
        Km = (1. - V_f) * E_m
        Kf = V_f * E_f
        Kc = Km + Kf
        
        # double sided debonding
        #q0 = self.crackbridge(w, l, T, Kr, Km, Lmin, Lmax)
        q0 = self.crackbridge(w, l, T, Kf, Km, V_f)

        # displacement at which the debonding to the closer clamp is finished
        w0 = (Lmin+l)*Lmin*Kc*T/Kf/Km
        
        # debonding of one side; the other side is clamped
        q1 = self.pullout(w, l, T, Kf, Km, V_f, Lmin, Lmax)
        
        # displacement at which the debonding is finished at both sides
        e1 = Lmax*Kc*T/Km/Kf
        w1 = e1*(l+Lmax/2.)+(e1+e1*Lmin/Lmax)*Lmin/2.
        
        # debonding completed at both sides, response becomes linear elastic
        q2 = self.linel(w , l, T, Kf, Km, V_f, Lmax, Lmin)

        # cut out definition ranges 
        q0 = H(w) * (q0 + 1e-15) * H(w0 - w)
        q1 = H(w - w0) * q1 * H(w1 - w)
        q2 = H(w - w1) * q2
        
        #add all parts
        q = q0 + q1 + q2

        # include breaking strain
        q = q * H(E_f * xi - q)
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
            T = 2. * tau * V_f / r        
            q = super(CBEMClampedFiberStressSP, self).__call__(w, tau, l, E_f, E_m, theta, xi, phi, Ll, Lr, V_f, r)            
            #stress in the free length
            q_l = q * H(l / 2 - abs(x))
            
            #stress in the part, where fiber transmits stress to the matrix
            q_e = (q - T/V_f * (abs(x) - l / 2.)) * H(abs(x) - l / 2.)
            #q_e = q_e * H(x + Ll) * H (Lr - x)
            
            #far field stress
            E_c = E_m * (1-V_f) + E_f * V_f
            q_const = q * V_f * E_f / E_c
            
            #putting all parts together
            q_x = q_l + q_e
            q_x = np.maximum(q_x, q_const)
            return q_x

if __name__ == '__main__':

    r = 0.00345
    V_f = 0.0103
    t = .1
    Ef = 200e3
    Em = 25e3
    l = 100.
    theta = 0.0
    xi = 0.017
    phi = 1.
    Ll = 40.
    Lr = 20.
    
    def Pw():
        plt.figure()
        w = np.linspace(0, 1, 300)
        P = CBEMClampedFiberStress()
        q = P(w, t, l, Ef, Em, theta, xi, phi, Ll, Lr, V_f, r) 
        plt.plot(w, q, lw=2, ls='-', color='black', label='CB_emtrx_stress')
        #plt.legend(loc='best')
        plt.ylim(0,)
        plt.ylabel('stress', fontsize=14)
        plt.xlabel('w', fontsize=14)
        plt.title('Pullout Resp Func Clamped Fiber EMTRX')

    def SP():
        plt.figure()
        cbcsp = CBEMClampedFiberStressSP()
        x = np.linspace(-40, 40, 300)
        w = np.linspace(0, 1, 300)[113]
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
