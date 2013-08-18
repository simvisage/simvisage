'''
Created on 17 Aug 2013

@author: Q
'''

from etsproxy.traits.api import HasTraits, Float, Property, cached_property
import numpy as np
from scipy.special import expi
from scipy.integrate import cumtrapz
from scipy.integrate import ode

class SFC_Hui(HasTraits):

    l0 = Float
    d = Float
    tau = Float
    sigma0 = Float
    rho = Float

    def delta(self, sigma):
        '''stress transfer lenth'''
        return self.d * sigma / 2. / self.tau

    def Lambda(self, sigma):
        '''mean number of flaws with strength less than sigma'''
        return (sigma/self.sigma0)**self.rho / self.l0
    
    def F(self, sigma, l):
        '''probability of failure at stress sigma and length l'''
        return 1. - np.exp(-l*self.Lambda(sigma))

    lambd = Property(Float, depends_on='rho')
    @cached_property
    def _get_lambd(self):
        return self.rho/(1 + self.rho)

    sigmac = Property(Float, depends_on='sigma0, l0, tau, d, rho')
    @cached_property
    def _get_sigmac(self):
        '''dimensionless stress'''
        return self.sigma0 * (2. * self.l0 * self.tau / self.d / self.sigma0) ** (self.lambd / self.rho)

    deltac = Property(Float, depends_on='sigma0, l0, tau, d, rho')
    @cached_property
    def _get_deltac(self):
        '''dimensionless length'''
        return self.l0 * (self.d * self.sigma0 / 2. / self.l0 / self.tau) ** self.lambd

    def Lambda_tilde(self, s):
        '''mean number of flaws with strength less than s'''
        return s ** self.rho

    def h(self, s):
        '''hazard rate at s'''
        return self.rho * s ** (self.rho - 1.)
    
    def A0(self, s):
        Y = s**(self.rho+1)/2.
        def Phi(x):
            return self.lambd * (0.577215664901532 + np.log(x) - expi(-x))
        def Psi(x):
            return np.exp(-2.*Phi(x))
        return s**(2.*self.rho) * np.exp(self.lambd * s ** (self.rho + 1)) * Psi(Y)

    def p1(self, s, x):
        return self.A0(s) * np.exp(-s**self.rho * x)
    
    def p2(self, s, x):
        def integ_scalar(s, x):
            t = np.linspace(x, s, 200)
            return np.trapz(self.A0(t)/t * np.exp(-t**self.rho * (x + t/2.)), t)
        integ_vect = np.vectorize(integ_scalar)
        return 2. * self.rho * integ_vect(s, x) + np.nan_to_num(self.p1(x, x))

    def p3(self, x):
        return self.p2(2 * x, x)

    def p11(self, s, x):
        return s**(2.*self.rho) * np.exp(self.lambd * s ** (self.rho + 1) -2.*self.lambd * (0.577215664901532 + np.log(s**(self.rho+1)/2.) - expi(-s**(self.rho+1)/2.)) -s**self.rho * x) 

    def p22(self, s, x):
        def integ_scalar(s, x):
            t = np.linspace(x, s, 200)
            def integrant(t):
                return t**(2.*self.rho) * np.exp(self.lambd * t ** (self.rho + 1) -2.*self.lambd * (0.577215664901532 + np.log(t**(self.rho+1)/2.) - expi(-t**(self.rho+1)/2.)) -t**self.rho * (x + t/2.)) /t
            return np.trapz(integrant(t), t)
        integ_vect = np.vectorize(integ_scalar)
        return 2. * self.rho * integ_vect(s, x) + np.nan_to_num(self.p11(x, x))

    def p33(self, x):
        return self.p22(2 * x, x)

    def p_x(self, s):
        x = self.x
        p1 = np.nan_to_num(self.p11(s, x)) * (x >= s)
        p2 = np.nan_to_num(self.p22(s, x)) * (x < s) * (x >= s / 2.)
        p3 = np.nan_to_num(self.p33(x)) * (x < s / 2.)
        p = p1 + p2 + p3
        #p_inf = np.trapz(p, x)
        return p# / p_inf

    def dp(self, si, p_i):
        x = self.x
        integ = np.trapz(p_i * (x >= si / 2.), x)
        integs = integ - np.hstack((0.0, cumtrapz(p_i * (x >= si / 2.), x)))
        integs = np.hstack((integs[x > si/2.], np.repeat(integs[-1], np.sum(x < si/2.))))
        dp2 = 2 *self.h(si) * integs * (x >= si / 2.) * (x < si)
        dp3 = (2* self.h(si) * integs - (x - si) * p_i * self.h(si)) * (x >= si)  
        return dp2 + dp3

    def p_num(self, s):
        p_i = self.p_x(1.)
        s_arr = np.linspace(1., s, 500)
        ds = s_arr[1] - s_arr[0]
        for si in s_arr:
            p_i += self.dp(si, p_i) * ds
        return p_i
    
    def p_num_scipy(self, s):
        p_0 = self.p_x(1.)
        r = ode(self.dp).set_integrator('dopri5')
        r.set_initial_value(p_0, 1.0)
        t1 = s
        dt = 0.1
        while r.successful() and r.t < t1:
            print 'step'
            r.integrate(r.t+dt)
        return r.y
        

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    import time
    sfc = SFC_Hui(l0=1., d=0.007, tau=0.1, sigma0=2200., rho=5.0)
    
    for rho in np.array([10.]):
        sfc.rho = rho
        sfc.x = np.linspace(0.01, 7., 300)
        #pdf_x = sfc.p_x(2., x)
        #print np.trapz(pdf_x, x)
        #print np.trapz(x*pdf_x_num, x)
        #cdf_x = np.hstack((0., cumtrapz(pdf_x, x)))
        plt.plot(sfc.x, sfc.p_x(1.), label='analytical 1.0')
        plt.plot(sfc.x, sfc.p_x(1.3), label='analytical 1.3')
        t = time.clock()
        euler = sfc.p_num(2.8)
        print 'euler', time.clock() - t
        t = time.clock()
        rk = sfc.p_num_scipy(2.8)
        print 'RK', time.clock() - t
        #plt.plot(sfc.x, euler, label='numerical 1.3')
        plt.plot(sfc.x, rk, label='ru-ku 1.3')
        #plt.ylim(0,4.0)
    plt.legend()
    plt.show()
    