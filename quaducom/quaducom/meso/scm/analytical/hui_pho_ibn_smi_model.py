'''
Created on Apr 20, 2012

model of single fiber composite by Hui, Phoenix, Ibnabeljalil and Smith

@author: rostar
'''

from etsproxy.traits.api import HasTraits, Float, Property, cached_property
import numpy as np
from scipy.special import expi
from scipy.integrate import cumtrapz

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
        Y = s**(self.rho+1.)/2.
        def Phi(x):
            return self.lambd * (0.577215664901532 + np.log(x) - expi(-x))
        def Psi(x):
            return np.exp(-2.*Phi(x))
        return s**(2.*self.rho) * np.exp(self.lambd * s ** (self.rho + 1)) * Psi(Y)

    def p1(self, s, x):
        return self.A0(s) * np.exp(-s**self.rho * x)
    
    def p2(self, s, x):
        def integ_scalar(s, x):
            t = np.linspace(x, s, 300)
            return np.trapz(self.A0(t)/t * np.exp(-t**self.rho * (x + t/2.)), t)
        integ_vect = np.vectorize(integ_scalar)
        return 2. * self.rho * integ_vect(s, x) + np.nan_to_num(self.p1(x, x))

    def p3(self, x):
        return self.p2(2 * x, x)

    def p11(self, s, x):
        return s**(2.*self.rho) * np.exp(self.lambd * s ** (self.rho + 1) - 2.*self.lambd * (0.577215664901532 + np.log(s**(self.rho+1)/2.) - expi(-s**(self.rho+1)/2.)) -s**self.rho * x) 

    def p22(self, s, x):
        def integ_scalar(s, x):
            t = np.linspace(x, s, 300)
            def integrant(t):
                return t**(2.*self.rho) * np.exp(self.lambd * t ** (self.rho + 1) -2.*self.lambd * (0.577215664901532 + np.log(t**(self.rho+1)/2.) - expi(-t**(self.rho+1)/2.)) -t**self.rho * (x + t/2.)) /t
            return np.trapz(integrant(t), t)
        integ_vect = np.vectorize(integ_scalar)
        return 2. * self.rho * integ_vect(s, x) + np.nan_to_num(self.p11(x, x))

    def p33(self, x):
        return self.p22(2 * x, x)

    def p_x(self, s, x):
        p1 = np.nan_to_num(self.p11(s, x)) * (x >= s)
        p2 = np.nan_to_num(self.p22(s, x)) * (x < s) * (x >= s / 2.)
        p3 = np.nan_to_num(self.p33(x)) * (x < s / 2.)
        p = p1 + p2 + p3
        return p

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    sfc = SFC_Hui(l0=1., d=0.007, tau=0.1, sigma0=2200., rho=5.0)
    
    for rho in np.array([500.]):
        sfc.rho = rho
        x = np.linspace(0.4, 1.2, 1000)
        pdf_x = sfc.p_x(4., x)
        pdf_x = pdf_x / np.trapz(pdf_x, x)
        print 2 * np.trapz(pdf_x * x, x)
        # Widom limiting case
        # print 2. / np.trapz(pdf_x, x)
        #cdf_x = np.hstack((0., cumtrapz(pdf_x * x, x)))
        plt.plot(x, pdf_x, label=str(rho))
        #plt.plot(x, cdf_x, label=str(rho))
    #plt.legend()
    plt.show()
    
