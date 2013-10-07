'''
Created on Jun 6, 2013

@author: rostar
'''

from etsproxy.traits.api import HasTraits, Float
import numpy as np
from math import pi
from scipy.special import gamma, gammainc


def H(x):
    return x >= 0.0


class WeibullFibers(HasTraits):
    '''meta class for the distribution of strength of Weibull fibers'''
    m = Float(5.)
    sV0 = Float(1.)
    V0 = Float(1.)


class fibers_dry(WeibullFibers):
    ''' distribution of strength of brittle fibers with flaws following
        the compound Poisson process.
        input parameters:   m = rate of the Poisson process / Weibull modulus
                            sV0 = scale parameter of flaw strength (in strain per volume unit)
                            V0 = reference volume (default 1.0)
                            r = fiber radius
                            L = fiber length'''  

    def cdf(self, e, r, L):
        s = (self.sV0 ** self.m * self.V0 / (L * r ** 2 * pi)) ** (1./self.m)
        return 1. - np.exp(- (e/s)**self.m)
    
    def pdf(self, e, r, L):
        s = (self.sV0 ** self.m * self.V0 / (L * r ** 2 * pi)) ** (1./self.m)
        return self.m / s * (e / s)**(self.m - 1.) * np.exp(- (e/s)**self.m)
    

class fibers_CB_rigid(WeibullFibers):
    ''' distribution of strength of brittle fibers-in-composite with flaws following
        the compound Poisson process.
        Single crack bridge with rigid matrix = fibers have strain e at the crack plane
        and 0 at the end of the debonded length a.
        input parameters:   m = rate of the Poisson process / Weibull modulus
                            sV0 = scale parameter of flaw strength (in strain per volume unit)
                            V0 = reference volume (default 1.0)
                            r = fiber radius
                            tau = bond strength
                            Ef = fiber modulus of elasticity'''

    def cdf(self, e, r, tau, Ef):
        m = self.m
        T = 2. * tau / r
        s = ((T*(m+1)*self.sV0**m*self.V0)/(2.*Ef*pi*r**2))**(1./(m+1))
        return 1. - np.exp(-(e/s) ** (m + 1))

    def pdf(self, e, r, tau, Ef):
        m = self.m
        T = 2. * tau / r
        s = ((T*(m+1.)*self.sV0**m*self.V0)/(2.*Ef*pi*r**2.))**(1./(m+1.))
        return (m+1.)/s * (e/s)**m * np.exp(- (e/s)**(m+1.))
    
class fibers_CB_elast(WeibullFibers):
    ''' distribution of strength of brittle fibers-in-composite with flaws following
        the compound Poisson process.
        Single crack bridge with elastic matrix = fibers have strain e at the crack plane
        and debond up to the distance a where the fiber strain equals the matrix strain.
        input parameters:   m = rate of the Poisson process / Weibull modulus
                            sV0 = scale parameter of flaw strength (in strain per volume unit)
                            V0 = reference volume (default 1.0)
                            r = fiber radius
                            tau = bond strength
                            Ef = fiber modulus of elasticity
                            a = debonded length'''

    def cdf(self, e, r, tau, Ef, a):
        m = self.m
        T = 2. * tau / r
        s = ((T*(m+1.)*self.sV0**m*self.V0)/(2.*Ef*pi*r**2.))**(1./(m+1.))
        a0 = e/(T/Ef)
        a = a0 * 0.1
        return 1. - np.exp(-(e/s) ** (m + 1) * (1.-(1.-a/a0)**(m+1.)))

    def pdf(self, e, r, tau, Ef, a):
        m = self.m
        T = 2. * tau / r
        s = ((T*(m+1.)*self.sV0**m*self.V0)/(2.*Ef*pi*r**2.))**(1./(m+1.))
        a0 = e/(T/Ef)
        a = a0 * 0.9
        part1 = ((m+1.)*(e/s) ** (m + 1) * (1.-(1.-a/a0)**(m+1.))/e - a * T *(m+1.)*(e/s) ** (m + 1) * (1.-a/a0)**(m+1.)/e**2/Ef )
        part2 = np.exp(-(e/s) ** (m + 1) * (1.-(1.-a/a0)**(m+1.)))
        return part1 * part2


class fibers_MC(WeibullFibers):
    ''' distribution of strength of brittle fibers-in-composite with flaws following
        the compound Poisson process.
        Multiple cracking is considered. In the analysed crack bridge, fibers have symmetry
        condition from left Ll and right Lr. Fiber strain beyond these conditions is assumed
        to be periodical.
        input parameters:   m = rate of the Poisson process / Weibull modulus
                            sV0 = scale parameter of flaw strength (in strain per volume unit)
                            V0 = reference volume (default 1.0)
                            r = fiber radius
                            tau = bond strength
                            Ef = fiber modulus of elasticity
                            Ll = symmetry condition left
                            Lr = symmetry condition right'''

    def cdf(self, e, r, tau, Ef, Ll, Lr):
        print 'if a0 < Ll'
        m = self.m
        T = 2. * tau / r
        s = ((T*(m+1.)*self.sV0**m*self.V0)/(Ef*pi*r**2.))**(1./(m+1.))
        a0 = e/(T/Ef)
        expL = a0 / Ll * (e/s) ** (m + 1) * (1.-(1.-Ll/a0)**(m+1.))
        expR = a0 / Lr * (e/s) ** (m + 1) * (1.-(1.-Lr/a0)**(m+1.))
        return 1. - np.exp(- expL - expR)

#        Ll = shorter_boundary
#        Lr = longer_boundary
#        if np.any(Lr > ef0 / strain_slope + 1e-15):
#            print 'wrong input for Pf (WeibullFibers)'
#        s = (strain_slope ** 2 * (m + 1) * self.sV0 ** m * self.V0 / pi / r ** 2) ** (1./(m+2))
#        a0 = ef0 / strain_slope + 1e-15
#        Pf = 1. - np.exp(- (ef0 / s ) ** (m + 2) *
#            ((1. - (1. - Ll / a0) ** (m + 1))/Ll   +  (1. - (1. - Lr / a0) ** (m + 1))/Lr))
#        return Pf * H(ef0)
#        a = 2 * max_strain / strain_slope
#        return weibull_min(m, scale=(self.V0 * self.sV0 ** m / fiber_radius ** 2 / pi / a)**(1./m)).cdf(max_strain)
        
    def weibull_fibers_cdf(self, ef0, strain_slope,
                           shorter_boundary, longer_boundary, r):
        m = self.shape
        s = ((strain_slope * (m+1) * self.sV0**m)/(2. * pi * r ** 2))**(1./(m+1))
        Pf_int = 1 - np.exp(-(ef0/s)**(m+1))
        I = s * gamma(1 + 1./(m+1)) * gammainc(1 + 1./(m+1), (ef0/s)**(m+1))
        Pf_broken = I / (m + 1) / ef0
        return Pf_int# - Pf_broken

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    from scipy.integrate import cumtrapz
    e = np.linspace(0.001, 0.04, 100)
    wfcbe = fibers_CB_elast(m=4.0, sV0=0.0026)
    CDF = wfcbe.cdf(e, 0.00345, 0.5, 200e3, 1.)
    plt.plot(e, CDF, label='CB')
    wfcbr = fibers_CB_rigid(m=4.0, sV0=0.0026)
    CDF = wfcbr.cdf(e, 0.00345, 0.5, 200e3)
    plt.plot(e, CDF, label='CB rigid')
    wfmc = fibers_MC(m=4.0, sV0=0.0026)
    CDF = wfmc.cdf(e, 0.00345, 0.5, 200e3)
    plt.plot(e, CDF, label='MC')
    plt.legend()
    plt.show()
    
    
    

