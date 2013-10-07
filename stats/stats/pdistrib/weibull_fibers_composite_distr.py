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

    def cdf(self, e, depsf, r):
        '''weibull_fibers_cdf_cb_rigid'''
        m = self.m
        s = ((depsf*(m+1)*self.sV0**m*self.V0)/(2.*pi*r**2))**(1./(m+1))
        return 1. - np.exp(-(e/s) ** (m + 1))


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

    def cdf(self, e, depsf, r, al, ar):
        '''weibull_fibers_cdf_cb_elast'''
        m = self.m
        s = ((depsf*(m+1.)*self.sV0**m*self.V0)/(pi*r**2.))**(1./(m+1.))
        a0 = e/depsf
        expL = (e/s) ** (m + 1) * (1.-(1.-al/a0)**(m+1.))
        expR = (e/s) ** (m + 1) * (1.-(1.-ar/a0)**(m+1.))
        return 1. - np.exp(-expL-expR)


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

    Ll = Float(1e10)
    Lr = Float(1e10)

    def cdf(self, e, depsf, r, al, ar):
        '''weibull_fibers_cdf_mc'''
        Ll, Lr, m = self.Ll, self.Lr, self.m
        s = ((depsf*(m+1.)*self.sV0**m*self.V0)/(pi*r**2.))**(1./(m+1.))
        a0 = e/depsf
        expLfree = (e/s) ** (m + 1) * (1.-(1.-al/a0)**(m+1.))
        expLfixed = a0 / Ll * expLfree
        maskL = al < Ll
        expL = expLfree * maskL + expLfixed * (maskL == False)
        expRfree = (e/s) ** (m + 1) * (1.-(1.-ar/a0)**(m+1.))
        expRfixed = a0 / Lr * expRfree
        maskR = ar < Lr
        expR = expRfree * maskR + expRfixed * (maskR == False)
        return 1. - np.exp(- expL - expR)
        
#    def weibull_fibers_cdf(self, ef0, strain_slope,
#                           shorter_boundary, longer_boundary, r):
#        m = self.shape
#        s = ((strain_slope * (m+1) * self.sV0**m)/(2. * pi * r ** 2))**(1./(m+1))
#        Pf_int = 1 - np.exp(-(ef0/s)**(m+1))
#        I = s * gamma(1 + 1./(m+1)) * gammainc(1 + 1./(m+1), (ef0/s)**(m+1))
#        Pf_broken = I / (m + 1) / ef0
#        return Pf_int - Pf_broken

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    from scipy.integrate import cumtrapz
    r = 0.00345
    tau = 0.5
    Ef = 200e3
    m = 4.0
    sV0 = 0.0026
    e = np.linspace(0.001, 0.04, 100)
    a = e * Ef / (2. * tau / r)
    wfcbe = fibers_CB_elast(m=m, sV0=sV0)
    CDF = wfcbe.cdf(e, 2*tau/Ef/r, r, 0.1 * a, 0.1 * a)
    plt.plot(e, CDF, label='CB')
    wfcbr = fibers_CB_rigid(m=m, sV0=sV0)
    CDF = wfcbr.cdf(e, 2*tau/Ef/r, r)
    plt.plot(e, CDF, label='CB rigid')
    wfmc = fibers_MC(m=m, sV0=sV0, Ll=14.0, Lr=14.0)
    CDF = wfmc.cdf(e, 2*tau/Ef/r, r, a, a)
    plt.plot(e, CDF, label='MC')
    plt.legend()
    plt.show()


