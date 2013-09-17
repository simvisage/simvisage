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
    '''the class contains a method for the evaluation of
    failure probability of continuous fibers in a composite
    with Weibull local strength distribution'''
    shape = Float(5.0)
    sV0 = Float
    V0 = Float(1.)

    def weibull_fibers_cdf(self, ef0, strain_slope,
                           shorter_boundary, longer_boundary, r):
        m = self.shape
        s = ((strain_slope * (m+1) * self.sV0**m)/(2. * pi * r ** 2))**(1./(m+1))
        Pf_int = 1 - np.exp(-(ef0/s)**(m+1))
        I = s * gamma(1 + 1./(m+1)) * gammainc(1 + 1./(m+1), (ef0/s)**(m+1))
        Pf_broken = I / (m + 1) / ef0
        return Pf_int - Pf_broken

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

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    wf = WeibullFibers(shape=4.0, sV0=0.5)
    CDF = wf.weibull_fibers_cdf(0.03, 0.0015, 20., 20., 0.00345)
    print CDF
    
    
    

