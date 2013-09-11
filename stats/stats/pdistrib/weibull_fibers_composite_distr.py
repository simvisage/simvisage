'''
Created on Jun 6, 2013

@author: rostar
'''

from etsproxy.traits.api import HasTraits, Float
from scipy.stats import weibull_min
import numpy as np
from math import pi


def H(x):
    return x >= 0.0 


class WeibullFibers(HasTraits):
    '''the class contains a method for the evaluation of
    failure probability of continuous fibers in a composite
    with Weibull local strength distribution'''
    shape = Float(5.0)
    sV0 = Float
    V0 = Float(1.)

    def weibull_fibers_cdf(self, max_strain, strain_slope,
                           shorter_boundary, longer_boundary, fiber_radius):
        m = self.shape
#        Ll = shorter_boundary
#        Lr = longer_boundary
#        if np.any(Lr > max_strain / strain_slope + 1e-15):
#            print 'wrong input for Pf (WeibullFibers)'
#        s = (strain_slope ** 2 * (m + 1) * self.sV0 ** m * self.V0 / pi / fiber_radius ** 2) ** (1./(m+2))
#        a0 = max_strain / strain_slope + 1e-15
#        Pf = 1. - np.exp(- (max_strain / s ) ** (m + 2) *
#            ((1. - (1. - Ll / a0) ** (m + 1))/Ll   +  (1. - (1. - Lr / a0) ** (m + 1))/Lr))
#        return Pf * H(max_strain)
        return weibull_min(m, scale=(self.V0 * self.sV0 ** m / fiber_radius ** 2 / pi / 10.)**(1./m)).cdf(max_strain)

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    wf = WeibullFibers(shape=4.0, sV0=0.5)
    CDF = wf.weibull_fibers_cdf(0.03, 0.0015, 20., 20., 0.00345)
    print CDF
    
    
    

