'''
Created on Jun 6, 2013

@author: rostar
'''

from etsproxy.traits.api import HasTraits, Float
from scipy.special import gamma
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

    def mean(self, depsf, r):
        m = self.shape
        s = (depsf * (m + 1) * self.sV0 ** m * self.V0 / 2 / pi / r **2)**(1./(m+1))
        return s * gamma(1. + 1/(m+1.))

    def weibull_fibers_Pf(self, epsy_arr, depsf, x_short, x_long, r_arr):
        m = self.shape
        x_short = np.hstack((x_short[1:], np.repeat(x_short[-1], len(epsy_arr)-len(x_short[1:]))))
        x_long = np.hstack((x_long[1:], np.repeat(x_long[-1], len(epsy_arr)-len(x_long[1:]))))
        s = depsf * (m + 1) * self.sV0 ** m * self.V0 / pi / r_arr **2
        a0 = epsy_arr / depsf
        Pf = 1. - np.exp( - epsy_arr ** (m+1)/s * (2. - (1-x_short/a0)**(m+1) - (1-x_long/a0)**(m+1)) )
        return Pf * H(epsy_arr)