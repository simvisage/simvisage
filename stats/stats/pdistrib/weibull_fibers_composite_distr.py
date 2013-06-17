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
        s = (depsf * (m + 1) * self.sV0 ** m * self.V0 / 2 / pi / r ** 2) ** (1. / (m + 1))
        return s * gamma(1. + 1 / (m + 1.))

    def cdf(self, max_strain, strain_slope, shorter_boundary, longer_boundary, fiber_radius):
        m = self.shape
        x_short = np.hstack((shorter_boundary[1:], np.repeat(shorter_boundary[-1], len(max_strain)-len(shorter_boundary[1:]))))
        x_long = np.hstack((longer_boundary[1:], np.repeat(longer_boundary[-1], len(max_strain)-len(longer_boundary[1:]))))
        s = strain_slope * (m + 1) * self.sV0 ** m * self.V0 / pi / fiber_radius ** 2
        a0 = max_strain / strain_slope
        Pf = 1. - np.exp(- max_strain ** (m + 1) / s *
            (2. - (1 - x_short / a0) ** (m + 1) - (1 - x_long / a0) ** (m + 1)))
        return Pf * H(max_strain)
