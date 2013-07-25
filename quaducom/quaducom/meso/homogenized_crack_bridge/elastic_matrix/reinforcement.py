
'''
Created on 23.10.2012

An instance or a list of instances of the Reinforcement class
can be used by the composite crack bridge model.

@author: Q
'''

import numpy as np
from spirrid.rv import RV
from etsproxy.traits.api import HasTraits, cached_property, \
    Float, Property, Int, Str, on_trait_change
from types import FloatType
from util.traits.either_type import EitherType
from math import pi
from scipy.special import gamma


def H(x):
    return x >= 0.0

class WeibullFibers(HasTraits):
    '''class evaluating damage for Weibull fibers with linearly decreasing stress'''
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

class Reinforcement(HasTraits):

    label = Str('reinforcement')
    r = EitherType(klasses=[FloatType, RV])
    V_f = Float
    E_f = Float
    xi = EitherType(klasses=[FloatType, RV, WeibullFibers])
    tau = EitherType(klasses=[FloatType, RV])
    n_int = Int

    @on_trait_change('n_int')
    def check(self):
        if self.n_int < 50:
            print 'Warning: integration with', self.n_int, 'points might not be precise enough.'
            print 'A minimum of 50 integration points is recommended'


class ContinuousFibers(Reinforcement):
    '''implements continuous reinforcement'''
    results = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    @cached_property
    def _get_results(self):
        stat_weights = 1.0
        if isinstance(self.tau, RV):
            tau = self.tau.ppf(
                np.linspace(.5 / self.n_int, 1. - .5 / self.n_int, self.n_int))
            stat_weights *= 1. / self.n_int
            nu_r_tau = np.ones_like(tau)
        else:
            tau = self.tau
            nu_r_tau = 1.0
        if isinstance(self.r, RV):
            r = self.r.ppf(
                np.linspace(.5 / self.n_int, 1. - .5 / self.n_int, self.n_int))
            stat_weights *= 1. / self.n_int
            r2 = r ** 2
            nu_r = r2 / np.mean(r2)
        else:
            r = self.r
            r2 = r ** 2
            nu_r = nu_r_tau * 1.0
        if isinstance(tau, np.ndarray) and isinstance(r, np.ndarray):
            r = r.reshape(1, self.n_int)
            tau = tau.reshape(self.n_int, 1)
            nu_r_r = (r2 / np.mean(r2)).reshape(1, self.n_int)
            nu_r_tau = np.ones(self.n_int).reshape(self.n_int, 1)
            nu_r = nu_r_r * nu_r_tau
            r_arr = (nu_r * np.mean(r2))**0.5
            return (2. * tau / r / self.E_f).flatten(), stat_weights, nu_r.flatten(), r_arr.flatten()
        else:
            r_arr = (nu_r * np.mean(r2))**0.5
            return 2. * tau / r / self.E_f, stat_weights, nu_r, r_arr

    depsf_arr = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    @cached_property
    def _get_depsf_arr(self):
        return self.results[0]

    stat_weights = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    @cached_property
    def _get_stat_weights(self):
        return self.results[1]
    
    nu_r = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    @cached_property
    def _get_nu_r(self):
        return self.results[2]
    
    r_arr = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    @cached_property
    def _get_r_arr(self):
        return self.results[3]  
    

