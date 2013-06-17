
'''
Created on 23.10.2012

An instance or a list of instances of the Reinforcement class
can be used by the composite crack bridge model.

@author: Q
'''

import numpy as np
from spirrid.rv import RV
from etsproxy.traits.api import HasTraits, cached_property, \
    Float, Property, Int, Str
from types import FloatType
from util.traits.either_type import EitherType
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers


class Reinforcement(HasTraits):
    '''common class for all reinforcement types'''
    label = Str('reinforcement')
    r = EitherType(klasses=[FloatType, RV])
    V_f = Float
    E_f = Float
    xi = EitherType(klasses=[FloatType, RV, WeibullFibers])
    tau = EitherType(klasses=[FloatType, RV])
    n_int = Int


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


class ShortFibers(Reinforcement):
    '''implements short fiber reinforcement'''
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
