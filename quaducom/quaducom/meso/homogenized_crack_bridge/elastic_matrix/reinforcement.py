

'''
Created on 23.10.2012

An instance or a list of instances of the Reinforcement class
can be used by the composite crack bridge model.

@author: Q
'''

from bmcs.utils.extra_traits.either_type import EitherType
from etsproxy.traits.api import HasTraits, cached_property, \
    Float, Property, Int, Str, on_trait_change, Enum
import numpy as np
from spirrid.rv import RV
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers


class Reinforcement(HasTraits):
    '''common class for all reinforcement types'''
    label = Str('reinforcement')
    r = EitherType(klasses=[float, RV])
    V_f = Float
    E_f = Float
    xi = EitherType(klasses=[float, RV, WeibullFibers])
    tau = EitherType(klasses=[float, RV])
    n_int = Int

    @on_trait_change('n_int')
    def check(self):
        if self.n_int < 50:
            print('Warning: integration with', self.n_int,
                  'points might not be accurate enough')
            print('a minimum of 50 integration points is recommended')


class ContinuousFibers(Reinforcement):
    '''implements continuous reinforcement'''
    results = Property(depends_on='r, V_f, E_f, xi, tau, n_int')

    @cached_property
    def _get_results(self):
        stat_weights = 1.0
        if isinstance(self.tau, RV):
            p_arr = np.linspace(0.5 / self.n_int, 1 -
                                0.5 / self.n_int, self.n_int)
            tau = self.tau.ppf(p_arr)
            mu_tau = np.mean(tau)
            # to eliminate zeros and equal numbers
            tau += np.linspace(mu_tau / 1e6, 2 * mu_tau / 1e6, len(tau))
            nu_r_tau = np.ones_like(tau)
            stat_weights = np.ones_like(tau) / self.n_int
        else:
            tau = self.tau
            nu_r_tau = 1.0
        if isinstance(self.r, RV):
            r = []
            p_arr = np.linspace(.005, 0.995, self.n_int + 1)
            for i, p in enumerate(p_arr[1:]):
                r_arr = self.r.ppf(np.linspace(p_arr[i], p, 500))
                pdf = self.r.pdf(r_arr)
                r.append(np.trapz(r_arr * pdf, r_arr) / (p - p_arr[i]))
            r = np.array(r)

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
            r_arr = (nu_r * np.mean(r2)) ** 0.5
            return (2. * tau / r / self.E_f).flatten(), stat_weights, nu_r.flatten(), r_arr.flatten()
        else:
            r_arr = (nu_r * np.mean(r2)) ** 0.5
            return 2. * tau / r / self.E_f, stat_weights, nu_r, r_arr

    depsf_arr = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    '''Derivatives of fiber strains along z - longitudinal axis.
    Eq. (3)
    '''
    @cached_property
    def _get_depsf_arr(self):
        return self.results[0]

    stat_weights = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    '''Statistical weights corresponding to the randomization.
    '''
    @cached_property
    def _get_stat_weights(self):
        return self.results[1]

    nu_r = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    '''Dimensionless fiber cross section
    Eq. (19)
    '''
    @cached_property
    def _get_nu_r(self):
        return self.results[2]

    r_arr = Property(depends_on='r, V_f, E_f, xi, tau, n_int')
    '''Random filament radius.
    '''
    @cached_property
    def _get_r_arr(self):
        return self.results[3]


class ShortFibers(Reinforcement):
    '''implements short fiber reinforcement'''
    phi = EitherType(
        klasses=[float, RV])  # inclination of short fibers to the crack plane normal
    Lf = Float  # length of short fibers
    bond_law = Enum('plastic', 'elasto-plastic')
    tau_fr = Float  # frictional bond at debonded interface
    k_tau = Float  # stiffness of the elastic adhesive bond
    beta = Float  # slip hardening coefficient
    snub = Float  # snubbing coefficient

    le = Property(depends_on='Lf')

    @cached_property
    def _get_le(self):
        return RV('uniform', loc=0.0, scale=self.Lf / 2.)


if __name__ == '__main__':
    pass
