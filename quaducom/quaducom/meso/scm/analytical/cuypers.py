# -------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Sep 21, 2009 by: rch

from etsproxy.traits.api import HasTraits, Float, Property, \
    cached_property, Range, Instance
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import weibull_min
from scipy.special import gamma
from .material import Material
from .reinf_cross_section import SimplyRatio, GridReinforcement
from util.traits.either_type import EitherType
from math import cos, pi as Pi, e


class Cuypers(HasTraits):

    '''
    Stochastic cracking model due to H. Cuypers (W. Curtin)
    matrix tensile strength is assumed to follow the Weibull distribution
    '''

    material = Instance(Material)

    def _material_default(self):
        return Material()

    m = Float(5.3, auto_set=False, enter_set=True,  # [-]
              desc='Weibull shape parameter for the matrix tensional strength [-]',
              modified=True)

    orientation = Range(low=0.0, high=Pi / 2, value=0.0,
                        auto_set=False, enter_set=True,  # [-]
                        desc='fiber orientation [rad]',
                        modified=True)

    sigma_mu = Float(12.0, auto_set=False, enter_set=True,  # [N/mm^2]
                     desc='mean matrix tensional strength [MPa]',
                     modified=True)

    sigma_fu = Float(720.0, auto_set=False, enter_set=True,  # [N/mm^2]
                     desc='fiber tensional strength [MPa]',
                     modified=True)

    reinf_ratio = EitherType(names=['grid fiber layout', 'explicit value'],
                             klasses=[GridReinforcement, SimplyRatio],
                             modified=True)

    Pf = Range(low=0.0, high=1.0 - 1e-15, value=0.9,
               auto_set=False, enter_set=True,  # [-]
               desc='probability of crack spacing to be of final range',
               modified=True)

    sigma_ref = Range(low=1e-10, high=80.0, value=12.0,
                      auto_set=False, enter_set=True,  # [MPa]
                      desc='reference matrix stress', modified=True)

    sigma_ref_2 = Property(Float, depends_on='m, sigma_mu')

    def _get_sigma_ref_2(self):
        return e ** (0.0991 / self.m) * self.scale_sigma_m

    rho = Property(Float, depends_on='reinf_ratio.rho,orientation')

    def _get_rho(self):
        if self.reinf_ratio.rho * cos(self.orientation) == 0:
            return 1e-15
        else:
            return self.reinf_ratio.rho * cos(self.orientation)

    V_f = Property(Float, depends_on='reinf_ratio.rho, orientation')

    @cached_property
    def _get_V_f(self):
        return self.rho

    V_m = Property(Float, depends_on='reinf_ratio.rho, orientation')

    @cached_property
    def _get_V_m(self):
        return 1 - self.rho

    alpha = Property(Float, depends_on='E_m,E_f,reinf_ratio.rho, orientation')

    @cached_property
    def _get_alpha(self):
        return (self.material.E_m * self.V_m) / (self.material.E_f * self.V_f)

    E_c = Property(Float, depends_on='E_m,E_f,reinf_ratio.rho, orientation')

    @cached_property
    def _get_E_c(self):
        return self.material.E_f * self.V_f + self.material.E_m * self.V_m

    delta_final = Property(Float, depends_on='E_m,E_f,rho,r,sigma_ref,tau,m')

    @cached_property
    def _get_delta_final(self):
        return self.sigma_ref * (self.V_m * self.material.r) / (self.V_f * 2 * self.material.tau)

    cs_final = Property(Float)

    def _get_cs_final(self):
        return 1.337 * self.delta_final

    def _get_delta(self, sigma_c):
        return sigma_c * (self.V_m * self.material.r * self.material.E_m) / (self.V_f * 2 * self.material.tau * self.E_c)

    # matrix strength scale parameter for the Weibull distribution with sigma_mu as
    # mean and m as shape parameter
    scale_sigma_m = Property(Float, depends_on='sigma_mu, m')

    @cached_property
    def _get_scale_sigma_m(self):
        return self.sigma_mu / gamma(1. + 1. / self.m)

    # composite scale parameter for the Weibull distribution with sigma_mu as
    # mean and m as shape parameter
    scale_sigma_c = Property(
        Float, depends_on='sigma_mu, m, E_m, E_f, reinf_ratio.rho, orientation')

    @cached_property
    def _get_scale_sigma_c(self):
        return self.scale_sigma_m / self.material.E_m * self.E_c

    def _get_cs(self, sigma_c):
        Pf = weibull_min.cdf(sigma_c, self.m, scale=self.scale_sigma_c)
        if Pf == 0:
            Pf = 1e-15
        return self.cs_final * 1.0 / Pf

    def eps_c(self, sigma_c):
        cs = self._get_cs(sigma_c)
        delta = self._get_delta(sigma_c)
        if cs > 2 * delta:
            return sigma_c / self.E_c * (1 + self.alpha * delta / cs)
        else:
            return sigma_c * (1. / (self.material.E_f * self.V_f) -
                              (self.alpha * cs) / (4. * delta * self.E_c))

    def _get_epsilon_c(self, sigma_c):
        get_epsilon_c = np.vectorize(self.eps_c)
        return get_epsilon_c(sigma_c)

    sigma_cu = Property(depends_on='rho, orientation, sigma_fu')

    @cached_property
    def _get_sigma_cu(self):
        '''Ultimate composite strength.
        The strength is given by the fiber strength related to the composite
        cross section by the reinforcement ratio rho and projected
        by the cosine of the fiber inclination into the loading direction.
        '''
        # 0.05 quantile strength of the matrix as a criterion for matrix failure
        # when this value is higher than the composite strength governed by the
        # reinforcement stress
        quantile = weibull_min.ppf(0.05, self.m, scale=self.scale_sigma_m)
        # composite failure due to matrix failure
        if self.sigma_fu * self.V_f * cos(self.orientation) < quantile / self.V_m:
            return quantile / self.V_m
        else:
            # composite failure due to fiber failure
            return self.sigma_fu * self.V_f * cos(self.orientation)

    csf = Property(
        depends_on='m, sigma_mu, sigma_ref, Pf, sigma_fu, orientation, E_f, E_m')

    @cached_property
    def _get_csf(self):
        # composite stress at Pf probability for CS to be of final range
        # scale parameter for composite stress
        scale_sigma_c = self.scale_sigma_m / self.E_m * self.E_c

        sigma = weibull_min.ppf(self.Pf, self.m, scale=scale_sigma_c)

        # point of reaching final crack spacing
        epsilon_csf = np.hstack((0,
                                 self._get_epsilon_c(sigma),
                                 self._get_epsilon_c(sigma)))
        sigma_csf = np.hstack((sigma, sigma, 0))
        return epsilon_csf, sigma_csf

    sig_eps_fn = Property(depends_on='+modified, reinf_ratio.+modified')

    @cached_property
    def _get_sig_eps_fn(self):
        '''Get the stress and strain arrays'''
        n_points = 100
        sigma_c_arr = np.linspace(0, self.sigma_cu, n_points)
        if self.sigma_cu == self.sigma_fu * self.V_f * cos(self.orientation):
            epsilon_c_arr = self._get_epsilon_c(sigma_c_arr)
        else:
            epsilon_c_arr = sigma_c_arr / self.E_c

        # stress with respect to reinforcement
        sigma_f_arr = sigma_c_arr / self.rho

        # stress of reinforcement with no matrix interaction
        sigma_fiber = epsilon_c_arr[[0, -1]] * self.material.E_f * self.rho

        return epsilon_c_arr, sigma_c_arr, sigma_f_arr, sigma_fiber

if __name__ == '__main__':
    c = Cuypers()
    plt.plot(c.sig_eps_fn[0], c.sig_eps_fn[1], lw=2)
    plt.show()
