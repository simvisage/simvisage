'''
Created on Sep 20, 2012

The CompositeCrackBridgeLoop class has a method for evaluating fibers and matrix
strain in the vicinity of a crack bridge.
Fiber diameter and bond coefficient can be set as random variables.
Reinforcement types can be combined by creating a list of Reinforcement
instances and defining it as the reinforcement_lst Trait in the
CompositeCrackBridgeLoop class.
The evaluation uses a python loop over the discretization point.

@author: rostar
'''
import numpy as np
from spirrid.rv import RV
from etsproxy.traits.api import HasTraits, cached_property, \
    Float, Property, Instance, List, Array
from types import FloatType
from reinforcement import Reinforcement, WeibullFibers
from scipy.optimize import fsolve, broyden2
import time as t
import time

class CompositeCrackBridgeLoop(HasTraits):

    reinforcement_lst = List(Instance(Reinforcement))
    w = Float
    E_m = Float
    Ll = Float
    Lr = Float

    V_f_tot = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_V_f_tot(self):
        V_f_tot = 0.0
        for reinf in self.reinforcement_lst:
            V_f_tot += reinf.V_f
        return V_f_tot

    E_c = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_E_c(self):
        E_fibers = 0.0
        for reinf in self.reinforcement_lst:
            E_fibers += reinf.V_f * reinf.E_f
        return self.E_m * (1. - self.V_f_tot) + E_fibers

    sorted_theta = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_sorted_theta(self):
        '''sorts the integral points by bond in descending order'''
        depsf_arr = np.array([])
        V_f_arr = np.array([])
        E_f_arr = np.array([])
        xi_arr = np.array([])
        stat_weights_arr = np.array([])
        nu_r_arr = np.array([])
        for reinf in self.reinforcement_lst:
            n_int = len(np.hstack((np.array([]), reinf.depsf_arr)))
            depsf_arr = np.hstack((depsf_arr, reinf.depsf_arr))
            V_f_arr = np.hstack((V_f_arr, np.repeat(reinf.V_f, n_int)))
            E_f_arr = np.hstack((E_f_arr, np.repeat(reinf.E_f, n_int)))
            xi_arr = np.hstack((xi_arr, np.repeat(reinf.xi, n_int)))
            stat_weights_arr = np.hstack((stat_weights_arr,
                                          np.repeat(reinf.stat_weights, n_int)))
            nu_r_arr = np.hstack((nu_r_arr, reinf.nu_r))
        argsort = np.argsort(depsf_arr)[::-1]
        return depsf_arr[argsort], V_f_arr[argsort], E_f_arr[argsort], \
                xi_arr[argsort],  stat_weights_arr[argsort], \
                nu_r_arr[argsort]

    sorted_depsf = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_sorted_depsf(self):
        return self.sorted_theta[0]

    sorted_V_f = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_sorted_V_f(self):
        return self.sorted_theta[1]

    sorted_E_f = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_sorted_E_f(self):
        return self.sorted_theta[2]

    sorted_xi = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_sorted_xi(self):
        return self.sorted_theta[3]

    sorted_stats_weights = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_sorted_stats_weights(self):
        return self.sorted_theta[4]

    sorted_nu_r = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_sorted_nu_r(self):
        return self.sorted_theta[5]

    sorted_xi_cdf = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_sorted_xi_cdf(self):
        '''breaking strain: CDF for random and Heaviside for discrete values'''
        # TODO: does not work for reinforcement types with the same xi
        methods = []
        masks = []
        for reinf in self.reinforcement_lst:
            masks.append(self.sorted_xi == reinf.xi)
            if isinstance(reinf.xi, FloatType):
                methods.append(lambda x: 1.0 * (reinf.xi <= x))
            elif isinstance(reinf.xi, RV):
                methods.append(reinf.xi._distr.cdf)
            elif isinstance(reinf.xi, WeibullFibers):
                methods.append(reinf.xi.weibull_fibers_Pf)
        return methods, masks

    def vect_xi_cdf(self, epsy, x_short, x_long):
        Pf = np.zeros_like(self.sorted_depsf)
        methods, masks = self.sorted_xi_cdf
        for i, method in enumerate(methods):
            if method.__name__ == 'weibull_fibers_Pf':
                Pf += method(epsy * masks[i], self.sorted_depsf,
                             x_short=x_short, x_long=x_long)
            else:
                Pf += method(epsy * masks[i])
        return Pf

    def dem_depsf(self, depsf, damage):
        '''evaluates the deps_m given deps_f
        at that point and the damage array'''
        Kf = self.sorted_V_f * self.sorted_nu_r * \
            self.sorted_stats_weights * self.sorted_E_f
        Kf_intact_bonded = np.sum(Kf * (depsf <= self.sorted_depsf)
                                         * (1. - damage))
        Kf_broken = np.sum(Kf * damage)
        Kf_add = Kf_intact_bonded + Kf_broken
        Km = (1. - self.V_f_tot) * self.E_m
        E_mtrx = Km + Kf_add
        mean_acting_T = np.sum(self.sorted_depsf * (self.sorted_depsf < depsf) *
                                   Kf * (1. - damage))
        return mean_acting_T / E_mtrx

    def double_sided(self, defi, x0, demi, em0, um0, damage):
        dxi = (-defi * x0 - demi * x0 + (defi * x0 ** 2 * demi
            + demi ** 2 * x0 ** 2 - 2 * defi * em0 * x0 + 2 *
            defi * um0 + defi * self.w - 2 * demi * em0 * x0 +
            2 * demi * um0 + demi * self.w) ** (.5)) / (defi + demi)
        dem = self.dem_depsf(defi, damage)
        emi = em0 + demi * dxi
        umi = um0 + (em0 + emi) * dxi / 2.
        return dxi, dem, emi, umi

    def one_sided(self, defi, x0, demi, em0, um0, clamped, damage):
        w = self.w
        xs = clamped[0]
        ums = clamped[1]
        dxi = (-xs * demi - demi * x0 - defi * xs - defi * x0 + (2 *
                demi * x0 * defi * xs + demi * x0 ** 2 * defi + 2 *
                demi ** 2 * x0 * xs + 3 * defi * xs ** 2 * demi - 2 *
                demi * xs * em0 - 2 * demi * em0 * x0 - 2 * defi *
                xs * em0 - 2 * defi * em0 * x0 + demi ** 2 * x0 ** 2 +
                2 * defi ** 2 * xs ** 2 + xs ** 2 * demi ** 2 + 2 *
                demi * um0 + 2 * demi * ums + 2 * demi * w + 2 * defi *
                um0 + 2 * defi * ums + 2 * defi * w) ** (0.5)) / (demi + defi)
        dem = self.dem_depsf(defi, damage)
        emi = em0 + demi * dxi
        umi = um0 + (em0 + emi) * dxi / 2.
        return dxi, dem, emi, umi

    def clamped(self, defi, xs, xl, ems, eml, ums, uml):
        c1 = eml * xl - uml
        c2 = ems * xs - ums
        c3 = defi * xl ** 2 / 2.
        c4 = defi * xs ** 2 / 2.
        c5 = (defi * (xl - xs) + (eml - ems)) * xs
        h = (self.w - c1 - c2 - c3 - c4 - c5) / (xl + xs)
        return defi * xl + eml + h

    def damage_residuum(self, iter_damage):
        um_short, em_short, x_short = [0.0], [0.0], [0.0]
        um_long, em_long, x_long = [0.0], [0.0], [0.0]
        init_dem = self.dem_depsf(np.infty, iter_damage)
        dem_short = [init_dem]
        dem_long = [init_dem]
        epsf0 = np.zeros_like(self.sorted_depsf)
        Lmin = min(self.Ll, self.Lr)
        Lmax = max(self.Ll, self.Lr)
        for i, defi in enumerate(self.sorted_depsf):
            if x_short[-1] < Lmin and x_long[-1] < Lmax:
                '''double sided pullout'''
                dxi, dem, emi, umi = self.double_sided(defi,
                                    x_short[-1], dem_short[-1],
                                    em_short[-1], um_short[-1], iter_damage)
                if x_short[-1] + dxi < Lmin:
                    # dx increment does not reach the boundary
                    dem_short.append(dem)
                    dem_long.append(dem)
                    x_short.append(x_short[-1] + dxi)
                    x_long.append(x_long[-1] + dxi)
                    em_short.append(emi)
                    em_long.append(emi)
                    um_short.append(umi)
                    um_long.append(umi)
                    epsf0[i] = (em_short[-1] + x_short[-1] * defi)
                else:
                    # boundary reached at shorter side
                    deltax = Lmin - x_short[-1]
                    x_short.append(Lmin)
                    em_short.append(em_short[-1] + dem_short[-1] * deltax)
                    um_short.append(um_short[-1] + (em_short[-2] + em_short[-1]) * deltax / 2.)
                    short_side = [x_short[-1], um_short[-1]]
                    dxi, dem, emi, umi = self.one_sided(defi, x_long[-1], dem_long[-1],
                                            em_long[-1], um_long[-1], short_side, iter_damage)

                    if x_long[-1] + dxi >= Lmax:
                        # boundary reached at longer side
                        deltax = Lmax - x_long[-1]
                        x_long.append(Lmax)
                        em_long.append(em_long[-1] + dem_long[-1] * deltax)
                        um_long.append(um_long[-1] + (em_long[-2] + em_long[-1]) * deltax / 2.)
                        epsf0_clamped = self.clamped(defi, x_short[-1], x_long[-1], em_short[-1],
                             em_long[-1], um_short[-1], um_long[-1])
                        epsf0[i] = epsf0_clamped
                    else:
                        dem_long.append(dem)
                        x_long.append(x_long[-1] + dxi)
                        em_long.append(emi)
                        um_long.append(umi)
                        epsf0[i] = (em_long[-1] + x_long[-1] * defi)

            elif x_short[-1] == Lmin and x_long[-1] < Lmax:
                #one sided pullout
                clamped = [x_short[-1], um_short[-1]]
                dxi, dem, emi, umi = self.one_sided(defi, x_long[-1], dem_long[-1],
                                    em_long[-1], um_long[-1], clamped, iter_damage)
                if x_long[-1] + dxi < Lmax:
                    dem_long.append(dem)
                    x_long.append(x_long[-1] + dxi)
                    em_long.append(emi)
                    um_long.append(umi)
                    epsf0[i] = (em_long[-1] + x_long[-1] * defi)
                else:
                    dxi = Lmax - x_long[-1]
                    x_long.append(Lmax)
                    em_long.append(em_long[-1] + dem_long[-1] * dxi)
                    um_long.append(um_long[-1] + (em_long[-2] + em_long[-1]) * dxi / 2.)
                    epsf0_clamped = self.clamped(defi, x_short[-1], x_long[-1], em_short[-1],
                                 em_long[-1], um_short[-1], um_long[-1])
                    epsf0[i] = epsf0_clamped

            elif x_short[-1] == Lmin and x_long[-1] == Lmax:
                #clamped fibers
                epsf0_clamped = self.clamped(defi, x_short[-1], x_long[-1], em_short[-1],
                             em_long[-1], um_short[-1], um_long[-1])
                epsf0[i] = epsf0_clamped
        self._x_arr = np.hstack((-np.array(x_short)[::-1][:-1], np.array(x_long)))
        self._epsm_arr = np.hstack((np.array(em_short)[::-1][:-1], np.array(em_long)))
        self._epsf0_arr = epsf0
        residuum = self.vect_xi_cdf(epsf0, x_short=x_short, x_long=x_long) - iter_damage
        return residuum

    _x_arr = Array
    def __x_arr_default(self):
        return np.repeat(1e-10, len(self.sorted_depsf))

    _epsm_arr = Array
    def __epsm_arr_default(self):
        return np.repeat(1e-10, len(self.sorted_depsf))

    _epsf0_arr = Array
    def __epsf0_arr_default(self):
        return np.repeat(1e-10, len(self.sorted_depsf))

    damage = Property(depends_on='w, Ll, Lr, reinforcement+')
    @cached_property
    def _get_damage(self):
        ff = time.clock()
        if self.w == 0.:
            damage = np.zeros_like(self.sorted_depsf)
        else:
            ff = t.clock()
            try:
                damage = broyden2(self.damage_residuum, 0.2 * np.ones_like(self.sorted_depsf), maxiter=20)
            except:
                print 'broyden2 does not converge fast enough: switched to fsolve for this step'
                damage = fsolve(self.damage_residuum, 0.2 * np.ones_like(self.sorted_depsf))
            print 'damage =', np.sum(damage) / len(damage), 'iteration time =', time.clock() - ff, 'sec'
        return damage

if __name__ == '__main__':
    from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
    from matplotlib import pyplot as plt

    reinf1 = Reinforcement(r=0.00345,#RV('uniform', loc=0.001, scale=0.005),
                          tau=RV('uniform', loc=4., scale=2.),
                          V_f=0.2,
                          E_f=70e3,
                          xi=RV('weibull_min', shape=5., scale=0.04),
                          n_int=100,
                          label='AR glass')

    reinf2 = Reinforcement(r=0.003,#RV('uniform', loc=0.002, scale=0.002),
                          tau=RV('uniform', loc=.3, scale=.05),
                          V_f=0.1,
                          E_f=200e3,
                          xi=WeibullFibers(shape=5., scale=0.02),
                          n_int=100,
                          label='carbon')

    ccb = CompositeCrackBridgeLoop(E_m=25e3,
                                 reinforcement_lst=[reinf1, reinf2],
                                 Ll=5.,
                                 Lr=8.,
                                 w=0.028)

    ccb.damage
    plt.plot(ccb._x_arr, ccb._epsm_arr, label='loop')
    plt.plot(np.zeros_like(ccb._epsf0_arr), ccb._epsf0_arr, 'ro')
    for i, depsf in enumerate(ccb.sorted_depsf):
        plt.plot(ccb._x_arr, np.maximum(ccb._epsf0_arr[i] - depsf*np.abs(ccb._x_arr),ccb._epsm_arr))
    plt.legend(loc='best')
    plt.xlim(-5,5)
    plt.show()
