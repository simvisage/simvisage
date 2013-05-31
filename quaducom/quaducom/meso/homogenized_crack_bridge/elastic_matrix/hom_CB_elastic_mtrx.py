'''
Created on Sep 20, 2012

The CompositeCrackBridge class has a method for evaluating fibers and matrix
strain in the vicinity of a crack bridge.
Fiber diameter and bond coefficient can be set as random variables.
Reinforcement types can be combined by creating a list of Reinforcement
instances and defining it as the reinforcement_lst Trait in the
CompositeCrackBridge class.
The evaluation is array based.

@author: rostar
'''
import numpy as np
from spirrid.rv import RV
from etsproxy.traits.api import HasTraits, cached_property, \
    Float, Property, Instance, List, Array
from types import FloatType
from reinforcement import Reinforcement, WeibullFibers
from scipy.optimize import fsolve, broyden2, root
import time as t
from scipy.integrate import cumtrapz


class CompositeCrackBridge(HasTraits):

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
        E_c = self.E_m * (1. - self.V_f_tot) + E_fibers
        return E_c * (1. + 1e-15)

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
        r_arr = np.array([])
        for reinf in self.reinforcement_lst:
            n_int = len(np.hstack((np.array([]), reinf.depsf_arr)))
            depsf_arr = np.hstack((depsf_arr, reinf.depsf_arr))
            V_f_arr = np.hstack((V_f_arr, np.repeat(reinf.V_f, n_int)))
            E_f_arr = np.hstack((E_f_arr, np.repeat(reinf.E_f, n_int)))
            xi_arr = np.hstack((xi_arr, np.repeat(reinf.xi, n_int)))
            stat_weights_arr = np.hstack((stat_weights_arr,
                                          np.repeat(reinf.stat_weights, n_int)))
            nu_r_arr = np.hstack((nu_r_arr, reinf.nu_r))
            r_arr = np.hstack((r_arr, reinf.r_arr))
        argsort = np.argsort(depsf_arr)[::-1]
        idxs = np.array([])
        for i, reinf in enumerate(self.reinforcement_lst):
            idxs = np.hstack((idxs, i * np.ones_like(reinf.depsf_arr)))
        masks = []
        for i, reinf in enumerate(self.reinforcement_lst):
            masks.append((idxs == i)[argsort])
        return depsf_arr[argsort], V_f_arr[argsort], E_f_arr[argsort], \
                xi_arr[argsort],  stat_weights_arr[argsort], \
                nu_r_arr[argsort], masks, r_arr[argsort]

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

    sorted_masks = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_sorted_masks(self):
        return self.sorted_theta[6]

    sorted_r = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_sorted_r(self):
        return self.sorted_theta[7]

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
                             x_short, x_long, self.sorted_r)
            else:
                Pf += method(epsy * masks[i])
        return Pf

    def dem_depsf_vect(self, depsf, damage):
        '''evaluates the deps_m given deps_f
        at that point and the damage array'''
        Kf = self.sorted_V_f * self.sorted_nu_r * \
            self.sorted_stats_weights * self.sorted_E_f
        Kf_intact_bonded = np.hstack((0.0, np.cumsum((Kf * (1. - damage)))))[:-1]
        Kf_broken = np.sum(Kf * damage)
        Kf_add = Kf_intact_bonded + Kf_broken
        Km = (1. - self.V_f_tot) * self.E_m
        E_mtrx = Km + Kf_add
        mu_T = np.cumsum((self.sorted_depsf * Kf * (1. - damage))[::-1])[::-1]
        return mu_T / E_mtrx

    def F(self, dems, amin):
        F = np.zeros_like(self.sorted_depsf)
        for i, mask in enumerate(self.sorted_masks):
            depsfi = self.sorted_depsf[mask]
            demsi = dems[mask]
            fi = 1. / (depsfi + demsi)
            F[mask] = np.hstack((np.array([0.0]), cumtrapz(fi, -depsfi)))
            if i==0:
                C = 0.0
            else:
                depsf0 = self.sorted_depsf[self.sorted_masks[i-1]]
                depsf1 = depsfi[0]
                idx = np.sum(depsf0 > depsf1) - 1
                depsf2 = depsf0[idx]
                a1 = np.exp(F[self.sorted_masks[i-1]][idx]/2. + np.log(amin))
                p = depsf2 - depsf1
                q = depsf1 + demsi[0]
                amin_i = np.sqrt(a1**2 + p/q*a1**2)
                C = np.log(amin_i/amin)
            F[mask] += 2 * C
        #plt.plot(-self.sorted_depsf, F)
        #plt.show()
        return F

    def profile(self, iter_damage, Lmin, Lmax):
        if np.any(iter_damage < 0.0) or np.any(iter_damage > 1.0):
            print 'heereee'
            return np.ones_like(iter_damage) * 0.5, np.ones_like(self.sorted_depsf), np.ones_like(self.sorted_depsf)
        # matrix strain derivative with resp. to z as a function of T
        dems = self.dem_depsf_vect(self.sorted_depsf, iter_damage)
        # initial matrix strain derivative
        init_dem = dems[0]
        # debonded length of fibers with Tmax
        amin = (self.w / (np.abs(init_dem) + np.abs(self.sorted_depsf[0])))**0.5
        # integrated f(depsf) - see article
        F = self.F(dems, amin)
        # a(T) for double sided pullout
        a1 = np.exp(F/2. + np.log(amin))
        if Lmin < a1[0] and Lmax < a1[0]:
            # all fibers debonded up to Lmin and Lmax
            a = np.hstack((-Lmin, 0.0, Lmax))
            em = np.hstack((init_dem * Lmin, 0.0, init_dem * Lmax))
            epsf0 = (self.sorted_depsf/2. * (Lmin**2 + Lmax**2) +
                     self.w + em[0] * Lmin / 2. + em[-1] * Lmax / 2.) / (Lmin + Lmax)
        elif Lmin < a1[0] and Lmax >= a1[0]:
            # all fibers debonded up to Lmin but not up to Lmax
            amin = -Lmin + np.sqrt(2 * Lmin**2 + 2*self.w / (self.sorted_depsf[0] + init_dem))
            C = np.log(amin**2 + 2*Lmin*amin - Lmin**2)
            a2 = np.sqrt(2*Lmin**2 + np.exp((F + C))) - Lmin
            if Lmax <= a2[-1]:
                idx = np.sum(a2 < Lmax) - 1
                a = np.hstack((-Lmin, 0.0, a2[:idx + 1], Lmax))
                em2 = np.cumsum(np.diff(np.hstack((0.0, a2)))*dems)
                em = np.hstack((init_dem*Lmin, 0.0, em2[:idx+1], em2[idx] + (Lmax-a2[idx])*dems[idx]))
                um = np.trapz(em, a)
                epsf01 = em2[:idx + 1] + a2[:idx + 1] * self.sorted_depsf[:idx+1]
                epsf02 = (self.w + um + self.sorted_depsf[idx + 1:] / 2. * (Lmin**2 + Lmax**2)) / (Lmin +Lmax)
                epsf0 = np.hstack((epsf01, epsf02))
            else:
                a = np.hstack((-Lmin, 0.0, a2, Lmax))
                em2 = np.cumsum(np.diff(np.hstack((0.0, a2)))*dems)
                em = np.hstack((init_dem * Lmin, 0.0, em2, em2[-1]))
                epsf0 = em2 + self.sorted_depsf * a2
        elif a1[0] < Lmin and a1[-1] > Lmin:
            # boundary condition position
            idx1 = np.sum(a1 <= Lmin)
            # a(T) for one sided pullout
            # first debonded length amin for one sided PO
            depsfLmin = self.sorted_depsf[idx1]
            p = (depsfLmin + dems[idx1])
            a_short = np.hstack((a1[:idx1], Lmin))
            em_short = np.cumsum(np.diff(np.hstack((0.0, a_short)))*dems[:idx1+1])
            emLmin = em_short[-1]
            umLmin = np.trapz(np.hstack((0.0, em_short)), np.hstack((0.0, a_short)))
            amin = - Lmin + np.sqrt(4 * Lmin**2 * p**2 - 4 * p * emLmin * Lmin + 4 * p * umLmin - 2 * p * Lmin**2 * depsfLmin + 2 * p * self.w)/p
            C = np.log(amin**2 + 2 * amin * Lmin - Lmin**2)
            a2 = (np.sqrt(2*Lmin**2 + np.exp(F + C - F[idx1])) - Lmin)[idx1:]
            # matrix strain profiles - shorter side
            a_short = np.hstack((-Lmin, -a1[:idx1][::-1], 0.0))
            dems_short = np.hstack((dems[:idx1], dems[idx1]))
            em_short = np.hstack((0.0,np.cumsum(np.diff(-a_short[::-1])*dems_short)))[::-1]
            if a2[-1] > Lmax:
                idx2 = np.sum(a2 <= Lmax)
                # matrix strain profiles - longer side
                a_long = np.hstack((a1[:idx1], a2[:idx2]))
                em_long = np.cumsum(np.diff(np.hstack((0.0, a_long)))*dems[:idx1+idx2])
                a = np.hstack((a_short, a_long, Lmax))
                em = np.hstack((em_short, em_long, em_long[-1] + (Lmax-a_long[-1])*dems[idx1+idx2]))
                um = np.trapz(em, a)
                epsf01 = em_long + a_long * self.sorted_depsf[:idx1+idx2]
                epsf02 = (self.w + um + self.sorted_depsf [idx1+idx2:] / 2. * (Lmin**2 + Lmax**2)) / (Lmin +Lmax)
                epsf0 = np.hstack((epsf01, epsf02))
            else:
                a_long = np.hstack((0.0, a1[:idx1], a2, Lmax))
                a = np.hstack((a_short, a_long[1:]))
                dems_long = dems
                em_long = np.hstack((np.cumsum(np.diff(a_long[:-1])*dems_long)))
                em_long = np.hstack((em_long, em_long[-1]))
                em = np.hstack((em_short, em_long))
                epsf0 = em_long[:-1] + self.sorted_depsf * a_long[1:-1]
        elif a1[-1] <= Lmin:
            #double sided pullout
            a = np.hstack((-Lmin, -a1[::-1], 0.0, a1, Lmax))
            em1 = np.cumsum(np.diff(np.hstack((0.0, a1)))*dems)
            em = np.hstack((em1[-1], em1[::-1], 0.0, em1, em1[-1]))
            epsf0 = em1 + self.sorted_depsf * a1
        self._x_arr = a
        self._epsm_arr = em
        self._epsf0_arr = epsf0
        a_short = -a[a<0.0][1:][::-1]
        if len(a_short) < len(self.sorted_depsf):
            a_short = np.hstack((a_short, Lmin * np.ones(len(self.sorted_depsf) - len(a_short))))
        a_long = a[a>0.0][:-1]
        if len(a_long) < len(self.sorted_depsf):
            a_long = np.hstack((a_long, Lmax * np.ones(len(self.sorted_depsf) - len(a_long))))
        return epsf0, a_short, a_long

    def damage_residuum(self, iter_damage):
        Lmin = min(self.Ll, self.Lr)
        Lmax = max(self.Ll, self.Lr)
        epsf0, x_short, x_long = self.profile(iter_damage, Lmin, Lmax)
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
        if self.w == 0.:
            damage = np.zeros_like(self.sorted_depsf)
        else:
            ff = t.clock()
            try:
                damage = root(self.damage_residuum, np.zeros_like(self.sorted_depsf), method='excitingmixing')
                if np.any(damage.x < 0.0) or np.any(damage.x > 1.0):
                    raise ValueError
                damage = damage.x
            except:
                print 'fast opt method does not converge: switched to a slower, robust method for this step'
                damage = root(self.damage_residuum, np.zeros_like(self.sorted_depsf), method='krylov')
                damage = damage.x
            #print 'damage =', np.sum(damage) / len(damage), 'iteration time =', t.clock() - ff, 'sec'
        return damage

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    reinf1 = Reinforcement(r=0.00345,#RV('uniform', loc=0.001, scale=0.005),
                          tau=RV('uniform', loc=1., scale=20.),
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

    ccb = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf1],
                                 Ll=.7,
                                 Lr=1.,
                                 w=0.028)

    reinf = Reinforcement(r=0.01,
                          tau=RV('uniform', loc=0.01, scale=.5),
                          V_f=0.05,
                          E_f=200e3,
                          xi=WeibullFibers(shape=5., sV0=0.00618983207723),
                          n_int=50,
                          label='carbon')

    ccb = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf],
                                 Ll=4.,
                                 Lr=87.,
                                 w=0.004)

    ccb.damage
    plt.plot(ccb._x_arr, ccb._epsm_arr, lw=2, color='red', ls='dashed', label='analytical')
    plt.plot(np.zeros_like(ccb._epsf0_arr), ccb._epsf0_arr, 'ro')
    for i, depsf in enumerate(ccb.sorted_depsf):
        plt.plot(ccb._x_arr, np.maximum(ccb._epsf0_arr[i] - depsf*np.abs(ccb._x_arr),ccb._epsm_arr))
    plt.legend(loc='best')
    plt.xlim(-5,109)
    plt.show()
