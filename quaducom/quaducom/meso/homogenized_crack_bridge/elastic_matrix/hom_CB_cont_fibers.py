'''
Created on Jan 21, 2015

@author: rostislavrypl
'''

from traits.api import HasTraits, List, Property, cached_property, Instance, Array, Float, Bool
from types import FloatType
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import Reinforcement
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers
from scipy.optimize import root
import time as t
from scipy.integrate import cumtrapz
import numpy as np
from spirrid.rv import RV


class CrackBridgeContFibers(HasTraits):
    
    cont_reinf_lst = List(Instance(Reinforcement))
    w = Float
    Ll = Float
    Lr = Float
    E_m = Float
    E_c = Float
    epsm_softening = Float
    w_unld = Float
    damage_switch = Bool(True)
    
    V_f_tot_cont = Property(depends_on='con_reinf_lst+')
    @cached_property
    def _get_V_f_tot_cont(self):
        V_f_tot = 0.0
        for reinf in self.cont_reinf_lst:
            V_f_tot += reinf.V_f
        return V_f_tot

    sorted_theta = Property(depends_on='cont_reinf_lst+')
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
        for reinf in self.cont_reinf_lst:
            n_int = len(np.hstack((np.array([]), reinf.depsf_arr)))
            depsf_arr = np.hstack((depsf_arr, reinf.depsf_arr))
            V_f_arr = np.hstack((V_f_arr, np.repeat(reinf.V_f, n_int)))
            E_f_arr = np.hstack((E_f_arr, np.repeat(reinf.E_f, n_int)))
            xi_arr = np.hstack((xi_arr, np.repeat(reinf.xi, n_int)))
            stat_weights_arr = np.hstack((stat_weights_arr, reinf.stat_weights))
            nu_r_arr = np.hstack((nu_r_arr, reinf.nu_r))
            r_arr = np.hstack((r_arr, reinf.r_arr))
        argsort = np.argsort(depsf_arr)[::-1]
        # sorting the masks for the evaluation of F
        idxs = np.array([])
        for i, reinf in enumerate(self.cont_reinf_lst):
            idxs = np.hstack((idxs, i * np.ones_like(reinf.depsf_arr)))
        masks = []
        for i, reinf in enumerate(self.cont_reinf_lst):
            masks.append((idxs == i)[argsort])
        max_depsf = [np.max(reinf.depsf_arr) for reinf in self.cont_reinf_lst]
        masks = [masks[i] for i in np.argsort(max_depsf)[::-1]]
        return depsf_arr[argsort], V_f_arr[argsort], E_f_arr[argsort], \
                xi_arr[argsort], stat_weights_arr[argsort], \
                nu_r_arr[argsort], masks, r_arr[argsort]

    sorted_depsf = Property(depends_on='cont_reinf_lst+')
    @cached_property
    def _get_sorted_depsf(self):
        return self.sorted_theta[0]

    sorted_V_f = Property(depends_on='cont_reinf_lst+')
    @cached_property
    def _get_sorted_V_f(self):
        return self.sorted_theta[1]

    sorted_E_f = Property(depends_on='cont_reinf_lst+')
    @cached_property
    def _get_sorted_E_f(self):
        return self.sorted_theta[2]

    sorted_xi = Property(depends_on='cont_reinf_lst+')
    @cached_property
    def _get_sorted_xi(self):
        return self.sorted_theta[3]

    sorted_stats_weights = Property(depends_on='cont_reinf_lst+')
    @cached_property
    def _get_sorted_stats_weights(self):
        return self.sorted_theta[4]

    sorted_nu_r = Property(depends_on='cont_reinf_lst+')
    @cached_property
    def _get_sorted_nu_r(self):
        return self.sorted_theta[5]

    sorted_masks = Property(depends_on='cont_reinf_lst+')
    @cached_property
    def _get_sorted_masks(self):
        return self.sorted_theta[6]

    sorted_r = Property(depends_on='cont_reinf_lst+')
    @cached_property
    def _get_sorted_r(self):
        return self.sorted_theta[7]

    sorted_xi_cdf = Property(depends_on='cont_reinf_lst+,Ll,Lr')
    @cached_property
    def _get_sorted_xi_cdf(self):
        '''breaking strain: CDF for random and Heaviside for discrete values'''
        # TODO: does not work for reinforcement types with the same xi
        methods = []
        masks = []
        for reinf in self.cont_reinf_lst:
            masks.append(self.sorted_xi == reinf.xi)
            if isinstance(reinf.xi, FloatType):
                methods.append(lambda x: 1.0 * (reinf.xi <= x))
            elif isinstance(reinf.xi, RV):
                methods.append(reinf.xi._distr.cdf)
            elif isinstance(reinf.xi, WeibullFibers):
                reinf.xi.Ll = self.Ll
                reinf.xi.Lr = self.Lr
                methods.append(reinf.xi.cdf)
        return methods, masks

    Kf = Property(depends_on='cont_reinf_lst+')
    @cached_property
    def _get_Kf(self):
        return self.sorted_V_f * self.sorted_nu_r * \
                self.sorted_stats_weights * self.sorted_E_f

    def vect_xi_cdf(self, epsy, x_short, x_long):
        Pf = np.zeros_like(self.sorted_depsf)
        methods, masks = self.sorted_xi_cdf
        for i, method in enumerate(methods):
            if method.__doc__ == 'weibull_fibers_cdf_mc':
                Pf[masks[i]] += method(epsy[masks[i]],
                                       self.sorted_depsf[masks[i]],
                                       self.sorted_r[masks[i]],
                                       x_short[masks[i]],
                                       x_long[masks[i]])
            elif method.__doc__ == 'weibull_fibers_cdf_cb_elast':
                Pf[masks[i]] += method(epsy[masks[i]],
                                       self.sorted_depsf[masks[i]],
                                       self.sorted_r[masks[i]],
                                       x_short[masks[i]],
                                       x_long[masks[i]])
            else:
                Pf[masks[i]] += method(epsy[masks[i]])
        return Pf

    def dem_depsf_vect(self, damage):
        '''evaluates the deps_m given deps_f
        at that point and the damage array'''
        Kf_intact = self.Kf * (1. - damage)
        Kf_intact_bonded = np.hstack((0.0, np.cumsum((Kf_intact))))[:-1]
        Kf_broken = np.sum(self.Kf - Kf_intact)
        Kf_add = Kf_intact_bonded + Kf_broken
        Km = self.E_c - np.sum(self.Kf)
        E_mtrx = Km + Kf_add
        mu_T = np.cumsum((self.sorted_depsf * Kf_intact)[::-1])[::-1]
        return mu_T / E_mtrx

    def F(self, dems, amin):
        '''Eq.(D.21)'''
        f = 1. / (self.sorted_depsf + dems)
        F = np.hstack((0., cumtrapz(f, -self.sorted_depsf)))
        return F

#     def F(self, dems, amin):
#         '''Auxiliary function (see Part II, appendix B)
#         '''
#         F = np.zeros_like(self.sorted_depsf)
#         for i, mask in enumerate(self.sorted_masks):
#             depsfi = self.sorted_depsf[mask]
#             demsi = dems[mask]
#             fi = 1. / (depsfi + demsi)
#             F[mask] = np.hstack((np.array([0.0]), cumtrapz(fi, -depsfi)))
#             if i == 0:
#                 C = 0.0
#             else:
#                 depsf0 = self.sorted_depsf[self.sorted_masks[i - 1]]
#                 depsf1 = depsfi[0]
#                 idx = np.sum(depsf0 > depsf1) - 1
#                 depsf2 = depsf0[idx]
#                 a1 = np.exp(F[self.sorted_masks[i - 1]][idx] / 2. + np.log(amin))
#                 p = depsf2 - depsf1
#                 q = depsf1 + demsi[0]
#                 amin_i = np.sqrt(a1 ** 2 + p / q * a1 ** 2)
#                 C = np.log(amin_i / amin)
#             F[mask] += 2 * C
#         return F

    def clamped(self, Lmin, Lmax, init_dem):
        a = np.hstack((-Lmin, 0.0, Lmax))
        em = np.hstack((init_dem * Lmin, 0.0, init_dem * Lmax))
        epsf0 = (self.sorted_depsf / 2. * (Lmin ** 2 + Lmax ** 2) +
                     self.w + em[0] * Lmin / 2. + em[-1] * Lmax / 2.) / (Lmin + Lmax)
        return a, em, epsf0

    def profile(self, iter_damage):
        '''
        Evaluates the maximum fiber strain, debonded lengths and matrix strain profile
        given the boundaries, crack opening and damage vector.
        '''
        Lmin = min(self.Ll, self.Lr)
        Lmax = max(self.Ll, self.Lr)
        # matrix strain derivative with resp. to z as a function of T
        dems = self.dem_depsf_vect(iter_damage)
        # initial matrix strain derivative
        init_dem = dems[0]
        # debonded length of fibers with Tmax
        amin = (self.w / (np.abs(init_dem) + np.abs(self.sorted_depsf[0]))) ** 0.5
        # integrated f(depsf) - see article
        F = self.F(dems, amin)
        # a1 is a(depsf) for double sided pullout
        a1 = amin * np.exp(F / 2.)
        #aX = np.exp((-np.log(np.abs(self.sorted_depsf) + dems) + np.log(self.w)) / 2.)
        if Lmin < a1[0] and Lmax < a1[0]:
            # all fibers debonded up to Lmin and Lmax
            a, em, epsf0 = self.clamped(Lmin, Lmax, init_dem)

        elif Lmin < a1[0] and Lmax >= a1[0]:
            # all fibers debonded up to Lmin but not up to Lmax
            amin = -Lmin + np.sqrt(2 * Lmin ** 2 + 2 * self.w / (self.sorted_depsf[0] + init_dem))
            C = np.log(amin ** 2 + 2 * Lmin * amin - Lmin ** 2)
            a2 = np.sqrt(2 * Lmin ** 2 + np.exp((F + C))) - Lmin
            if Lmax < a2[0]:
                a, em, epsf0 = self.clamped(Lmin, Lmax, init_dem)
            else:
                if Lmax <= a2[-1]:
                    idx = np.sum(a2 < Lmax) - 1
                    a = np.hstack((-Lmin, 0.0, a2[:idx + 1], Lmax))
                    em2 = np.cumsum(np.diff(np.hstack((0.0, a2))) * dems)
                    em = np.hstack((init_dem * Lmin, 0.0, em2[:idx + 1], em2[idx] + (Lmax - a2[idx]) * dems[idx]))
                    um = np.trapz(em, a)
                    epsf01 = em2[:idx + 1] + a2[:idx + 1] * self.sorted_depsf[:idx + 1]
                    epsf02 = (self.w + um + self.sorted_depsf[idx + 1:] / 2. * (Lmin ** 2 + Lmax ** 2)) / (Lmin + Lmax)
                    epsf0 = np.hstack((epsf01, epsf02))
                else:
                    a = np.hstack((-Lmin, 0.0, a2, Lmax))
                    em2 = np.cumsum(np.diff(np.hstack((0.0, a2))) * dems)
                    em = np.hstack((init_dem * Lmin, 0.0, em2, em2[-1]))
                    epsf0 = em2 + self.sorted_depsf * a2
        elif a1[0] < Lmin and a1[-1] > Lmin:
            # some fibers are debonded up to Lmin, some are not
            # boundary condition position
            idx1 = np.sum(a1 <= Lmin)
            # a(T) for one sided pullout
            # first debonded length amin for one sided pull out
#             depsfLmin = self.sorted_depsf[idx1]
#             p = (depsfLmin + dems[idx1])
#             a_short = np.hstack((a1[:idx1], Lmin))
#             em_short = np.cumsum(np.diff(np.hstack((0.0, a_short))) * dems[:idx1 + 1]) 
#             emLmin = em_short[-1]
#             umLmin = np.trapz(np.hstack((0.0, em_short)), np.hstack((0.0, a_short)))
#             amin = -Lmin + np.sqrt(4 * Lmin ** 2 * p ** 2 - 4 * p * emLmin * Lmin + 4 * p * umLmin - 2 * p * Lmin ** 2 * depsfLmin + 2 * p * self.w) / p
#             C = np.log(amin ** 2 + 2 * amin * Lmin - Lmin ** 2)
#             a2 = (np.sqrt(2 * Lmin ** 2 + np.exp(F + C - F[idx1])) - Lmin)[idx1:]
            a2 = np.sqrt(2 * Lmin ** 2 + np.exp(F[idx1:])*2*self.w/(self.sorted_depsf[0] + init_dem))-Lmin
            # matrix strain profiles - shorter side
            a_short = np.hstack((-Lmin, -a1[:idx1][::-1], 0.0))
            dems_short = np.hstack((dems[:idx1], dems[idx1]))
            em_short = np.hstack((0.0, np.cumsum(np.diff(-a_short[::-1]) * dems_short)))[::-1]
            if a2[-1] > Lmax:
                idx2 = np.sum(a2 <= Lmax)
                # matrix strain profiles - longer side
                a_long = np.hstack((a1[:idx1], a2[:idx2]))
                em_long = np.cumsum(np.diff(np.hstack((0.0, a_long))) * dems[:idx1 + idx2])
                a = np.hstack((a_short, a_long, Lmax))
                em = np.hstack((em_short, em_long, em_long[-1] + (Lmax - a_long[-1]) * dems[idx1 + idx2]))
                um = np.trapz(em, a)
                epsf01 = em_long + a_long * self.sorted_depsf[:idx1 + idx2]
                epsf02 = (self.w + um + self.sorted_depsf [idx1 + idx2:] / 2. * (Lmin ** 2 + Lmax ** 2)) / (Lmin + Lmax)
                epsf0 = np.hstack((epsf01, epsf02))
            else:
                a_long = np.hstack((0.0, a1[:idx1], a2, Lmax))
                a = np.hstack((a_short, a_long[1:]))
                dems_long = dems
                em_long = np.hstack((np.cumsum(np.diff(a_long[:-1]) * dems_long)))
                em_long = np.hstack((em_long, em_long[-1]))
                em = np.hstack((em_short, em_long))
                epsf0 = em_long[:-1] + self.sorted_depsf * a_long[1:-1]
        elif a1[-1] <= Lmin:
            # double sided pullout
            a = np.hstack((-Lmin, -a1[::-1], 0.0, a1, Lmax))
            em1 = np.cumsum(np.diff(np.hstack((0.0, a1))) * dems)
            em = np.hstack((em1[-1], em1[::-1], 0.0, em1, em1[-1]))
            epsf0 = em1 + self.sorted_depsf * a1
        a_short = -a[a < 0.0][1:][::-1]
        if len(a_short) < len(self.sorted_depsf):
            a_short = np.hstack((a_short, Lmin * np.ones(len(self.sorted_depsf) - len(a_short))))
        a_long = a[a > 0.0][:-1]
        if len(a_long) < len(self.sorted_depsf):
            a_long = np.hstack((a_long, Lmax * np.ones(len(self.sorted_depsf) - len(a_long))))
        return epsf0 + self.epsm_softening, a_short, a_long, em + self.epsm_softening, a

    def damage_residuum(self, iter_damage):
        if np.any(iter_damage < 0.0) or np.any(iter_damage > 1.0):
            return np.ones_like(iter_damage) * 2.0
        else:
            epsf0, x_short, x_long, epsm_arr, x_arr = self.profile(iter_damage)
            residuum = self.vect_xi_cdf(epsf0, x_short=x_short, x_long=x_long) - iter_damage
            return residuum

    x_arr = Property(Array, depends_on='reinforcement_lst,w,Ll,Lr,E_m')
    @cached_property
    def _get_x_arr(self):
        return self.profile(self.damage)[4] 

    epsm_arr = Property(Array, depends_on='reinforcement_lst,w,Ll,Lr,E_m')
    @cached_property
    def _get_epsm_arr(self):
        return self.profile(self.damage)[3]
    
    epsf_arr = Property(Array, depends_on='reinforcement_lst,w,Ll,Lr,E_m')
    @cached_property
    def _get_epsf_arr(self):
        epsf_xi = self.sorted_depsf.reshape(len(self.sorted_depsf),1) * np.abs(self.x_arr.reshape(1,len(self.x_arr)))
        epsf_x = np.maximum(self.epsf0_arr.reshape(len(self.epsf0_arr),1) - epsf_xi, self.epsm_arr.reshape(1,len(self.epsm_arr)))
        epsf_x = np.mean(epsf_x, axis=0)
        return epsf_x

    epsf0_arr = Property(Array, depends_on='reinforcement_lst,w,Ll,Lr,E_m')
    @cached_property
    def _get_epsf0_arr(self):
        return self.profile(self.damage)[0]

    maximum_damage = Array

    damage = Property(depends_on='w, Ll, Lr, reinforcement+')
    @cached_property
    def _get_damage(self):
        if len(self.maximum_damage) == 0:
                self.maximum_damage = np.zeros_like(self.sorted_depsf)
        if self.w == 0.:
            damage = np.zeros_like(self.sorted_depsf)
        else:
            ff = t.clock()
            try:
                damage = root(self.damage_residuum, self.maximum_damage,
                              method='excitingmixing', options={'maxiter':100})
                if np.any(damage.x < 0.0) or np.any(damage.x > 1.0):
                    raise ValueError
                damage = damage.x
            except:
                print 'fast opt method does not converge: switched to a slower, robust method for this step'
                damage = root(self.damage_residuum, np.ones_like(self.sorted_depsf) * 0.2,
                              method='krylov')
                damage = damage.x
            # print 'damage =', np.sum(damage) / len(damage), 'iteration time =', t.clock() - ff, 'sec'
        if self.damage_switch == False:
            return np.maximum(damage, self.maximum_damage)
        elif self.damage_switch == True:
            self.maximum_damage = np.maximum(damage, self.maximum_damage)
            return self.maximum_damage

if __name__ == '__main__':
    pass
    