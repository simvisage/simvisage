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
from scipy.integrate import cumtrapz
from scipy.interpolate.interpolate import interp2d
from scipy.optimize import fminbound, brentq
from scipy.optimize import root
from traits.api import HasTraits, cached_property, \
    Float, Property, Instance, List, Array, Tuple, Bool, provides
from traits.has_traits import on_trait_change
from traitsui.api import ModelView

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
import numpy as np
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import \
    Reinforcement, ContinuousFibers, ShortFibers
from spirrid import SPIRRID
from spirrid.i_rf import IRF
from spirrid.rf import RF
from spirrid.rv import RV
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers
import time as t


class CrackBridgeContFibersGf(HasTraits):

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
            stat_weights_arr = np.hstack(
                (stat_weights_arr, reinf.stat_weights))
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
            if isinstance(reinf.xi, float):
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
        amin = (self.w / (np.abs(init_dem) +
                          np.abs(self.sorted_depsf[0]))) ** 0.5
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
            amin = -Lmin + np.sqrt(2 * Lmin ** 2 + 2 *
                                   self.w / (self.sorted_depsf[0] + init_dem))
            C = np.log(amin ** 2 + 2 * Lmin * amin - Lmin ** 2)
            a2 = np.sqrt(2 * Lmin ** 2 + np.exp((F + C))) - Lmin
            if Lmax < a2[0]:
                a, em, epsf0 = self.clamped(Lmin, Lmax, init_dem)
            else:
                if Lmax <= a2[-1]:
                    idx = np.sum(a2 < Lmax) - 1
                    a = np.hstack((-Lmin, 0.0, a2[:idx + 1], Lmax))
                    em2 = np.cumsum(np.diff(np.hstack((0.0, a2))) * dems)
                    em = np.hstack(
                        (init_dem * Lmin, 0.0, em2[:idx + 1], em2[idx] + (Lmax - a2[idx]) * dems[idx]))
                    um = np.trapz(em, a)
                    epsf01 = em2[:idx + 1] + a2[:idx + 1] * \
                        self.sorted_depsf[:idx + 1]
                    epsf02 = (
                        self.w + um + self.sorted_depsf[idx + 1:] / 2. * (Lmin ** 2 + Lmax ** 2)) / (Lmin + Lmax)
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
            a2 = np.sqrt(
                2 * Lmin ** 2 + np.exp(F[idx1:]) * 2 * self.w / (self.sorted_depsf[0] + init_dem)) - Lmin
            # matrix strain profiles - shorter side
            a_short = np.hstack((-Lmin, -a1[:idx1][::-1], 0.0))
            dems_short = np.hstack((dems[:idx1], dems[idx1]))
            em_short = np.hstack(
                (0.0, np.cumsum(np.diff(-a_short[::-1]) * dems_short)))[::-1]
            if a2[-1] > Lmax:
                idx2 = np.sum(a2 <= Lmax)
                # matrix strain profiles - longer side
                a_long = np.hstack((a1[:idx1], a2[:idx2]))
                em_long = np.cumsum(
                    np.diff(np.hstack((0.0, a_long))) * dems[:idx1 + idx2])
                a = np.hstack((a_short, a_long, Lmax))
                em = np.hstack(
                    (em_short, em_long, em_long[-1] + (Lmax - a_long[-1]) * dems[idx1 + idx2]))
                um = np.trapz(em, a)
                epsf01 = em_long + a_long * self.sorted_depsf[:idx1 + idx2]
                epsf02 = (self.w + um + self.sorted_depsf[idx1 + idx2:] / 2. * (
                    Lmin ** 2 + Lmax ** 2)) / (Lmin + Lmax)
                epsf0 = np.hstack((epsf01, epsf02))
            else:
                a_long = np.hstack((0.0, a1[:idx1], a2, Lmax))
                a = np.hstack((a_short, a_long[1:]))
                dems_long = dems
                em_long = np.hstack(
                    (np.cumsum(np.diff(a_long[:-1]) * dems_long)))
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
            a_short = np.hstack(
                (a_short, Lmin * np.ones(len(self.sorted_depsf) - len(a_short))))
        a_long = a[a > 0.0][:-1]
        if len(a_long) < len(self.sorted_depsf):
            a_long = np.hstack(
                (a_long, Lmax * np.ones(len(self.sorted_depsf) - len(a_long))))
        return epsf0 + self.epsm_softening, a_short, a_long, em + self.epsm_softening, a

    def damage_residuum(self, iter_damage):
        if np.any(iter_damage < 0.0) or np.any(iter_damage > 1.0):
            return np.ones_like(iter_damage) * 2.0
        else:
            epsf0, x_short, x_long, epsm_arr, x_arr = self.profile(iter_damage)
            residuum = self.vect_xi_cdf(
                epsf0, x_short=x_short, x_long=x_long) - iter_damage
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
        epsf_xi = self.sorted_depsf.reshape(
            len(self.sorted_depsf), 1) * np.abs(self.x_arr.reshape(1, len(self.x_arr)))
        epsf_x = np.maximum(self.epsf0_arr.reshape(
            len(self.epsf0_arr), 1) - epsf_xi, self.epsm_arr.reshape(1, len(self.epsm_arr)))
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
                              method='excitingmixing', options={'maxiter': 100})
                if np.any(damage.x < 0.0) or np.any(damage.x > 1.0):
                    raise ValueError
                damage = damage.x
            except:
                print(
                    'fast opt method does not converge: switched to a slower, robust method for this step')
                damage = root(self.damage_residuum, np.ones_like(self.sorted_depsf) * 0.2,
                              method='krylov')
                damage = damage.x
            # print 'damage =', np.sum(damage) / len(damage), 'iteration time
            # =', t.clock() - ff, 'sec'
        if self.damage_switch == False:
            return np.maximum(damage, self.maximum_damage)
        elif self.damage_switch == True:
            self.maximum_damage = np.maximum(damage, self.maximum_damage)
            return self.maximum_damage


@provides(IRF)
class CBShortFiber(RF):
    '''
    Micromechanical response of a short fiber bridging a crack
    '''
    xi = Float(distr=['weibull_min', 'uniform'])
    E_f = Float(distr=['uniform', 'norm'])
    r = Float(distr=['uniform', 'norm'])
    le = Float(distr=['uniform'])
    tau = Float(distr=['norm', 'uniform', 'weibull_min'])
    snub = Float(distr=['uniform', 'norm'])
    phi = Float(distr=['sin2x', 'uniform'])
    w = Float
    C_code = ''

    def __call__(self, w, tau, r, E_f, le, phi, snub, xi, epsm_softening):
        T = 2. * tau / r
        # debonding stage
        ef0_deb = np.sqrt(T * w / E_f)
        # crack opening at which debonding is finished
        w0 = le ** 2 * T / E_f
        # pulling out stage - the fiber is pulled out from the
        # side with the shorter embedded length only
        ef0_pull = (le + w0 - w) * T / E_f
        ef0 = (ef0_deb * (w < w0) + ef0_pull * (w > w0)) * \
            np.exp(phi * snub) + epsm_softening
        # include breaking strain
        ef0 = ef0 * (ef0 < xi) * (ef0 > 0.0)
        return ef0


class CBShortFiberSP(CBShortFiber):
    '''
    stress profile for a crack bridged by a short fiber
    '''
    x = Float(distr=['uniform'])
    C_code = ''

    def __call__(self, w, x, tau, r, E_f, le, phi, snub, xi, epsm_softening):
        epsf0 = super(CBShortFiberSP, self).__call__(
            w, tau, r, E_f, le, phi, snub, xi, epsm_softening)
        T = 2. * tau / r
        epsf_x = epsf0 / np.exp(snub * phi) - np.abs(x) * T / E_f
        epsf_x = epsf_x * (epsf_x > 0.0) * np.exp(snub * phi)
        return epsf_x * E_f


class CrackBridgeShortFibersGf(HasTraits):

    short_reinf_lst = List(Instance(Reinforcement))
    w = Float
    E_c = Float
    E_m = Float
    epsm_softening = Float
    damage_switch = Bool(True)

    sorted_V_f = Property(depends_on='short_reinf_lst+')

    @cached_property
    def _get_sorted_V_f(self):
        return np.array([reinf.V_f for reinf in self.short_reinf_lst])

    sorted_E_f = Property(depends_on='short_reinf_lst+')

    @cached_property
    def _get_sorted_E_f(self):
        return np.array([reinf.E_f for reinf in self.short_reinf_lst])

    x_arr = Property(Array, depends_on='short_reinf_lst+')

    @cached_property
    def _get_x_arr(self):
        Lf_lst = []
        for reinf in self.short_reinf_lst:
            Lf_lst.append(reinf.Lf)
        max_Lf = np.max(np.array(Lf_lst))
        # !!! an even number has to be set as step for the zero position to be in the linspace !!!
        x_arr = np.linspace(-max_Lf / 2., max_Lf / 2., 61)
        return x_arr

    spirrid_lst = Property(List(Instance(SPIRRID)),
                           depends_on='short_reinf_lst+')

    @cached_property
    def _get_spirrid_lst(self):
        spirrid_epsm_list = []
        for reinf in self.short_reinf_lst:
            cb = CBShortFiberSP()
            spirrid = SPIRRID(q=cb,
                              sampling_type='LHS',
                              theta_vars=dict(epsm_softening=self.epsm_softening,
                                              tau=reinf.tau,
                                              E_f=reinf.E_f,
                                              r=reinf.r,
                                              xi=reinf.xi,
                                              snub=reinf.snub,
                                              le=reinf.le,
                                              phi=reinf.phi),
                              n_int=reinf.n_int)
            spirrid_epsm_list.append(spirrid)
        return spirrid_epsm_list

    spirrid_evaluation_cached = Property(Array, depends_on='short_reinf_lst+')

    @cached_property
    def _get_spirrid_evaluation_cached(self):
        interpolators_lst = []
        for i, spirr in enumerate(self.spirrid_lst):
            Lfi = self.short_reinf_lst[i].Lf

            def minfunc_short_fibers(w):
                spirr.eps_vars = dict(w=np.array([w]),
                                      x=np.array([0.0]))
                return -spirr.mu_q_arr.flatten()
            w_maxi = fminbound(minfunc_short_fibers, 0.0,
                               Lfi / 3., maxfun=20, disp=0)
            w_arri = np.hstack((np.linspace(0.0, w_maxi, 15),
                                np.linspace(w_maxi + 1e-10, Lfi / 2., 15)))
            spirr.eps_vars = dict(w=w_arri,
                                  x=self.x_arr)
            interpolators_lst.append(
                interp2d(self.x_arr, w_arri, spirr.mu_q_arr, fill_value=0.0))
        return interpolators_lst

    epsm_arr = Property(
        Array, depends_on='short_reinf_lst+,w,E_m,epsm_softening')

    @cached_property
    def _get_epsm_arr(self):
        epsm_x_arr = np.zeros(len(self.x_arr))
        for i, interpolator in enumerate(self.spirrid_evaluation_cached):
            sigf_x_i = interpolator(self.x_arr, self.w)
            Ff_x_i = sigf_x_i * self.sorted_V_f[i]
            Fmax_i = np.max(Ff_x_i)
            epsm_x_i = (Fmax_i - Ff_x_i) / self.E_c
            epsm_x_arr += epsm_x_i.flatten() + self.epsm_softening
        return epsm_x_arr

    epsf0_arr = Property(Array, depends_on='short_reinf_lst+,w')

    @cached_property
    def _get_epsf0_arr(self):
        return np.array([interpolator(0., self.w) / self.sorted_E_f[i] for i, interpolator
                         in enumerate(self.spirrid_evaluation_cached)]).flatten()


class CompositeCrackBridgeGf(HasTraits):

    reinforcement_lst = List(Instance(Reinforcement))
    w = Float
    E_m = Float
    Ll = Float
    Lr = Float
    ft = Float
    Gf = Float(1.0)
    w_unld = Float(0.0)
    damage_switch = Bool(True)

    @on_trait_change('damage_switch')
    def switch_damage(self):
        '''freezes the loading history'''
        self.cont_fibers_instance.damage_switch = self.damage_switch
        self.short_fibers_instance.damage_switch = self.damage_switch

    epsm_softening = Property(depends_on='w,ft,Gf,E_m')

    @cached_property
    def _get_epsm_softening(self):
        if self.w >= self.w_unld:
            if self.damage_switch == True:
                self.w_unld += self.w - self.w_unld
            return self.ft * np.exp(-self.ft / self.Gf * self.w) / self.E_m
        else:
            return self.w / self.w_unld * self.ft * np.exp(-self.ft / self.Gf * self.w_unld) / self.E_m

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

    sorted_reinf_lst = Property(
        Tuple(List, List), depends_on='reinforcement_lst')

    @cached_property
    def _get_sorted_reinf_lst(self):
        cont_reinf_lst = []
        short_reinf_lst = []
        for reinf in self.reinforcement_lst:
            if reinf.__class__ == ContinuousFibers:
                cont_reinf_lst.append(reinf)
            elif reinf.__class__ == ShortFibers:
                short_reinf_lst.append(reinf)
        return cont_reinf_lst, short_reinf_lst

    cont_fibers_instance = Instance(CrackBridgeContFibersGf)

    def _cont_fibers_instance_default(self):
        return CrackBridgeContFibersGf()

    cont_fibers = Property(Instance(CrackBridgeContFibersGf),
                           depends_on='reinforcement_lst+,Ll,Lr,E_m,w')

    @cached_property
    def _get_cont_fibers(self):
        cbcf = self.cont_fibers_instance
        cbcf.w = self.w
        cbcf.Ll = self.Ll
        cbcf.Lr = self.Lr
        cbcf.E_m = self.E_m
        cbcf.E_c = self.E_c
        cbcf.w_unld = self.w_unld
        cbcf.cont_reinf_lst = self.sorted_reinf_lst[0]
        cbcf.epsm_softening = self.epsm_softening
        # print self.w_unld, self.w, self.epsm_softening, self.ft, self.Gf
        return cbcf

    short_fibers_instance = Instance(CrackBridgeShortFibersGf)

    def _short_fibers_instance_default(self):
        return CrackBridgeShortFibersGf()

    short_fibers = Property(
        Instance(CrackBridgeShortFibersGf), depends_on='reinforcement_lst+,E_m,w')

    @cached_property
    def _get_short_fibers(self):
        cbsf = self.short_fibers_instance
        cbsf.w = self.w
        cbsf.E_m = self.E_m
        cbsf.E_c = self.E_c
        cbsf.short_reinf_lst = self.sorted_reinf_lst[1]
        cbsf.epsm_softening = self.epsm_softening
        return cbsf

    _x_arr = Property(Array, depends_on='w,E_m,Ll,Lr,reinforcement_lst+')

    @cached_property
    def _get__x_arr(self):
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) != 0:
            added_x = np.hstack(
                (self.cont_fibers.x_arr, self.short_fibers.x_arr))
            sorted_unique_x = np.unique(added_x)
            return sorted_unique_x
        elif len(self.sorted_reinf_lst[0]) != 0:
            return self.cont_fibers.x_arr
        elif len(self.sorted_reinf_lst[1]) != 0:
            return self.short_fibers.x_arr

    _epsm_arr = Property(Array, depends_on='w,E_m,Ll,Lr,reinforcement_lst+')

    @cached_property
    def _get__epsm_arr(self):
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) != 0:
            epsm_cont_interp = MFnLineArray(
                xdata=self.cont_fibers.x_arr, ydata=self.cont_fibers.epsm_arr)
            epsm_short_interp = MFnLineArray(
                xdata=self.short_fibers.x_arr, ydata=self.short_fibers.epsm_arr)
            added_epsm_cont = self.cont_fibers.epsm_arr + \
                epsm_short_interp.get_values(self.cont_fibers.x_arr)
            added_epsm_short = self.short_fibers.epsm_arr + \
                epsm_cont_interp.get_values(self.short_fibers.x_arr)
            sorted_unique_idx = np.unique(np.hstack(
                (self.cont_fibers.x_arr, self.short_fibers.x_arr)), return_index=True)[1]
            return np.hstack((added_epsm_cont, added_epsm_short))[sorted_unique_idx]
        elif len(self.sorted_reinf_lst[0]) != 0:
            return self.cont_fibers.epsm_arr
        elif len(self.sorted_reinf_lst[1]) != 0:
            self.short_fibers.w = self.w
            return self.short_fibers.epsm_arr

    _epsf_arr = Property(Array, depends_on='w,E_m,Ll,Lr,reinforcement_lst+')

    @cached_property
    def _get__epsf_arr(self):
        ''' only for continuous reinforcement '''
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) == 0:
            self.cont_fibers.w = self.w
            return self.cont_fibers.epsf_arr
        else:
            raise ValueError('epsf can only be computed for continuous fibers')

    _epsf0_arr = Property(Array, depends_on='w,E_m,Ll,Lr,reinforcement_lst+')

    @cached_property
    def _get__epsf0_arr(self):
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) != 0:
            epsf0_cont = self.cont_fibers.epsf0_arr
            epsf0_short = self.short_fibers.epsf0_arr
        elif len(self.sorted_reinf_lst[0]) != 0:
            epsf0_cont = self.cont_fibers.epsf0_arr
            epsf0_short = np.array([])
        elif len(self.sorted_reinf_lst[1]) != 0:
            epsf0_cont = np.array([])
            epsf0_short = self.short_fibers.epsf0_arr
        return epsf0_cont, epsf0_short

    _epsf0_arr_cont = Property(
        Array, depends_on='w,E_m,Ll,Lr,reinforcement_lst+')

    @cached_property
    def _get__epsf0_arr_cont(self):
        return self._epsf0_arr[0]

    _epsf0_arr_short = Property(
        Array, depends_on='w,E_m,Ll,Lr,reinforcement_lst+')

    @cached_property
    def _get__epsf0_arr_short(self):
        return self._epsf0_arr[1]

    sigma_c = Property(depends_on='w,E_m,Ll,Lr,reinforcement_lst+')

    @cached_property
    def _get_sigma_c(self):
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) != 0:
            sigma_c_cont = np.sum(self._epsf0_arr_cont * self.cont_fibers.sorted_stats_weights *
                                  self.cont_fibers.sorted_V_f * self.cont_fibers.sorted_nu_r *
                                  self.cont_fibers.sorted_E_f * (1. - self.cont_fibers.damage))
            sigma_c_short = np.sum(self._epsf0_arr_short * self.short_fibers.sorted_V_f *
                                   self.short_fibers.sorted_E_f)
        elif len(self.sorted_reinf_lst[0]) != 0:
            sigma_c_cont = np.sum(self._epsf0_arr_cont * self.cont_fibers.sorted_stats_weights *
                                  self.cont_fibers.sorted_V_f * self.cont_fibers.sorted_nu_r *
                                  self.cont_fibers.sorted_E_f * (1. - self.cont_fibers.damage))
            sigma_c_short = 0.0
        elif len(self.sorted_reinf_lst[1]) != 0:
            sigma_c_cont = 0.0
            sigma_c_short = np.sum(self._epsf0_arr_short * self.short_fibers.sorted_V_f *
                                   self.short_fibers.sorted_E_f)
        return sigma_c_cont + sigma_c_short + self.epsm_softening * self.E_m * (1. - self.V_f_tot)

    secant_K = Property(depends_on='w,E_m,Ll,Lr,reinforcement_lst+')

    @cached_property
    def _get_secant_K(self):
        ''' secant stiffness at given w '''
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) == 0:
            self.cont_fibers.w = self.w
            ef0, a_short, a_long, em, a = self.cont_fibers.profile(
                self.cont_fibers.damage)
            K_cont = np.sum(self.cont_fibers.sorted_stats_weights *
                            self.cont_fibers.sorted_V_f * self.cont_fibers.sorted_nu_r *
                            self.cont_fibers.sorted_E_f * (1. - self.cont_fibers.damage) /
                            (a_short + a_long))
            return K_cont
        else:
            raise ValueError(
                'secant stiffness not yet implemented for short fibers')


class CompositeCrackBridgeViewGf(ModelView):

    model = Instance(CompositeCrackBridgeGf)
    results = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_results(self):
        if self.model.w <= 0.0:
            self.model.w = 1e-15
        sigma_c = self.model.sigma_c
#         Kf_broken = np.sum(self.model.cont_fibers.sorted_V_f * self.model.cont_fibers.sorted_nu_r *
#                            self.model.cont_fibers.sorted_stats_weights * self.model.cont_fibers.sorted_E_f *
#                            self.model.cont_fibers.damage)
        if self.model.Ll > self.model.Lr:
            return -self.model._x_arr[::-1], self.model._epsm_arr[::-1], sigma_c, self.model._epsf_arr[::-1]
        else:
            return self.model._x_arr, self.model._epsm_arr, sigma_c, self.model._epsf_arr

    x_arr = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_x_arr(self):
        return self.results[0]

    epsm_arr = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_epsm_arr(self):
        return self.results[1]

    epsf_arr = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_epsf_arr(self):
        return self.results[3]

    sigma_c = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_sigma_c(self):
        return self.results[2]

    def sigma_c_arr(self, w_arr, u=False, damage=False):
        sigma_c_lst = []
        u_lst = []
        damage_lst = []
        for i, w in enumerate(w_arr):
            self.model.w = w
            sigma_c_lst.append(self.sigma_c)
            if u == True:
                u_lst.append(self.u_evaluated)
            if damage == True:
                damage_lst.append(np.sum(self.model.cont_fibers.damage *
                                         self.model.cont_fibers.sorted_stats_weights *
                                         self.model.cont_fibers.sorted_nu_r))
        if u == True or damage == True:
            return np.array(sigma_c_lst), np.array(u_lst), np.array(damage_lst)
        else:
            return np.array(sigma_c_lst)

    def secant_K(self, w_arr):
        secant_K_lst = []
        for w_i in w_arr:
            self.model.w = w_i
            secant_K_lst.append(self.model.secant_K)
        return np.array(secant_K_lst)

    u_evaluated = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_u_evaluated(self):
        return self.model.w + np.trapz(self.epsm_arr, self.x_arr)

    sigma_c_max = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_sigma_c_max(self):
        frozen_w_unld = self.model.w_unld
        frozen_damage = self.model.cont_fibers.maximum_damage

        def minfunc_sigma(w):
            self.model.w_unld = frozen_w_unld
            self.model.cont_fibers.maximum_damage = frozen_damage
            self.model.w = w
            stiffness_loss = np.sum(
                self.model.cont_fibers.Kf * self.model.cont_fibers.damage) / np.sum(self.model.cont_fibers.Kf)
            if stiffness_loss > 0.90:
                return 1. + w
            #plt.plot(w, self.sigma_c, 'ro')
            return -self.sigma_c

        def residuum_stiffness(w):
            self.model.w_unld = frozen_w_unld
            self.model.cont_fibers.maximum_damage = frozen_damage
            self.model.w = w
            stiffness_loss = np.sum(
                self.model.Kf * self.model.damage) / np.sum(self.model.Kf)
            if stiffness_loss > 0.90:
                return 1. + w
            if stiffness_loss < 0.65 and stiffness_loss > 0.45:
                residuum = 0.0
            else:
                residuum = stiffness_loss - 0.5
            return residuum

        if len(self.model.sorted_reinf_lst[0]) == 0:
            # there are only short fibers
            def minfunc_short_fibers(w):
                self.model.w_unld = frozen_w_unld
                self.model.cont_fibers.maximum_damage = frozen_damage
                self.model.w = w
                return -self.sigma_c
            w_max = fminbound(minfunc_short_fibers, 0.0,
                              3.0, maxfun=10, disp=0)
            return self.sigma_c, w_max
        else:
            # continuous or mixed fibers
            try:
                w_max = brentq(residuum_stiffness, 0.0, min(
                    0.1 * (self.model.Ll + self.model.Lr), 20.))
            except:
                w_max = 0.03 * (self.model.Ll + self.model.Lr)
            w_points = np.linspace(0, w_max, len(
                self.model.reinforcement_lst) + 1)
            w_maxima = []
            sigma_maxima = []
            for i, w in enumerate(w_points[1:]):
                w_maxima.append(
                    fminbound(minfunc_sigma, w_points[i], w_points[i + 1], maxfun=10, disp=0))
                sigma_maxima.append(self.sigma_c)
            return sigma_maxima[np.argmax(np.array(sigma_maxima))], w_maxima[np.argmax(np.array(sigma_maxima))]

    def apply_load(self, sigma):
        if sigma > self.sigma_c_max[0]:
            raise ValueError(
                'applied load ', sigma, 'MPa is larger than composite strength ', self.sigma_c_max[0], 'MPa')
        else:
            def residuum(w):
                self.model.w = float(w)
                return sigma - self.sigma_c
            brentq(residuum, 0.0, min(self.sigma_c_max[1], 20.))

    def sigma_f_lst(self, w_arr):
        sigma_f_arr = np.zeros(len(w_arr) *
                               len(self.model.reinforcement_lst)).reshape(len(w_arr),
                                                                          len(self.model.reinforcement_lst))
        masks = [((self.model.sorted_xi == reinf.xi) *
                  (self.model.sorted_E_f == reinf.E_f) *
                  (self.model.sorted_V_f == reinf.V_f))
                 for reinf in self.model.reinforcement_lst]
        for i, w in enumerate(w_arr):
            if w == 0.0:
                self.model.w = 1e-15
            else:
                self.model.w = w
            self.model.damage
            for j, reinf in enumerate(self.model.reinforcement_lst):
                sigma_fi = np.sum(self.model._epsf0_arr * self.model.sorted_stats_weights * self.model.sorted_nu_r *
                                  self.model.sorted_E_f * (1. - self.model.damage) * masks[j])
                sigma_f_arr[i, j] = sigma_fi
        return sigma_f_arr

    Welm = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_Welm(self):
        Km = self.results[4]
        bonded_l = self.epsm_arr[0] ** 2 * Km * \
            (self.model.Ll - np.abs(self.x_arr[0]))
        bonded_r = self.epsm_arr[-1] ** 2 * Km * \
            (self.model.Lr - np.abs(self.x_arr[-1]))
        return 0.5 * (np.trapz(self.epsm_arr ** 2 * Km, self.x_arr) + bonded_l + bonded_r)

    Welf = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_Welf(self):
        Kf = self.model.E_c - self.results[4]
        bonded_l = self.mu_epsf_arr[0] ** 2 * Kf * \
            (self.model.Ll - np.abs(self.x_arr[0]))
        bonded_r = self.mu_epsf_arr[-1] ** 2 * Kf * \
            (self.model.Lr - np.abs(self.x_arr[-1]))
        return 0.5 * (np.trapz(self.mu_epsf_arr ** 2 * Kf, self.x_arr) + bonded_l + bonded_r)

    W_el_tot = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_W_el_tot(self):
        '''total elastic energy stored in the specimen'''
        return self.Welf + self.Welm

    W_inel_tot = Property(
        depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')

    @cached_property
    def _get_W_inel_tot(self):
        '''total inelastic energy dissipated during loading up to w'''
        return self.U - self.W_el_tot

    U_line = Property(
        depends_on='model.E_m, model.Ll, model.Lr, model.reinforcement_lst+, w_arr_energy')

    @cached_property
    def _get_U_line(self):
        '''work done by external force - mfn_line'''
        w_arr = self.w_arr_energy
        u_lst = []
        F_lst = []
        for w in w_arr:
            self.model.w = w
            u_lst.append(self.u_evaluated)
            F_lst.append(self.sigma_c)
        u_arr = np.array(u_lst)
        F_arr = np.array(F_lst)
        U_line = MFnLineArray(xdata=w_arr, ydata=np.hstack(
            (0, cumtrapz(F_arr, u_arr))))
        return U_line

    U = Property(
        depends_on='model.E_m, model.Ll, model.Lr, model.reinforcement_lst+, model.w')

    @cached_property
    def _get_U(self):
        '''work done by external force U(w)'''
        return self.U_line.get_values(self.model.w)

    w_arr_energy = Array

    def get_sigma_m_x_input(self, sigma):
        self.apply_load(sigma)
        line = MFnLineArray(xdata=self.x_arr,
                            ydata=self.epsm_arr)
        return line.get_values(self.x_input)


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    from stats.pdistrib.weibull_fibers_composite_distr import fibers_MC
    from spirrid.rv import RV

    tau_scale = 1.53419049
    tau_shape = 0.90615
    tau_loc = 0.00
    xi_shape = 8.6
    xi_scale = 1.0114
    ft = 6.0
    Gf = 3.0

    reinf1 = ContinuousFibers(r=3.5e-3,
                              tau=RV('weibull_min', loc=0.01,
                                     scale=.1, shape=2.),
                              V_f=0.005,
                              E_f=200e3,
                              xi=fibers_MC(m=7., sV0=0.005),
                              label='carbon',
                              n_int=100)

    reinf_cont = ContinuousFibers(r=3.5e-3,
                                  tau=RV('gamma', loc=tau_loc,
                                         scale=tau_scale, shape=tau_shape),
                                  V_f=0.001,
                                  E_f=181e3,
                                  xi=fibers_MC(m=xi_shape, sV0=xi_scale),
                                  label='carbon',
                                  n_int=20)

    reinf_short = ShortFibers(bond_law='plastic',
                              r=3.5e-3,
                              tau=2.1,
                              V_f=0.03,
                              E_f=200e3,
                              xi=20.,
                              snub=0.5,
                              phi=RV('sin2x', scale=1.0, shape=0.0),
                              Lf=5.,
                              label='carbon',
                              n_int=201)

    ccb = CompositeCrackBridgeGf(E_m=25e3,
                                 ft=ft,
                                 Gf=Gf,
                                 reinforcement_lst=[reinf_cont],
                                 Ll=50.,
                                 Lr=50.,
                                 w=.1)

    for i, depsf in enumerate(ccb.cont_fibers.sorted_depsf):
        epsf_x = np.maximum(
            ccb._epsf0_arr[0][i] - depsf * np.abs(ccb._x_arr), ccb._epsm_arr)
#             if i == 0:
#                 plt.plot(ccb._x_arr, epsf_x, color='blue', label='fibers')
#             else:
        #plt.plot(ccb._x_arr, epsf_x, color='black', alpha=1-0.5*ccb.cont_fibers.damage[i])

    epsf0_combined = np.hstack((ccb._epsf0_arr[0], ccb._epsf0_arr[1]))
    #plt.plot(np.zeros_like(epsf0_combined), epsf0_combined, 'ro', label='maximum')
    plt.plot(ccb._x_arr, ccb._epsm_arr, lw=2, color='blue', label='matrix')
    #plt.plot(ccb._x_arr, ccb._epsf_arr, lw=2, color='red', label='mean fiber')
    plt.legend(loc='best')
    plt.ylabel('matrix and fiber strain [-]')
#     plt.ylabel('long. position [mm]')
    plt.show()
