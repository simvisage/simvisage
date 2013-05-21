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
# Created on Jun 14, 2010 by: rch
import time
from etsproxy.traits.api import \
    Float, Str, implements, Range, List, Property, cached_property, Tuple, Array

import numpy as np

from spirrid.i_rf import \
    IRF

from spirrid.rf import \
    RF

from matplotlib import pyplot as plt

from scipy.stats import weibull_min
from scipy.optimize import fsolve


def H(x):
    return x >= 0.0


class CBEMClampedFiberStressResidual(RF):
    '''
    Crack bridged by a fiber with constant frictional interface
    to an elastic matrix; clamped fiber end; residual stress
    carried by broken filaments.
    '''

    implements(IRF)

    title = Str('crack bridge - clamped fiber with constant friction')

    Pf = Range(0, 1, auto_set=False, enter_set=True, input=True,
                distr=['uniform'])

    m = Float(5.0, auto_set=False, enter_set=True, input=True,
                distr=['weibull_min', 'uniform'],
                desc='filament Weibull shape parameter')

    s0 = Float(.02, auto_set=False, enter_set=True, input=True,
                distr=['weibull_min', 'uniform'],
                desc='filament Weibull scale parameter at l = 10 mm')

    tau = Float(2.5, auto_set=False, enter_set=True, input=True,
                distr=['uniform', 'norm'])

    l = Float(10.0, auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='free length')

    r = Float(0.013, auto_set=False, input=True,
              enter_set=True, desc='fiber radius in mm')

    E_f = Float(72e3, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    E_m = Float(30e3, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'])

    V_f = Float(0.0175, auto_set=False, enter_set=True, input=True,
              distr=['uniform'])

    theta = Float(0.01, auto_set=False, enter_set=True, input=True,
                  distr=['uniform', 'norm'], desc='slack')

    phi = Float(1., auto_set=False, enter_set=True, input=True,
                  distr=['uniform', 'norm'], desc='bond quality')

    Ll = Float(1., auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='embedded length - left',
               ctrl_range=(0.0, 1.0, 10))

    Lr = Float(.5, auto_set=False, enter_set=True, input=True,
              distr=['uniform'], desc='embedded length - right',
               ctrl_range=(0.0, 1.0, 10))

    w = Float(auto_set=False, enter_set=True, input=True,
               distr=['uniform'], desc='crack width',
               ctrl_range=(0.0, 1.0, 10))

    x = Float(auto_set=False, enter_set=True, input=True,
               distr=['uniform'], desc='crack width',
               ctrl_range=(0.0, 1.0, 10))

    x_label = Str('crack opening [mm]')
    y_label = Str('force [N]')

    C_code = Str('')

    def crackbridge(self, w, l, T, Kf, Km, V_f):
        # Phase A : Both sides debonding .
        Kc = Kf + Km
        c = Kc * T * l / 2.
        q0 = (np.sqrt(c ** 2 + w * Kf * Km * Kc * T) - c) / Km
        return q0 / V_f

    def pullout(self, w, l, T, Kf, Km, V_f, Lmin, Lmax):
        # Phase B : Debonding of shorter side is finished
        Kc = Kf + Km
        c = Kc * T * (Lmin + l)
        f = T ** 2 * Lmin ** 2 * Kc ** 2
        q1 = (np.sqrt(c ** 2. + f + 2 * w * Kc * T * Kf * Km) - c) / Km
        return q1 / V_f

    def linel(self, w, l, T, Kf, Km, V_f, Lmax, Lmin):
        # Phase C: Both sides debonded - linear elastic behavior.
        Kc = Kf + Km
        q2 = (2. * w * Kf * Km + T * Kc * (Lmin ** 2 + Lmax ** 2)) \
                / (2. * Km * (Lmax + l + Lmin))
        return q2 / V_f

    # reshaped parameters for the evaluation of L_bond_mean
    params = List

    # cached array of mean bonded lengths of broken filaments
    L_bond_mean = Property(depends_on='CB_params')

    @cached_property
    def _get_L_bond_mean(self):
        return self.L_bond_mean_func(*self.params)

    def L_bond_mean_func(self, q, tau, x_n, l, E_f, E_m,
                         theta, Pf, phi, Ll, Lr, V_f, r, s0, m):
        '''numerical solver for filament breaks with evaluation
        of the mean remaining bonded lengths'''
        L_bond_x = np.abs(x_n) - l / 2.
        L_bond_x = L_bond_x * H(L_bond_x)
        # initial guess for fsolve - solution for a bond free filament
        q_init_guess = weibull_min(m, scale=s0).ppf(0.5) \
                                * E_f * np.ones_like(q)
        left_side = np.log(1. - Pf) * s0 ** m
        q_break = fsolve(self.opt_q0, q_init_guess.flatten(),
                         args=(q_init_guess.shape, tau, x_n, l,
                               V_f, E_f, E_m, theta, Pf, r, m, left_side),
                         xtol=0.01,
                         col_deriv=True)
        # load q at filament breakage
        q_break = q_break.reshape(q_init_guess.shape)
        # evaluation of the mean bonded length - weighted mean of all
        # all possible bonded lengths
        eps_n0 = self.eps_n(q_break, tau, x_n, l, V_f, E_f, E_m, theta, r)
        CDF_n0 = weibull_min(m, scale=s0).cdf(eps_n0)
        L_bond_mean = np.sum(CDF_n0 * L_bond_x, axis= -1) / (1e-15 +
                            np.sum(CDF_n0, axis= -1))
        return L_bond_mean

    def opt_q0(self, q0, q_shape, tau, x_n, l, V_f, E_f,
               E_m, theta, Pf, r, m, left_side):
        '''filament breaks are found by fsolve as roots of this function '''
        right = np.sum(self.eps_n(q0.reshape(q_shape), tau, x_n, l,
                    V_f, E_f, E_m, theta, r) ** m, axis= -1)
        value = left_side + right
        return value.flatten()

    def eps_n(self, q, tau, x_n, l, V_f, E_f, E_m, theta, r):
        '''evaluates the strain profiles along x_n'''
        # stress in the free length
        l = l * (1 + theta)
        q_l = q * H(l / 2. - abs(x_n))
        # stress in the part, where fiber transmits stress to the matrix
        Tf = 2. * tau / r
        q_e = (q - Tf * (abs(x_n) - l / 2.)) * H(abs(x_n) - l / 2.)
        # far field stress
        E_c = E_m * (1 - V_f) + E_f * V_f
        q_const = q * V_f * E_f / E_c
        # putting all parts together
        q_tau = q_l + q_e
        q_x = np.maximum(q_tau, q_const)
        # cutting out the part, where broken filaments can be pulled out
        eps_n = q_x / E_f * H(q_tau)
        return eps_n

    def fil_break(self, q, x_n, tau, l, E_f, E_m, theta, Pf,
                  phi, Ll, Lr, V_f, r, s0, m):
        '''evaluates the filament breakage probabilities chain_sf
        and calls the method for evaluating the mean bonded length'''
        # add another dimension (of x_n) to the parameters
        varlist = [q, tau, l, E_f, E_m, theta, Pf, phi, Ll, Lr, V_f, r, s0, m]
        reshaped = []
        for var in varlist:
            if isinstance(var, np.ndarray):
                shape = tuple(list(var.shape) + [1])
                var = var.reshape(shape)
            reshaped.append(var)
        q, tau, l, E_f, E_m, theta, Pf, phi, Ll, Lr, V_f, r, s0, m = reshaped
        # evaluate survival probability of the fiber along a crack bridge
        eps_n = self.eps_n(q, tau, x_n, l, V_f, E_f, E_m, theta, r)
        CDF_n = weibull_min(m, scale=s0).cdf(eps_n)
        # survival probability of the whole system along x_n
        sf = 1 - CDF_n
        chain_sf = sf.prod(axis= -1)
        # store the parameters for evaluation of the breaking load
        self.params = [q, tau, x_n, l, E_f, E_m, theta, Pf,
                       phi, Ll, Lr, V_f, r, s0, m]
        # residual bonded length
        L_bond_mean = self.L_bond_mean
        return chain_sf, L_bond_mean

    # params to be saved for traits notification on their change
    CB_params = Tuple

    def __call__(self, w, tau, l, E_f, E_m, theta, Pf,
                 phi, Ll, Lr, V_f, r, s0, m):
        self.CB_params = (tau, l, E_f, E_m, theta, Pf,
                          phi, Ll, Lr, V_f, r, s0, m)
        # assigning short and long embedded length
        Lmin = np.minimum(Ll, Lr)
        Lmax = np.maximum(Ll, Lr)
        Lmin = np.maximum(Lmin - l / 2., 1e-15)
        Lmax = np.maximum(Lmax - l / 2., 1e-15)
        # maximum condition for free length
        l = np.minimum(l / 2., Lr) + np.minimum(l / 2., Ll)
        # defining variables
        w = w - theta * l
        l = l * (1 + theta)
        w = H(w) * w
        T = 2. * tau * V_f / r
        Km = (1. - V_f) * E_m
        Kf = V_f * E_f
        Kc = Km + Kf

        # double sided debonding
        # q0 = self.crackbridge(w, l, T, Kr, Km, Lmin, Lmax)
        q0 = self.crackbridge(w, l, T, Kf, Km, V_f)

        # displacement at which the debonding to the closer clamp is finished
        w0 = (Lmin + l) * Lmin * Kc * T / Kf / Km

        # debonding of one side; the other side is clamped
        q1 = self.pullout(w, l, T, Kf, Km, V_f, Lmin, Lmax)

        # displacement at which the debonding is finished at both sides
        e1 = Lmax * Kc * T / Km / Kf
        w1 = e1 * (l + Lmax / 2.) + (e1 + e1 * Lmin / Lmax) * Lmin / 2.

        # debonding completed at both sides, response becomes linear elastic
        q2 = self.linel(w, l, T, Kf, Km, V_f, Lmax, Lmin)

        # cut out definition ranges
        q0 = H(w) * (q0 + 1e-15) * H(w0 - w)
        q1 = H(w - w0) * q1 * H(w1 - w)
        q2 = H(w - w1) * q2

        # add all parts
        q = q0 + q1 + q2
        # include breaking strain
        n = (Ll + Lr) / 10.
        if n < 1.0:
            n = 1
        else:
            n = int(n)
        # within a crack bridge divide the filament
        # into parts with constant strengths
        x_n = np.linspace(-Ll, Lr, n).reshape(tuple([1] * len(q.shape) + [n]))
        # create an array of strains along these parts
        chain_sf, L_bond_mean = self.fil_break(q, x_n, tau,
            l, E_f, E_m, theta, Pf, phi, Ll, Lr, V_f, r, s0, m)
        # print np.unravel_index(a.argmin(axis = -1), a.shape)
        surv = q * H(Pf - (1 - chain_sf))
        broken = L_bond_mean * T / V_f * H((1 - chain_sf) - Pf)
        q = surv + broken
        return q


class CBEMClampedFiberStressResidualSP(CBEMClampedFiberStressResidual):
        '''
        Stress profile for a crack bridged by a fiber with constant
        frictional interface to an elastic matrix; clamped fiber end;
        residual stress carried by broken filaments.
        '''
        x = Float(0.0, auto_set=False, enter_set=True, input=True,
                  distr=['uniform'], desc='distance from crack')

        x_label = Str('position [mm]')
        y_label = Str('force [N]')

        C_code = Str('')

        def __call__(self, w, x, tau, l, E_f, E_m, theta, Pf,
                     phi, Ll, Lr, V_f, r, s0, m):
            Tf = 2. * tau / r
            q = super(CBEMClampedFiberStressResidualSP, self).__call__(w,
                    tau, l, E_f, E_m, theta, Pf, phi, Ll, Lr, V_f, r, s0, m)
            self.q = q
            # stress in the free length
            l = l * (1 + theta)
            q_l = q * H(l / 2 - abs(x))
            # stress in the part, where fiber transmits stress to the matrix
            q_e = (q - Tf * (abs(x) - l / 2.)) * H(abs(x) - l / 2.)
            # q_e = q_e * H(x + Ll) * H (Lr - x)

            # far field stress
            E_c = E_m * (1 - V_f) + E_f * V_f
            q_const = q * V_f * E_f / E_c

            # putting all parts together
            q_x = q_l + q_e
            q_x = np.maximum(q_x, q_const)
            return q_x

        q = Array(Float)

if __name__ == '__main__':

    r = 0.00345
    V_f = 0.0103
    t = .1
    Ef = 200e3
    Em = 25e3
    l = 10.
    theta = 0.
    Pf = 0.5
    phi = 1.
    Ll = 40.
    Lr = 20.
    s0 = 0.02
    m = 5.0

    def Pw():
        plt.figure()
        w = np.linspace(0, 1, 500)
        P = CBEMClampedFiberStressResidual()
        q = P(w, t, l, Ef, Em, theta, Pf, phi, Ll, Lr, V_f, r, s0, m)
        plt.plot(w, q, lw=2, ls='-', color='black', label='CB_emtrx_stress')
        plt.legend(loc='best')
        plt.ylim(0,)
        plt.ylabel('stress', fontsize=14)
        plt.xlabel('w', fontsize=14)
        plt.title('Pullout Resp Func Clamped Fiber EMTRX Residual Stress')

    def SP(w, x, t, l, Ef, Em, theta, Pf, phi, Ll, Lr, V_f, r, s0, m):
        cbcsp = CBEMClampedFiberStressResidualSP()
        x = np.linspace(-40, 40, 300)
        w = .4
        q = cbcsp(w, x, t, l, Ef, Em, theta, Pf, phi, Ll, Lr, V_f, r, s0, m)
        plt.plot(x, q, lw=2, color='black', label='stress along filament')
        plt.ylabel('stress', fontsize=14)
        plt.xlabel('position', fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.title('Stress Along Filament at w = %.1f' % w)
        # plt.legend(loc='best')
        plt.ylim(0,)


    Pw()
    #SP()
    plt.show()
