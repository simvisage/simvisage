'''
Created on Apr 20, 2012
stress - strain relation for composites according to Ahn and Curtin 1996
@author: rostar
'''
from material import Material
from etsproxy.traits.api import Float, Property, cached_property, Instance, HasTraits
from scipy.stats import weibull_min
import numpy as np
from math import e
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from scipy.integrate import odeint
import pickle


def dirac(x):
    return np.sign(x == x[np.argmin(np.abs(x))]) / (x[1] - x[0])


class Curtin(HasTraits):

    material = Instance(Material)

    L = Property()

    def _get_S(self):
        '''distribution of the short fragments'''
        return weibull_min(8, scale=10., loc=1.)

    M = Property()

    def _get_M(self):
        '''distribution of the medium fragments'''
        return weibull_min(10, scale=15., loc=3.)

    S = Property()

    def _get_L(self):
        '''distribution of the long fragments'''
        return weibull_min(12, scale=20., loc=5.)

    alpha = Property()

    def _get_alpha(self):
        return self.material.V_m * self.material.E_m / \
            self.material.V_f / self.material.E_c

    delta = Property()

    def _get_delta(self):
        return self.material.r / 2. / \
            self.material.tau * self.alpha * self.material.sigma

    diff_delta = Property()

    def _get_diff_delta(self):
        return self.delta / self.material.sigma

    sigma_fL = Property()

    def _get_sigma_fL(self):
        return self.delta * self.alpha * self.material.sigma * self.L.moment(1) + \
            self.material.sigma * self.material.E_f / \
            self.material.E_c * self.L.moment(2)

    sigma_fM = Property()

    def _get_sigma_fM(self):
        return self.material.sigma / self.material.V_f * self.M.moment(1) - \
            self.alpha / 4. / self.delta * \
            self.material.sigma * self.M.moment(2)

    sigma_fS = Property()

    def _get_sigma_fS(self):
        return self.material.sigma / self.material.V_f * self.S.moment(1) - \
            self.alpha / 4. / self.delta * \
            self.material.sigma * self.S.moment(2)

    eps = Property()

    def _get_eps(self):
        return 1. / self.Lc / self.material.E_f * (self.sigma_fL + self.sigma_fM + self.sigma_fS)

# ===============================
#     auxiliary function
# ===============================

    def ints(self, t):
        s = np.linspace(1e-10, t, 500)
        return np.trapz((1. - e ** (-s)) / s, s)

    def intt(self, psi):
        t = np.linspace(1e-10, psi, 500)
        vect_ints = np.vectorize(self.ints)
        return np.trapz(e ** (-2 * vect_ints(t)), t)

    psi_line = Property()

    @cached_property
    def _get_psi_line(self):
        try:
            psi = open('psi_Curtin', 'r')
            line = pickle.load(psi)
        except:
            eta_psi_vect = np.vectorize(self.intt)
            psi = np.linspace(1e-3, 4.99, 300)
            eta_psi_func = eta_psi_vect(psi)
            eta2 = np.linspace(
                0.7476 - e ** (-2. * 0.577216) / 5, 0.74759, 500)
            psi2 = e ** (-2. * 0.577216) / (0.7476 - eta2)
            line = MFnLineArray(xdata=np.hstack((eta_psi_func, eta2)),
                                ydata=np.hstack((psi, psi2)))
        return line

    def psi(self, eta):
        return self.psi_line.get_values(eta, k=1)

    def diff_psi(self, eta):
        return self.psi_line.get_diffs(eta, k=1)

    def eta(self, psi):
        psi_line = MFnLineArray(xdata=self.psi_line.ydata,
                                ydata=self.psi_line.xdata)
        return psi_line.get_values(psi, k=1)

# ===============================


# distribution of long segments

    def PL(self, x, eta):
        L = self.psi(eta) ** 2 / self.diff_psi(eta) * \
            e ** (-(x / self.delta - 2.) * self.psi(eta))
        return L / eta / self.delta

# distribution of medium segments
    def PM(self, x, eta):
        eta_linsp = np.linspace(1e-10, eta, 100)
        if type(x) == type(1.0):
            M = 2. * np.trapz(self.psi(eta_linsp) *
                              e ** (-(x / self.delta - 1.) * self.psi(eta_linsp)), eta_linsp)
        else:
            M = 2. * np.trapz(self.psi(eta_linsp) *
                              e ** (-(x[:, np.newaxis] / self.delta - 1.) * self.psi(eta_linsp)), eta_linsp)
        return M.flatten() / eta / self.delta

    eta_comp = Float

    def NLPs_ode(self, U, sigma):
        xL = np.linspace(2 * self.delta, self.material.L_c, 100)
        eta = U[0] * self.delta / U[1]
        self.eta_comp = eta
        Lstar = U[0] * np.trapz((xL - 2 * self.delta) * self.PL(xL, eta), xL)
        print self.PL(self.delta, eta)
        dN = -U[0] * self.PM(np.array([self.delta]), eta) * \
            self.diff_delta + self.sigma_mu_distr(Lstar).pdf(sigma)
        dL = -U[0] * self.delta * self.PM(self.delta, eta) * self.diff_delta
        res = np.array([dN, dL]).flatten()
        print res
        return res

    def solveODE(self):
        return odeint(self.NLPs_ode, y0=np.array([0.01, self.material.L_c]), t=np.linspace(0.1, self.material.sigma, 50))

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    material = Instance(Material)
    material = Material(sigma=10.1,
                        E_f=72e3,
                        E_m=30e3,
                        tau=.1,
                        r=5e-4,
                        V_f=0.0175,
                        l_0=10,
                        L_c=1000
                        )

    #eta = np.linspace(0.001,0.74, 100)
    #psi = np.linspace(1e-10, 20., 100)
    c = Curtin(material=material)
    c.solveODE()
    x2 = np.linspace(c.delta, 2 * c.delta, 100)
    x3 = np.linspace(2 * c.delta, 10 * c.delta, 100)
#    plt.plot(eta, c.psi(eta), label = 'psi')
#    plt.plot(eta, e**(-2.*0.577216)/(0.7476-eta), label = 'psi_1')
#    plt.plot(eta, 1./(1-eta/0.7476)**0.7476, label = 'psi_2')
#    plt.plot(eta, eta + eta**2 + 7./6.*eta**3 + 13./9.* eta**4, label = 'psi_expansion')
#    plt.plot(eta, c.diff_psi(eta), label = 'diff psi')
#    plt.plot(psi, c.eta(psi), label = 'eta')
    plt.plot(x2, c.PM(x2, .2), label='M fragments')
    plt.plot(x3, c.PL(x3, .2), label='L fragments')
    plt.legend(loc='best')
    plt.show()
