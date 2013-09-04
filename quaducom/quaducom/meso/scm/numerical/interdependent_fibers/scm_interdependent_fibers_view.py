'''
Created on Jul 26, 2012

@author: rostar
'''

from etsproxy.traits.api import \
    Instance, Array, List, cached_property, Property
from matplotlib import pyplot as plt
from etsproxy.traits.ui.api import ModelView
from spirrid.rv import RV
from stats.misc.random_field.random_field_1D import RandomField
import numpy as np
import copy
from scm_interdependent_fibers_model import SCM
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import Reinforcement, ContinuousFibers
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx_view import CompositeCrackBridgeView


class SCMView(ModelView):

    model = Instance(SCM)

    def crack_widths(self, sigma_c):
        # find the index of the nearest value in the load range
        idx = np.abs(self.model.load_sigma_c_arr - sigma_c).argmin()
        # evaluate the relative strain e_rel between fibers
        # and matrix for the given load
        e_rel = self.mu_epsf_x[idx, :] - self.eps_m_x[idx, :]
        # pick the cracks that emerged at the given load
        cb_load = self.model.cb_list(sigma_c)
        if cb_load[0] is not None:
            # find the symmetry points between cracks as
            # the 0 element of their x range
            idxs = []
            for cb in cb_load:
                idxs.append(np.where(cb.position +
                                cb.x[0] == self.model.x_arr)[0])
            # add the index of the last point
            idxs.append(self.model.nx - 1)
            # list of crack widths to be filled in a loop with integrated e_rel
            crack_widths = [np.trapz(e_rel[idx:idxs[i + 1]],
                            self.model.x_arr[idx:idxs[i + 1]])
                            for i, idx in enumerate(idxs[:-1])]
            return np.array(crack_widths, ndmin=1)
        else:
            return np.array(0.0, ndmin=1)

    eval_w = Property(List, depends_on='model')
    @cached_property
    def _get_eval_w(self):
        return [self.crack_widths(load) for load in self.model.load_sigma_c_arr]

    w_mean = Property(Array, depends_on='model')
    @cached_property
    def _get_w_mean(self):
        return np.array([np.mean(w) for w in self.eval_w])

    w_median = Property(Array, depends_on='model')
    @cached_property
    def _get_w_median(self):
        return np.array([np.median(w) for w in self.eval_w])

    w_stdev = Property(Array, depends_on='model')
    @cached_property
    def _get_w_stdev(self):
        return np.array([np.std(w) for w in self.eval_w])

    w_max = Property(Array, depends_on='model')
    @cached_property
    def _get_w_max(self):
        return np.array([np.max(w) for w in self.eval_w])

    x_area = Property(depends_on='model.')
    @cached_property
    def _get_x_area(self):
        return  np.ones_like(self.model.load_sigma_c_arr)[:, np.newaxis] \
            * self.model.x_arr[np.newaxis, :]

    sigma_m_x = Property(depends_on='model.')
    @cached_property
    def _get_sigma_m_x(self):
        sigma_m_x = np.zeros_like(self.model.load_sigma_c_arr[:, np.newaxis]
                                  * self.model.x_arr[np.newaxis, :])
        for i, q in enumerate(self.model.load_sigma_c_arr):
            sigma_m_x[i, :] = self.model.sigma_m(q)
        return sigma_m_x

    eps_m_x = Property(Array, depends_on='model.')
    @cached_property
    def _get_eps_m_x(self):
        return self.sigma_m_x / self.model.CB_model.E_m

    mu_epsf_x = Property(depends_on='model.')
    @cached_property
    def _get_mu_epsf_x(self):
        mu_epsf_x = np.zeros_like(self.model.load_sigma_c_arr[:, np.newaxis]
                                  * self.model.x_arr[np.newaxis, :])
        for i, q in enumerate(self.model.load_sigma_c_arr):
            mu_epsf_x[i, :] = self.model.epsf_x(q)
        return mu_epsf_x

    eps_sigma = Property(depends_on='model.')
    @cached_property
    def _get_eps_sigma(self):
        eps = np.trapz(self.mu_epsf_x, self.x_area, axis=1) / self.model.length
        eps = eps[np.isnan(eps) == False]
        if len(eps) != len(self.model.load_sigma_c_arr):
            eps = list(eps) + [list(eps)[-1]]
            sigma = copy.copy(self.model.load_sigma_c_arr[:len(eps)])
            sigma[-1] = 0.0
            return eps, sigma
        else:
            return eps, self.model.load_sigma_c_arr

if __name__ == '__main__':
    length = 1000.
    nx = 1000
    random_field = RandomField(seed=True,
                               lacor=5.,
                                xgrid=np.linspace(0., length, 400),
                                nsim=1,
                                loc=.0,
                                shape=15.,
                                scale=2.5,
                                non_negative_check=True,
                                distribution='Weibull'
                               )

    reinf = ContinuousFibers(r=0.0035,
                          tau=RV('weibull_min', loc=0.006, shape=1.2, scale=.03),
                          V_f=0.011,
                          E_f=240e3,
                          xi=WeibullFibers(shape=5.0, sV0=0.0026),
                          label='carbon')

    CB_model = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf],
                                 )

    scm = SCM(length=length,
              nx=nx,
              random_field=random_field,
              CB_model=CB_model,
              load_sigma_c_arr=np.linspace(0.01, 15., 100),
              )

    scm_view = SCMView(model=scm)
    scm_view.model.evaluate()

    def plot():
        eps, sigma = scm_view.eps_sigma
        plt.figure()
        plt.plot(eps, sigma, color='black', lw=2, label='model')
        plt.legend(loc='best')
        plt.xlabel('composite strain [-]')
        plt.ylabel('composite stress [MPa]')
#         plt.figure()
#         plt.hist(scm_view.crack_widths(16.), bins=20, label='load = 20 MPa')
#         plt.hist(scm_view.crack_widths(13.), bins=20, label='load = 15 MPa')
#         plt.hist(scm_view.crack_widths(10.), bins=20, label='load = 10 MPa')
#         plt.legend(loc='best')
        plt.figure()
        plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_mean,
                 color='green', lw=2, label='mean crack width')
        plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_median,
                 color='blue', lw=2, label='median crack width')
        plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_mean + scm_view.w_stdev,
                 color='black', label='stdev')
        plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_mean - scm_view.w_stdev,
                 color='black')
        plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_max,
                 ls='dashed', color='red', label='max crack width')
        plt.legend(loc='best')
        plt.show()
    plot()
