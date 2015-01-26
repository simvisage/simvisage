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
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import ContinuousFibers
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge


class SCMView(ModelView):

    model = Instance(SCM)

    def crack_widths(self, sigma_c):
        # find the index of the nearest value in the load range
        idx = np.abs(self.model.load_sigma_c_arr - sigma_c).argmin()
        # pick the cracks that emerged at the given load
        current_cb_list = self.model.cb_list(sigma_c)
        if current_cb_list[0] is not None:
            w_lst = []
            for CB in current_cb_list:
                w_lst.append(CB.get_w_sigma_c(sigma_c))
            return np.array(w_lst)
        else:
            return np.array(0.0, ndmin=1)
    
    def crack_density(self, sigma_c):
        # pick the cracks that emerged at the given load
        cb_load = self.model.cb_list(sigma_c)
        if cb_load[0] is not None:
            return float(len(cb_load)) / self.model.length
        else:
            return np.array(0.0, ndmin=1) / self.model.length

    w_density = Property(List, depends_on='model')
    @cached_property
    def _get_w_density(self):
        return [self.crack_density(load) for load in self.model.load_sigma_c_arr]

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

    eps_sigma = Property(depends_on='model.')
    @cached_property
    def _get_eps_sigma(self):
        eps_lst = []
        u_m_tot = np.trapz(self.eps_m_x, self.x_area, axis=1)
        for i, sig in enumerate(self.model.load_sigma_c_arr):
            w_arr_i = np.array(self.eval_w[i])
            eps_lst.append((np.sum(w_arr_i) + u_m_tot[i])/self.model.length)
        eps = np.array(eps_lst)
        eps = eps[np.isnan(eps) == False]
        if len(eps) != len(self.model.load_sigma_c_arr):
            eps = list(eps) + [list(eps)[-1]]
            sigma = copy.copy(self.model.load_sigma_c_arr[:len(eps)])
            sigma[-1] = 0.0
            return eps, sigma
        else:
            return eps, self.model.load_sigma_c_arr

if __name__ == '__main__':
    from stats.pdistrib.weibull_fibers_composite_distr import fibers_MC
    length = 500.
    nx = 2000
    random_field = RandomField(seed=False,
                               lacor=1.,
                               length=length,
                               nx=500,
                               nsim=1,
                               loc=.0,
                               shape=15.,
                               scale=4.0,
                               distr_type='Weibull'
                               )

    reinf1 = ContinuousFibers(r=3.5e-3,
                              tau=RV('weibull_min', loc=0.01, scale=.1, shape=2.),
                              V_f=0.005,
                              E_f=200e3,
                              xi=fibers_MC(m=7., sV0=0.005),
                              label='carbon',
                              n_int=500)

    reinf2 = ContinuousFibers(r=3.5e-3,
                              tau=RV('weibull_min', loc=0.01, scale=.1, shape=2.),
                              V_f=0.005,
                              E_f=200e3,
                              xi=fibers_MC(m=7., sV0=0.005),
                              label='carbon',
                              n_int=500)

    CB_model = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf1, reinf2],
                                 )

    scm = SCM(length=length,
              nx=nx,
              random_field=random_field,
              CB_model=CB_model,
              load_sigma_c_arr=np.linspace(0.01, 20., 100),
              )

    scm_view = SCMView(model=scm)
    scm_view.model.evaluate()

    def plot():
        eps, sigma = scm_view.eps_sigma
        plt.plot(eps, sigma, color='black', lw=2, label='model')
#         plt.legend(loc='best')
#         plt.xlabel('composite strain [-]')
#         plt.ylabel('composite stress [MPa]')
#         plt.figure()
#         plt.hist(scm_view.crack_widths(15.), bins=20, label='load = 15 MPa')
#         plt.hist(scm_view.crack_widths(10.), bins=20, label='load = 10 MPa')
#         plt.hist(scm_view.crack_widths(5.), bins=20, label='load = 5 MPa')
#         plt.ylabel('frequency [-]')
#         plt.xlabel('crack width [mm]') 
#         plt.legend(loc='best')
#         plt.xlim(0)
#         plt.figure()
#         plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_mean,
#                  color='green', lw=2, label='mean crack width')
#         plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_median,
#                 color='blue', lw=2, label='median crack width')
#         plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_mean + scm_view.w_stdev,
#                 color='black', label='stdev')
#         plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_mean - scm_view.w_stdev,
#                 color='black')
#         plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_max,
#                  ls='dashed', color='red', label='max crack width')
#         plt.ylabel('crack width [mm]')
#         plt.xlabel('composite stress [MPa]')
#         plt.legend(loc='best')
#         plt.figure()
#         plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_density,
#                  color='black', lw=2, label='crack density')
#         plt.legend(loc='best')
#         plt.ylabel('crack density [1/mm]')
#         plt.xlabel('composite stress [MPa]')
        plt.show()

    plot()