'''
Created on Jul 26, 2012

@author: rostar
'''

from etsproxy.traits.api import \
    Instance, Array, List, cached_property, Property
from etsproxy.traits.ui.api import ModelView
from spirrid.rv import RV
from stats.misc.random_field.random_field_1D import RandomField
import numpy as np
import copy
from .scm_model import SCM
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge
import matplotlib.pyplot as plt


class SCMView(ModelView):

    model = Instance(SCM)

    def crack_widths(self, load):
        # pick the cracks that emerged at the given load
        if load > self.model.cracking_stress_lst[0]:
            cb_lst = self.model.get_current_cracking_state(load)
            w_lst = [cb_i.get_w(load) for cb_i in cb_lst]
            return np.array(w_lst)
        else:
            return np.array(0.0, ndmin=1)
    
    def crack_density(self, load):
        # pick the cracks that emerged at the given load
        if load > self.model.cracking_stress_lst[0]:
            return float(len(self.model.cracking_stress_lst)) / self.model.length
        else:
            return np.array([1. / self.model.length])

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
        for i, load in enumerate(self.model.load_sigma_c_arr): 
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

    eps_sigma_altern = Property(depends_on='model.')
    @cached_property
    def _get_eps_sigma_altern(self):
        eps_lst = []
        for load in self.model.load_sigma_c_arr: 
            if self.model.cracking_stress_lst[0] < load:
                cb_lst = self.model.get_current_cracking_state(load)
                u_i = 0.0
                for cb_i in cb_lst:
                    u_i += np.trapz(cb_i.get_epsf_x(load), cb_i.x)
                eps_lst.append(u_i/self.model.length)
            else:
                eps_lst.append(load / self.model.CB_model.E_c)
        eps = np.array(eps_lst)
        if len(eps) != len(self.model.load_sigma_c_arr):
            eps = list(eps) + [list(eps)[-1]]
            sigma = copy.copy(self.model.load_sigma_c_arr[:len(eps)])
            sigma[-1] = 0.0
            return eps, sigma
        else:
            return eps, self.model.load_sigma_c_arr
             
if __name__ == '__main__':
    from stats.pdistrib.weibull_fibers_composite_distr import fibers_MC
    from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import ContinuousFibers, ShortFibers
    from matplotlib import pyplot as plt
    length = 200.
    nx = 1000
    random_field = RandomField(seed=True,
                           lacor=1.,
                           length=length,
                           nx=1000,
                           nsim=1,
                           loc=.0,
                           shape=50.,
                           scale=3.4,
                           distr_type='Weibull')

    reinf_cont = ContinuousFibers(r=3.5e-3,
                              tau=RV('weibull_min', loc=0.01, scale=.1, shape=2.),
                              V_f=0.01,
                              E_f=200e3,
                              xi=fibers_MC(m=7., sV0=0.005),
                              label='carbon',
                              n_int=100)

    reinf_short = ShortFibers(bond_law = 'plastic',
                        r=.2,
                        tau=1.,
                        V_f=0.01,
                        E_f=200e3,
                        xi=10.,
                        snub=0.5,
                        phi=RV('sin2x', scale=1.0, shape=0.0),
                        Lf=20.,
                        label='short steel fibers',
                        n_int=50)

    CB_model = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf_cont, reinf_short],
                                 )
    scm = SCM(length=length,
              nx=nx,
              random_field=random_field,
              CB_model=CB_model,
              load_sigma_c_arr=np.linspace(0.01, 20., 200),
              n_BC_CB = 10)

    scm_view = SCMView(model=scm)
    scm_view.model.evaluate() 

    def plot():
        eps, sigma = scm_view.eps_sigma
        plt.plot(eps, sigma, color='black', lw=2, label='model')
        plt.legend(loc='best')
        plt.xlabel('composite strain [-]')
        plt.ylabel('composite stress [MPa]')
        plt.figure()
        plt.hist(scm_view.crack_widths(15.), bins=20, label='load = 15 MPa')
        plt.hist(scm_view.crack_widths(10.), bins=20, label='load = 10 MPa')
        plt.hist(scm_view.crack_widths(5.), bins=20, label='load = 5 MPa')
        plt.ylabel('frequency [-]')
        plt.xlabel('crack width [mm]') 
        plt.legend(loc='best')
        plt.xlim(0)
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
        plt.ylabel('crack width [mm]')
        plt.xlabel('composite stress [MPa]')
        plt.legend(loc='best')
        plt.figure()
        plt.plot(scm_view.model.load_sigma_c_arr, scm_view.w_density,
                 color='black', lw=2, label='crack density')
        plt.legend(loc='best')
        plt.ylabel('crack density [1/mm]')
        plt.xlabel('composite stress [MPa]')
        plt.show()

    plot()