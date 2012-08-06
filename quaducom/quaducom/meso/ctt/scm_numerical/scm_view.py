'''
Created on Jul 26, 2012

@author: rostar
'''

from etsproxy.traits.api import \
    Instance, Array, List, cached_property, Property
from etsproxy.traits.ui.api import ModelView
from stats.spirrid.spirrid import FunctionRandomization
from stats.spirrid.rv import RV
from stats.misc.random_field.random_field_1D import RandomField
import numpy as np
from matplotlib import pyplot as plt
from quaducom.meso.ctt.scm_numerical.scm_model import SCM

class SCMView(ModelView):
    
    model = Instance(SCM)

    def crack_widths(self, sigma_c):
        # find the index of the nearest value in the load range
        idx = np.abs(self.model.load_sigma_c - sigma_c).argmin()
        # evaluate the relative strain e_rel between fibers and matrix for the given load
        e_rel = self.eps_f[idx,:] - self.eps_m[idx,:]
        # pick the cracks that emerged at the given load
        cb_load = []
        for cb in self.model.cb_list:
            if cb.crack_load_sigma_c <= sigma_c:
                cb_load.append(cb)
        if len(cb_load) != 0:
            # find the symmetry points between cracks as the 0 element of their x range
            indexes = []
            for cb in cb_load:
                indexes.append(np.where(cb.position + cb.x[0] == self.model.x_arr)[0])
            # add the index of the last point
            indexes.append(self.model.nx-1)
            # list of crack widths to be filled in a loop with integrated e_rel
            crack_widths = [np.trapz(e_rel[idx:indexes[i+1]], self.model.x_arr[idx:indexes[i+1]])
                            for i, idx in enumerate(indexes[:-1])]  
            return np.array(crack_widths)
        else:
            return 0.0
    
    eval_w = Property(List, depends_on = 'model')
    @cached_property
    def _get_eval_w(self):
        return [self.crack_widths(load) for load in self.model.load_sigma_c]

    w_mean = Property(Array, depends_on = 'model')
    @cached_property
    def _get_w_mean(self):
        return np.array([np.mean(w) for w in self.eval_w])

    w_median = Property(Array, depends_on = 'model')
    @cached_property
    def _get_w_median(self):
        return np.array([np.median(w) for w in self.eval_w])
    
    w_stdev = Property(Array, depends_on = 'model')
    @cached_property
    def _get_w_stdev(self):
        return np.array([np.std(w) for w in self.eval_w])
    
    w_max = Property(Array, depends_on = 'model')
    @cached_property
    def _get_w_max(self):
        return np.array([np.max(w) for w in self.eval_w])
    
    x_area = Property(depends_on = 'model.')
    def _get_x_area(self):
        return  np.ones_like(self.model.load_sigma_c)[:, np.newaxis] * self.model.x_arr[np.newaxis, :]

    eps_m = Property(Array, depends_on = 'model.')
    @cached_property
    def _get_eps_m(self):
        return self.model.sigma_m_x / self.model.cb_randomization.tvars['E_m']

    sigma_f = Property(Array, depends_on = 'model.')
    @cached_property
    def _get_sigma_f(self):
        Vf = self.model.cb_randomization.tvars['V_f']
        return (self.model.load_sigma_c[:, np.newaxis] - self.model.sigma_m_x * (1.-Vf))/Vf

    eps_f = Property(Array, depends_on = 'model.')
    @cached_property
    def _get_eps_f(self):
        return self.sigma_f / self.model.cb_randomization.tvars['E_f']

    eps_sigma = Property(depends_on = 'model.')
    @cached_property
    def _get_eps_sigma(self):
        eps = np.trapz(self.eps_f, self.x_area, axis = 1) / self.model.length
        if np.sum(np.isnan(eps))==0:
            sigma = self.model.load_sigma_c
        else:
            idx = np.sum(np.isnan(eps) == False)
            eps = eps[:idx]
            eps[-1] = eps[-2]
            sigma = self.model.load_sigma_c[:idx]
            sigma[-1] = 0.0
        return eps, sigma
        
if __name__ == '__main__':
    from quaducom.micro.resp_func.cb_emtrx_clamped_fiber_stress import \
    CBEMClampedFiberStressSP
    from quaducom.micro.resp_func.cb_emtrx_clamped_fiber_stress_residual import \
    CBEMClampedFiberStressResidualSP
    
    # filaments
    r = 0.00345
    Vf = 0.0103
    tau = 0.1#RV('uniform', loc = 0.02, scale = .2) # 0.5
    Ef = 200e3
    Em = 25e3
    l = RV( 'uniform', scale = 20., loc = 2. )
    theta = 0.0
    xi = RV( 'weibull_min', scale = 0.02, shape = 5 ) # 0.017
    phi = 1.
    Pf = RV( 'uniform', scale = 1., loc = 0. )
    m = 5.0
    s0 = 10.02

    length = 2000.
    nx = 2000
    random_field = RandomField(seed = True,
                               lacor = 4.,
                                xgrid = np.linspace(0., length, 600),
                                nsim = 1,
                                loc = .0,
                                shape = 7.5,
                                scale = 7.8,
                                non_negative_check = True,
                                distribution = 'Weibull'
                               )
#
#    rf = CBEMClampedFiberStressSP()
#    rand = FunctionRandomization(   q = rf,
#                                    tvars = dict(tau = tau,
#                                                 l = l,
#                                                 E_f = Ef,
#                                                 theta = theta,
#                                                 xi = xi,
#                                                 phi = phi,
#                                                 E_m = Em,
#                                                 r = r,
#                                                 V_f = Vf
#                                                 ),
#                                    n_int = 20
#                                    )

    rf = CBEMClampedFiberStressResidualSP()
    rand = FunctionRandomization(   q = rf,
                                    tvars = dict(tau = tau,
                                                 l = l,
                                                 E_f = Ef,
                                                 theta = theta,
                                                 Pf = Pf,
                                                 phi = phi,
                                                 E_m = Em,
                                                 r = r,
                                                 V_f = Vf,
                                                 m = m,
                                                 s0 = s0
                                                 ),
                                    n_int = 20
                                    )

    scm = SCM(length = length,
              nx = nx,
              random_field = random_field,
              cb_randomization = rand,
              cb_type = 'mean',
              load_sigma_c_min = 2.1,
              load_sigma_c_max = 20.,
              load_n_sigma_c = 80,
              n_w = 60,
              n_x = 101,
              n_BC = 10
              )
    
    scm_view = SCMView(model = scm)
    scm_view.model.evaluate()

    def plot():
        eps, sigma = scm_view.eps_sigma
        plt.figure()
        plt.plot(eps, sigma, color = 'black', lw = 2, label = 'model')
        plt.legend(loc = 'best')
        plt.xlabel('composite strain [-]')
        plt.ylabel('composite stress [MPa]')
        plt.figure()
        plt.hist(scm_view.crack_widths(20.), bins = 20, label = 'load = 20 MPa')
        plt.hist(scm_view.crack_widths(15.), bins = 20, label = 'load = 15 MPa')
        plt.hist(scm_view.crack_widths(10.), bins = 20, label = 'load = 10 MPa')
        plt.legend(loc = 'best')
        plt.figure()
        plt.plot(scm.load_sigma_c, scm_view.w_mean, color = 'green', lw = 2, label = 'mean crack width')
        plt.plot(scm.load_sigma_c, scm_view.w_median, color = 'blue', lw = 2, label = 'median crack width')
        plt.plot(scm.load_sigma_c, scm_view.w_mean + scm_view.w_stdev, color = 'black', label = 'stdev')
        plt.plot(scm.load_sigma_c, scm_view.w_mean - scm_view.w_stdev, color = 'black')
        plt.plot(scm.load_sigma_c, scm_view.w_max, ls = 'dashed', color = 'red', label = 'max crack width')
        plt.legend(loc = 'best')
        plt.show()
    plot()

