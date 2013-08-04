'''
Created on 11 May 2013

the Interpolator class evaluates crack bridges along x and for a range of load and BCs
the Interpolator2 class evaluates crack bridges along x and for a range of load and for given fixed BCs

@author: Q
'''
from etsproxy.traits.api import HasTraits, Property, cached_property, \
    Instance, Array, Float, Int
import numpy as np
from spirrid.rv import RV
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import ContinuousFibers
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx_view import CompositeCrackBridgeView
from scipy.optimize import brentq
import time
from scipy.interpolate import griddata, interp2d
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

class Interpolator(HasTraits):

    CB_model = Instance(CompositeCrackBridge)
    n_w = Int
    n_BC = Int
    n_x = Int
    load_sigma_c_arr = Array
    length = Float

    CB_model_view = Property(Instance(CompositeCrackBridgeView), depends_on='CB_model')
    @cached_property
    def _get_CB_model_view(self):
        for i, reinf in enumerate(self.CB_model.reinforcement_lst):
            self.CB_model.reinforcement_lst[i].n_int = self.n_x
        return CompositeCrackBridgeView(model=self.CB_model)

    def max_sigma_w(self, Ll, Lr):
        self.CB_model_view.model.Ll = Ll
        self.CB_model_view.model.Lr = Lr
        max_sigma_c, max_w = self.CB_model_view.sigma_c_max
        if max_sigma_c < self.load_sigma_c_arr[-1]:
            return max_sigma_c, max_w
        else:
            def residuum(w):
                self.CB_model_view.model.w = float(w)
                sigma_c = self.CB_model_view.sigma_c
                return self.load_sigma_c_arr[-1] - sigma_c
            max_w = brentq(residuum, 0.0, max_w)
            return self.load_sigma_c_arr[-1], max_w

    BC_range = Property(depends_on = 'n_BC, CB_model')
    @cached_property
    def _get_BC_range(self):
        self.max_sigma_w(np.inf, np.inf)
        Lmax = min(self.CB_model_view.x_arr[-2], self.length)
        return np.linspace(1.0, Lmax, self.n_BC)

    result_values = Property(Array, depends_on='CB_model, load_sigma_c_arr, n_w, n_x, n_BC')
    @cached_property
    def _get_result_values(self):
        L_arr = self.BC_range
        Ll_arr = np.array([])
        Lr_arr = np.array([])
        x_arr = np.array([])
        sigma_c_arr = np.array([])
        mu_epsf_arr = np.array([])
        epsm_arr = np.array([])
        loops_tot = self.n_BC ** 2
        max_sigma_c_arr = np.zeros((self.n_BC, self.n_BC))
        for i, ll in enumerate(L_arr):
            for j, lr in enumerate(L_arr):
                if j >= i:
                    # find maximum
                    sigma_c_max, wmax = self.max_sigma_w(ll, lr)
                    max_sigma_c_arr[i, j] = max_sigma_c_arr[j, i] = sigma_c_max
                    w_arr = np.linspace(0.0, wmax, self.n_w)
                    mu_sigma_c, x, mu_epsf, epsm, ll_arr, lr_arr = self.CB_model_view.w_x_res(w_arr, ll, lr, self.length)
                    # store the particular result for BC ll and lr into the result array
                    Ll_arr = np.hstack((Ll_arr, ll_arr))
                    Lr_arr = np.hstack((Lr_arr, lr_arr))
                    x_arr = np.hstack((x_arr, x))
                    sigma_c_arr = np.hstack((sigma_c_arr, mu_sigma_c))
                    mu_epsf_arr = np.hstack((mu_epsf_arr, mu_epsf))
                    epsm_arr = np.hstack((epsm_arr, epsm))
                current_loop = i * len(L_arr) + j + 1
                print 'progress: %2.1f %%' % \
                (current_loop / float(loops_tot) * 100.)
        points = np.array([Ll_arr, Lr_arr, x_arr, sigma_c_arr])
        return points, mu_epsf_arr, epsm_arr, max_sigma_c_arr

    def interpolate(self, load, x, Ll, Lr, idx):
        strength = self.interpolate_max_sigma_c(Ll, Lr)
        if load > strength:
            raise ValueError('applied load', load, 'MPa higher then strength ', strength, 'MPa')
        else:
            L_l, L_r = self.get_L(Ll, Lr)
            L_l = self.BC_range[np.argwhere(L_l <= self.BC_range)[0]]
            L_r = self.BC_range[np.argwhere(L_r <= self.BC_range)[0]]
            
            maskBC = (self.result_values[0][0] == L_l) * (self.result_values[0][1] == L_r) == 1
            
            
            sigma_c_masked = self.result_values[0][-1][maskBC]
            sigma_c_low = np.max(sigma_c_masked[np.argwhere(sigma_c_masked < load)])
            sigma_c_high = np.min(sigma_c_masked[np.argwhere(sigma_c_masked > load)])
            mask_sigma_c_low = (self.result_values[0][-1] == sigma_c_low)
            mask_sigma_c_high = (self.result_values[0][-1] == sigma_c_high)
            line_low = MFnLineArray(xdata=self.result_values[0][2][maskBC * mask_sigma_c_low == 1],
                                    ydata=self.result_values[idx][maskBC * mask_sigma_c_low == 1])
            values_low = line_low.get_values(x) * (sigma_c_high - load) / (sigma_c_high - sigma_c_low)
            line_high = MFnLineArray(xdata=self.result_values[0][2][maskBC * mask_sigma_c_high == 1],
                                    ydata=self.result_values[idx][maskBC * mask_sigma_c_high == 1])
            values_high = line_high.get_values(x) * (load - sigma_c_low) / (sigma_c_high - sigma_c_low) 
            return values_low + values_high

    def interpolate_mu_epsf(self, load, x, Ll, Lr):
        return self.interpolate(load, x, Ll, Lr, 1)

    def interpolate_epsm(self, load, x, Ll, Lr):
        return self.interpolate(load, x, Ll, Lr, 2)

    def interpolate_max_sigma_c(self, Ll, Lr):
        L_l, L_r = self.get_L(Ll, Lr)
        L_l = self.BC_range[np.argwhere(L_l <= self.BC_range)[0]]
        L_r = self.BC_range[np.argwhere(L_r <= self.BC_range)[0]]
        interp = interp2d(self.BC_range, self.BC_range, self.result_values[3])
        return interp(L_l, L_r)

    def get_L(self, Ll, Lr):
        self.result_values
        BC_line = MFnLineArray(xdata=self.BC_range, ydata=self.BC_range, extrapolate='constant')
        return BC_line.get_values([Ll, Lr])

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    reinf = ContinuousFibers(r=0.00345,
                          tau=RV('weibull_min', loc=0.007, shape=1.22, scale=.04),
                          V_f=0.0103,
                          E_f=180e3,
                          xi=WeibullFibers(shape=5.0, sV0=0.0032),
                          n_int=200,
                          label='carbon')

    CB_model = CompositeCrackBridge(E_m=25e3, reinforcement_lst=[reinf])

    ir = Interpolator(CB_model=CB_model,
                             load_sigma_c_arr=np.linspace(0.0, 25., 50),
                             n_w=50,
                             n_BC=3,
                             n_x=100,
                             length=500.
                             )

    Ll = 2.6
    Lr = 3.
    x_arr = np.linspace(-Ll, Lr, 500)
    Ll = 5.
    Lr = 2.
    x_arr = np.linspace(-Ll, Lr, 200)
    sigma_c = np.linspace(1., 11., 7)
    for i, s in enumerate(sigma_c):
        plt.plot(x_arr, ir.interpolate_epsm(s, x_arr, Ll, Lr))
        plt.plot(x_arr, ir.interpolate_mu_epsf(s, x_arr, Ll, Lr))
    plt.show()
