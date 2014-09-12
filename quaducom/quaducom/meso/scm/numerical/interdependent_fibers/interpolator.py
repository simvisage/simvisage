'''
Created on 11 May 2013

the Interpolator class evaluates crack bridges along x and for a range of load and BCs
the Interpolator2 class evaluates crack bridges along x and for a range of load and for given fixed BCs

@author: Q
'''
from etsproxy.traits.api import HasTraits, Property, cached_property, \
    Instance, Array, Float, Int
import numpy as np
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import Reinforcement, ContinuousFibers
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx_view import CompositeCrackBridgeView
from scipy.interpolate import interp2d
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
import pickle
import os

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
            self.CB_model_view.apply_load(self.load_sigma_c_arr[-1])
            return self.load_sigma_c_arr[-1], self.CB_model_view.model.w

    BC_range = Property(depends_on='n_BC, CB_model')
    @cached_property
    def _get_BC_range(self):
        self.max_sigma_w(1e5, 1e5)
        Lmax = min(self.CB_model_view.x_arr[-2], self.length)
        bc_range = np.logspace(np.log10(1.0), np.log10(Lmax), self.n_BC)
        return bc_range

    def w_x_res(self, w_arr, ll, lr, maxBC):
        self.CB_model_view.model.Ll = ll
        self.CB_model_view.model.Lr = lr
        epsm = np.array([0.0])
        mu_epsf = np.array([0.0])
        x = np.array([0.0])
        sigma_c = np.array([0.0])
        ll_arr = np.array([0.0])
        lr_arr = np.array([0.0])
        for w in w_arr:
            self.CB_model_view.model.w = w
            if self.CB_model_view.sigma_c > sigma_c[-1]:
                sigma_c = np.hstack((sigma_c, np.ones(len(self.CB_model_view.x_arr) + 2) * self.CB_model_view.sigma_c))
                ll_arr = np.hstack((ll_arr, ll * np.ones(len(self.CB_model_view.x_arr) + 2)))
                lr_arr = np.hstack((lr_arr, lr * np.ones(len(self.CB_model_view.x_arr) + 2)))
                epsm = np.hstack((epsm, self.CB_model_view.epsm_arr[0], self.CB_model_view.epsm_arr, self.CB_model_view.epsm_arr[-1]))
                mu_epsf = np.hstack((mu_epsf, self.CB_model_view.mu_epsf_arr[0], self.CB_model_view.mu_epsf_arr, self.CB_model_view.mu_epsf_arr[-1]))
                x = np.hstack((x, -maxBC - 1e-10, self.CB_model_view.x_arr, maxBC + 1e-10))
                if ll != lr:
                    ll_arr = np.hstack((ll_arr, lr * np.ones(len(self.CB_model_view.x_arr) + 2)))
                    lr_arr = np.hstack((lr_arr, ll * np.ones(len(self.CB_model_view.x_arr) + 2)))
                    epsm = np.hstack((epsm, self.CB_model_view.epsm_arr[-1], self.CB_model_view.epsm_arr[::-1], self.CB_model_view.epsm_arr[0]))
                    mu_epsf = np.hstack((mu_epsf, self.CB_model_view.mu_epsf_arr[-1], self.CB_model_view.mu_epsf_arr[::-1], self.CB_model_view.mu_epsf_arr[0]))
                    x = np.hstack((x, -maxBC - 1e-10, -self.CB_model_view.x_arr[::-1], maxBC + 1e-10))
                    sigma_c = np.hstack((sigma_c, np.ones(len(self.CB_model_view.x_arr) + 2) * self.CB_model_view.sigma_c))
        return sigma_c, x, mu_epsf, epsm, ll_arr, lr_arr

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
                    mu_sigma_c, x, mu_epsf, epsm, ll_arr, lr_arr = self.w_x_res(w_arr, ll, lr, self.length)
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
        interp_arr_points = open('interp_arr_points.pkl', 'wb')
        pickle.dump([points, max_sigma_c_arr], interp_arr_points, -1)
        interp_arr_points.close()

        interp_arr_mu_epsf_arr = open('interp_arr_mu_epsf_arr.pkl', 'wb')
        pickle.dump([mu_epsf_arr], interp_arr_mu_epsf_arr, -1)
        interp_arr_mu_epsf_arr.close()

        interp_arr_epsm_arr = open('interp_arr_epsm_arr.pkl', 'wb')
        pickle.dump([epsm_arr], interp_arr_epsm_arr, -1)
        interp_arr_epsm_arr.close()
        os.chdir(os.pardir)
        return [points, mu_epsf_arr, epsm_arr, max_sigma_c_arr]

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

    def get_strain_profiles(self, Ll, Lr):
        L_l, L_r = self.get_L(Ll, Lr)
        L_l = self.BC_range[np.argwhere(L_l <= self.BC_range)[0]]
        L_r = self.BC_range[np.argwhere(L_r <= self.BC_range)[0]]
        maskBC = (self.result_values[0][0] == L_l) * (self.result_values[0][1] == L_r) == 1

        points = np.array([self.result_values[0][2][maskBC],
                           self.result_values[0][3][maskBC]])
        mu_epsf_arr = self.result_values[1][maskBC]
        epsm_arr = self.result_values[2][maskBC]
        return points, mu_epsf_arr, epsm_arr

if __name__ == '__main__':
    pass
