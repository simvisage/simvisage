
'''
Created on 11 May 2013

the Interpolator class evaluates crack bridges along x and for a range of load and BCs
the Interpolator2 class evaluates crack bridges along x and for a range of load and for given fixed BCs

@author: Q
'''
from etsproxy.traits.api import HasTraits, Property, cached_property, \
    Instance, Array, Float, Int
import numpy as np
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx_view import CompositeCrackBridgeView
from scipy.interpolate import interp2d, griddata
from matplotlib import pyplot as plt
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
import time as t


class RepresentativeCB(HasTraits):
    CB_model = Instance(CompositeCrackBridge)
    n_w = Int
    n_BC = Int
    load_sigma_c_arr = Array
    length = Float
    CB_model_view = Property(Instance(CompositeCrackBridgeView), depends_on='CB_model')
    @cached_property
    def _get_CB_model_view(self):
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

    def w_x_res(self, w_arr, ll, lr):
        self.CB_model_view.model.Ll = ll
        self.CB_model_view.model.Lr = lr
        part_epsm_interp_lst = [MFnLineArray(xdata=np.linspace(-1e10, 1e10, 10), ydata=np.zeros(10))]
        sigma_c = [0.0]
        for w in w_arr:
            self.CB_model_view.model.w = w
            if self.CB_model_view.sigma_c > sigma_c[-1]:
                sigma_c.append(self.CB_model_view.sigma_c)
                x_i = np.hstack((-self.length - 1e-1, self.CB_model_view.x_arr, self.length + 1e-1))
                epsm_i = np.hstack((self.CB_model_view.epsm_arr[0], self.CB_model_view.epsm_arr, self.CB_model_view.epsm_arr[-1]))
                part_epsm_interp_lst.append(MFnLineArray(xdata=x_i, ydata=epsm_i))
            else:
                break                
        return np.array(sigma_c), part_epsm_interp_lst

    interpolator_lists = Property(Array, depends_on='CB_model, load_sigma_c_arr, n_w, n_x, n_BC')
    @cached_property
    def _get_interpolator_lists(self):
            epsm_interpolators = np.zeros((self.n_BC, self.n_BC), dtype = np.object)
            w_interpolators = np.zeros((self.n_BC, self.n_BC), dtype = np.object)
            loops_tot = self.n_BC ** 2
            max_sigma_c_arr = np.zeros((self.n_BC, self.n_BC))
            for i, ll in enumerate(self.BC_range):
                for j, lr in enumerate(self.BC_range):
                    if j >= i:
                        # find maximum
                        sigma_c_max, wmax = self.max_sigma_w(ll, lr)
                        max_sigma_c_arr[i, j] = max_sigma_c_arr[j, i] = sigma_c_max
    
                        w_arr = np.linspace(0.0, wmax, self.n_w)
                        mu_sigma_c, part_epsm_interp_lst = self.w_x_res(w_arr, ll, lr)
                        epsm_interp_lst = [mu_sigma_c, part_epsm_interp_lst]
                        epsm_interpolators[i,j] = epsm_interpolators[j,i] = epsm_interp_lst
                        w_interp_lst = [mu_sigma_c, w_arr]
                        w_interpolators[i,j] = w_interpolators[j,i] = w_interp_lst
                    current_loop = i * len(self.BC_range) + j + 1
                    print 'progress: %2.1f %%' % \
                    (current_loop / float(loops_tot) * 100.)
            interp_max_sigma_c = interp2d(self.BC_range, self.BC_range, max_sigma_c_arr, fill_value = None)
            return interp_max_sigma_c, epsm_interpolators, w_interpolators

    def get_BC_idxs(self, Ll, Lr):
        if Ll > self.BC_range[-1]:
            ll_idx = -1
        else: 
            ll_idx = np.argwhere(Ll <= self.BC_range)[0][0]
        if Lr > self.BC_range[-1]:
            lr_idx = -1
        else: 
            lr_idx = np.argwhere(Lr <= self.BC_range)[0][0]
        return ll_idx, lr_idx
        
    def interpolate_max_sigma_c(self, Ll, Lr):
        return self.interpolator_lists[0](Ll, Lr)
    
    def interpolate_epsm(self, Ll, Lr, sigma_c, x_arr):
        ll_idx, lr_idx = self.get_BC_idxs(Ll, Lr)
        epsm_interpolator_lst = self.interpolator_lists[1][ll_idx, lr_idx]
        sigc = epsm_interpolator_lst[0]
        sigc_high = np.argwhere(sigc > sigma_c)[0][0]
        sigc_low = np.argwhere(sigc < sigma_c)[-1][0]
        coeff_low = (sigc[sigc_high] - sigma_c) / (sigc[sigc_high] - sigc[sigc_low])
        coeff_high = (sigma_c - sigc[sigc_low]) / (sigc[sigc_high] - sigc[sigc_low])
        if Lr >= Ll:
            epsm = epsm_interpolator_lst[1][sigc_low].get_values(x_arr) * coeff_low + \
                   epsm_interpolator_lst[1][sigc_high].get_values(x_arr) * coeff_high
        else:
            epsm = epsm_interpolator_lst[1][sigc_low].get_values(-x_arr[::-1]) * coeff_low + \
                    epsm_interpolator_lst[1][sigc_high].get_values(-x_arr[::-1]) * coeff_high
            epsm = epsm[::-1]
        return epsm
        
    def interpolate_w(self, Ll, Lr, sigma_c):
        ll_idx, lr_idx = self.get_BC_idxs(Ll, Lr)
        w_interpolator_data = self.interpolator_lists[2][ll_idx, lr_idx]
        w_interpolator = MFnLineArray(xdata=w_interpolator_data[0][1:], ydata=w_interpolator_data[1])
        return w_interpolator.get_values(np.array([sigma_c]))

if __name__ == '__main__':
    pass