
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
from scipy.interpolate import interp2d
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray


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
            self.CB_model_view.apply_load(self.load_sigma_c_arr[-1] - 1e-10)
            return self.load_sigma_c_arr[-1] - 1e-10, self.CB_model_view.model.w

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
        cb_epsm_interpolators_lst = [MFnLineArray(xdata=np.linspace(-1e5,1e5,5), ydata=np.zeros(5))]
        cb_epsf_interpolators_lst = [MFnLineArray(xdata=np.linspace(-1e5,1e5,5), ydata=np.zeros(5))]
        sigma_c_lst = [0.0]
        w_lst = [0.0]
        for w in w_arr:
            self.CB_model_view.model.w = w
            if self.CB_model_view.sigma_c > sigma_c_lst[-1]:
                w_lst.append(w)
                sigma_c_lst.append(self.CB_model_view.sigma_c)
                x_i = np.hstack((-self.length - 1e-1, self.CB_model_view.x_arr, self.length + 1e-1))
                epsm_i = np.hstack((self.CB_model_view.epsm_arr[0], self.CB_model_view.epsm_arr, self.CB_model_view.epsm_arr[-1]))
                epsf_i = np.hstack((self.CB_model_view.epsf_arr[0], self.CB_model_view.epsf_arr, self.CB_model_view.epsf_arr[-1]))
                cb_epsm_interpolators_lst.append(MFnLineArray(xdata=x_i, ydata=epsm_i))
                cb_epsf_interpolators_lst.append(MFnLineArray(xdata=x_i, ydata=epsf_i))   
        w_interpolator = MFnLineArray(xdata=np.array(sigma_c_lst), ydata=np.array(w_lst))
        return w_interpolator, [sigma_c_lst, cb_epsm_interpolators_lst], [sigma_c_lst, cb_epsf_interpolators_lst]


    interpolator_lists = Property(Array, depends_on='CB_model, load_sigma_c_arr, n_w, n_x, n_BC')
    @cached_property
    def _get_interpolator_lists(self):
            epsm_interpolators = np.zeros((self.n_BC, self.n_BC), dtype = np.object)
            epsf_interpolators = np.zeros((self.n_BC, self.n_BC), dtype = np.object)
            w_interpolators = np.zeros((self.n_BC, self.n_BC), dtype = np.object)
            loops_tot = self.n_BC ** 2
            max_sigma_c_arr = np.zeros((self.n_BC, self.n_BC))
            for i, ll in enumerate(self.BC_range):
                for j, lr in enumerate(self.BC_range):
                    if j >= i:
                        # find maximum
                        sigma_c_max, wmax = self.max_sigma_w(ll, lr)
                        max_sigma_c_arr[i, j] = max_sigma_c_arr[j, i] = sigma_c_max
                        w_arr0 = np.linspace(1e-10, wmax, self.n_w)
                        w_interpolator, epsm_interp_lst, epsf_interp_lst = self.w_x_res(w_arr0, ll, lr)
                        epsm_interpolators[i,j] = epsm_interpolators[j,i] = epsm_interp_lst
                        epsf_interpolators[i,j] = epsf_interpolators[j,i] = epsf_interp_lst
                        w_interpolators[i,j] = w_interpolators[j,i] = w_interpolator
                    current_loop = i * len(self.BC_range) + j + 1
                    print 'progress: %2.1f %%' % \
                    (current_loop / float(loops_tot) * 100.)
            interp_max_sigma_c = interp2d(self.BC_range, self.BC_range, max_sigma_c_arr, fill_value = None)
            return interp_max_sigma_c, epsm_interpolators, w_interpolators, epsf_interpolators

    def get_BC_idxs(self, Ll, Lr):
        if Ll >= self.BC_range[-1]:
            ll_idx_high = -1
            ll_idx_low = -1
        elif Ll <= self.BC_range[0]:
            ll_idx_high = 0
            ll_idx_low = 0
        else:
            ll_idx_high = np.argwhere(Ll <= self.BC_range)[0][0]
            ll_idx_low = np.argwhere(Ll >= self.BC_range)[-1][0]
        if Lr > self.BC_range[-1]:
            lr_idx_high = -1
            lr_idx_low = -1
        elif Lr <= self.BC_range[0]:
            lr_idx_high = 0
            lr_idx_low = 0
        else:
            lr_idx_high = np.argwhere(Lr <= self.BC_range)[0][0]
            lr_idx_low = np.argwhere(Lr >= self.BC_range)[-1][0]
        return ll_idx_high, lr_idx_high, ll_idx_low, lr_idx_low
        
    def interpolate_max_sigma_c(self, Ll, Lr):
        return self.interpolator_lists[0](Ll, Lr)
    
    def interpolate_epsm(self, Ll, Lr, sigma_c, x_arr):
        ll_idx_high, lr_idx_high, ll_idx_low, lr_idx_low = self.get_BC_idxs(Ll, Lr)
        epsm_interpolator_lst = self.interpolator_lists[1][ll_idx_high, lr_idx_high]
        sigc = np.array(epsm_interpolator_lst[0])
        if sigma_c > sigc[-1]:
            # applied stress is higher than crack bridge strength
            return np.repeat(np.nan, len(x_arr)) 
        else:
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

    def interpolate_epsf(self, Ll, Lr, sigma_c, x_arr):
        ll_idx_high, lr_idx_high, ll_idx_low, lr_idx_low = self.get_BC_idxs(Ll, Lr)
        epsf_interpolator_lst = self.interpolator_lists[3][ll_idx_high, lr_idx_high]
        sigc = np.array(epsf_interpolator_lst[0])
        if sigma_c > sigc[-1]:
            # applied stress is higher than crack bridge strength
            return np.repeat(np.nan, len(x_arr)) 
        else:
            sigc_high = np.argwhere(sigc > sigma_c)[0][0]
            sigc_low = np.argwhere(sigc < sigma_c)[-1][0]
            coeff_low = (sigc[sigc_high] - sigma_c) / (sigc[sigc_high] - sigc[sigc_low])
            coeff_high = (sigma_c - sigc[sigc_low]) / (sigc[sigc_high] - sigc[sigc_low])
            if Lr >= Ll:
                epsf = epsf_interpolator_lst[1][sigc_low].get_values(x_arr) * coeff_low + \
                       epsf_interpolator_lst[1][sigc_high].get_values(x_arr) * coeff_high
            else:
                epsf = epsf_interpolator_lst[1][sigc_low].get_values(-x_arr[::-1]) * coeff_low + \
                        epsf_interpolator_lst[1][sigc_high].get_values(-x_arr[::-1]) * coeff_high
                epsf = epsf[::-1]
            return epsf
        
    def interpolate_w(self, Ll, Lr, sigma_c):
        '''
        interpolation of w using the approach of interpolation on a 4-node rectangular finite element
        '''
        ll_idx_high, lr_idx_high, ll_idx_low, lr_idx_low = self.get_BC_idxs(Ll, Lr)
        ll_low, ll_high, lr_low, lr_high = self.BC_range[ll_idx_low], self.BC_range[ll_idx_high], self.BC_range[lr_idx_low], self.BC_range[lr_idx_high]
        # evaluating nodal values / note that w2 = w3
        w1 = self.interpolator_lists[2][ll_idx_low, lr_idx_low].get_values(np.array([sigma_c]))
        w2 = self.interpolator_lists[2][ll_idx_low, lr_idx_high].get_values(np.array([sigma_c]))
        w3 = self.interpolator_lists[2][ll_idx_high, lr_idx_low].get_values(np.array([sigma_c]))
        w4 = self.interpolator_lists[2][ll_idx_high, lr_idx_high].get_values(np.array([sigma_c]))
        nodal_values = [w1, w2, w3, w4]
        # shape functions
        if ll_idx_low == ll_idx_high:
            if lr_idx_low == lr_idx_high:
                # no interpolation
                return w1
            else:
                # 1D interpolation
                a = lr_high - lr_low
                N1 = lambda x: - (x - lr_high)/a
                N2 = lambda x: (x - lr_low)/a
                return w1 * N1(Lr) + w2*N2(Lr)
        else:
            if lr_idx_low == lr_idx_high:
                # 1D interpolation
                a = ll_high - ll_low
                N1 = lambda x: - (x - ll_high)/a
                N2 = lambda x: (x - ll_low)/a
                return w2 * N1(Ll) + w3*N2(Ll)
            else:
                #2D interpolation 
                ab = (ll_high-ll_low) * (lr_high - lr_low)
                N1 = lambda x,y: (x - ll_high) * (y - lr_high) / ab 
                N2 = lambda x,y: - (x - ll_low) * (y - lr_high) / ab
                N3 = lambda x,y: (x - ll_low) * (y - lr_low) / ab
                N4 = lambda x,y: - (x - ll_high) * (y - lr_low) / ab
                shape_functions = [N1, N2, N3, N4]
                # interpolate w
                w_interpolated = 0.0
                for i, Ni in enumerate(shape_functions):
                    w_interpolated += nodal_values[i] * Ni(Ll, Lr)
                return w_interpolated

if __name__ == '__main__':
    pass