'''
Created on 11 May 2013

the Interpolator class evaluates crack bridges along x and for a range of load and BCs
the Interpolator2 class evaluates crack bridges along x and for a range of load and for given fixed BCs

@author: Q
'''
from stats.spirrid import make_ogrid as orthogonalize
from etsproxy.traits.api import HasTraits, Property, cached_property, \
    Instance, Array, List, Float, Int
import numpy as np
from scipy import ndimage
import types
from etsproxy.mayavi import mlab as m
from spirrid.rv import RV
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import Reinforcement, ContinuousFibers
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx_view import CompositeCrackBridgeView
from scipy.optimize import minimize, fminbound, brentq
import time

def H(x):
    return x >= 0.0

def orthogonalize_filled(args):
    '''creates meshgrid up to third dimension
    given a list of 1D arrays and floats
    '''

    array_list = []
    array_args = []
    for arg in args:
        if isinstance(arg, np.ndarray):
            array_args.append(arg)

    if len(array_args) == 0:
        meshgrid = np.array([1.])
    elif len(array_args) == 1:
        meshgrid = [array_args[0]]
    elif len(array_args) == 2 or len(array_args) == 2:
        meshgrid = np.meshgrid(*array_args)
    else:
        raise NotImplementedError('''max number of arrays as
                                input is currently 3''')

    i = 0
    if len(meshgrid) == 0:
        meshgrid = np.array([1.])
    for arg in args:
        if isinstance(arg, np.ndarray):
            array_list.append(meshgrid[i])
            i += 1
        elif isinstance(arg, types.FloatType):
            array_list.append(np.ones_like(meshgrid[0]) * arg)
    return array_list


class NDIdxInterp(HasTraits):

    # nd array of values (measured, computed..) of
    # size orthogonalize(axes_values)
    data = Array

    # list of control input parameter values
    axes_values = List(Array(float))

    def __call__(self, *gcoords, **kw):
        '''kw: dictionary of values to interpolate for;
        len(kw) has to be equal the data dimension
        '''
        order = kw.get('order', 1)
        mode = kw.get('mode', 'nearest')

        # check if the number of dimensions to interpolate
        # in equals the number of given coordinates
        if len(self.axes_values) != len(gcoords):
            raise TypeError('''method takes {req} arguments
            ({given} given)'''.format(req=len(self.axes_values),
                                      given=len(gcoords)))
        icoords = self.get_icoords(gcoords)

        # create a meshgrid for the interpolation
        icoords = orthogonalize_filled(icoords)
        data = self.data
        # interpolate the value (linear)
        # a, b, c = [0.5, 0.5, 0.5], [0, 0, 0], [0, 1, 2]
        # icoords = [a, b, c]
        val = ndimage.map_coordinates(data, icoords, order=order, mode=mode)
        return val

    def get_icoords(self, gcoords):
        '''
        gcoords: values to be interpolated for
        this method transforms the global coords to "index" coords
        '''
        icoords = [np.interp(gcoord, axis_values, np.arange(len(axis_values)))
                   for gcoord, axis_values in zip(gcoords, self.axes_values)]
        return icoords


class Interpolator2(HasTraits):

    CB_model = Instance(CompositeCrackBridgeView)
    n_w = Int
    n_x = Int
    n_BC = Int
    load_sigma_c_arr = Array
    length = Float
    
    def max_sigma_w(self, Ll, Lr):
        self.CB_model.model.Ll = Ll
        self.CB_model.model.Lr = Lr
        def min_full_sigma_c(w):
            self.CB_model.model.w = float(w)
            sigma_c = self.CB_model.sigma_c
            return - sigma_c
        wmax_full = minimize(min_full_sigma_c, 0.0, options = dict(maxiter = 10)).x
        if self.CB_model.sigma_c < self.load_sigma_c_arr[-1]:
            return self.CB_model.sigma_c, wmax_full
        else:
            def residuum(w):
                self.CB_model.model.w = float(w)
                sigma_c = self.CB_model.sigma_c
                return self.load_sigma_c_arr[-1] - sigma_c
            result = brentq(residuum, 0.0, wmax_full)
            return self.CB_model.sigma_c, result
    
    BC_range = Property(depends_on = 'n_BC, CB_model')
    @cached_property
    def _get_BC_range(self):
        self.max_sigma_w(np.inf, np.inf)
        Lmax = min(self.CB_model.x_arr[-2], self.length)
        return np.linspace(1.0, Lmax, self.n_BC)
    
    x_arr = Property(depends_on = 'n_BC, CB_model, n_x')
    @cached_property
    def _get_x_arr(self):
        return np.linspace(-self.BC_range[-1], self.BC_range[-1], self.n_x) 

    def preinterpolate(self, sigma_c_w_x, sigma_c_cutoff, x_range):
        # values to create array grid
        axes_values = [sigma_c_cutoff, x_range]
        preinterp = NDIdxInterp(data=sigma_c_w_x, axes_values=axes_values)
        # values to interpolate for
        interp_coords = [self.load_sigma_c_arr, self.x_arr]
        return preinterp(*interp_coords, mode='constant')

    result_values = Property(Array)
    @cached_property
    def _get_result_values(self):
        L_arr = self.BC_range
        epsm_w_x = np.zeros((len(self.load_sigma_c_arr), self.n_x,
                           self.n_BC, self.n_BC))
        mu_epsf_w_x = np.zeros((len(self.load_sigma_c_arr), self.n_x,
                           self.n_BC, self.n_BC))
        mu_sigma_c_w = np.zeros((len(self.load_sigma_c_arr), self.n_x,
                           self.n_BC, self.n_BC))
        loops_tot = len(L_arr) ** 2
        for i, ll in enumerate(L_arr):
            for j, lr in enumerate(L_arr):
                if j >= i:
                    # adapt w range
                    sigma_c_max, wmax = self.max_sigma_w(ll, lr)
                    w_arr = np.linspace(0.0, wmax, self.n_w)
                    # evaluate the result (2D (w,x) SPIRRID with adapted ranges x and w
                    CB_epsm_w_x, CB_mu_epsf_w_x, sigma_c = self.CB_model.w_x_results(w_arr, self.x_arr)
                    eps_vars = orthogonalize([np.arange(len(w_arr))/5., np.arange(len(self.x_arr))])
                    m.surf(eps_vars[0], eps_vars[1], CB_epsm_w_x*500000)
                    m.surf(eps_vars[0], eps_vars[1], CB_mu_epsf_w_x*5000)
                    m.show()
                    # preinterpolate particular result for the given x and sigma ranges
                    epsm_w_x_preinterp = \
                    self.preinterpolate(CB_epsm_w_x, sigma_c, self.x_arr).T
                    mu_epsf_w_x_preinterp = \
                    self.preinterpolate(CB_mu_epsf_w_x, sigma_c, self.x_arr).T
                    mask = np.where(self.load_sigma_c_arr
                                    <= sigma_c_max, 1, np.NaN)[:, np.newaxis]
                    epsm_w_x_preinterp = epsm_w_x_preinterp * mask
                    mu_epsf_w_x_preinterp = mu_epsf_w_x_preinterp * mask
                    mu_sigma_c_w_preinterp = np.ones_like(epsm_w_x_preinterp) * self.load_sigma_c_arr[:, np.newaxis] * mask
#                     eps_vars = orthogonalize([np.arange(len(w_arr))/5., self.x_arr])
#                     m.surf(eps_vars[0], eps_vars[1], epsm_w_x_preinterp*500)
#                     m.surf(eps_vars[0], eps_vars[1], mu_epsf_w_x_preinterp*500)
#                     m.surf(eps_vars[0], eps_vars[1], epsm_w_x_interp*1000)
#                     m.show()
                    # store the particular result for BC ll and lr into the result array
                    epsm_w_x[:, :, i, j] = epsm_w_x_preinterp
                    epsm_w_x[:, :, j, i] = epsm_w_x_preinterp[:,::-1]
                    mu_epsf_w_x[:, :, i, j] = mu_epsf_w_x_preinterp
                    mu_epsf_w_x[:, :, j, i] = mu_epsf_w_x_preinterp[:,::-1]
                    mu_sigma_c_w[:, :, i, j] = mu_sigma_c_w_preinterp
                    mu_sigma_c_w[:, :, j, i] = mu_sigma_c_w_preinterp[:,::-1]
                    current_loop = i * len(L_arr) + j + 1
                    print 'progress: %2.1f %%' % \
                    (current_loop / float(loops_tot) * 100.)
        axes_values = [self.load_sigma_c_arr, self.x_arr, self.BC_range, self.BC_range]
        interp_epsm = NDIdxInterp(data=epsm_w_x, axes_values=axes_values)
        interp_epsf = NDIdxInterp(data=mu_epsf_w_x, axes_values=axes_values) 
        interp_sigmac = NDIdxInterp(data=mu_sigma_c_w, axes_values=axes_values) 
        return interp_epsm, interp_epsf, interp_sigmac

    interpolator_epsm = Property(depends_on = 'CB_model, load_sigma_c_max, load_n_sigma_c, n_w, n_x, n_BC')
    @cached_property
    def _get_interpolator_epsm(self):
        return self.result_values[0]
    
    interpolator_mu_epsf = Property(depends_on = 'CB_model, load_sigma_c_max, load_n_sigma_c, n_w, n_x, n_BC')
    @cached_property
    def _get_interpolator_mu_epsf(self):
        return self.result_values[1]
    
    interpolator_mu_sigma_c = Property(depends_on = 'CB_model, load_sigma_c_max, load_n_sigma_c, n_w, n_x, n_BC')
    @cached_property
    def _get_interpolator_mu_sigma_c(self):
        return self.result_values[2]


class Interpolator(HasTraits):

    CB_model = Instance(CompositeCrackBridgeView)
    n_w = Int
    n_BC = Int
    load_sigma_c_arr = Array
    length = Float
    
    n_x = Property(depends_on='CB_model')
    @cached_property
    def _get_n_x(self):
        return len(self.CB_model.x_arr) - 2
    
    def max_sigma_w(self, Ll, Lr):
        self.CB_model.model.Ll = Ll
        self.CB_model.model.Lr = Lr
        def min_full_sigma_c(w):
            self.CB_model.model.w = float(w)
            sigma_c = self.CB_model.sigma_c
            return - sigma_c
        wmax_full = minimize(min_full_sigma_c, 0.0, options = dict(maxiter = 10)).x
        if self.CB_model.sigma_c < self.load_sigma_c_arr[-1]:
            return self.CB_model.sigma_c, wmax_full
        else:
            def residuum(w):
                self.CB_model.model.w = float(w)
                sigma_c = self.CB_model.sigma_c
                return self.load_sigma_c_arr[-1] - sigma_c
            result = brentq(residuum, 0.0, wmax_full)
            return self.CB_model.sigma_c, result
    
    BC_range = Property(depends_on = 'n_BC, CB_model')
    @cached_property
    def _get_BC_range(self):
        self.max_sigma_w(np.inf, np.inf)
        Lmax = min(self.CB_model.x_arr[-2], self.length)
        return np.linspace(1.0, Lmax, self.n_BC)

    result_values = Property(Array)
    @cached_property
    def _get_result_values(self):
        L_arr = self.BC_range
        epsm_w_x = np.zeros((self.n_w, self.n_x,
                           self.n_BC, self.n_BC))
        mu_epsf_w_x = np.zeros((self.n_w, self.n_x,
                           self.n_BC, self.n_BC))
        mu_sigma_c_w = np.zeros((self.n_w, self.n_x,
                           self.n_BC, self.n_BC))
        x_ndim = np.zeros((self.n_w, self.n_x,
                           self.n_BC, self.n_BC))
        L_BC = np.zeros((self.n_w, self.n_x,
                           self.n_BC, self.n_BC))
        R_BC = np.zeros((self.n_w, self.n_x,
                           self.n_BC, self.n_BC))
        loops_tot = self.n_BC ** 2
        for i, ll in enumerate(L_arr):
            for j, lr in enumerate(L_arr):
                if j >= i:
                    # find maximum
                    sigma_c_max, wmax = self.max_sigma_w(ll, lr)
                    w_arr = np.linspace(0.0, wmax, self.n_w)
                    # 
                    epsm, mu_epsf, x, mu_sigma_c = self.CB_model.w_x_res(w_arr)
                    mask = np.where(mu_sigma_c <= sigma_c_max, 1, np.NaN)
                    epsm = epsm * mask
                    mu_epsf = mu_epsf * mask
                    # store the particular result for BC ll and lr into the result array
                    print mask
                    epsm_w_x[:, :, i, j] = epsm
                    epsm_w_x[:, :, j, i] = epsm[:,::-1]
                    mu_epsf_w_x[:, :, i, j] = mu_epsf
                    mu_epsf_w_x[:, :, j, i] = mu_epsf[:,::-1]
                    mu_sigma_c_w[:, :, i, j] = mu_sigma_c
                    mu_sigma_c_w[:, :, j, i] = mu_sigma_c[:,::-1]
                    x_ndim[:, :, i, j] = x
                    x_ndim[:, :, j, i] = x[:,::-1]
                    current_loop = i * len(L_arr) + j + 1
                    print 'progress: %2.1f %%' % \
                    (current_loop / float(loops_tot) * 100.)
        axes_values = [self.load_sigma_c_arr, self.x_arr, self.BC_range, self.BC_range]
        interp_epsm = NDIdxInterp(data=epsm_w_x, axes_values=axes_values)
        interp_epsf = NDIdxInterp(data=mu_epsf_w_x, axes_values=axes_values) 
        interp_sigmac = NDIdxInterp(data=mu_sigma_c_w, axes_values=axes_values) 
        return interp_epsm, interp_epsf, interp_sigmac

    interpolator_epsm = Property(depends_on = 'CB_model, load_sigma_c_max, load_n_sigma_c, n_w, n_x, n_BC')
    @cached_property
    def _get_interpolator_epsm(self):
        return self.result_values[0]
    
    interpolator_mu_epsf = Property(depends_on = 'CB_model, load_sigma_c_max, load_n_sigma_c, n_w, n_x, n_BC')
    @cached_property
    def _get_interpolator_mu_epsf(self):
        return self.result_values[1]
    
    interpolator_mu_sigma_c = Property(depends_on = 'CB_model, load_sigma_c_max, load_n_sigma_c, n_w, n_x, n_BC')
    @cached_property
    def _get_interpolator_mu_sigma_c(self):
        return self.result_values[2]

if __name__ == '__main__':
    from stats.spirrid import make_ogrid as orthogonalize
    from matplotlib import pyplot as plt

    reinf = ContinuousFibers(r=0.00345,#RV('uniform', loc=0.002, scale=0.002),
                          tau=RV('uniform', loc=0.02, scale=20.),
                          V_f=0.15,
                          E_f=200e3,
                          xi=WeibullFibers(shape=5., sV0=0.00618983207723),
                          n_int=50,
                          label='carbon')

    model = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf],
                                 Ll=100.,
                                 Lr=100.)

    ccb_post = CompositeCrackBridgeView(model=model)
    
    ir = Interpolator(CB_model = ccb_post,
                             load_sigma_c_arr = np.linspace(0.0, 25., 50),
                             n_w = 30,
                             n_BC = 3
                             )

    def plot():
        
        w_arr = np.linspace(0.0, 0.1, 20)
        ll = 2.0
        lr = 5.0
        x = np.linspace(-ll, lr, 51)
        #resm, resf = ccb_post.w_x_results(w_arr, x)
        #eps_vars = orthogonalize([w_arr * 100., x])
        #m.surf(eps_vars[0], eps_vars[1], resm * 500.)
        #m.surf(eps_vars[0], eps_vars[1], resf * 500.)
        #m.show()
        sigma = ir.load_sigma_c_arr

        eps_vars = orthogonalize([np.linspace(0.0, 0.0254229525299, len(sigma)) * 100., x])
        # mu_q_nisp = nisp(P, x, Ll, Lr)[0]
        mu_epsm = ir.interpolator_epsm(sigma, x, 100., 100.)
        mu_epsf = ir.interpolator_mu_epsf(sigma, x, 100., 100.)
        sig_c = ir.interpolator_mu_sigma_c(sigma, x, 100., 100.)
        #mu_q_isp2 = ir(sigma, x, Ll, Lr)

#        plt.plot(np.arange(len(sigma)), sigma/0.0103)
#        plt.plot(np.arange(len(sigma)), np.max(mu_q_isp,axis = 0))
#        plt.show()
        # n_mu_q_arr = mu_q_nisp / np.max(np.fabs(mu_q_nisp))
        m.surf(eps_vars[0], eps_vars[1], mu_epsm * 500.)
        m.surf(eps_vars[0], eps_vars[1], mu_epsf * 500.)
        m.surf(eps_vars[0], eps_vars[1], sig_c /(25e3*0.85 + 200e3*0.15)*500)
        m.show()

    plot()
