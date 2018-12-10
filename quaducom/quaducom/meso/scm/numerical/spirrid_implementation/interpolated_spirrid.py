'''
Created on Aug 17, 2011

@author: rostar
'''

from etsproxy.traits.api import HasTraits, Property, cached_property, \
    Instance, Array, List, Float, Int, Dict
from spirrid import SPIRRID
from spirrid.rv import RV
import numpy as np
from scipy import ndimage
import types
import copy


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
        elif isinstance(arg, float):
            array_list.append(np.ones_like(meshgrid[0]) * arg)
    return array_list


def find_closest_higher(array, scalar):
    idx = np.argwhere((array - scalar) > 0.0).flat[0]
    return array.flat[idx]


def arg_find_closest_higher(array, scalar):
    idx = np.argwhere((array - scalar) > 0.0).flat[0]
    return idx


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


class RangeAdaption(HasTraits):
    ''' incorporates methods for defining the ranges of BCs, w and x'''

    spirrid = Instance(SPIRRID)
    load_sigma_c_max = Float
    load_n_sigma_c = Int
    n_w = Int
    n_x = Int
    n_BC = Int

    tvar_dict = Property(Dict, depends_on='spirrid.theta_vars')

    @cached_property
    def _get_tvar_dict(self):
        tvar_dict = copy.deepcopy(self.spirrid.theta_vars)
        for key, value in zip(list(tvar_dict.keys()), list(tvar_dict.values())):
            if isinstance(value, RV):
                tvar_dict[key] = value._distr.mean
        return tvar_dict

    load_sigma_c = Property(Array,
                    depends_on='load_n_sigma_c, load_sigma_c_max')

    @cached_property
    def _get_load_sigma_c(self):
        return np.linspace(1e-6, self.load_sigma_c_max, self.load_n_sigma_c)

    load_sigma_f = Property(Array,
                    depends_on='load_n_sigma_c, load_sigma_c_max')

    @cached_property
    def _get_load_sigma_f(self):
        Vf = self.spirrid.theta_vars['V_f']
        return self.load_sigma_c / Vf

    delta0 = Property(Float, depends_on='spirrid.theta_vars')

    @cached_property
    def _get_delta0(self):
        ''' length of the ascending part of the stress
        profile for mean material parameters '''
        r = self.tvar_dict['r']
        V_f = self.tvar_dict['V_f']
        E_m = self.tvar_dict['E_m']
        E_f = self.tvar_dict['E_f']
        tau = self.tvar_dict['tau']
        delta0 = 0.5 * self.load_n_sigma_c / V_f / tau \
        / (V_f * E_f + (1 - V_f) * E_m) * r * (1 - V_f) * E_m
        return delta0

    def w_guess(self, ll, lr):
        l = self.tvar_dict['l']
        theta = self.tvar_dict['theta']
        E_f = self.tvar_dict['E_f']
        l = l * (1. + theta)
        epsmax = self.load_sigma_f[-1] / E_f
        if ll < l / 2.:
            left = ll * epsmax
            l_left = ll
        elif ll > l / 2. and ll - l / 2. < self.delta0:
            left = l / 2. * epsmax + (ll - (l / 2.)) * epsmax / 2.
            l_left = l / 2.
        elif ll > l / 2. + self.delta0:
            left = (l + self.delta0) / 2. * epsmax
            l_left = l / 2.
        if lr < l / 2.:
            right = lr * epsmax
            l_right = lr
        elif lr > l / 2. and lr - l / 2. < self.delta0:
            right = l / 2. * epsmax + (lr - (l / 2.)) * epsmax / 2.
            l_right = l / 2.
        elif lr > l / 2. + self.delta0:
            right = (l + self.delta0) / 2. * epsmax
            l_right = l / 2.
        wmax_init = left + right + (l_left + l_right) * theta
        return np.linspace(0.0, wmax_init, self.n_w)

    x_init = Property(Array)

    @cached_property
    def _get_x_init(self):
        l = self.tvar_dict['l']
        theta = self.tvar_dict['theta']
        x_init = np.linspace(-self.delta0 - l * (1 + theta) / 2.,
                             self.delta0 + l * (1 + theta) / 2., self.n_x)
        x_opt = self.adapt_x_range(self.w_opt_init[-1], x_init, 1000., 1000.)
        return x_opt

    w_opt_init = Property(Array)

    @cached_property
    def _get_w_opt_init(self):
        return self.adapt_w_range(1000., 1000.)[0]

    BC_range = Property(Array)

    @cached_property
    def _get_BC_range(self):
        return np.linspace(self.x_init[-1] \
                / self.n_BC, self.x_init[-1], self.n_BC)

    def adapt_w_range(self, ll, lr):
        ''' method for adapting the w range for given BC '''
        # three cases of the peak response are distinguished:
        # 1) response includes a higher value than the applied load
        # 2) response peak is lower and lies within the given w range
        # 3) response peak is lower and lies beyond the given w range
        w_guess = self.w_guess(ll, lr)
        self.spirrid.eps_vars = dict(w=w_guess)
        self.spirrid.theta_vars['x'] = 0.0
        self.spirrid.theta_vars['Ll'] = ll
        self.spirrid.theta_vars['Lr'] = lr
        sigma_f = self.spirrid.mu_q_arr
        # finding appropriate range for w
        load_sigma_f = self.load_sigma_f
        # is the evaluated stress higher than the applied stress? case 1)
        # if so, the CB will not break during the loading process
        if np.max(sigma_f) >= np.max(load_sigma_f):
            w_opt, sigma_opt = self.adapt_w_higher_stress(sigma_f, load_sigma_f)
        # the evaluated stress has a peak lower
        # than the maximum of the applied stress
        else:
            w_opt, sigma_opt = self.adapt_w_lower_stress(sigma_f, load_sigma_f)
        # cut off the part after peak value and invert w and sigma_f
        idxmax = np.argmax(sigma_opt)
        sigma_f_cutoff = sigma_opt[:idxmax + 1]
        w_cutoff = w_opt[:idxmax + 1]
        del self.spirrid.theta_vars['x']
        del self.spirrid.theta_vars['Ll']
        del self.spirrid.theta_vars['Lr']
        return w_cutoff, sigma_f_cutoff

    def adapt_w_higher_stress(self, sigma_f, load_sigma_f):
        # case 1) stretch the w range
        count = 0
        while np.max(sigma_f) / np.max(load_sigma_f) > 1.1 or \
            float(np.argmax(sigma_f)) / float(len(sigma_f) - 1) < 0.9:
            # find the first nonzero value
            min_idx = np.where((sigma_f == sigma_f[0]) == False)[0][0] - 1
            wmin = self.spirrid.eps_vars['w'][min_idx]
            # find the closest higher value to the max applied stress
            idx_closest = arg_find_closest_higher(sigma_f, np.max(load_sigma_f))
            wmax = self.spirrid.eps_vars['w'][idx_closest]
            # adapt the w range and evaluate spirrid with the adapted range
            self.spirrid.eps_vars = dict(w=np.linspace(wmin, wmax, self.n_w))
            sigma_f = self.spirrid.mu_q_arr
            count += 1
            if count > 3:
                raise ValueError('''got stuck in a loop adapting w
                - try to change the w range''')
        return self.spirrid.eps_vars['w'], sigma_f

    def adapt_w_lower_stress(self, sigma_f, load_sigma_f):
        # the peak is within the w range - stretch the w range
        if np.argmax(sigma_f) != len(sigma_f) - 1:
            count = 0
            while np.argmax(sigma_f) / float(len(sigma_f) - 1) < 0.9:
                # find the first nonzero value
                min_idx = np.where((sigma_f == sigma_f[0]) == False)[0][0] - 1
                wmin = self.spirrid.eps_vars['w'][min_idx]
                wmax = self.spirrid.eps_vars['w'][np.argmax(sigma_f)] * 1.05
                self.spirrid.eps_vars = dict(w=np.linspace(wmin, wmax, self.n_w))
                sigma_f = self.spirrid.mu_q_arr
                count += 1
                if count > 3:
                    raise ValueError('''got stuck in a loop adapting w
                    - try to change the w range''')
            return self.spirrid.eps_vars['w'], sigma_f
        # the peak is beyond the w range
        else:
            # stretch the w range until case 1) or 2) is attained
            count = 0
            while np.argmax(sigma_f) == len(sigma_f) - 1 and \
             np.max(sigma_f) < np.max(load_sigma_f):
                factor = np.max(load_sigma_f) / np.max(sigma_f)
                wmax = self.spirrid.eps_vars['w'][-1] * factor * 1.2
                self.spirrid.eps_vars['w'] = np.linspace(0.0, wmax, self.n_w)
                sigma_f = self.spirrid.mu_q_arr
                count += 1
                if count > 3:
                    raise ValueError('''got stuck in a loop adapting w
                    - try to change the w range''')
            # case 1)
            if np.argmax(sigma_f) == len(sigma_f) - 1 and \
            np.max(sigma_f) > np.max(load_sigma_f):
                w_opt, sigma_opt = \
                self.adapt_w_higher_stress(sigma_f, load_sigma_f)
            elif np.argmax(sigma_f) != len(sigma_f) - 1 and \
            np.max(sigma_f) > np.max(load_sigma_f):
                w_opt, sigma_opt = \
                self.adapt_w_higher_stress(sigma_f, load_sigma_f)
            # case 2)
            elif np.argmax(sigma_f) != len(sigma_f) - 1 and \
            np.max(sigma_f) < np.max(load_sigma_f):
                w_opt, sigma_opt = \
                self.adapt_w_lower_stress(sigma_f, load_sigma_f)
            return w_opt, sigma_opt

    def adapt_x_range(self, wmax, x_init, ll, lr):
        ''' adapts the x range for a given crack opening '''
        self.spirrid.theta_vars['w'] = wmax
        self.spirrid.theta_vars['Ll'] = ll
        self.spirrid.theta_vars['Lr'] = lr
        self.spirrid.eps_vars = dict(x=x_init)
        qx = self.spirrid.mu_q_arr
        x = x_init
        while np.min(qx) == np.max(qx) or \
        len(np.where(np.min(qx) == qx)[0]) < 3:
            # TODO smarter adaption of the x range
            x = np.linspace(x[0] - 0.5 * self.delta0,
                            x[-1] + 0.5 * self.delta0,
                            self.n_x)
            self.spirrid.eps_vars['x'] = x
            qx = self.spirrid.mu_q_arr
        # look for the index in q_x_arr where the constant
        # part changes to non linear stress profile
        idxs_const_x = np.where(qx == np.min(qx))[0]
        idx_const_x = np.max(idxs_const_x[len(idxs_const_x)
                                / 2 - 2], idxs_const_x[0])
        x_opt = np.linspace(x[idx_const_x], -x[idx_const_x], self.n_x)
        del self.spirrid.theta_vars['w']
        del self.spirrid.theta_vars['Ll']
        del self.spirrid.theta_vars['Lr']
        return x_opt


class InterpolatedSPIRRID(HasTraits):

    adaption = Instance(RangeAdaption)

    def preinterpolate(self, mu_w_x, sigma_f_cutoff, x_adapt):
        # values to create array grid
        axes_values = [sigma_f_cutoff, x_adapt]
        preinterp = NDIdxInterp(data=mu_w_x, axes_values=axes_values)
        # values to interpolate for
        load_sigma_f = self.initial_eps_vars[0]
        x = self.initial_eps_vars[1]
        interp_coords = [load_sigma_f, x]
        return preinterp(*interp_coords, mode='constant')

    interp_grid = Property()

    @cached_property
    def _get_interp_grid(self):
        print('evaluating mean response and adapting ranges...')
        spirrid_result = self.spirrid_result
        print('complete')
        axes_values = self.initial_eps_vars
        ni = NDIdxInterp(data=spirrid_result, axes_values=axes_values)
        return ni

    spirrid_result = Property(Array)

    @cached_property
    def _get_spirrid_result(self):
        Ll = self.adaption.BC_range
        Lr = Ll
        result = np.zeros((self.adaption.load_n_sigma_c, self.adaption.n_x,
                           self.adaption.n_BC, self.adaption.n_BC))
        loops_tot = len(Ll) * len(Lr)
        for i, ll in enumerate(Ll):
            for j, lr in enumerate(Lr):
                # adapt w range
                w_opt, sigma_f_opt = self.adaption.adapt_w_range(ll, lr)
                # adapt x range
                x_opt = np.linspace(-ll, lr, self.adaption.n_x)
                self.adaption.spirrid.eps_vars = dict(w=w_opt, x=x_opt)
                self.adaption.spirrid.theta_vars['Ll'] = ll
                self.adaption.spirrid.theta_vars['Lr'] = lr
                # evaluate 2D (w,x) SPIRRID with adapted ranges x and w
                mu_w_x = self.adaption.spirrid.mu_q_arr
                # preinterpolate particular result for the given x and sigma ranges
                mu_w_x_interp = \
                self.preinterpolate(mu_w_x, sigma_f_opt, x_opt).T
                mask = np.where(self.adaption.load_sigma_f
                                <= sigma_f_opt[-1], 1, np.NaN)[:, np.newaxis]
                mu_w_x_interp = mu_w_x_interp * mask
                # eps_vars = orthogonalize([np.arange(len(w_opt)), np.arange(len(x_opt))])
                # m.surf(eps_vars[0], eps_vars[1], mu_w_x * 0.05)
                # m.surf(eps_vars[0], eps_vars[1], mu_w_x_interp * 0.05)
                # m.show()
                # store the particular result for BC ll and lr into the result array
                result[:, :, i, j] = mu_w_x_interp
                current_loop = i * len(Lr) + j + 1
                print('progress: %2.1f %%' % \
                (current_loop / float(loops_tot) * 100.))
        return result

    initial_eps_vars = Property(List(Array))

    @cached_property
    def _get_initial_eps_vars(self):
        return [self.adaption.load_sigma_f,
                self.adaption.x_init,
                self.adaption.BC_range,
                self.adaption.BC_range]

    def __call__(self, *args):
        '''
        evaluation of force profile in the vicinity of a crack bridge
        '''
        args = list(args)
        load_sigma_f = args[0] / self.adaption.spirrid.theta_vars['V_f']
        # fiber stress
        args[0] = load_sigma_f
        self.load_sigma_f = load_sigma_f
        return self.interp_grid(*args)

if __name__ == '__main__':
    from etsproxy.mayavi import mlab
    from stats.spirrid import make_ogrid as orthogonalize
    from matplotlib import pyplot as plt
    from quaducom.micro.resp_func.cb_emtrx_clamped_fiber_stress import \
        CBEMClampedFiberStressSP
    from quaducom.micro.resp_func.cb_emtrx_clamped_fiber_stress_residual \
        import CBEMClampedFiberStressResidualSP

    # filaments
    r = 0.00345
    V_f = 0.0103
    tau = 0.1  # RV('uniform', loc=0.02, scale=.2)
    E_f = 200e3
    E_m = 25e3
    l = RV('uniform', scale=20., loc=2.)
    theta = 0.0
    xi = RV('weibull_min', scale=0.01, shape=5)
    phi = 1.
    Ll = 40.
    Lr = 20.
    s0 = 0.01
    m = 5.0
    Pf = RV('uniform', loc=0., scale=1.0)

    rf = CBEMClampedFiberStressSP()
    ra = RangeAdaption(load_sigma_c_max=22.0,
                       load_n_sigma_c=30,
                       n_w=100,
                       n_x=101,
                       n_BC=2,
                       spirrid=SPIRRID(q=rf,
                                       sampling_type='LHS',
                                       theta_vars=dict(tau=tau,
                                                  l=l,
                                                  E_f=E_f,
                                                  theta=theta,
                                                  xi=xi,
                                                  phi=phi,
                                                  E_m=E_m,
                                                  r=r,
                                                  V_f=V_f
                                                        ),
                                        n_int=20),
                                    )

    isp = InterpolatedSPIRRID(adaption=ra)
    print('no res', isp(isp.adaption.load_sigma_c, np.linspace(-50, 50, 20), 50., 50.))

    rf = CBEMClampedFiberStressResidualSP()
    ra = RangeAdaption(load_sigma_c_max=22.0,
                       load_n_sigma_c=30,
                       n_w=100,
                       n_x=101,
                       n_BC=2,
                       spirrid=SPIRRID(q=rf,
                                       sampling_type='LHS',
                                       theta_vars=dict(tau=tau,
                                                   l=l,
                                                   E_f=E_f,
                                                   theta=theta,
                                                   Pf=Pf,
                                                   phi=phi,
                                                   E_m=E_m,
                                                   r=r,
                                                   V_f=V_f,
                                                   m=m,
                                                   s0=s0
                                                        ),
                                        n_int=20),
                                    )

    isp = InterpolatedSPIRRID(adaption=ra)
    print('res', isp(isp.adaption.load_sigma_c, np.linspace(-50, 50, 20), 50., 50.))

    def plot():
        sigma = isp.adaption.load_sigma_c
        Ll = 50.
        Lr = 50.
        x = np.linspace(-Ll, Lr, 201)

        eps_vars = orthogonalize([np.arange(len(sigma)), np.arange(len(x))])
        # mu_q_nisp = nisp(P, x, Ll, Lr)[0]
        mu_q_isp = isp(sigma, x, 1., 14.)
        mu_q_isp2 = isp(sigma, x, Ll, Lr)

#        plt.plot(np.arange(len(sigma)), sigma/0.0103)
#        plt.plot(np.arange(len(sigma)), np.max(mu_q_isp,axis = 0))
#        plt.show()
        # n_mu_q_arr = mu_q_nisp / np.max(np.fabs(mu_q_nisp))
        mlab.surf(eps_vars[0], eps_vars[1], mu_q_isp / 10.)
        mlab.surf(eps_vars[0], eps_vars[1], mu_q_isp2 / 10.)
        mlab.show()

    # plot()
