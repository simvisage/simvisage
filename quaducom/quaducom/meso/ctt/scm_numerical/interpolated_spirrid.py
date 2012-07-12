'''
Created on Aug 17, 2011

@author: rostar
'''

from etsproxy.traits.api import HasTraits, Property, cached_property, \
    Instance, Array, List, Float, Int
from stats.spirrid.spirrid import SPIRRID
from stats.spirrid.rv import RV
from quaducom.micro.resp_func.cb_emtrx_clamped_fiber_stress import \
    CBEMClampedFiberStressSP
import numpy as np
from scipy import ndimage
from matplotlib import pyplot as plt
import types

def orthogonalize_filled(args):
    '''
    creates meshgrid up to third dimension
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
        raise NotImplementedError, 'max number of arrays as input is currently 3'

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

def find_closest_higher(array, scalar):
    idx = np.argwhere((array - scalar) > 0.0).flat[0]
    return array.flat[idx]

def arg_find_closest_higher(array, scalar):
    idx = np.argwhere((array - scalar) > 0.0).flat[0]
    return idx

class NDIdxInterp(HasTraits):

    data = Array # nd array of values (measured, computed..) of size orthogonalize(axes_values)
    axes_values = List(Array(float)) # list of control input parameter values

    def __call__(self, *gcoords, **kw):
        ''' 
        kw: dictionary of values to interpolate for;
        len(kw) has to be equal the data dimension
        '''
        order = kw.get('order', 1)
        mode = kw.get('mode', 'nearest')

        # check if the number of dimensions to interpolate in equals the number of given coordinates
        if len(self.axes_values) != len(gcoords):
            raise TypeError, 'method takes {req} arguments ({given} given)'.format(req = \
                len(self.axes_values), given = len(gcoords))
        icoords = self.get_icoords(gcoords)

        # create a meshgrid for the interpolation
        icoords = orthogonalize_filled(icoords)
        data = self.data
        # interpolate the value (linear)
        # a, b, c = [0.5, 0.5, 0.5], [0, 0, 0], [0, 1, 2]
        # icoords = [a, b, c]
        val = ndimage.map_coordinates(data, icoords, order = 1)
        return val

    def get_icoords(self, gcoords):
        '''
        gcoords: values to be interpolated for
        this method transforms the global coords to "index" coords
        '''
        icoords = [np.interp(gcoord, axis_values, np.arange(len(axis_values)))
                   for gcoord, axis_values in zip(gcoords, self.axes_values) ]
        return icoords

class RangeAdaption(HasTraits):
    ''' incorporates methods for defining the ranges of BCs, w and x'''
    
    spirrid = Instance(SPIRRID)
    max_load_sigma = Float
    n_load_sigma = Int
    n_w = Int
    n_x = Int
    n_BC = Int
    
    load_sigma = Property(Array, depends_on = 'n_load_sigma, max_load_sigma')
    @cached_property
    def _get_load_sigma(self):
        return np.linspace(1e-6, self.max_load_sigma, self.n_load_sigma)

    load_sigma_f = Property(Array, depends_on = 'n_load_sigma, max_load_sigma')
    @cached_property
    def _get_load_sigma_f(self):
        Vf = self.spirrid.tvars['V_f']
        return self.load_sigma / Vf

    stress_transm_props = Property(Float, depends_on = 'spirrid.tvars')
    @cached_property
    def _get_stress_transm_props(self):
        ''' evaluates the half shielded length at maximum load'''
        r = self.spirrid.tvars['r']
        if isinstance(r, RV):
            r = r._distr.mean
        tau = self.spirrid.tvars['tau']
        if isinstance(tau, RV):
            tau = tau._distr.mean
        l = self.spirrid.tvars['l']
        if isinstance(l, RV):
            l = l._distr.mean
        Vf = self.spirrid.tvars['V_f']
        w_approx = (r*self.max_load_sigma/2./tau + l/2.) * self.max_load_sigma * Vf
        return r*self.max_load_sigma/2./tau, l/2., w_approx
    
    delta = Property(Float, depends_on = 'spirrid.tvars')
    @cached_property
    def _get_delta(self):
        return self.stress_transm_props[0] + self.stress_transm_props[1]

    def adapt_BC_range(self):
        ''' evaluates the necessary range for BCs and delivers a linspace of Lr and Ll '''
        # set the BCs to infinity and evaluate max crack opening wmax
        ll, lr = 10e15, 10e15
        # start with an approximate w range
        w_init = np.linspace(0, 1.5 * self.stress_transm_props[2], self.n_w)
        # adapt the w range to the optimum
        w_opt, sigma_f_opt_w = self.adapt_w_range(ll, lr, w_init)
        # get the crack opening at peak stress and initial guess for the x range
        x_init = np.linspace(-self.delta, self.delta, self.n_x)
        plt.plot(sigma_f_opt_w)
        plt.show()
        print np.argmax(sigma_f_opt_w)
        wmax = w_opt[np.argmax(sigma_f_opt_w)[0]]
        # evaluate the max BC value
        x_opt, sigma_f_opt_x = self.adapt_x_range(ll, lr, x_init, wmax)   
        
    def adapt_w_range(self, ll, lr, w):
        ''' method for adapting the w range for given BC '''
        # three cases of the peak response are distinguished:
        # 1) response includes a higher value than the applied load
        # 2) response peak is lower and lies within the given w range
        # 3) response peak is lower and lies beyond the given w range
        self.spirrid.evars = dict(w = w)
        self.spirrid.tvars['x'] = 0.0
        self.spirrid.tvars['Ll'] = ll
        self.spirrid.tvars['Lr'] = lr
        sigma_f = self.spirrid.mu_q_arr
        # finding appropriate range for w
        load_sigma_f = self.load_sigma_f
        # is the evaluated stress higher than the applied stress? case 1)
        # if so, the CB will not break during the loading process
        if np.max(sigma_f) >= np.max(load_sigma_f):
            w_opt, sigma_opt = self.adapt_w_higher_stress(sigma_f, load_sigma_f, self.n_w)
        # the evaluated stress has a peak lower than the maximum of the applied stress
        else:
            w_opt, sigma_opt = self.adapt_w_lower_stress(sigma_f, load_sigma_f, self.n_w)
        # cut off the part after peak value and invert w and sigma_f
        idxmax = np.argmax(sigma_opt)
        sigma_f_cutoff = sigma_opt[:idxmax+1]
        w_cutoff = w_opt[:idxmax]
        return sigma_f_cutoff, w_cutoff
    
    def adapt_w_higher_stress(self, sigma_f, load_sigma_f, len_w):
        # case 1) stretch the w range
        count = 0
        while np.max(sigma_f)/np.max(load_sigma_f) > 1.1 or \
            float(np.argmax(sigma_f))/float(len(sigma_f) - 1) < 0.9:
            # find the closest higher value to the max applied stress
            idx_closest = arg_find_closest_higher(sigma_f, np.max(load_sigma_f))
            wmax = self.spirrid.evars['w'][idx_closest]
            #adapt the w range and evaluate spirrid with the adapted range
            self.spirrid.evars['w'] = np.linspace(0.0, wmax, len_w)
            sigma_f = self.spirrid.mu_q_arr
            count += 1
            if count > 3:
                raise ValueError('got stuck in a loop adapting w - try to change the w range')
        return self.spirrid.evars['w'], sigma_f

    def adapt_w_lower_stress(self, sigma_f, load_sigma_f, len_w):
        # the peak is within the w range - stretch the w range
        if np.argmax(sigma_f) != len(sigma_f) -1:
            count = 0
            while np.argmax(sigma_f)/float(len(sigma_f)-1) < 0.9:
                wmax = self.spirrid.evars['w'][np.argmax(sigma_f)] * 1.05
                self.spirrid.evars['w'] = np.linspace(0.0, wmax, len_w)
                sigma_f = self.spirrid.mu_q_arr
                count += 1
                if count > 3:
                    raise ValueError('got stuck in a loop adapting w - try to change the w range')
            return self.spirrid.evars['w'], sigma_f
        # the peak is beyond the w range
        else:
            # stretch the w range until case 1) or 2) is attained
            count = 0
            while np.argmax(sigma_f) == len(sigma_f) -1 and np.max(sigma_f) < np.max(load_sigma_f):
                factor = np.max(load_sigma_f)/np.max(sigma_f)
                wmax = self.spirrid.evars['w'][-1] * factor * 1.2
                self.spirrid.evars['w'] = np.linspace(0.0, wmax, len_w)
                sigma_f = self.spirrid.mu_q_arr
                count += 1
                if count > 3:
                    raise ValueError('got stuck in a loop adapting w - try to change the w range')
            # case 1)
            if np.argmax(sigma_f) == len(sigma_f) -1 and np.max(sigma_f) > np.max(load_sigma_f):
                w_opt, sigma_opt = self.adapt_w_higher_stress(sigma_f, load_sigma_f, len_w)
            elif np.argmax(sigma_f) != len(sigma_f) -1 and np.max(sigma_f) > np.max(load_sigma_f):
                w_opt, sigma_opt = self.adapt_w_higher_stress(sigma_f, load_sigma_f, len_w)
            # case 2)
            elif np.argmax(sigma_f) != len(sigma_f) -1 and np.max(sigma_f) < np.max(load_sigma_f):
                w_opt, sigma_opt = self.adapt_w_lower_stress(sigma_f, load_sigma_f, len_w)
            return w_opt, sigma_opt

    def adapt_x_range(self, ll, lr, x_init, wmax):
        ''' adapts the x range for a given crack opening '''
        self.spirrid.tvars['w'] = wmax
        self.spirrid.tvars['Ll'] = ll
        self.spirrid.tvars['Lr'] = lr
        self.spirrid.evars = dict(x = x_init)
        del self.spirrid.tvars['x']
        Em = self.spirrid.tvars['E_m']
        Ef = self.spirrid.tvars['E_f']
        Vf = self.spirrid.tvars['V_f']
        qx = self.spirrid.mu_q_arr
        qmax = np.max(qx)
        Ec = Em * (1.-Vf) + Ef * Vf
        q_const = qmax * Vf * Ef/Ec
        x = x_init
        print self.spirrid.tvars
        while np.abs(np.min(qx) - q_const) > q_const/10e6:
            print 'looping x range opt'
            x = np.linspace(x[0] - 0.5 * self.delta,
                            x[-1] + 0.5 * self.delta,
                            self.n_x)
            self.spirrid.evars['x'] = x
            qx = self.spirrid.mu_q_arr
            print q_const, np.min(qx)
            plt.plot(qx)
            plt.show()   
        x_bound = qx[np.where((qx - q_const) < q_const/10e6)[0][0]]
        x_opt = np.linspace(x_bound, -x_bound, self.n_x)
        self.spirrid.evars['x'] = x_opt
        sigma_f_opt_x = self.spirrid.mu_q_arr
        plt.plot(np.repeat(q_const,len(qx)))
        plt.plot(sigma_f_opt_x, color = 'blue')
        plt.plot(qx, lw = 2, color = 'red')
        plt.show()
        return x_opt, sigma_f_opt_x

class InterpolatedSPIRRID(HasTraits):

    adaption = Instance(RangeAdaption)
    
    def adapt_w_higher_stress(self, sigma_f, load_sigma_f, len_w):
        # case 1) stretch the w range
        count = 0
        while np.max(sigma_f)/np.max(load_sigma_f) > 1.1 or \
            float(np.argmax(sigma_f))/float(len(sigma_f) - 1) < 0.9:
            # find the closest higher value to the max applied stress
            idx_closest = arg_find_closest_higher(sigma_f, np.max(load_sigma_f))
            wmax = self.spirrid.evars['w'][idx_closest]
            #adapt the w range and evaluate spirrid with the adapted range
            self.spirrid.evars['w'] = np.linspace(0.0, wmax, len_w)
            sigma_f = self.spirrid.mu_q_arr
            count += 1
            if count > 3:
                raise ValueError('got stuck in a loop adapting w - try to change the w range')

    def adapt_w_lower_stress(self, sigma_f, load_sigma_f, len_w):
        # the peak is within the w range - stretch the w range
        if np.argmax(sigma_f) != len(sigma_f) -1:
            count = 0
            while np.argmax(sigma_f)/float(len(sigma_f)-1) < 0.9:
                wmax = self.spirrid.evars['w'][np.argmax(sigma_f)] * 1.05
                self.spirrid.evars['w'] = np.linspace(0.0, wmax, len_w)
                sigma_f = self.spirrid.mu_q_arr
                count += 1
                if count > 3:
                    raise ValueError('got stuck in a loop adapting w - try to change the w range')
        # the peak is beyond the w range
        else:
            # stretch the w range until case 1) or 2) is attained
            count = 0
            while np.argmax(sigma_f) == len(sigma_f) -1 and np.max(sigma_f) < np.max(load_sigma_f):
                factor = np.max(load_sigma_f)/np.max(sigma_f)
                wmax = self.spirrid.evars['w'][-1] * factor * 1.2
                self.spirrid.evars['w'] = np.linspace(0.0, wmax, len_w)
                sigma_f = self.spirrid.mu_q_arr
                count += 1
                if count > 3:
                    raise ValueError('got stuck in a loop adapting w - try to change the w range')
            # case 1)
            if np.argmax(sigma_f) == len(sigma_f) -1 and np.max(sigma_f) > np.max(load_sigma_f):
                self.adapt_w_higher_stress(sigma_f, load_sigma_f, len_w)
            elif np.argmax(sigma_f) != len(sigma_f) -1 and np.max(sigma_f) > np.max(load_sigma_f):
                self.adapt_w_higher_stress(sigma_f, load_sigma_f, len_w)
            # case 2)
            elif np.argmax(sigma_f) != len(sigma_f) -1 and np.max(sigma_f) < np.max(load_sigma_f):
                self.adapt_w_lower_stress(sigma_f, load_sigma_f, len_w)

    def adapt_w_range(self, w, ll, lr):
        # three cases of the peak response are distinguished:
        # 1) response includes a higher value than the applied load
        # 2) response peak is lower and lies within the given w range
        # 3) response peak is lower and lies beyond the given w range
        self.spirrid.evars = dict(w = self.spirrid.evars['w'])
        self.spirrid.tvars['x'] = 0.0
        self.spirrid.tvars['Ll'] = ll
        self.spirrid.tvars['Lr'] = lr
        sigma_f = self.spirrid.mu_q_arr
        # finding appropriate range for w
        load_sigma_f = self.initial_evars[0]
        # is the evaluated stress higher than the applied stress? case 1)
        # if so, the CB will not break during the loading process
        if np.max(sigma_f) >= np.max(load_sigma_f):
            self.adapt_w_higher_stress(sigma_f, load_sigma_f, len(w))
        # the evaluated stress has a peak lower than the maximum of the applied stress
        else:
            self.adapt_w_lower_stress(sigma_f, load_sigma_f, len(w))
        # cut off the part after peak value and invert w and sigma_f
        idxmax = np.argmax(self.spirrid.mu_q_arr)
        sigma_f_cutoff = self.spirrid.mu_q_arr[:idxmax+1]
        return sigma_f_cutoff
    
    delta = Property(Float)
    @cached_property
    def _get_delta(self):
        r = self.spirrid.tvars['r']
        if isinstance(r, RV):
            r = r._distr.mean
        tau = self.spirrid.tvars['tau']
        if isinstance(tau, RV):
            tau = tau._distr.mean
        l = self.spirrid.tvars['l']
        if isinstance(l, RV):
            l = l._distr.mean
        sigma = self.load_sigma_f[-1]
        return r*sigma/2./tau + l/2.

    def adapt_x_range(self, x, wmax, i, Ll, j, Lr):
        adapt = RangeAdaption(spirrid = self.spirrid)
        adapt.adapt_BC_range()
        
        try:
            Ll[i+1]
            l_bound = np.min([Ll[i+1], 2*self.delta])
        except IndexError:
            l_bound = np.min([Ll[i], 2*self.delta])        
        try:
            Lr[j+1]
            r_bound = np.min([Lr[j+1], 2*self.delta])
        except IndexError:
            r_bound = np.min([Lr[j], 2*self.delta])
        return np.linspace(-l_bound, r_bound, len(x))

    def preinterpolate(self, mu_w_x, sigma_f_cutoff, x_adapt):
        # values to create array grid
        axes_values = [sigma_f_cutoff, x_adapt]
        preinterp = NDIdxInterp(data = mu_w_x, axes_values = axes_values)
        # values to interpolate for
        load_sigma_f = self.initial_evars[0]
        x = self.initial_evars[1]
        interp_coords = [load_sigma_f, x]
        return preinterp(*interp_coords, mode = 'constant')

    interp_grid = Property()
    @cached_property
    def _get_interp_grid(self):
        print 'evaluating mean response and adapting ranges...'        
        spirrid_result = self.spirrid_result
        print 'complete'
        axes_values = self.initial_evars
        ni = NDIdxInterp(data = spirrid_result, axes_values = axes_values)
        return ni

    spirrid_result = Property(Array)
    @cached_property
    def _get_spirrid_result(self):
        self.adaption.adapt_BC_range()
        w = self.adaption.spirrid.evars['w']
        x = self.adaption.spirrid.evars['x']
        Ll = self.adaption.spirrid.evars['Ll']
        Lr = self.adaption.spirrid.evars['Lr']      
        result = np.zeros((len(self.initial_evars[0]),len(x),len(Ll),len(Lr)))
        loops_tot = len(Ll)*len(Lr)
        for i, ll in enumerate(Ll):
            for j, lr in enumerate(Lr):
                # adapt w range
                sigma_f_cutoff = self.adapt_w_range(w, ll, lr)
                w_adapted = self.adaption.spirrid.evars['w'][:len(sigma_f_cutoff)]
                # adapt x range
                x_adapted = self.adapt_x_range(x, w_adapted[-1], i, Ll, j, Lr)
                self.spirrid.evars['w'] = w_adapted
                del self.spirrid.tvars['w']
                self.spirrid.evars['x'] = x_adapted
                # evaluate 2D (w,x) SPIRRID with adapted ranges x and w
                mu_w_x = self.spirrid.mu_q_arr
                # preinterpolate particular result for the given x and sigma ranges
                mu_w_x_interp = self.preinterpolate(mu_w_x, sigma_f_cutoff, self.spirrid.evars['x']).T
                mask = np.where(self.load_sigma_f <= sigma_f_cutoff[-1], 1, np.NaN)[:,np.newaxis]               
                mu_w_x_interp = mu_w_x_interp*mask
#                e_arr = orthogonalize([np.arange(len(self.initial_evars[0])), np.arange(len(x))])
#                m.surf(e_arr[0], e_arr[1], mu_w_x_interp/np.max(mu_w_x)*50.)
#                m.show()
                # store the particular result for BC ll and lr into the result array 
                result[:,:,i,j] = mu_w_x_interp
                current_loop = i*len(Lr)+j+1
                print 'progress: %2.1f %%' %(current_loop/float(loops_tot)*100.)
        return result

    load_sigma_f = Array    
    initial_evars = Property(List(Array))
    @cached_property
    def _get_initial_evars(self):
        return [self.load_sigma_f,
                self.spirrid.evars['x'],
                self.spirrid.evars['Lr'],
                self.spirrid.evars['Ll']]

    def __call__(self, *args):
        '''
        evaluation of force profile in the vicinity of a crack bridge
        '''
        args = list(args)
        load_sigma_f = args[0]/self.adaption.spirrid.tvars['V_f']
        # fiber stress
        args[0] = load_sigma_f
        self.load_sigma_f = load_sigma_f
        return self.interp_grid(*args)

if __name__ == '__main__':
    import etsproxy.mayavi.mlab as m
    from stats.spirrid import make_ogrid as orthogonalize

    # filaments
    r = 0.00345
    Vf = 0.0103
    tau = 0.5 #RV('uniform', loc = 0.02, scale = .01) # 0.5
    Ef = 200e3
    Em = 25e3
    l = RV( 'uniform', scale = 10., loc = 0. )
    theta = 0.0
    xi = 0.01#RV( 'weibull_min', scale = 0.12, shape = 5 ) # 0.017
    phi = 1.
    w = np.linspace(0.0, .8, 121)
    #x = np.linspace(-50., 50., 51)
    #Ll = np.linspace(0.5,50,3)
    #Lr = np.linspace(0.5,50,3)

    rf = CBEMClampedFiberStressSP()
    ra = RangeAdaption(max_load_sigma = 20.0,
                       n_load_sigma = 100,
                       n_w = 50,
                       n_x = 30,
                       n_BC = 10,
                       spirrid = SPIRRID(q = rf,
                                          sampling_type = 'LHS',
                                       evars = dict(w = w,
                                                   #x = x,
                                                   #Ll = Ll,
                                                   #Lr = Lr,
                                                    ),
                                     tvars = dict(x = 0.0,
                                                  Ll = 20000000.,
                                                  Lr = 20000000.,
                                                  tau = tau,
                                                   l = l,
                                                   E_f = Ef,
                                                   theta = theta,
                                                   xi = xi,
                                                   phi = phi,
                                                   E_m = Em,
                                                   r = r,
                                                   V_f = Vf
                                                        ),
                                        n_int = 30),
                                    )
    
    isp = InterpolatedSPIRRID(adaption = ra)
    def plot():
        sigma = np.linspace(0, 30, 120)
        Ll = 20.
        Lr = 30.
        x = np.linspace(-Ll, Lr, 100)
    
        e_arr = orthogonalize([np.arange(len(sigma)), np.arange(len(x))])
        #mu_q_nisp = nisp(P, x, Ll, Lr)[0]
        mu_q_isp = isp(sigma, x, Ll, Lr)
        #n_mu_q_arr = mu_q_nisp / np.max(np.fabs(mu_q_nisp))
        #m.surf(e_arr[0], e_arr[1], mu_q_nisp)
        m.surf(e_arr[0], e_arr[1], mu_q_isp/50.)
        m.show()
        
    plot()