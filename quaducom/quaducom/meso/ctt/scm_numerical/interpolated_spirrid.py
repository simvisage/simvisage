'''
Created on Aug 17, 2011

@author: rostar
'''

from etsproxy.traits.api import HasTraits, Property, cached_property, \
    Instance, Array, List
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
        val = ndimage.map_coordinates(data, icoords, order = 1, mode = mode)
        return val

    def get_icoords(self, gcoords):
        '''
        gcoords: values to be interpolated for
        this method transforms the global coords to "index" coords
        '''
        icoords = [np.interp(gcoord, axis_values, np.arange(len(axis_values)))
                   for gcoord, axis_values in zip(gcoords, self.axes_values) ]
        return icoords

class InterpolatedSPIRRID(HasTraits):

    spirrid = Instance(SPIRRID)
    
    def adapt_w_higher_stress(self, sigma_f, load_sigma_f, len_w):
        # case 1) stretch the w range
        while np.max(sigma_f)/np.max(load_sigma_f) > 1.1 or \
            float(np.argmax(sigma_f))/float(len(sigma_f) - 1) < 0.9:
            # find the closest higher value to the max applied stress
            idx_closest = arg_find_closest_higher(sigma_f, np.max(load_sigma_f))
            wmax = self.spirrid.evars['w'][idx_closest]
            #adapt the w range and evaluate spirrid with the adapted range
            self.spirrid.evars['w'] = np.linspace(0.0, wmax, len_w)
            sigma_f = self.spirrid.mu_q_arr

    def adapt_w_lower_stress(self, sigma_f, load_sigma_f, len_w):
        # the peak is within the w range - stretch the w range
        if np.argmax(sigma_f) != len(sigma_f) -1:
            while np.argmax(sigma_f)/float(len(sigma_f)-1) < 0.9:
                wmax = self.spirrid.evars['w'][np.argmax(sigma_f)] * 1.05
                self.spirrid.evars['w'] = np.linspace(0.0, wmax, len_w)
                sigma_f = self.spirrid.mu_q_arr
        # the peak is beyond the w range
        else:
            # stretch the w range until case 1) or 2) is attained
            while np.argmax(sigma_f) == len(sigma_f) -1 and np.max(sigma_f) < np.max(load_sigma_f):
                factor = np.max(load_sigma_f)/np.max(sigma_f)
                wmax = self.spirrid.evars['w'][-1] * factor * 1.2
                self.spirrid.evars['w'] = np.linspace(0.0, wmax, len_w)
                sigma_f = self.spirrid.mu_q_arr
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

    def adapt_x_range(self, x, i, Ll, j, Lr):
        try:
            Ll[i+1]
            l_bound = np.min([Ll[i+1], -x[0]])
        except IndexError:
            l_bound = np.min([Ll[i], -x[0]])        
        try:
            Lr[j+1]
            r_bound = np.min([Lr[j+1], x[-1]])
        except IndexError:
            r_bound = np.min([Lr[j], x[-1]])
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

    interp_grid = Property(depends_on = 'spirrid.evars')
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
        w = self.spirrid.evars['w']
        x = self.spirrid.evars['x']
        Ll = self.spirrid.evars['Ll']
        Lr = self.spirrid.evars['Lr']      
        result = np.zeros((len(self.initial_evars[0]),len(x),len(Ll),len(Lr)))
        loops_tot = len(Ll)*len(Lr)
        for i, ll in enumerate(Ll):
            for j, lr in enumerate(Lr):
                # adapt w range
                sigma_f_cutoff = self.adapt_w_range(w, ll, lr)
                self.spirrid.evars['w'] = self.spirrid.evars['w'][:len(sigma_f_cutoff)]
                # adapt x range
                adapted_x = self.adapt_x_range(x, i, Ll, j, Lr)
                self.spirrid.evars['x'] = adapted_x
                del self.spirrid.tvars['x']
                # evaluate 2D (w,x) SPIRRID with adapted ranges x and w
                mu_w_x = self.spirrid.mu_q_arr
#                e_arr = orthogonalize([np.arange(len(w)), np.arange(len(x))])
#                m.surf(e_arr[0], e_arr[1], mu_w_x/np.max(mu_w_x)*50.)
#                m.show()
                # preinterpolate particular result for the given x and sigma ranges
                mu_w_x_interp = self.preinterpolate(mu_w_x, sigma_f_cutoff, self.spirrid.evars['x']).T
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
        load_sigma_f = args[0]/self.spirrid.tvars['V_f']
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
    tau = 0.1 #RV('uniform', loc = 0.02, scale = .01) # 0.5
    Ef = 200e3
    Em = 25e3
    l = RV( 'uniform', scale = 10., loc = 0. )
    theta = 0.0
    xi = 0.01#RV( 'weibull_min', scale = 0.12, shape = 5 ) # 0.017
    phi = 1.
    w = np.linspace(0.0, .1, 51)
    x = np.linspace(-50., 50., 51)
    Ll = np.linspace(0.5,50,10)
    Lr = np.linspace(0.5,50,10)

    rf = CBEMClampedFiberStressSP()
    isp = InterpolatedSPIRRID(spirrid = SPIRRID(q = rf,
                                          sampling_type = 'LHS',
                                       evars = dict(w = w,
                                                   x = x,
                                                   Ll = Ll,
                                                   Lr = Lr,
                                                    ),
                                     tvars = dict(tau = tau,
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

    def plot():
        sigma = np.linspace(0, 20, 120)
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
#    nisp = NonInterpolatedSPIRRID(spirrid = SPIRRID(q = rf,
#                                          sampling_type = 'LHS',
#                                       evars = dict(w = np.linspace(0.0, .7, 31),
#                                                   x = np.linspace(-10., 10., 11),
#                                                   Ll = Ll,
#                                                   Lr = Lr,
#                                                    ),
#                                     tvars = dict(tau = tau,
#                                                   l = l,
#                                                   E_f = Ef,
#                                                   theta = theta,
#                                                   xi = xi,
#                                                   phi = phi,
#                                                   E_m = Em,
#                                                   r = r,
#                                                   V_f = Vf
#                                                        ),
#                                        n_int = 20),
#                                    )