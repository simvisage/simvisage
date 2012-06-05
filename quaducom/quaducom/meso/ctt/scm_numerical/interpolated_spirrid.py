'''
Created on Aug 17, 2011

@author: rostar
'''

from etsproxy.traits.api import HasTraits, Property, cached_property, \
    implements, Instance, Float, Array, List, Int
from stats.spirrid.spirrid import SPIRRID
from stats.spirrid.rv import RV
from quaducom.micro.resp_func.cb_emtrx_clamped_fiber import \
    CBEMClampedFiberSP
from etsproxy.traits.ui.api import View, Item, VGroup
import numpy as np
from stats.spirrid.rf import \
    RF
from math import pi
from scipy import ndimage
from mathkit.mfn import MFnLineArray
import types
from matplotlib import pyplot as plt

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


class NDIdxInterp(HasTraits):

    data = Array # nd array of values (measured, computed..) of size orthogonalize(axes_values)
    axes_values = List(Array(float)) # list of control input parameter values

    def __call__(self, *gcoords, **kw):
        ''' 
        kw: dictionary of values to interpolate for;
        len(kw) has to be equal the data dimension
        '''
        order = kw.get('order', 1)

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

        val = ndimage.map_coordinates(data, icoords, order = 1, mode = 'nearest')
        return val

    def get_icoords(self, gcoords):
        '''
        gcoords: values to be interpolated for
        this method transforms the global coords to "index" coords
        '''
        icoords = [np.interp(gcoord, axis_values, np.arange(len(axis_values)))
                   for gcoord, axis_values in zip(gcoords, self.axes_values) ]
        return icoords

import copy

class InterpolatedSPIRRID(HasTraits):

    ctrl_arr = Array

    spirrid = Instance(SPIRRID)
    
    def P2w(self, *args):

        # set the x position to 0.0
        print 'args', args
        args = list(args)
        P_given = args[0]
        args[0] = self.spirrid.evars['w']
        npoints = len(args[0])
        args[1] = 0.0
        P = self.interp_grid(*args)
        if np.argmax(P) != len(P)-1:
            raise ValueError, 'maximum force %.1f reached' % np.max(P)
        else:
            if np.max(P_given)/np.max(P) > 1.0 and np.max(P_given)/np.max(P) < 1.2:
                w = self.spirrid.evars['w']
            else:             
                while np.max(P_given)/np.max(P) < 1.0 or np.max(P_given)/np.max(P) > 1.3:
                    print 'ratio before iteration', np.max(P_given)/np.max(P)
                    print 'iterating the w interval'
                    Pw_iter_line = MFnLineArray(extrapolate = 'diff', xdata = P, ydata = self.spirrid.evars['w'])
                    wmax = Pw_iter_line.get_values(np.max(P_given)* 1.2)
                    w = np.linspace(0,wmax,npoints)
                    self.spirrid.evars['w'] = w
                    args[0] = w
                    P = self.interp_grid(*args)
                    print 'ratio after iteration',np.max(P_given)/np.max(P)
            Pw_line = MFnLineArray(xdata = P, ydata = w)       
            return Pw_line

    interp_grid = Property(depends_on = 'spirrid.evars')
    @cached_property
    def _get_interp_grid(self):
        print 'RECALCULATING'
        print 'evaluating spirrid...'
        data = self.spirrid.mu_q_arr
        print 'complete'
        axes_values = self.spirrid.evar_lst
        ni = NDIdxInterp(data = data, axes_values = axes_values)
        return ni

    def __call__(self, *args):
        '''
        evaluation of force profile in the vicinity of a crack bridge
        '''
        P = args[0]
        P2w_line = self.P2w(*args)
        w = P2w_line.get_values(P)
        return self.interp_grid(w, *args[1:]), w


class NonInterpolatedSPIRRID(HasTraits):

    ctrl_arr = Property(Array)
    def _get_ctrl_arr(self):
        return self.spirrid.evars['w']

    spirrid = Instance(SPIRRID)

    def mu_P_w(self, Ll, Lr):
        sp = copy.copy(self.spirrid)
        sp.evars.pop('x')
        sp.evars.pop('Ll')
        sp.evars.pop('Lr')
        sp.tvars['x'] = 0.0
        sp.tvars['Ll'] = Ll
        sp.tvars['Lr'] = Lr
        return sp.mu_q_arr
    
    def P2w(self, *args):
        P = self.mu_P_w(args[2], args[3])
        id_max = np.argmax(P)
        Pw_line = MFnLineArray(xdata = P[:id_max], ydata = self.ctrl_arr[:id_max])
#        plt.plot(Pw_line.ydata, Pw_line.xdata)
#        plt.show()
        return Pw_line

    def spirrid_response(self,w,*args):
        sp = copy.copy(self.spirrid)
        sp.evars = {'w':w, 'x':args[0]}
        sp.tvars['Ll'] = args[1]
        sp.tvars['Lr'] = args[2]
        return sp.mu_q_arr.T

    def __call__(self, *args):
        '''
        evaluation of force profile in the vicinity of a crack bridge
        '''
        P = args[0]
        P2w_line = self.P2w(*args)
        if np.max(P) > np.max(P2w_line.xdata):
            raise ValueError, 'maximum force %.1f reached' % np.max(P2w_line.xdata)
        w = P2w_line.get_values(P)
        return self.spirrid_response(w, *args[1:]), w


if __name__ == '__main__':


    import etsproxy.mayavi.mlab as m
    from stats.spirrid import make_ogrid as orthogonalize

    # filaments
    tau = 2.01
    Af = 5.31e-4
    Ef = 72e3
    Am = 50.
    Em = 30e3
    l = 0.0#RV('uniform', 5.0, 20.0)
    theta = 0.0
    phi = 1.
    Ll = np.linspace(0.01, 100., 3)
    Lr = np.linspace(0.01, 100., 3)
    Nf = 1700.
    xi = 50.#RV( 'weibull_min', scale = 0.017, shape = 5 )

    rf = CBEMClampedFiberSP()
    nisp = NonInterpolatedSPIRRID(spirrid = SPIRRID(q = rf,
                                          sampling_type = 'LHS',
                                       evars = dict(w = np.linspace(0.0, .1, 100),
                                                   x = np.linspace(-10., 10., 100),
                                                   Ll = Ll,
                                                   Lr = Lr,
                                                    ),
                                     tvars = dict(tau = tau,
                                                   l = l,
                                                   A_r = Af,
                                                   E_r = Ef,
                                                   theta = theta,
                                                   xi = xi,
                                                   phi = phi,
                                                   E_m = Em,
                                                   A_m = Am,
                                                   Nf = Nf,
                                                        ),
                                        n_int = 20),
                                    )

    isp = InterpolatedSPIRRID(spirrid = SPIRRID(q = rf,
                                          sampling_type = 'LHS',
                                       evars = dict(w = np.linspace(0.0, .1, 100),
                                                   x = np.linspace(-10., 10., 61),
                                                   Ll = Ll,
                                                   Lr = Lr,
                                                    ),
                                     tvars = dict(tau = tau,
                                                   l = l,
                                                   A_r = Af,
                                                   E_r = Ef,
                                                   theta = theta,
                                                   xi = xi,
                                                   phi = phi,
                                                   E_m = Em,
                                                   A_m = Am,
                                                   Nf = Nf,
                                                        ),
                                        n_int = 20),
                                    )

    def plot():
        
        P = np.linspace(0, 500, 100)
        x = np.linspace(-10., 10., 61)
        Ll = 100.
        Lr = 100.
    
        e_arr = orthogonalize([P, x])
        mu_q_nisp = nisp(P, x, Ll, Lr)[0]
        mu_q_isp = isp(P, x, Ll, Lr)[0]
        #n_mu_q_arr = mu_q_nisp / np.max(np.fabs(mu_q_nisp))
        m.surf(e_arr[0], e_arr[1], mu_q_nisp)
        m.surf(e_arr[0], e_arr[1], mu_q_isp)
        m.show()
        
    plot()