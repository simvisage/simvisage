'''
Created on Aug 17, 2011

@author: rostar
'''

from enthought.traits.api import HasTraits, Property, cached_property, \
    implements, Instance, Float, Array, List, Int
from stats.spirrid.spirrid import SPIRRID
from stats.spirrid.rv import RV
from quaducom.micro.resp_func.cb_emtrx_clamped_fiber import \
    CBEMClampedFiberSP
from enthought.traits.ui.api import View, Item, VGroup
import numpy as np
from stats.spirrid.rf import \
    RF
from math import pi
from scipy import ndimage
from mathkit.mfn import MFnLineArray
#import scitools.numpytools as st
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

        val = ndimage.map_coordinates(data, icoords, order = order, mode = 'nearest')
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

    ctrl_arr = Property(Array)
    def _get_ctrl_arr(self):
        return self.spirrid.evars['w']

    spirrid = Instance(SPIRRID)

    def P2w(self, *args):

        # set the x position to 0.0
        args = list(args)
        args[0] = self.ctrl_arr
        args[1] = 0.0
        P = self.interp_grid(*args)

        id_max = np.argmax(P)
        Pw_line = MFnLineArray(xdata = P[:id_max], ydata = self.ctrl_arr[:id_max])
        return Pw_line

    interp_grid = Property()
    @cached_property
    def _get_interp_grid(self):
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
        if np.any(P > np.max(P2w_line.xdata)):
            raise ValueError, 'maximum force %.1f reached' % np.max(P2w_line.xdata)
        w = P2w_line.get_values(P)
        return self.interp_grid(w, *args[1:]), w

if __name__ == '__main__':


    import enthought.mayavi.mlab as m
    from stats.spirrid import make_ogrid as orthogonalize
    from matplotlib import pyplot as plt

    # filaments
    tau = 2.0
    Af = 5.31e-4
    Ef = 72e3
    Am = 50.
    Em = 30e3
    l = RV('uniform', 5.0, 20.0)
    theta = 0.0
    phi = 1.
    Ll = np.linspace(0.01, 100., 5)
    Lr = np.linspace(0.01, 100., 5)
    Nf = 1700.
    xi = 50.#RV( 'weibull_min', scale = 0.017, shape = 5 )

    rf = CBEMClampedFiberSP()
    isp = InterpolatedSPIRRID(spirrid = SPIRRID(q = rf,
                                          sampling_type = 'LHS',
                                       evars = dict(w = np.linspace(0.0, .25, 40),
                                                   x = np.linspace(-100., 100., 100),
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

#    e = isp.spirrid.evars['x']
#    print isp.spirrid.mu_q_arr.shape
#    q = isp.spirrid.mu_q_arr[25, :, 4, 4]
#    print q
#    print 'qm =', np.max(q) / rf.Kc * rf.Kr
#    plt.plot(e, q)
#    plt.show()

    data = isp.spirrid.mu_q_arr
    axes_values = isp.spirrid.evar_lst

    P = np.linspace(0, 250, 500)
    x = np.linspace(-200., 200., 500)
    Ll = 6.67779632721
    Lr = 5.00834724541

    ni = NDIdxInterp(data = data, axes_values = axes_values)

    pr = (P / rf.Kc * rf.Kr)[np.newaxis, :] * np.ones_like(x)[:, np.newaxis]


    e_arr = orthogonalize([P, x])
    n_e_arr = [ e / np.max(np.fabs(e)) for e in e_arr ]

    mu_q_arr = isp(P, x, Ll, Lr)[0]
    #pr = isp(P, x, 10., Lr)
    n_mu_q_arr = mu_q_arr / np.max(np.fabs(mu_q_arr))
    #n_pr = pr / np.max(np.fabs(mu_q_arr))
    m.surf(n_e_arr[0], n_e_arr[1], n_mu_q_arr)
    #m.surf(n_e_arr[0], n_e_arr[1], n_pr)

    m.show()
