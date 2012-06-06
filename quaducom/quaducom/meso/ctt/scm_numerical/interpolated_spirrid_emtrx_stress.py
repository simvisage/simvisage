'''
Created on Aug 17, 2011

@author: rostar
'''

from enthought.traits.api import HasTraits, Property, cached_property, \
    implements, Instance, Float, Array, List, Int
from stats.spirrid.spirrid import SPIRRID
from stats.spirrid.rv import RV
from quaducom.resp_func.cb_emtrx_clamped_fiber_stress import \
    CBEMClampedFiberStressSP
from enthought.traits.ui.api import View, Item, VGroup
from quaducom.ctt.homogenized_crack_bridges.i_homogenized_cb import ICB
import numpy as np
from stats.spirrid.rf import \
    RF
from math import pi
from scipy import ndimage
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

class NDIdxInterp(HasTraits):

    data = Array # nd array of values (measured, computed..) of size orthogonalize(axes_values)
    axes_values = List(Array(float)) # list of control input parameter values

    def __call__(self, *gcoords, **kw):
        ''' 
        kw: dictionary of values to interpolate for;
        len(kw) has to be equal the data dimension
        '''
        order = kw.get('order', 1)
        if len(self.axes_values) != len(gcoords):
            raise TypeError, 'method takes {req} arguments ({given} given)'.format(req=\
                len(self.axes_values), given=len(gcoords))

        # the value to be interpolated is defined by coordinates assigning
        # the data grid; e.g. if coords = array([[1.5, 2.5]]) then coords.T
        # point to the half distance between the 1. and 2. entry in the first
        # dimension and the half distance between the 2. and 3. entry of the 2. dim
        icoords = self.get_icoords(gcoords)
        # specify the data for the interpolation
        data = self.data
        # interpolate the value (linear)
#        a, b, c = [0.5, 0.5, 0.5], [0, 0, 0], [0, 1, 2]
#        icoords = [a, b, c]
        val = ndimage.map_coordinates(data, icoords, order=order, mode='nearest')
        return val.flatten()

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

        arg_list = [np.repeat(arg, len(self.ctrl_arr))
                    for arg in args]
        args = [self.ctrl_arr, np.zeros(len(self.ctrl_arr))] + arg_list
        P = self.interp_grid(*args)
        id_max = np.argmax(P)
        Pw_line = MFnLineArray(ydata=self.ctrl_arr[:id_max], xdata=P[:id_max])
        return Pw_line

    interp_grid = Property()
    @cached_property
    def _get_interp_grid(self):
        print 'evaluating spirrid...'
        data = self.spirrid.mu_q_arr
        print 'complete'
        axes_values = self.spirrid.evar_list
        ni = NDIdxInterp(data=data, axes_values=axes_values)
        return ni

    def __call__(self, P, x, *args):
        '''
        evaluation of force profile in the vicinity of a crack bridge
        '''
        if P > np.max(self.P2w(*args).xdata):
            raise ValueError, 'maximum force reached'
        ctrl = self.P2w(*args).get_values(np.array(P))
        n = len(x)
        arg_list = [np.repeat(arg, n) for arg in args]
        args = [np.repeat(ctrl, n), x] + arg_list
        return self.interp_grid(*args)

    P_w = Property(depends_on='+modified')
    @cached_property
    def _get_P_w(self):
        '''
        evaluation of force-crack width relation for a crack bridge
        '''
        return MFnLineArray(xdata=self.P2w.ydata, ydata=self.P2w.xdata)

if __name__ == '__main__':

    import enthought.mayavi.mlab as m
    from stats.spirrid import orthogonalize
    
    rf_test = CBEMClampedFiberStressSP()
    w_guess = np.linspace(0.0, 1.0, 20)
    isp = InterpolatedSPIRRID(spirrid=SPIRRID(q=rf_test,
                                          sampling_type='LHS',
                                          evars=dict(w=w_guess,
                                                       x=np.linspace(-60., 60., 30),
                                                       Ll=np.linspace(10., 11., 2),
                                                       Lr=np.linspace(20., 21., 2),
                                                        ),
                                          tvars=dict(tau=0.15, #RV( 'uniform', 0.7, 1.0 ),
                                                       l=RV('uniform', 2.0, 10.0),
                                                       E_r=72e3,
                                                       theta=RV('uniform', 0.0, .10),
                                                       xi=0.5, #RV( 'weibull_min', scale = 0.017, shape = 5, n_int = 10 ),
                                                       phi=1.0,
                                                       E_m=30e3,
                                                       V_f=0.0175,
                                                       r=0.013
                                                        ),
                                        n_int=20),
                                    )

    x = np.linspace(-30, 30, 300)
    P = np.linspace(0., 580, 200)
    forces = []
    for p in P:
        forces.append(isp(p, x, 10., 20.))
#        plt.plot( x, mfy.get_force_x_reinf( x ) )
#    plt.show()

    e_arr = orthogonalize([x, P])
    n_e_arr = [ e / np.max(np.fabs(e)) for e in e_arr ]

    mu_q_arr = np.array(forces)

    n_mu_q_arr = mu_q_arr / np.max(np.fabs(mu_q_arr))
    m.surf(n_e_arr[0], n_e_arr[1], n_mu_q_arr)
    
    m.show()
