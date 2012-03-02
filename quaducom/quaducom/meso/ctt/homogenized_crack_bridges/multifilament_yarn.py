'''
Created on Aug 17, 2011

@author: rostar
'''

from enthought.traits.api import HasTraits, Property, cached_property, \
    implements, Instance, Float, Array, List, Int
from stats.spirrid.spirrid import SPIRRID
from stats.spirrid.rv import RV
from quaducom.resp_func.cb_clamped_fiber import \
    CBClampedFiberSP
from quaducom.resp_func.cb_emtrx_clamped_fiber import \
    CBEMClampedFiberSP
from enthought.traits.ui.api import View, Item, VGroup
from quaducom.ctt.homogenized_crack_bridges.i_homogenized_cb import ICB
import numpy as np
from stats.spirrid.rf import \
    RF
from math import pi
from scipy import ndimage
from mathkit.mfn import MFnLineArray

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
            raise TypeError, 'method takes {req} arguments ({given} given)'.format(req = \
                len(self.axes_values), given = len(gcoords))

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
        val = ndimage.map_coordinates(data, icoords, order = order, mode = 'nearest')
        return val.flatten()

    def get_icoords(self, gcoords):
        '''
        gcoords: values to be interpolated for
        this method transforms the global coords to "index" coords
        '''
        icoords = [np.interp(gcoord, axis_values, np.arange(len(axis_values)))
                   for gcoord, axis_values in zip(gcoords, self.axes_values) ]
        return icoords


class MultifilamentYarn(HasTraits):

    implements(ICB)

    # composite length
    length = Float(modified = True) # [mm]

    # applied force
    P = Float(modified = True) # [N]

    # closest crack from left
    Ll = Float(modified = True) # [mm]

    # closest crack from right
    Lr = Float(modified = True) # [mm]

    rf = Instance(RF)
    def _rf_default(self):
        #return CBClampedFiberSP()
        return CBEMClampedFiberSP()

    Er = Float(72e3, auto_set = False, enter_set = True,
               desc = 'AR-glass modulus of elasticity [N/mm2]', modified = True, param = True)

    Em = Float(30e3, auto_set = False, enter_set = True,
               desc = 'matrix modulus of elasticity [N/mm2]', modified = True, param = True)

    Ar = Float(0.89, auto_set = False, enter_set = True,
           desc = 'yarn cross sectional area in [mm2]', modified = True, param = True)

    Nf = Float(1700., auto_set = False, enter_set = True,
           desc = 'number of parallel fibers', modified = True, param = True)

    Ac = Float(20., auto_set = False, enter_set = True,
             desc = 'composite cross section [mm2]', modified = True, param = True)

    tau = Float(0.15, auto_set = False, enter_set = True,
         desc = 'shear force per unit area [N/mm2]', modified = True, param = True)

    theta = Float(0.0, auto_set = False, enter_set = True,
         desc = 'slack', modified = True, param = True)

    xi = Float(0.0179, auto_set = False, enter_set = True,
         desc = 'breaking strain', modified = True, param = True)

    phi = Float(1.0, auto_set = False, enter_set = True,
         desc = 'bond quality', modified = True, param = True)

    l = Float(2.0, auto_set = False, enter_set = True,
         desc = 'free length [mm]', modified = True, param = True)


    w_arr = np.linspace(0.0, 1.0, 60)
    x_arr = np.linspace(-60., 60., 100)
    Ll_arr = np.linspace(10., 11., 2)
    Lr_arr = np.linspace(20., 21., 2)


    Am = Property(depends_on = 'Ac, Ar')
    @cached_property
    def _get_Am(self):
        return self.Ac - self.Ar

    Kr = Property(depends_on = 'Ar, Er')
    @cached_property
    def _get_Kr(self):
        return self.Ar * self.Er

    Km = Property(depends_on = 'Ac, Ar, Em')
    @cached_property
    def _get_Km(self):
        return self.Am * self.Em

    Kc = Property(depends_on = 'Ac, Ar, Er, Em')
    @cached_property
    def _get_Kc(self):
        return self.Kr + self.Km

    spirrid_grid = Property(depends_on = 'rf')
    @cached_property
    def _get_spirrid_grid(self):
        '''
        delivers the full grid of mean forces;
        grid size and dimensions depend on the specified control variables 
        '''

        # specify the single filament crack bridge parameters
        rf = self.rf
        s = SPIRRID(q = rf,
             sampling_type = 'LHS',
             evars = dict(w = self.w_arr,
                           x = self.x_arr,
                           Ll = self.Ll_arr,
                           Lr = self.Lr_arr,
                            ),
             tvars = dict(tau = self.tau, #RV( 'uniform', 0.7, 1.0 ),
                           l = RV('uniform', 2.0, 10.0),
                           A_r = self.Ar, #26e-3,
                           E_f = self.Er, #72e3,
                           theta = RV('uniform', 0.0, .10),
                           xi = 0.5, #RV( 'weibull_min', scale = 0.017, shape = 5, n_int = 10 ),
                           phi = self.phi,
                           E_m = self.Em,
                           A_m = self.Am,
                           Nf = self.Nf
                            ),
             n_int = 20)

        return s.evar_list, s.mu_q_arr

    P2w = Property(depends_on = 'rf, Ll, Lr')
    @cached_property
    def _get_P2w(self):
        w, x, ll, lr = self.w_arr, np.zeros(len(self.w_arr)), \
        np.repeat(self.Ll, len(self.w_arr)), np.repeat(self.Lr, len(self.w_arr))
        P = self.interp_grid(w, x, ll, lr)
        id_max = np.argmax(P)
        Pw_line = MFnLineArray(ydata = self.w_arr[:id_max], xdata = P[:id_max])
        return Pw_line

    interp_grid = Property(depends_on = 'rf, length')
    @cached_property
    def _get_interp_grid(self):
        print 'evaluating spirrid...'
        data = self.spirrid_grid[1]
        print 'complete'
        axes_values = self.spirrid_grid[0]
        ni = NDIdxInterp(data = data, axes_values = axes_values)
        return ni

    def get_force_x_reinf(self, x):
        '''
        evaluation of force profile in the vicinity of a crack bridge
        '''
        if self.P > np.max(self.P2w.xdata):
            raise ValueError, 'maximum force reached'
        w = self.P2w.get_values(np.array(self.P))
        n = len(x)
        w, x, ll, lr = np.repeat(w, n), x, np.repeat(self.Ll, n), np.repeat(self.Ll, n)
        return self.interp_grid(w, x, ll, lr)

    def get_eps_x_reinf(self, x):
        '''
        evaluation of strain profile in the vicinity of a crack bridge
        '''
        return self.get_force_x_reinf(x) / self.Kr

    def get_sigma_x_reinf(self, x):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        return self.get_force_x_reinf(x) / self.Ar

    P_w = Property(depends_on = '+modified')
    @cached_property
    def _get_P_w(self):
        '''
        evaluation of force-crack width relation for a crack bridge
        '''
        return MFnLineArray(xdata = self.P2w.ydata, ydata = self.P2w.xdata)

    traits_view = View(
                       VGroup(
                           Item('dr', resizable = False, springy = True),
                           Item('Er', resizable = False, springy = False),
                           Item('Em', resizable = False, springy = False),
                           Item('tau', resizable = False, springy = False),
                           Item('Ac', resizable = False, springy = False),
                           springy = True,
                           label = 'CB parameters',
                           dock = 'tab',
                           id = 'cb.multifilament_yarn.params',
                        ),
                            id = 'cb.multifilament_yarn',
                            dock = 'fixed',
                            scrollable = True,
                            resizable = True,
                            height = 0.8, width = 0.8
                                   )

class MultifilamentYarnRM(MultifilamentYarn):

    rf = Instance(RF)
    def _rf_default(self):
        return CBClampedFiberSP()

    spirrid_grid = Property(depends_on = 'rf')
    @cached_property
    def _get_spirrid_grid(self):
        '''
        delivers the full grid of mean forces;
        grid size and dimensions depend on the specified control variables 
        '''

        # specify the single filament crack bridge parameters
        rf = self.rf
        rf.tau = self.tau * np.sqrt(self.Nf)
        rf.l = self.l
        rf.D_f = np.sqrt(self.Ar / pi) * 2
        rf.E_f = self.Er
        rf.theta = self.theta
        rf.xi = self.xi
        rf.phi = self.phi
        rf.Ll = self.Ll
        rf.Lr = self.Lr

        s = SPIRRID(q = rf,
             sampling_type = 'LHS',
             evars = dict(w = self.w_arr,
                           x = self.x_arr,
                           Ll = self.Ll_arr,
                           Lr = self.Lr_arr,
                            ),
             tvars = dict(tau = RV('uniform', 7., 15.0),
                           l = RV('uniform', 2.0, 10.0),
                           D_f = rf.D_f, #26e-3,
                           E_f = rf.E_f, #72e3,
                           theta = rf.theta, #0.0,
                           xi = 0.5, #RV( 'weibull_min', scale = 0.017, shape = 5, n_int = 10 ),
                           phi = rf.phi,
                            ),
             n_int = 20)

        return s.evar_list, s.mu_q_arr

if __name__ == '__main__':

    import enthought.mayavi.mlab as m
    from stats.spirrid import orthogonalize
    from matplotlib import pyplot as plt

#    mfy = MultifilamentYarnRM()
#    mfy.Ll = 10.
#    mfy.Lr = 25.
#
#    x = np.linspace(-30, 30, 300)
#    P = np.linspace(0., 580, 200)
#    forces = []
#    for p in P:
#        mfy.P = p
#        forces.append(mfy.get_force_x_reinf(x))
##        plt.plot( x, mfy.get_force_x_reinf( x ) )
##    plt.show()
#    e_arr = orthogonalize([x, P])
#    n_e_arr = [ e / np.max(np.fabs(e)) for e in e_arr ]
#
#    mu_q_arr = np.array(forces)
#
#    n_mu_q_arr = mu_q_arr / np.max(np.fabs(mu_q_arr))
#    m.surf(n_e_arr[0], n_e_arr[1], n_mu_q_arr)

    mfy = MultifilamentYarn()
    mfy.Ll = 10.
    mfy.Lr = 20.

    x = np.linspace(-30, 30, 300)
    P = np.linspace(0., 580, 200)
    forces = []
    for p in P:
        mfy.P = p
        forces.append(mfy.get_force_x_reinf(x))
#        plt.plot( x, mfy.get_force_x_reinf( x ) )
#    plt.show()
    e_arr = orthogonalize([x, P])
    n_e_arr = [ e / np.max(np.fabs(e)) for e in e_arr ]

    mu_q_arr = np.array(forces)

    n_mu_q_arr = mu_q_arr / np.max(np.fabs(mu_q_arr))
    m.surf(n_e_arr[0], n_e_arr[1], n_mu_q_arr)

    m.show()
