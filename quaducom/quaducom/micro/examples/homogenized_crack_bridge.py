'''
Created on Sep 30, 2011

@author: rostar
'''

from quaducom.resp_func.cb_clamped_fiber import CBClampedFiberSP
from stats.spirrid.spirrid import SPIRRID, orthogonalize
from stats.spirrid.rv import RV
from matplotlib import pyplot as plt
import numpy as np
from mayavi import mlab as m

if __name__ == '__main__':

    w = np.linspace( 0, 2., 100 )
    x = np.linspace( -60, 110, 100 )

    cbcsp = CBClampedFiberSP()

    s = SPIRRID( q = cbcsp,
         sampling_type = 'LHS',
         implicit_var_eval = False,
         evars = dict( w = w,
                       x = x,
                        ),
         tvars = dict( tau = 0.1,
                         l = 3.5,
                         D_f = 26e-3,
                         E_f = RV( 'uniform', loc = 65e3, scale = 10e3 ),
                         theta = RV( 'uniform', loc = 0.0, scale = 0.1 ),
                         xi = RV( 'weibull_min', scale = 0.02, shape = 5 ),
                         phi = 1.,
                         Ll = 20.,
                         Lr = 80. ),
         n_int = 20 )

    mu_q_arr = np.array( s.mu_q_arr )

    def mlab3Dview():

        e_arr = orthogonalize( [x, w] )
        n_e_arr = [ e / np.max( np.fabs( e ) ) for e in e_arr ]

        n_mu_q_arr = mu_q_arr / np.max( np.fabs( mu_q_arr ) )
        m.surf( n_e_arr[0], n_e_arr[1], n_mu_q_arr )
        m.show()

    def w_P():
        plt.plot( w, mu_q_arr[:, np.argwhere( x >= 0 )[0]], color = 'black', lw = 2 )
        plt.title( 'force in reinforcement at the crack plane: x = 0' )
        plt.show()

    def x_P():
        plt.plot( x, mu_q_arr[20, :], lw = 2, label = 'w = %0.2f mm' % w[20] )
        plt.plot( x, mu_q_arr[60, :], lw = 2, label = 'w = %0.2f mm' % w[40] )
        plt.plot( x, mu_q_arr[80, :], lw = 2, label = 'w = %0.2f mm' % w[80] )
        plt.title( 'force profile along reinforcement at various crack openings' )
        plt.legend()
        plt.show()

    mlab3Dview()
    w_P()
    x_P()
