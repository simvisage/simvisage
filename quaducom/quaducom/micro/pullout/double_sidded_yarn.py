'''
Created on 3.12.2010

Response function for a double sided pullout of a continuous filament
with different friction at both sides.

@author: Q
'''

from etsproxy.traits.api import \
    Float, Str, implements

from numpy import sign, linspace

from matplotlib import pyplot as plt

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf import \
    RF

def H( x ):
    return ( sign( x ) + 1.0 ) / 2.0


class DoublePullout( RF ):

    implements( IRF )

    title = Str( 'double yarn pullout' )

    xi = Float( 0.014, auto_set = False, enter_set = True,
                distr = ['weibull_min'] )
    tau1 = Float( 6., auto_set = False, enter_set = True,
                distr = ['uniform', 'norm'] )
    tau2 = Float( 6., auto_set = False, enter_set = True,
                distr = ['uniform', 'norm'] )
    L = Float( 1., auto_set = False, enter_set = True,
              distr = ['uniform'] )
    l = Float( 0.01, auto_set = False, enter_set = True,
              distr = ['uniform'] )
    A = Float( 1., auto_set = False,
              enter_set = True, distr = ['uniform', 'weibull_min'] )
    E_mod = Float( 1., auto_set = False, enter_set = True,
                  distr = ['uniform'] )
    theta = Float( 0.01, auto_set = False, enter_set = True,
                  distr = ['uniform', 'norm'] )

    w = Float( auto_set = False, enter_set = True,
               ctrl_range = ( 0.0, 1.0, 10 ) )

    C_code = '''
        '''

    def __call__( self, w, tau1, tau2, L, l, A, E_mod, theta, xi ):
        l = l * ( 1 + theta )
        w = w - theta * l
        P_ = ( -tau1 * tau2 * l + ( tau1 * tau2 * ( tau2 * tau1 * l ** 2 + 2 * w * H( w ) * E_mod * A * tau2 + 2 * w * H( w ) * E_mod * A * tau1 ) ) ** ( 1. / 2. ) ) / ( tau2 + tau1 )
        P = P_ * H( A * E_mod * xi - P_ )
        return P

if __name__ == '__main__':
    dp = DoublePullout()
    X = linspace( 0, 0.0025, 100 )
    Y = dp( X, 6.0, 6.0, 1., 0.01, .89e-6 / 1600, 72.0e9, 0.01, 0.014 )
    plt.plot( X, Y, linewidth = 2 )
    plt.show()
