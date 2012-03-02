'''
Created on Aug 5, 2010

@author: rostislav
'''

from enthought.traits.api import \
    HasTraits, Float, Str, implements, Range

from math import pi, e

from numpy import sign, linspace, array, cos, sqrt, argmax, hstack, max

from matplotlib import pyplot as plt

from stats.spirrid.i_rf import \
    IRF

from stats.spirrid.rf import \
    RF

def Heaviside( x ):
    return sign( sign( x ) + 1.0 )

class ConstantFrictionFreeLengthFiniteFiber( RF ):
    '''
    '''

    implements( IRF )

    title = Str( 'pull-out with constant friction and free length ' )

    tau = Float( 0.5e6, auto_set = False, enter_set = True,
                distr = ['uniform', 'weibull_min', 'norm'] )

    # free length
    l = Float( 0.01, auto_set = False, enter_set = True,
                distr = ['uniform', 'norm'] )

    E = Float( 70e9, auto_set = False, enter_set = True,
               distr = ['uniform'] )

    A = Float( 5.30929158457e-10, auto_set = False, enter_set = True,
               distr = ['uniform', 'weibull_min'] )


    # waviness in strains
    slack = Float( 0.0, auto_set = False, enter_set = True,
               distr = ['uniform'] )

    # embedded length
    L = Float( 0.03, auto_set = False, enter_set = True,
               distr = ['uniform'] )

    # breking stress
    sigma_fu = Float( 1200.e6, auto_set = False, enter_set = True,
               distr = ['uniform'] )

    u = Float( auto_set = False, enter_set = True,
               ctrl_range = ( 0.0, 1.0, 10 ) )

    # 

    def __call__( self, u, tau, l, E, A, slack, L, sigma_fu ):
        ''' method for evaluating the u-F diagram for
        pullout with constant friction at the interface and
        free fiber length sticking out of the matrix '''
        # defining tau as length dependent
        tau = tau * sqrt( 4 * A * pi )
        # constitutive law for a pullout without free length 
        q = ( -l * ( 1 + slack ) * tau * Heaviside( u / 2 - l * ( slack ) ) + \
                + sqrt( ( l * ( 1 + slack ) * tau ) ** 2 * Heaviside( u / 2 - l * ( slack ) ) \
                + 2 * E * A * ( u / 2 - ( l * slack ) ) * Heaviside( u / 2 - l * ( slack ) ) ) )

        # deformation of the free length added
        continuous = q * Heaviside( L - l - q / tau )\
                + ( L - l ) * tau * Heaviside( l + q / tau - L )

        # a check whether the fiber has reached the breaking stress        
        return continuous * Heaviside( sigma_fu - continuous / A )


if __name__ == '__main__':
    cf = ConstantFrictionFreeLengthFiniteFiber()
    u = linspace( 0, 0.8, 200 )
    q = cf( u, 0.5e6, 0.01, 70.e9, 0.89e-6, 0.0, 0.1, 300.e6 )
    plt.plot( u, q, linewidth = 1 )
    plt.show()
