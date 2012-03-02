'''
Created on Aug 6, 2010

@author: rostislav
'''

from enthought.traits.api import \
    Float #@UnresolvedImport
from numpy import sqrt, exp, sum
from math import pi
#from enthought.mayavi import mlab

from stats.spirrid import SPIRRID
from matplotlib import pyplot as plt
from const_frict_free_length_finite_fiber import ConstantFrictionFreeLengthFiniteFiber

class CorSPIRRID( SPIRRID ):
    '''class evaluating the mean value of a function with
    two normally distributed variables correlated with
    the factor rho '''

    # correlation factor
    rho = Float( 0.99 )

    def _get_pdf_theta_grid( self ):
        ''' returns the joint correlated pdf '''
        mu_x = self.rv_dict['l'].pd.distr_type.mean
        sig_x = sqrt( self.rv_dict['l'].pd.distr_type.variance )

        mu_y = self.rv_dict['tau'].pd.distr_type.mean
        sig_y = sqrt( self.rv_dict['tau'].pd.distr_type.variance )

        x_theta, y_theta = self.theta_ogrid

        # evaluating the theta difference
        dx = x_theta[1, 0] - x_theta[0, 0]
        dy = y_theta[0, 1] - y_theta[0, 0]

        # normalization
        x = ( x_theta - mu_x ) / sig_x
        y = ( y_theta - mu_y ) / sig_y

        # dx and dy should be halved at around the grid region 
        return self.stdPDF_NOR_bivar( x, y ) * dx * dy / sig_x / sig_y

    def stdPDF_NOR_bivar( self, x, y ):
        rho = self.rho
        # substitution 1
        Z = -2. * rho * x * y + x * x + y * y
        # substitution 2
        X = 1. - ( rho ** 2 )
        return exp( -0.5 * Z / X ) / ( sqrt( X ) * 2. * pi );


def run():
    # Quantities for the response function
    # and randomization

    # construct a default response function for a single filament

    rf = ConstantFrictionFreeLengthFiniteFiber( tau = 0.5e6,
                                               l = 0.01,
                                               E = 70e9,
                                               A = 0.89e-6,
                                               slack = 0.0,
                                               L = 0.3,
                                               sigma_fu = 200e6 )

    # construct the integrator and provide
    # it with the response function.

    max_u = 1.8

    s = SPIRRID( rf = rf,
                 min_eps = 0.00, max_eps = max_u, n_eps = 80,
                 cached_dG = True, compiled_QdG_loop = False,
                compiled_eps_loop = False )

    # construct the correlated integrator and
    # provide it with the response function.

    sc = CorSPIRRID( rf = rf, rho = -0.99,
                 min_eps = 0.00, max_eps = max_u, n_eps = 80,
                 cached_dG = True, compiled_QdG_loop = False,
                compiled_eps_loop = False )

    n_int = 100

    s.add_rv( 'tau', distribution = 'norm', loc = 0.3e6, scale = 0.1e6, n_int = n_int )
    s.add_rv( 'l', distribution = 'norm', loc = 0.03, scale = 0.01, n_int = n_int )

    sc.add_rv( 'tau', distribution = 'norm', loc = 0.3e6, scale = 0.1e6, n_int = n_int )
    sc.add_rv( 'l', distribution = 'norm', loc = 0.03, scale = 0.01, n_int = n_int )

    # define a tables with the run configurations to start in a batch
    print 'sum', sum( sc.dG_grid )
    print 'sum', sum( s.dG_grid )

#    z = sc.pdf_theta_grid
#    x, y = mgrid[-4:4:50j, -4:4:50j]
#
#    mlab.surf(x, y, z)
#    mlab.show()

    s.mean_curve.plot( plt, color = 'red', linewidth = 2, label = 'spirrid' )
    sc.mean_curve.plot( plt, color = 'blue', linewidth = 2, label = 'cor spirrid' )

    plt.xlabel( 'strain [-]' )
    plt.ylabel( 'stress' )
    plt.legend( loc = 'lower right' )

    plt.title( s.rf.title )
    plt.show()

if __name__ == '__main__':
    run()
