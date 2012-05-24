'''
Created on Jun 16, 2010

@author: rostislav
'''
'''
Created on May 27, 2010

@author: rostislav
'''

from etsproxy.traits.api import HasTraits, Float, Str, implements, \
    Float, cached_property, Property, Array, Instance
from matplotlib import pyplot as plt
from numpy import linspace, array, frompyfunc
from stats.spirrid import SPIRRID
from stats.spirrid.i_rf import IRF
from stats.spirrid.rf import RF

from stats.spirrid.rf_filament import \
    Filament

class Yarn( RF ):
    ''' rf for the outer integration as a spirrid call'''
    implements( IRF )

    title = Str( 'yarn statistics' )

    loc_theta = Float( auto_set = False,
                       enter_set = True,
                       distr = ['uniform'] )

    rf = Instance( Filament )
    def _rf_default( self ):
        return Filament( E_mod = 70.e9, xi = 0.02,
                         A = 5.30929158457e-10, theta = 0, lambd = 0 )

    # construct the spirrid object outside the called method so that
    # it is constructed only once
    s_inner = Property( Instance( SPIRRID ), depends_on = 'rf, rf.changed', )
    @cached_property
    def _get_s_inner( self ):

        s = SPIRRID( rf = self.rf,
                     cached_dG = True,
                     compiled_QdG_loop = False,
                     compiled_eps_loop = False )

        # construct the random variables
        n_int = 10

        s.add_rv( 'xi', distribution = 'weibull_min', scale = 0.02,
                  shape = 10., n_int = n_int )
        s.add_rv( 'theta', distribution = 'uniform', loc = self.loc_theta,
                  scale = 0.01, n_int = n_int )
        #.add_rv('lambd', distribution = 'uniform', loc = 0.0, scale = 0.01, n_int = n_int)

        return s

    def __call__( self, eps, loc_theta ):
        ''' the inner integral method must be vectorized using
        frompyfunc because the pdistrib instances used in spirrid
        don't work with arrays as distribution parameters '''
        func = frompyfunc( self._resp, 2, 1 )
        return func( eps, loc_theta )

    def _resp( self, eps, loc_theta ):
        '''private method implementing the response
        '''
        # integrator instance
        s = self.s_inner

        # add the control variable array cv - CRIMINAL HACK [rch] !!!!
        s.cv = array( [eps] ).flatten()

        r = s.results[0]
        print 'input', s.cv, 'output', r

        return r

if __name__ == '__main__':

    rf = Yarn()

    # control variable - strains
    # can be a float number or array
    cv = linspace( 0.0, 0.06, 10 )

    # construct the outer integrator and provide it with the response function.
    s_outer = SPIRRID( rf = rf,
                 cv = cv,
                 cached_dG = True,
                 compiled_QdG_loop = False,
                 compiled_eps_loop = False )

    # construct the random variables

    n_int = 3

    s_outer.add_rv( 'loc_theta', distribution = 'uniform', loc = 0.00,
                    scale = 0.01, n_int = n_int )

    # plot the mean filament and mean yarn response
    plt.plot( cv, rf( cv, 0.01 ), linewidth = 2,
              label = 'random filament properties' )
    plt.plot( cv, s_outer.results[0], linewidth = 2,
              label = 'random bundle properties' )

    plt.legend( loc = 'best' )
    plt.title( rf.title )
    plt.show()
