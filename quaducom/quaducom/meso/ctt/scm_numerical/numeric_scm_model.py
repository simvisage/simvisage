'''
Created on 10.03.2011

@author: rrypl
'''
from enthought.traits.api import \
    Instance, Enum, Bool, on_trait_change, Int, Event, Array, Tuple, List
from enthought.traits.ui.api import \
    View, Item, VGroup, HGroup, ModelView, HSplit, VSplit
from enthought.traits.ui.menu import OKButton

from enthought.traits.api import HasTraits, Float, Property, \
                                cached_property, Range, Button
from enthought.traits.ui.api import View, Item, Tabbed, VGroup, \
                                VSplit, Group
from enthought.traits.ui.menu import OKButton, CancelButton
from matplotlib import pyplot as plt

from numpy import linspace, frompyfunc, max, abs, array, hstack, sign, argmax, mean, ones
from numpy.random import rand
from scipy.stats import weibull_min
from scipy.special import gamma
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline

from quaducom.meso.ctt.scm_cuypers.reinf_cross_section import \
    SimplyRatio, GridReinforcement
from util.traits.either_type import EitherType
from quaducom.ctt.scm_cuypers.analyt_scm_model import SCM

from math import cos, pi as Pi, e

from stats.spirrid.spirrid_nd import SPIRRID
from quaducom.crackbridge.crack_bridge import StressInFiberWithConstantFriction

def H( x ):
    return sign( sign( x ) + 1. )

class CrackSimulation( SCM ):

    length = Float( 1000., auto_set = False, enter_set = True, # [mm]
                 desc = 'total specimen length', modified = True )

    nx = Int( 1000, auto_set = False, enter_set = True,
                 desc = 'number of length discretization points', modified = True )

    n_cdf = Int( 50, auto_set = False, enter_set = True,
                 desc = 'number of CDF discretization points', modified = True )

    applied_stress = Range( low = 1e-10, high = 50.0, value = 1.2,
                         auto_set = False, enter_set = True, # [MPa]
                 desc = 'current stress', ctrl_param = True )

    # minimal length at which the variable values are not correlated any more
    l_rho = Float( 3., auto_set = False, enter_set = True, # [mm]
                 desc = 'autocorrelation length', modified = True )

    cracks = List

    result_arr = Tuple( Array, Array, Array )

    current_stress = Float

    x_arr = Property( Array, depends_on = 'length,nx' )
    @cached_property
    def _get_x_arr( self ):
        '''discretizes the specimen length'''
        return linspace( 0, self.length, self.nx )

    CDF_arr = Property( Array, depends_on = 'n_cdf' )
    @cached_property
    def _get_CDF_arr( self ):
        '''discretizes the CDF for equal spaced stress values'''
        return weibull_min( self.m, scale = self.scale_sigma_m ).cdf( self.Pf_sigma_arr )

    # TODO linspace stresses
    Pf_sigma_arr = Property( Array, depends_on = 'n_cdf' )
    @cached_property
    def _get_Pf_sigma_arr( self ):
        '''equivalently discretizes the matrix stresses'''
        low = weibull_min( self.m, scale = self.scale_sigma_m ).ppf( 1e-10 )
        high = weibull_min( self.m, scale = self.scale_sigma_m ).ppf( 1. - 1e-10 )
        return linspace( low, high, self.n_cdf )

    sigma_m_ff = Property( depends_on = '+ctrl_param,+modified, reinf_ratio.+modified' )
    @cached_property
    def _get_sigma_m_ff( self ):
        '''stress in the matrix in an uncracked composite'''
        return self.applied_stress * self.E_m / self.E_c

    # TODO three parameter Weibull?
    def surv_arr( self, sigma ):
        '''survival probabilities for stresses along the specimen'''
        return 1 - weibull_min( self.m, scale = self.scale_sigma_m ).cdf( sigma )

    sigma_m_x = Property( Array, depends_on = 'cracks, +modified, reinf_ratio.+modified, +ctrl_param' )
    @cached_property
    def _get_sigma_m_x( self ):
        '''creates the matrix stress array given the current loading and an array of crack positions'''
        sigma_x = array( len( self.x_arr ) * [self.sigma_m_ff] )

        if len( self.result_arr ) > 0:
            for crack in self.cracks:
                cr_arr = self.x_sig( abs( self.x_arr - crack ), self.applied_stress )
                mask = sigma_x <= cr_arr
                sigma_x = sigma_x * mask + cr_arr * ( mask == False )
        return sigma_x

    spirrid_response = Property( depends_on = 'tau' )
    @cached_property
    def _get_spirrid_response( self ):

        rf = StressInFiberWithConstantFriction( u = 0.05 )
        s = SPIRRID( rf = rf,
                 cached_dG = True,
                 compiled_QdG_loop = False,
                 compiled_eps_loop = False
                 )

        n_int = 50
        s.add_rv( 'tau', distribution = 'uniform', loc = 0.1, scale = 10.0, n_int = n_int )
        #s.add_rv( 'l', distribution = 'uniform', loc = 0.0, scale = 5.0, n_int = n_int )
        #s.add_rv( 'D_f', distribution = 'uniform', loc = 25.0e-3, scale = 5.0e-3, n_int = n_int )

        eps_list = s.eps_list
        mu_q = s.mu_q_grid
        return eps_list, mu_q

    def x_sig( self, x, sigma ):
        eps_list, mu_q = self.spirrid_response
        eps_sig = InterpolatedUnivariateSpline( mu_q[0, :], eps_list[1] )
        if max( mu_q ) > sigma:
            pass
        else:
            raise ValueError( 'applied stress higher than the maximum in micromechanical evaluation of a CB' )
        eps = eps_sig( sigma )
        spline = RectBivariateSpline( eps_list[0], eps_list[1], mu_q )
        sigma_f = spline.ev( x, ones( len( x ) ) * eps ) / self.V_f
        sigma_m = ( sigma - sigma_f * self.V_f ) / self.V_m
        return sigma_m

    random_field = Property( Array, depends_on = 'l_rho, length, m, sigma_mu' )
    @cached_property
    def _get_random_field( self ):
        '''generates an array of random matrix strength'''

        np = int( self.length / self.l_rho )
        if np <= 3:
            np = 4
        print np
        x_not_corr_points = linspace( 0, np * self.l_rho, np)
        y_not_corr_points = weibull_min( self.m, scale = self.scale_sigma_m ).ppf( rand( np ) )
        spline = InterpolatedUnivariateSpline( x_not_corr_points, y_not_corr_points )
        return spline( self.x_arr )

    def crack_devel( self ):
        '''checks the matrix stress profile and compares it to the matrix CDF; adds new cracks'''
        x = self.x_arr

        # loading
        if self.sigma_m_ff > self.current_stress:
            self.current_stress = self.sigma_m_ff
            # find the points, where the strength is lower than the stress
            while sum( self.sigma_m_x >= self.random_field ) > 0.:
                cr_pos = argmax( self.sigma_m_x - self.random_field )
                self.cracks.append( self.x_arr[cr_pos] )

        # unloading           
        else:
            self.current_stress = self.sigma_m_ff

        eps_m = self.sigma_m_x / self.E_m
        eps_f = ( ( 1 + self.alpha ) * array( len( x ) * [self.applied_stress / self.E_c] ) - eps_m * self.alpha )

        self.sig_eps[0].append( mean( eps_f ) )
        self.sig_eps[1].append( self.applied_stress )
        self.result_arr = ( x, eps_m, eps_f )

    @on_trait_change( '+modified, reinf_ratio.+modified' )
    def reset_history( self ):
        ''' if some params are changed, the sigma - eps history is not valid any more and is reseted'''
        self.cracks = []
        self.no_of_cracks = ( [0.], [0.] )
        self.sig_eps = ( [0.], [0.] )
        self.crack_devel()

    # for plotting the current Pf value in the CDF
    CDF_p = Property( depends_on = '+ctrl_param,+modified, reinf_ratio.+modified' )
    def _get_CDF_p( self ):
        return self.applied_stress * self.E_m / self.E_c, \
                weibull_min( self.m, scale = self.scale_sigma_m ).cdf( self.applied_stress * self.E_m / self.E_c )

    # for plotting the current Pf value in the CDF
    CDF = Property( depends_on = '+ctrl_param,+modified, reinf_ratio.+modified' )
    def _get_CDF( self ):
        sigma = linspace( 0, 30, 200 )
        return sigma, weibull_min( self.m, scale = self.scale_sigma_m ).cdf( sigma )

    sig_eps = Tuple( List( [0.0] ), List( [0.0] ) )


    no_of_cracks = Tuple( List, List )
    def get_no_of_cracks( self ):
        self.no_of_cracks[0].append( self.current_stress )
        self.no_of_cracks[1].append( len( self.cracks ) )


    def launch( self, max, points ):
        self.reset_history()
        for stress in linspace( 0.01, max, points ):
            self.applied_stress = stress


    traits_view = View( 
                          Group( 
                            Tabbed( 
                              VGroup( 
                                   Group( Item( 'E_f', resizable = False,
                                               label = 'E-modulus',
                                               tooltip = "Young's modulus of the fiber" ),

                                         Item( 'sigma_fu', resizable = False,
                                               label = 'strength',
                                               help = "Strength of the fibers"
                                         ),
                                         label = 'fibers',
                                    ),

                                   Group( 
                                        Item( 'l_rho', resizable = False,
                                               label = 'autocorrelation length' ),
                                       Item( 'E_m', resizable = False,
                                             label = 'E-modulus',
                                             help = "Young's modulus of the matrix" ),
                                       Item( 'sigma_mu', resizable = False,
                                             label = 'strength',
                                             help = "Scale parameter of the matrix strength"
                                             'roughly corresponding to the mean strength' ),
                                       Item( 'm', resizable = False,
                                             label = 'Weibull-modulus',
                                             help = "Weibull modulus of the matrix strength distribution"
                                             'defining the scatter of the strength' ),
                                        label = 'matrix',
                                        scrollable = False,
                                        ),
                                   label = 'Components',
                                   dock = 'tab',
                                   id = 'scm.model.component_params',
                               ),
                                   VGroup( 
                                       Item( 'tau', resizable = False, springy = True ),
                                       Item( 'r', resizable = False, springy = False ),
                                       springy = True,
                                       label = 'Bond',
                                       dock = 'tab',
                                       id = 'scm.model.params',
                                    ),
                                    id = 'scm.model.allparams',
                                  ),
                               VGroup( 
                                   Item( 'reinf_ratio@', show_label = False, resizable = True ),
                                   label = 'Cross section parameters',
                                   dock = 'tab',
                                   id = 'scm.model.reinf_ratio',
                               ),
                               VGroup( 
                                    Item( 'orientation', label = 'fiber orientation' ),
                                    ),
                               id = 'scm.model.splitter',
                               springy = False,
                               layout = 'split',
                            ),
                            id = 'scm.model',
                            dock = 'fixed',
                            scrollable = True,
                            resizable = True,
                            buttons = [OKButton],
                            height = 0.8, width = 0.8
                                   )
