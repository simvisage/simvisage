#### \TODO -  adapt this to the current state and add comments

from ibvpy.api import IBVPSolve
#from mfn_line import MFn2DDataGrid
from numpy import frompyfunc, linspace, meshgrid, array, float_, zeros
from math import log, sin, cos
from scipy.linalg import solve, norm
from enthought.traits.api import HasTraits, implements
#from enthought.traits.ui.api import *
#from enthought.traits.ui.menu import *
from ibvpy.api import RTraceGraph, TStepper as TS
from ibvpy.core.tstepper import TStepper
from ibvpy.core.i_tstepper_eval import ITStepperEval
from ibvpy.core.tstepper_eval import TStepperEval as TSE
from time import sleep

from matplotlib import pylab as p

# Include the distinction between init, update, iter calls to the
# corrector/predictor
#
if 0:
    F = array( [ 0.5, 0.1] )
    fn = lambda x, y: sin( x ) * cos( y )
    fn1 = lambda x, y: cos( x ) * cos( y )
    fn2 = lambda x, y: sin( x ) * -sin( y )
    fn11 = lambda x, y:-sin( x ) * cos( y )
    fn12 = lambda x, y: cos( x ) * -sin( y )
    fn21 = lambda x, y: cos( x ) * -sin( y )
    fn22 = lambda x, y: sin( x ) * -cos( y )
    U_0 = array( [8.3, 4.0] )

if 1:
    F = array( [20., 20.] )
    fn = lambda x, y: pow( x, 4 ) * pow( y, 4 ) + 2 * pow( x, 2 ) + x * y + pow( y, 2 )
    fn1 = lambda x, y: 4 * pow( x, 3 ) * pow( y, 4 ) + 2 * 2 * pow( x, 1 ) + y
    fn2 = lambda x, y: pow( x, 4 ) * 4 * pow( y, 3 ) + 2 * pow( y, 1 ) + x
    fn11 = lambda x, y: 12 * pow( x, 2 ) * pow( y, 4 ) + 2 * 2
    fn12 = lambda x, y: 4 * pow( x, 3 ) * 4 * pow( y, 3 ) + 1
    fn21 = lambda x, y: 4 * pow( x, 3 ) * 4 * pow( y, 3 ) + 1
    fn22 = lambda x, y: pow( x, 4 ) * 12 * pow( y, 2 ) + 2
    U_0 = array( [0., 0.] )

if 0:
    F = array( [20., 10.] )
    fn = lambda x, y: 2 * pow( x, 2 ) + x * y + pow( y, 2 )
    fn1 = lambda x, y: 2 * 2 * pow( x, 1 ) + y
    fn2 = lambda x, y: 2 * pow( y, 1 ) + x
    fn11 = lambda x, y: 2 * 2
    fn12 = lambda x, y: 4 * pow( y, 3 ) + 1
    fn21 = lambda x, y: 4 * pow( y, 3 ) + 1
    fn22 = lambda x, y: 2
    U_0 = [0., 0.]

##################
# Equation system
#
dG_n = lambda u: array( [[ fn11( u[0], u[1] ), fn12( u[0], u[1] ) ],
                         [ fn21( u[0], u[1] ), fn22( u[0], u[1] ) ]] )

G_n = lambda u: array( [ fn1( u[0], u[1] ) , fn2( u[0], u[1] )  ] )


class TSFn2D( TSE ):
    '''
    Class managing the single step.
    '''
    implements( ITStepperEval )

    def x_default( self ):
        potential = MFn2DDataGrid( x_min = 0., x_max = 10., n_x = 10,
                                  y_min = 0., y_max = 10., n_y = 10,
                                  ufunc = fn )
        fn_x = MFn2DDataGrid( x_min = 0., x_max = 10., n_x = 10,
                                  y_min = 0., y_max = 10., n_y = 10,
                                  ufunc = fn1 )
        fn_y = MFn2DDataGrid( x_min = 0., x_max = 10., n_x = 10,
                                  y_min = 0., y_max = 10., n_y = 10,
                                  ufunc = fn2 )
        return {'potential' : potential,
                'fn_x' : fn_x,
                'fn_y' : fn_y}

    def setup( self, sctx ):
        pass

    def new_cntl_var( self ):
        global U_0
        return U_0

    def new_resp_var( self ):
        return zeros( 2, float_ )

    def get_corr_pred( self, sctx, u_k, d_u, tn, tn1 ):
        '''
        Additionally handle the constraints
        '''
        dG = dG_n( u_k )
        G = G_n( u_k )
        return G, dG

    def update_state( self, sctx, U ):
        pass

if __name__ == '__main__':

    from ibvpy.api import TLoop, TLine, BCDof

    tl = TLoop( tstepper = TS( tse = TSFn2D() ),
                tline = TLine( min = 0.0, step = 1.0, max = 20.0 ) )
    tl.tstepper.bcond_list = [ BCDof( var = 'f', dof = 0, value = 20.0,
                             time_function = lambda t: 1.0 + 0.33 * t ** 2 ),
                      BCDof( var = 'f', dof = 1, value = 20.0,
                             time_function = lambda t: 1.0 + 0.3 * t ** 2 ) ]
    tl.tstepper.rtrace_list = [ RTraceGraph( name = 'u01 vs. u1',
                                var_x = 'U_k', idx_x = 0,
                                var_y = 'U_k', idx_y = 1,
                                update_on = 'iteration' ),
                      RTraceGraph( name = 'u01 vs. u1',
                                var_x = 'U_k', idx_x = 0,
                                var_y = 'U_k', idx_y = 1,
                                update_on = 'update' )
                      ]

    tl.eval()

    for rt in tl.tstepper.rtrace_list:
        rt.refresh()
        rt.trace.plot( p, 'o-' )
    p.show()
