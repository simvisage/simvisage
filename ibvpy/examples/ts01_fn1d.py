
from numpy import frompyfunc, linspace, meshgrid, array, float_, zeros
from math import log, sin, cos
from scipy.linalg import solve, norm

from traits.api import Array, Bool, Callable, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, Callable, List, \
    TraitDict, Any, implements, WeakRef
from traitsui.api import Item, View, HGroup, ListEditor, VGroup

from ibvpy.api import RTraceGraph
from ibvpy.core.tstepper import TStepper as TS
from ibvpy.core.i_tstepper_eval import ITStepperEval
from ibvpy.core.tstepper_eval import TStepperEval
from time import sleep
from math import pi
from matplotlib import pylab as p
from envisage.ui.workbench.api import WorkbenchApplication


class TSFn1D(TStepperEval):

    '''
    Class managing the single step.
    '''
    implements(ITStepperEval)

    def new_cntl_var(self):
        return zeros(1, float_)

    def new_resp_var(self):
        return zeros(1, float_)

#    def get_corr_pred( self, sctx, u, d_u, tn, tn1 ):
#        '''
#        Get the corrector and predictor.
#        '''
#        sleep( 0.1 )
#
#        x = u[0]
#        G_n = array( [ x * sin( x ) ] )
#        dG_n = array( [[ sin( x ) + x * cos( x ) ]] )
#
#        return G_n, dG_n

    def get_corr_pred(self, sctx, u, d_u, tn, tn1):
        '''
        Get the corrector and predictor.

        Implicit tracing of the quadratic curve:

        x = linspace( 0, 2.0, 100 )
        y = -( x - 1 ) ** 2 + 1
        p.plot( x, y )
        p.show()
        '''
        x = u[0]
        G_n = array([-(x - 1) ** 2 + 1], dtype=float)
        # dG_n = array( [[ -2 * ( 0 - 1 )  ]], dtype = float )
        dG_n = array([[-2 * (x - 1)]], dtype=float)

        return G_n, dG_n

if __name__ == '__main__':

    from ibvpy.api import TLoop, TLine, BCDof, \
        IBVPSolve as IS

    tl = TLoop(tstepper=TS(tse=TSFn1D()),
               tolerance=1e-8,
               # debug = True,
               KMAX=4,
               verbose_iteration=True,
               tline=TLine(min=0.0, step=0.25, max=1.0))
    # tl.tstepper.bcond_list = [ BCDof( var = 'u', dof = 0, value = 1.0 ) ]
    tl.tstepper.bcond_list = [BCDof(var='f', dof=0, value=0.9)]
    tl.tstepper.rtrace_list = [
        RTraceGraph(name='u_update',
                    var_x='U_k', idx_x=0,
                    var_y='F_int', idx_y=0,
                    record_on='update'),
        RTraceGraph(name='u_iter',
                    var_x='U_k', idx_x=0,
                    var_y='F_int', idx_y=0,
                    record_on='iteration')
    ]

    tl.eval()

    for rt in tl.tstepper.rtrace_list:
        rt.refresh()
        rt.trace.plot(p, 'o-')
    p.show()
