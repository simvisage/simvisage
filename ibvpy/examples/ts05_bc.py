
from etsproxy.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

from etsproxy.traits.ui.api import \
     Item, View, HGroup, ListEditor, VGroup, Group

from etsproxy.traits.ui.menu import \
     NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
     MenuBar, Separator

from math  import \
     pow, fabs

from numpy import \
     array, zeros, int_, float_, ix_, dot, linspace, hstack, vstack, arange

from scipy.linalg import \
     inv, det

import time

from ibvpy.fets.fets_eval import FETSEval

#-----------------------------------------------------------------------------
# FETS1D2L
#-----------------------------------------------------------------------------

class FETS1D2L(FETSEval):

    debug_on = True
    
    # Dimensional mapping
    dim_slice = slice(0,1)    
    
    n_e_dofs = Int(2)
    n_nodal_dofs = Int(1)
        
    E = Float( 1.0, label = "Young's modulus" )
    A = Float( 1.0, label = "Cross-sectional area" )
    field_entity_type = 'line'
    vtk_r = [[-1.],[1.]]
    
    def get_N_geo_mtx( self, r_pnt ):
        r = r_pnt[0]
        N_mtx=array([[0.5 -r/2.,0.5 + r/2.]])
        return N_mtx
    
    def get_dNr_geo_mtx(self, r_pnt ):
        '''
        Return the matrix of shape function derivatives.
        Used for the conrcution of the Jacobi matrix.

        @TODO - the B matrix is used
        just for uniaxial bar here with a trivial differential
        operator.
        '''
        return array([[-1./2,1./2]])

    def get_N_mtx(r_pnt):
        return self.get_N_geo_mtx()

    def get_dNr_mtx(self,r_pnt):
        '''
        Return the derivatives of the shape functions
        '''
        return self.get_dNr_geo_mtx(r_pnt)
    
    def get_B_mtx( self, r_pnt, X_mtx ):
        J_mtx = self.get_J_mtx(r_pnt,X_mtx)
        dNr_mtx = self.get_dNr_mtx( r_pnt )
        dNx_mtx = dot( inv( J_mtx ), dNr_mtx  )
        return dNx_mtx

#    # @todo: (alex): is this method used?
#    def get_X_mtx( self, X_mtx ):
#        X_mtx = self._get_X_mtx()
#        # make sure that the supplied global coordinate is within the
#        # element domain
#        return X_mtx

    def get_mtrl_corr_pred(self, sctx, eps_mtx):
        D_mtx = array([[ self.E * self.A ]])
        sig_mtx = dot( D_mtx, eps_mtx )
        return sig_mtx, D_mtx

#----------------------- example --------------------

if __name__ == '__main__':
    
    from ibvpy.api import \
        TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval

    # Tseval for a discretized line domain
    #
    tseval  = DOTSEval( fets_eval = FETS1D2L() )
    
    # Discretization
    #
    shape = 3
    domain = MGridDomain( lengths = (1.,0.,0.), 
                             shape = (shape,0,0),
                             n_e_nodes_geo = (1,0,0), 
                             n_e_nodes_dof = (1,0,0), 
                             node_map_geo = [0,1], 
                             node_map_dof = [0,1] )
    
    # Put the tseval (time-stepper) into the spatial context of the
    # discretization and specify the response tracers to evaluate there.
    #
    right_dof = shape
    ts = TS( tse = tseval,
	     sdomain = domain,
	     bcond_list = [ BCDof(var='u', dof = 0, value = 0.),
                     BCDof(var='u', dof = right_dof, value = 1 ),
## 			 BCDof(var='f', dof = right_dof - 1, value = 1 ),
##                          BCDof(var='f', dof = right_dof - 1, value = 0,
##                                link_dofs = [right_dof],
##                                link_coeffs = [1.]),
                     BCDof(var='u', dof = right_dof - 2, value = 0.8,
                           link_dofs = [right_dof-1],
                           link_coeffs = [-1.]),
                         ],
	     rtrace_list = [ RTraceGraph(name = 'Fi,right over u_right' ,
				                        var_y = 'F_int', idx_y = right_dof,
				                        var_x = 'U_k', idx_x = right_dof,
				                        update_on = 'update'),
                     RTraceGraph(name = 'Fi,right over u_right - 1' ,
				              var_y = 'F_int', idx_y = right_dof-1,
				              var_x = 'U_k', idx_x = right_dof-1,
				              update_on = 'update'),
                     RTraceGraph(name = 'Fi,right over u_right - 2' ,
				              var_y = 'F_int', idx_y = right_dof-2,
				              var_x = 'U_k', idx_x = right_dof-2,
				              update_on = 'update'),
                     RTraceDomainField(name = 'Strain Field' ,
                              var = 'eps', idx = 0,
                              update_on = 'update'),
			 ]		     
	     )

    #print ts.elcoord_array
    # Add the time-loop control
    #
    tl = TLoop( tstepper = ts,
	     DT = 1,
	     tline  = TLine( min = 0.0,  max = 1.0 ))

    tl.eval()
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    sim = IS( tloop = tl )
    sim.configure_traits()
    
