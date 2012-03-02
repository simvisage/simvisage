
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

from enthought.traits.ui.api import \
     Item, View, HGroup, ListEditor, VGroup, Group

from enthought.traits.ui.menu import \
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
    
    def get_B_mtx( self, r, X ):
        J_mtx = self.get_J_mtx(r,X)
        dNr_mtx = self.get_dNr_mtx( r )
        dNx_mtx = dot( inv( J_mtx ), dNr_mtx  )
        return dNx_mtx

        
    def get_x_mtx( self, X ):
        X = self._get_X_mtx()
        # make sure that the supplied global coordinate is within the
        # element domain

    def get_mtrl_corr_pred(self, sctx, eps_mtx):
        D_mtx = array([[ self.E * self.A ]])
        sig_mtx = dot( D_mtx, eps_mtx )
        return sig_mtx, D_mtx


#----------------------- example --------------------

if __name__ == '__main__':
    from ibvpy.api import \
        TStepper as TS, MGridDomain, RTraceGraph, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
    
    # Tseval for a discretized line domain
    #
    tseval  = DOTSEval( fets_eval = FETS1D2L() )
    
    from ibvpy.rtrace.rt_domain_field import MeshGridAdaptor

    # Define a mesh domain adaptor as a cached property to 
    # be constracted on demand
    mgrid_adaptor = MeshGridAdaptor( n_nodal_dofs = 1,
                                     n_e_nodes_geo = (1,0,0), 
                                     n_e_nodes_dof = (1,0,0), 
                                     node_map_geo = [0,1], 
                                     node_map_dof = [0,1] )
    
    # Discretization
    #
    line_domain = MGridDomain( lengths = (1.,0.,0.),
                                  shape = (10,0,0),
                                  adaptor = mgrid_adaptor )
    
    # Put the tseval (time-stepper) into the spatial context of the
    # discretization and specify the response tracers to evaluate there.
    #
    right_dof = 10
    ts = TS( tse = tseval,
	     sdomain = line_domain,
	     bcond_list = [ BCDof(var='u', dof = 0, value = 0.),
                     BCDof(var='f', dof = 10, value = 20 ) ],
	     rtrace_list = [ RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
				              var_y = 'F_int', idx_y = right_dof,
				              var_x = 'U_k', idx_x = right_dof,
				              update_on = 'update')  			 ]		     
	     )

    # Add the time-loop control
    #
    tl = TLoop( tstepper = ts,
	     DT = 3.7,
	     tline  = TLine( min = 0.0,  max = 60.0 ))
    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    sim = IS( tloop = tl )
    sim.configure_traits()
    
