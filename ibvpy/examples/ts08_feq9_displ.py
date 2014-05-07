
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
     array, zeros, int_, float_, ix_, dot, linspace, hstack, vstack, arange, \
     identity

from scipy.linalg import \
     inv, det

import time

from ibvpy.api import FETSEval

#-----------------------------------------------------------------------------
# FEBar1D
#-----------------------------------------------------------------------------

class FEQ8U(FETSEval):
    debug_on = True

    # Dimensional mapping
    dim_slice = slice(0, 2)

    n_e_dofs = Int(8 * 2)
    t = Float(1.0, label='thickness')
    E = Float(1.0, label="Young's modulus")
    nu = Float(0.2, label="Poisson's ratio")

    # Integration parameters
    #
    ngp_r = 3
    ngp_s = 3

    field_entity_type = 'quad'
    # 4 corner nodes, 4 edge nodes and 1 interior nodes
    vtk_r = [[-1., -1.],
                        [ 0., -1.],
                        [ 1., -1.],
                        [-1., 0.],
                        [ 0., 0.],
                        [ 1., 0.],
                        [-1., 1.],
                        [ 0., 1.],
                        [ 1., 1.]]
    field_faces = [[0, 1, 4, 3],
                   [1, 2, 5, 4],
                   [3, 4, 7, 6],
                   [4, 5, 8, 7]]

    n_nodal_dofs = Int(2)

    # Order of node positions for the formulation of shape function
    #
    _node_coord_map = Array(Float, (8, 2),
                             [[-1., -1.],
                              [ 1., -1.],
                              [ 1., 1.],
                              [-1., 1.],
                              [ 0., -1.],
                              [ 1., 0.],
                              [ 0., 1.],
                              [-1., 0.]])

    # Order of node positions for the geometry approximation
    #
    # @TODO is it possible to use 'vtk_r' instead ?
    _node_coord_map_geo = Array('float_', (4, 2),
                             [[-1., -1.],
                              [ 1., -1.],
                              [ 1., 1.],
                              [-1., 1.]])


    #---------------------------------------------------------------------
    # Method required to represent the element geometry
    #---------------------------------------------------------------------
    def get_N_geo_mtx(self, r):
        '''
        Return the value of shape functions for the specified local coordinate r
        '''
        cx = self._node_coord_map_geo
        N_geo = array([[ 1 / 4.*(1 + r[0] * cx[i, 0]) * (1 + r[1] * cx[i, 1]) for i in range(0, 4) ]])
        return N_geo

    def get_dNr_geo_mtx(self, r):
        '''
        Return the matrix of shape function derivatives.
        Used for the conrcution of the Jacobi matrix.

        @TODO - the B matrix is used
        just for uniaxial bar here with a trivial differential
        operator.
        '''
        cx = self._node_coord_map_geo
        dNr_geo = array([[ 1 / 4.*cx[i, 0] * (1 + r[1] * cx[i, 1]) for i in range(0, 4) ],
                          [ 1 / 4.*cx[i, 1] * (1 + r[0] * cx[i, 0]) for i in range(0, 4) ]])
        return dNr_geo

    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------
    def get_N_mtx(self, r):
        '''
        Returns the matrix of the shape functions used for the field approximation
        containing zero entries. The number of rows corresponds to the number of nodal
        dofs. The matrix is evaluated for the specified local coordinate r.
        '''
        cx = self._node_coord_map

        N_c = [ 1 / 4.*(1 + r[0] * cx[i, 0]) * (1 + r[1] * cx[i, 1]) * (r[0] * cx[i, 0] + r[1] * cx[i, 1] - 1) for i in range(0, 4) ]

        N_0m1 = 1 / 2.*(1 - r[0] ** 2) * (1 + cx[4, 1] * r[1])
        N_10 = 1 / 2.*(1 + cx[5, 0] * r[0]) * (1 - r[1] ** 2)
        N_01 = 1 / 2.*(1 - r[0] ** 2) * (1 + cx[6, 1] * r[1])
        N_m10 = 1 / 2.*(1 + cx[7, 0] * r[0]) * (1 - r[1] ** 2)

        N = array([ N_c + [ N_0m1, N_10, N_01, N_m10 ]])

        I_mtx = identity(self.n_nodal_dofs, float)
        N_mtx_list = [I_mtx * N[0, i] for i in range(0, N.shape[1])]
        N_mtx = hstack(N_mtx_list)
        return N_mtx

    def get_dNr_mtx(self, r):
        '''
        Return the derivatives of the shape functions used for the field approximation
        '''
        return zeros((2, 8))

    def get_B_mtx(self, r, X):
        J_mtx = self.get_J_mtx(r, X)
        dNr_mtx = self.get_dNr_mtx(r)
        dNx_mtx = dot(inv(J_mtx), dNr_mtx)
        Bx_mtx = zeros((3, 16), dtype='float_')
        for i in range(0, 8):
            Bx_mtx[0, i * 2] = dNx_mtx[0, i]
            Bx_mtx[1, i * 2 + 1] = dNx_mtx[1, i]
            Bx_mtx[2, i * 2] = dNx_mtx[1, i]
            Bx_mtx[2, i * 2 + 1] = dNx_mtx[0, i]
        return Bx_mtx

    def get_mtrl_corr_pred(self, sctx, eps_mtx):
        D_mtx = zeros((3, 3), dtype='float_')
        E = self.E
        nu = self.nu
        D_factor = E / (1 - nu * nu)
        D_mtx[0, 0] = D_factor * 1
        D_mtx[1, 0] = D_factor * nu
        D_mtx[0, 1] = D_factor * nu
        D_mtx[1, 1] = D_factor * 1
        D_mtx[2, 2] = D_factor * (1 - nu) / 2.
        sig_mtx = dot(D_mtx, eps_mtx)
        return sig_mtx, D_mtx

    def get_x_mtx(self, X):
        X = self._get_X_mtx()
        # make sure that the supplied global coordinate is within the
        # element domain


#----------------------- example --------------------

if __name__ == '__main__':
    from core.api import \
        TStepper as TS, MGridDomain, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval

    fets_eval = FEQ8U()

    single_elem = True
    if single_elem:
        from core.melem_domain import MElemDomain, RTraceElemField
        quad_elem = MElemDomain(X=[[0, 0],
                                      [10, 0],
                                      [10, 10],
                                      [0, 10]])

        ts = TS(tse=fets_eval,
                 sdomain=quad_elem,
                 bcond_list=[ BCDof(var='u', dof=i, value=0.) for i in [0, 1, 3] ] +
                           [ BCDof(var='f', dof=i, value=10) for i in [5, 7] ],
                 rtrace_list=[ RTraceElemField(name='Displacement' ,
                                          var='u', idx=0,
                                          record_on='update'),
                             RTraceElemField(name='Deformation' ,
                                          var='eps', idx=0,
                                          record_on='update'),
                             RTraceElemField(name='Jacobi determinant' ,
                                          var='J_det', idx=0,
                                          record_on='update'),
                             RTraceElemField(name='N0' ,
                                          var='N_mtx', idx=0,
                                          record_on='update'),
                             RTraceElemField(name='B0' ,
                                          var='B_mtx0', idx=0,
                                          record_on='update'),
                             RTraceElemField(name='B1' ,
                                          var='B_mtx1', idx=0,
                                          record_on='update'),
                             RTraceElemField(name='B2' ,
                                          var='B_mtx2', idx=0,
                                          record_on='update')
                                                                                                                              ]
                 )

        # Add the time-loop control
        #
        tl = TLoop(tstepper=ts,
                 DT=0.5,
                 tline=TLine(min=0.0, max=0.0))

        tl.eval()
        # Put the whole stuff into the simulation-framework to map the
        # individual pieces of definition into the user interface.
        #
        sim = IS(tloop=tl)
        sim.configure_traits()

    else:
        # Tseval for a discretized line domain
        #
        tseval = DOTSEval(fets_eval=fets_eval)

        # Discretization
        #
        line_domain = MGridDomain(lengths=(3., 1., 0.), shape=(8, 8, 0), n_nodal_dofs=2)

        # Put the tseval (time-stepper) into the spatial context of the
        # discretization and specify the response tracers to evaluate there.
        #
        right_dof = 2
        ts = TS(tse=tseval,
             sdomain=line_domain,
             bcond_list=[ BCDof(var='u', dof=i, value=0.) for i in line_domain.get_left_dofs()[:, 0] ] +
                        [ BCDof(var='u', dof=i, value=0.) for i in line_domain.get_left_dofs()[:, 1] ] +
                        [ BCDof(var='u', dof=i, value=20) for i in line_domain.get_right_dofs()[:, 1] ],
             rtrace_list=[ RTraceGraph(name='Fi,right over u_right (iteration)' ,
                      var_y='F_int', idx_y=right_dof,
                      var_x='U_k', idx_x=right_dof,
                      record_on='update'),
#                      RTraceDomainField(name = 'Flux field' ,
#                      var = 'eps', idx = 0,
#                      record_on = 'update'),
                      RTraceDomainField(name='Temperature' ,
                      var='u', idx=0,
                      record_on='update'),
#                             RTraceDomainField(name = 'N0' ,
#                                          var = 'N_mtx', idx = 0,
#                                          record_on = 'update')

                 ]
             )

        # Add the time-loop control
        #
        tl = TLoop(tstepper=ts,
             DT=0.5,
             tline=TLine(min=0.0, max=1.0))

        tl.eval()
        # Put the whole stuff into the simulation-framework to map the
        # individual pieces of definition into the user interface.
        #
        sim = IS(tloop=tl)
        sim.configure_traits()
