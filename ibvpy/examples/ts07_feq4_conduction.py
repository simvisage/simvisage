

from numpy import \
    array, zeros, int_, float_, ix_, dot, linspace, hstack, vstack, arange
from scipy.linalg import \
    inv
from traits.api import \
    Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
    Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
    on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, \
    Dict

from ibvpy.fets.fets_eval import FETSEval, RTraceEvalElemFieldVar


#-----------------------------------------------------------------------------
# FEQ4
#-----------------------------------------------------------------------------
class FEQ4T(FETSEval):

    debug_on = True

    # Material parameters
    #
    k = Float(1.0, label='conductivity')

    # Dimensional mapping
    #
    dim_slice = slice(0, 2)

    # System mapping parameters
    #
    n_e_dofs = Int(4)
    n_nodal_dofs = Int(1)

    # Order of node positions for the formulation of shape function
    # (isoparametric formulation)
    # Order of node positions for the formulation of shape function
    #
    dof_r = [[-1, -1], [1, -1], [1, 1], [-1, 1]]
    geo_r = [[-1, -1], [1, -1], [1, 1], [-1, 1]]

    # Integration parameters
    #
    ngp_r = 2
    ngp_s = 2

    # Field visualization attributes
    #
    field_entity_type = 'quad'

    # 4 corner nodes, 4 edge nodes and 1 interior nodes
    vtk_r = [[-1., -1.],
             [0., -1.],
             [1., -1.],
             [-1., 0.],
             [0., 0.],
             [1., 0.],
             [-1., 1.],
             [0., 1.],
             [1., 1.]]
    vtk_cells = [[0, 1, 4, 3],
                 [1, 2, 5, 4],
                 [3, 4, 7, 6],
                 [4, 5, 8, 7]]
    vtk_cell_types = 'Quad'

    #---------------------------------------------------------------------
    # Method required to represent the element geometry
    #---------------------------------------------------------------------
    def get_N_geo_mtx(self, r_pnt):
        '''
        Return the value of shape functions for the specified local coordinate r
        '''
        cx = array(self.geo_r, dtype='float_')
        Nr = array([[1 / 4. * (1 + r_pnt[0] * cx[i, 0])
                     * (1 + r_pnt[1] * cx[i, 1]) for i in range(0, 4)]])
        return Nr

    def get_dNr_geo_mtx(self, r_pnt):
        '''
        Return the matrix of shape function derivatives.
        Used for the conrcution of the Jacobi matrix.

        @TODO - the B matrix is used
        just for uniaxial bar here with a trivial differential
        operator.
        '''
        cx = array(self.geo_r, dtype='float_')
        dNr_geo = array([[1 / 4. * cx[i, 0] * (1 + r_pnt[1] * cx[i, 1]) for i in range(0, 4)],
                         [1 / 4. * cx[i, 1] * (1 + r_pnt[0] * cx[i, 0]) for i in range(0, 4)]])
        return dNr_geo

    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------
    def get_N_mtx(self, r_pnt):
        return self.get_N_geo_mtx(r_pnt)

    def get_dNr_mtx(self, r_pnt):
        '''
        Return the derivatives of the shape functions
        '''
        return self.get_dNr_geo_mtx(r_pnt)

    def get_B_mtx(self, r_pnt, X_mtx):
        J_mtx = self.get_J_mtx(r_pnt, X_mtx)
        dNr_mtx = self.get_dNr_mtx(r_pnt)
        dNx_mtx = dot(inv(J_mtx), dNr_mtx)
        return dNx_mtx

    def get_mp_state_array_size(self, sctx):
        return 0

    def setup(self, sctx):
        return

    def get_mtrl_corr_pred(self, sctx, eps_eng, d_eps_eng, tn, tn1, eps_avg=None):
        D_mtx = array([[self.k, 0], [0, self.k]])
        sig_mtx = dot(D_mtx, eps_eng)
        return sig_mtx, D_mtx

    rte_dict = Trait(Dict)

    def _rte_dict_default(self):
        '''
        RTraceEval dictionary with standard field variables.
        '''
        rte_dict = self._debug_rte_dict()
        rte_dict.update({'eps_app': RTraceEvalElemFieldVar(eval=self.get_eps_mtx33),
                         'eps0_app': RTraceEvalElemFieldVar(eval=self.get_eps0_mtx33),
                         'eps1t_app': RTraceEvalElemFieldVar(eval=self.get_eps1t_mtx33),
                         'u': RTraceEvalElemFieldVar(eval=self.get_u)})
        return rte_dict

#----------------------- example --------------------

if __name__ == '__main__':
    from ibvpy.api import \
        TStepper as TS, FEGrid, RTraceGraph, RTraceDomainField, TLoop, \
        TLine, BCDof, DOTSEval, RTraceDomainListField

    fets_eval = FEQ4T()

    single_elem = True
    if single_elem:
        quad_elem = FEGrid(coord_min=(-5, -5),
                           coord_max=(10, 10),
                           fets_eval=fets_eval,
                           shape=(1, 1))
        print quad_elem.dof_grid_spec.node_coords
        print 'xxxx'
        ts = TS(
            sdomain=quad_elem,
            bcond_list=[BCDof(var='u', dof=i, value=0.) for i in [0, 3]] +
            [BCDof(var='f', dof=i, value=10) for i in [1, 2]],
            rtrace_list=[RTraceDomainListField(name='Temperature',
                                                    var='u', warp=False, idx=0,
                                                    record_on='update'),
                         RTraceDomainListField(name='Flux',
                                               var='eps', warp=False, idx=0,
                                               record_on='update'),
                         RTraceDomainListField(name='Jacobi determinant',
                                               var='J_det', warp=False, idx=0,
                                               record_on='update'),
                         #                                RTraceDomainListField(name = 'Shape functions' ,
                         #                                   var = 'N_mtx', warp = False, idx = 0,
                         # record_on = 'update'),
                         ]
        )

        # Add the time-loop control
        #
        tl = TLoop(tstepper=ts,
                   tline=TLine(min=0.0, step=0.5, max=1.0))

        tl.eval()

    else:
        # Tseval for a discretized line domain
        #
        tseval = DOTSEval(fets_eval=FEQ4T())

        # Discretization
        #
        line_domain = MGridDomain(lengths=(3., 1., 0.),
                                  shape=(8, 8, 0),
                                  n_e_nodes_geo=(1, 1, 0),
                                  n_e_nodes_dof=(1, 1, 0),
                                  node_map_geo=[0, 1, 3, 2],
                                  node_map_dof=[0, 1, 3, 2])

        # Put the tseval (time-stepper) into the spatial context of the
        # discretization and specify the response tracers to evaluate there.
        #
        right_dof = 2
        ts = TS(tse=tseval,
                sdomain=line_domain,
                bcond_list=[BCDof(var='u', dof=i, value=0.) for i in line_domain.get_left_dofs()[:, 0]] +
                [BCDof(var='u', dof=i, value=20)
                 for i in line_domain.get_right_dofs()[:, 0]],
                rtrace_list=[RTraceGraph(name='Fi,right over u_right (iteration)',
                                         var_y='F_int', idx_y=right_dof,
                                         var_x='U_k', idx_x=right_dof,
                                         record_on='update'),
                             #                      RTraceDomainField(name = 'Flux field' ,
                             #                      var = 'eps', idx = 0,
                             #                      record_on = 'update'),
                             RTraceDomainField(name='Temperature',
                                               var='u', idx=0,
                                               record_on='update'),
                             #                             RTraceDomainField(name = 'Shape functions' ,
                             #                                          var = 'N_mtx', idx = 0,
                             # record_on = 'update')

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

    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource=tl)
    app.main()
