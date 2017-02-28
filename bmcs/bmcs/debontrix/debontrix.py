'''
Created on Dec 20, 2016

This script demonstrates the looples implementation
of the finite element code for multilayer composite.

The current version demonstrating a uni-axial two-layer
continuum discretized with linear shape functions

Planned extensions:
   * Based on the element template generate the shape functions 
     based on the partition of unity concept. Generate symbolic 
     shape functions using sympy
   * Use sympy to construct the shape function derivatives up to
     the required order
   * Introduce the differential operator putting the differential
     terms together to represent the weak form
   * Introduce the integration template within an element template
   * Introduce the visualization template for an element.
   * Single layer works
   * Put the object into a view window.
   * Use a non-linear material model to demonstrate bond.
   * Use a non-linear material model to demonstrate strain-localization
   * Provide visualization of shear flow and slip.
   * Questions to be answered using this app

   * What is the relation to multi-dimensional case?
   * Generalization 
@author: rch
'''

from timeit import Timer

from ibvpy.api import \
    TStepperEval, IFETSEval, FETSEval, FEGrid
from ibvpy.mats.mats1D import MATS1DElastic
from mathkit.matrix_la.sys_mtx_assembly import SysMtxAssembly
from traits.api import \
    Int, implements, Array, \
    List, Property, cached_property, \
    Instance, Float

import numpy as np
import pylab as p
import sympy as sp


n_C = 2

ONE = np.ones((1,), dtype=np.float_)
DELTA_cd = np.identity(n_C)
c1 = np.arange(n_C) + 1
SWITCH_C = np.power(-1.0, c1)
SWITCH_CD = np.power(-1.0, c1[np.newaxis, :] + c1[:, np.newaxis])

r_ = sp.symbols('r')


class FETS1D2L(FETSEval):
    '''Example of a finite element definition.
    '''

    implements(IFETSEval)

    dim_slice = slice(0, 1)
    n_e_dofs = Int(2 * n_C)
    n_nodal_dofs = Int(n_C)
    dof_r = Array(value=[[-1], [1]])
    geo_r = Array(value=[[-1], [1]])
    vtk_r = Array(value=[[-1.], [1.]])
    vtk_cells = [[0, 1]]
    vtk_cell_types = 'Line'

    r_m = Array(value=[[-1], [1]], dtype=np.float_)
    w_m = Array(value=[1, 1], dtype=np.float_)

    Nr_i_geo = List([(1 - r_) / 2.0,
                     (1 + r_) / 2.0, ])

    dNr_i_geo = List([- 1.0 / 2.0,
                      1.0 / 2.0, ])

    Nr_i = Nr_i_geo
    dNr_i = dNr_i_geo

    N_mi_geo = Property()

    @cached_property
    def _get_N_mi_geo(self):
        return self.get_N_mi(sp.Matrix(self.Nr_i_geo, dtype=np.float_))

    dN_mid_geo = Property()

    @cached_property
    def _get_dN_mid_geo(self):
        return self.get_dN_mid(sp.Matrix(self.dNr_i_geo, dtype=np.float_))

    N_mi = Property()

    @cached_property
    def _get_N_mi(self):
        return self.get_N_mi(sp.Matrix(self.Nr_i, dtype=np.float_))

    dN_mid = Property()

    @cached_property
    def _get_dN_mid(self):
        return self.get_dN_mid(sp.Matrix(self.dNr_i, dtype=np.float_))

    def get_N_mi(self, Nr_i):
        return np.array([Nr_i.subs(r_, r)
                         for r in self.r_m], dtype=np.float_)

    def get_dN_mid(self, dNr_i):
        dN_mdi = np.array([[dNr_i.subs(r_, r)]
                           for r in self.r_m], dtype=np.float_)
        return np.einsum('mdi->mid', dN_mdi)


class DOTSEval(TStepperEval):

    fet = Instance(FETSEval)

    def _fet_default(self):
        return FETS1D2L()

    L_x = Float(1, input=True)

    eta = Float(0.1, input=True)

    n_E = Int(10, input=True)

    fe_grid = Property(Instance(FEGrid), depends_on='+input')

    @cached_property
    def _get_fe_grid(self):
        return FEGrid(coord_min=(0., ),
                      coord_max=(self.L_x, ),
                      shape=(self.n_E, ),
                      fets_eval=FETS1D2L())

    #=========================================================================
    # index maps
    #=========================================================================

    dof_ECid = Property(depends_on='+input')
    '''For a given element, layer, node number and dimension
    return the dof number
    '''
    @cached_property
    def _get_dof_ECid(self):
        dof_EiCd = self.fe_grid.dof_grid.cell_dof_map[..., np.newaxis]
        return np.einsum('EiCd->ECid', dof_EiCd)

    I_Ei = Property(depends_on='+input')
    '''For a given element and its node number return the global index
    of the node'''
    @cached_property
    def _get_I_Ei(self):
        return self.fe_grid.geo_grid.cell_grid.cell_node_map

    dof_E = Property(depends_on='+input')
    '''Get ordered array of degrees of freedom corresponding to each element.
    '''
    @cached_property
    def _get_dof_E(self):
        return self.dof_ECid.reshape(-1, self.fet.n_e_dofs)

    #=========================================================================
    # Coordinate arrays
    #=========================================================================

    X_Id = Property(depends_on='+input')
    'Coordinate of the node `I` in dimension `d`'
    @cached_property
    def _get_X_Id(self):
        return self.fe_grid.geo_grid.cell_grid.point_x_arr

    X_Eid = Property(depends_on='+input')
    'Coordinate of the node `i` in  element `E` in dimension `d`'
    @cached_property
    def _get_X_Eid(self):
        return self.X_Id[self.I_Ei, :]

    X_Emd = Property(depends_on='+input')
    'Coordinate of the integration point `m` of an element `E` in dimension `d`'
    @cached_property
    def _get_X_Emd(self):
        N_mi_geo = self.fet.N_mi_geo
        return np.einsum('mi,Eid->Emd', N_mi_geo, self.X_Eid)

    X_J = Property(depends_on='+input')
    '''Return ordered vector of nodal coordinates respecting the the order
    of the flattened array of elements, nodes and spatial dimensions.'''
    @cached_property
    def _get_X_J(self):
        return self.X_Eid.flatten()

    X_M = Property(depends_on='+input')
    '''Return ordered vector of global coordinates of integration points
    respecting the the order of the flattened array of elements, 
    nodes and spatial dimensions. Can be used for point-value visualization
    of response variables.'''
    @cached_property
    def _get_X_M(self):
        return self.X_Emd.flatten()

    #=========================================================================
    # cached time-independent terms
    #=========================================================================
    BB_ECidDjf = Property(Array)
    '''Product of shape function derivatives  mappings 
    in every integration point
    '''

    def _get_BB_ECidDjf(self):
        return self.constant_terms[0]

    NN_ECidDjf = Property(Array)
    '''Product of shape functions in every integration point
    '''

    def _get_NN_ECidDjf(self):
        return self.constant_terms[1]

    dN_Eimd = Property
    '''Shape function derivatives in every integration point
    '''

    def _get_dN_Eimd(self):
        return self.constant_terms[2]

    sN_Cim = Property
    '''Slip operator between the layers C = 0,1
    '''

    def _get_sN_Cim(self):
        return self.constant_terms[3]

    constant_terms = Property(depends_on='+input')
    '''Procedure calculating all constant terms of the finite element
    algorithm including the geometry mapping (Jacobi), shape functions and the kinematics needed
    for the integration of stresses and stifnesses in every material point.
    '''
    @cached_property
    def _get_constant_terms(self):
        fet = self.fet
        dN_mid_geo = fet.dN_mid_geo
        N_mi = fet.N_mi
        dN_mid = fet.dN_mid
        w_m = fet.w_m

        A_C = np.ones((n_C,), dtype=np.float_)
        A_C[0] *= self.eta
        # Geometry approximation / Jacobi transformation
        J_Emde = np.einsum('mid,Eie->Emde', dN_mid_geo, self.X_Eid)
        J_det_Em = np.linalg.det(J_Emde)
        J_inv_Emed = np.linalg.inv(J_Emde)
        # Quadratic forms
        dN_Eimd = np.einsum('mid,Eied->Eime', dN_mid, J_inv_Emed)
        sN_Cim = np.einsum('C,mi->Cim', SWITCH_C, N_mi)
        BB_ECidDjf = np.einsum('m, CD, C, Eimd,Ejmf,Em->ECidDjf',
                               w_m, DELTA_cd, A_C, dN_Eimd, dN_Eimd, J_det_Em)
        NN_ECidDjf = np.einsum('m, d,f,CD,mi,mj,Em->ECidDjf',
                               w_m, ONE, ONE, SWITCH_CD, N_mi, N_mi, J_det_Em)
        return BB_ECidDjf, NN_ECidDjf, dN_Eimd, sN_Cim

    G = Float(1.0, bc_changed=True)

    w = Float(0.01, bc_changed=True)

    d = Property(depends_on='+input,+bc_changed')

    @cached_property
    def _get_d(self):

        fet = self.fet
        n_dof_tot = self.fe_grid.n_dofs

        K_ECidDjf = self.BB_ECidDjf + self.NN_ECidDjf * self.G
        # Stapled matrices
        K_Eij = K_ECidDjf.reshape(-1, fet.n_e_dofs, fet.n_e_dofs)
        # System matrix
        K = SysMtxAssembly()
        K.add_mtx_array(K_Eij, self.dof_E)
        K.register_constraint(0, 0.0)
        K.register_constraint(n_dof_tot - 1, self.w)
        F_ext = np.zeros((n_dof_tot,), np.float_)
        K.apply_constraints(F_ext)
        d = K.solve(F_ext)
        return d

    d_C = Property

    def _get_d_C(self):
        d_ECid = self.d[self.dof_ECid]
        return np.einsum('ECid->EidC', d_ECid).reshape(-1, n_C)

    eps_C = Property

    def _get_eps_C(self):
        d_ECid = self.d[self.dof_ECid]
        eps_EmdC = np.einsum('Eimd,ECid->EmdC', self.dN_Eimd, d_ECid)
        return eps_EmdC.reshape(-1, n_C)

    u_C = Property
    '''Displacement field
    '''

    def _get_u_C(self):
        d_ECid = self.d[self.dof_ECid]
        N_mi = self.fet.N_mi
        u_EmdC = np.einsum('mi,ECid->EmdC', N_mi, d_ECid)
        return u_EmdC.reshape(-1, n_C)

    s = Property

    def _get_s(self):
        d_ECid = self.d[self.dof_ECid]
        s_Emd = np.einsum('Cim,ECid->Emd', self.sN_Cim, d_ECid)
        return s_Emd.flatten()

    f_IC = Property

    def _get_f_IC(self):
        K_ECidDjf = self.BB_ECidDjf + self.NN_ECidDjf * self.G
        d_ECid = self.d[self.dof_ECid]
        f_ECid = np.einsum('ECidDjf,EDjf->ECid', K_ECidDjf, d_ECid)
        f_Ei = f_ECid.reshape(-1, self.fet.n_e_dofs)
        f_I = np.bincount(self.dof_E.flatten(), weights=f_Ei.flatten())
        return f_I.reshape(-1, n_C)

    def plot_f_C(self, ax):
        ax.plot(self.X_Id.flatten(), self.f_IC)

    def plot_u_C(self, ax):
        ax.plot(self.X_J, self.u_C)

    def plot_eps_C(self, ax):
        ax.plot(self.X_M, self.eps_C)

    def plot_s(self, ax):
        ax.plot(self.X_J, self.s)

    def plot(self, fig):
        ax = fig.add_subplot(221)
        self.plot_f_C(ax)
        ax = fig.add_subplot(222)
        self.plot_eps_C(ax)
        ax = fig.add_subplot(223)
        self.plot_s(ax)
        ax = fig.add_subplot(224)
        self.plot_u_C(ax)


def plot_mlab(X, I_Ei, fet):
    subcell_offsets, subcell_lengths, subcells, subcell_types = \
        fet.vtk_node_cell_data
    print subcell_offsets
    print subcell_lengths
    print subcells
    print subcell_types
    print I_Ei

if __name__ == '__main__':
    dots = DOTSEval(n_E=5,
                    L_x=1.0,
                    G=1.0)

    print 's_C', dots.X_J, dots.u_C
    fig = p.figure()
    dots.plot(fig)
    p.show()
