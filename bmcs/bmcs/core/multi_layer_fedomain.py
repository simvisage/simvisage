'''
Created on Dec 20, 2016

This script demonstrates the looples implementation
of the finite element code for multilayer composite.

The current version demonstrating a uni-axial two-layer
continuum discritized with linear shape functions

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

@author: rch
'''

from traits.api import \
    Int, implements, Array, \
    List, Property, cached_property
from ibvpy.fets.fets_eval import IFETSEval, FETSEval
from ibvpy.mats.mats1D import MATS1DElastic
from ibvpy.mesh.fe_grid import FEGrid
from mathkit.matrix_la.sys_mtx_assembly import SysMtxAssembly
import numpy as np
import pylab as p
import sympy as sp

n_c = 2
n_e = 1000

DELTA_cd = np.identity(n_c)
c1 = np.arange(n_c) + 1
SWITCH_cd = np.power(-1.0, c1[np.newaxis, :] + c1[:, np.newaxis])

r_ = sp.symbols('r')

L_x, L_y, L_z = n_e, 1, 1


class FETS1D2L(FETSEval):
    '''Two layer bar element.
    '''

    implements(IFETSEval)

    dim_slice = slice(0, 1)
    n_e_dofs = Int(2 * n_c)
    n_nodal_dofs = Int(n_c)
    dof_r = Array(value=[[-1], [1]])
    geo_r = Array(value=[[-1], [1]])
    vtk_r = Array(value=[[-1.], [1.]])
    vtk_cells = [[0, 1]]
    vtk_cell_types = 'Line'

    r_m = Array(value=[[-1], [1]], dtype=np.float_)
    w_m = Array(value=[1, 1], dtype=np.float_)

    # Integration parameters
    ngp_r = 2

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


fet = FETS1D2L()

fe_grid = FEGrid(coord_min=(0., ),
                 coord_max=(L_x, ),
                 shape=(n_e, ),
                 fets_eval=fet)
n_dof_tot = fe_grid.n_dofs

N_mi_geo = fet.N_mi_geo
dN_mid_geo = fet.dN_mid_geo
N_mi = fet.N_mi
dN_mid = fet.dN_mid
w_m = fet.w_m

dof_Eic = fe_grid.dof_grid.cell_dof_map
I_Ei = fe_grid.geo_grid.cell_grid.cell_node_map
X_Id = fe_grid.geo_grid.cell_grid.point_x_arr
X_Eid = X_Id[I_Ei]

A_c = np.ones((n_c,), dtype=np.float_)
# Geometry approximation / Jacobi transformation
J_Emde = np.einsum('mid,Eie->Emde', dN_mid_geo, X_Eid)
J_det_Em = np.linalg.det(J_Emde)
J_inv_Emed = np.linalg.inv(J_Emde)
# Quadratic forms
dN_Eimd = np.einsum('mid,Eied->Eime', dN_mid, J_inv_Emed)
BB_ECidDjf = np.einsum('m, CD, C, Eimd,Ejmf,Em->EiCdjDf',
                       w_m, DELTA_cd, A_c, dN_Eimd, dN_Eimd, J_det_Em)
NN_ECidDjf = np.einsum('m, CD,mi,mj,Em->EiCjD',
                       w_m, SWITCH_cd, N_mi, N_mi, J_det_Em) * 0.1
# Stapled matrices
BB_Eij = BB_ECidDjf.reshape(-1, fet.n_e_dofs, fet.n_e_dofs)
NN_Eij = NN_ECidDjf.reshape(-1, fet.n_e_dofs, fet.n_e_dofs)
# System matrix
dof_E = dof_Eic.reshape(-1, fet.n_e_dofs)
K = SysMtxAssembly()
K.add_mtx_array(BB_Eij + NN_Eij, dof_E)
K.register_constraint(0, 0.0)
K.register_constraint(n_dof_tot - 1, 0.01)
F_ext = np.zeros((n_dof_tot,), np.float_)
K.apply_constraints(F_ext)
d = K.solve(F_ext)
# Postprocessing
d_Eic = d[dof_Eic]
d_c = d_Eic.reshape(-1, n_c)
X = X_Eid.flatten()


def plot_pylab(X, d_c):
    p.plot(X, d_c)
    p.show()


def plot_mlab(X, I_Ei):
    subcell_offsets, subcell_lengths, subcells, subcell_types = \
        fet.vtk_node_cell_data
    print subcell_offsets
    print subcell_lengths
    print subcells
    print subcell_types
    print fe_grid.geo_grid.cell_grid.point_X_arr.shape
    print fe_grid.geo_grid.cell_grid.point_x_arr.shape
    print I_Ei

plot_pylab(X, d_c)
