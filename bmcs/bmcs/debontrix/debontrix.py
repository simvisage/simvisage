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

from ibvpy.fets.fets_eval import IFETSEval, FETSEval
from ibvpy.mats.mats1D import MATS1DElastic
from ibvpy.mesh.fe_grid import FEGrid
from mathkit.matrix_la.sys_mtx_assembly import SysMtxAssembly
from traits.api import \
    Int, implements, Array, \
    List, Property, cached_property

import numpy as np
import pylab as p
import sympy as sp


n_C = 2
n_E = 10

ONE = np.ones((1,), dtype=np.float_)
DELTA_cd = np.identity(n_C)
c1 = np.arange(n_C) + 1
SWITCH_cd = np.power(-1.0, c1[np.newaxis, :] + c1[:, np.newaxis])

r_ = sp.symbols('r')

L_x, L_y, L_z = 1, 1, 1


class FETS1D2L(FETSEval):
    '''Two layer bar element.
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
#     r_m = Array(value=[[0.0]], dtype=np.float_)
#     w_m = Array(value=[2.0], dtype=np.float_)

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
                 shape=(n_E, ),
                 fets_eval=fet)
n_dof_tot = fe_grid.n_dofs

N_mi_geo = fet.N_mi_geo
dN_mid_geo = fet.dN_mid_geo
N_mi = fet.N_mi
dN_mid = fet.dN_mid
w_m = fet.w_m

dof_EiCd = fe_grid.dof_grid.cell_dof_map[..., np.newaxis]
dof_ECid = np.einsum('EiCd->ECid', dof_EiCd)
I_Ei = fe_grid.geo_grid.cell_grid.cell_node_map
X_Id = fe_grid.geo_grid.cell_grid.point_x_arr
X_Eid = X_Id[I_Ei, :]
X_Emd = np.einsum('mi,Eid->Emd', N_mi_geo, X_Eid)

A_C = np.ones((n_C,), dtype=np.float_)
A_C[0] *= 10.0
G = 1.0
# Geometry approximation / Jacobi transformation
J_Emde = np.einsum('mid,Eie->Emde', dN_mid_geo, X_Eid)
J_det_Em = np.linalg.det(J_Emde)
J_inv_Emed = np.linalg.inv(J_Emde)
# Quadratic forms
dN_Eimd = np.einsum('mid,Eied->Eime', dN_mid, J_inv_Emed)
BB_ECidDjf = np.einsum('m, CD, C, Eimd,Ejmf,Em->ECidDjf',
                       w_m, DELTA_cd, A_C, dN_Eimd, dN_Eimd, J_det_Em)
NN_ECidDjf = np.einsum('m, d,f,CD,mi,mj,Em->ECidDjf',
                       w_m, ONE, ONE, SWITCH_cd, N_mi, N_mi, J_det_Em) * G
K_ECidDjf = BB_ECidDjf + NN_ECidDjf
# Stapled matrices
K_Eij = K_ECidDjf.reshape(-1, fet.n_e_dofs, fet.n_e_dofs)
# System matrix
dof_E = dof_ECid.reshape(-1, fet.n_e_dofs)
K = SysMtxAssembly()
K.add_mtx_array(K_Eij, dof_E)
K.register_constraint(0, 0.0)
K.register_constraint(n_dof_tot - 1, 0.01)
F_ext = np.zeros((n_dof_tot,), np.float_)
K.apply_constraints(F_ext)
d = K.solve(F_ext)
# Postprocessing
X_J = X_Eid.flatten()
X_M = X_Emd.flatten()
d_ECid = d[dof_ECid]
d_C = np.einsum('ECid->EidC', d_ECid).reshape(-1, n_C)
eps_EmdC = np.einsum('Eimd,ECid->EmdC', dN_Eimd, d_ECid)
eps_C = eps_EmdC.reshape(-1, n_C)
sig_EmdC = eps_EmdC
f_ECid = np.einsum('ECidDjf,EDjf->ECid', K_ECidDjf, d_ECid)
f_Ei = f_ECid.reshape(-1, fet.n_e_dofs)
f_I = np.bincount(dof_E.flatten(), weights=f_Ei.flatten())
f_IC = f_I.reshape(-1, n_C)
P = f_IC[-1, -1]


def plot_pylab(X, field_C):
    p.plot(X, field_C)


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

p.subplot(211)
plot_pylab(X_Id.flatten(), f_IC)
p.subplot(212)
plot_pylab(X_J, eps_C)
#plot_pylab(X_Id.flatten(), f_IC)
p.show()
