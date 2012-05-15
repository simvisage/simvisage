

#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Dec 9, 2009 by: rch

'''
MUSHROOF - Demostrator SFB532

@todos 
1) geo_transform for the bottom bar [Done]
2) quadratic  boundary condition not working

'''
from enthought.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool, HasTraits

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel

#from ibvpy.rtrace.rt_domain_list_field import \
#    RTraceDomainListField

#from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
#    MATS2DElastic
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic
#from ibvpy.mats.mats3D.mats3D_sdamage.mats3D_sdamage import \
#    MATS3DScalarDamage
#from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import \
#    MATS3DMicroplaneDamage
#from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
#    MATS2D5MicroplaneDamage, PhiFnGeneral, PhiFnStrainHardening

from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U
#from ibvpy.fets.fets2D5.fets2D58h import \
#    FETS2D58H
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U

from ibvpy.mesh.fe_grid import \
    FEGrid

from mathkit.mfn import MFnLineArray
from numpy import \
    array, tensordot, dot, zeros, c_, ix_, shape, \
    cos, sin, arctan, where, abs, all, any, diag
from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos, pi as Pi
from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

from geo_column import GEOColumn

# Interpolation
from scipy.interpolate import Rbf

from simiter.sim_pstudy import ISimModel, SimOut, SimPStudy

from hp_shell import HPShell

from time import time

def temperature_strain(X_pnt, x_pnt):

    alpha = 1.3e-5
    t_up = -100
    t_lo = -100
    delta_t = t_lo + (t_up - t_lo) * x_pnt[2]
    epsilon_0 = alpha * delta_t
    # return the initial volumetric strain tensor with n_dims  
    return diag([ epsilon_0 for i in range(3) ])

#class Column( HasTraits ):
#
#    X0 = Array( float, value = [ 4., 4., -3. ] )
#
#    width_top = Float( 0.45, unit = 'm' )
#    width_bottom = Float( 0.3, unit = 'm' )
#    height = Float( 3.0, unit = 'm' )
#    r_pipe = Float( 0.1, unit = 'm' )
#
#    def __call__( self, points ):
#
#        xi, yi, zi = points[:, 0], points[:, 1], points[:, 2]
#
#        # grid from 0 to 1, shift origin for rotation
#        # 
#        xi -= 0.5
#        yi -= 0.5
#
#        # setting of global coordinates different width over height 
#        #
#        xi = ( ( self.width_top - self.width_bottom ) * zi + self.width_bottom ) * xi
#        yi = ( ( self.width_top - self.width_bottom ) * zi + self.width_bottom ) * yi
#        zi *= self.height
#
#        # rotation of 45 with global coordinates
#        #
#        x = cos( Pi / 4.0 ) * xi + sin( Pi / 4.0 ) * yi
#        y = -sin( Pi / 4.0 ) * xi + cos( Pi / 4.0 ) * yi
#        z = zi + self.X0[2]
#
#        #shift of internal elements
#        #
#        r_0 = where( ( xi ** 2 + yi ** 2 ) ** 0.5 >= self.r_pipe, self.r_pipe , ( xi ** 2 + yi ** 2 ) ** 0.5 )
##        r_0 = where( ( abs( xi ) == 0.5 or abs( yi ) == 0.5 ).all() , self.r_pipe , ( xi ** 2 + yi ** 2 ) ** 0.5 )
#
#        scale_r = where ( r_0 == 0, 0, self.r_pipe / r_0 )
#
#        x = scale_r * x + self.X0[0]
#        y = scale_r * y + self.X0[1]
#        return c_[x, y, z]

class SFBMushRoofModel(IBVModel):
    '''SFB - Demonstrator model specification.
    '''
    implements(ISimModel)

    # dimensions of one quarter of the shell structure [m]
    #
    length_xy = Float(4.)
    length_z = Float(1.062)
    t_shell = Float(0.06)

    # choose model and discretization:
    #
    mushroof_part = 'one'
    n_elems_xy_quarter = Int(4, ps_levels = (3, 5, 2))
    n_elems_z = Int(1, input = True)

    # grid parameters:
    #
    shift_elems_column = Bool(False, input = True)
    const_edge_elem = Bool (False, input = True)
    width_column = Float(0.45, input = True) # [m]
    n_elems_col = Int(1, input = True)
    delta_h_scalefactor = Float(1.00)

    # Material properties: Youngs-modulus, poision ratio 
    #
    E_roof = Float(28700) # [MPa]
    E_column = Float(37000) # [MPa] Annahme C60
    nu = Float(0.2) # [-]

    n_elems_xy = Property(Int , depends_on = '+ps_levels, +input')
    def _get_n_elems_xy(self):
        return int (self.n_elems_xy_quarter * 2)


    # variable type of the finite element
    fets = Instance(FETSEval , ps_levels = ['fe_linear',
#                                            'fe2d5_quad_serendipity',
                                            'fe_quad_serendipity' ])
    def _fets_default(self):
#        return self.fe_linear
        return self.fe_quad_serendipity

    mats_roof = Instance(MATS3DElastic)
    def _mats_roof_default(self):
        return MATS3DElastic(E = self.E_roof, nu = self.nu)#, initial_strain = temperature_strain )

    mats_column = Instance(MATS3DElastic)
    def _mats_column_default(self):
        return MATS3DElastic(E = self.E_column, nu = self.nu)#, initial_strain = temperature_strain )

    fe_linear = Instance(FETSEval, transient = True)
    def _fe_linear_default(self):
        return FETS3D8H(mats_eval = self.mats_roof)

    fe_quad_serendipity = Instance(FETSEval, transient = True)
    def _fe_quad_serendipity_default(self):
        return FETS3D8H20U(mats_eval = self.mats_roof)

    fe2d5_quad_serendipity = Instance(FETSEval, transient = True)
    def _fe2d5_quad_serendipity_default(self):
        return FETS2D58H20U(mats_eval = self.mats_roof)

    fe_quad_lagrange = Instance(FETSEval, transient = True)
    def _fe_quad_lagrange_default(self):
        return FETS3D8H27U(mats_eval = self.mats_roof)

    fe_column = Instance(FETSEval, transient = True)
    def _fe_column_default(self):
        fets = FETS3D8H20U(mats_eval = self.mats_column)
        fets.vtk_r *= 0.95
        return fets

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name = 'u_z_free_corner', unit = 'm')]
#                 SimOut( name = 'maximum principle stress', unit = 'MPa' ), ]

    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()

        U_z = U[self.corner_dof_r00][0, 0, 2]

#        max_princ_stress = max( self.max_princ_stress._get_field_data().flatten() )

        return array([ U_z ],
#                       max_princ_stress ],
                        dtype = 'float_')

    tline = Instance(TLine)
    def _tline_default(self):
        return TLine(min = 0.0, step = 1.0, max = 1.0)

    max_princ_stress = Instance(RTraceDomainListField)
    def _max_princ_stress_default(self):
        return RTraceDomainListField(name = 'max principle stress' , idx = 0,
                                      var = 'max_principle_sig', warp = True,
                                      record_on = 'update',)

    rtrace_list = List
    def _rtrace_list_default(self):
        return [  self.max_princ_stress, self.sig_app, self.u  ]

    sig_app = Property(Instance(RTraceDomainListField), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_sig_app(self):
        return RTraceDomainListField(name = 'sig_app' ,
#                                      position = 'int_pnts',
                                      var = 'sig_app',
                                      record_on = 'update',)

    u = Property(Instance(RTraceDomainListField), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_u(self):
        return RTraceDomainListField(name = 'displacement' ,
                                      var = 'u', warp = True,
                                      record_on = 'update',)



    #----------------------------------------------------------------------------------
    # HP-Shell specification
    #----------------------------------------------------------------------------------

    hp_shells = Property(List(HPShell) , depends_on = '+ps_levels, +input')
    @cached_property
    def _get_hp_shells(self):
        X_list = [ [0, 0, 0], [8, 0, 0], [0, 8, 0], [8, 8, 0] ]
        return [
                 HPShell(length_x = self.length_xy,
                    length_y = self.length_xy,
                    length_z = self.length_z,
                    t_shell = self.t_shell,
                    X0 = X,
                    n_elems_xy = self.n_elems_xy,
                    n_elems_z = self.n_elems_z,
                    shift_elems_column = self.shift_elems_column,
                    delta_h_scalefactor = self.delta_h_scalefactor,
                    const_edge_elem = self.const_edge_elem,
                    width_column = self.width_column,
                    n_elems_col = self.n_elems_col,
                    mushroof_part = self.mushroof_part)
                 for X in X_list ]

    fe_grid_roofs = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_grid_roofs(self):
        return [
                FEGrid(coord_min = (0.0, 0.0, 0.0),
                       coord_max = (1.0, 1.0, 1.0),
                       geo_transform = hp_shell,
                       shape = (self.n_elems_xy, self.n_elems_xy, self.n_elems_z),
                       fets_eval = self.fets)
                for hp_shell in self.hp_shells ]


    #----------------------------------------------------------------------------------
    # Column specification
    #----------------------------------------------------------------------------------

    delta_h_scale_factor = Float(1.0, input = True)
    column_width_top = Float(0.45, unit = 'm')
    column_width_bottom = Float(0.3, unit = 'm')
    column_height = Float(3.0, unit = 'm')

    columns = Property(Instance(HPShell) , depends_on = '+ps_levels')
    @cached_property
    def _get_columns(self):
        ch = self.column_height / self.delta_h_scale_factor
        X_list = [[ 4., 4., -ch ],
                  [ 12., 4., -ch ],
                  [ 4., 12., -ch ],
                  [ 12., 12., -ch ]]
        return [ GEOColumn(width_top = self.column_width_top,
                            X0 = X,
                            width_bottom = self.column_width_bottom,
                            height = ch)
                            for X in X_list ]

    fe_grid_columns = Property(Instance(FEGrid), depends_on = '+ps_levels')
    @cached_property
    def _get_fe_grid_columns(self):
        fe_grid_columns = [ FEGrid(coord_min = (0.0, 0.0, 0.0),
                             coord_max = (1.0, 1.0, 1.0),
                             geo_transform = column,
                             shape = (3, 3, 3),
                             fets_eval = self.fe_column)
                             for column in self.columns ]
#        for i, fe_grid_column in enumerate( fe_grid_columns ):
#            interior_elems = fe_grid_column[ 1:-1, 1:-1, :, :, :, : ].elems
#            fe_grid_column.inactive_elems = list( interior_elems )
#            fe_grid_columns[i] = fe_grid_column
        return fe_grid_columns


    #----------------------------------------------------
    # loading cases (LC):
    #----------------------------------------------------

    #--- LC0: verifying load for ps_study
    # g = 1.0 kN/m^2 
    # orientation: global z-direction; 
    f_z = Float(-1e-3, unit = 'MN/m^3')

    #--- LC1: dead load
    # g = 22.4 kN/m^3 
    # orientation: global z-direction; 
    material_density_roof = Float(-0.0224, unit = 'MN/m^3')
    material_density_column = Float(-0.0224, unit = 'MN/m^3')

    #--- LC2 additional dead load 
    # gA = 0,20 kN/m^2 
    # orientation: global z-direction (following the curved structure); 
    surface_load_gA = Float(-0.20e-3, unit = 'MN/m^2')

    #--- LC3 snow
    # s = 0,79 kN/m^2 
    # orientation: global z-direction (projection); 
    surface_load_s = Float(-0.84e-3, unit = 'MN/m^2')

    #--- LC4 wind (pressure) 
    # w = 0,13 kN/m^2 
    # orientation: local t-direction (surface normal); 
    surface_load_w = Float(-0.13e-3, unit = 'MN/m^2')

    force_bc_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_force_bc_list(self):
        # NOTE: additional line-loads at the edge of the roof need to be considered!  

        force_bc_list = []
        for roof, column in zip(self.fe_grid_roofs, self.fe_grid_columns):
            upper_surface = roof[:, :, -1, :, :, -1]
#            whole_roof = roof[:, :, :, :, :, :]
#            whole_column = column[:, :, :, :, :, :]

            force_bc = [
#                         # LC0: dead load
                         BCSlice(var = 'f', value = self.f_z, dims = [2],
                                  integ_domain = 'global',
                                  slice = upper_surface),

#                         # LC1: dead load
#                         BCSlice( var = 'f', value = self.material_density_roof, dims = [2],
#                                  integ_domain = 'global',
#                                  slice = whole_roof ),
#                         BCSlice( var = 'f', value = self.material_density_column, dims = [2],
#                                  integ_domain = 'global',
#                                  slice = whole_column ),
                         # LC2: additional dead load
#                         BCSlice( var = 'f', value = self.surface_load_gA, dims = [2],
#                                  integ_domain = 'global',
#                                  slice = upper_surface ),
#                         # LC3: snow load         
#                         BCSlice( var = 'f', value = self.surface_load_s, dims = [2],
#                                  integ_domain = 'global',
#                                  slice = upper_surface )
                       ]
            force_bc_list += force_bc

        return force_bc_list

    displ_bc_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_displ_bc_list(self):

        displ_bc_list = []
        for roof, column in zip(self.fe_grid_roofs, self.fe_grid_columns):
            bc_list = [
                       BCSlice(var = 'u'  , dims = [0, 1, 2],
                                slice = roof[self.n_elems_xy_quarter - 1,
                                               self.n_elems_xy_quarter - 1,
                                               0,
                                               0, -1, 0 ],
                                link_slice = column[ 0, 0, -1, 0, 0, -1],
                                link_coeffs = [1.0],
                                value = 0.),
                        BCSlice(var = 'u'  , dims = [0, 1, 2],
                                slice = roof[self.n_elems_xy_quarter - 1,
                                               self.n_elems_xy_quarter - 1 ,
                                               0,
                                               - 1, 0, 0],
                                link_slice = column[ -1, 0, -1, -1, 0, -1],
                                link_coeffs = [1.0],
                                value = 0.),
                        BCSlice(var = 'u'  , dims = [0, 1, 2],
                                slice = roof[self.n_elems_xy_quarter,
                                               self.n_elems_xy_quarter,
                                               0,
                                               - 1, 0, 0 ],
                                link_slice = column[ -1, -1, -1, -1, -1, -1],
                                link_coeffs = [1.0],
                                value = 0.),
                        BCSlice(var = 'u'  , dims = [0, 1, 2],
                                    slice = roof[self.n_elems_xy_quarter,
                                                   self.n_elems_xy_quarter,
                                                   0,
                                                   0, -1, 0 ],
                                    link_slice = column[ 0, -1, -1, 0, -1, -1],
                                    link_coeffs = [1.0],
                                    value = 0.),
                        BCSlice(var = 'u', dims = [0, 1, 2],
                                 slice = column[ :, :, 0, :, :, 0 ],
#                                 slice = column[ 0, 0, 0, -1, -1, 0 ],
                                 value = 0.0)
                        ]
            displ_bc_list += bc_list
        return displ_bc_list


    wind_bc_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_wind_bc_list(self):

        w_left_edge = 0.91e-3 # [MN/m]
        w_right_edge = 0.39e-3 # [MN/m]
        w_bottom_edge_0 = 1.14e-3 # [Mn/m]
        w_bottom_edge_1 = 0.65e-3 # [Mn/m]
        w_top_edge_0 = -w_bottom_edge_0 # [Mn/m]
        w_top_edge_1 = -w_bottom_edge_1 # [Mn/m]
        w_face_0_left = 0.89e-3 # MN/m^2
        w_face_0_right = -0.13e-3 # MN/m^2
        w_face_1_left = 0.72e-3 # MN/m^2
        w_face_1_right = w_face_0_right

        r00, r10, r01, r11 = self.fe_grid_roofs

        left_edge_0 = r00[ 0, :, -1, 0, :, -1 ]
        left_edge_1 = r01[ 0, :, -1, 0, :, -1 ]
        right_edge_0 = r10[ -1, :, -1, 0, :, -1 ]
        right_edge_1 = r11[ -1, :, -1, 0, :, -1 ]
        bottom_edge_0 = r00[ :, 0, -1, :, 0, -1 ]
        bottom_edge_1 = r10[ :, 0, -1, :, 0, -1 ]
        top_edge_0 = r01[ :, -1, -1, :, -1, -1 ]
        top_edge_1 = r11[ :, -1, -1, :, -1, -1 ]

        n_e_q = self.n_elems_xy_quarter
        face_00_left = r00[ :n_e_q, :, -1, :, :, -1 ]
        face_00_right = r00[ n_e_q:, :, -1, :, :, -1 ]
        face_01_left = r01[ :n_e_q, :, -1, :, :, -1 ]
        face_01_right = r01[ n_e_q:, :, -1, :, :, -1 ]
        face_10_left = r10[ :n_e_q, :, -1, :, :, -1 ]
        face_10_right = r10[ n_e_q:, :, -1, :, :, -1 ]
        face_11_left = r11[ :n_e_q, :, -1, :, :, -1 ]
        face_11_right = r11[ n_e_q:, :, -1, :, :, -1 ]

        wind_list = [
                   # left edge - x direction
                   BCSlice(var = 'f'  , dims = [0],
                            slice = left_edge_0,
                            value = w_left_edge),
                   BCSlice(var = 'f'  , dims = [0],
                            slice = left_edge_1,
                            value = w_left_edge),
                   # right edge - x direction
                   BCSlice(var = 'f'  , dims = [0],
                            slice = right_edge_0,
                            value = w_right_edge),
                   BCSlice(var = 'f'  , dims = [0],
                            slice = right_edge_1,
                            value = w_right_edge),
                   # bottom edge - y direction
                   BCSlice(var = 'f'  , dims = [1],
                            slice = bottom_edge_0,
                            value = w_bottom_edge_0),
                   BCSlice(var = 'f'  , dims = [1],
                            slice = bottom_edge_1,
                            value = w_bottom_edge_1),
                   # top edge - y direction
                   BCSlice(var = 'f'  , dims = [1],
                            slice = top_edge_0,
                            value = w_top_edge_0),
                   BCSlice(var = 'f'  , dims = [1],
                            slice = top_edge_1,
                            value = w_top_edge_1),
                   # upper face left - left 
                   BCSlice(var = 'f'  , dims = [2],
                            slice = face_00_left,
                            value = w_face_0_left),
                   BCSlice(var = 'f'  , dims = [2],
                            slice = face_01_left,
                            value = w_face_0_left),
                   # upper face left - right 
                   BCSlice(var = 'f'  , dims = [2],
                            slice = face_00_right,
                            value = w_face_0_right),
                   BCSlice(var = 'f'  , dims = [2],
                            slice = face_01_right,
                            value = w_face_0_right),
                   # upper face right - left 
                   BCSlice(var = 'f'  , dims = [2],
                            slice = face_10_left,
                            value = w_face_1_left),
                   BCSlice(var = 'f'  , dims = [2],
                            slice = face_11_left,
                            value = w_face_1_left),
                   # upper face right - right 
                   BCSlice(var = 'f'  , dims = [2],
                            slice = face_10_right,
                            value = w_face_1_right),
                   BCSlice(var = 'f'  , dims = [2],
                            slice = face_11_right,
                            value = w_face_1_right),
                        ]
        return wind_list

    link_bc_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_link_bc_list(self):
        r00, r10, r01, r11 = self.fe_grid_roofs
        link_bc_list = [
                       BCSlice(var = 'u'  , dims = [0, 1, 2],
                                slice = r00[-1, :, -1, -1, :-1, -1 ],
                                link_slice = r10[ 0, :, -1, 0, :-1, -1],
                                link_coeffs = [1.0],
                                value = 0.),
                        BCSlice(var = 'u'  , dims = [0, 1, 2],
                                slice = r01[:, 0, -1, :-1, 0, -1],
                                link_slice = r00[ :, -1, -1, :-1, -1, -1],
                                link_coeffs = [1.0],
                                value = 0.),
                        BCSlice(var = 'u'  , dims = [0, 1, 2],
                                slice = r11[ 0, :, -1, 0, 1:, -1 ],
                                link_slice = r01[ -1, :, -1, -1, 1:, -1],
                                link_coeffs = [1.0],
                                value = 0.),
                        BCSlice(var = 'u'  , dims = [0, 1, 2],
                                slice = r10[ :, -1, -1, 1:, -1, -1 ],
                                link_slice = r11[ :, 0, -1, 1:, 0, -1],
                                link_coeffs = [1.0],
                                value = 0.),
                        ]
        return link_bc_list


    # time loop
    tloop = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_tloop(self):
        self.fets.vtk_r *= 0.95

        roofs = self.fe_grid_roofs

        columns = self.fe_grid_columns

        r00, r10, r01, r11 = roofs

        self.corner_dof_r00 = r00[0, 0, 0, 0, 0, 0].dofs

        rtrace_list = self.rtrace_list

        ts = TS(sdomain = roofs + columns,
                 dof_resultants = True,
#                 bcond_list = self.wind_bc_list + self.displ_bc_list + self.link_bc_list,
                 bcond_list = self.force_bc_list + self.displ_bc_list + self.link_bc_list,
#                 bcond_list = self.displ_bc_list + self.link_bc_list,
                 rtrace_list = rtrace_list
               )

        # Add the time-loop control
        tloop = TLoop(tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline)

        return tloop

#        inactivate_elems = Property(depends_on='+ps_levels, +input' )
#        def _get_inactivate_elems( self ):
#            for i in range( 0, len( self.fe_grid_columns ) ):
#                interior_elems = self.fe_grid_columns[i][ 1:-1, 1:-1, :, :, :, : ].elems
#                self.fe_grid_columns[i].inactivate_elems = list( interior_elems )


if __name__ == '__main__':

    sim_model = SFBMushRoofModel(n_elems_xy_quarter = 4,
                                  n_elems_z = 1,
                                  const_edge_elem = False,
                                  shift_elems_column = True,
                                  width_column = 0.45,
                                  n_elems_col = 1,
                                  inactivate_elems_column = True)


    do = 'ui'

    if do == 'eval':
        print 'eval', sim_model.peval()

    if do == 'ui':

        sim_model.peval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource = sim_model)
        app.main()

    elif do == 'ps':

        sim_ps = SimPStudy(sim_model = sim_model)
        sim_ps.configure_traits()

    elif do == 'pickle':

        import pickle
        filename = '/tmp/sim.pickle'
        file = open(filename, 'w')
        pickle.dump(sim_model, file)
        file.close()
        file = open(filename, 'r')
        sm = pickle.load(file)
        file.close()
