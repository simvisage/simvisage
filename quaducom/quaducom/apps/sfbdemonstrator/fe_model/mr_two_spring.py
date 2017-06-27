'''
MUSHROOF - For double symmetric loading cases, one Roof is sufficient 

TODO: @ Andreas
     - split dead load cases, include additional dead load case
'''

from etsproxy.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool, HasTraits, Enum, Dict, Str

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel


from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.mats.mats1D import \
    MATS1DElastic

from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets1D.fets1D2l import \
    FETS1D2L
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U
from ibvpy.mesh.fe_grid import \
    FEGrid
from ibvpy.mesh.fe_refinement_grid import \
    FERefinementGrid
from ibvpy.mesh.fe_subdomain import \
    FESubDomain
from ibvpy.mesh.fe_domain import \
    FEDomain

from ibvpy.mesh.fe_spring_array import FESpringArray
from mathkit.mfn import MFnLineArray

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, shape, \
    cos, sin, arctan, where, abs, all, any, diag, \
    argsort, sum, transpose, vstack, append, arange, \
    hstack, reshape, shape, size, zeros_like, cumsum, \
    dsplit, copy, dtype, sort, ones_like, unique

from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos, pi as Pi

from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

# Interpolation
from scipy.interpolate import Rbf

from geo_column import GEOColumn

from simiter.sim_pstudy import\
    ISimModel, SimOut, SimPStudy, SimArray#, SimArrayView

from hp_shell import HPShell

from mush_roof_model import MushRoofModel

from matplotlib.pyplot import bar, show, axhline, ion, ioff, xlabel, ylabel, title, figure, savefig

import csv

def construct_fe_domain( fe_domain_list ):
    _sdomain = FEDomain()
    for d in fe_domain_list:
        if isinstance( d, FEGrid ):
            fe_rgrid = FERefinementGrid( domain = _sdomain,
                                         fets_eval = d.fets_eval )
            d.level = fe_rgrid
        elif isinstance( d, FESubDomain ):
            d.domain = _sdomain
        else:
            raise TypeError, 'The list can contain only FEGrid or FERefinementGrid'
    return _sdomain

class MRtwo( MushRoofModel ):

    implements( ISimModel )
    mushroof_part = 'one'
    #----------------------------------------------------
    # elements
    #----------------------------------------------------

    n_elems_xy_quarter = Int( 20, input = True )#, ps_levels = [4, 16, 5] )
    n_elems_z = Int( 1, input = True )#, ps_levels = [1, 2, 1] )
    n_elems_col_z = Int( 10  , input = True )#, ps_levels = [5, 20, 3 ] )
    n_elems_col_xy = Int( 2 , input = True )#, ps_levels = [2,4,1]   )
    shift_elems = True


    vtk_r = Float( 0.95, input = True )
    #default column
    #default roof
    fe_roof = Instance( ( FETSEval ), depends_on = '+ps_levels, +input' )
    def _fe_roof_default( self ):
        fets = self.fe_quad_serendipity_roof
        fets.vtk_r *= self.vtk_r
        return fets

    #default plate
    fe_plate = Instance( ( FETSEval ), depends_on = '+ps_levels, +input' )
    def _fe_plate_default ( self ):
        fets = self.fe_quad_serendipity_plate
        fets.ngp_r = 3
        fets.ngp_s = 3
        fets.ngp_t = 3
        fets.vtk_r *= self.vtk_r
        return fets


    #default column
    fe_column = Instance( ( FETSEval ), transient = True , depends_on = '+ps_levels, +input' )
    def _fe_column_default( self ):
        fets = self.fe_quad_serendipity_column
        fets.vtk_r *= self.vtk_r
        return fets


    #----------------------------------------------------
    # geometric transformation
    #----------------------------------------------------

    hp_shells = Property( depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_hp_shells( self ):
        X_list = [ [0, 0, 0], [8, 0, 0] ]
        return [
                 HPShell( t_shell = self.t_shell,
                         X0 = X,
                         delta_h_scalefactor = self.delta_h_scalefactor,
                         n_elems_xy = self.n_elems_xy,
                         n_elems_z = self.n_elems_z,
                         shift_elems = self.shift_elems,
                         const_edge_elem = self.const_edge_elem,
                         shift_array = self.shift_array,
                         mushroof_part = self.mushroof_part )
                 for X in X_list ]

    columns = Property( Instance( GEOColumn ) , depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_columns( self ):
        # column heights for a total heights of constant 4m
        #
        ch = self.h_col / self.delta_h_scalefactor
        X_list = [[ 4., 4., -ch ],
                  [ 12., 4., -ch ], ]
        return [ GEOColumn( width_top_col = self.width_top_col,
                            X0 = X,
                            width_bottom_col = self.width_bottom_col,
                            h_col = ch - self.t_plate )
                            for X in X_list ]

    plates = Property( Instance( GEOColumn ) , depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_plates( self ):
        ch = self.t_plate
        X_list = [[ 4., 4., -ch ],
                  [ 12., 4., -ch ], ]
        return [GEOColumn( width_top = self.width_top_col,
                          width_bottom = self.width_top_col,
                          X0 = X,
                          h_col = self.t_plate )
                          for X in X_list ]
    #----------------------------------------------------
    # grid 
    #----------------------------------------------------

    fe_grid_roofs = Property( Instance( FEGrid ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_fe_grid_roofs( self ):
        return [
                FEGrid( coord_min = ( 0.0, 0.0, 0.0 ),
                       coord_max = ( 1.0, 1.0, 1.0 ),
                       geo_transform = hp_shell,
                       shape = ( self.n_elems_xy, self.n_elems_xy, self.n_elems_z ),
                       fets_eval = self.fe_roof )
                for hp_shell in self.hp_shells ]

    fe_grid_columns = Property( Instance( FEGrid ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_fe_grid_columns( self ):
        return [ FEGrid( coord_min = ( 0.0, 0.0, 0.0 ),
                         coord_max = ( 1.0, 1.0, 1.0 ),
                         geo_transform = column,
                         shape = ( self.n_elems_col_xy, self.n_elems_col_xy, self.n_elems_col_z ),
                         fets_eval = self.fe_column )
                         for column in self.columns ]

    fe_grid_plates = Property( Instance( FEGrid ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_fe_grid_plates( self ):
        return [
                FEGrid( coord_min = ( 0.0, 0.0, 0.0 ),
                       coord_max = ( 1.0, 1.0, 1.0 ),
                       geo_transform = plates,
                       shape = ( self.n_elems_col_xy, self.n_elems_col_xy, 2 ),
                       fets_eval = self.fe_plate )
                for plates in self.plates ]

    # grid where all dofs are zero to apply boundary conditions at symmetric axes
    #
    fe_grid_sym_axes = Property( Instance( FEGrid ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_fe_grid_sym_axes( self ):
        return [FEGrid( coord_min = ( 0.0, 8.0, 1.06 ),
                      coord_max = ( 16.0, 8.02, 1.062, ),
                      shape = ( self.dofs_r0_link_top[0].shape[0], 1, 1 ),
                      fets_eval = FETS3D8H( mats_eval = MATS3DElastic( E = 1, nu = 0 ) ) )]

    #----------------------------------------------------
    # evaluation
    #----------------------------------------------------


    def peval( self ):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''

        # displacement from sorted slices
        #
        U = self.tloop.eval()

        if self.link_type == 'exc_V_ip':

            N_ip = hstack( ( self.r0_r1_N_ip.get_spring_forces( U ),
                             self.sym_N_ip.get_spring_forces( U ) ) )
            V_ip = zeros_like( N_ip )
            V_op = hstack( ( self.r0_r1_V_op.get_spring_forces( U ),
                             zeros_like( self.sym_N_ip.get_spring_forces( U ) ) ) )

            F_export = vstack( ( N_ip, V_ip, V_op ) ).transpose()


            X_F_export = vstack( ( self.dofs_r0_link_right[1],
                                 self.dofs_link_top[1] ) )

        if self.link_type == 'inc_V_ip':
            N_ip = hstack( ( self.r0_r1_N_ip.get_spring_forces( U ),
                             self.sym_N_ip.get_spring_forces( U ) ) )
            V_ip = hstack( ( self.r0_r1_V_ip.get_spring_forces( U ),
                             zeros_like( self.sym_N_ip.get_spring_forces( U ) ) ) )
            V_op = hstack( ( self.r0_r1_V_op.get_spring_forces( U ),
                             zeros_like( self.sym_N_ip.get_spring_forces( U ) ) ) )

            F_export = vstack( ( N_ip, V_ip, V_op ) ).transpose()


            X_F_export = vstack( ( self.dofs_r0_link_right[1],
                                 self.dofs_link_top[1] ) )


        self.X_F_export = X_F_export
        self.F_export = F_export

        self.U_export = U[self.dofs_u[0]].reshape( -1, 3 )
        self.X_U_export = self.dofs_u[1].reshape( -1, 3 )

        return array( [ ] )


    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [  SimOut( name = 'U_z', unit = 'm' ), ]


    X_F_export = Array  # X, Y, Z
    F_export = Array  # Nip, Vip, Vop order defined in peval


    def export_int_force_data( self, filename = 'F_int_data.csv' ):
        '''exports X_F_export and F_export data to csv - worksheet
        '''

        print '*** writing hinge force data data to file,', filename, ' ***'

        X_data = self.X_F_export.reshape( -1, 3 )
        F_data = self.F_export.reshape( -1, 3 )

        data = hstack( ( X_data, F_data ) )
        file = open( filename, 'w' )

        writer = csv.writer( file, delimiter = ";", lineterminator = "\n" )
        writer.writerow( ['x[m]', 'y[m]', 'z[m]', 'N_ip[MN]', 'V_ip[MN]', 'V_op[MN]' ] )
        writer.writerows( data )

        file = file.close()
        return

    X_U_export = Array
    U_export = Array


    def export_edge_u_data( self, filename = 'U_data.csv' ):
        '''exports X_U_export and U_export data to csv - worksheet
        '''
        print '*** writing displacement data to file,', filename, ' ***'

        X_data = self.X_U_export.reshape( -1, 3 )
        U_data = self.U_export.reshape( -1, 3 )

        data = hstack( ( X_data, U_data ) )
        file = open( filename, 'w' )

        writer = csv.writer( file, delimiter = ";", lineterminator = "\n" )
        writer.writerow( ['X[m]', 'Y[m]', 'Z[m]', 'u-x[m]', 'u-y[m]', 'u-z[m]' ] )
        writer.writerows( data )

        file = file.close()
        return


    X_array_linked = Property()
    def _get_X_array_linked( self ):
        """gets the X coordinates of all linked elements depending on shift_array
        """

        linked_array = vstack( ( self.shift_array[:self.not_linked_elem],
                               self.shift_array[self.not_linked_elem + 1:] ) )

        #global coordinates from local coordinates of one quarter
        #
        x_array_linked = hstack( ( 4.0 - linked_array[:, 0],
                                4.0 + linked_array[:, 0],
                                12.0 - linked_array[:, 0],
                                12.0 + linked_array[:, 0] ) )
        return sort( x_array_linked )

    Y_array_linked = Property()
    def _get_Y_array_linked( self ):
        """gets the Y coordinates of all linked elements depending on shift_array
        """

        linked_array = vstack( ( self.shift_array[:self.not_linked_elem],
                               self.shift_array[self.not_linked_elem + 1:] ) )

        #global coordinates from local coordinates of one quarter
        #
        y_array_linked = hstack( ( 4 - linked_array[:, 1],
                                 4 + linked_array[:, 1], ) )
        return sort( y_array_linked )


    #----------------------------------------------------
    # boundaries
    #----------------------------------------------------

    bc_plate_roof_link_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_bc_plate_roof_link_list( self ):
        '''
        links all plate corner nodes of each elements to the adjacent elements of the roof
        '''

        bc_plate_roof_link_list = []
        n_from_center = 1 + self.not_linked_elem

        for roof, plate in zip( self.fe_grid_roofs, self.fe_grid_plates ):
            slice_1 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = roof[self.n_elems_xy_quarter - n_from_center ,
                                               self.n_elems_xy_quarter, 0,
                                               0, 0, 0 ],
                                  link_slice = plate[ 0 , 0 , -1, 0, 0, -1], link_coeffs = [1.0],
                                  value = 0. )]
            slice_2 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = roof[self.n_elems_xy_quarter,
                                               self.n_elems_xy_quarter - n_from_center, 0,
                                               0, 0, 0 ],
                                  link_slice = plate[ -1, 0, -1, -1, 0, -1], link_coeffs = [1.0],
                                  value = 0. )]

            slice_3 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = roof[self.n_elems_xy_quarter + n_from_center,
                                               self.n_elems_xy_quarter, 0,
                                               0, 0, 0 ],
                                  link_slice = plate[ -1 , -1 , -1, -1, -1, -1], link_coeffs = [1.0],
                                  value = 0. )]

            slice_4 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = roof[self.n_elems_xy_quarter ,
                                               self.n_elems_xy_quarter + n_from_center, 0,
                                               0, 0, 0 ],
                                  link_slice = plate[ 0 , -1 , -1, 0, -1, -1], link_coeffs = [1.0],
                                  value = 0. )]
            bc_plate_roof_link_list = bc_plate_roof_link_list + \
                              slice_1 + slice_2 + slice_3 + slice_4

        return bc_plate_roof_link_list


    bc_roof_top_roof_low_link_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_bc_roof_top_roof_low_link_list( self ):
        '''
        links all plate corner nodes of each elements to the adjacent elements of the roof
        '''
        n_from_center = 1 + self.not_linked_elem
        bc_roof_top_roof_low_link_list = []
        for roof, plate in zip( self.fe_grid_roofs, self.fe_grid_plates ):
            slice_1 = [BCSlice( var = 'u'  , dims = [ 2],
                                  link_slice = roof[self.n_elems_xy_quarter - n_from_center ,
                                               self.n_elems_xy_quarter, 0,
                                               0, 0, 0 ],
                                  slice = roof[self.n_elems_xy_quarter - n_from_center ,
                                                    self.n_elems_xy_quarter, -1,
                                                    0, 0, -1 ], link_coeffs = [1.0],
                                  value = 0. )]
            slice_2 = [BCSlice( var = 'u'  , dims = [ 2],
                                  link_slice = roof[self.n_elems_xy_quarter,
                                               self.n_elems_xy_quarter - n_from_center, 0,
                                               0, 0, 0 ],
                                  slice = roof[self.n_elems_xy_quarter ,
                                                    self.n_elems_xy_quarter - n_from_center, -1,
                                                    0, 0, -1 ], link_coeffs = [1.0],
                                  value = 0. )]

            slice_3 = [BCSlice( var = 'u'  , dims = [ 2],
                                  link_slice = roof[self.n_elems_xy_quarter + n_from_center,
                                               self.n_elems_xy_quarter, 0,
                                               0, 0, 0 ],
                                  slice = roof[self.n_elems_xy_quarter + n_from_center ,
                                                    self.n_elems_xy_quarter, -1,
                                                    0, 0, -1 ], link_coeffs = [1.0],
                                  value = 0. )]

            slice_4 = [BCSlice( var = 'u'  , dims = [ 2],
                                link_slice = roof[self.n_elems_xy_quarter ,
                                             self.n_elems_xy_quarter + n_from_center, 0,
                                             0, 0, 0 ],
                                slice = roof[self.n_elems_xy_quarter  ,
                                                    self.n_elems_xy_quarter + n_from_center, -1,
                                                    0, 0, -1 ], link_coeffs = [1.0],
                                value = 0. )]

            slice_5 = [BCSlice( var = 'u'  , dims = [ 2],
                                link_slice = roof[self.n_elems_xy_quarter ,
                                             self.n_elems_xy_quarter , 0,
                                             0, 0, 0 ],
                                slice = roof[self.n_elems_xy_quarter  ,
                                                    self.n_elems_xy_quarter, -1,
                                                    0, 0, -1 ], link_coeffs = [1.0],
                                value = 0. )]


            bc_roof_top_roof_low_link_list = bc_roof_top_roof_low_link_list + \
                                             slice_1 + slice_2 + slice_3 + slice_4 + slice_5

        return bc_roof_top_roof_low_link_list



    bc_plate_column_link_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_bc_plate_column_link_list( self ):
        '''
        links all column nodes to plate nodes but only corner nodes of each element
        for 2x2 elements 9 elements are linked
        '''
        bc_plate_column_link_list = []
        for column, plate in zip( self.fe_grid_columns, self.fe_grid_plates ):

            slice_1 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = plate[:, :, 0, -1, -1, 0 ],
                                  link_slice = column[ :, :, -1 , -1, -1, -1], link_coeffs = [1.0],
                                  value = 0. )]

            slice_2 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = plate[:, :, 0, 0, 0, 0 ],
                                  link_slice = column[ :, :, -1 , 0, 0, -1], link_coeffs = [1.0],
                                  value = 0. )]

            slice_3 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = plate[:, :, 0, 0, -1, 0 ],
                                  link_slice = column[ :, :, -1 , 0, -1, -1], link_coeffs = [1.0],
                                  value = 0. )]

            slice_4 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = plate[:, :, 0, -1, 0, 0 ],
                                  link_slice = column[ :, :, -1 , -1, 0, -1], link_coeffs = [1.0],
                                  value = 0. )]

            bc_plate_column_link_list = bc_plate_column_link_list + \
                                        slice_1 + slice_2 + slice_3 + slice_4

        return bc_plate_column_link_list
#        return [BCSlice( var = 'u'  , dims = [0, 1, 2],
#                         slice = plate[:,:,0,:,:, 0 ],
#                         link_slice = column[ :,:,-1 ,:,:,-1], link_coeffs = [1.0], value = 0. )]


    bc_col_clamped_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_bc_col_clamped_list( self ):
        '''clamped support of the columns at the bottom
        '''
        constraint = []
        for column in self.fe_grid_columns:
            constraint = constraint + [ BCSlice( var = 'u', dims = [0, 1, 2],
                                                slice = column[ :, :, 0, :, :, 0 ],
                                                value = 0.0 ) ]
        return constraint

    bc_col_link_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_bc_col_link_list( self ):
        '''links column top to roof
        ( all column corner nodes of each elements to the adjacent elements of the roof )
        '''
        bc_col_link_list = []

        n_from_center = 1 + self.not_linked_elem


        for roof, column in zip( self.fe_grid_roofs, self.fe_grid_columns ):
            slice_1 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = roof[self.n_elems_xy_quarter - n_from_center ,
                                               self.n_elems_xy_quarter, 0,
                                               0, 0, 0 ],
                                  link_slice = column[ 0 , 0 , -1, 0, 0, -1], link_coeffs = [1.0],
                                  value = 0. )]
            slice_2 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = roof[self.n_elems_xy_quarter,
                                               self.n_elems_xy_quarter - n_from_center, 0,
                                               0, 0, 0 ],
                                  link_slice = column[ -1, 0, -1, -1, 0, -1], link_coeffs = [1.0],
                                  value = 0. )]

            slice_3 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = roof[self.n_elems_xy_quarter + n_from_center,
                                               self.n_elems_xy_quarter, 0,
                                               0, 0, 0 ],
                                  link_slice = column[ -1 , -1 , -1, -1, -1, -1], link_coeffs = [1.0],
                                  value = 0. )]

            slice_4 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                                  slice = roof[self.n_elems_xy_quarter ,
                                               self.n_elems_xy_quarter + n_from_center, 0,
                                               0, 0, 0 ],
                                  link_slice = column[ 0 , -1 , -1, 0, -1, -1], link_coeffs = [1.0],
                                  value = 0. )]


            bc_col_link_list = bc_col_link_list + slice_1 + slice_2 + slice_3 + slice_4

        return bc_col_link_list


    # column support as hinge at the bottom
    # 
    bc_col_hinge_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_bc_col_hinge_list( self ):
        '''column support as hinge at the bottom
        '''
        constraint = []
        n_el_col_xy_half = self.n_elems_col_xy / 2.0
        for column in self.fe_grid_columns:
            for i in range( 0, n_el_col_xy_half ):
                dof_const = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                             slice = column[i , 0 , 0, 0, 0, 0 ],
                             link_slice = column[-1 - i , -1, 0, -1 , -1, 0], link_coeffs = [-1.0],
                             value = 0.0 )]
                constraint = constraint + dof_const
            for i in range( 0, n_el_col_xy_half ):
                dof_const = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                             slice = column[0 , -1 - i , 0, 0, -1, 0 ],
                             link_slice = column[-1, i , 0, -1, 0 , 0], link_coeffs = [-1.0],
                             value = 0.0 )]
                constraint = constraint + dof_const

        return constraint


    #-------------------------------------
    # arbitrary linking 
    #-------------------------------------


    shift_array = Array( value = [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1], ] , input = True )

    link_type = Enum( 'exc_V_ip', 'inc_V_ip' )

    linked_dims = {'exc_V_ip':[0, 2],
                   'inc_V_ip':[0, 1, 2]}#dims of link, no in plane shear forces for exc V_ip

    not_linked_elem = Int( 1 )  # not linked element depends on shift array     

    select_hinge_pos = Bool( True )
    link_elem_array = Property( Array, depends_on = '+ps_levels, +input' )
    def _get_link_elem_array( self ):
        if self.select_hinge_pos == False:
            return arange( 0, self.n_elems_xy_quarter, 1, dtype = int )
        if self.select_hinge_pos == True:
            e_shifted = cumsum( self.shift_array[:, 2] )
            # create link list of elements 
            # always [:,:,0,0,0] not [:,:,:,-1,-1,-1] will be linked
            # without not_linked elem of shift array because this is only needed for column link
            #
            idx_e_shifted = arange( size( e_shifted ) )

            e_linked = append( e_shifted[:self.not_linked_elem], e_shifted[self.not_linked_elem + 1:] )
            return array( append( self.n_elems_xy_quarter - e_linked,
                           self.n_elems_xy_quarter + e_linked ), dtype = int )


    bc_sym_axes_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_bc_sym_axes_list( self ):
        ''' elements needed for sym axes spring are constrained in all directions
        '''
        fe_sym = self.fe_grid_sym_axes[0]
        return [BCSlice( var = 'u', dims = [0, 1, 2],
                         slice = fe_sym[ : , :, :, :, :, :],
                         value = 0.0 )]


#    bc_roof_link_list = Property( List, depends_on = '+ps_levels, +input' )
#    @cached_property
#    def _get_bc_roof_link_list( self ):
#        '''link roof to roof
#        and boundary conditions for the symmetry assumption
#        '''
#        roof_0, roof_1 = self.fe_grid_roofs
#        if self.select_hinge_pos == False:
#
#
#            # constraint at symmetric axes for roof 0
#            #
#            constraint_0 = [BCSlice( var = 'u', dims = [1],
#                                     # get only every corner node
#                                     slice = roof_0[ : , -1, -1, 0, -1, -1],
#                                     value = 0.0 )]
#
#            # constraint at symmetric axes for roof 1
#            #
#            constraint_1 = [BCSlice( var = 'u', dims = [1],
#                                     slice = roof_1[ : , -1, -1, -1, -1, -1],
#                                     value = 0.0 )]
#
#            # constraint for center node connecting roof 0 and 1 at synmmmetric axes
#            #
#            constraint_0_1 = [BCSlice( var = 'u', dims = [1],
#                                     slice = roof_0[-1 , -1, -1, -1, -1, -1],
#                                     value = 0.0 )]
#
#
#            # link of roof 0 and 1
#            #
#
#            link_0_1 = [BCSlice( var = 'u'  , dims = self.linked_dims[self.link_type],
#                                         slice = roof_0[-1, :, -1, -1, 0, -1 ],
#                                         link_slice = roof_1[ 0, :, -1, 0, 0, -1],
#                                         link_coeffs = [1.0],
#                                         value = 0. )]
#            link_0_1_last = [BCSlice( var = 'u'  , dims = self.linked_dims[self.link_type],
#                                         slice = roof_0[-1, -1, -1, -1, -1, -1 ],
#                                         link_slice = roof_1[ 0, -1, -1, 0, -1, -1],
#                                         link_coeffs = [1.0],
#                                         value = 0. )]
#
#            bc_roof_link_list = constraint_0 + constraint_1 + constraint_0_1 + \
#                                link_0_1 + link_0_1_last
#
#        if self.select_hinge_pos == True:
#            constraint_list = []
#            for elem in self.link_elem_array:
#                # link
#                #     
#                link_0_1 = [BCSlice( var = 'u'  , dims = self.linked_dims[self.link_type],
#                                         slice = roof_0[-1, elem, -1, -1, 0, -1 ],
#                                         link_slice = roof_1[ 0, elem, -1, 0, 0, -1],
#                                         link_coeffs = [1.0],
#                                         value = 0. )]
#                # symmetric constraint
#                #
#                const_0 = [BCSlice( var = 'u'  , dims = [1],
#                                         slice = roof_0[elem, -1, -1, 0, -1, -1 ],
#                                         value = 0. )]
#                const_1 = [BCSlice( var = 'u'  , dims = [1],
#                                         slice = roof_1[elem, -1, -1, 0, -1, -1 ],
#                                         value = 0. )]
#
#                constraint_list = constraint_list + link_0_1 + const_0 + const_1
#
#            bc_roof_link_list = constraint_list
#
#        return bc_roof_link_list


    lc = Enum( 'lc_g',
               'lc_s_sym', 'lc_s_asym',
               'lc_w_pos', 'lc_w_neg', 'lc_w_asym',
               'lc_shrink', input = True )

    lc_dict = Property( Dict )
    def _get_lc_dict( self ):
        return {'lc_g' : self.lc_g_list,
                'lc_s_sym' : self.lc_s_sym_list,
                'lc_s_asym' : self.lc_s_asym_list,
                'lc_w_pos' : self.lc_w_pos_list,
                'lc_w_neg' : self.lc_w_neg_list,
                'lc_w_asym' : self.lc_w_asym_list,
                'lc_shrink' : [] }


    lc_list = Property( List )
    def _get_lc_list( self ):
        return self.lc_dict[self.lc]

    lc_g_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_lc_g_list( self ):
        '''loading case: dead loads ( self - weight and 2cm add thickness )
        '''
        #slices   
        #
        roof_0 , roof_1 = self.fe_grid_roofs
        col_0, col_1 = self.fe_grid_columns
        upper_surf_0 = roof_0[:, :, -1, :, :, -1]
        upper_surf_1 = roof_1[:, :, -1, :, :, -1]
        bottom_edge_roof_0 = roof_0[:, 0, -1, :, 0, -1]
        bottom_edge_roof_1 = roof_1[:, 0, -1, :, 0, -1]
        left_edge_roof_0 = roof_0[0, :, -1, 0, :, -1]
        right_edge_roof_1 = roof_1[-1, :, -1, -1, :, -1]

        # loads in global z- direction
        #
        material_density_roof = -22.4e-3    # [MN/m^3]
        material_density_column = -26e-3    # [MN/m^3]
        add_t = 0.02 * material_density_roof# [MN/m^2]
        add_surf_load = -0.2e-3             # [MN/m^2]
        edge_load = -0.35e-3                # [MN/m]

        return [ BCSlice( var = 'f', value = material_density_roof, dims = [2],
                          integ_domain = 'global',
                          slice = roof_0[:, :, :, :, :, :] ),
                 BCSlice( var = 'f', value = material_density_roof, dims = [2],
                          integ_domain = 'global',
                          slice = roof_1[:, :, :, :, :, :] ),

                 BCSlice( var = 'f', value = material_density_column, dims = [2],
                          integ_domain = 'global',
                          slice = col_0[:, :, :, :, :, :] ),
                 BCSlice( var = 'f', value = material_density_column, dims = [2],
                          integ_domain = 'global',
                          slice = col_1[:, :, :, :, :, :] ),

                 BCSlice( var = 'f', value = add_surf_load + add_t, dims = [2],
                          integ_domain = 'global',
                          slice = upper_surf_0 ),
                 BCSlice( var = 'f', value = add_surf_load + add_t, dims = [2],
                          integ_domain = 'global',
                          slice = upper_surf_1 ),

                 BCSlice( var = 'f', value = edge_load, dims = [2],
                          integ_domain = 'global',
                          slice = bottom_edge_roof_0 ),
                 BCSlice( var = 'f', value = edge_load, dims = [2],
                          integ_domain = 'global',
                          slice = bottom_edge_roof_1 ),

                 BCSlice( var = 'f', value = edge_load, dims = [2],
                          integ_domain = 'global',
                          slice = left_edge_roof_0 ),
                 BCSlice( var = 'f', value = edge_load, dims = [2],
                          integ_domain = 'global',
                          slice = right_edge_roof_1 )]

    lc_s_sym_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_lc_s_sym_list( self ):
        '''loading case: snow symmetric
        '''
        #slices   
        roof_0 , roof_1 = self.fe_grid_roofs
        upper_surf_0 = roof_0[:, :, -1, :, :, -1]
        upper_surf_1 = roof_1[:, :, -1, :, :, -1]
        #loads
        s_sysm = -0.85e-3 # [MN/m^2] 
        return [ BCSlice( var = 'f', value = s_sysm, dims = [2],
                          integ_domain = 'global',
                          slice = upper_surf_0 ),
                 BCSlice( var = 'f', value = s_sysm, dims = [2],
                          integ_domain = 'global',
                          slice = upper_surf_1 ), ]

    lc_s_asym_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_lc_s_asym_list( self ):
        '''loading case: snow asymmetric
        '''
        #slices   
        roof_0 , roof_1 = self.fe_grid_roofs
        upper_surf_0_l = roof_0[:self.n_elems_xy_quarter, :, -1, :, :, -1]
        upper_surf_0_r = roof_0[self.n_elems_xy_quarter:, :, -1, :, :, -1]

        upper_surf_1_l = roof_1[:self.n_elems_xy_quarter, :, -1, :, :, -1]
        upper_surf_1_r = roof_1[self.n_elems_xy_quarter:, :, -1, :, :, -1]

        #loads
        s_left = -0.425e-3 # [MN/m^2] 
        s_right = -0.85e-3 # [MN/m^2] 

        return [# left side in z-direction 
                BCSlice( var = 'f', value = s_left, dims = [2],
                          integ_domain = 'global',
                          slice = upper_surf_0_l ),
                BCSlice( var = 'f', value = s_left, dims = [2],
                          integ_domain = 'global',
                          slice = upper_surf_1_l ),
                # right side in z-direction            
                BCSlice( var = 'f', value = s_right, dims = [2],
                          integ_domain = 'global',
                          slice = upper_surf_0_r ),
                BCSlice( var = 'f', value = s_right, dims = [2],
                          integ_domain = 'global',
                          slice = upper_surf_1_r ), ]

    lc_w_asym_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_lc_w_asym_list( self ):
        '''loading case: wind asymmetric
        '''
        # slices
        r0, r1 = self.fe_grid_roofs
        n_e_q = self.n_elems_xy_quarter
        left_edge = r0[ 0, :, -1, 0, :, -1 ]
        right_edge = r1[-1, :, -1, 0, :, -1 ]
        bottom_edge_0_left = r0[ :n_e_q, 0, -1, :, 0, -1 ]
        bottom_edge_0_right = r0[n_e_q:, 0, -1, :, 0, -1]
        bottom_edge_1 = r1[ :, 0, -1, :, 0, -1 ]
        face_0_left = r0[ :n_e_q, :, -1, :, :, -1 ]
        face_0_right = r0[ n_e_q:, :, -1, :, :, -1 ]
        face_1_left = r1[ :n_e_q, :, -1, :, :, -1 ]
        face_1_right = r1[ n_e_q:, :, -1, :, :, -1 ]

        # loads
        w_left_edge = 0.91e-3 # [MN/m]
        w_right_edge = 0.39e-3 # [MN/m]
        w_bottom_edge_0_left = -1.56e-3 # [Mn/m] tension
        w_bottom_edge_0_right = -1.04e-3 # [Mn/m] tension
        w_bottom_edge_1 = -0.65e-3 # [Mn/m]tension
        w_face_0_left = 0.89e-3 # MN/m^2
        w_face_0_right = -0.13e-3 # MN/m^2
        w_face_1_left = 0.72e-3 # MN/m^2
        w_face_1_right = w_face_0_right

        return [   # left edge - x direction
                   BCSlice( var = 'f'  , dims = [0],
                            slice = left_edge,
                            value = w_left_edge ),
                   # right edge - x direction
                   BCSlice( var = 'f'  , dims = [0],
                            slice = right_edge,
                            value = w_right_edge ),
                   # bottom edge - y direction
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_0_left,
                            value = w_bottom_edge_0_left ),
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_0_right,
                            value = w_bottom_edge_0_right ),
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_1,
                            value = w_bottom_edge_1 ),
                   # upper face left - left 
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_0_left,
                            value = w_face_0_left ),
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_0_right,
                            value = w_face_0_right ),
                   # upper face left - left 
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_1_left,
                            value = w_face_1_left ),
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_1_right,
                            value = w_face_1_right ),
                ]

    lc_w_pos_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_lc_w_pos_list( self ):
        '''loading case: wind symmetric pressure
        '''
        # slices
        r0, r1 = self.fe_grid_roofs
        n_e_q = self.n_elems_xy_quarter
        left_edge = r0[ 0, :, -1, 0, :, -1 ]
        right_edge = r1[-1, :, -1, 0, :, -1 ]
        bottom_edge_0_left = r0[ :n_e_q, 0, -1, :, 0, -1 ]
        bottom_edge_0_right = r0[n_e_q:, 0, -1, :, 0, -1]
        bottom_edge_1 = r1[ :, 0, -1, :, 0, -1 ]
        face_0_left = r0[ :n_e_q, :, -1, :, :, -1 ]
        face_0_right = r0[ n_e_q:, :, -1, :, :, -1 ]
        face_1_left = r1[ :n_e_q, :, -1, :, :, -1 ]
        face_1_right = r1[ n_e_q:, :, -1, :, :, -1 ]

        # loads
        w_left_edge = 0.91e-3 # [MN/m] pressure
        w_right_edge = 0.39e-3 # [MN/m] tension
        w_bottom_edge_0_left = -1.56e-3 # [Mn/m] tension
        w_bottom_edge_0_rigth = -1.04e-3 # [Mn/m] tension
        w_bottom_edge_1 = -0.65e-3 # [Mn/m]tension
        w_face_0_left = 0.89e-3 # MN/m^2 tension
        w_face_0_right = -0.13e-3 # MN/m^2 pressure
        w_face_1_left = w_face_0_right # MN/m^2 pressure
        w_face_1_right = w_face_0_right # MN/m^2 pressure

        return [   # left edge - x direction
                   BCSlice( var = 'f'  , dims = [0],
                            slice = left_edge,
                            value = w_left_edge ),
                   # right edge - x direction
                   BCSlice( var = 'f'  , dims = [0],
                            slice = right_edge,
                            value = w_right_edge ),
                   # bottom edge - y direction
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_0_left,
                            value = w_bottom_edge_0_left ),
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_0_right,
                            value = w_bottom_edge_0_rigth ),
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_1,
                            value = w_bottom_edge_1 ),
                   # upper face left - left 
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_0_left,
                            value = w_face_0_left ),
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_0_right,
                            value = w_face_0_right ),
                   # upper face left - left 
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_1_left,
                            value = w_face_1_left ),
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_1_right,
                            value = w_face_1_right ),
                ]

    lc_w_neg_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_lc_w_neg_list( self ):
        '''loading case: wind symmetric suction
        '''
        # slices
        r0, r1 = self.fe_grid_roofs
        n_e_q = self.n_elems_xy_quarter
        left_edge = r0[ 0, :, -1, 0, :, -1 ]
        right_edge = r1[-1, :, -1, 0, :, -1 ]
        bottom_edge_0_left = r0[ :n_e_q, 0, -1, :, 0, -1 ]
        bottom_edge_0_right = r0[n_e_q:, 0, -1, :, 0, -1]
        bottom_edge_1 = r1[ :, 0, -1, :, 0, -1 ]
        face_0_left = r0[ :n_e_q, :, -1, :, :, -1 ]
        face_0_right = r0[ n_e_q:, :, -1, :, :, -1 ]
        face_1_left = r1[ :n_e_q, :, -1, :, :, -1 ]
        face_1_right = r1[ n_e_q:, :, -1, :, :, -1 ]

        # loads
        w_left_edge = 0.91e-3 # [MN/m] pressure
        w_right_edge = 0.39e-3 # [MN/m] tension
        w_bottom_edge_0_left = -1.56e-3 # [Mn/m] tension
        w_bottom_edge_0_rigth = -1.04e-3 # [Mn/m] tension
        w_bottom_edge_1 = -0.65e-3 # [Mn/m]tension
        w_face_0_left = 0.89e-3 # [MN/m^2] tension
        w_face_0_right = 0.39e-3 # [MN/m^2] tension
        w_face_1_left = 0.72e-3 # [MN/m^2] tension
        w_face_1_right = 0.39e-3 # [MN/m^2]tension

        return [   # left edge - x direction
                   BCSlice( var = 'f'  , dims = [0],
                            slice = left_edge,
                            value = w_left_edge ),
                   # right edge - x direction
                   BCSlice( var = 'f'  , dims = [0],
                            slice = right_edge,
                            value = w_right_edge ),
                   # bottom edge - y direction
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_0_left,
                            value = w_bottom_edge_0_left ),
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_0_right,
                            value = w_bottom_edge_0_rigth ),
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_1,
                            value = w_bottom_edge_1 ),
                   # upper face left - left 
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_0_left,
                            value = w_face_0_left ),
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_0_right,
                            value = w_face_0_right ),
                   # upper face left - left 
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_1_left,
                            value = w_face_1_left ),
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_1_right,
                            value = w_face_1_right ),
                ]

    #===============================================================================
    # response tracer
    #===============================================================================

    rtrace_list = List
    def _rtrace_list_default( self ):
        return [  self.u]

#        return [  self.max_princ_stress, self.sig_app, self.u]


    # linked dofs

    dofs_r0_link_right = Property( List, depends_on = '+ps_levels, +input' )
    def _get_dofs_r0_link_right( self ):
        """linked_dofs at r0_right
            [:,0] for x, [:,1] for y, [:,2] for z,"""
        r0, r1 = self.fe_grid_roofs
        r_slice = r0[ -1, self.link_elem_array, -1, -1, 0, -1]
        return r_slice.dofs[:, 0, :], r_slice.dof_X[ :, 0, : ]

    dofs_r1_link_left = Property( List, depends_on = '+ps_levels, +input' )
    def _get_dofs_r1_link_left( self ):
        """linked_dofs at r1_left
            [:,0] for x, [:,1] for y, [:,2] for z,"""
        r0, r1 = self.fe_grid_roofs
        r_slice = r1[ 0, self.link_elem_array, -1, 0, 0, -1]
        return r_slice.dofs[:, 0, :], r_slice.dof_X[ :, 0, : ]

    dofs_r0_link_top = Property( List, depends_on = '+ps_levels, +input' )
    def _get_dofs_r0_link_top( self ):
        """ linked_dofs at r0_top , shape=(n,3) 
            [:,0] for x, [:,1] for y, [:,2] for z,"""
        r0, r1 = self.fe_grid_roofs
        r_slice = r0[ self.link_elem_array, -1, -1, 0, -1, -1]
        return r_slice.dofs[:, 0, :], r_slice.dof_X[ :, 0, : ]

    dofs_r1_link_top = Property( List, depends_on = '+ps_levels, +input' )
    def _get_dofs_r1_link_top( self ):
        """ linked_dofs at r0_top , shape=(n,3) 
            [:,0] for x, [:,1] for y, [:,2] for z,"""
        r0, r1 = self.fe_grid_roofs
        r_slice = r1[ self.link_elem_array, -1, -1, 0, -1, -1]
        return r_slice.dofs[:, 0, :], r_slice.dof_X[ :, 0, : ]

    dofs_link_top = Property( List, depends_on = '+ps_levels, +input' )
    def _get_dofs_link_top( self ):
        """ linked_dofs at r0_top , shape=(n,3) 
            [:,0] for x, [:,1] for y, [:,2] for z,"""
        r0, r1 = self.fe_grid_roofs
        r_slice_0 = r0[ self.link_elem_array, -1, -1, 0, -1, -1]
        r_slice_1 = r1[ self.link_elem_array, -1, -1, 0, -1, -1]
        return vstack( ( r_slice_0.dofs[:, 0, :], r_slice_1.dofs[:, 0, :] ) ), \
               vstack( ( r_slice_0.dof_X[ :, 0, : ], r_slice_1.dof_X[ :, 0, : ] ) )


    # outout dofs u
    #
    dofs_u = Property ( List, depends_on = '+ps_levels, +input' )
    def _get_dofs_u( self ):
        """get slice for u output = all edges"""
        r0, r1 = self.fe_grid_roofs
        r0_right = r0[ -1, :, -1, -1, :, -1]
        r0_left = r0[ 0, :, -1, 0, :, -1]
        r0_top = r0[ :, -1, -1, :, -1, -1]
        r0_bottom = r0[ :, 0, -1, :, -1, -1]

        r1_right = r1[ -1, :, -1, -1, :, -1]
        r1_left = r1[ 0, :, -1, 0, :, -1]
        r1_top = r1[ :, -1, -1, :, -1, -1]
        r1_bottom = r1[ :, 0, -1, :, -1, -1]
        return vstack( ( r0_right.dofs,
                         r0_left.dofs,
                         r0_top.dofs,
                         r0_bottom.dofs,
                         r1_right.dofs,
                         r1_left.dofs,
                         r1_top.dofs,
                         r1_bottom.dofs ) ), \
               vstack( ( r0_right.dof_X,
                         r0_left.dof_X,
                         r0_top.dof_X,
                         r0_bottom.dof_X,
                         r1_right.dof_X,
                         r1_left.dof_X,
                         r1_top.dof_X,
                         r1_bottom.dof_X ) )
        return



    # symmetric axes
    # 
    dofs_r0_link_sym = Property( List, depends_on = '+ps_levels, +input' )
    def _get_dofs_r0_link_sym( self ):
        fe_sym = self.fe_grid_sym_axes[0]
        return unique( fe_sym[:, :, :, :, 0, 0].dofs[:, 0, 0].reshape( -1 ) )

    dofs_r1_link_sym = Property( List, depends_on = '+ps_levels, +input' )
    def _get_dofs_r1_link_sym( self ):
        fe_sym = self.fe_grid_sym_axes[0]
        return unique( fe_sym[:, :, :, :, 0, 0].dofs[:, 0, 1].reshape( -1 ) )

    dofs_link_sym = Property( List, depends_on = '+ps_levels, +input' )
    def _get_dofs_link_sym( self ):
        fe_sym = self.fe_grid_sym_axes[0]
        return unique( fe_sym[:, :, :, :, 0, 0].dofs[:, 0, 1:] )


    #----------------------------------------------------
    # k_values
    #----------------------------------------------------


   # plate dimensions
   #

#    d_e_b = Float( 0.15, unit = 'm' ) # distance outer bolt to edge
#    d_i_b = Float( 0.3, unit = 'm' ) # distance between bolts
#    w_pl = Float( 0.15, unit = 'm' ) # width of plate
#    h_pl = Float( 0.03, unit = 'm' ) # distance between bolts
#    E_steel = Float( 210000.0, unit = 'MN/m^2' )



    k_n_ip = Property( Float )
    def _get_k_n_ip( self ):
#        A = self.w_pl * self.h_pl
#        d_11 = ( self.d_e_b + 0.5 ** 2 * self.d_i_b ) / ( self.E_steel * A )
#        return 1 / d_11 * 2.0
        return 20.0

    k_v_ip = Property( Float )
    def _get_k_v_ip( self ):
#        I_ip = self.h_pl * self.w_pl ** 3 / 12
#        d_11 = ( self.d_e_b ** 3 + self.d_e_b ** 2 * self.d_i_b ) / ( self.E_steel * I_ip )
#        return 1 / d_11 * 2.0
        return 4.0

    k_v_op = Property( Float )
    def _get_k_v_op( self ):
#        I_op = self.h_pl ** 3 * self.w_pl / 12
#        d_11 = ( self.d_e_b ** 3 + self.d_e_b ** 2 * self.d_i_b ) / ( self.E_steel * I_op )
#        return 1 / d_11 * 2.0
        return 14.0

    fe_domain = Property( depends_on = '+input' )
    def _get_fe_domain( self ):
        # Hack to get all domains together nedds to be done before linking for
        # boundary conditions plroblems
        # 
        roofs = self.fe_grid_roofs
        columns = self.fe_grid_columns
        plates = self.fe_grid_plates
        sym_axes = self.fe_grid_sym_axes
        return construct_fe_domain( roofs + columns + plates + sym_axes )


#  TO GET LINKING AUTOMATICALLY DIDNT WORK BEACUASE OF INPUT PROBLEMS


#    r0_r1_N_ip = Property( depends_on ='+input')
#    def _get_r0_r1_N_ip(self):
#        return  FESpringArray( domain = self.fe_domain,
#                               dofs_1 = self.dofs_r0_link_right[0][:, 0],
#                               dofs_2 = self.dofs_r1_link_left[0][:, 0], k_value = self.k_n_ip )
#    
#    r0_r1_V_ip = Property (depends_on ='+input')
#    def _get_r0_r1_V_ip(self):
#        return FESpringArray( domain = self.fe_domain,
#                              dofs_1 = self.dofs_r0_link_right[0][:, 1],
#                              dofs_2 = self.dofs_r1_link_left[0][:, 1], k_value = self.k_v_ip )
#
#    r0_r1_V_op = Property (depends_on ='+input')
#    def _get_r0_r1_V_op(self):
#            return FESpringArray( domain = self.fe_domain,
#                                  dofs_1 = self.dofs_r0_link_right[0][:, 2],
#                                  dofs_2 = self.dofs_r1_link_left[0][:, 2], k_value = self.k_v_op )
#      
#    sym_N_ip = Property( depends_on ='+input')
#    def _get_sym_N_ip(self):
#        return  FESpringArray( domain = self.fe_domain,
#                               dofs_1 = self.dofs_r0_link_top[0][:, 1],
#                               dofs_2 = self.dofs_r0_link_sym, k_value = self.k_n_ip / 2.0 )

#    r0_r1_N_ip = Property(FESpringArray, depends_on ='+input,+ps_levels')
#    def _get_r0_r1_N_ip(self):
#        return FESpringArray( domain = fe_domain,
#                       dofs_1 = self.dofs_r0_link_right[0][:, 0],
#                       dofs_2 = self.dofs_r1_link_left[0][:, 0], k_value = self.k_n_ip )

    #----------------------------------------------------
    # time loop
    #----------------------------------------------------

    tloop = Property( depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_tloop( self ):
        roofs = self.fe_grid_roofs
        fe_domain = self.fe_domain

        ts = TS( sdomain = fe_domain,
                 dof_resultants = True,
                 bcond_list =

                            # boundary conditions of the selected loading case
                            # e.g "lc_g_list", "lc_s_sym_list", etc.
                              self.lc_list +

                            # links (kinematic constraints and symmetric assumptions)
                              self.bc_sym_axes_list +
                              self.bc_plate_column_link_list +
                              self.bc_plate_roof_link_list +
                              self.bc_roof_top_roof_low_link_list +

                            # column at bottom as clamped support (compare only: hinged optional)                              
                              self.bc_col_clamped_list,
#                              self.bc_col_hinge_list ,

                 rtrace_list = self.rtrace_list

               )

        # Add the time-loop control
        #

        tloop = TLoop( tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline )
        # between roofs
        #



        self.r0_r1_N_ip = FESpringArray( domain = fe_domain,
                       dofs_1 = self.dofs_r0_link_right[0][:, 0],
                       dofs_2 = self.dofs_r1_link_left[0][:, 0], k_value = self.k_n_ip )
        if self.link_type == 'inc_V_ip':
            self.r0_r1_V_ip = FESpringArray( domain = fe_domain,
                           dofs_1 = self.dofs_r0_link_right[0][:, 1],
                           dofs_2 = self.dofs_r1_link_left[0][:, 1], k_value = self.k_v_ip )

        self.r0_r1_V_op = FESpringArray( domain = fe_domain,
                                         dofs_1 = self.dofs_r0_link_right[0][:, 2],
                                         dofs_2 = self.dofs_r1_link_left[0][:, 2], k_value = self.k_v_op )

        #sym_axes
        #
        self.sym_N_ip = FESpringArray( domain = fe_domain,
                                       dofs_1 = self.dofs_link_top[0][:, 1],
                                       dofs_2 = self.dofs_link_sym, k_value = self.k_n_ip * 2 )


        # output slices
#        #
#        r0, r1 = roofs
#        #slices r0
#        #
#        self.dofs_r0_right = r0[-1, :, -1, -1, :, -1].dofs
#        self.dofs_r0_left = r0[0, :, -1, 0, :, -1].dofs
#        self.dofs_r0_top = r0[:, -1, -1, :, -1, -1].dofs
#        self.dofs_r0_bottom = r0[:, 0, -1, :, 0, -1].dofs
#
#        #get all the dofs that are subjedcted to be linked
#        #
#
#        #slices r0 global coordinates
#        #
#        self.dof_X_r0_right = r0[-1, :, -1, -1, :, -1].dof_X
#        self.dof_X_r0_left = r0[0, :, -1, 0, :, -1].dof_X
#        self.dof_X_r0_top = r0[:, -1, -1, :, -1, -1].dof_X
#        self.dof_X_r0_bottom = r0[:, 0, -1, :, 0, -1].dof_X
#
#        #slices r1 dofs
#        #
#        self.dofs_r1_right = r1[-1, :, -1, -1, :, -1].dofs
#        self.dofs_r1_left = r1[0, :, -1, 0, :, -1].dofs
#        self.dofs_r1_top = r1[:, -1, -1, :, -1, -1].dofs
#        self.dofs_r1_bottom = r1[:, 0, -1, :, 0, -1].dofs
#
#        #slices r1 global coordinates
#        #
#        self.dof_X_r1_right = r1[-1, :, -1, -1, :, -1].dof_X
#        self.dof_X_r1_left = r1[0, :, -1, 0, :, -1].dof_X
#        self.dof_X_r1_top = r1[:, -1, -1, :, -1, -1].dof_X
#        self.dof_X_r1_bottom = r1[:, 0, -1, :, 0, -1].dof_X
#
#        #long slices over both domains
#        #
#        self.dofs_sym = vstack( ( self.dofs_r0_top, self.dofs_r1_top ) )
#        self.dofs_bottom = vstack( ( self.dofs_r0_bottom, self.dofs_r1_bottom ) )
#
#        self.dof_X_sym = vstack( ( self.dof_X_r0_top, self.dof_X_r1_top ) )
#        self.dof_X_bottom = vstack( ( self.dof_X_r0_bottom , self.dof_X_r1_bottom ) )

        return tloop



if __name__ == '__main__':


    sim_model = MRtwo()
    do = 'export_hinge'
    print sim_model.k_v_op
    print sim_model.k_v_ip
    print sim_model.k_n_ip


    if do == 'ui':

        print 'ui'

        sim_model.lc = 'lc_g'
        dofs_u = sim_model.dofs_u[0].shape

#        sim_model.shift_array = array( [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1],
#                                 [0.5, 0.5 , 2],
#                                 [1.5 , 1.5, 4],
#                                 [2.5, 2.5, 4],
#                                 [3.5, 3.5, 4],
#                                    ] )
#        print sim_model.shift_array

        sim_model.not_linked_elem = 0
#        sim_model.export_int_force_data()
#        print sim_model.dofs_link_sym
#        print sim_model.dofs_r0_link_right

        sim_model.peval()
#        from ibvpy.plugins.ibvpy_app import IBVPyApp
##
#        app = IBVPyApp( ibv_resource = sim_model )
#        app.main()

    if do == 'export_hinge':

        print 'export_hinge'


        #loading cases
        #
        lc_list = [
                   'lc_shrink' ,
#                   'lc_g',
#                   'lc_s_sym',
#                   'lc_s_asym',
#                   'lc_w_pos',
#                   'lc_w_neg',
#                   'lc_w_asym'
                   ]

        # link_cases
        #
        link_case_list = [
#                          'equal_25cm',
#                          'equal_50cm',
#                          'equal_100cm',
#                          'equal_200cm',
#                          'gap_100cm',
#                          'corner_2to4_25cm',
#                          'middle_0to3_25cm',
                          'middle_0to2_25cm'
                          ]

        # V_ip within hinge exists if  freedom is linked or no linked
        #
        link_type_list = [
#                          'exc_V_ip',
                          'inc_V_ip'
                          ]


        # shift_array_dict for the different types of links
        #       
        shift_array_dict = {'equal_25cm':
                           array( [[0.125 , 0.125, 1],
                                  [0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1],
                                  [0.375, 0.375, 1],
                                  [0.625, 0.625, 1],
                                  [0.875, 0.875, 1],
                                  [1.125, 1.125, 1],
                                  [1.375, 1.375, 1],
                                  [1.625, 1.625, 1],
                                  [1.875, 1.875, 1],
                                  [2.125, 2.125, 1],
                                  [2.375, 2.375, 1],
                                  [2.625, 2.625, 1],
                                  [2.875, 2.875, 1],
                                  [3.125, 3.125, 1],
                                  [3.375, 3.375, 1],
                                  [3.625, 3.625, 1],
                                  [3.875, 3.875, 1], ] ),
                           'equal_50cm':
                           array( [[0.25 , 0.25, 1],
                                 [0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1],
                                 [0.75, 0.75 , 2],
                                 [1.25 , 1.25, 2],
                                 [1.75 , 1.75, 2],
                                 [2.25, 2.25, 2],
                                 [2.75, 2.75, 2],
                                 [3.25, 3.25, 2],
                                 [3.75, 3.75, 2]] ),
                           'equal_100cm':
                           array( [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1],
                                 [0.5, 0.5 , 2],
                                 [1.5 , 1.5, 4],
                                 [2.5, 2.5, 4],
                                 [3.5, 3.5, 4],
                                 ] ),
                           'equal_200cm':
                           array( [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1],
                                 [1.0, 1.0 , 7],
                                 [3.0, 3.0, 7],
                                 ] ),
                           'gap_100cm':
                           array( [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1],
                                 [1.0, 1.0 , 5],
                                 [2.0, 2.0 , 5],
                                 [3.0, 3.0, 5], ] ),
                           'corner_2to4_25cm':
                           array( [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1],
                                  [2.125, 2.125, 5],
                                  [2.375, 2.375, 1],
                                  [2.625, 2.625, 1],
                                  [2.875, 2.875, 1],
                                  [3.125, 3.125, 1],
                                  [3.375, 3.375, 1],
                                  [3.625, 3.625, 1],
                                  [3.875, 3.875, 1], ] ),
                          'middle_0to3_25cm':
                           array( [[0.125 , 0.125, 1],
                                  [0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1],
                                  [0.375, 0.375, 1],
                                  [0.625, 0.625, 1],
                                  [0.875, 0.875, 1],
                                  [1.125, 1.125, 1],
                                  [1.375, 1.375, 1],
                                  [1.625, 1.625, 1],
                                  [1.875, 1.875, 1],
                                  [2.125, 2.125, 1],
                                  [2.375, 2.375, 1],
                                  [2.625, 2.625, 1],
                                  [2.875, 2.875, 1], ] ),
                           'middle_0to2_25cm':
                           array( [[0.125 , 0.125, 1],
                                  [0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1],
                                  [0.375, 0.375, 1],
                                  [0.625, 0.625, 1],
                                  [0.875, 0.875, 1],
                                  [1.125, 1.125, 1],
                                  [1.375, 1.375, 1],
                                  [1.625, 1.625, 1],
                                  [1.875, 1.875, 1], ] )
                           }

        n_elems_quarter_xy_dict = {'equal_25cm': 19,
                                   'equal_50cm': 18,
                                   'equal_100cm': 17,
                                   'equal_200cm': 17,
                                   'gap_100cm': 18,
                                   'corner_2to4_25cm': 15,
                                   'middle_0to3_25cm':18,
                                   'middle_0to2_25cm':18}

        not_linked_dict = {'equal_25cm': 1,
                           'equal_50cm': 1,
                           'equal_100cm': 0,
                           'equal_200cm': 0,
                           'gap_100cm': 0,
                           'corner_2to4_25cm': 0,
                           'middle_0to3_25cm':1,
                           'middle_0to2_25cm':1}

        for link_type in link_type_list:

            for link_case in link_case_list:
                del sim_model
                sim_model = MRtwo()
                sim_model.select_hinge_pos = True


                for lc in lc_list:

                    sim_model.link_type = link_type
                    sim_model.n_elems_quarter_xy = n_elems_quarter_xy_dict[link_case]
                    sim_model.shift_array = shift_array_dict[link_case]
                    sim_model.not_linked_elem = not_linked_dict[link_case]
                    sim_model.lc = lc
                    print  'link_type', sim_model.link_type
                    print 'loading_case:', sim_model.lc
                    print "n_elems_quarter_xy", sim_model.n_elems_quarter_xy
                    sim_model.peval()
                    sim_model.export_int_force_data( filename = 'spring_1_' + link_type + '_' + link_case + '_hf_' + lc + '.csv' )
                    sim_model.export_edge_u_data( filename = 'spring_1_' + link_type + '_' + link_case + '_u_' + lc + '.csv' )
#                    sim_model.F_int_bar_plot( filename = 'spring_1_' + link_type + '_' + link_case + '_plt_' + lc + '_' )


#        
#        for lc in lc_list:
#            for         


    elif do == 'ps':

        sim_ps = SimPStudy( sim_model = sim_model )
        sim_ps.configure_traits()

    elif do == 'pickle':

        import pickle
        filename = '/tmp/sim.pickle'
        file = open( filename, 'w' )
        pickle.dump( sim_model, file )
        file.close()
        file = open( filename, 'r' )
        sm = pickle.load( file )
        file.close()



    elif do == 'cs':

#        sim_ps = SimPStudy( sim_model = sim_model )
        pstudy = SimArray( sim_model = sim_model )
        pstudy_view = SimArrayView( model = pstudy )
#        pstudy_view._start_study()
        pstudy_view.configure_traits()


        factor_list = pstudy.factor_list


        n_levels_list = [ factor.get_n_levels() for factor in factor_list ]
#        print 'get array0', pstudy.output_array.shape
#        print 'get array1', pstudy.output_array
        from matplotlib.pyplot import *
        output_array = pstudy.output_array



        cs = 'column'
        choosen_x , choosen_z = 10, 1


        if cs == 'shell':
            print output_array.shape
            y_0 = output_array[:, 0, :]
            y_1 = output_array[:, 1, :]
#            y_2 = output_array[:,2,:]
        if cs == 'column':

            y_0 = output_array
#            y_1 = output_array[1,:,:]
#            y_2 = output_array[2,:,:]


        # shape of input [kombinations , number of ps_level]
        #   
        input_array = pstudy.input_table

        #choosen output
        #

#        cs = Enum('column','shell') 


        idx_x_axes = Int
        idx_z_axes = Int

        for i, name in enumerate( pstudy.factor_names ):
            if name == 'n_elems_col_xy':
               idx_z_axes = i
            if name == 'n_elems_col_z':
               idx_x_axes = i
            if name == 'n_elems_xy_quarter':
               idx_x_axes = i
            if name == 'n_elems_z':
               idx_z_axes = i


        def n_dof_serendipity_grid( x, y, z ):
            return ( ( 2 * x + 1 ) * ( 2 * y + 1 ) * ( 2 * z + 1 ) - 6 * x * y * z \
                   + ( ( x - 1 ) * y * z + ( y - 1 ) * x * z + ( z - 1 ) * x * y ) - x * y * z ) * 3.0

        if cs == 'shell':
            n_dofs = n_dof_serendipity_grid( input_array[:, 0] * 2,
                                            input_array[:, 0] * 2,
                                            input_array[:, 1] )
            n_dofs_0 = n_dofs[where( input_array[:, idx_z_axes] == 1 )]
            n_dofs_1 = n_dofs[where( input_array[:, idx_z_axes] == 2 )]
            n_dofs_2 = n_dofs[where( input_array[:, idx_z_axes] == 3 )]

        if cs == 'column':
            n_dofs = n_dof_serendipity_grid( ones_like( input_array ) * 2,
                                            ones_like( input_array ) * 2,
                                            input_array )
            n_dofs_0 = n_dofs
#            n_dofs_0 = n_dofs[where(input_array[:,idx_z_axes]==2)]
#            n_dofs_1 = n_dofs[where(input_array[:,idx_z_axes]==4)]
#            n_dofs_2 = n_dofs[where(input_array[:,idx_z_axes]==4)]





        idx_1 = where( input_array[:, idx_x_axes] == choosen_x )
#        idx_2 = where(input_array[:,idx_z_axes]== choosen_z)

        for i in  idx_1[0]:
#            for j in idx_2[0]:
#            if i==j:
            idx_choosen = i

#        for i in  idx_1[0]:
#            idx_choosen = i

        n_dofs_choosen = n_dofs[idx_choosen]

        y_choosen = c_[y_0,
#                       y_1,
#                       y_2
                       ][where( c_[n_dofs_0,
#                       n_dofs_1,
#                       n_dofs_2
                       ] == n_dofs_choosen )]

        fig = figure( facecolor = 'white' )
        ax1 = fig.add_subplot( 1, 1, 1 )

#        from matplotlib import rc
#        rc( 'text', usetex = True )
#        rc( 'font', **{'family':'serif', 'serif':'Times'} )
#        rc( 'text.latex',  )
#        ax1.set_xticks((50000,100000,300000))
#        ax1.set_xticklabels(("50000","100000", '300000'))

        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize( 20 )

        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize( 20 )



        if cs == 'shell':

            ax1.plot( n_dofs_0 , y_0, color = 'b', label = '1', linewidth = 1.5 )
            ax1.plot( n_dofs_1 , y_1, color = 'g', label = '2', linewidth = 1.5 )
#            ax1.plot(n_dofs_2 ,y_2, color = 'r', label = '3', linewidth=1.5)
            ax1.plot( n_dofs_choosen.reshape( -1, 1 ), y_choosen.reshape( -1, 1 ), marker = 'o', color = 'b', markerfacecolor = 'b', markersize = 8 )#
            xlim ( 0, 15000 )

        if cs == 'column':
            xlim ( 0, 2000 )
            ax1.plot( n_dofs_0 , y_0, color = 'b', label = '2', linewidth = 1.5 )
#            ax1.plot(n_dofs_1 ,y_1, color = 'g', label = '4', linewidth=1.5)
#            ax1.plot(n_dofs_2 ,y_2, color = 'r', label = '4', linewidth=1.5)
            ax1.plot( n_dofs_choosen.reshape( -1, 1 ), y_choosen.reshape( -1, 1 ), marker = 'o', color = 'b', markerfacecolor = 'b', markersize = 8 )#

        ylim ( -0.0040, -0.0038 )

#        ax1.plot( [0.25,0.75],[0.45/2**0.5,2], 'bo')
        ax1.set_xlabel( 'Freiheitsgrade ', fontsize = 24 )
        ax1.set_ylabel( 'U [m]', fontsize = 24 )

        legend()
        show()

        filename = 'cs_' + cs + '.csv'
#        
        print '*** writing study data to file,', filename, ' ***'

        if cs == 'column':
            X_data = vstack( ( hstack( ( n_dofs_0.reshape( -1, 1 ), ones_like( n_dofs_0.reshape( -1, 1 ) ) * 2, y_0 ) ) ) )
        if cs == 'shell':
            X_data = vstack( ( hstack( ( n_dofs_0.reshape( -1, 1 ), ones_like( n_dofs_0.reshape( -1, 1 ) ) * 1, y_0 ) ),
                            hstack( ( n_dofs_1.reshape( -1, 1 ), ones_like( n_dofs_1.reshape( -1, 1 ) ) * 2, y_1 ) ) ) )



        file = open( filename, 'w' )
#        
        writer = csv.writer( file, delimiter = ";", lineterminator = "\n" )
        writer.writerow( ['nodes' + cs , 'n_elems_t', 'U_z'] )
        writer.writerows( X_data )

        file = file.close()





#        pyplot.plot
#        sim_array.clear_cache()
#        print 'array_content', sim_array.output_array
#        sim_array.configure_traits()

#        sim_ps.configure_traits()

