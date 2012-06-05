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
from ibvpy.mats.mats3D.mats3D_tensor import \
    map3d_sig_eng_to_mtx

from ibvpy.fets.fets_eval import \
    FETSEval
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

from ibvpy.mesh.fe_spring_array import \
     FESpringArray

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, shape, \
    cos, sin, arctan, where, abs, all, any, diag, \
    argsort, sum, transpose, vstack, append, arange, \
    hstack, reshape, shape, size, zeros_like, cumsum, \
    dsplit, copy, dtype, sort, ones_like, max, argmax

# Interpolation
from scipy.interpolate import Rbf

from math import \
    sqrt, asin, acos, pi as Pi

from matplotlib.pyplot import \
    bar, show, axhline, ion, ioff, xlabel, ylabel, title, figure, savefig

import csv


from simiter.sim_pstudy import\
    ISimModel, SimOut, SimPStudy, SimArray, SimArrayView


from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

from geo_column import GEOColumn

from hp_shell import HPShell

from mush_roof_model import MushRoofModel


class MRtwo( MushRoofModel ):

    implements( ISimModel )

    # parameter used by hp_shell
    # if set to "one" hp shell uses symmetry of one roof to construct geometry
    # otherwise only one quarter is assumed.
    #
    mushroof_part = 'one'

    # shift of elements in hp_shell
    #
    shift_elems = True

    #----------------------------------------------------
    # elements
    #----------------------------------------------------

    n_elems_xy_quarter = Int( 10, input = True )#, ps_levels = [4, 16, 5] )
    n_elems_z = Int( 1, input = True )#, ps_levels = [1, 2, 1] )

    n_elems_col_z = Int( 10 , input = True )#, ps_levels = [5, 20, 3 ] )
    n_elems_col_xy = Int( 2 , input = True )#, ps_levels = [2,4,1]   )

    # params are defined in the base class "MushRoofModel"
    # uncomment if pstudy is to be run:
    #
#    delta_h_scalefactor = Float( 1.0 , ps_levels = [1.0, 1.3, 3] )
#    width_bottom_col = Float( .35, ps_levels = [0.35, 0.6, 6] )


    vtk_r = Float( 0.9, input = True )

    # fets used for roof
    #
    fe_roof = Instance( ( FETSEval ), depends_on = '+ps_levels, +input' )
    def _fe_roof_default( self ):
        # fe_quad_serendipity_roof is defined in base class 
        # connected with material properties of the roof
        #
        fets = self.fe_quad_serendipity_roof
        fets.vtk_r *= self.vtk_r
        return fets

    # fets used for steel plate
    #
    fe_plate = Instance( ( FETSEval ), depends_on = '+ps_levels, +input' )
    def _fe_plate_default ( self ):
        fets = self.fe_quad_serendipity_plate
        # default integration scheme defined as 2x2x2
        # use higher order integration needed to get better 
        # results at connection points with the column
        #
        fets.ngp_r = 3
        fets.ngp_s = 3
        fets.ngp_t = 3
        fets.vtk_r *= self.vtk_r
        return fets

    # fets used for columns
    #
    fe_column = Instance( ( FETSEval ), transient = True , depends_on = '+ps_levels, +input' )
    def _fe_column_default( self ):
        # default integration (2x2x2) sufficient in case of 2x2 discretization
        # (NOTE: for 1x1 discetrization hourglassing occured)
        #
        fets = self.fe_quad_serendipity_column
        fets.vtk_r *= self.vtk_r
        return fets

    #----------------------------------------------------
    # geometric transformation
    #----------------------------------------------------

    hp_shells = Property( depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_hp_shells( self ):
        # list of local origins of each roof defined in global coordinates
        # (for MRtwo two roofs are considered due to symmetry)
        #
        X_list = [ [0, 0, 0], [8.0, 0, 0] ]
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
        # translation to move comumns in roof center
        #
        X_list = [[ 4., 4., -ch ],
                  [ 12, 4., -ch ], ]
        return [ GEOColumn( width_top = self.width_top_col,
                            X0 = X,
                            width_bottom = self.width_bottom_col,
                            h_col = ch - self.t_plate )
                            for X in X_list ]

    plates = Property( Instance( GEOColumn ) , depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_plates( self ):
        # column transformation can be used to obtain the geometry of the plate
        #
        ch = self.t_plate
        X_list = [[ 4., 4., -ch ],
                  [ 12, 4., -ch ], ]
        return [GEOColumn( width_top = self.width_top_col,
                           width_bottom = self.width_top_col,
                           X0 = X,
                           h_col = self.t_plate )
                           for X in X_list ]
    #----------------------------------------------------
    # fe-grids 
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

    #----------------------------------------------------
    # ps_study
    #----------------------------------------------------
    def peval( self ):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''

        #----------------------------------------------------
        # SORTING
        #----------------------------------------------------

        # the data will be evaluated and for the output the dofs need to be sorted
        # numbering of dofs within one element corresponds to numbering in Belytschko...

        # @TODO: 
        #        - one array with coordinates and dof number for all necessary slices
        #          therefore no sorting necessary

        # sorted dofs, shape remains unchanged
        #
        dofs_r0_right = self.sort_by_dofs( self.dofs_r0_right, self.dofs_r0_right )
        dofs_r0_left = self.sort_by_dofs( self.dofs_r0_left, self.dofs_r0_left )
        dofs_r1_right = self.sort_by_dofs( self.dofs_r1_right, self.dofs_r1_right )
        dofs_r1_left = self.sort_by_dofs( self.dofs_r1_left, self.dofs_r1_left )
        dofs_bottom = self.sort_by_dofs( self.dofs_bottom, self.dofs_bottom )
        dofs_sym = self.sort_by_dofs( self.dofs_sym, self.dofs_sym )

        # sorted coordinates, shape remains unchanged
        #
        dof_X_r0_right = self.sort_by_dofs( self.dofs_r0_right, self.dof_X_r0_right )
        dof_X_r0_left = self.sort_by_dofs( self.dofs_r0_left, self.dof_X_r0_left )
        dof_X_r1_right = self.sort_by_dofs( self.dofs_r1_right, self.dof_X_r1_right )
        dof_X_r1_left = self.sort_by_dofs( self.dofs_r1_left, self.dof_X_r1_left )
        dof_X_bottom = self.sort_by_dofs( self.dofs_bottom, self.dof_X_bottom )
        dof_X_sym = self.sort_by_dofs( self.dofs_sym, self.dof_X_sym )


        # corner nodes are doubled therefore undo slices
        #
        dofs_r0_right = append( dofs_r0_right[:, :-1, :],
                                dofs_r0_right[-1, -1, :] ).reshape( -1, 3 )
        dofs_r0_left = append( dofs_r0_left[:, :-1, :],
                                dofs_r0_left[-1, -1, :] ).reshape( -1, 3 )
        dofs_bottom = append( dofs_bottom[:, :-1, :],
                           dofs_bottom[-1, -1, :] ).reshape( -1, 3 )
        dofs_r1_right = append( dofs_r1_right[:, :-1, :],
                                dofs_r1_right[-1, -1, :] ).reshape( -1, 3 )
        dofs_r1_left = append( dofs_r1_left[:, :-1, :],
                                dofs_r1_left[-1, -1, :] ).reshape( -1, 3 )
        dofs_sym = append( dofs_sym[:, :-1, :],
                           dofs_sym[-1, -1, :] ).reshape( -1, 3 )


        dof_X_r0_right = append( dof_X_r0_right[:, :-1, :],
                                dof_X_r0_right[-1, -1, :] ).reshape( -1, 3 )
        dof_X_r0_left = append( dof_X_r0_left[:, :-1, :],
                                dof_X_r0_left[-1, -1, :] ).reshape( -1, 3 )
        dof_X_bottom = append( dof_X_bottom[:, :-1, :],
                           dof_X_bottom[-1, -1, :] ).reshape( -1, 3 )
        dof_X_r1_right = append( dof_X_r1_right[:, :-1, :],
                                dof_X_r1_right[-1, -1, :] ).reshape( -1, 3 )
        dof_X_r1_left = append( dof_X_r1_left[:, :-1, :],
                                dof_X_r1_left[-1, -1, :] ).reshape( -1, 3 )
        dof_X_sym = append( dof_X_sym[:, :-1, :],
                           dof_X_sym[-1, -1, :] ).reshape( -1, 3 )

        #----------------------------------------------------
        # EVALUATE DATA
        #----------------------------------------------------


        # DISPLACEMENT
        #
        U = self.tloop.eval()

        # slice r0
        # 
        U_r0_right = U[dofs_r0_right]
        U_r0_left = U[dofs_r0_left]

        # slice r0
        #
        U_r1_right = U[dofs_r1_right]
        U_r1_left = U[dofs_r1_left]

        # slice over both domains
        #
        U_bottom = U[dofs_bottom]
        U_sym = U[dofs_sym]

        # INTERNAL FORCES
        #
        F_int = self.tloop.tstepper.F_int

        # slice r0
        # 
        F_int_r0_right = F_int[dofs_r0_right]
        F_int_r0_left = F_int[dofs_r0_left]

        # slice r1
        #
        F_int_r1_right = F_int[dofs_r1_right]
        F_int_r1_left = F_int[dofs_r1_left]

        # slice over both domains
        #
        F_int_bottom = F_int[dofs_bottom]
        F_int_sym = F_int[dofs_sym]




        #----------------------------------------------------
        # EXPORT DATA
        #----------------------------------------------------


        # export of FORCE and corresponding COORDINATES 

        # coordinate of internal forces at all nodes at slices 
        #
        X_F_export = vstack( ( dof_X_r0_right,
                                  dof_X_sym ) )

        # export only linked dofs 
        # position is defined by X-array_linked and Y_array_linked
        #
        if self.select_hinge_pos == True:
            # if selected linking is choosen export only selected nodes
            #
            idx = []
            # X value corresponds to X array linked
            #
            for X in self.X_array_linked:
                # due to round of errors small epsilon introduced to find the necessary values
                #
                idx = idx + list( where( abs( X_F_export[:, 0] - X ) <= 0.00001 )[0] )

            # Y value corresponds to Y array linked
            #
            for Y in self.Y_array_linked:
                # due to round of errors small epsilon introduced to find the necessary values
                #
                idx = idx + list( where( abs( X_F_export[:, 1] - Y ) <= 0.00001 )[0] )

            X_F_export = X_F_export[idx, :]

        # coordinate of internal forces for export
        #
        self.X_F_export = X_F_export


        # sort F for export with in the following order
        # [:,:,0] = N_ip
        # [:,:,1] = V_ip
        # [:,:,2] = V_op
        # order of internal forces
        # [:,:,0] = F_x
        # [:,:,1] = F_y
        # [:,:,2] = F_z
        #
        F_int_r0_right_export = copy( F_int_r0_right )
        F_sym_export = copy( F_int_sym )
        F_sym_export[:, 1] = F_int_sym[:, 0]
        F_sym_export[:, 0] = F_int_sym[:, 1]

        # order needs to be identical with X_F_export
        #
        F_export = vstack( ( F_int_r0_right_export, F_sym_export ) )
        if self.select_hinge_pos == True:
            # if selected linking is choosen export only selected nodes
            # same index as before
            #
            F_export = F_export[idx, :]

        # forces of internal forces for export
        #
        self.F_export = F_export


        # export of DISPLACEMENT and corresponding COORDINATES 


        # coordinates at all nodes at slices, here all edges of each roof
        #
        self.X_U_export = vstack( ( dof_X_sym,
                                  dof_X_r0_left,
                                  dof_X_bottom,
                                  dof_X_r1_right,
                                  dof_X_r0_right ) )

        # delta values for displacement (VERSATZ)
        #
        dU_sym = zeros_like( U_sym[:, 1:] )
        dU_r0_left = U_r0_left[:, 1:] - U_r0_left[:, 1:]  # dUx can be neglected
        dU_bottom = zeros_like( U_bottom[:, 1:] )
        dU_r1_right = zeros_like( U_r1_right[:, 1:] )
        dU_r0_right = U_r0_right[:, 1:] - U_r1_left[:, 1:] # dUx can be neglected

        # some bad DATA arrangement to get all together 
        #
        data_U_sym = hstack( ( U_sym, dU_sym ) )
        data_U_r0_left = hstack( ( U_r0_left, dU_r0_left ) )
        data_U_bottom = hstack( ( U_bottom, dU_bottom ) )
        data_U_r1_right = hstack( ( U_r1_right, dU_r1_right ) )
        data_U_r0_right = hstack( ( U_r0_right, dU_r0_right ) )

        # displacement at all nodes at slices, here all edges of each roof
        #    
        self.U_export = vstack( ( data_U_sym,
                                data_U_r0_left,
                                data_U_bottom,
                                data_U_r1_right,
                                data_U_r0_right, ) )

        U_r0_right = U_bottom[0, 2]

        return array( [ U_bottom[argmax( abs( U_bottom[:, 2].reshape( -1 ) ) )][2]] )


    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [  SimOut( name = '$u_z$', unit = '[mm]' ), ]


    X_F_export = Array  # X, Y, Z
    F_export = Array    # Nip, Vip, Vop order defined in peval


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


    X_U_export = Array  # X, Y, Z
    U_export = Array    # 'u-x[m]', 'u-y[m]', 'u-z[m]', 'du-y[m]', 'du-z[m]'


    def export_edge_u_data( self, filename = 'U_data.csv' ):
        '''exports X_U_export and U_export data to csv - worksheet
        '''
        print '*** writing displacement data to file,', filename, ' ***'

        X_data = self.X_U_export.reshape( -1, 3 )
        U_data = self.U_export.reshape( -1, 5 )

        data = hstack( ( X_data, U_data ) )
        file = open( filename, 'w' )

        writer = csv.writer( file, delimiter = ";", lineterminator = "\n" )
        writer.writerow( ['X[m]', 'Y[m]', 'Z[m]', 'u-x[m]', 'u-y[m]', 'u-z[m]', 'du-y[m]', 'du-z[m]' ] )
        writer.writerows( data )

        file = file.close()
        return


    X_array_linked = Property()
    def _get_X_array_linked( self ):
        """gets the X coordinates of all linked elements depending on shift_array
        """
        # @TODO: 
        # make linked_array a new method of MRTwo used in multiple methods
        #
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

        # @TODO: 
        # make linked_array a new method of MRTwo used in multiple methods
        #
        linked_array = vstack( ( self.shift_array[:self.not_linked_elem],
                               self.shift_array[self.not_linked_elem + 1:] ) )

        #global coordinates from local coordinates of one quarter
        #
        y_array_linked = hstack( ( 4 - linked_array[:, 1],
                                   4 + linked_array[:, 1], ) )
        return sort( y_array_linked )


    # internal Forces plot no more need for that, plot now done in lcc_tbale_hf
    #
#    def F_int_bar_plot( self, filename ):
#        """ creates a bar plot of the internal forces and saves them
#        """
#        from matplotlib import pyplot
#        print '*** plot of hinge force data data to plot file,', filename, ' ***'
#
#        X_data = self.X_F_export.reshape( -1, 3 )        # 1st x 2nd y 3rd z
#        F_data = self.F_export.reshape( -1, 3 ) * 1000    # 1st N_ip 2nd V_ip 3rd V_op
#
#        def save_bar_plot( x, y, filename = 'bla', title = 'Title',
#                          xlabel = 'xlabel', ylabel = 'ylavel',
#                          width = 0.1, xmin = 0, xmax = 1000 ):
#            fig = pyplot.figure( facecolor = "white" )
#            ax1 = fig.add_subplot( 1, 1, 1 )
#            ax1.bar( x , y , width = width, align = 'center', color = 'blue' )
#            ax1.set_xlim( xmin, xmax )
#            ax1.set_xlabel( xlabel, fontsize = 22 )
#            ax1.set_ylabel( ylabel, fontsize = 22 )
#            fig.savefig( filename, orientation = 'portrait', bbox_inches = 'tight' )
#            pyplot.clf()
#
#        # symmetric axes
#        #
#        idx_sym = where( abs( X_data[:, 1] - 8 ) <= 0.0001 )
#        X_sym = X_data[idx_sym, 0].reshape( -1 )
#        F_sym = F_data[idx_sym]
#
#        save_bar_plot( X_sym, F_sym[:, 0],
#                      xlabel = '$X$ [m]', ylabel = '$N_{ip}$ [kN]',
#                      filename = filename + 'N_ip' + '_sym',
#                      xmin = 0.0, xmax = 16.0 )
#        save_bar_plot( X_sym, F_sym[:, 1],
#                      xlabel = '$X$ [m]', ylabel = '$V_{ip}$ [kN]',
#                      filename = filename + 'V_ip' + '_sym',
#                      xmin = 0.0, xmax = 16.0 )
#        save_bar_plot( X_sym, F_sym[:, 2],
#                      xlabel = '$X$ [m]', ylabel = '$V_{op}$ [kN]',
#                      filename = filename + 'V_op' + '_sym',
#                      xmin = 0.0, xmax = 16.0 )
#
#        # r0_r1
#        #
#        idx_r0_r1 = where( abs( X_data[:, 0] - 8 ) <= 0.0001 )
#        X_r0_r1 = X_data[idx_r0_r1, 1].reshape( -1 )
#        F_r0_r1 = F_data[idx_r0_r1]
#
#        save_bar_plot( X_r0_r1, F_r0_r1[:, 0],
#                      xlabel = '$X$ [m]', ylabel = '$N_{ip}$ [kN]',
#                      filename = filename + '_' + 'N_ip' + '_r0_r1',
#                      xmin = 0.0, xmax = 8.0 )
#        save_bar_plot( X_r0_r1, F_r0_r1[:, 1],
#                      xlabel = '$X$ [m]', ylabel = '$V_{ip}$ [kN]',
#                      filename = filename + '_' + 'V_ip' + '_r0_r1',
#                      xmin = 0.0, xmax = 8.0 )
#        save_bar_plot( X_r0_r1, F_r0_r1[:, 2],
#                      xlabel = '$X$ [m]', ylabel = '$V_{op}$ [kN]',
#                      filename = filename + '_' + 'V_op' + '_r0_r1',
#                      xmin = 0.0, xmax = 8.0 )


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

#    bc_col_link_list = Property( List, depends_on = '+ps_levels, +input' )
#    @cached_property
#    def _get_bc_col_link_list( self ):
#        '''links column top to roof
#        ( all column corner nodes of each elements to the adjacent elements of the roof )
#        '''
#        bc_col_link_list = []
#
#        n_from_center = 1 + self.not_linked_elem
#
#
#        for roof, column in zip( self.fe_grid_roofs, self.fe_grid_columns ):
#            slice_1 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
#                                  slice = roof[self.n_elems_xy_quarter - n_from_center ,
#                                               self.n_elems_xy_quarter, 0,
#                                               0, 0, 0 ],
#                                  link_slice = column[ 0 , 0 , -1, 0, 0, -1], link_coeffs = [1.0],
#                                  value = 0. )]
#            slice_2 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
#                                  slice = roof[self.n_elems_xy_quarter,
#                                               self.n_elems_xy_quarter - n_from_center, 0,
#                                               0, 0, 0 ],
#                                  link_slice = column[ -1, 0, -1, -1, 0, -1], link_coeffs = [1.0],
#                                  value = 0. )]
#
#            slice_3 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
#                                  slice = roof[self.n_elems_xy_quarter + n_from_center,
#                                               self.n_elems_xy_quarter, 0,
#                                               0, 0, 0 ],
#                                  link_slice = column[ -1 , -1 , -1, -1, -1, -1], link_coeffs = [1.0],
#                                  value = 0. )]
#
#            slice_4 = [BCSlice( var = 'u'  , dims = [0, 1, 2],
#                                  slice = roof[self.n_elems_xy_quarter ,
#                                               self.n_elems_xy_quarter + n_from_center, 0,
#                                               0, 0, 0 ],
#                                  link_slice = column[ 0 , -1 , -1, 0, -1, -1], link_coeffs = [1.0],
#                                  value = 0. )]
#        return bc_col_link_list


    # column support as hinge at bottom of column
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
    # arbitrary linking of roof domains
    #-------------------------------------

    # array of shape (n,3) used for shifting of shell elements
    # example shift_array = array([2.0, 2.0, 3]) 
    # use one row for only two sections within the quarter with varying element finess
    # (varying discretization).
    # first value in row of shift array := x-coordinate (e.g. 2m) of the quarter shell
    # 2nd value in row of  shift array := y-coordinate  (e.g. 2m) of the quarter shell
    # 3rd value number of elements to be used for this section of the quarter shell
    # (e.g. 3 elements) between 0 and 2.0m.
    #
    # NOTE: the number of elements between 2.0m and 4.0m is derived based on 
    # the total number of elements used for one quarter (n_elems_quarter_xy)
    # minus the number of elements used for the previous sections. 
    # Make sure to define the sum of all 3rd values in shift array smaller
    # then (n_elems_quarter_xy -1)
    #
    shift_array = Array( value = [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1], ], input = True )

    # link_type of roof connections (links)
    # in-plane degree of freedom between the roofs either linked or not
    #
    link_type = Enum( 'exc_V_ip', 'inc_V_ip' )
    linked_dims = {'exc_V_ip':[0, 2],
                   'inc_V_ip':[0, 1, 2]}

    # Set to true if the links between the roof are defined by shift_array 
    # only at the defined set of nodes
    #
    # Set to "false" if all nodes between the roofs are to be linked at edge  
    #
    select_hinge_pos = Bool( True )

    # defines element of shift_array which is not linked 
    # predefined in dictionary of shift_array_dict - either 0 or 1
    # If set to one the node used to link the roof and the column is the first 
    # node defined in shift array 
    #
    not_linked_elem = Int( 0 ) #not linked element depends on shift array     


    link_elem_array = Property( Array )
    def _get_link_elem_array( self ):
        """ from shift_array get element number which have to be linked
            shift_array is here by only defined for one quarter
        """
        # abreviation "e" = elements
        e_shifted = cumsum( self.shift_array[:, 2] )
        # create link list of elements 

        # NOTE: the definition depends directly on the definition of BCSlice
        # used for linking. 
        # In BCSlice the following syntax is applied:
        # [:,:,:,0,0,0] not [:,:,:,-1,-1,-1] will be linked
        # without unselected node defined by "not_linked_elem" and the
        # corresponding coordinate in "shift_array" which is only used for link
        # of column with the roof
        #

        # not_linked element neglected
        #
        # @TODO: replace with (varify):
        # e_shifted = delete(e_shifted,[self.not_linked_elem])
        #
        e_linked = append( e_shifted[:self.not_linked_elem],
                           e_shifted[self.not_linked_elem + 1:] )

        # for linking the global coordinate system is used 
        # respectively the coordinate system used for the slicing
        # within the BCSlice
        #
        # element list for the whole roof instead of one quarter
        #
        e_linked_roof = append( self.n_elems_xy_quarter - e_linked,
                                self.n_elems_xy_quarter + e_linked )
        return e_linked_roof

    bc_roof_link_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_bc_roof_link_list( self ):
        '''link roof to roof
        and boundary conditions for the symmetry assumption
        '''
        roof_0, roof_1 = self.fe_grid_roofs
        if self.select_hinge_pos == False:

            # constraint at symmetric axes for roof 0
            #
            constraint_0 = [BCSlice( var = 'u', dims = [1],
                                     # get only every corner node
                                     slice = roof_0[ : , -1, -1, 0, -1, -1],
                                     value = 0.0 )]

            # constraint at symmetric axes for roof 1
            #
            constraint_1 = [BCSlice( var = 'u', dims = [1],
                                     slice = roof_1[ : , -1, -1, -1, -1, -1],
                                     value = 0.0 )]

            # constraint for center node connecting roof 0 and 1 at synmmmetric axes
            #
            constraint_0_1 = [BCSlice( var = 'u', dims = [1],
                                     slice = roof_0[-1 , -1, -1, -1, -1, -1],
                                     value = 0.0 )]

            # link of roof 0 and 1
            #
            link_0_1 = [BCSlice( var = 'u'  , dims = self.linked_dims[self.link_type],
                                         slice = roof_0[-1, :, -1, -1, 0, -1 ],
                                         link_slice = roof_1[ 0, :, -1, 0, 0, -1],
                                         link_coeffs = [1.0],
                                         value = 0. )]

            # center point of all roofs at symmetry axis and at 01-axis
            # is only linked between the roofs r0 and r1
            #
            link_0_1_last = [BCSlice( var = 'u'  , dims = self.linked_dims[self.link_type],
                                         slice = roof_0[-1, -1, -1, -1, -1, -1 ],
                                         link_slice = roof_1[ 0, -1, -1, 0, -1, -1],
                                         link_coeffs = [1.0],
                                         value = 0. )]

            bc_roof_link_list = constraint_0 + constraint_1 + constraint_0_1 + \
                                link_0_1 + link_0_1_last

        # if a predefined positioning for the roof linking is to be used:
        #
        if self.select_hinge_pos == True:
            constraint_list = []
            for elem in self.link_elem_array:
                # link between roofs r0 and r1
                #     
                link_0_1 = [BCSlice( var = 'u'  , dims = self.linked_dims[self.link_type],
                                         # slice selects only cornernodes 
                                         slice = roof_0[-1, elem, -1, -1, 0, -1 ],
                                         link_slice = roof_1[ 0, elem, -1, 0, 0, -1],
                                         link_coeffs = [1.0],
                                         value = 0. )]
                # symmetric constraint (symmetry axis)
                #
                const_0 = [BCSlice( var = 'u'  , dims = [1],
                                         # slice selects only cornernodes 
                                         slice = roof_0[elem, -1, -1, 0, -1, -1 ],
                                         value = 0. )]
                const_1 = [BCSlice( var = 'u'  , dims = [1],
                                         # slice selects only cornernodes 
                                         slice = roof_1[elem, -1, -1, 0, -1, -1 ],
                                         value = 0. )]

                constraint_list = constraint_list + link_0_1 + const_0 + const_1

            bc_roof_link_list = constraint_list

        return bc_roof_link_list


    #-------------------------------------
    # loading cases
    #-------------------------------------
    # @TODO: - lc_shrink not defined just empty list only possible by Bool in mushroof_model
    #          so far would be better to be done automatically, name necessary for export 
    #

    # load cases
    #
    lc = Enum( 'lc_g',
               'lc_s_sym', 'lc_s_asym',
               'lc_w_pos', 'lc_w_neg', 'lc_w_asym',
               'lc_shrink', input = True )

    # Dictionary for corresponding methods
    #
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

#    
#    lc_shrinkage = Property( List, input=True, depends_on = '+ps_levels, +input' )
#    @cached_property
#    def _get_lc_shrinkage(self):
#        self.mats_roof = MATS3DElastic( E = self.E_roof, 
#                                        nu = self.nu, 
#                                        initial_strain = self.temperature_strain_z )
#        
#        self.mats_column = MATS3DElastic( E = self.E_column, 
#                                          nu = self.nu, 
#                                          initial_strain = self.temperature_strain_Z )



    #----------------------------------------------------
    # response tracer
    #----------------------------------------------------

    rtrace_list = List
    def _rtrace_list_default( self ):
        return [  self.max_princ_stress, self.sig_app, self.u]

    #----------------------------------------------------
    # time loop
    #----------------------------------------------------

    tloop = Property( depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_tloop( self ):

        roofs = self.fe_grid_roofs
        columns = self.fe_grid_columns
        plates = self.fe_grid_plates

        ts = TS( sdomain = roofs + columns + plates,
                 dof_resultants = True,
                 bcond_list =

                            # boundary conditions of the selected loading case
                            # e.g "lc_g_list", "lc_s_sym_list", etc.
                              self.lc_list +



                            # links (kinematic constraints and symmetric assumptions)
                              self.bc_roof_link_list +
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

        # output slices
        #
        r0, r1 = roofs


        #output slices 
        # @todo: did not manage to create external methodwould be 
        #        better for input problems and long method loading before evaluation
        # 
        #slices r0
        #
        self.dofs_r0_right = r0[-1, :, -1, -1, :, -1].dofs
        self.dofs_r0_left = r0[0, :, -1, 0, :, -1].dofs
        self.dofs_r0_top = r0[:, -1, -1, :, -1, -1].dofs
        self.dofs_r0_bottom = r0[:, 0, -1, :, 0, -1].dofs

        #get all the dofs that are subjedcted to be linked
        #
        #slices r0 global coordinates
        #
        self.dof_X_r0_right = r0[-1, :, -1, -1, :, -1].dof_X
        self.dof_X_r0_left = r0[0, :, -1, 0, :, -1].dof_X
        self.dof_X_r0_top = r0[:, -1, -1, :, -1, -1].dof_X
        self.dof_X_r0_bottom = r0[:, 0, -1, :, 0, -1].dof_X

        #slices r1 dofs
        #
        self.dofs_r1_right = r1[-1, :, -1, -1, :, -1].dofs
        self.dofs_r1_left = r1[0, :, -1, 0, :, -1].dofs
        self.dofs_r1_top = r1[:, -1, -1, :, -1, -1].dofs
        self.dofs_r1_bottom = r1[:, 0, -1, :, 0, -1].dofs

        #slices r1 global coordinates
        #
        self.dof_X_r1_right = r1[-1, :, -1, -1, :, -1].dof_X
        self.dof_X_r1_left = r1[0, :, -1, 0, :, -1].dof_X
        self.dof_X_r1_top = r1[:, -1, -1, :, -1, -1].dof_X
        self.dof_X_r1_bottom = r1[:, 0, -1, :, 0, -1].dof_X

        #long slices over both domains
        #
        self.dofs_sym = vstack( ( self.dofs_r0_top, self.dofs_r1_top ) )
        self.dofs_bottom = vstack( ( self.dofs_r0_bottom, self.dofs_r1_bottom ) )
        self.dof_X_sym = vstack( ( self.dof_X_r0_top, self.dof_X_r1_top ) )
        self.dof_X_bottom = vstack( ( self.dof_X_r0_bottom , self.dof_X_r1_bottom ) )

        return tloop



if __name__ == '__main__':

    #sim_model = MRtwo()

#    interior_elems = sim_model.fe_grid_column[ 1:-1, 1:-1, :, :, :, : ].elems
#    sim_model.fe_grid_column.inactive_elems = list( interior_elems )

#    do = 'export_hinge'
#    do = 'eval' 
    do = 'ui'
#    do = 'ps' 
#    do = 'cs' 

    if do == 'eval':
        print '*** eval ***'
        sim_model.peval()

    if do == 'ui':
        print '*** ui ***'

        # input options
        #
        sim_model = MRtwo( link_type = 'exc_V_ip',
    #        sim_model.link_type = 'inc_V_ip'
            lc = 'lc_w_asym',
    #        sim_model.lc = 'lc_g'

            n_elems_xy_quarter = 17, # corresponds to 'equal_100cm'
    #        sim_model.n_elems_xy_quarter = 10 # needs to be defined before shift array otherwise 
    #                                          # linear function not evaluable

            select_hinge_pos = True, # needs to be true for shift

    #        sim_model.shift_array = array( [[2, 2, 1]] ) # shift for one quarter

            # corresponds to 'equal_100cm'
            #
            shift_array = array( [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1],
                                            [0.5, 0.5 , 2],
                                            [1.5 , 1.5, 4],
                                            [2.5, 2.5, 4],
                                            [3.5, 3.5, 4]] ),

            # element of shif_array that should not be linked
            # (corresponds to 'equal_100cm')
            #
            not_linked_elem = 0
            )
        # evaluation
        #
        sim_model.peval()

        # export
        #
#        sim_model.export_edge_u_data( filename = 'U_data.csv' )
#        sim_model.export_int_force_data( filename = 'F_int_data.csv' )

        # visualisation
        # 
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp( ibv_resource = sim_model )
        app.main()


    if do == 'export_hinge':

        #loading cases
        #         !!!! lc_shrink not working automatically needs to be done over mushroof_model class!!!!!
        #
        lc_list = [
#                   'lc_shrink' ,
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
                          'equal_100cm',
#                          'equal_200cm',
#                          'gap_100cm',
#                          'corner_2to4_25cm',
#                          'middle_0to3_25cm',
#                          'middle_0to2_25cm'
                          ]

        # V_ip within hinge exists if  freedom is linked or no linked
        #
        link_type_list = [
                          'exc_V_ip',
#                          'inc_V_ip'
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

        sim_model.select_hinge_pos = True


        for link_type in link_type_list:

            for link_case in link_case_list:

                for lc in lc_list:

                        print 'link_type:', link_type
                        print 'link_case:', link_case
                        print 'loading_case:', lc

                        # load case
                        #
                        sim_model.lc = lc

                        # shifting
                        #
                        sim_model.n_elems_xy_quarter = n_elems_quarter_xy_dict[link_case]

                        sim_model.not_linked_elem = not_linked_dict[link_case]

                        sim_model.shift_array = shift_array_dict[link_case]

                        # linking type
                        #
                        sim_model.link_type = link_type


                        # evaluation
                        #
                        sim_model.peval()

                        # export
                        #
                        sim_model.export_int_force_data( filename = link_type + '_' + link_case + '_hf_' + lc + '.csv' )
                        sim_model.export_edge_u_data( filename = link_type + '_' + link_case + '_u_' + lc + '.csv' )

    # parametric studies
    #
    elif do == 'ps':
        print "*** ps ***"
        sim_ps = SimPStudy( sim_model = sim_model )
        sim_ps.configure_traits()

    # pickle never used from me
    #
    elif do == 'pickle':
        print "*** pickle ***"

        import pickle
        filename = '/tmp/sim.pickle'
        file = open( filename, 'w' )
        pickle.dump( sim_model, file )
        file.close()
        file = open( filename, 'r' )
        sm = pickle.load( file )
        file.close()



    elif do == 'cs':

        # not working so well really bad is basically just extracting the ps study parameters
        # 
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



        cs = 'shell'
        choosen_x , choosen_z = 10, 1


        if cs == 'shell':
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
#            n_dofs_2 = n_dofs[where( input_array[:, idx_z_axes] == 3 )]

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
                       y_1,
#                       y_2
                       ][where( c_[n_dofs_0,
                       n_dofs_1,
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
            ax1.plot( n_dofs_choosen, y_choosen, marker = 'o', color = 'b', markerfacecolor = 'b', markersize = 8 )#
            xlim ( 0, 15000 )

        if cs == 'column':
            xlim ( 0, 2000 )
            ax1.plot( n_dofs_0 , y_0, color = 'b', label = '2', linewidth = 1.5 )
#            ax1.plot(n_dofs_1 ,y_1, color = 'g', label = '4', linewidth=1.5)
#            ax1.plot(n_dofs_2 ,y_2, color = 'r', label = '4', linewidth=1.5)
            ax1.plot( n_dofs_choosen.reshape( -1, 1 ), y_choosen.reshape( -1, 1 ), marker = 'o', color = 'b', markerfacecolor = 'b', markersize = 8 )#

        ylim ( -0.030, -0.025 )

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

