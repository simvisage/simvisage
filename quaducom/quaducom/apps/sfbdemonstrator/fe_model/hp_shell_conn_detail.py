'''
Created on Jun 16, 2010

@author: alexander
'''


from etsproxy.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, Enum, \
    Dict, Bool, Int

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, mgrid, arange, \
    where, sum, sin, cos, vstack, hstack, argmax, newaxis, size, \
    shape, sqrt, frompyfunc, ones_like, zeros_like, ones, any, all, \
    sort, argsort, concatenate, add

from .rsurface_reader import \
    read_rsurface, normalize_rsurfaces

from os.path import join

from math import pi
from matplotlib.pyplot import *

# Interpolation
from scipy.interpolate import Rbf

def delete_second_rows( arr, nx = 20, ny = 20 ):
    '''Remove the second and second to last column and row from the regular grid of points.
    for lowerface_4x4m.robj and upperface_4x4m.robj data file
    '''
    va = arr
    va = va.reshape( ( nx, ny, 3 ) )
    va_a = vstack( [ va[ :1, :, : ], va[ 2:-2, :, : ], va[-1:, :, : ] ] )
    va_b = hstack( [ va_a[ :, :1, : ], va_a[ :, 2:-2, : ], va_a[ :, -1:, : ] ] )
    n_size = va_b.shape[0] * va_b.shape[1]
    va = va_b.reshape( ( n_size, 3 ) )
    return va

class HPShell( HasTraits ):
    '''Geometry definition.
    '''

    #===========================================================================
    # geometric variables and values
    #===========================================================================

    # part of mushroof
    #
    # @todo: "four" is not supported by "n_elems_xy_dict"in mushroff_model
    mushroof_part = Enum( 'detail', 'quarter', 'one', 'four' )

    # origin
    # @todo: define as "_default"
    #
    X0 = Array( float, value = [ 0., 0., 0. ] )

    # element properties of grid
    #
    n_elems_xy = Int( 6 )
    n_elems_z = Int( 3 )

    n_elems_xy_quarter = Property( Int )
    @cached_property
    def _get_n_elems_xy_quarter( self ):
        return  self.n_elems_xy / self.scale_size

    # standard array for column shift
    # shift array shifts the elements towards the defined coordinates
    # [x_shift,y_shift,element_between]
    # array needs to have the right shape (:,3)!!
    #
    # @todo: define as "_default"
    #
    shift_array = Array ( float,
                         shape = ( None, 3 ) ,
                         value = [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1]] )

    # dimensions of the shell for one quarter of mush_roof
    # @todo: add "_quarter"
    #
    length_x = Float( 4.0 ) # [m]
    length_y = Float( 4.0 ) # [m]
    length_z = Float( 1.062 ) # [m]
    t_shell = Float( 0.06 ) # [m]
    width_column = Float( 0.45 ) # [m]

    length_x_detail = Float( 1.5 ) # [m]
    length_y_detail = Float( 1.5 ) # [m]

    scale_size_detail_x = Property( Float )
    def _get_scale_size_detail_x( self ):
        return self.length_x_detail / self.length_x * 2.

    scale_size_detail_y = Property( Float )
    def _get_scale_size_detail_y( self ):
        return self.length_y_detail / self.length_y * 2.


    # scale factor for different mushroof parts
    # Defines the proportion between the lenght of the model
    # with respect to the length of a quarter shell as the 
    # basic substructure of which the model consists of.
    # @todo: add "depend_on" or remove "cached_property"
    # @todo: move to class definition of "mushroof_model" not in "HPShell"
    #        (also see there "n_elems_dict" with implicit "scale_factor")
    #
    scale_size = Property( Float )
#    @cached_property
    def _get_scale_size( self ):
#        scale_size_detail = self.lenght_x_detail / self.length_x
        scale_dict = {'detail' : self.scale_size_detail_x, 'quarter': 1.0 , 'one' : 2.0, 'four' : 4.0 }
        return scale_dict[self.mushroof_part]

    # factor to scale delta_h (inclination of the shell)
    # The thickness remains unchanged by this scaling, e.g. self.thickness = 0.06 [m]
    #
    delta_h_scalefactor = Float( 1.00 ) # [-] 

    # shift of column elements
    #
    shift_elems = Bool( True )

    # const_edge element operator
    # (for non-linear analysis only, where an element layer of constant 
    # thickness is needed to simulate the reinforced behaviour of the 
    # concrete.
    #
    const_edge_elem = Bool( False )
    t_edge = Float ( 0.03 ) # [m]
    n_elems_edge = Int( 1 ) #number of dofs used for edge refinement  

    #===========================================================================
    # reading options
    #===========================================================================

    # "lowerface_cut_off" - option replaces constant height for the coordinates 
    # which connect to the column (this cuts of the shell geometry horizontally 
    # at the bottom of the lower face of the shell geometry.
    # Option should be used for the robj-file with 4x4m geometry
    #
    cut_off_lowerface = Bool( True )

    # corresponds to the delta in the geometry .obj-file with name '4x4m' as a cut off
    #
    delta_h = Float( 1.00 ) # [m]

    # choose geometric file (obj-data file)
    #
    geo_input_name = Enum( '4x4m', '02' )

    # filter for '4x4m' file needs to be done to have regular grid
    # in order to rbf-function leading to stable solution without oscilation
    #
    geo_filter = Dict( {'4x4m' : delete_second_rows } )

    def _read_arr( self, side = 'lowerface_' ):
        file_name = side + self.geo_input_name + '.robj'
        file_path = join( 'geometry_files', file_name )
        v_arr = read_rsurface( file_path )
        filter = self.geo_filter.get( self.geo_input_name, None )
        if filter != None:
            v_arr = list(filter( v_arr ))
        return v_arr

    # array of the vertex positions in global 
    # x,y,z-coordinates defining the lower surface of the shell 
    vl_arr = Property( Array( float ) )
    @cached_property
    def _get_vl_arr( self ):
        vl_arr = self._read_arr( 'lowerface_' )
        if self.cut_off_lowerface == True:
            print('--- lower face z-coords cut off ---')

            # z-values of the coords from the lower face are cut off. 
            # From the highest z-coordinate of the lower face the vertical
            # distance is 1 m (=delta h). At this limit the lower face is
            # cut off. Global z coordinate is assumed to point up.
            #
            vl_z_max = max( vl_arr[:, 2] )
            vl_z_min = vl_z_max - self.delta_h
            vl_arr_z = where( vl_arr[:, 2] < vl_z_min, vl_z_min, vl_arr[:, 2] )
            vl_arr = c_[vl_arr[:, 0:2], vl_arr_z]
        return vl_arr

    # array of the vertex positions in global 
    # x,y,z-coordinates defining the upper surface of the shell
    vu_arr = Property( Array( float ) )
    @cached_property
    def _get_vu_arr( self ):
        return self._read_arr( 'upperface_' )

    #------------------------------------------------------------------------------ 
    # hp_shell geometric transformation 
    #------------------------------------------------------------------------------ 

    def __call__( self, points ):
        '''Return the global coordinates of the supplied local points.
        '''

        # number of local grid points for each coordinate direction
        # values must range between 0 and 1
        #
        xi, yi, zi = points[:, 0], points[:, 1], points[:, 2]
        print("xi", xi)
        print("xi.shape", xi.shape)

        # size of total structure
        #
        # @todo: move to class definition of "mushroof_model" and send to "__call__"
        scale_size = self.scale_size
        print("scale_size", scale_size)

        # @todo: add "_quarter" (see above)
        length_x_tot = self.length_x * scale_size
        length_y_tot = self.length_y * scale_size
        n_elems_xy_quarter = self.n_elems_xy_quarter

        # distance from origin for each mushroof_part
        #
        def d_origin_fn( self, coords ):
#            if self.mushroof_part == 'quarter':
#                return coords
#            if self.mushroof_part == 'one':
#                return abs( 2.0 * coords - 1.0 )

            if self.mushroof_part == 'detail':
                print('in d_origin_fn')
                return abs( 1.0 * coords - 0.5 ) * scale_size

#            # @todo: corresponding "scale_factor" needs to be added 
#            #        in order for this to work
#            if self.mushroof_part == 'four':
#                return  where( coords < 0.5, abs( 4 * coords - 1 ), abs( -4 * coords + 3 ) )


        # values are used to calculate the z-coordinate using RBF-function of the quarter
        # (= values of the distance to the origin as absolute value)
        #
        xi_rbf = d_origin_fn( self, xi )
        print('xi_rbf', xi_rbf)
        yi_rbf = d_origin_fn( self, yi )

        # normalized coordinates of the vertices for lower- and upperface
        # NOTE: the underline character indicates a normalized value 
        #
        vl_arr_, vu_arr_ = normalize_rsurfaces( self.vl_arr,
                                                self.vu_arr )

        # use a radial basis function approximation (rbf) (i.e. interpolation of
        # scattered data) based on the normalized vertex points of the lower face
        #
        x_ = vl_arr_[:, 0]

        # flip the orientation of the local coordinate system in the
        # corresponding y-direction depending on the data file
        #
        geo_input_name = self.geo_input_name
        if geo_input_name == '4x4m':
            y_ = vl_arr_[:, 1]
        else:
            y_ = 1 - vl_arr_[:, 1]

        z_ = vl_arr_[:, 2]
        rbf = Rbf( x_, y_, z_, function = 'cubic' )

        # get the z-value at the supplied local grid points
        # of the lower face
        #
        zi_lower_ = rbf( xi_rbf, yi_rbf )

        # use a radial basis function approximation (rbf) (i.e. interpolation of
        # scattered data) based on the normalized vertex points of the upper face
        #
        x_ = vu_arr_[:, 0]

        # flip the orientation of the local coordinate system in the
        # corresponding y-direction depending on the data file
        #
        geo_input_name = self.geo_input_name
        if geo_input_name == '4x4m':
            y_ = vu_arr_[:, 1]
        else:
            y_ = 1 - vu_arr_[:, 1]

        z_ = vu_arr_[:, 2]
        rbf = Rbf( x_, y_, z_, function = 'cubic' )

        # get the z-value at the supplied local grid points
        # of the upper face
        # 
        # note that zi_upper_ is a normalized coordinate!
        #
        zi_upper_ = rbf( xi_rbf, yi_rbf )

        # thickness is multiplied by the supplied zi coordinate
        #
        z_ = ( zi_lower_ + ( zi_upper_ - zi_lower_ ) * zi / self.delta_h_scalefactor ) * self.delta_h_scalefactor

        # coordinates of origin
        #
        X, Y, Z = self.X0

        print('--- geometric transformation done ---')

        # multiply the local grid points with the real dimensions in order to obtain the 
        # global coordinates of the mushroof_part:
        #
        return c_[ X + xi * length_x_tot, Y + yi * length_y_tot, Z + z_ * self.length_z ]



if __name__ == '__main__':

    from numpy import mgrid, c_, hstack, vstack, shape
    from etsproxy.mayavi import mlab
    hp = HPShell()
    hp.n_elems_xy = 10
    hp.n_elems_z = 1
    hp.n_elems_col = 1
    hp.n_elems_edge = 1

    hp.t_shell = 0.06
    hp.t_edge = 0.03
    hp.width_column = 0.45

    hp.const_edge_elem = False
    hp.shift_elems_column = True
    hp.cut_off_lowerface = True
    hp.geo_input_name = '4x4m'
    hp.mushroof_part = 'one'

    X, Y, Z = mgrid[0:1:11j, 0:1:11j , 0:1:2j]
    gpoints = c_[ X.flatten(), Y.flatten(), Z.flatten() ]

    fp1 = hp( gpoints )
    print(shape( fp1 ))
    mlab.points3d( fp1[:, 0], fp1[:, 1], fp1[:, 2], scale_factor = 0.05 , resolution = 8 )

    mlab.show()
