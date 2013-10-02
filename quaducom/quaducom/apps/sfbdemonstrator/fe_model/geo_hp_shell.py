'''
Created on Jun 16, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, Bool, File

from mathkit.mfn import MFnLineArray

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, mgrid, arange, \
    where, sign, fabs

from math import sqrt, asin, acos

from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

# Interpolation
from scipy.interpolate import Rbf

from simiter.sim_pstudy import ISimModel, SimOut, SimPStudy

class GeoHPShell( HasTraits ):
    '''Geometry definition.
    '''
    # dimensions of the shell structure [m]
    # (only one quart of the shell structure)
    #
    # NOTE: lenth_z = 1.0 m + 0.062 m  
    length_x = Float( 4.0 )
    length_y = Float( 4.0 )
    length_z = Float( 1.062 )
    thickness = Float( 0.062 )

    center = Array( float, value = [0, 0, 0] )

    # cut of the z-coordinates of the lowerface if set to True
    #
    cut_off_lowerface = Bool( True )

    #geo_file_upper = File( 'upperface_02.robj' )
    geo_file_upper = File( 'geometry_files/upperface_4x4m.robj' )

    # array of the vertex positions in global 
    # x,y,z-coordinates defining the upper surface of the shell
    #
    vu_arr = Property( Array( float ), depends_on = 'geo_file_upper' )
    @cached_property
    def _get_vu_arr( self ):
        return read_rsurface( self.geo_file_upper )

    #geo_file_lower = File( 'lowerface_02.robj' )
    geo_file_lower = File( 'geometry_files/lowerface_4x4m.robj' )

    # array of the vertex positions in global 
    # x,y,z-coordinates defining the lower surface of the shell 
    #
    vl_arr = Property( Array( float ), depends_on = 'geo_file_lower' )
    @cached_property
    def _get_vl_arr( self ):
        vl_arr = read_rsurface( self.geo_file_lower )

        if self.cut_off_lowerface == True:
            print '--- lower face z-coords cut off ---'

            # z-values of the coords from the lower face are cut off. 
            # From the highest z-coordinate of the lower face the vertical
            # distance is 1 m (=delta h). At this limit the lower face is
            # cut off. Global z coordinate is assumed to point up.
            #
            delta_h = 1.0
            vl_z_max = max( vl_arr[:, 2] )
            vl_z_min = vl_z_max - delta_h
            vl_arr_z = where( vl_arr[:, 2] < vl_z_min, vl_z_min, vl_arr[:, 2] )
            vl_arr = c_[vl_arr[:, 0:2], vl_arr_z]

        return vl_arr

    def __call__( self, points ):
        '''Return the global coordinates of the supplied local points.
        '''

        # number of local grid points for each coordinate direction
        # values must range between 0 and 1

        Xc, Yc, Zc = self.center

        Xi, Yi, Zi = points[:, 0] - Xc, points[:, 1] - Yc, points[:, 2] - Zc

        sxi, syi = sign( Xi ), sign( Yi )
        xi, yi, zi = fabs( Xi ), fabs( Yi ), Zi

        # normalized coordinates of the vertices for lower- and upperface
        # NOTE: the underline character indicates a normalized value 
        vl_arr_, vu_arr_ = normalize_rsurfaces( self.vl_arr,
                                                self.vu_arr )

        # use a radial basis function approximation (rbf) (i.e. interpolation of
        # scattered data) based on the normalized vertex points of the lower face
        x_ = vl_arr_[:, 0] * self.length_x

        # for file 'upperface02.robj' and 'lowerface02.robj' the 
        # orientation of the local coordinate system in the
        # corresponding y-direction needs to be fliped
        # e.g. y_ = 1 - vl_arr_[:, 1]
        # for the file 'upperface_4x4m.robj' and 'lowerface_4x4m.robj'
        # this is not the case!
        #
        y_ = vl_arr_[:, 1] * self.length_y

        z_ = vl_arr_[:, 2] * self.length_z

        rbf = Rbf( x_, y_, z_, function = 'cubic' )
        # get the z-value at the supplied local grid points
        # of the lower face
        zi_lower_ = rbf( xi, yi )

        # use a radial basis function approximation (rbf) (i.e. interpolation of
        # scattered data) based on the normalized vertex points of the upper face
        x_ = vu_arr_[:, 0] * self.length_x

        # flip the orientation of the local coordinate system in the
        # corresponding y-direction 
        # for file 'upperface02.robj' and 'lowerface02.robj' the 
        # orientation of the local coordinate system in the
        # corresponding y-direction needs to be fliped
        # e.g. y_ = 1 - vl_arr_[:, 1]
        # for the file 'upperface_4x4m.robj' and 'lowerface_4x4m.robj'
        # this is not the case!
        #
        y_ = vu_arr_[:, 1] * self.length_y

        z_ = vu_arr_[:, 2] * self.length_z

        rbf = Rbf( x_, y_, z_, function = 'cubic' )
        # get the z-value at the supplied local grid points
        # of the upper face
        zi_upper_ = rbf( xi, yi )

        # thickness is multiplied by the supplied zi coordinate
        z_ = zi_lower_ + ( zi_upper_ - zi_lower_ ) * zi

        # multiply the local grid points with the real dimensions in order to obtain the 
        # global coordinates of the structure:
        #
        Xi, Yi, Zi = sxi * xi, syi * yi, z_ # * self.length_z
        return c_[ Xi + Xc, Yi + Yc , Zi + Zc  ]


if __name__ == '__main__':
    hp = GeoHPShell()

    points = arange( 21 ).reshape( 7, 3 ) / 20.
    print points
    print hp( points )
    
    
    from numpy import mgrid, c_, hstack, vstack, shape
    from etsproxy.mayavi import mlab
    fp1 = hp( points )
    mlab.points3d( fp1[:, 0], fp1[:, 1], fp1[:, 2], scale_factor = 0.05 , resolution = 8 )
    mlab.show()