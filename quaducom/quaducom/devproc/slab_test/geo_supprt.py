'''
Created on Jun 16, 2010

@author: alexander
'''
from etsproxy.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, Enum, \
    Dict, Bool, Int, List

import numpy as np

from numpy import c_

import scipy as sp

from os.path import join

from math import pi


class GeoSPT(HasTraits):
    '''Geometry definition of the tappered support added from four hp-surfaces and a planar top-level.
    used for linear reduction of the stiffness of the support.
    '''

    #-----------------------------------------------------------------
    # geometric parameters of the support
    #-----------------------------------------------------------------

    #-----------------
    # geometry:
    #-----------------

    # x and y-direction 
    #
    width_supprt = Float(0.1, input = True)

    # z-direction
    #
    thickness_supprt = Float(0.02, input = True)
    
    # specify offset (translation) of the elastomer patch
    # in global coordinates
    #
    xyoffset = Float(0.0, input = True)
    zoffset = Float(0.0, input = True)

    def __call__(self, pts):
        print '*** geo_slab_test called ***' 
        
        L = self.width_supprt 
        t = self.thickness_supprt
        
        x_, y_, z_ = pts.T
                
        x = np.zeros_like(x_)
        y = np.zeros_like(y_)
        z = np.zeros_like(y_)

        # 1. quadrant     
        #
        bool_x = x_ >= 0.5
        bool_y = y_ >= 0.5 
        bool_xy = bool_x * bool_y
        idx_xy = np.where( bool_xy == 1. )[0]
        rx = (x_[ idx_xy ]-0.5)/0.5
        ry = (y_[ idx_xy ]-0.5)/0.5
        dz_red_ = 0.9 * (rx + ry - rx * ry )
        x[ idx_xy ] = L / 2. + rx * L / 2.
        y[ idx_xy ] = L / 2. + ry * L / 2.
        z[ idx_xy ] = dz_red_ * t

        # 2. quadrant     
        #
        bool_x = x_ >= 0.5
        bool_y = y_ < 0.5 
        bool_xy = bool_x * bool_y
        idx_xy = np.where( bool_xy == 1. )[0]
        rx = (x_[ idx_xy ]-0.5)/0.5
        ry = y_[ idx_xy ]/0.5
        dz_red_ = 0.9 * (rx + (1-ry) - rx * (1-ry) )
        x[ idx_xy ] = L / 2. + rx * L / 2.
        y[ idx_xy ] = ry * L / 2.
        z[ idx_xy ] = dz_red_ * t

        # 3. quadrant     
        #
        bool_x = x_ < 0.5
        bool_y = y_ < 0.5 
        bool_xy = bool_x * bool_y
        idx_xy = np.where( bool_xy == 1. )[0]
        rx = x_[ idx_xy ]/0.5
        ry = y_[ idx_xy ]/0.5
        dz_red_ = 0.9 * ((1-rx) + (1-ry) - (1-rx) * (1-ry) )
        x[ idx_xy ] = rx * L / 2.
        y[ idx_xy ] = ry * L / 2.
        z[ idx_xy ] = dz_red_ * t

        # 4. quadrant     
        #
        bool_x = x_ < 0.5
        bool_y = y_ >= 0.5 
        bool_xy = bool_x * bool_y
        idx_xy = np.where( bool_xy == 1. )[0]
        rx = x_[ idx_xy ]/0.5
        ry = (y_[ idx_xy ]-0.5)/0.5
        dz_red_ = 0.9 * ((1-rx) + ry - (1-rx) * ry )
        x[ idx_xy ] = rx * L / 2.
        y[ idx_xy ] = L / 2. + ry * L / 2.
        z[ idx_xy ] = dz_red_ * t

        # apply z-reduction only at lower row     
        #
        idx_z = np.where( z_ == 1. )[0]
        z[ idx_z ] = z_[ idx_z ] * t
        
        pts = np.c_[x, y, z]
        
        # specify the offset of the position of the support with respect to the global coordinate system
        #
        xyoffset = self.xyoffset
        zoffset = self.zoffset
        offset_arr = np.array([xyoffset, xyoffset, zoffset]) 
        pts = pts + offset_arr

        return pts


if __name__ == '__main__':

    from numpy import mgrid, c_, hstack, vstack, shape
    from etsproxy.mayavi import mlab

    geo_supprt = GeoSPT()

    shape_xy = 2
    shape_z = 1

    grid = mgrid[0:1:complex(0, shape_xy + 1),
                 0:1:complex(0, shape_xy + 1 ),
                 0:1:complex(0, shape_z + 1)]

    X, Y, Z = grid

    gpoints = c_[ X.flatten(), Y.flatten(), Z.flatten() ]
    
    mlab.figure(bgcolor=(1.,1.,1.,))
    fp1 = geo_supprt(gpoints)
    mlab.points3d(fp1[:, 0], fp1[:, 1], fp1[:, 2],
                   scale_factor = 0.003 ,
                   resolution = 8)
    mlab.axes()
    mlab.show()