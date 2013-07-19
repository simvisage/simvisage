'''
Created on Jun 16, 2010

@author: alexander
'''
from etsproxy.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, Enum, \
    Dict, Bool, Int, List

import numpy as np

import scipy as sp

from os.path import join

from math import pi


class GeoST(HasTraits):
    '''Geometry definition of the slab test with round load introduction area 
    corresponding to steel plate in the test setup.
    '''

    #-----------------------------------------------------------------
    # geometric parameters of the slab
    #-----------------------------------------------------------------
    # NOTE: coordinate system is placed where the symmetry planes cut each other, 
    # i.e the center of the load introduction area (=middle of steel plate)

    # discretization of total slab in x- and y-direction (region 'L')
    #
    shape_xy = Int(10, input = True)

    # discretization of the load introduction plate (region 'R')
    #
    shape_R = Int(2, input = True)
        
    # discretization in z-direction 
    # (thickness direction):
    #
    shape_z = Int(1, input = True)

    #-----------------
    # geometry:
    #-----------------

    # x and y-direction 
    #
    length_quarter = Float(0.625, input = True)

    # Radius of load introduction plate 
    #
    radius_plate = Float(0.10, input = True)
    
    # z-direction
    #
    thickness = Float(0.06, input = True)

#    # used regular discretization up to y = L1
#    # (by default use regular discretization up to support)
#    #
#    L1 = Float(0.30, input = True)

    def __call__(self, pts):
        print '*** geo_slab_test called ***' 
        
        x_, y_, z_ = pts.T
                
        R = self.radius_plate
        L = self.length_quarter
        t = self.thickness
        
        #-------------------------------------------
        # transformation to global coordinates 
        #-------------------------------------------

        x = np.zeros_like(x_)
        y = np.zeros_like(y_)
        z = z_ * t

        # ratio of the discretization, i.e. number of elements for each region     
        #
        r_ = 1.* self.shape_R / self.shape_xy

        # 1. quadrant     
        #
        bool_x = x_ > r_
        bool_y = y_ > r_ 
        bool_xy = bool_x * bool_y
        idx_xy = np.where( bool_xy == 1. )[0]
        x[ idx_xy ] = R + (x_[ idx_xy ]-r_)/(1-r_) * (L - R)
        y[ idx_xy ] = R + (y_[ idx_xy ]-r_)/(1-r_) * (L - R)

        # 2. quadrant     
        #
        bool_x = x_ > r_
        bool_y = y_ <= r_ 
        bool_xy = bool_x * bool_y
        idx_xy = np.where( bool_xy == 1. )[0]
        xR = R * np.cos( y_[ idx_xy ]/r_ * np.pi/4. ) 
        x[ idx_xy ] = xR + (x_[ idx_xy ]-r_)/(1-r_) * (L - xR)
        y[ idx_xy ] = y_[ idx_xy ]/r_ * R

        # 3. quadrant (mesh in load introduction area)
        #
        bool_x = x_ <= r_
        bool_y = y_ <= r_
        bool_xy = bool_x * bool_y 
        idx_xy = np.where( bool_xy == 1. )[0]
        xR = R * np.cos( y_[ idx_xy ]/r_ * np.pi/4. ) 
        yR = R * np.cos( x_[ idx_xy ]/r_ * np.pi/4. ) 
        x[ idx_xy ] = x_[ idx_xy ]/r_ * xR 
        y[ idx_xy ] = y_[ idx_xy ]/r_ * yR
        
        # 4. quadrant     
        #
        bool_x = x_ <= r_
        bool_y = y_ > r_ 
        bool_xy = bool_x * bool_y
        idx_xy = np.where( bool_xy == 1. )[0]
        x[ idx_xy ] = x_[ idx_xy ]/r_ * R
        yR = R * np.cos( x_[ idx_xy ]/r_ * np.pi/4. ) 
        y[ idx_xy ] = yR + (y_[ idx_xy ]-r_)/(1-r_) * (L - yR )
  
        pts = c_[x, y, z]
        return pts




if __name__ == '__main__':

    from numpy import mgrid, c_, hstack, vstack, shape
    from etsproxy.mayavi import mlab

    st = GeoST()

    grid = mgrid[0:1:complex(0, st.shape_xy + 1),
                 0:1:complex(0, st.shape_xy + 1),
                 0:1:complex(0, st.shape_z )]

    X, Y, Z = grid

    gpoints = c_[ X.flatten(), Y.flatten(), Z.flatten() ]
    
    mlab.figure(bgcolor=(1.,1.,1.,))
    fp1 = st(gpoints)
    mlab.points3d(fp1[:, 0], fp1[:, 1], fp1[:, 2],
                   scale_factor = 0.02 ,
                   resolution = 8)
    mlab.axes()
    mlab.show()