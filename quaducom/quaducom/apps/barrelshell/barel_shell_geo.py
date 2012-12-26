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


class BarelShellGeo(HasTraits):
    '''Geometry definition of the barrel shell.
    '''

    #-----------------------------------------------------------------
    # geometric parameters of the shell
    #-----------------------------------------------------------------
    # NOTE: coordinate system is placed where the symmetry planes cut each other:

    # discretization in length-direction (y-direction)
    #
    shape_l = Int(20, input = True)

    # discretization in arc-direction 
    # (use aquidistant spacing along the arc):
    #
    shape_s = Int(4, input = True)

    # discretization in z-direction 
    # (thickness direction):
    #
    shape_z = Int(1, input = True)

    #-----------------
    # geometry:
    #-----------------

    # y-direction 
    # NOTE: (length_quarter = half of the total length (due to symmetry!)) 
    #
    length_quarter = Float(2.20, input = True)

    # x-direction 
    # NOTE: (width_quarter = half of the total arc width (due to symmetry!)) 
    width_quarter = Float(1.07, input = True)

    # z-direction
    #
    arc_height = Float(0.50, input = True)
    thickness = Float(0.02, input = True)

    # specify the length where the cutting of the barrel shell starts
    # (by default this corresponds to 'arc_height' 
    # yielding a 45 deg angle of the rounded edge in the y-z-plane)
    #
    Lr = Float(0.50, input = True)

    # used regular discretization up to y = L1
    # (by default use regular discretization up to support)
    #
    L1 = Float(1.30, input = True)

    def __call__(self, pts):
        print '*** geo_barrel_shell called ***'

        x_, y_, z_ = pts.T

        L = self.length_quarter
        b = self.width_quarter
        f = self.arc_height
        t = self.thickness

        #-------------------------------------------
        # transformation for 'cylinder coordinates' 
        #-------------------------------------------

        # calculate the arc radius:
        # 
        R = f / 2. + b ** 2 / (2.*f)

        # calculate the arc angle [rad]
        beta = sp.arctan(b / (R - f))

        # cylinder coordinates of the barrel shell
        #
        y = y_ * L
        x = (R - z_ * t) * np.sin(x_ * beta)
        z = (R - z_ * t) * np.cos(x_ * beta) - R + f

        #-------------------------------------------
        # cut of free edge by 45 deg
        #-------------------------------------------

        # rounded length
        Lr = self.Lr

        # length to be substracted in y-direction (linear relation with respect to the z-axis)
        #
        delta_yr = (1. - z / f) * Lr

        # used regular discretization up to y = L1
        L1 = self.L1

        # substract 'yr' for y_ = 1.0 (edge) and substract 0. for y_ = L1/L
        # and interpolate linearly within 'L' and 'L1'
        #
        idx_r = np.where(y_ > L1 / L)[0]
        y[ idx_r ] -= ((y_[ idx_r ] - L1 / L) / (1.0 - L1 / L) * delta_yr[ idx_r ])

        pts = np.c_[x, y, z]

        return pts

if __name__ == '__main__':

    from etsproxy.mayavi import mlab

    bs = BarelShellGeo()

    grid = np.mgrid[0:1:complex(0, bs.shape_s + 1),
                 0:1:complex(0, bs.shape_l + 1),
                 0:1:complex(0, bs.shape_z + 1)]

    X, Y, Z = grid

    gpoints = np.c_[ X.flatten(), Y.flatten(), Z.flatten() ]
    fp1 = bs(gpoints)
    mlab.points3d(fp1[:, 0], fp1[:, 1], fp1[:, 2],
                   scale_factor = 0.05 ,
                   resolution = 8)

    mlab.show()
