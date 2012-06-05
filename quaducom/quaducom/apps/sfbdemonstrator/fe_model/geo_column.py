'''
Created on Jul 30, 2010

@author: abach
'''

from etsproxy.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, Enum, \
    Dict, Bool, Int

from numpy import \
    sin, cos, mgrid, where, c_

from math import \
    pi

class GEOColumn( HasTraits ):

    h_col = Float( 3.00, unit = 'm' )
    X0 = Array( float, value = [ 4., 4., -3.0 ] )
    width_top = Float( 0.45, unit = 'm' )
    width_bottom = Float( 0.35, unit = 'm' )
    r_pipe = Float( 0.1, unit = 'm' )

    def __call__( self, points ):

        xi, yi, zi = points[:, 0], points[:, 1], points[:, 2]

        # grid from 0 to 1, shift origin for rotation
        # 
        xi -= 0.5
        yi -= 0.5

        # setting of global coordinates different width over h_col 
        # f_QS:
        #
        xi = ( ( self.width_top - self.width_bottom ) * zi + self.width_bottom ) * xi
        yi = ( ( self.width_top - self.width_bottom ) * zi + self.width_bottom ) * yi
        zi *= self.h_col

        # rotation of 45 with global coordinates
        #
        x = cos( pi / 4.0 ) * xi + sin( pi / 4.0 ) * yi
        y = -sin( pi / 4.0 ) * xi + cos( pi / 4.0 ) * yi
        z = zi + self.X0[2]

        # @TODO: kill element in the center of the columns:
        #shift of internal elements
        #
#        r_0 = where( ( xi ** 2 + yi ** 2 ) ** 0.5 >= self.r_pipe, self.r_pipe , ( xi ** 2 + yi ** 2 ) ** 0.5 )
#        scale_r = where ( r_0 == 0, 0, self.r_pipe / r_0 )
#        x = scale_r * x + self.X0[0]
#        y = scale_r * y + self.X0[1]

        x = x + self.X0[0]
        y = y + self.X0[1]
        return c_[x, y, z]
#        return c_[xi + self.X0[0], yi+ self.X0[1], z]
if __name__ == '__main__':

    col = GEOColumn()

    X, Y, Z = mgrid[0:1:11j, 0:1:11j , 0:1:4j]
    gpoints = c_[ X.flatten(), Y.flatten(), Z.flatten() ]

    print col( gpoints )
