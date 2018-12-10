'''
Created on Dec 22, 2009

@author: alexander
'''

from numpy import array, loadtxt, vstack, c_


def read_rsurface( filename ):
    '''Read the rhino3D data format. 
    Extract solely the vertex coordinates.
    NOTE: the file ending .obj already exists in eclipse
          therefore the ending of the input data file 
          has been changed to .robj
    '''
    # file containing the geometry data
    geo_file = file( filename )

    vertices_list = []
    for line in geo_file:
        # the coordinates of the vertices are indicated in 
        # the rhino3D-input file by the first character of 
        # the line given being "v"
        if line[0] == "v" and line[0:2] != "vp":
            # split each line and cut of the "v"-key
            coords = line.split()
            x, y, z = list(map( float, coords[1:] ))
            vertices_list.append( [x, y, z] )

    # convert the ordered list of x,y,z-coordinats in an 
    # array with 3 columns
    vertices_arr = array( vertices_list )
    return vertices_arr


def normalize_rsurfaces( vl_arr, vu_arr ):
    '''normalize the geometry data by the division by the 
    maximum ranges in x, y, and z-direction, respectively.
    '''
    # index "l" indicates lower face of the shell geometry
    # index "u" indicates upper face of the shell geometry
    # vl_arr is the array of vertices of the lower face

    # derive the number of vertices used for the lowerface
    vl_shape = vl_arr.shape[0]

    v_arr = vstack( [vl_arr, vu_arr] )

    vx_arr, vy_arr, vz_arr = v_arr[:, 0], v_arr[:, 1], v_arr[:, 2]
    xmin, xmax = min( vx_arr ), max( vx_arr )
    ymin, ymax = min( vy_arr ), max( vy_arr )
    zmin, zmax = min( vz_arr ), max( vz_arr )

    xrange = xmax - xmin
    yrange = ymax - ymin
    zrange = zmax - zmin

    vx_arr_ = ( vx_arr - xmin ) * 1 / xrange
    vy_arr_ = ( vy_arr - ymin ) * 1 / yrange
    vz_arr_ = ( vz_arr - zmin ) * 1 / zrange

    v_arr_ = c_[ vx_arr_, vy_arr_, vz_arr_ ]

    # return the normalized vertex-coordinates 
    # indicated by he underline character ("_")
    # for the lowerface and upperface
    return v_arr_[:vl_shape, :], v_arr_[vl_shape:, :]


if __name__ == '__main__':

    from os.path import join

    dir = 'geometry_files'

    # read in the global coordinates of the vertices of the
    # lower and upper shell surface 
    geo_lowerfile = 'lowerface_350x350cm_neu2.robj'
    geo_upperfile = 'upperface_350x350cm_neu2.robj'
    vl_arr = read_rsurface( join( dir, geo_lowerfile ) )
    vu_arr = read_rsurface( join( dir, geo_upperfile ) )

    # normalize the vertex-coordinates
    vl_arr_, vu_arr_ = normalize_rsurfaces( vl_arr, vu_arr )

    # print the normalized coordinates 
    print('vl_arr - normalized')
    print(vl_arr)
    print('\n')

    print('vu_arr - normalized')
    print(vu_arr_)


