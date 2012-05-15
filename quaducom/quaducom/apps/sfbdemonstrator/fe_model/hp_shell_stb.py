'''
Created on Jun 16, 2010

@author: andreas
'''

from enthought.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, Enum, \
    Dict, Bool, Int

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, mgrid, arange, \
    where, sum, sin, cos, vstack, hstack, argmax, newaxis, size, \
    shape, sqrt, frompyfunc, ones_like, loadtxt, arange, sqrt, \
    zeros, arctan, sin, cos, ones_like, \
    vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
    copy, c_, where, mean, arctan, cos, min, argmin

from os.path import join

from math import pi

from enthought.traits.api import \
    HasTraits, Float, Array, Bool, Enum, Dict

from numpy import \
    c_, ix_, mgrid, transpose, shape

from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

# Interpolation
from scipy.interpolate import Rbf

from math import pi

import csv

from enthought.mayavi.mlab import colorbar, show, points3d
from enthought.mayavi.api import Engine

def delete_second_rows( arr, nx = 20, ny = 20 ):
    '''Remove the second and second to last column and row from the regular grid of points.
    for lowerface_4x4m.robj and upperface_4x4m.robj data file
    '''
    print '--- deleting second rows ---'
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
    # dimensions of the shell structure [m]
    # (only one quart of the shell structure)
    #
    # NOTE: lenth_z = 1.0 m + 0.062 m = 1.062  
    # NOTE: lenth_z = 0.865 m + 0.062 m = 0.927  
    length_x = Float( 3.50 )
    length_y = Float( 3.50 )
    length_z = Float( 0.927 )

    # corresponds to the delta in the geometry obj file '4x4m'
    delta_h = Float( 0.865 ) # [m]

    # factor to scale height of lower surface 
    # thickness remains unchanged as 0.06 m
    #
    delta_h_scalefactor = Float( 1.00 ) # [-] 

    # cut of the z-coordinates of the lowerface if set to True
    #
    cut_off_lowerface = Bool( True )

    geo_input_name = Enum( '350x350cm' )
    geo_filter = Dict( {'4x4m' : delete_second_rows } )

    def _read_arr( self, side = 'lowerface_' ):
        file_name = side + self.geo_input_name + '.robj'
        file_path = join( 'geometry_files', file_name )
        v_arr = read_rsurface( file_path )

        filter = self.geo_filter.get( self.geo_input_name, None )
        if filter != None:
            v_arr = filter( v_arr )
        return v_arr

    # array of the vertex positions in global 
    # x,y,z-coordinates defining the lower surface of the shell 
    vl_arr = Property( Array( float ) )
    @cached_property
    def _get_vl_arr( self ):
        vl_arr = self._read_arr( 'lowerface_' )
        if self.cut_off_lowerface == True:
            print '--- lower face z-coords cut off ---'

            # z-values of the coords from the lower face are cut off. 
            # From the highest z-coordinate of the lower face the vertical
            # distance is 1 m (=delta h). At this limit the lower face is
            # cut off. Global z coordinate is assumed to point up.
            #
            vl_z_max = max( vl_arr[:, 2] )
            if self.geo_input_name == '4x4m':
                # NOTE: the global z-coordinates are given in the geo data file in [m]
                # no conversion of unites necessary (self.delta_h is given in [m])
                delta_h = self.delta_h
            elif self.geo_input_name == '350x350cm':
                # NOTE: the global z-coordinates are given in the geo data file in [cm]
                # convert delta_h from [m] to [cm]
                #
                delta_h = self.delta_h * 100.
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

    def get_mid_surface_and_thickness( self, points, perpendicular_t = True ):
        '''Return the global coordinates of the supplied local points.
        '''
        print '*** get mid surface and thickness ***'

        #-----------------------------------------------
        # get the global coordinates as defined in the 
        # input file and transform them to the coordinate
        # system of the master quarter
        #-----------------------------------------------

        # 
        if self.geo_input_name == '350x350cm':
            X0 = [3.50, 3.50, 0.]
        else:
            X0 = [0., 0., 0.]

        # number of global grid points for each coordinate direction
        #
        xi, yi = points[:, 0] - X0[0], points[:, 1] - X0[1]

        # NOTE:
        # -- The available rbf-function is only defined for a quarter of one shell.
        # in order to get z and t values for an entire shell the abs-function 
        # is used. The coordinate system of the quarter must be defined in the 
        # lower left corner; the coordinate systemn of the entire one shell must 
        # be defined in the center of the shell so that the coordinate system
        # for the master quarter remains unchanged. 
        # -- The transformation is performed based on the defined class attributes
        # of hp_shell_stb: length_x, length_y, length_z, delta_h, delta_h_scalefactor
        # characterizing the properties of the master quarter 

        # number of local grid points for each coordinate direction
        # values must range between 0 and 1
        #         
        points_tilde_list = []
        for i_row in range( points.shape[0] ):
            # get the x, y coordinate pair defined in the input 
            # file in global coordinates
            #
            x = xi[i_row]
            y = yi[i_row]

            # transform values to local coordinate system, 
            # i.e. move point to the 'master roof' containing the 
            # global coordinate system:
            #
            if x <= self.length_x and y <= self.length_y:
                # point lays in first (master) roof 
                #
                x_tilde = x
                y_tilde = y

            elif x >= self.length_x and y <= self.length_y:
                # point lays in second roof:
                #
                # roof length = 2* length of the master quarter
                # (e.g. 2*4,0m = 8,00m for obj-file "4x4m")
                x_tilde = x - 2 * self.length_x
                y_tilde = y

            elif x <= self.length_x and y >= self.length_y:
                # point lays in third roof:
                #
                x_tilde = x
                y_tilde = y - 2 * self.length_y

            elif x >= self.length_x and y >= self.length_y:
                # point lays in fourth roof:
                #
                x_tilde = x - 2 * self.length_x
                y_tilde = y - 2 * self.length_y

            points_tilde_list.append( [x_tilde, y_tilde] )

        points_tilde_arr = array( points_tilde_list, dtype = 'float_' )
        xi_tilde = points_tilde_arr[:, 0]
        yi_tilde = points_tilde_arr[:, 1]
#        print 'points_tilde_arr', points_tilde_arr


        xi_ = abs( xi_tilde ) / self.length_x
        yi_ = abs( yi_tilde ) / self.length_y

        #-----------------------------------------------
        # get the normalized rbf-function for the upper 
        # and lower face of the master quarter
        #-----------------------------------------------
        # NOTE: the underline character indicates a normalized value 

        # normalized coordinates of the vertices for lower- and upper-face
        #
        vl_arr_, vu_arr_ = normalize_rsurfaces( self.vl_arr,
                                                self.vu_arr )

        # use a radial basis function approximation (rbf) (i.e. interpolation of
        # scattered data) based on the normalized vertex points of the lower face
        #
        x_ = vl_arr_[:, 0]
        y_ = vl_arr_[:, 1]

        if self.geo_input_name == '350x350cm':
            x_ = 1 - vl_arr_[:, 0]

        z_l_ = vl_arr_[:, 2]
        rbf_l = Rbf( x_, y_, z_l_, function = 'cubic' )

        # get the z-value at the supplied local grid points
        # of the lower face
        #
        zi_lower_ = rbf_l( xi_, yi_ )

        # use a radial basis function approximation (rbf) (i.e. interpolation of
        # scattered data) based on the normalized vertex points of the upper face
        #
        x_ = vu_arr_[:, 0]
        y_ = vu_arr_[:, 1]

        if self.geo_input_name == '350x350cm':
            x_ = 1 - vu_arr_[:, 0]

        z_u_ = vu_arr_[:, 2]
        rbf_u = Rbf( x_, y_, z_u_ , function = 'cubic' )

        # get the z-value at the supplied local grid points
        # of the upper face 
        #
        zi_upper_ = rbf_u( xi_, yi_ )

        # approach of the slope to get thickness perpendicular to slope
        #
        # thickness is multiplied by the supplied zi coordinate
        # and z value of mid plane
        #
        t_ = zi_upper_ - zi_lower_

        z_middle_ = ( zi_lower_ + ( zi_upper_ - zi_lower_ ) * 0.5 / self.delta_h_scalefactor ) * self.delta_h_scalefactor

        if perpendicular_t == True:
            # delta shift of x and y for estimation of slope will be done in 4 direction
            # 0, 45, 90 and 135 degrees 
            print "--- perpendicular ---"
            delta = 0.000001

            # shift in x

            dz_x_p_ = ( rbf_u( xi_ + delta, yi_ ) + rbf_l( xi_ + delta, yi_ ) ) / 2.0
            dz_x_m_ = ( rbf_u( xi_ - delta, yi_ ) + rbf_l( xi_ - delta, yi_ ) ) / 2.0

            slope_x_ = ( dz_x_p_ - dz_x_m_ ) / ( 2.0 * delta )
            angle_x = arctan( slope_x_ * self.length_z / self.length_x )
            f_1 = cos( angle_x )

            # shift in y 

            dz_y_p_ = ( rbf_u( xi_, yi_ + delta ) + rbf_l( xi_, yi_ + delta ) ) / 2.0
            dz_y_m_ = ( rbf_u( xi_, yi_ - delta ) + rbf_l( xi_, yi_ - delta ) ) / 2.0

            slope_y_ = ( dz_y_p_ - dz_y_m_ ) / ( 2.0 * delta )
            angle_y = arctan( slope_y_ * self.length_z / self.length_x )
            f_2 = cos( angle_y )

            #shift +x +y; -x -y

            dz_x_p_y_p_ = ( rbf_u( xi_ + delta, yi_ + delta ) + rbf_l( xi_ + delta, yi_ + delta ) ) / 2.0
            dz_x_m_y_m_ = ( rbf_u( xi_ - delta, yi_ - delta ) + rbf_l( xi_ - delta, yi_ - delta ) ) / 2.0

            slope_x_p_y_p_ = ( dz_x_p_y_p_ - dz_x_m_y_m_ ) / ( 2.0 * sqrt( 2 ) * delta )
            angle_x_p_y_p = arctan( slope_x_p_y_p_ * self.length_z /
                                       ( self.length_x ** 2 + self.length_y ** 2 ) ** 0.5 )
            f_3 = cos( angle_x_p_y_p )

            # shift in +x,-y ; -x and +y

            dz_x_p_y_m_ = ( rbf_u( xi_ + delta, yi_ - delta ) + rbf_l( xi_ + delta, yi_ - delta ) ) / 2.0
            dz_x_m_y_p_ = ( rbf_u( xi_ - delta, yi_ + delta ) + rbf_l( xi_ - delta, yi_ + delta ) ) / 2.0

            slope_x_p_y_m_ = ( dz_x_p_y_m_ - dz_x_m_y_p_ ) / ( sqrt( 2 ) * 2.0 * delta )
            angle_x_p_y_m = arctan( slope_x_p_y_m_ * self.length_z /
                                   ( self.length_x ** 2 + self.length_y ** 2 ) ** 0.5 )
            f_4 = cos( angle_x_p_y_m )

            # obtain minimum factor for good estimate of maximum slope

            factor = min( [f_1, f_2, f_3, f_4], axis = 0 )
            t_ = t_ * factor

        return  xi, yi, z_middle_ * self.length_z, t_ * self.length_z

    def _read_thickness_data( self, file_name ):
        '''to read the stb - X and Y coordinates ( m ) save the xls - worksheet
        to a csv - file using ';' as filed delimiter and ' ' ( blank )
        as text delimiter.
        Stb Data needs to have same range of values in X and Y direction and same unit [m],
        as defined as length_x and length_y
        '''
        print '*** reading thickness data from file: ', file_name, ' ***'

        # get the column headings defined in the second row 
        # of the csv thickness input file
        # "Nr.;X;Y;Z;[mm]"
        # 
        file = open( file_name, 'r' )
        lines = file.readlines()
        column_headings = array( lines[1].split( ';' ) )
        elem_no_idx = where( 'Nr.' == column_headings )[0]
        X_idx = where( 'X' == column_headings )[0]
        Y_idx = where( 'Y' == column_headings )[0]
        Z_idx = where( 'Z' == column_headings )[0]
        thickness_idx = where( '[mm]\n' == column_headings )[0]

        input_arr = loadtxt( file_name, delimiter = ';', skiprows = 2 )

        # elem number:
        #
        elem_no = input_arr[:, elem_no_idx]

        # coordinates [m]:
        #
        X = input_arr[:, X_idx][:, 0]
        Y = input_arr[:, Y_idx][:, 0]

#        print 'thickness_idx', thickness_idx
        if thickness_idx != []:
            thickness_stb = input_arr[:, thickness_idx][:, 0] / 1000.
            return elem_no, X, Y, thickness_stb
        else:
            thickness_stb = ones_like( elem_no )
            return elem_no, X, Y, thickness_stb


    def _read_elem_coords( self, file_name ):
        '''x,y -coordinates must be read from old file
        '''
        input_arr = loadtxt( file_name, delimiter = ';', skiprows = 2 )

        elem_no = input_arr[:, 0]
        X = input_arr[:, 2]
        Y = input_arr[:, 3]

        return elem_no, X, Y

    def _read_nodal_coords( self, file_name ):
        '''read the nodal coordinates of the mid - surface
        defined in a csv - file. To export the excel sheet
        to csv use ";" as a field delimiter and "" ( none )
        as a text delimiter.
        Note that some lines do not contain values !
        '''
        print '*** reading nodal coordinates from file: ', file_name, ' ***'

        file = open( file_name, 'r' )

        # read the column headings (first two lines)
        #
        first_line = file.readline()
        second_line = file.readline()
        column_headings = second_line.split( ';' )
        # remove '\n' from last string element in list
        column_headings[-1] = column_headings[-1][:-1]
        column_headings_arr = array( column_headings )

        # check in which column the node number and the 
        # carthesian coordinates can be found 
        #
        elem_no_idx = where( 'Nr.' == column_headings_arr )[0]
        X_idx = where( 'X [m]' == column_headings_arr )[0]
        Y_idx = where( 'Y [m]' == column_headings_arr )[0]
        Z_idx = where( 'Z [m]' == column_headings_arr )[0]

        lines = file.readlines()

        lines_list = [ line.split( ';' ) for line in lines ]

        empty_lines_idx = []
        ll = []
        for i_line, line in enumerate( lines_list ):

            # check if line contains values or only a node number!
            #
            if line[1] == 'Standard':
                ll.append( [line[elem_no_idx], line[X_idx], line[Y_idx], line[Z_idx]] )
            else:
                # NOTE: current number in file starts with 1, index in loop starts with 0
                # therefore add 1 in the index list
                #
                empty_lines_idx.append( i_line + 1 )

        input_arr = array( ll, dtype = 'float_' )

        node_no = input_arr[:, 0]
        X = input_arr[:, 1]
        Y = input_arr[:, 2]
        Z = input_arr[:, 2]

        return node_no, X, Y, Z, empty_lines_idx


    def compare_thickness_values( self, thickness, thickness_stb ):
        '''get relative difference between the calucated thickness
        read in from the obj file, cut of and projected with respect to
        the approximated data given from stb.
        '''
        thickness = thickness.reshape( shape( thickness_stb ) )
        error = abs( 1 - thickness / thickness_stb ) * 100
        return error


    def export_midsurface_data( self, node_no, x, y, z_middle, file_name, empty_lines_idx ):
        '''exports data to csv - worksheet
        '''
        print '*** writing middle surface data to file,', file_name, ' ***'

        data = c_[node_no, x, y, z_middle]
        file = open( file_name, 'w' )
        writer = csv.writer( file, delimiter = ";", lineterminator = "\n" )
        writer.writerow( ['node_number', 'x[m]', 'y[m]', 'z[m]'] )
        writer.writerows( data )

        file = file.close()

        # if file contains empty lines add them at the positions 
        # defined in 'empty_lines_idx'
        #
        if len( empty_lines_idx ) != 0:

            print '--- file contains ', len( empty_lines_idx ), ' empty_lines ---'

            # file without empty lines
            #
            file = open( file_name, 'r' )
            lines = file.readlines()

            # overwrite file including empty lines
            #
            file = open( file_name, 'w' )

            # index 'n' runs in the array without empty lines
            # index 'i' runs in the array with empty lines
            #
            n = 0
            for i in range( data.shape[0] + len( empty_lines_idx ) ):

                if i in empty_lines_idx:
                    file.writelines( str( i ) + ";;;;\n" )
                else:
                    file.writelines( lines[n] )
                    n += 1

            # add last line:
            #
            file.writelines( lines[-1 ] )

            file.close()
            print '--- empty lines added to file ---'

        return


    def export_thickness_data( self, elem_no, x, y, t, file_name ):
        '''exports data to csv - worksheet
        '''
        print '*** writing thickness data to file,', file_name, ' ***'

        data = c_[elem_no, x, y, t * 1000]
        print shape( data )
        writer = csv.writer( open( file_name, 'w' ), delimiter = ";", lineterminator = "\n" )
        writer.writerow( ['element_number', 'x[m]', 'y[m]', 't[mm]'] )
        writer.writerows( data )
        return

    @show
    def show( self, x, y, z_middle, displayed_value ):
        """Test contour_surf on regularly spaced co-ordinates like MayaVi.
        """
        print '*** plotting data***'
        s = points3d( X, Y, z_middle, displayed_value, colormap = "gist_rainbow", mode = "cube", scale_factor = 0.3 )

        sb = colorbar( s )
        # Recorded script from Mayavi2
        #try:
        #    engine = mayavi.engine
        #except NameError:
        #    from enthought.mayavi.api import Engine
        #    engine = Engine()
        #    engine.start()
        #if len(engine.scenes) == 0:
        #    engine.new_scene()
        # ------------------------------------------- 
        glyph = s#.pipeline.scenes[0].children[0].children[0].children[0]
        glyph.glyph.glyph_source.glyph_source.center = array( [ 0., 0., 0.] )
        glyph.glyph.glyph_source.glyph_source.progress = 1.0
        glyph.glyph.glyph_source.glyph_source.x_length = 0.6
        glyph.glyph.glyph_source.glyph_source.y_length = 0.6
        sb.scalar_bar.title = 'thickness [m]'
        #print s.pipeline
        #s.scene.background = (1.0, 1.0, 1.0)

        return s

if __name__ == '__main__':

    #------------------------------------------------------------------
    # NUMBER OF SHELLS
    #------------------------------------------------------------------
    #
    n_shells = 1

    # if delta_h = 1,30 m set scalefactor to  1.3. 
    #
    delta_h_scalefactor = 1.0

    # tags for automatic file name generation
    #
    if n_shells == 1:
        n_shells_str = '1shell'
    elif n_shells == 4:
        n_shells_str = '4shells'

    delta_h_scalefactor_str = 'delta_h_' + str( int( 1000 * delta_h_scalefactor ) ) + 'mm'

    # for 7x7m geometry:

    #
    n_shells_str = '1shell_7x7m'
    delta_h_scalefactor_str = 'delta_h_865mm'

    # NOTE:
    # Leave 'length_z' and 'delta_h' as they are. 
    # They correspond in both cases (i.e. delta_h = 1,00m or 1,30m) 
    # to the values specified in the original obj-file ('4x4m')
    #
    hp = HPShell( cut_off_lowerface = True,
                  delta_h_scalefactor = delta_h_scalefactor )

    print '*** INPUT ***'
    print '--- number of shells = ', n_shells, ' ---'
    print '--- delta h scale factor = ', delta_h_scalefactor, ' ---'
    print '---', delta_h_scalefactor_str, ' ---\n'

    #------------------------------------------------------------------
    # THICKNESS
    # get a csv file with the exact thickness of the shell elements
    # in the order of the supplied csv coordinates (x,y)
    # (NOTE: the thickness coordinates are defined in the middle of the elements)  
    # --> t_elem in [mm]
    #------------------------------------------------------------------
    print '*** CALCULATE THICKNESS *** '

    # 1) read x,y coord values from csv-file used to specify shell thickness  
    #
    input_data_thickness = 'RFEM_file_export/input_data_elem_coords_' + n_shells_str + '.csv'
    file_name = input_data_thickness

    if n_shells == 1:
        elem_no, X, Y, thickness_stb = hp._read_thickness_data( file_name )

    elif n_shells == 4:
        # element coordinates at taken from a state data output file given from stb.
        # This order has not been changed since an is used for the generation of 
        # the thickness output file. To read the file a separate read procedure has been defined.
        #
        elem_no, X, Y = hp._read_elem_coords( file_name )

    # 2) get thickness data calculated using 'get_thickness_and_middlesurface' 
    #
    xi , yi , z_middle , t = hp.get_mid_surface_and_thickness( c_[X, Y], perpendicular_t = True )
    print "thickness at corner must evaluate to 62 mm (almost no slope): ", t[ argmax( z_middle ) ]
    print "minimum thickness must evaluate to about 60 mm (if perpendicular = True)", min( t )

    # 3) export thickness to file
    #
    ouput_data_thickness = 'RFEM_file_export/output_data_thickness_' + n_shells_str + '_' + delta_h_scalefactor_str + '.csv'
    hp.export_thickness_data( elem_no, X, Y, t, ouput_data_thickness )

    # 4) plot midsurface and thickness in mayavi
    #
    hp.show( X, Y, z_middle, t )

    ### possible further step show error
    #
#    error = hp.compare_thickness_values( t, thickness_stb )
#    print " - - ->difference between stb - thickness ( approximation ) as read from file and calculated thickness ( exact ) ", error
#    hp.show( X, Y, z_middle, error )
#
    #------------------------------------------------------------------
    # MIDSURFACE
    # get a csv file with the nodal coords of the shell elements
    # in the order of the supplied csv coordinates (x,y) 
    # --> z_node in [m]
    #------------------------------------------------------------------
    print '\n'
    print '*** CALCULATE MIDSURFACE *** '

    # 1) read x,y coords from csv-file used to specify shell middle surface
    # (NOTE: the midsurface coordinates are defined at the nodes)  
    #
    file_name = 'RFEM_file_export/input_data_nodal_coords_' + n_shells_str + '.csv'
    node_no, X, Y, Z, empty_lines_idx = hp._read_nodal_coords( file_name )

    # 2) get midsurface data calculated using 'get_thickness_and_middlesurface' 
    #
    xi , yi , z_middle , t = hp.get_mid_surface_and_thickness( c_[X, Y], perpendicular_t = True )
    print "position and value of lowest shell thickness: ", t[argmin( t )], ' at node ', argmin( t )

    # shift coordinate system in z-direction to start with node_no = 0 
    # at position (x,y,z) = (0,0,0). Change the sign in order to point with 
    # the z-coordinate down.  

    print ' - - -z - coordinates before shifting:- - -'
    print 'z_middle[0] ( z value in origin: must evaluate to 31cm / 2 = 15, 5 cm; origin at lowerface cut off position )', z_middle[0]
    print 'z_middle[-1] ( edge value ) must evaluate to 1, 03 m = delta h + 6cm / 2', z_middle[-1]
    z_middle -= z_middle[0]
    z_middle *= -1

    # OPTIONAL:
#    # move the coordinate system in the support of the column if desired
#    #
#    if delta_h_scalefactor == 1.0:
#        z_middle -= 3.15
#        delta_h_scalefactor_str = delta_h_scalefactor_str + '_315'
#    elif delta_h_scalefactor == 1.3:
#        z_middle -= 2.85
#        delta_h_scalefactor_str = delta_h_scalefactor_str + '_285'

    print ' - - -z - coordinates after shifting:- - -'
    print 'z_middle[0] ( must evaluate to 0 )', z_middle[0]
    print 'z_middle[-1] ( must evaluate to 100 cm - 31cm / 2 + 6cm / 2 = 87, 5 cm )', z_middle[-1]

    # 3) export middle surface to file
    #    
    output_data_midsurface = 'RFEM_file_export/output_data_midsurface_' + n_shells_str + '_' + delta_h_scalefactor_str + '.csv'
    hp.export_midsurface_data( node_no, X, Y, z_middle, output_data_midsurface, empty_lines_idx )

    # 4) plot midsurface and thickness in mayavi
    #
    hp.show( X, Y, z_middle, t )



