#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Jul 1, 2010 by: rch

from etsproxy.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool

from etsproxy.util.home_directory import \
    get_home_directory

from etsproxy.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from etsproxy.mayavi import \
    mlab

from etsproxy.traits.ui.table_column import \
    ObjectColumn

from etsproxy.traits.ui.menu import \
    OKButton, CancelButton

from etsproxy.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, sqrt, frompyfunc

from math import pi
from string import split
import os

#from apps.projects.sfb532demo.assess_shell.ls_table import LSTable
from quaducom.assess.assess_shell.ls_table import LSTable

class InfoShellMngr( HasTraits ):
    '''Assessment tool
    '''
    #------------------------------------------
    # specify default data input files:
    #------------------------------------------

    # raw input file for thicknesses (and coordinates)
    #
    geo_data_file = Str
    def _geo_data_file_default( self ):
        return 'input_geo_data.csv'

    # raw input file for element numbers, coordinates, and stress_resultants
    #
    state_data_file = Str
    def _state_data_file_default( self ):
        return 'input_state_data.csv'

    info_shell_uls = Property( Instance( LSTable ),
                               depends_on = 'state_data_file' )
    @cached_property
    def _get_info_shell_uls( self ):
        return LSTable( geo_data = self.geo_data,
                          state_data = self.state_data,
                          ls = 'ULS' )

    info_shell_sls = Property( Instance( LSTable ),
                               depends_on = 'state_data_file' )
    @cached_property
    def _get_info_shell_sls( self ):
        return LSTable( geo_data = self.geo_data,
                          state_data = self.state_data,
                          ls = 'SLS' )

    #------------------------------------------
    # read the geometry data from file 
    # (corrds and thickness):
    #------------------------------------------

    def _read_geo_data( self, file_name ):
        '''to read the stb-thickness save the xls-worksheet 
        to a csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
        # get the column headings defined in the second row 
        # of the csv thickness input file
        # "Nr.;X;Y;Z;[mm]"
        #
        file = open( file_name, 'r' )
        first_line = file.readline()
        second_line = file.readline()
        column_headings = second_line.split( ';' )
        # remove '\n' from last string element in list
        column_headings[-1] = column_headings[-1][:-1]
        column_headings_arr = array( column_headings )
        elem_no_idx = where( 'Nr.' == column_headings_arr )[0]
        X_idx = where( 'X' == column_headings_arr )[0]
        Y_idx = where( 'Y' == column_headings_arr )[0]
        Z_idx = where( 'Z' == column_headings_arr )[0]
        thickness_idx = where( '[mm]' == column_headings_arr )[0]
        # read the float data:
        #
        input_arr = loadtxt( file_name, delimiter = ';', skiprows = 2 )

        # element number:
        #
        elem_no = input_arr[:, elem_no_idx]

        # coordinates [m]:
        #
        X = input_arr[:, X_idx]
        Y = input_arr[:, Y_idx]
        Z = input_arr[:, Z_idx]

        # element thickness [mm]:
        #
        thickness = input_arr[:, thickness_idx]

        return  {'elem_no':elem_no,
                 'X':X, 'Y':Y, 'Z':Z,
                 'thickness':thickness }

    # coordinates and element thickness read from file:
    # 
    geo_data = Property( Dict, depends_on = 'geo_data_file' )
    @cached_property
    def _get_geo_data( self ):
        return self._read_geo_data( self.geo_data_file )

    #------------------------------------------
    # read the state data from file 
    # (elem_no, coords, moments and normal forces):
    #------------------------------------------

    def _Xread_state_data( self, file_name ):
        '''to read the stb-stress_resultants save the xls-worksheet 
        to a csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
        ######################################################
        # this method returns only the MAXIMUM VALUES!!!
        # @todo: dublicate the elem number and coordinates and add also the minimum values
        ######################################################

        # get the column headings defined in the second row 
        # of the csv soliciotations input file
        #
#        column_headings = array(["Nr.","Punkt","X","Y","Z","mx","my","mxy","vx","vy","nx","ny","nxy"])
        file = open( file_name, 'r' )
        lines = file.readlines()
        column_headings = lines[1].split( ';' )
        # remove '\n' from last string element in list
        column_headings[-1] = column_headings[-1][:-1]
        column_headings_arr = array( column_headings )
        elem_no_idx = where( 'Nr.' == column_headings_arr )[0]
        X_idx = where( 'X' == column_headings_arr )[0]
        Y_idx = where( 'Y' == column_headings_arr )[0]
        Z_idx = where( 'Z' == column_headings_arr )[0]
        mx_idx = where( 'mx' == column_headings_arr )[0]
        my_idx = where( 'my' == column_headings_arr )[0]
        mxy_idx = where( 'mxy' == column_headings_arr )[0]
        nx_idx = where( 'nx' == column_headings_arr )[0]
        ny_idx = where( 'ny' == column_headings_arr )[0]
        nxy_idx = where( 'nxy' == column_headings_arr )[0]

        # define arrays containing the information from the raw input file
        #

        # @todo: check how loadtxt can be used directly
        #        instead of reading lines line by line?
        # input_arr = loadtxt( file_name , delimiter=';', skiprows = 2 )

        # read max-values (=first row of each double-line):
        # read stress_resultant csv-input file line by line in steps of two
        # starting in the third line returning an data array.
        #
        n_columns = len( lines[2].split( ';' ) )
        # auxiliary line in array sliced off in input array
        data_array = zeros( n_columns )
        for n in range( 2, len( lines ), 2 ):
            line_split = lines[n].split( ';' )
            line_array = array( line_split, dtype = float )
            data_array = vstack( [ data_array, line_array ] )
        input_arr = data_array[1:, :]

        # element number:
        #
        elem_no = input_arr[:, elem_no_idx]

        # coordinates [m]:
        #
        X = input_arr[:, X_idx]
        Y = input_arr[:, Y_idx]
        Z = input_arr[:, Z_idx]

        # moments [kNm/m]
        #
        mx = input_arr[:, mx_idx]
        my = input_arr[:, my_idx]
        mxy = input_arr[:, mxy_idx]

        # normal forces [kN/m]:
        #
        nx = input_arr[:, nx_idx]
        ny = input_arr[:, ny_idx]
        nxy = input_arr[:, nxy_idx]

        return { 'elem_no':elem_no, 'X':X, 'Y':Y, 'Z':Z,
                 'mx':mx, 'my':my, 'mxy':mxy,
                 'nx':nx, 'ny':ny, 'nxy':nxy }


    def _read_state_data( self, file_name ):
        '''to read the stb-stress_resultants save the xls-worksheet 
        to a csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
        ######################################################
        # this method returns only the MAXIMUM VALUES!!!
        # @todo: dublicate the elem number and coordinates and add also the minimum values
        ######################################################

        # get the column headings defined in the second row 
        # of the csv solicitations input file
        #
#        column_headings = array(["Nr.","Punkt","X","Y","Z","mx","my","mxy","vx","vy","nx","ny","nxy"])
        file = open( file_name, 'r' )
        lines = file.readlines()
        column_headings = lines[1].split( ';' )
        # remove '\n' from last string element in list
        column_headings[-1] = column_headings[-1][:-1]
        column_headings_arr = array( column_headings )
        elem_no_idx = where( 'Nr.' == column_headings_arr )[0]
        X_idx = where( 'X' == column_headings_arr )[0]
        Y_idx = where( 'Y' == column_headings_arr )[0]
        Z_idx = where( 'Z' == column_headings_arr )[0]
        m1_idx = where( 'm1' == column_headings_arr )[0]
        m2_idx = where( 'm2' == column_headings_arr )[0]
        n1_idx = where( 'n1' == column_headings_arr )[0]
        n2_idx = where( 'n2' == column_headings_arr )[0]

        # define arrays containing the information from the raw input file
        #

        # @todo: check how loadtxt can be used directly
        #        instead of reading lines line by line?
        # input_arr = loadtxt( file_name , delimiter=';', skiprows = 2 )

        # read max-values (=first row of each double-line):
        # read stress_resultant csv-input file line by line in steps of two
        # starting in the third line returning an data array.
        #
        n_columns = len( lines[2].split( ';' ) )
        # auxiliary line in array sliced off in input array
        data_array = zeros( n_columns )
        for n in range( 2, len( lines ), 2 ):
            line_split = lines[n].split( ';' )
            line_array = array( line_split, dtype = float )
            data_array = vstack( [ data_array, line_array ] )
        input_arr = data_array[1:, :]

        # element number:
        #
        elem_no = input_arr[:, elem_no_idx]

        # coordinates [m]:
        #
        X = input_arr[:, X_idx]
        Y = input_arr[:, Y_idx]
        Z = input_arr[:, Z_idx]

        # moments [kNm/m]
        #
        mx = input_arr[:, m1_idx]
        my = input_arr[:, m2_idx]
        mxy = zeros_like( input_arr[:, m1_idx] )

        # normal forces [kN/m]:
        #
        nx = input_arr[:, n1_idx]
        ny = input_arr[:, n2_idx]
        nxy = zeros_like( input_arr[:, m1_idx] )

        return { 'elem_no':elem_no, 'X':X, 'Y':Y, 'Z':Z,
                 'mx':mx, 'my':my, 'mxy':mxy,
                 'nx':nx, 'ny':ny, 'nxy':nxy }

    # ------------------------------------------------------------
    # Read state data from file and assign attributes 
    # ------------------------------------------------------------

    # coordinates and stress_resultants read from file:
    # 
    state_data = Property( Dict, depends_on = 'state_data_file' )
    @cached_property
    def _get_state_data( self ):
        return self._read_state_data( self.state_data_file )


    # ------------------------------------------------------------
    # check input files for consistency
    # ------------------------------------------------------------

    @on_trait_change( 'geo_data_file, state_data_file' )
    def _check_input_files_for_consistency( self ):
        '''check if the element order of the thickness input file is 
        identical to the order in the stress_resultant input file
        '''
        if not all( self.geo_data['X'] ) == all( self.state_data['X'] ) or \
            not all( self.geo_data['Y'] ) == all( state_data['Y'] ) or \
            not all( self.geo_data['Z'] ) == all( state_data['Z'] ):
            raise ValueError, 'coordinates in file % s and file % s are not identical. Check input files for consistency ! ' \
                    % ( self.geo_data_file, self.state_data_file )
        else:
            print ' *** input files checked for consistency ( OK ) *** '
            return True


    # ------------------------------------------------------------
    # View 
    # ------------------------------------------------------------

    traits_view = View( VGroup( 
                              Item( 'state_data_file',
                              label = 'Evaluated input file for stress_resultants ',
                              style = 'readonly', emphasized = True ),

                        Item( 'geo_data_file',
                              label = 'Evaluated input file for thicknesses ',
                               style = 'readonly', emphasized = True ),
                        Tabbed( 
                        Group( 
                        Item( 'info_shell_sls@', show_label = False ),
                            label = 'SLS' ),
                        Group( 
                        Item( 'info_shell_uls@', show_label = False ),
                            label = 'ULS' ),
                            ),
                        ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )


if __name__ == '__main__':

    ifs = InfoShellMngr( state_data_file = 'input_state_data.csv',
                         geo_data_file = 'input_geo_data.csv' )

    print ifs.info_shell_uls.assess_value
#    print ifs.info_shell_sls.assess_value

    ifs.configure_traits()


#    print ifs.columns
#
#    ifs.selected_dir = 'y'
#    print ifs.columns
#
#    ifs.selected_sr = 'N'
#    print ifs.columns
#
#    ifs.selected_ls = ULS
#    print ifs.columns
#
#    print 'n1'
#    print ifs.n1
#    print ifs.my_M
#
#    print 'n_table'
#    print ifs.current_ls.ls_table

