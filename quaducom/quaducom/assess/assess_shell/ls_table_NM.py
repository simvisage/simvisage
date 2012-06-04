'''
Created on Jun 23, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, Trait

from etsproxy.util.home_directory import \
    get_home_directory

from etsproxy.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from etsproxy.mayavi import \
    mlab

from matplotlib import pyplot

#from etsproxy.mayavi.mlab import \
#    colorbar, show, points3d
#
#from etsproxy.mayavi.api import \
#    Engine


from etsproxy.traits.ui.table_column import \
    ObjectColumn

from etsproxy.traits.ui.menu import \
    OKButton, CancelButton

from etsproxy.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, frompyfunc, min, shape, \
                max as ndmax, min as ndmin

from math import pi
from string import split
import os

from promod.simdb import \
    SimDB

import os
import pickle
import string
from os.path import join

# Access to the top level directory of the database
#
simdb = SimDB()


class LSArrayAdapter ( TabularAdapter ):

    columns = Property
    def _get_columns( self ):
#        print 'GETTING COLUMNS', self.object.columns, self.object, self.object.__class__
        columns = self.object.columns
        return [ ( name, idx ) for idx, name in enumerate( columns ) ]

    font = 'Courier 10'
    alignment = 'right'
    format = '%5.3f'#'%g'
    even_bg_color = Color( 0xE0E0FF )
    width = Float( 80 )

class LS( HasTraits ):
    '''Limit state class
    '''

    # backward link to the info shell to access the
    # input data when calculating 
    # the limit-state-specific values
    #
    ls_table = WeakRef

    #-------------------------------
    # ls columns
    #-------------------------------
    # defined in the subclasses
    #
    ls_columns = List

    #-------------------------------
    # sr columns
    #-------------------------------

    # stress resultant columns - for ULS this is defined in the subclasses
    #
    sr_columns = List( [] )

    #-------------------------------
    # geo columns form info shell
    #-------------------------------
    geo_columns = List( [] )

    #-------------------------------
    # state columns form info shell
    #-------------------------------
    state_columns = List( ['Px', 'Py', 'Pz', 'Mx', 'My', 'Mz',
                           'Px_', 'Py_', 'Pz', 'Mx_', 'My_', 'Mz',
                           'Mx_tot_', 'My_tot_', 'ex_tot_', 'ey_tot_', 'exy_tot_'] )

    # global coords:
    #
    Px = Property( Float )
    def _get_Px( self ):
        return self.ls_table.Px

    Py = Property( Float )
    def _get_Py( self ):
        return self.ls_table.Py

    Pz = Property( Float )
    def _get_Pz( self ):
        return self.ls_table.Pz

    Mx = Property( Float )
    def _get_Mx( self ):
        return self.ls_table.Mx

    My = Property( Float )
    def _get_My( self ):
        return self.ls_table.My

    Mz = Property( Float )
    def _get_Mz( self ):
        return self.ls_table.Mz

    # local coords:
    #
    Px_ = Property( Float )
    def _get_Px_( self ):
        return self.ls_table.Px_

    Py_ = Property( Float )
    def _get_Py_( self ):
        return self.ls_table.Py_

    Mx_ = Property( Float )
    def _get_Mx_( self ):
        return self.ls_table.Mx_

    My_ = Property( Float )
    def _get_My_( self ):
        return self.ls_table.My_

    # total moments end excentricities
    #
    Mx_tot_ = Property( Float )
    def _get_Mx_tot_( self ):
        return self.ls_table.Mx_tot_

    My_tot_ = Property( Float )
    def _get_My_tot_( self ):
        return self.ls_table.My_tot_

    ex_tot_ = Property( Float )
    def _get_ex_tot_( self ):
        return self.ls_table.ex_tot_

    ey_tot_ = Property( Float )
    def _get_ey_tot_( self ):
        return self.ls_table.ey_tot_

    exy_tot_ = Property( Float )
    def _get_exy_tot_( self ):
        return self.ls_table.exy_tot_

    #-------------------------------
    # ls table
    #-------------------------------

    # all columns associated with the limit state including the corresponding
    # stress resultants
    #
    columns = Property( List )
    @cached_property
    def _get_columns( self ):
        columns = self.geo_columns \
                + self.state_columns \
                + self.sr_columns \
                + self.ls_columns

        return columns


    # select column used for sorting the data in selected sorting order 
    #
    sort_column = Enum( values = 'columns' )
    def _sort_column_default( self ):
        return self.columns[-1]

    sort_order = Enum( 'descending', 'ascending', 'unsorted' )

    #-------------------------------------------------------
    # get the maximum value of the selected variable 
    # 'max_in_column' of the current sheet (only one sheet)
    #-------------------------------------------------------

    # get the maximum value of the chosen column
    #
    max_in_column = Enum( values = 'columns' )
    def _max_in_column_default( self ):
        return self.columns[-1]

    max_value = Property( depends_on = 'max_in_column' )
    def _get_max_value( self ):
        col = getattr( self, self.max_in_column )[:, 0]
        return max( col )

    #-------------------------------------------------------
    # get ls_table for View
    #-------------------------------------------------------

    # stack columns together for table used by TabularEditor
    #
    ls_array = Property( Array, depends_on = 'sort_column, sort_order' )
    @cached_property
    def _get_ls_array( self ):

        arr_list = [ getattr( self, col ) for col in self.columns ]

        # get the array currently selected by the sort_column enumeration
        #
        sort_arr = getattr( self, self.sort_column )[:, 0]
        sort_idx = argsort( sort_arr )
        ls_array = hstack( arr_list )

        if self.sort_order == 'descending':
            return ls_array[ sort_idx[::-1] ]
        if self.sort_order == 'ascending':
            return ls_array[ sort_idx ]
        if self.sort_order == 'unsorted':
            return ls_array


    # name of the trait that is used to assess the evaluated design
    #
    assess_name = Str( '' )

    #-------------------------------
    # ls group
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    ls_group = VGroup( 
                        HGroup( #Item( 'assess_name' ),
                                Item( 'max_in_column' ),
                                Item( 'max_value', style = 'readonly', format_str = '%6.2f' ),
                              ),
                        HGroup( Item( 'sort_column' ),
                                Item( 'sort_order' ),
                              ),
                     )


class ULS( LS ):
    '''Ultimate limit state
    '''

    #--------------------------------------------------------
    # ULS: material parameters (Inputs)
    #--------------------------------------------------------
    # tensile strength of the textile reinforcement [kN/m]
#    f_Rdtex = Float( 18.0, input = True )

    # ------------------------------------------------------------
    # ULS: outputs
    # ------------------------------------------------------------

    ls_columns = List( [] )

#    assess_name = 'max_My_tot_'
#    max_My_tot_ = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_My_tot_( self ):
#        return ndmax( abs( self.My_tot_) )

#    assess_name = 'max_ey_tot_'
#    max_ey_tot_ = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_ey_tot_( self ):
#        return ndmax( self.ey_tot_ )

#    assess_name = 'max_Pz'
#    max_Pz = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_Pz( self ):
#        return ndmax( self.Pz )

#    assess_name = 'max_ex_tot_'
#    max_ex_tot_ = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_ex_tot_( self ):
#        return ndmax( self.ex_tot_ )

#    assess_name = 'max_Mx_tot_'
#    max_Mx_tot_ = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_Mx_tot_( self ):
#        return ndmax( abs( self.Mx_tot_ ) )

    # essential for geotechnical stability; fullfilled if smaller then 1/9
    # (= gap smaller then center of foundation)
    # 
    assess_name = 'max_exy_tot_'
    max_exy_tot_ = Property( depends_on = '+input' )
    @cached_property
    def _get_max_exy_tot_( self ):
        return ndmax( self.exy_tot_ )


    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    traits_view = View( 
                       VGroup( 
                        HGroup( 
                            VGroup( 
                                Item( name = 'max_My', label = 'max_My' ),
                                label = 'resistance'
                                  ),
                             ),

                        VGroup( 
                            Include( 'ls_group' ),
                            Item( 'ls_array', show_label = False,
                                  editor = TabularEditor( adapter = LSArrayAdapter() ) )
                              ),
                            ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )



LSLIST = [ ULS]

class LSTable( HasTraits ):
    '''Assessment tool
    '''

    is_id = Int( 0 )

#    elem_no = Property( Array )
#    def _get_elem_no( self ):
#        return self.state_data['elem_no']

    # state data: stress resultants 
    # 
    state_data = Dict
    geo_data = Dict

    # global
    #
    Px = Property( Array )
    def _get_Px( self ):
        return self.state_data['Px']

    Py = Property( Array )
    def _get_Py( self ):
        return self.state_data['Py']

    Pz = Property( Array )
    def _get_Pz( self ):
        return self.state_data['Pz']

    Mx = Property( Array )
    def _get_Mx( self ):
        return self.state_data['Mx']

    My = Property( Array )
    def _get_My( self ):
        return self.state_data['My']

    Mz = Property( Array )
    def _get_Mz( self ):
        return self.state_data['Mz']

    # local
    #
    Px_ = Property( Array )
    def _get_Px_( self ):
        return sqrt( 2 ) / 2 * ( self.Py + self.Px )

    Py_ = Property( Array )
    def _get_Py_( self ):
        return sqrt( 2 ) / 2 * ( self.Py - self.Px )

    Mx_ = Property( Array )
    def _get_Mx_( self ):
        return sqrt( 2 ) / 2 * ( self.My + self.Mx )

    My_ = Property( Array )
    def _get_My_( self ):
        return sqrt( 2 ) / 2 * ( self.My - self.Mx )

    # total moments and excentricities at 
    # bottom side of the foundation (used for 
    # geotechnical stability "Kippnachweis (2. Kenrweite)" )
    #
    Mx_tot_ = Property( Array )
    def _get_Mx_tot_( self ):
        # positive Py_ creates a positive (!) My 
        # (height foundation = 0.60 m)
        #
        return self.Mx_ + 0.60 * self.Py_

    My_tot_ = Property( Array )
    def _get_My_tot_( self ):
        # positive Px_ creates a negative (!) My 
        # (height foundation = 0.60 m)
        #
        return self.My_ - 0.60 * self.Px_

    ex_tot_ = Property( Array )
    def _get_ex_tot_( self ):
        # positive values for Pz are compressive forces
        # therefore add (!) own weight foundation (=100 kN)
        #
        return abs( self.My_tot_ / ( self.Pz + 100. ) )

    ey_tot_ = Property( Array )
    def _get_ey_tot_( self ):
        # positive values for Pz are compressive forces
        # therefore add (!) own weight foundation
        #
        return abs( self.Mx_tot_ / ( self.Pz + 100. ) )

    exy_tot_ = Property( Array )
    def _get_exy_tot_( self ):
        return ( self.ex_tot_ / 2.5 ) ** 2 + ( self.ey_tot_ / 2.5 ) ** 2


    #------------------------------------------
    # combinations of limit states, stress resultants and directions
    #------------------------------------------

    ls = Trait( 'ULS',
                {'ULS' : ULS} )

    ls_class = Instance( LS )
    def _ls_class_default( self ):
        '''ls instances, e.g. ULS()
        '''
        ls_class = self.ls_
        return ls_class( ls_table = self )


    #------------------------------------------
    # get arrays for the TabularEditor:
    #------------------------------------------

#    ls_class = Property( Instance( LS ) )
#    def _get_ls_class( self ):
#        ls_class = self.ls_
#        return ls_class( ls_table = self )

    assess_value = Property
    def _get_assess_value( self ):
        ls = self.ls_class
        return getattr( ls, ls.assess_name )

    # ------------------------------------------------------------
    # View 
    # ------------------------------------------------------------

    traits_view = View( Tabbed( 
                            Item( 'ls_class@' , label = "ls", show_label = False ),
                            scrollable = False,
                         ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )

if __name__ == '__main__':

    ls = LS()
    print ls.ls_table
