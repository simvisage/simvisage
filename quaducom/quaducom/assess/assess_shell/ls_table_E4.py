'''
Created on Jun 23, 2010

@author: alexander
'''

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, Trait

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from enthought.mayavi import \
    mlab

from matplotlib import pyplot

#from enthought.mayavi.mlab import \
#    colorbar, show, points3d
#
#from enthought.mayavi.api import \
#    Engine


from enthought.traits.ui.table_column import \
    ObjectColumn

from enthought.traits.ui.menu import \
    OKButton, CancelButton

from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, frompyfunc, min, shape, \
                max as ndmax, min as ndmin

from math import pi
from string import split
import os

#from scipy.io import read_array

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
#    geo_columns = List( [ 'elem_no' ] )
    geo_columns = List( ['b_elem', 'h_elem'] )

    # width [m]        
    #
    b_elem = Property( Float )
    def _get_b_elem( self ):
        return self.ls_table.b_elem

    # heigth [m]
    #
    h_elem = Property( Float )
    def _get_h_elem( self ):
        return self.ls_table.h_elem

    # ------------------------------------------------------------
    # ULS - derived params:
    # ------------------------------------------------------------

    # area [m2]
    #
    A = Property( Float )
    def _get_A( self ):
        return self.h_elem * self.b_elem

    # moment of inertia (strong axis) [m3]
    #
    Wy = Property( Float )
    def _get_Wy( self ):
        return self.b_elem * self.h_elem ** 2 / 6.

    # moment of inertia (weak axis) [m3]
    #
    Wz = Property( Float )
    def _get_Wz( self ):
        return self.h_elem * self.b_elem ** 2 / 6.

    #-------------------------------
    # state columns form info shell
    #-------------------------------
    state_columns = List( ['N', 'Vy', 'Vz', 'MT', 'My', 'Mz', 'sig_N', 'sig_My', 'sig_Mz', ] )

    N = Property( Float )
    def _get_N( self ):
        return self.ls_table.N

    Vy = Property( Float )
    def _get_Vy( self ):
        return self.ls_table.Vy

    Vz = Property( Float )
    def _get_Vz( self ):
        return self.ls_table.Vz

    MT = Property( Float )
    def _get_MT( self ):
        return self.ls_table.MT

    My = Property( Float )
    def _get_My( self ):
        return self.ls_table.My

    Mz = Property( Float )
    def _get_Mz( self ):
        return self.ls_table.Mz

    #-------------------------------
    # derived state data - stresses
    #-------------------------------

    sig_N = Property( Float )
    def _get_sig_N( self ):
        return ( self.N / 1000. ) / self.A # convert [kN] to [MN]

    sig_My = Property( Float )
    def _get_sig_My( self ):
        return ( self.My / 1000. ) / self.Wy # convert [kNm] to [MNm]

    sig_Mz = Property( Float )
    def _get_sig_Mz( self ):
        return ( self.Mz / 1000. ) / self.Wz # convert [kNm] to [MNm]

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

    #---------------------------------
    # plot outputs in mlab-window 
    #---------------------------------

    plot_column = Enum( values = 'columns' )
    plot_3D = Button
    def _plot_3D_fired( self ):
        X = self.X_hf[:, 0]
        Y = self.Y_hf[:, 0]
        Z = self.Z_hf[:, 0]

        plot_col = getattr( self, self.plot_column )[:, 0]
        scale = 1 / max( plot_col )

        mlab.figure( figure = "SFB532Demo",
                     bgcolor = ( 1.0, 1.0, 1.0 ),
                     fgcolor = ( 0.0, 0.0, 0.0 ) )

        mlab.points3d( X, Y, ( -1.0 ) * Z, plot_col,
                       colormap = "copper",
                       mode = "cube",
                       scale_factor = scale )
        mlab.outline()

        mlab.scalarbar( title = self.plot_column, orientation = 'vertical' )

        mlab.show

        plot_column = Enum( values = 'columns' )



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
                                Item( 'plot_column' ),
                                Item( 'plot_3D' ),
                              ),
                     )


class SLS( LS ):
    '''Serviceability limit state
    '''

    #--------------------------------------------------------
    # SLS: material parameters (Inputs)
    #--------------------------------------------------------
    # tensile strength of the concrete [MPa]
    f_ctm = Float( 4.0, input = True )

    # ------------------------------------------------------------
    # SLS: outputs
    # ------------------------------------------------------------

    ls_columns = List( ['sig_tl', 'sig_tr', 'sig_bl', 'sig_br', ] )

    ls_values = Property( depends_on = '+input' )
    @cached_property
    def _get_ls_values( self ):
        '''get the outputs for SLS:
        stresses in the four corners of the bar
        '''
        sig_tl = self.sig_N - self.sig_My + self.sig_Mz
        sig_tr = self.sig_N - self.sig_My - self.sig_Mz
        sig_bl = self.sig_N + self.sig_My + self.sig_Mz
        sig_br = self.sig_N + self.sig_My - self.sig_Mz

        return { 'sig_tl':sig_tl, 'sig_tr':sig_tr,
                 'sig_bl':sig_bl, 'sig_br':sig_br }

    # stresses at corners of the bar
    #
    # (top-left)
    sig_tl = Property
    def _get_sig_tl( self ):
        return self.ls_values['sig_tl']

    # (top-right)
    sig_tr = Property
    def _get_sig_tr( self ):
        return self.ls_values['sig_tr']

    # (bottom-left)
    sig_bl = Property
    def _get_sig_bl( self ):
        return self.ls_values['sig_bl']

    # (bottom-right)
    sig_br = Property
    def _get_sig_br( self ):
        return self.ls_values['sig_br']

    assess_name = 'sig_ol'

    max_sig_tl = Property( depends_on = '+input' )
    @cached_property
    def _get_max_sig_tl( self ):
        return ndmax( self.sig_tl )


    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    traits_view = View( 
                       VGroup( 
                        HGroup( 
                            VGroup( 
                                Item( name = 'f_ctm', label = 'tensile strength [MPa]:  f_ctm ' ),
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



class ULS( LS ):
    '''Ultimate limit state
    '''

    #--------------------------------------------------------
    # ULS: material parameters (Inputs)
    #--------------------------------------------------------
    # tensile strength of the textile reinforcement [kN/m]
    f_Rdtex = Float( 18.0, input = True )

    # ------------------------------------------------------------
    # ULS: outputs
    # ------------------------------------------------------------

    ls_columns = List( ['sig_tot', 'N_equiv', 'n_equiv', 'n_tex', ] )

    ls_values = Property( depends_on = '+input' )
    @cached_property
    def _get_ls_values( self ):
        '''get the outputs for ULS
        '''
        sig_tot = self.sig_N + abs( self.sig_My ) + abs( self.sig_Mz )
        n_equiv = sig_tot * self.b_elem * 1000. # convert from [MN/m] to [kN/m] 
        N_equiv = sig_tot * self.A * 1000. # convert [MN] to [kN]
        n_tex = n_equiv / self.f_Rdtex

        return { 'sig_tot':sig_tot, 'N_equiv':N_equiv,
                 'n_equiv':n_equiv, 'n_tex':n_tex }

    sig_tot = Property
    def _get_sig_tot( self ):
        return self.ls_values['sig_tot']

    N_equiv = Property
    def _get_N_equiv( self ):
        return self.ls_values['N_equiv']

    n_equiv = Property
    def _get_n_equiv( self ):
        return self.ls_values['n_equiv']

    n_tex = Property
    def _get_n_tex( self ):
        return self.ls_values['n_tex']

    assess_name = 'max_n_tex'

    max_n_tex = Property( depends_on = '+input' )
    @cached_property
    def _get_max_n_tex( self ):
        return ndmax( self.n_tex )

#    assess_name = 'min_N'
#
#    min_N = Property( depends_on = '+input' )
#    @cached_property
#    def _get_min_N( self ):
#        return ndmin( self.N )



    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    traits_view = View( 
                       VGroup( 
                        HGroup( 
                            VGroup( 
                                Item( name = 'f_Rdtex', label = 'resistance per layer [kN/m]:  f_Rdtex ' ),
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



LSLIST = [ ULS, SLS ]

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

    N = Property( Array )
    def _get_N( self ):
        return self.state_data['N']

    Vy = Property( Array )
    def _get_Vy( self ):
        return self.state_data['Vy']

    Vz = Property( Array )
    def _get_Vz( self ):
        return self.state_data['Vz']

    MT = Property( Array )
    def _get_MT( self ):
        return self.state_data['MT']

    My = Property( Array )
    def _get_My( self ):
        return self.state_data['My']

    Mz = Property( Array )
    def _get_Mz( self ):
        return self.state_data['Mz']

    b_elem = Property( Array )
    def _get_b_elem( self ):
        return self.geo_data['b_elem']

    h_elem = Property( Array )
    def _get_h_elem( self ):
        return self.geo_data['h_elem']


    #------------------------------------------
    # combinations of limit states, stress resultants and directions
    #------------------------------------------

    ls = Trait( 'ULS',
                {'ULS' : ULS,
                 'SLS' : SLS } )

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
