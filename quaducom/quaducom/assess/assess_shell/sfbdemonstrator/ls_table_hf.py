'''
Created on Jun 23, 2010

@author: alexander
'''

from traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, Trait

from traitsui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from mayavi import \
    mlab

from matplotlib import pyplot

#from etsproxy.mayavi.mlab import \
#    colorbar, show, points3d
#
#from etsproxy.mayavi.api import \
#    Engine


from traitsui.table_column import \
    ObjectColumn

from traitsui.menu import \
    OKButton, CancelButton

from traitsui.tabular_adapter \
    import TabularAdapter

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, frompyfunc, min, shape, \
                max as ndmax

from math import pi
from string import split
import os

from matresdev.db.simdb import \
    SimDB
    
import os
import pickle
import string
from os.path import join

# Access to the top level directory of the database
#
simdb = SimDB()


class LSArrayAdapter (TabularAdapter):

    columns = Property
    def _get_columns(self):
#        print 'GETTING COLUMNS', self.object.columns, self.object, self.object.__class__
        columns = self.object.columns
        return [ (name, idx) for idx, name in enumerate(columns) ]

    font = 'Courier 10'
    alignment = 'right'
    format = '%5.2f'#'%g'
    even_bg_color = Color(0xE0E0FF)
    width = Float(80)


class LS( HasTraits ):
    '''Limit state class:
    Evaluation of the hinge forces linking the roof domains. 
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
    geo_columns = List( [] )



#    elem_no = Property( Float )
#        return self.ls_table.elem_no

    #-------------------------------
    # state columns form info shell
    #-------------------------------
    geo_columns = List( ['X_hf', 'Y_hf', 'Z_hf'] )
    state_columns = List( ['N_ip', 'V_ip', 'V_op'] )

    N_ip = Property( Float )
    def _get_N_ip( self ):
        return self.ls_table.N_ip

    V_ip = Property( Float )
    def _get_V_ip( self ):
        return self.ls_table.V_ip

    V_op = Property( Float )
    def _get_V_op( self ):
        return self.ls_table.V_op

    X_hf = Property( Array )
    def _get_X_hf( self ):
        return self.ls_table.X_hf

    Y_hf = Property( Array )
    def _get_Y_hf( self ):
        return self.ls_table.Y_hf

    Z_hf = Property( Array )
    def _get_Z_hf( self ):
        return self.ls_table.Z_hf
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

    hf_columns = Property( List )
    @cached_property
    def _get_hf_columns( self ):
        return ['sym', 'r0_r1']

    sym = Property
    def _get_sym( self ):
        return 'sym'

    r0_r1 = Property
    def _r0_r1( self ):
        return 'r0_r1'


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
#        if self.plot_column == 'n_tex':
#            plot_col = where( plot_col < 0, 0, plot_col )

        mlab.figure( figure = "SFB532Demo",
                     bgcolor = ( 1.0, 1.0, 1.0 ),
                     fgcolor = ( 0.0, 0.0, 0.0 ) )

        mlab.points3d( X, Y, ( -1.0 ) * Z, plot_col,
#                       colormap = "gist_rainbow",
#                       colormap = "Reds",
                       colormap = "copper",
                       mode = "cube",
                       scale_factor = scale )
        mlab.outline()

        mlab.scalarbar( title = self.plot_column, orientation = 'vertical' )

        mlab.show

        plot_column = Enum( values = 'columns' )


    plot_max_Value = Button
    def _plot_max_Value( self ):
        X = self.X_hf[:, 0]
        Y = self.Y_hf[:, 0]
        Z = self.Z_hf[:, 0]


        plot_col = getattr( self, self.plot_column )[:, 0]
        scale = 1 / max( plot_col )
#        if self.plot_column == 'n_tex':
#            plot_col = where( plot_col < 0, 0, plot_col )

        mlab.figure( figure = "SFB532Demo",
                     bgcolor = ( 1.0, 1.0, 1.0 ),
                     fgcolor = ( 0.0, 0.0, 0.0 ) )

        mlab.points3d( X, Y, ( -1.0 ) * Z, plot_col,
#                       colormap = "gist_rainbow",
#                       colormap = "Reds",
                       colormap = "copper",
                       mode = "cube",
                       scale_factor = scale )
        mlab.outline()

        mlab.scalarbar( title = self.plot_column, orientation = 'vertical' )

        mlab.show


    plot_hf_columns = Enum ( values = 'hf_columns' )
    plot_hf = Button
    def _plot_hf_fired( self ):
        X = self.X_hf[:, 0]
        Y = self.Y_hf[:, 0]
        Z = self.Z_hf[:, 0]
        pyplot.ioff()
        plot_col = self.plot_hf_columns
        if plot_col == 'sym':
            idx = where( Y == 8.0 )
            x_bar = X[idx]
            x_label = 'X [m]'

        if plot_col == 'r0_r1':
            idx = where( X == 8.0 )
            x_bar = Y[idx]
            x_label = 'Y [m]'

        N_ip = self.N_ip[idx]
        V_op = self.V_op[idx]

#        pyplot.clf()
        fig = pyplot.figure( facecolor = "white" )
        ax1 = fig.add_subplot( 1, 1, 1 )
        width = 0.05
        ax1.bar( x_bar - width / 2, N_ip, width, align = 'center', color = 'blue', label = '$N_{ip}$' )
        ax1.bar( x_bar + width / 2, V_op, width, align = 'center', color = 'red', label = '$V_{op}$' )
        ax1.set_xlabel( x_label, fontsize = 22 )
        ax1.set_ylabel( '$F_{int}$ [kN]', fontsize = 22 )
#        ax1.set_title( 'At axis' + plot_col )
        ax1.set_xlim( 0, max( x_bar ) )


        ax1.legend()
        pyplot.show()

#        fig.savefig(plot_col)
#        fig.savefig(facecolor = 'w', edgecolor = 'w',
#                       orientation = 'portrait')
        pyplot.clf()

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
                                Item( 'plot_hf_columns' ),
                                Item( 'plot_hf' ),
                              ),
                     )


class SLS( LS ):
    '''Serviceability limit state
    '''

    # ------------------------------------------------------------
    # SLS: material parameters (Inputs)
    # ------------------------------------------------------------

    # tensile strength [MPa]
    f_ctk = Float( 4.0, input = True )

    # flexural tensile strength [MPa]
    f_m = Float( 5.0, input = True )

    # ------------------------------------------------------------
    # SLS - derived params:
    # ------------------------------------------------------------

    # area
    #
    A = Property( Float )
    def _get_A( self ):
        return self.ls_table.D_elem * 1.

    # moment of inertia
    #
    W = Property( Float )
    def _get_W( self ):
        return 1. * self.ls_table.D_elem ** 2 / 6.

    # ------------------------------------------------------------
    # SLS: outputs
    # ------------------------------------------------------------

    ls_columns = List( ['sig_n', 'sig_m', 'eta_n', 'eta_m', 'eta_tot', ] )

    ls_values = Property( depends_on = '+input' )
    @cached_property
    def _get_ls_values( self ):
        '''get the outputs for SLS
        '''
        n = self.n
        m = self.m

        A = self.A
        W = self.W
        f_ctk = self.f_ctk
        f_m = self.f_m

        sig_n = n / A / 1000.
        sig_m = abs( m / W ) / 1000.
        eta_n = sig_n / f_ctk
        eta_m = sig_m / f_m
        eta_tot = eta_n + eta_m

        return { 'sig_n':sig_n, 'sig_m':sig_m,
                 'eta_n':eta_n, 'eta_m':eta_m,
                 'eta_tot':eta_tot }

    sig_n = Property
    def _get_sig_n( self ):
        return self.ls_values['sig_n']

    sig_m = Property
    def _get_sig_m( self ):
        return self.ls_values['sig_m']

    eta_n = Property
    def _get_eta_n( self ):
        return self.ls_values['eta_n']

    eta_m = Property
    def _get_eta_m( self ):
        return self.ls_values['eta_m']

    eta_tot = Property
    def _get_eta_tot( self ):
        return self.ls_values['eta_tot']

    assess_name = 'max_eta_tot'

#    assess_name = 'max_sig1_up'
#
#    max_sig1_up = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_sig1_up( self ):
#        return ndmax( self.sig1_up )

    # @todo: make it possible to select the assess value:
    #
#    assess_name = Enum( values = 'columns' )
#    def _assess_name_default( self ):
#        return self.columns[-1]

    max_eta_tot = Property( depends_on = '+input' )
    @cached_property
    def _get_max_eta_tot( self ):
        return ndmax( self.eta_tot )


    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    traits_view = View( VGroup( 
                            HGroup( Item( name = 'f_ctk', label = 'Tensile strength concrete [MPa]: f_ctk ' ),
                                    Item( name = 'f_m', label = 'Flexural tensile trength concrete [MPa]: f_m ' )
                                   ),
                            VGroup( 
                                Include( 'ls_group' ),

                                # @todo: currently LSArrayAdapter must be called both 
                                #        in SLS and ULS separately to configure columns 
                                #        arrangement individually
                                #
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
    Hinge_forces
    '''

    #--------------------------------------------------------
    # ULS: material parameters (Inputs)
    #--------------------------------------------------------

    # gamma-factor 
    gamma_N = Float( 1.0, input = True )
    gamma_V = Float( 1.0, input = True )

    # characteristic shear resistance 
    V_Rk = Float( 20.0, input = True )

    # characteristic pull-out resistance 
    N_Rk = Float( 11.3, input = True )

    # ------------------------------------------------------------
    # ULS - derived params:
    # ------------------------------------------------------------

    V_Rd = Property( Float )
    def _get_V_Rd( self ):
        return self.V_Rk / self.gamma_N

    N_Rd = Property( Float )
    def _get_N_Rd( self ):
        return self.N_Rk / self.gamma_V

    # ------------------------------------------------------------
    # ULS: outputs
    # ------------------------------------------------------------

    ls_columns = List( [] )

    sr_columns = [ 'eta_N', 'eta_V', 'eta_tot']

    ls_case = []

    conection_type = 'glued_bar'

    ls_values = Property( depends_on = '+input' )
    @cached_property
    def _get_ls_values( self ):
        '''get the outputs for ULS
        '''
        V_ed = self.N_ip / 2
        N_ed = self.V_op * 2
        print("shape", self.N_ip.shape)
#        N = sqrt( N_ip ** 2 + V_ip ** 2 )
#        N = N_ip
#        V = V_op


        # eta
        eta_N = abs( N_ed / self.N_Rd )
        eta_V = abs( V_ed / self.V_Rd )


        eta_tot = ( eta_V + eta_N ) / 1.2
        print("eta_tot", eta_tot.shape)
        eta_tot = where( eta_tot < eta_V, eta_V, eta_tot )
        eta_tot = where( eta_tot < eta_N, eta_N, eta_tot )


        return { 'eta_N':eta_N, 'eta_V':eta_V, 'eta_tot':eta_tot }

    eta_N = Property
    def _get_eta_N( self ):
        return self.ls_values['eta_N']

    eta_V = Property
    def _get_eta_V( self ):
        return self.ls_values['eta_V']

    eta_tot = Property
    def _get_eta_tot( self ):
        return self.ls_values['eta_tot']

    assess_name = 'max_eta_tot'

    max_eta_tot = Property( depends_on = '+input' )
    @cached_property
    def _get_max_eta_tot( self ):
        return ndmax( self.eta_tot )


    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    traits_view = View( 
                       VGroup( 
                        HGroup( 
                            VGroup( 
                                Item( name = 'gamma_N', label = 'safety factor normal force [-]:  gamma_N ' ),
                                Item( name = 'gamma_V', label = 'safety factor shear force [-]:  gamma_V ' ),
                                label = 'safety factors'
                                  ),
                            VGroup( 
                                Item( name = 'V_Rk', label = 'characteristic shear resistance [MPa]:  V_Rk ', format_str = "%.1f" ),
                                Item( name = 'N_Rk', label = 'characteristic axial resistance [MPa]:  N_Rk ', format_str = "%.1f" ),
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

LSLIST = [ SLS, ULS ]

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

    N_ip = Property( Array )
    def _get_N_ip( self ):
        return self.state_data['N_ip']

    V_ip = Property( Array )
    def _get_V_ip( self ):
        return self.state_data['V_ip']

    V_op = Property( Array )
    def _get_V_op( self ):
        return self.state_data['V_op']

    X_hf = Property( Array )
    def _get_X_hf( self ):
        return self.geo_data['X_hf']

    Y_hf = Property( Array )
    def _get_Y_hf( self ):
        return self.geo_data['Y_hf']

    Z_hf = Property( Array )
    def _get_Z_hf( self ):
        return self.geo_data['Z_hf']


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
    print(ls.ls_table)
