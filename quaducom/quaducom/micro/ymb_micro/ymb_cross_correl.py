'''
Created on Nov 19, 2010

@author: kelidas
'''
from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, on_trait_change, List
from enthought.traits.ui.api import Group
from matplotlib.figure import Figure
from util.traits.editors.mpl_figure_editor import MPLFigureEditor
from numpy import arange
from ymb_data import IYMBData, YMBSegmentData, YMBSource, var_dict
import numpy.ma as ma
from ymb_auto_correl import MatSpearman

from enthought.traits.ui.api \
    import View, Item, TabularEditor

from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

class ArrayAdapter ( TabularAdapter ):

    columns = Property
    def _get_columns( self ):
        names = self.object.var_name_list
        cols = [ ( name, i ) for i, name in enumerate( names ) ]
#        n_columns = getattr( self.object, self.name ).shape[1]
#        cols = [ ( str( i ), i ) for i in range( n_columns ) ]
        return [ ( 'i', 'index' ) ] + cols

    font = 'Courier 15'
    alignment = 'right'
    format = '%6.2f'
    index_text = Property

    def _get_index_text ( self ):
        return self.object.var_name_list[ self.row ]

tabular_editor = TabularEditor( adapter = ArrayAdapter() )

class YMBCrossCorrel( HasTraits ):

    data = Instance( IYMBData )

    # convert the dictionary keys to an ordered list.
    var_name_list = Property( List )
    @cached_property
    def _get_var_name_list( self ):
        return sorted( var_dict.keys() )

    # list of data arrays in the order of the var_name_list
    var_arr_list = Property( List, depends_on = 'data.input_change' )
    @cached_property
    def _get_var_arr_list( self ):
        return [ getattr( self.data, var_dict[ var_name ] ).flatten()[:, None]
                 for var_name in self.var_name_list ]

    corr_arr = Property( Array, depends_on = 'data.input_change' )
    @cached_property
    def _get_corr_arr( self ):
        print 'redrawing cross correl'
        # get the list of names and sort them alphabetically
        corr_data = ma.hstack( self.var_arr_list )
        # @kelidas: return small differences between ma and numpy corrcoef
        #return ma.corrcoef( corr_data, rowvar = False, allow_masked = True )
        return MatSpearman( corr_data )

    figure = Instance( Figure )

    def _figure_default( self ):
        figure = Figure()
        figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        return figure

    data_changed = Event( True )
    @on_trait_change( 'data, data.input_change' )
    def _redraw( self ):
        figure = self.figure
        figure.clear()
        var_data = self.corr_arr

        figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        axes = figure.axes[0]
        axes.clear()
        x_coor = arange( var_data.shape[1] )
        axes.grid()
        for i in range( 0, var_data.shape[1] ):
            axes.plot( x_coor[i:] - x_coor[i] ,
                       var_data[i, ( i ):], '-x' )
        axes.set_xlabel( '$\mathrm{x}\, [\mu\mathrm{m}]$', fontsize = 16 )
        axes.set_ylabel( '$\mathrm{correlation}$', fontsize = 16 )
        axes.set_ylim( -1, 1 )

        self.data_changed = True

    traits_view_mpl = View( Group( 
#                       Group( Item( 'figure', style = 'custom',
#                              editor = MPLFigureEditor(),
#                              show_label = False )
#                              , id = 'figure.view' ),
                        Item( 'corr_arr',
                              show_label = False,
                              style = 'readonly',
                              editor = TabularEditor( adapter = ArrayAdapter() ) )
                       ),
                       resizable = True,
                        )

    traits_view = View( Item( 'corr_arr', editor = tabular_editor, show_label = False ),
                 resizable = True,
                 scrollable = True,
                 buttons = ['OK', 'Cancel' ],
                 width = 1.0, height = 0.5 )

if __name__ == '__main__':

    yarn_type = 'VET'
    data = YMBSegmentData( source = YMBSource( yarn_type = yarn_type ) )
    corr = YMBCrossCorrel( data = data )
    corr.configure_traits()
