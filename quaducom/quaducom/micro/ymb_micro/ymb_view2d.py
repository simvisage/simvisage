'''
Created on Nov 20, 2010

@author: kelidas
'''

from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Event, Array, Instance, Range, on_trait_change, Bool, Trait, DelegatesTo, \
    Constant
from enthought.traits.ui.api import View, Item, Group
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from numpy import hstack
from util.traits.editors.mpl_figure_editor import MPLFigureEditor
from numpy import ceil, max
from ymb_data import YMBData, YMBSource, YMBSegmentData, var_dict


cdict = {
    'red':   ( ( 0.0, 1.0, 1.0 ), ( 0.5, 0.0, 0.0 ), ( 1.0, 0.0, 0.0 ) ),
    'green': ( ( 0.0, 0.0, 0.0 ), ( 0.5, 1.0, 1.0 ), ( 1.0, 0.0, 0.0 ) ),
    'blue':  ( ( 0.0, 0.0, 0.0 ), ( 0.5, 0.0, 0.0 ), ( 1.0, 1.0, 0.5 ) )
    }
my_cmap_lin = LinearSegmentedColormap( 'my_colormap_lin', cdict, 256 )


class YMBView2D( HasTraits ):

    data = Instance( YMBData )

    zero = Constant( 0 )
    slider_max = Property()
    def _get_slider_max( self ):
        return self.data.n_cuts - 1

    var_enum = Trait( 'radius',
                      var_dict, modified = True )

    cut_slider = Range( 'zero', 'slider_max', mode = 'slider', auto_set = False,
                         enter_set = True, modified = True )

    circle_diameter = Float( 20, enter_set = True, auto_set = False, modified = True )

    underlay = Bool( False, modified = True )

    variable = Property( Array, depends_on = 'var_enum' )
    @cached_property
    def _get_variable( self ):
        return getattr( self.data, self.var_enum_ )

    figure = Instance( Figure )

    def _figure_default( self ):
        figure = Figure()
        figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        return figure

    data_changed = Event( True )
    @on_trait_change( '+modified, data.input_changed' )
    def _redraw( self ):
        # TODO: set correct ranges, fix axis range (axes.xlim)
        self.figure.clear()
        self.figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        figure = self.figure
        axes = figure.axes[0]
        axes.clear()
        y_arr, z_arr = self.data.cut_data[1:3]
        y_raw_arr, z_raw_arr = self.data.cut_raw_data[0:2]
        offset = hstack( [0, self.data.cut_raw_data[5]] )
        scalar_arr = self.variable
        mask = y_arr[:, self.cut_slider] > -1

        axes.scatter( y_raw_arr[offset[self.cut_slider]:offset[self.cut_slider + 1]],
                      z_raw_arr[offset[self.cut_slider]:offset[self.cut_slider + 1]],
                      s = self.circle_diameter, color = 'k', marker = 'x', label = 'identified filament in cut' )
        scat = axes.scatter( y_arr[:, self.cut_slider][mask], z_arr[:, self.cut_slider][mask],
                      s = self.circle_diameter, c = scalar_arr[:, self.cut_slider][mask], cmap = my_cmap_lin, label = 'connected filaments' )
        axes.set_xlabel( '$y\, [\mathrm{mm}]$', fontsize = 16 )
        axes.set_ylabel( '$z\, [\mathrm{mm}]$', fontsize = 16 )

        axes.set_xlim( [0, ceil( max( y_arr ) )] )
        axes.set_ylim( [0, ceil( max( z_arr ) )] )
        axes.legend()
        figure.colorbar( scat )

        if self.underlay == True:
            axes.text( axes.get_xlim()[0], axes.get_ylim()[0],
                        'That\'s all at this moment :-)', color = 'red', fontsize = 20 )
        self.data_changed = True


    traits_view = View( Group( Item( 'var_enum' ),
                       Item( 'cut_slider', springy = True ),
                       Item( 'circle_diameter', springy = True ),
                       Item( 'underlay', springy = True ),
                        ),
                        Item( 'figure', style = 'custom',
                              editor = MPLFigureEditor(),
                              show_label = False ),
                       resizable = True,
                        )




if __name__ == '__main__':
    yarn_type = 'VET'
    cut = YMBView2D( data = YMBSegmentData( source = YMBSource( yarn_type = yarn_type ) ) )
    cut.configure_traits()
