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
# Created on Dec 21, 2010 by: kelidas


from enthought.traits.api import HasTraits, Instance, on_trait_change, \
    Trait, Property, Event, Bool
from enthought.traits.ui.api import View, Item
from matplotlib.figure import Figure
from matplotlib.mlab import griddata
from numpy import ones_like, array, zeros_like, linspace, max, diff, min
from util.traits.editors.mpl_figure_editor import MPLFigureEditor
from ymb_data import \
    IYMBData, YMBSegmentData, YMBCutData, YMBSource, var_dict
import numpy.ma as ma
from matplotlib.colors import LinearSegmentedColormap

cdict = {
    'red':   ( ( 0.0, 1.0, 1. ), ( 0.5, 0.0, 0.0 ), ( 1., 0., 0. ) ),
    'green': ( ( 0.0, 0.0, 0.0 ), ( 0.5, 1.0, 1.0 ), ( 1, 0, 0 ) ),
    'blue':  ( ( 0.0, 0.0, 0.0 ), ( 0.5, 0.0, 0.0 ), ( 1., 1.0, 0.5 ) )
    }
my_cmap_lin = LinearSegmentedColormap( 'my_colormap_lin', cdict, 256 )


class YMBFieldVar( HasTraits ):
    data = Instance( IYMBData )

    n_cols = Property()
    def _get_n_cols( self ):
        return self.data.n_cuts

    var_enum = Trait( 'radius', var_dict, modified = True )
    scalar_arr = Property( depends_on = 'var_enum' )
    def _get_scalar_arr( self ):
        return getattr( self.data, self.var_enum_ )

    sorted_on = Bool( False, modified = True )

    scalar_arr_sorted = Property( depends_on = 'var_enum' )
    def _get_scalar_arr_sorted( self ):
        ''' Return scalar array sorted by the shortest distance from the edge 
        '''
        scalar_arr = zeros_like( getattr( self.data, self.var_enum_ ) )
        scalar_mask_arr = zeros_like( getattr( self.data, self.var_enum_ ) )
        distance_arr = self.data.edge_distance.filled()
        for i in range( 0, self.n_cols ):
            scalar_mask_arr[:, i] = zip( *sorted( zip( distance_arr[:, i], getattr( self.data, self.var_enum_ ).mask[:, i] ), reverse = True ) )[1]
            scalar_arr[:, i] = zip( *sorted( zip( distance_arr[:, i], getattr( self.data, self.var_enum_ ).filled()[:, i] ), reverse = True ) )[1]
        return ma.array( scalar_arr, mask = array( scalar_mask_arr, dtype = bool ) )

    figure = Instance( Figure, () )
    def _figure_default( self ):
        figure = Figure()
        figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        return figure

    data_changed = Event( True )
    @on_trait_change( '+modified, data' )
    def _redraw( self ):
        self.figure.clear()
        self.figure.add_axes( [0.1, 0.1, 0.8, 0.8] )
        figure = self.figure
        axes = figure.axes[0]
        axes.clear()

        if self.sorted_on == True:
            scalar_arr = self.scalar_arr_sorted
        else:
            scalar_arr = self.scalar_arr

        xi = linspace( min( self.data.cut_x ), max( self.data.cut_x ), 100 )

        x = ( ones_like( scalar_arr ) * self.data.cut_x ).flatten()
        ny_row = scalar_arr.shape[0]
        dy = max( diff( self.data.cut_x ) )
        yi = linspace( 0, ny_row * dy, ny_row )
        y = ( ones_like( scalar_arr ).T * linspace( 0, ny_row * dy, ny_row ) ).T.flatten()
        z = scalar_arr.flatten()
        zi = griddata( x, y, z, xi, yi, interp = 'nn' )

        # contour the gridded data, plotting dots at the nonuniform data points
        #axes.contour( xi, yi, zi, 20, linewidths = .5, colors = 'k' )
        # plotting filled contour
        axes.contourf( xi, yi, zi, 200, cmap = my_cmap_lin ) # my_cmap_lin
        scat = axes.scatter( x, y, marker = 'o', c = z, s = 20, linewidths = 0, cmap = my_cmap_lin )
        figure.colorbar( scat )

        self.data_changed = True

    view = View( 
                'var_enum',
                'sorted_on',
                Item( 'figure', style = 'custom',
                          editor = MPLFigureEditor(),
                          show_label = False ),
                id = 'yarn_structure_view',
                resizable = True,
                scrollable = True,
                dock = 'tab',
                width = 0.8,
                height = 0.4
                        )

if __name__ == '__main__':
    yarn_type = 'VET'
    my_model = YMBFieldVar( data = YMBCutData( source = YMBSource( yarn_type = yarn_type ) ) )
    my_model.configure_traits()
