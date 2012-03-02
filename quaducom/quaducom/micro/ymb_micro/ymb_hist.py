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
# Created on Dec 13, 2010 by: rch

from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Instance, Int, on_trait_change, Bool, Str, Tuple, List, Float
from enthought.traits.ui.api import Item, View, Group, HGroup, HSplit
from util.traits.editors.mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
from pylab import setp
from numpy import histogram, prod, invert, linspace, sqrt, mean, var
from ymb_data import YMBCutData, YMBSegmentData, YMBSlider, YMBSource


class YMBHist( HasTraits ):

    slider = Instance( YMBSlider )

    figure = Instance( Figure )

    bins = Int( 20, auto_set = False, enter_set = True, modified = True )
    xlimit_on = Bool( False, modified = True )
    ylimit_on = Bool( False, modified = True )
    xlimit = Float( 100, auto_set = False, enter_set = True, modified = True )
    ylimit = Float( 100, auto_set = False, enter_set = True, modified = True )
    multi_hist_on = Bool( False, modified = True )
    stats_on = Bool( False, modified = True )
    normed_on = Bool( False, modified = True )

    normed_hist = Property( depends_on = '+modified, slider.input_change' )
    @cached_property
    def _get_normed_hist( self ):
        data = self.slider.stat_data
        h, b = histogram( data, bins = self.bins, normed = True )
        return h, b

    range = Property
    def _get_range( self ):
        h, b = self.normed_hist
        return ( min( b ), max( b ) )

    bin_width = Property
    def _get_bin_width( self ):
        return ( self.range[1] - self.range[0] ) / self.bins

    def _figure_default( self ):
        figure = Figure()
        figure.add_axes( self.axes_adjust )
        return figure

    edge_color = Str( None )
    face_color = Str( None )
    axes_adjust = List( [ 0.1, 0.1, 0.8, 0.8 ] )

    data_changed = Event( True )
    @on_trait_change( '+modified, slider.input_change' )
    def _redraw( self ):
        figure = self.figure
        axes = figure.axes[0]
        axes.clear()
        if self.multi_hist_on == True:
            histtype = 'step'
            lw = 3
            plot_data = getattr( self.slider.data, self.slider.var_enum_ )
            for i in range( 0, plot_data.shape[1] ):
                axes.hist( plot_data[:, i].compressed(), bins = self.bins, histtype = histtype, color = 'gray' )
        if self.multi_hist_on == False:
            histtype = 'bar'
            lw = 1
        var_data = self.slider.stat_data

        axes.hist( var_data, bins = self.bins, histtype = histtype, linewidth = lw, \
                   normed = self.normed_on, edgecolor = self.edge_color, facecolor = self.face_color )
        if self.stats_on == True:
            xint = axes.xaxis.get_view_interval()
            yint = axes.yaxis.get_view_interval()
            axes.text( xint[0], yint[0], 'mean = %e, std = %e' % ( mean( var_data ), sqrt( var( var_data ) ) ) )

        # redefine xticks labels
        #inter = axes.xaxis.get_view_interval()
        #axes.set_xticks( linspace( inter[0], inter[1], 5 ) )
        axes.set_xlabel( self.slider.var_enum )
        axes.set_ylabel( 'frequency' )#, fontsize = 16
        setp( axes.get_xticklabels(), position = ( 0, -.025 ) )
        if self.xlimit_on == True:
            axes.set_xlim( 0, self.xlimit )
        if self.ylimit_on == True:
            axes.set_ylim( 0, self.ylimit )
        self.data_changed = True

    view = View( 
                   Group( Item( 'figure', style = 'custom',
                                  editor = MPLFigureEditor(),
                                  show_label = False, id = 'figure.view' ),
                    HGroup( 
                             Item( 'bins' ),
                            Item( 'ylimit_on', label = 'Y limit' ),
                            Item( 'ylimit', enabled_when = 'ylimit_on == True',
                                          show_label = False ),
                            Item( 'stats_on', label = 'stats' ),
                            Item( 'normed_on', label = 'norm' ),
                            Item( 'multi_hist_on', label = 'multi' ) ),
                    label = 'histogram',
                    dock = 'horizontal',
                    id = 'yarn_hist.figure',
                ),
                Group( Item( 'slider', style = 'custom', show_label = False ),
                    label = 'yarn data',
                    dock = 'horizontal',
                    id = 'yarn_hist.config',
                ),
                id = 'yarn_structure_view',
                resizable = True,
                scrollable = True,
                #width = 0.8,
                #height = 0.4
                        )


if __name__ == '__main__':

    yarn_type = 'VET'
    source = YMBSource( yarn_type = yarn_type )

    data = YMBCutData( source = source, cf_limit = 0.5 )

    slider = YMBSlider( var_enum = 'slack', data = data )
    yarn = YMBHist( slider = slider )

    yarn.configure_traits()

