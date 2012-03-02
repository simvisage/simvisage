'''
Created on Nov 19, 2010

@author: kelidas
'''
from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Range, on_trait_change, Bool, Trait, Constant, Int, Str, List
from enthought.traits.ui.api import Group, HGroup, View, Item
from matplotlib.figure import Figure
from numpy import min, array, histogram2d, vstack, hstack, corrcoef, prod, ones, invert
from scipy.optimize import leastsq
from scipy.stats import spearmanr
from util.traits.editors.mpl_figure_editor import MPLFigureEditor
from ymb_data import IYMBData, YMBSegmentData, YMBSource, var_dict, YMBCutData
from matplotlib.pyplot  import setp
import numpy.ma as ma

def MatSpearman( mat ):
    ''' Return correlation matrix
    '''
    n = mat.shape[1]
    corr = ones( ( n, n ) )
    for i in range( 0, n ):
        for j in range( i + 1, n ):
            x = mat[:, i]
            y = mat[:, j]
            maska = x.mask * y.mask
            x = ma.array( x.filled(), mask = maska )
            y = ma.array( y.filled(), mask = maska )
            corr[i, j] = spearmanr( x.compressed(), y.compressed() )[0]
            corr[j, i] = spearmanr( x.compressed(), y.compressed() )[0]
    return corr

class YMBAutoCorrel( HasTraits ):

    data = Instance( IYMBData )

    var_enum = Trait( 'radius',
                      var_dict )

    input_change = Event
    @on_trait_change( 'var_enum, data.input_change' )
    def _set_input_change( self ):
        print 'YMBAutoCorrel input change'
        self.input_change = True

    corr_arr = Property( Array, depends_on = 'var_enum, data.input_change' )
    @cached_property
    def _get_corr_arr( self ):
        corr_data = getattr( self.data, self.var_enum_ )
        # @kelidas: return small differences between ma and numpy corrcoef
        #print MatSpearman( corr_data )
        #return ma.corrcoef( corr_data, rowvar = False, allow_masked = True )
        return MatSpearman( corr_data )


    fit_correl = Property()
    def _get_fit_correl( self ):
        x_coor = self.data.x_coord
        var_data = self.corr_arr
        x = []
        y = []
        for i in range( 0, var_data.shape[1] ):
            x.append( x_coor[i:] - x_coor[i] )
            y.append( var_data[i, ( i ):] )
        x = hstack( x )
        y = hstack( y )
        p0 = [1., 1., 1., 1.]
        plsq = leastsq( self.residual_ls, p0, args = ( y, x ) )
        return plsq[0]

    def residual_ls( self, p, y, x ):
        err = y - self.peval( x, p )
        return err

    def peval( self, x, p ):
        return  p[0] * x ** 3 + p[1] * x ** 2 + p[2] * x + p[3]

    traits_view = View( Item( 'var_enum', label = 'Variable' ) )


class YMBAutoCorrelView( HasTraits ):

    correl_data = Instance( YMBAutoCorrel )

    axes_adjust = List( [0.1, 0.1, 0.8, 0.8] )

    data = Property
    def _get_data( self ):
        return self.correl_data.data

    zero = Constant( 0 )
    slider_max = Property()
    def _get_slider_max( self ):
        return self.data.n_cuts - 1

    cut_slider = Range( 'zero', 'slider_max', mode = 'slider', auto_set = False, enter_set = True, modified = True )
    vcut_slider = Range( 'zero', 'slider_max', mode = 'slider', auto_set = False, enter_set = True, modified = True )

    cut_slider_on = Bool( False, modified = True )

    color = Str( 'blue' )

    figure = Instance( Figure )

    def _figure_default( self ):
        figure = Figure()
        figure.add_axes( self.axes_adjust )
        return figure

    data_changed = Event( True )
    @on_trait_change( 'correl_data.input_change, +modified' )
    def _redraw( self ):
        # TODO: set correct ranges, fix axis range (axes.xlim)
        print 'redrawing xxxx'
        figure = self.figure
        figure.clear()
        var_data = self.correl_data.corr_arr
        id = self.cut_slider
        if self.cut_slider_on == True:
            i = self.cut_slider
            j = self.vcut_slider
            plot_data = getattr( self.data, self.correl_data.var_enum_ )
            #plot_data = vstack( [plot_data[:, i], plot_data[:, j]] ).T
            # plot only values > -1
            #plot_data = plot_data[prod( plot_data >= 0, axis = 1, dtype = bool )]
            plot_data_x = plot_data[:, i]
            plot_data_y = plot_data[:, j]
            plot_data_corr = min( corrcoef( plot_data_x, plot_data_y ) )
            plot_data_corr_spear = spearmanr( plot_data_x, plot_data_y )[0]

            left, width = 0.1, 0.65
            bottom, height = 0.1, 0.65
            bottom_h = left_h = left + width + 0.02

            rect_scatter = [left, bottom, width, height]
            rect_histx = [left, bottom_h, width, 0.2]
            rect_histy = [left_h, bottom, 0.2, height]

            axScatter = figure.add_axes( rect_scatter )
            axHistx = figure.add_axes( rect_histx )
            axHisty = figure.add_axes( rect_histy )
            axScatter.clear()
            axHistx.clear()
            axHisty.clear()

            from matplotlib.ticker import NullFormatter
            axHistx.xaxis.set_major_formatter( NullFormatter() )
            axHisty.yaxis.set_major_formatter( NullFormatter() )

            axScatter.scatter( plot_data_x,
                               plot_data_y )

            #binwidth = 0.25
            #xymax = max( [max( abs( self.data.cf[:, j] ) ), max( abs( self.data.cf[:, i] ) )] )
            #lim = ( int( xymax / binwidth ) + 1 ) * binwidth

            #axScatter.set_xlim( ( -lim, lim ) )
            #axScatter.set_ylim( ( -lim, lim ) )

            #bins = arange( -lim, lim + binwidth, binwidth )
            axHistx.hist( plot_data_x.compressed(), bins = 40 )
            axHisty.hist( plot_data_y.compressed(), bins = 40, orientation = 'horizontal' )
            axHistx.set_xlim( axScatter.get_xlim() )
            axHisty.set_ylim( axScatter.get_ylim() )

            axScatter.set_xlabel( '$\mathrm{cut\, %i}$' % self.cut_slider, fontsize = 16 )
            axScatter.set_ylabel( '$\mathrm{cut\, %i}$' % self.vcut_slider, fontsize = 16 )
            axScatter.text( axScatter.get_xlim()[0], axScatter.get_ylim()[0],
                             'actual set correlation %.3f (Pearson), %.3f (Spearman)' % ( plot_data_corr, plot_data_corr_spear ), color = 'r' )

        if self.cut_slider_on == False:
            figure.add_axes( self.axes_adjust )
            axes = figure.axes[0]
            axes.clear()
            x_coor = self.data.x_coord
            axes.grid()
            for i in range( 0, var_data.shape[1] ):
                axes.plot( x_coor[i:] - x_coor[i],
                           var_data[i, ( i ):], '-x', color = self.color )
            # approximate by the polynomial (of the i-th order)
            #axes.plot( x_coor, self.correl_data.peval( x_coor, self.correl_data.fit_correl ), 'b', linewidth = 3 )
            setp( axes.get_xticklabels(), position = ( 0, -.025 ) )
            axes.set_xlabel( '$x \, [\mathrm{mm}]$', fontsize = 15 )
            axes.set_ylabel( '$\mathrm{correlation}$', fontsize = 15 )
            axes.set_ylim( -1, 1 )

        self.data_changed = True


    traits_view = View( Group( Item( 'correl_data', show_label = False, style = 'custom' ),
                              HGroup( 
                       Item( 'cut_slider_on', label = 'Scatter' ),
                       Item( 'cut_slider', show_label = False, springy = True, enabled_when = 'cut_slider_on == True' ),
                       Item( 'vcut_slider', show_label = False, springy = True, enabled_when = 'cut_slider_on == True' ),
                       ),
                       Group( Item( 'figure', style = 'custom',
                              editor = MPLFigureEditor(),
                              show_label = False )
                              , id = 'figure.view' ),
                       ),
                       resizable = True,
                        )



if __name__ == '__main__':

    yarn_type = 'VET'
    source = YMBSource( yarn_type = yarn_type )

    data = YMBCutData( source = source, cf_limit = .2 )
    corr = YMBAutoCorrelView( correl_data = YMBAutoCorrel( data = data ) )
    corr.configure_traits()

#    from matplotlib.colors import LinearSegmentedColormap
#
#    from numpy import mgrid, max
#    mg = mgrid[0:100:36j, 0:100:36j]
#    d = corr.data.cf
#    d2 = d[:, 0:2][prod( d[:, 0:2] >= 0, axis = 1, dtype = bool )]
#    a, b, c = histogram2d( d2[:, 0], d2[:, 1], bins = 36 )
#
#    # my own colormap
#    cdict = {
#    'red':   ( ( 0.0, 1.0, 1.0 ),
#                ( 1 / max( a ) , 1.0, 1.0 ),
#                 ( 1 / max( a ), 0.0, 0.0 ),
#                  ( 1.0, 0.0, 0.0 ) ),
#    'green': ( ( 0.0, 0.0, 0.0 ),
#                ( 1 / max( a ) , 0.0, 0.0 ),
#                 ( 1 / max( a ), 1.0, 1.0 ),
#                  ( 1.0, 0.0, 0.0 ) ),
#    'blue':  ( ( 0.0, 0.0, 0.0 ),
#                ( 1 / max( a ) , 0.0, 0.0 ),
#                 ( 1 / max( a ), 0.0, 0.0 ),
#                  ( 1.0, 1.0, 0.5 ) ) }
#
#    my_cmap_lin = LinearSegmentedColormap( 'my_colormap_lin', cdict, 256 )
#    plt.title( 'Contact fraction scatter -- cut0 x cut1' )
#    scat = plt.scatter( mg[0], mg[1], c = a, s = 40, cmap = my_cmap_lin )
#    plt.colorbar( scat )
#    plt.show()
#

