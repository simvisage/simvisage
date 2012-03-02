'''
Created on Dec 13, 2010

@author: kelidas
'''
from enthought.traits.api import Float, Property, cached_property, Int, \
    on_trait_change, Interface, implements
from stats.pdistrib.pdistrib import IPDistrib
from numpy import mean, hstack
from ymb_hist import YMBHist
from ymb_data import YMBCutData


class YMBDistrib( YMBHist ):

    implements( IPDistrib )

    n_segments = Property
    def _get_n_segments( self ):
        return self.bins

    def _set_n_segments( self, value ):
        self.bins = value

    x_array = Property( depends_on = '+modified, slider.input_change' )
    @cached_property
    def _get_x_array( self ):
        h, b = self.normed_hist
        return b

    pdf_array = Property( depends_on = '+modified, slider.input_change' )
    @cached_property
    def _get_pdf_array( self ):
        # get the normed histogram implemented in the base class YMBHist
        h, b = self.normed_hist
        return h

    dx = Property()
    def _get_dx( self ):
        return self.bin_width

    def get_pdf_array( self, x_array ):
        # @todo - the supplied x_array must considered.
        #
        return self.pdf_array

    data_mean = Property( Float, depends_on = '+modified, slider.input_change' )
    @cached_property
    def _get_data_mean( self ):
        data = self.slider.stat_data
        return mean( data[data >= 0] )

    @on_trait_change( '+modified, slider.input_change' )
    def _redraw( self ):
        figure = self.figure
        axes = figure.axes[0]
        axes.clear()

        var_data = self.slider.stat_data

        axes.hist( var_data, bins = self.bins, normed = True,
                   edgecolor = self.edge_color, facecolor = self.face_color )

        x_arr = self.x_array[:-1] + self.dx / 2.0

        axes.plot( x_arr, self.pdf_array, 'ro' )
        axes.set_xlabel( self.slider.var_enum )

        if self.ylimit_on == True:
            axes.set_ylim( 0, self.ylimit )
        self.data_changed = True

if __name__ == '__main__':

    from ymb_data import YMBData, YMBSlider, YMBSource

    data = YMBCutData( source = YMBSource() )

    slider = YMBSlider( data = data, var_enum = 'contact fraction' )

    ymb_pd = YMBDistrib( slider = slider, n_int = 20 )

    print ymb_pd.x_array.shape
    print ymb_pd.pdf_array.shape

    ymb_pd.configure_traits()


