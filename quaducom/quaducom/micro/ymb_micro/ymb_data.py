'''
Created on Nov 19, 2010

@author: kelidas
'''
from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Tuple, Interface, implements, Enum, Str
from enthought.traits.trait_types import DelegatesTo
from enthought.traits.ui.api import Item, View, HGroup, RangeEditor
from numpy import loadtxt, min, array, arange, ones_like, cumsum, vstack, \
    hstack, sum, zeros_like, zeros, ones, where, unique, pi, invert, \
    prod
from os.path import join
from promod.simdb import SimDB
from scipy.sparse import csr_matrix
import numpy.ma as ma
import os
import re

simdb = SimDB()
#data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'VET', 'raw_data' )
#data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'MAG', 'raw_data' )
#data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'TEST' )

var_dict = {'radius' : 'radius',
                       'cross sectional area':'cs_area',
                       'shortest distance from the edge' : 'edge_distance',
                       'contact fraction' : 'cf',
                       'slack' : 'slack',
                       'bond free length': 'bond_free_length',
                       'bond free length x': 'bond_free_length_x'}

yarn_list = ['MAG', 'VET', 'VET_110225', 'VET3c']

class YMBSource( HasTraits ):
    '''
        Read raw data files and convert length units px to micrometers
    '''
    yarn_type = Enum( yarn_list, changed_source = True )

    root_dir = Str( join( simdb.exdata_dir, 'trc', 'yarn_structure' ), changed_source = True )

    data_dir = Property( Directory, depends_on = '+changed_source' ) # Directory( data_dir, entries = 10, changed_source = True, auto_set = False, enter_set = True )
    @cached_property
    def _get_data_dir( self ):
        return join( self.root_dir, self.yarn_type, 'raw_data' )

    # constant -- convert pixels to milimeters
    px2mum = Constant( 1200, input = True )

    cut_raw_data = Property( Tuple, depends_on = 'data_dir' )
    @cached_property
    def _get_cut_raw_data( self ):
        '''
            Given a directory name, returns Tuple of NumPy arrays of variables 
            that can be obtained from raw data (y, z, r, d, contact fraction).
            Cut filenames must satisfy the following format '*Schnitt(num).txt'
            Return:
                y,z -- 1D arrays of coordinates
                r -- 1D array of diameters 
                d -- 1D array of distances from the edge
                cf -- 2D array of contact fraction components
                offset -- 1D array of number of filaments in the cut
        '''
        files = []
        num = []
        paths = os.listdir( self.data_dir )  # list of paths in that dir
        for fname in paths:
            match = re.search( r'\w+Schnitt(\d+).txt', fname ) # find files of given pattern
            if match:
                files.append( os.path.abspath( os.path.join( self.data_dir, fname ) ) )
                num.append( int( match.group( 1 ) ) )
        num, files = zip( *sorted( zip( num, files ) ) )
        y_arr = array( [] )
        z_arr = array( [] )
        r_arr = array( [] )
        d_arr = array( [] )
        cf_component_arr = array( [] ).reshape( 0, 36 )
        offset = array( [] )

        for n, file in zip( num, files ):
            print 'reading cut data_file %s -- slice %i' % ( file, n )
            data = loadtxt( file , delimiter = ' ',
                              skiprows = 1,
                              usecols = None )

            y = data[ :, 1] / self.px2mum
            z = data[ :, 2] / self.px2mum
            r = data[ :, 3] / self.px2mum
            d = data[ :, 4] / self.px2mum
            cf_component = data[ :, 5:41]

            y_arr = hstack( [y_arr, y] )
            z_arr = hstack( [ z_arr, z] )
            r_arr = hstack( [ r_arr, r] )
            d_arr = hstack( [ d_arr, d] )
            cf_component_arr = vstack( [cf_component_arr, cf_component] )
            offset = hstack( [offset, len( y )] )

        y_shift = min( y_arr )
        z_shift = min( z_arr )
        if y_shift < 0:
            print 'y-coordinites shifted by %f' % y_shift
            y_arr = y_arr - y_shift
        if z_shift < 0:
            print 'z-coordinites shifted by %f' % z_shift
            z_arr = z_arr - z_shift
        return y_arr, z_arr, r_arr, d_arr, cf_component_arr, cumsum( offset )

    # number of cuts in the raw data (numbering from zero)
    n_cuts = Property( Int, depends_on = '+changed_source' )
    @cached_property
    def _get_n_cuts( self ):
        n_cuts = self.connectivity.shape[1]
        return n_cuts

    # number of filaments in the raw data (numbering from zero)
    n_filaments = Property( Int, depends_on = '+changed_source' )
    @cached_property
    def _get_n_filaments( self ):
        #y_list = self.cut_raw_data[0]
        n_filaments = self.connectivity.shape[0] # max( [ cut.shape[0] for cut in y_list ] )
        return n_filaments

    # yarn length
    yarn_length = Property( Int, depends_on = '+changed_source' )
    @cached_property
    def _get_yarn_length( self ):
        length = self.cut_x[-1]
        return length

    # yarn cut coordinates
    # TODO: change with new data sets? probably yes
    cut_x = Property( Array, depends_on = '+changed_source' )
    @cached_property
    def _get_cut_x( self ):
        '''
            Read data with cut x-coordinates
            First cut has x-coordinate = 0
            Return NumPy array
        '''
        print 'reading data_file -- slice_height.txt'
        cut_x_file = join( self.data_dir, 'cut_x.txt' )
        cut_x = loadtxt( cut_x_file, dtype = 'float',
                                skiprows = 1, usecols = ( 1, 3 ) ) / self.px2mum
        cut_x = unique( cut_x.flatten() )
        return cut_x

    # filament connections between cuts
    connectivity = Property( Array, depends_on = '+changed_source' )
    @cached_property
    def _get_connectivity( self ):
        '''
            Read connectivity data.
            Return NumPy array
        '''
        print 'reading data_file -- connectivity.txt'
        connectivity_file = join( self.data_dir, 'connectivity.txt' )
        connect = loadtxt( connectivity_file, delimiter = ' ', dtype = 'float', skiprows = 1 )
        # replace 0 by -1 (-1 = not identified)
        # TODO: will be changed in raw data
        #connect[ equal( 0, connect )] = -1
        return connect.astype( int )

    traits_view = View( Item( 'yarn_type' ) )

class YMBData( HasTraits ):
    '''
        Preprocessing raw data
    '''

    source = Instance( YMBSource )

    cut_raw_data = Property( depends_on = 'source.+changed_source' )
    def _get_cut_raw_data( self ):
        return self.source.cut_raw_data

    n_cuts = Property( Int, depends_on = 'source.+changed_source' )
    def _get_n_cuts( self ):
        return self.source.n_cuts

    n_filaments = Property( Int, depends_on = 'source.+changed_source' )
    def _get_n_filaments( self ):
        return self.source.n_filaments

    yarn_length = DelegatesTo( 'source' )
    cut_x = DelegatesTo( 'source' )
    connectivity = DelegatesTo( 'source' )

    #------------------------------------------------------------------------
    # Change triggers
    #------------------------------------------------------------------------
    # the following events are introduced in order to react in the
    # clients of the YMBData to the changes in the input source
    # and in the configuration in a distinguished way.
    # 
    source_change = Event
    @on_trait_change( 'source.+changed_source' )
    def _set_source_change( self ):
        self.source_change = True

    input_change = Event
    @on_trait_change( 'source.+changed_source, +changed_config' )
    def _set_input_change( self ):
        self.input_change = True

    #-------------------------------------------------------------------------
    # Config option
    #-------------------------------------------------------------------------
    # limiting value for considering a filament as embedded
    #
    cf_limit = Range( value = .2, low = 0.0,
                      changed_config = True,
                      high = 1.0, enter_set = True, auto_set = False )

    # identify the filament points that have no  
    mask_arr = Property( Array, depends_on = 'source.+changed_source' )
    @cached_property
    def _get_mask_arr( self ):
        return self.connectivity < 0

    cut_data = Property( Tuple, depends_on = 'source.+changed_source' )
    @cached_property
    def _get_cut_data( self ):
        '''
            Prepare grid from raw data
            Output: 2D NumPy arrays (x, y, z, r, d, cf) 
            (cut (columns) x filament (rows))
        '''
        y_arr, z_arr, r_arr, d_arr, cf_arr, offset = self.cut_raw_data
        slice_height_arr = self.cut_x
        filam_connect_arr = self.connectivity
        mask_arr = self.mask_arr

        offset_arr = zeros_like( offset )
        offset_arr[1:] = offset[:-1]

        offset_arr = ones_like( filam_connect_arr ) * offset_arr
        map_arr = zeros_like( filam_connect_arr )
        map_arr[invert( mask_arr )] = filam_connect_arr[invert( mask_arr )] + offset_arr[invert( mask_arr )]
        map_arr[mask_arr] = -1
        map_arr = array( map_arr, dtype = 'int' )

        # stack 1d arrays in list into single array with -1 at the end
        y_arr = hstack( [y_arr, [-1]] )
        z_arr = hstack( [z_arr, [-1]] )
        r_arr = hstack( [r_arr, [-1]] )
        d_arr = hstack( [d_arr, [-1]] )

        # solve cf from contact fraction components (0-contact with matrix,
        # 1-contact with fiber, 2-no contact) -- 36 values around fiber diameter
        cf_arr = sum( cf_arr == 0, axis = 1 ) / 36.
        cf_arr = hstack( [ cf_arr , [-1]] )

        # map data arrays according to connection data
        x_arr = ones_like( map_arr ) * slice_height_arr

        x_arr = ma.array( x_arr, mask = mask_arr, fill_value = -1 )
        y_arr = ma.array( y_arr[map_arr], mask = mask_arr, fill_value = -1 )
        z_arr = ma.array( z_arr[map_arr], mask = mask_arr, fill_value = -1 )
        r_arr = ma.array( r_arr[map_arr], mask = mask_arr, fill_value = -1 )
        d_arr = ma.array( d_arr[map_arr], mask = mask_arr, fill_value = -1 )
        cf_arr = ma.array( cf_arr[map_arr], mask = mask_arr, fill_value = -1 )

        return x_arr, y_arr, z_arr, r_arr, d_arr, cf_arr

    #-----------------------------------------------------------------
    # Field variables
    #-----------------------------------------------------------------
    # naming:
    # fc - stands for var[ filament_idx, cut_idx ]
    # fs - stands for var[ filament_idx, segment_idx ]
    #  
    # cut-based fields come first
    #
    fc_radius = Property( Array )
    def _get_fc_radius( self ):
        return self.cut_data[3]

    fc_edge_distance = Property( Array )
    def _get_fc_edge_distance( self ):
        return self.cut_data[4]

    fc_cf = Property( Array )
    def _get_fc_cf( self ):
        return self.cut_data[5]

    # filament area
    fc_cs_area = Property( Array, depends_on = 'source.+changed_source' )
    @cached_property
    def _get_fc_cs_area( self ):
        area = pi * self.fc_radius ** 2
        return area

    # segment-based fields:
    # starting with fs_

    temp_filament = Array()
    temp_switch = Array()

    # real free filament length (sum of the filament length (linear interpolation) without matrix contact)
    fs_bond_free_length = Property( Array, depends_on = 'input_change' )
    @cached_property
    def _get_fs_bond_free_length( self ):

        x_data = hstack( [zeros( ( len( self.fs_length_between_cuts[:, 0] ), 1 ) ),
                          cumsum( self.fs_length_between_cuts, axis = 1 )] )
        cfl = self._free_length( x_data )
        return cfl

    fs_bond_free_length_x = Property( Array, depends_on = 'input_change' )
    @cached_property
    def _get_fs_bond_free_length_x( self ):

        x_data = ones( ( len( self.fs_length_between_cuts[:, 0] ), 1 ) ) * self.cut_x
        cfl = self._free_length( x_data )
        return cfl

    fs_slack = Property( Array, depends_on = 'input_change' )
    @cached_property
    def _get_fs_slack( self ):
        '''
            return slack in the cut of the free length part of the filament
        '''
        # slack from displacement that cause filament tension
        mask = self.fs_bond_free_length_x != 0
        slack = self._get_filament_dx().copy()
        slack[mask] /= self.fs_bond_free_length_x[mask]
        slack[self.fs_length_between_cuts.mask] = -1
        slack = ma.masked_array( slack, mask = self.fs_length_between_cuts.mask, fill_value = -1 )
        return slack

    # identify bond free parts of filaments 
    # and calculate distances between start and end point of these parts 
    def _free_length( self, data ):
        cf = self.fc_cf.copy()
        cf[cf <= self.cf_limit] = 0
        cf[cf != 0] = 1 # contact
        cf[cf == 0] = 0 #  no contact

        bound = cf[:, :-1] * cf[:, 1:]
        bound[bound > 0] = -1
        bound[bound == 0] = 1
        bound[self.fs_length_between_cuts.mask] = -1

        bound_data = hstack( [ -ones( ( bound.shape[ 0 ], 1 ) ), bound / abs( bound ),
                         - ones( ( bound.shape[ 0 ], 1 ) ) ] )

        filaments, switches = where( bound_data[:, 0:-1] * bound_data[:, 1:] < 0 )
        n_segments = len( filaments ) / 2.
        fil_idx_arr = filaments.reshape( ( n_segments, 2 ) )
        switch_idx_arr = switches.reshape( ( n_segments, 2 ) )

        segment_lengths = data[ fil_idx_arr[:, 1], switch_idx_arr[:, 1] ] - data[ fil_idx_arr[:, 0], switch_idx_arr[:, 0] ]

        # data for sparse matrix CSR
        data_row = segment_lengths.repeat( switch_idx_arr[:, 1] - switch_idx_arr[:, 0] )
        # row position of data_row
        row = fil_idx_arr[:, 0].repeat( switch_idx_arr[:, 1] - switch_idx_arr[:, 0] )

        # position assuming that bound_data.flatten()
        switch_idx_arr = ( switch_idx_arr + fil_idx_arr * bound_data.shape[1] ).flatten()
        switch_idx_arr2 = switch_idx_arr[1:] - switch_idx_arr[:-1]

        # create array  [1,0,1,0,......]
        aran = arange( 0, len( switch_idx_arr2 ) )
        mask = aran % 2 == 0
        aran[mask] = True
        aran[invert( mask )] = False
        # repeat values
        aran = ( aran.repeat( switch_idx_arr2 ) ).astype( bool )

        a = arange( min( switch_idx_arr ), max( switch_idx_arr ) )
        # column position of data_row
        col = a[aran] - row * bound_data.shape[1]

        cfl_real = array( csr_matrix( ( data_row, ( row, col ) ), shape = bound.shape ).todense() )
        cfl_real[self.fs_length_between_cuts.mask] = -1
        cfl_real = ma.array( cfl_real, mask = self.fs_length_between_cuts.mask )
        return cfl_real

    def _get_filament_dx( self ):
        ''' Displacement of filament end-point to become stretched. '''
        y_arr, z_arr = self.cut_data[1:3]
        dy_free = self._free_length( y_arr )
        dz_free = self._free_length( z_arr )
        dx = ma.sqrt( self.fs_bond_free_length ** 2
                - dy_free ** 2
                - dz_free ** 2 ) - self.fs_bond_free_length_x
        return dx

    # length obtained by linear interpolation between cuts
    fs_length_between_cuts = Property( Array, depends_on = 'source.+changed_source' )
    @cached_property
    def _get_fs_length_between_cuts( self ):
        '''
            Return linear filament length between cuts.
            Matrix contain mask information.
        '''
        x_arr = self.cut_data[0]
        y_arr = self.cut_data[1]
        z_arr = self.cut_data[2]
        length = ma.sqrt( ( x_arr[:, 1:] - x_arr[:, :-1] ) ** 2
                + ( y_arr[:, 1:] - y_arr[:, :-1] ) ** 2
                + ( z_arr[:, 1:] - z_arr[:, :-1] ) ** 2 )
        return length

    # yarn area in the cut
    cut_area = Property( Array, depends_on = 'source.+changed_source' )
    @cached_property
    def _get_cut_area( self ):
        area = sum( self.cs_area, axis = 0 )
        return area.reshape( 1, self.n_cuts )

    # number of filaments in the raw data (numbering from zero)
    # theoretical number of continuous filaments (approx)
    n_continuous_filaments_aver = Property( Int, depends_on = 'source.+changed_source' )
    @cached_property
    def _get_n_continuous_filaments_aver( self ):
        n_fil = sum( self.fs_length_between_cuts ) / self.yarn_length
        return int( n_fil )

    n_continuous_filaments = Property( Int, depends_on = 'source.+changed_source' )
    @cached_property
    def _get_n_continuous_filaments( self ):
        n_fil = sum( prod( invert( self.mask_arr ), axis = 1 ) )
        return int( n_fil )

    traits_view = View( HGroup( Item( 'source@', show_label = False, springy = True ),
                                Item( 'cf_limit', springy = True ),
                                Item( 'n_continuous_filaments', style = 'readonly', label = 'continuous', tooltip = 'Number of continuous filaments' ),
                                Item( 'n_filaments', style = 'readonly', label = 'out of', tooltip = 'Total number of filaments' ),
                                Item( 'n_continuous_filaments_aver', style = 'readonly', label = 'aver', tooltip = 'Average number of continuous filaments' ),
                                 ),
                                 )

class IYMBData( Interface ):
    '''Interface required from the data used by YMBSlider
    '''
    n_cols = Property

    x_coord = Property

    radius = Property()

    edge_distance = Property()

    cf = Property()

    cs_area = Property()

    bond_free_length = Property()

    bond_free_length_x = Property()

    slack = Property()

    length_between_cuts = Property()


class YMBCutData( YMBData ):
    '''
    '''
    implements( IYMBData )

    n_cols = Property
    def _get_n_cols( self ):
        return self.n_cuts

    x_coord = Property
    def _get_x_coord( self ):
        return self.source.cut_x

    radius = Property()
    def _get_radius( self ):
        return self.fc_radius

    edge_distance = Property()
    def _get_edge_distance( self ):
        return self.fc_edge_distance

    cf = Property()
    def _get_cf( self ):
        return self.fc_cf

    cs_area = Property()
    def _get_cs_area( self ):
        return self.fc_cs_area

    def extrapolate_segments( self, arr ):
        s = arr.shape
        fs_arr = zeros( ( s[0], s[1] + 1 ), dtype = float )
        fs_arr[:, :-1] += arr.filled()
        fs_arr[:, 1:] += arr.filled()
        fs_arr[:, 1:-1] /= 2
        mask = ones( ( s[0], s[1] + 1 ), dtype = bool )
        a = invert( arr.mask )
        mask[:, :-1] = a#[:, :-1]
        mask[:, 1:] *= a#[:, 1:]
        mask = invert( mask )
        return ma.array( fs_arr, mask = mask, fill_value = -1 )

    bond_free_length = Property()
    def _get_bond_free_length( self ):
        return self.extrapolate_segments( self.fs_bond_free_length )

    bond_free_length_x = Property()
    def _get_bond_free_length_x( self ):
        return self.extrapolate_segments( self.fs_bond_free_length_x )

    slack = Property()
    def _get_slack( self ):
        return self.extrapolate_segments( self.fs_slack )

    length_between_cuts = Property()
    def _get_length_between_cuts( self ):
        return self.extrapolate_segments( self.fs_length_between_cuts )

class YMBSegmentData( YMBData ):
    '''
    '''
    n_cols = Property
    def _get_n_cols( self ):
        return self.n_cuts - 1

    x_coord = Property
    def _get_x_coord( self ):
        cut_x = self.source.cut_x
        return ( cut_x[1:] + cut_x[:-1] ) / 2.0

    radius = Property()
    def _get_radius( self ):
        return ( self.fc_radius[:, 1:] + self.fc_radius[:, :-1] ) / 2.0

    edge_distance = Property()
    def _get_edge_distance( self ):
        return ( self.fc_edge_distance[:, 1:] + self.fc_edge_distance[:, :-1] ) / 2.0

    cf = Property()
    def _get_cf( self ):
        return ( self.fc_cf[:, 1:] + self.fc_cf[:, :-1] ) / 2.0

    cs_area = Property()
    def _get_cs_area( self ):
        return ( self.fc_cs_area[:, 1:] + self.fc_cs_area[:, :-1] ) / 2.0

    bond_free_length = Property()
    def _get_bond_free_length( self ):
        return self.fs_bond_free_length

    bond_free_length_x = Property()
    def _get_bond_free_length_x( self ):
        return self.fs_bond_free_length_x

    slack = Property()
    def _get_slack( self ):
        return self.fs_slack

    length_between_cuts = Property()
    def _get_length_between_cuts( self ):
        return self.fs_length_between_cuts

class YMBSlider( HasTraits ):
    '''
        Slicing data arrays (2d NumPy array)
        Return data for statistics (1d NumPy array)
    '''
    data = Instance( IYMBData )

    n_cols = Property
    def _get_n_cols( self ):
        return self.data.n_cols - 1

    n_rows = Property
    def _get_n_rows( self ):
        return self.data.n_filaments - 1

    var_enum = Trait( 'radius',
                      {'radius' : 'radius',
                       'cross sectional area':'cs_area',
                       'yarn area in the cut':'cut_area',
                       'shortest distance from the edge' : 'edge_distance',
                       'contact fraction' : 'cf',
                       'slack' : 'slack',
                       'bond free length': 'bond_free_length',
                       'bond free length x': 'bond_free_length_x'} )

    zero = Constant( 0 )
    cut_slider = Int( changed_range = True )
    cut_slider_on = Bool( True, changed_range = True )

    filament_slider = Int( changed_range = True )
    filament_slider_on = Bool( False, changed_range = True )

    # bundled change event to simplify the outside dependencies
    input_change = Event
    @on_trait_change( 'var_enum, +changed_range, data.+changed_source, data.+changed_config' )
    def _set_input_change( self ):
        self.input_change = True

    stat_data = Property( Array,
                          depends_on = 'input_change' )
    @cached_property
    def _get_stat_data( self ):
        if self.cut_slider_on == True:
            return getattr( self.data, self.var_enum_ )[:, self.cut_slider].compressed()
        if self.filament_slider_on == True:
            return  getattr( self.data, self.var_enum_ )[self.filament_slider, :].compressed()
        if self.cut_slider_on == False and self.filament_slider_on == False:
            return  getattr( self.data, self.var_enum_ )[ : , : ].compressed()

    traits_view = View( #VGroup( Item( 'data@', show_label=False ) ),
                        Item( 'var_enum', label = 'Variable' ),
                        HGroup( 
                        Item( 'cut_slider_on', label = 'Cut slider' ),
                        Item( 'cut_slider', show_label = False, springy = True,
                              enabled_when = 'cut_slider_on == True',
                              editor = RangeEditor( low_name = 'zero',
                                                    high_name = 'n_cols',
                                                    mode = 'slider',
                                                    auto_set = False,
                                                    enter_set = False,
                                                    ),
                               ),
                        ),
                        HGroup( 
                        Item( 'filament_slider_on', label = 'Filament slider' ),
                        Item( 'filament_slider', show_label = False, springy = True,
                              enabled_when = 'filament_slider_on == True',
                              editor = RangeEditor( low_name = 'zero',
                                                    high_name = 'n_rows',
                                                    mode = 'slider',
                                                    auto_set = False,
                                                    enter_set = False,
                                                    ),
                                )
                        ) )



if __name__ == '__main__':

    yarn_type = 'VET'
    data = YMBCutData( source = YMBSource( yarn_type = yarn_type ), cf_limit = .5 )
    yarn = YMBSlider( data = data )
    yarn.configure_traits()





