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
# Created on Nov 29, 2010 by: rch

from numpy import array, where, nonzero, array, hstack, ones, linspace, zeros_like, \
                arange, invert
from scipy.sparse import csr_matrix

def i_ident():
    data = array( [ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1 ] )

    switches = where( data[0:-1] * data[1:] < 0 )[0]

    print data
    print 'switches'
    switches = hstack( [array( [-1] ), switches, [len( data ) - 1]] )
    print switches

    print 'current cut position'
    cut_pos = 11

    print 'find the indices of the interval at 5-th cut'
    idx = where( switches >= cut_pos )[0][0]

    right_limit = switches[ idx ] + 1
    left_limit = switches[ idx - 1 ] + 1

    print 'right_limit', right_limit
    print 'left_limit', left_limit

    print data[ left_limit: right_limit ]

def i_ident_2d():
    from numpy import searchsorted, split, dtype
    from numpy.random import random
    data_bas = array( [[ 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1 ],
                   [ -1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1 ],
                   [ -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]] )
    #data_bas = random( ( 2000, 15 ) ) - .5
    x_data = linspace( 0, 10, data_bas.shape[1] )

    data = hstack( [ -ones( ( data_bas.shape[ 0 ], 1 ) ), data_bas / abs( data_bas ),
                     - ones( ( data_bas.shape[ 0 ], 1 ) ) ] )
    print 'data', data

    x_data = hstack( [ x_data[0], x_data, x_data[-1] ] )
    print 'x_data', x_data

    r = range( 0, data.shape[0] )
    filaments, switches = where( data[:, 0:-1] * data[:, 1:] < 0 )
    print 'filaments', filaments
    print 'switches', switches
    n_segments = len( filaments ) / 2
    fil_idx_arr = filaments.reshape( ( n_segments, 2 ) )
    switch_idx_arr = switches.reshape( ( n_segments, 2 ) )

    print fil_idx_arr
    print switch_idx_arr

    segment_lengths = x_data[ switch_idx_arr[:, 1] + 1 ] - x_data[ switch_idx_arr[:, 0] ]
    print 'segment lengths', segment_lengths

    # data for sparse matrix CSR
    data_row = segment_lengths.repeat( switch_idx_arr[:, 1] - switch_idx_arr[:, 0] )
    # row info for data
    row = fil_idx_arr[:, 0].repeat( switch_idx_arr[:, 1] - switch_idx_arr[:, 0] )

    print 'data_row', data_row
    print 'row', row

    # position assuming that data.flatten()
    switch_idx_arr = ( switch_idx_arr + fil_idx_arr * data.shape[1] ).flatten()
    print 'switch idx array', switch_idx_arr
    switch_idx_arr2 = switch_idx_arr[1:] - switch_idx_arr[:-1]
    print 'asdfdas', switch_idx_arr2


    # create array  [1,0,1,0,......]
    aran = arange( 0, len( switch_idx_arr2 ) )
    mask = aran % 2 == 0
    aran[mask] = True
    aran[invert( mask )] = False
    # repeat values
    aran = ( aran.repeat( switch_idx_arr2 ) ).astype( bool )
    print aran

    a = arange( min( switch_idx_arr ), max( switch_idx_arr ) )
    print a
    print 'fas', a[aran]
    col = a[aran] - row * data.shape[1]

    print csr_matrix( ( data_row, ( row, col ) ), shape = data_bas.shape ).todense()
    print data_bas

i_ident_2d()
