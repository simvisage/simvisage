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

from numpy import array, where, nonzero, array, hstack


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
    data = array( [[ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1 ],
                   [ -1, -1, 1, 1, 1, -1, -1, 1, 1, -1, -1, 1 ],
                   [ -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]] )
    #data = random( ( 2000, 15 ) ) - .5
    data = data / abs( data )
    print 'data', data

    r = range( 0, data.shape[0] )
    switches = where( data[:, 0:-1] * data[:, 1:] < 0 )
    print 'switches', switches
    switch_idx = searchsorted( switches[0], r )
    print 'switch_idx', switch_idx
    spl = split( switches[1], switch_idx[1:] )
    print 'spl', spl
    spl = [hstack( [array( [-1] ), x, [data.shape[1] - 1]] ) for x in spl]

    print data
    print 'switches', spl

    print 'current cut position'
    cut_pos = 5

    print 'find the indices of the interval at 5-th cut'
    idx = [where( x >= cut_pos )[0][0] for x in spl]
    #idx = array( idx )
    print idx

    right_limit = [spl[x][ id ] + 1 for x, id in zip( r, idx )]
    left_limit = [spl[x][ id - 1 ] + 1 for x, id in zip( r, idx )]

    print 'right_limit', right_limit
    print 'left_limit', left_limit

    res = array( [sum( x[ ll: rl ] ) for x, ll, rl in zip( data, left_limit, right_limit )] )
    res[res > 0] = 0
    print 'array of free lengths in segment -%i-' % cut_pos, abs( res )


#i_ident()
i_ident_2d()
