'''
Created on Aug 17, 2011

Copyright (c) 2009, IMB, RWTH Aachen.
All rights reserved.

This software is provided without warranty under the terms of the BSD
license included in simvisage/LICENSE.txt and may be redistributed only
under the conditions described in the aforementioned license.  The license
is also available online at http://www.simvisage.com/licenses/BSD.txt

Thanks for using Simvisage open source!

@author: rostar
'''

from enthought.traits.api import \
    HasTraits, List, Interface, Str, Float

class ICB( Interface ):
    """
    Abstract class representing a homogenized crack bridge.
    As a realization any function class with the get_eps_x_reinf and
    get_P_w members may be included.
    """

    # applied force
    P = Float

    # closest crack from left
    Ll = Float

    # closest crack from right
    Lr = Float

    def get_eps_x_reinf( self ):
        '''
        evaluation of strain profile in the vicinity of a crack bridge
        '''
        pass

    def get_P_w( self ):
        '''
        evaluation of force-crack width relation for a crack bridge
        '''
        pass
