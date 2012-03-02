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
# Created on Dec 20, 2010 by: rch

    #@cached_property
    def _get_p_s( self ):
        dict = {'phi':1.0,
                'ell':0.0,
                'theta':0.0}
        return dict[self.parameter]

    p_c = Property( depends_on = 'parameter, specimen' )
    #@cached_property
    def _get_p_c( self ):
        dict = {'phi' : 0.0,
                'ell' : self.p_c_ell,
                'theta' : self.p_c_theta}
        return dict[self.parameter]

    k_s = Property( depends_on = 'parameter' )
    #@cached_property
    def _get_k_s( self ):
        dict = {'phi':self.k_s_phi,
                'ell':self.k_s_phi,
                'theta':0.0}
        return dict[self.parameter]

    k_c = Property( depends_on = 'parameter, specimen' )
    #@cached_property
    def _get_k_c( self ):
        dict = {'phi' : self.k_c_phi,
                'ell' : self.k_c_phi,
                'theta' : self.k_c_theta}
        return dict[self.parameter]

    k_s_phi = Property( depends_on = 'specimen' )
    #@cached_property
    def _get_k_s_phi( self ):
        dict = {'1':1.308,
                '2':1.612,
                '3':1.492}
        return dict[self.specimen]

    k_c_phi = Property( depends_on = 'specimen' )
    #@cached_property
    def _get_k_c_phi( self ):
        dict = {'1':2.109,
                '2':2.233,
                '3':2.668}
        return dict[self.specimen]

    k_c_theta = Property( depends_on = 'specimen' )
    #@cached_property
    def _get_k_c_theta( self ):
        dict = {'1':3.336,
                '2':4.808,
                '3':2.708}
        return dict[self.specimen]

    p_c_ell = Property( depends_on = 'specimen' )
    #@cached_property
    def _get_p_c_ell( self ):
        dict = {'1' : 9.093,
                '2' : 10.579,
                '3' : 8.664}
        return dict[self.specimen]

    p_c_theta = Property( depends_on = 'specimen' )
    #@cached_property
    def _get_p_c_theta( self ):
        dict = {'1' : 0.04321,
                '2' : 0.04420,
                '3' : 0.00023}
        return dict[self.specimen]

