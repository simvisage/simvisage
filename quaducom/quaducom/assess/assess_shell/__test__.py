'''
Created on Aug 8, 2009

@author: alex
'''
import unittest

from etsproxy.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, DelegatesTo, Callable

from etsproxy.util.home_directory import \
    get_home_directory

from etsproxy.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from etsproxy.mayavi import \
    mlab

from etsproxy.traits.ui.table_column import \
    ObjectColumn

from etsproxy.traits.ui.menu import \
    OKButton, CancelButton

from etsproxy.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import \
    array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
    vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
    copy, c_, newaxis, argmax, where, sqrt, frompyfunc, sum, \
    ones, transpose, shape, append, argmin, argmax, fabs

from .ls_table import \
    LSTable, ULS, SLS

from .lcc_table import \
    LCCTable, LCCTableULS, LCCTableSLS, LCC, LC

from math import pi
from string import split
import os

from promod.simdb import \
    SimDB

import os
import pickle
import string
from os.path import join

# Access to the top level directory of the database
#
simdb = SimDB()


class TestCombiArr( unittest.TestCase ):
    '''test for a cycle of different show cases.
    '''

    def setUp( self ):
        '''Test combination ULS/SLS.
        '''
        lc_list = [
            LC( name = 'G1', category = 'dead-load' ),
            LC( name = 'G2', category = 'additional dead-load' ),
            LC( name = 'Q1', category = 'imposed-load',
               exclusive_to = ['Q2', 'Q3'], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.2 ),
            LC( name = 'Q2', category = 'imposed-load',
               exclusive_to = ['Q1', 'Q3'], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.2 ),
            LC( name = 'Q3', category = 'imposed-load',
               exclusive_to = ['Q1', 'Q2'], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.2 ),
            LC( name = 'Q4', category = 'imposed-load',
               exclusive_to = ['Q5'], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.2 ),
            LC( name = 'Q5', category = 'imposed-load',
               exclusive_to = ['Q4'], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.2 ),
            LC( name = 'Q6', category = 'imposed-load',
               exclusive_to = [], psi_0 = 0.7, psi_1 = 0.5, psi_2 = 0.2 ),
          ]

        self.lct_list = [ LCCTableULS( show_lc_characteristic = False, lc_list = lc_list ),
                          LCCTableSLS( show_lc_characteristic = False, lc_list = lc_list, combination_SLS = 'perm' ),
                          LCCTableSLS( show_lc_characteristic = False, lc_list = lc_list, combination_SLS = 'freq' ),
                          LCCTableSLS( show_lc_characteristic = False, lc_list = lc_list, combination_SLS = 'rare' )
                        ]

        self.n_combi_list = [ 188, 24, 47, 47 ]

    def test_combi_arr( self ):
        '''  
        '''
        for i in range( len( self.lct_list ) ):

            # calculate the combi_arr and retrieve its shape:
            #
            lct = self.lct_list[i]
            n_lcc = lct.combi_arr.shape[0]

            # number of combinations for the test 'lc_list' 
            # obtained using InfoGraph:
            #
            n_combi = self.n_combi_list[i]
            print('n_combi', n_combi)

            # check if values are equal:
            #
            self.assertEqual( n_lcc, n_combi )


    def tearDown( self ):
        '''
        '''
        pass

if __name__ == "__main__":
    unittest.main()


