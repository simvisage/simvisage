'''
Created on Dec 12, 2010

@author: kelidas
'''
from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Tuple, Interface, implements, Enum, Str, List
from enthought.traits.trait_types import DelegatesTo
from enthought.traits.ui.api import Item, View, HGroup, RangeEditor, EnumEditor
import unittest
from numpy import min, array, all, sum, round
from os.path import join
from promod.simdb import SimDB
from StringIO import StringIO
from ymb_data import YMBCutData, YMBSource, yarn_list
from ymb_micro import YMBMicro
import numpy.ma as ma

simdb = SimDB()
test_dir = join( '__TEST__', 'raw_data' )
yarn = '__TEST__'
PX2MUM = 1.2

class S( YMBSource ):

    yarn_type = Enum( '__TEST__', changed_source = True , editor = EnumEditor( values = ['__TEST__'] ) )
    root_dir = Str( '' )

    view = View( Item( 'yarn_type' ) )


class TestYarnSource_0( unittest.TestCase ):
    def setUp( self ):
        print 'setting up'
        source = S( yarn_type = yarn, root_dir = '' )
        self.yarn_data = YMBCutData( source = source, cf_limit = 0.0 )

    def test_cut_data( self ):
        mask = array( [[False, False, False, False, False, False, False, False, False, False],
                [False, False, False, False, False, False, False, False, False, False],
                [False, False, False, False, False, False, False, True, True, True],
                [False, False, False, False, False, False, False, False, False, False],
                [False, False, False, False, False, False, True, True, True, True]], dtype = bool )
        x_true = ma.masked_array( [[0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 4.0],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 4.0],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, -1, -1, -1],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 4.0],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, -1, -1, -1, -1]],
                                    mask = mask, fill_value = -1.0 )

        y_true = ma.masked_array( [[0.05, 0.07, 0.05, 0.01, 0.04, 0.09, 0.04, 0.0, 0.05, 0.09],
                                   [0.05, 0.08, 0.06, 0.1, 0.05, 0.07, 0.04, 0.07, 0.09, 0.12],
                                   [0.25, 0.28, 0.26, 0.24, 0.26, 0.28, 0.3, -1, -1, -1],
                                   [0.15, 0.17, 0.15, 0.13, 0.15, 0.17, 0.15, 0.13, 0.15, 0.17],
                                   [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, -1, -1, -1, -1]],
                                   mask = mask, fill_value = -1.0 )

        z_true = ma.masked_array( [[0.01, 0.03, 0.06, 0.01, 0.0, 0.05, 0.1, 0.06, 0.01, 0.05],
                                   [0.21, 0.18, 0.15, 0.19, 0.21, 0.18, 0.15, 0.18, 0.21, 0.18],
                                   [0.01, 0.04, 0.07, 0.04, 0.01, 0.04, 0.06, -1, -1, -1],
                                   [0.11, 0.13, 0.15, 0.13, 0.11, 0.13, 0.15, 0.13, 0.11, 0.13],
                                   [0.21, 0.21, 0.21, 0.21, 0.21, 0.21, -1, -1, -1, -1]],
                                   mask = mask, fill_value = -1.0 )


        cf_true = ma.masked_array( [[1.0, 0.75, 0.5, 0.25, 0.0, 0.111111111111, 0.25, 0.5, 0.75, 1.0],
                               [1.0, 0.75, 0.25, 0.5, 0.0, 0.5, 0.25, 0.5, 0.25, 0.0],
                               [0.0, 0.25, 0.5, 0.75, 1.0, 0.75, 0.5, -1, -1, -1],
                               [0.0, 0.111111111111, 0.25, 0.361111111111, 0.5, 0.611111111111, 0.75, 0.861111111111, 1.0, 1.0],
                               [0.5, 0.25, 0.0, 0.0, 0.25, 0.75, -1, -1, -1, -1]],
                               mask = mask, fill_value = -1.0 )
        bfl_x_true_0 = ma.masked_array( [[0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                                        [0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.6],
                                        [0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, -1],
                                        [0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                        [0.0, 1.3, 1.3, 1.3, 0.0, -1, -1, -1, -1]],
                                        mask = mask[:, 1:], fill_value = -1.0 )

        bfl_true_0 = ma.masked_array( [[0.0, 0.0, 0.0, 1.00597424891, 1.00597424891, 0.0, 0.0, 0.0, 0.0],
                                       [0.0, 0.0, 0.0, 1.00418995281, 1.00418995281, 0.0, 0.0, 0.0, 0.601498129673],
                                       [0.108627804912, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, -1],
                                       [0.103923048454, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                       [0.0, 1.3, 1.3, 1.3, 0.0, -1, -1, -1, -1]],
                                       mask = mask[:, 1:], fill_value = -1.0 )

        slack = array( [[0.0, 0.0, 0.0, 0.00199011446037, 0.00199011446037, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.00369191553907, 0.00369191553907, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, -1 ],
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, -1, -1 ]] )
        slack_true = ma.masked_array( slack,
                                      mask = mask[:, 1:], fill_value = -1.0 )

        length = ma.masked_array( 
                                [[0.103923048454, 0.302158898595, 0.504083326445, 0.500999001995, 0.504975246918, 0.504975246918, 0.503189825016, 0.504975246918, 0.602660766933],
                                 [0.108627804912, 0.302158898595, 0.503189825016, 0.502891638427, 0.50129831438, 0.501796771612, 0.501796771612, 0.50129831438, 0.601498129673],
                                 [0.108627804912, 0.302158898595, 0.50129831438, 0.50129831438, 0.50129831438, 0.500799361022, -1, -1, -1],
                                 [0.103923048454, 0.301330383466, 0.500799361022, 0.500799361022, 0.500799361022, 0.500799361022, 0.500799361022, 0.500799361022, 0.600666296707],
                                 [0.1, 0.3, 0.5, 0.5, 0.5, -1, -1, -1, -1]],
                                 mask = mask[:, 1:], fill_value = -1.0 )

        print 'x_array test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[0] - x_true ), decimals = 10 ) == 0, True )
        print 'y_array test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[1] - y_true ), decimals = 10 ) == 0, True )
        print 'z_array test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[2] - z_true ), decimals = 10 ) == 0, True )
        print 'contact fraction test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[5] - cf_true ), decimals = 10 ) == 0, True )
        print 'segment bond free length x test'
        self.assertEqual( round( sum( self.yarn_data.fs_bond_free_length_x - bfl_x_true_0 ), decimals = 10 ) == 0, True )
        print 'segment bond free length test'
        self.assertEqual( round( sum( self.yarn_data.fs_bond_free_length - bfl_true_0 ), decimals = 10 ) == 0, True )
        print 'segment length between cuts test'
        self.assertEqual( round( sum( self.yarn_data.fs_length_between_cuts - length ), decimals = 10 ) == 0, True )
        print 'slack test'
        self.assertEqual( round( sum( self.yarn_data.fs_slack - slack_true ), decimals = 10 ) == 0, True )


        print 'yarn cf = 0 tested'


class TestYarnSource_50( unittest.TestCase ):
    def setUp( self ):
        print 'setting up'
        self.yarn_data = YMBCutData( source = S( yarn_type = yarn ), cf_limit = 0.5 )

    def test_cut_data( self ):
        mask = array( [[False, False, False, False, False, False, False, False, False, False],
                [False, False, False, False, False, False, False, False, False, False],
                [False, False, False, False, False, False, False, True, True, True],
                [False, False, False, False, False, False, False, False, False, False],
                [False, False, False, False, False, False, True, True, True, True]], dtype = bool )
        x_true = ma.masked_array( [[0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 4.0],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 4.0],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, -1, -1, -1],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 4.0],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, -1, -1, -1, -1]],
                                    mask = mask, fill_value = -1.0 )

        y_true = ma.masked_array( [[0.05, 0.07, 0.05, 0.01, 0.04, 0.09, 0.04, 0.0, 0.05, 0.09],
                                   [0.05, 0.08, 0.06, 0.1, 0.05, 0.07, 0.04, 0.07, 0.09, 0.12],
                                   [0.25, 0.28, 0.26, 0.24, 0.26, 0.28, 0.3, -1, -1, -1],
                                   [0.15, 0.17, 0.15, 0.13, 0.15, 0.17, 0.15, 0.13, 0.15, 0.17],
                                   [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, -1, -1, -1, -1]],
                                   mask = mask, fill_value = -1.0 )

        z_true = ma.masked_array( [[0.01, 0.03, 0.06, 0.01, 0.0, 0.05, 0.1, 0.06, 0.01, 0.05],
                                   [0.21, 0.18, 0.15, 0.19, 0.21, 0.18, 0.15, 0.18, 0.21, 0.18],
                                   [0.01, 0.04, 0.07, 0.04, 0.01, 0.04, 0.06, -1, -1, -1],
                                   [0.11, 0.13, 0.15, 0.13, 0.11, 0.13, 0.15, 0.13, 0.11, 0.13],
                                   [0.21, 0.21, 0.21, 0.21, 0.21, 0.21, -1, -1, -1, -1]],
                                   mask = mask, fill_value = -1.0 )


        cf_true = ma.masked_array( [[1.0, 0.75, 0.5, 0.25, 0.0, 0.111111111111, 0.25, 0.5, 0.75, 1.0],
                               [1.0, 0.75, 0.25, 0.5, 0.0, 0.5, 0.25, 0.5, 0.25, 0.0],
                               [0.0, 0.25, 0.5, 0.75, 1.0, 0.75, 0.5, -1, -1, -1],
                               [0.0, 0.111111111111, 0.25, 0.361111111111, 0.5, 0.611111111111, 0.75, 0.861111111111, 1.0, 1.0],
                               [0.5, 0.25, 0.0, 0.0, 0.25, 0.75, -1, -1, -1, -1]],
                               mask = mask, fill_value = -1.0 )
        bfl_x_true_50 = ma.masked_array( [[0.0, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 0.0],
                                          [0.0, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9],
                                          [0.9, 0.9, 0.9, 0.0, 0.0, 0.5, -1, -1, -1],
                                          [1.9, 1.9, 1.9, 1.9, 1.9, 0.0, 0.0, 0.0, 0.0],
                                          [1.9, 1.9, 1.9, 1.9, 1.9, -1, -1, -1, -1]],
                                        mask = mask[:, 1:], fill_value = -1.0 )

        bfl_true_50 = ma.masked_array( [[0.0, 3.32535679281, 3.32535679281, 3.32535679281, 3.32535679281, 3.32535679281, 3.32535679281, 3.32535679281, 0.0],
                                        [0.0, 3.91592866369, 3.91592866369, 3.91592866369, 3.91592866369, 3.91592866369, 3.91592866369, 3.91592866369, 3.91592866369],
                                        [0.912085017887, 0.912085017887, 0.912085017887, 0.0, 0.0, 0.500799361022, -1, -1, -1],
                                        [1.90765151499, 1.90765151499, 1.90765151499, 1.90765151499, 1.90765151499, 0.0, 0.0, 0.0, 0.0],
                                        [1.9, 1.9, 1.9, 1.9, 1.9, -1, -1, -1, -1]],
                                       mask = mask[:, 1:], fill_value = -1.0 )

        slack = array( [[0.0, 0.00764742508657, 0.00764742508657, 0.00764742508657, 0.00764742508657, 0.00764742508657, 0.00764742508657, 0.00764742508657, 0.0],
                         [0.0, 0.00403188834656, 0.00403188834656, 0.00403188834656, 0.00403188834656, 0.00403188834656, 0.00403188834656, 0.00403188834656, 0.00403188834656],
                         [0.0128185094833, 0.0128185094833, 0.0128185094833, 0.0, 0.0, -1.11022302463e-16, -1, -1, -1],
                         [0.00391674818861, 0.00391674818861, 0.00391674818861, 0.00391674818861, 0.00391674818861, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, -1, -1]] )
        slack_true = ma.masked_array( slack,
                                      mask = mask[:, 1:], fill_value = -1.0 )

        length = ma.masked_array( 
                                [[0.103923048454, 0.302158898595, 0.504083326445, 0.500999001995, 0.504975246918, 0.504975246918, 0.503189825016, 0.504975246918, 0.602660766933],
                                 [0.108627804912, 0.302158898595, 0.503189825016, 0.502891638427, 0.50129831438, 0.501796771612, 0.501796771612, 0.50129831438, 0.601498129673],
                                 [0.108627804912, 0.302158898595, 0.50129831438, 0.50129831438, 0.50129831438, 0.500799361022, -1, -1, -1],
                                 [0.103923048454, 0.301330383466, 0.500799361022, 0.500799361022, 0.500799361022, 0.500799361022, 0.500799361022, 0.500799361022, 0.600666296707],
                                 [0.1, 0.3, 0.5, 0.5, 0.5, -1, -1, -1, -1]],
                                 mask = mask[:, 1:], fill_value = -1.0 )

        print 'x_array test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[0] - x_true ), decimals = 10 ) == 0, True )
        print 'y_array test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[1] - y_true ), decimals = 10 ) == 0, True )
        print 'z_array test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[2] - z_true ), decimals = 10 ) == 0, True )
        print 'contact fraction test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[5] - cf_true ), decimals = 10 ) == 0, True )
        print 'segment bond free length x test'
        self.assertEqual( round( sum( self.yarn_data.fs_bond_free_length_x - bfl_x_true_50 ), decimals = 10 ) == 0, True )
        print 'segment bond free length test'
        self.assertEqual( round( sum( self.yarn_data.fs_bond_free_length - bfl_true_50 ), decimals = 10 ) == 0, True )
        print 'segment length between cuts test'
        self.assertEqual( round( sum( self.yarn_data.fs_length_between_cuts - length ), decimals = 10 ) == 0, True )
        print 'slack test'
        self.assertEqual( round( sum( self.yarn_data.fs_slack - slack_true ), decimals = 10 ) == 0, True )


        print 'yarn cf = 50 tested'


class TestYarnSource_100( unittest.TestCase ):
    def setUp( self ):
        print 'setting up'
        self.yarn_data = YMBCutData( source = S( yarn_type = yarn ), cf_limit = 1.0 )

    def test_cut_data( self ):
        mask = array( [[False, False, False, False, False, False, False, False, False, False],
                [False, False, False, False, False, False, False, False, False, False],
                [False, False, False, False, False, False, False, True, True, True],
                [False, False, False, False, False, False, False, False, False, False],
                [False, False, False, False, False, False, True, True, True, True]], dtype = bool )
        x_true = ma.masked_array( [[0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 4.0],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 4.0],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, -1, -1, -1],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4, 4.0],
                                    [0.0, 0.1, 0.4, 0.9, 1.4, 1.9, -1, -1, -1, -1]],
                                    mask = mask, fill_value = -1.0 )

        y_true = ma.masked_array( [[0.05, 0.07, 0.05, 0.01, 0.04, 0.09, 0.04, 0.0, 0.05, 0.09],
                                   [0.05, 0.08, 0.06, 0.1, 0.05, 0.07, 0.04, 0.07, 0.09, 0.12],
                                   [0.25, 0.28, 0.26, 0.24, 0.26, 0.28, 0.3, -1, -1, -1],
                                   [0.15, 0.17, 0.15, 0.13, 0.15, 0.17, 0.15, 0.13, 0.15, 0.17],
                                   [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, -1, -1, -1, -1]],
                                   mask = mask, fill_value = -1.0 )

        z_true = ma.masked_array( [[0.01, 0.03, 0.06, 0.01, 0.0, 0.05, 0.1, 0.06, 0.01, 0.05],
                                   [0.21, 0.18, 0.15, 0.19, 0.21, 0.18, 0.15, 0.18, 0.21, 0.18],
                                   [0.01, 0.04, 0.07, 0.04, 0.01, 0.04, 0.06, -1, -1, -1],
                                   [0.11, 0.13, 0.15, 0.13, 0.11, 0.13, 0.15, 0.13, 0.11, 0.13],
                                   [0.21, 0.21, 0.21, 0.21, 0.21, 0.21, -1, -1, -1, -1]],
                                   mask = mask, fill_value = -1.0 )


        cf_true = ma.masked_array( [[1.0, 0.75, 0.5, 0.25, 0.0, 0.111111111111, 0.25, 0.5, 0.75, 1.0],
                               [1.0, 0.75, 0.25, 0.5, 0.0, 0.5, 0.25, 0.5, 0.25, 0.0],
                               [0.0, 0.25, 0.5, 0.75, 1.0, 0.75, 0.5, -1, -1, -1],
                               [0.0, 0.111111111111, 0.25, 0.361111111111, 0.5, 0.611111111111, 0.75, 0.861111111111, 1.0, 1.0],
                               [0.5, 0.25, 0.0, 0.0, 0.25, 0.75, -1, -1, -1, -1]],
                               mask = mask, fill_value = -1.0 )
        bfl_x_true_100 = ma.masked_array( [[4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
                                           [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
                                           [2.4, 2.4, 2.4, 2.4, 2.4, 2.4, -1, -1, -1],
                                           [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
                                           [1.9, 1.9, 1.9, 1.9, 1.9, -1, -1, -1, -1]],
                                        mask = mask[:, 1:], fill_value = -1.0 )

        bfl_true_100 = ma.masked_array( [[4.03194060819, 4.03194060819, 4.03194060819, 4.03194060819, 4.03194060819, 4.03194060819, 4.03194060819, 4.03194060819, 4.03194060819],
                                         [4.02455646861, 4.02455646861, 4.02455646861, 4.02455646861, 4.02455646861, 4.02455646861, 4.02455646861, 4.02455646861, 4.02455646861],
                                         [2.41548100767, 2.41548100767, 2.41548100767, 2.41548100767, 2.41548100767, 2.41548100767, -1, -1, -1],
                                         [4.01071589476, 4.01071589476, 4.01071589476, 4.01071589476, 4.01071589476, 4.01071589476, 4.01071589476, 4.01071589476, 4.01071589476],
                                         [1.9, 1.9, 1.9, 1.9, 1.9, -1, -1, -1, -1]],
                                       mask = mask[:, 1:], fill_value = -1.0 )

        slack = array( [[0.00788593935499, 0.00788593935499, 0.00788593935499, 0.00788593935499, 0.00788593935499, 0.00788593935499, 0.00788593935499, 0.00788593935499, 0.00788593935499],
                         [0.00595895694733, 0.00595895694733, 0.00595895694733, 0.00595895694733, 0.00595895694733, 0.00595895694733, 0.00595895694733, 0.00595895694733, 0.00595895694733],
                         [0.00601908137182, 0.00601908137182, 0.00601908137182, 0.00601908137182, 0.00601908137182, 0.00601908137182, -1, -1, -1],
                         [0.00265404017527, 0.00265404017527, 0.00265404017527, 0.00265404017527, 0.00265404017527, 0.00265404017527, 0.00265404017527, 0.00265404017527, 0.00265404017527],
                         [0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, -1, -1]]
                      )
        slack_true = ma.masked_array( slack,
                                      mask = mask[:, 1:], fill_value = -1.0 )

        length = ma.masked_array( 
                                [[0.103923048454, 0.302158898595, 0.504083326445, 0.500999001995, 0.504975246918, 0.504975246918, 0.503189825016, 0.504975246918, 0.602660766933],
                                 [0.108627804912, 0.302158898595, 0.503189825016, 0.502891638427, 0.50129831438, 0.501796771612, 0.501796771612, 0.50129831438, 0.601498129673],
                                 [0.108627804912, 0.302158898595, 0.50129831438, 0.50129831438, 0.50129831438, 0.500799361022, -1, -1, -1],
                                 [0.103923048454, 0.301330383466, 0.500799361022, 0.500799361022, 0.500799361022, 0.500799361022, 0.500799361022, 0.500799361022, 0.600666296707],
                                 [0.1, 0.3, 0.5, 0.5, 0.5, -1, -1, -1, -1]],
                                 mask = mask[:, 1:], fill_value = -1.0 )

        print 'x_array test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[0] - x_true ), decimals = 10 ) == 0, True )
        print 'y_array test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[1] - y_true ), decimals = 10 ) == 0, True )
        print 'z_array test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[2] - z_true ), decimals = 10 ) == 0, True )
        print 'contact fraction test'
        self.assertEqual( round( sum( self.yarn_data.cut_data[5] - cf_true ), decimals = 10 ) == 0, True )
        print 'segment bond free length x test'
        self.assertEqual( round( sum( self.yarn_data.fs_bond_free_length_x - bfl_x_true_100 ), decimals = 10 ) == 0, True )
        print 'segment bond free length test'
        self.assertEqual( round( sum( self.yarn_data.fs_bond_free_length - bfl_true_100 ), decimals = 10 ) == 0, True )
        print 'segment length between cuts test'
        self.assertEqual( round( sum( self.yarn_data.fs_length_between_cuts - length ), decimals = 10 ) == 0, True )
        print 'slack test'
        self.assertEqual( round( sum( self.yarn_data.fs_slack - slack_true ), decimals = 10 ) == 0, True )


        print 'yarn cf = 100 tested'

if __name__ == '__main__':
    unittest.main()

