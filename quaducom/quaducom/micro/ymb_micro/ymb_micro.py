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

#import wxversion
#wxversion.select( '2.8' )

from etsproxy.traits.api import \
    HasTraits, Instance, Property, cached_property, List, WeakRef

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, VSplit, VGroup, Tabbed, \
    TreeEditor, TreeNode, Handler, ObjectTreeNode

from ymb_data import YMBSource, IYMBData, YMBCutData, YMBSegmentData
from ymb_data import YMBSlider
from ymb_hist import YMBHist
from ymb_auto_correl import YMBAutoCorrelView, YMBAutoCorrel
from ymb_cross_correl import YMBCrossCorrel
from ymb_view3d import YMBView3D
from ymb_view2d import YMBView2D
from ymb_pullout import YMBPullOut
from ymb_pdistrib import YMBDistrib

import os
#os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'


class TitleHandler( Handler ):
    """ Change the title on the UI.
    """

    def object_yarn_type_changed( self, info ):
        """ Called whenever the "yarn" attribute changes on the handled
        object.
        """
        info.ui.title = '%s -- Yarn-Matrix-Bond Microstructure Lab' % info.object.data.source.yarn_type

class YMBNodeTree( HasTraits ):

    ymb_micro = WeakRef
    children = 'node_list'
    # empty view
    traits_view = View()

class YMBInspector( YMBNodeTree ):
    node_list = Property( List )
    @cached_property
    def _get_node_list( self ):
        return [ self.ymb_micro.view_3d, self.ymb_micro.cut_view ]

class YMBStatistics( YMBNodeTree ):
    node_list = Property( List )
    @cached_property
    def _get_node_list( self ):
        return [ self.ymb_micro.histogram,
                self.ymb_micro.auto_correl,
                self.ymb_micro.cross_correl
                ]

class YMBModel( YMBNodeTree ):
    node_list = Property( List )
    @cached_property
    def _get_node_list( self ):
        return [ self.ymb_micro.pullout ]

class YMBTree( YMBNodeTree ):
    node_list = Property( List )
    @cached_property
    def _get_node_list( self ):
        return [ YMBInspector( ymb_micro = self.ymb_micro ),
                 YMBStatistics( ymb_micro = self.ymb_micro ),
                 YMBModel( ymb_micro = self.ymb_micro ) ]

tree_editor = TreeEditor( 
    nodes = [
        TreeNode( node_for = [YMBTree],
                  auto_open = True,
                  label = '=root',
                  children = 'node_list' ),
        TreeNode( node_for = [YMBInspector],
                  auto_open = True,
                  label = '=Data inspection',
                  children = 'node_list' ),
        TreeNode( node_for = [ YMBStatistics ],
                  auto_open = True,
                  label = '=Statistics',
                  children = 'node_list',
                  ),
        TreeNode( node_for = [ YMBModel ],
                  auto_open = True,
                  label = '=Model',
                  children = 'node_list',
                  ),
        TreeNode( node_for = [ YMBPullOut ],
                  auto_open = False,
                  label = '=Pull-out',
                  children = '',
                  ),
        TreeNode( node_for = [ YMBView3D ],
                  auto_open = False,
                  label = '=view 3D',
                  children = '',
                  ),
        TreeNode( node_for = [ YMBView2D ],
                  auto_open = False,
                  label = '=view 2D',
                  children = '',
                  ),
        TreeNode( node_for = [ YMBHist ],
                  auto_open = False,
                  label = '=Histogram',
                  children = '',
                  ),
        TreeNode( node_for = [ YMBAutoCorrelView ],
                  auto_open = False,
                  label = '=Auto correlation',
                  children = '',
                  ),
        TreeNode( node_for = [ YMBCrossCorrel ],
                  auto_open = False,
                  label = '=Cross correlation',
                  children = '',
                  ),
                 ],
    hide_root = True
    )


class YMBMicro( HasTraits ):
    '''
    Yarn-Matrix-Bond Microstructure analysis
    '''

    data = Instance( IYMBData )

    view_3d = Property( Instance( YMBView3D ) )
    @cached_property
    def _get_view_3d( self ):
        return YMBView3D( data = self.data )

    cut_view = Property( Instance( YMBView2D ) )
    @cached_property
    def _get_cut_view( self ):
        return YMBView2D( data = self.data )

    slider = Property( Instance( YMBSlider ), depends_on = 'data.input_change' )
    @cached_property
    def _get_slider( self ):
        return YMBSlider( data = self.data )

    histogram = Property( Instance( YMBHist ), depends_on = 'data.input_change' )
    @cached_property
    def _get_histogram( self ):
        return YMBHist( slider = self.slider )

    auto_correl_data = Property( Instance( YMBAutoCorrel ), depends_on = 'data.input_change' )
    @cached_property
    def _get_auto_correl_data( self ):
        return YMBAutoCorrel( data = self.data )

    auto_correl = Property( Instance( YMBAutoCorrelView ), depends_on = 'data.input_change' )
    @cached_property
    def _get_auto_correl( self ):
        return YMBAutoCorrelView( correl_data = self.auto_correl_data )

    cross_correl = Property( Instance( YMBCrossCorrel ), depends_on = 'data.input_change' )
    @cached_property
    def _get_cross_correl( self ):
        return YMBCrossCorrel( data = self.data )

    pullout = Instance( YMBPullOut )
    def _pullout_default( self ):
        return YMBPullOut( data = self.data )

    ymb_tree = Property( Instance( YMBTree ) )
    @cached_property
    def _get_ymb_tree( self ):
        return YMBTree( ymb_micro = self )


    # The main view
    view = View( VGroup( 
               Group( 
                     Item( 'data@', show_label = False ),
                     label = 'Data source'
                     ),
               Group( 
                   Item( 
                        name = 'ymb_tree',
                        id = 'ymb_tree',
                        editor = tree_editor,
                        show_label = False,
                        resizable = True ),
                    orientation = 'vertical',
                    show_labels = True,
                    show_left = True, ),
                   dock = 'tab',
                   id = 'qdc.ymb.split',
                ),
                id = 'qdc.ymb',
                dock = 'horizontal',
                drop_class = HasTraits,
    #            handler = TreeHandler(),
                resizable = True,
                scrollable = True,
                title = 'Yarn-Matrix-Bond Microstructure Lab',
                handler = TitleHandler(),
                width = .8,
                height = .5 )


if __name__ == '__main__':
    from os.path import join
    from matresdev.db.simdb import SimDB
    simdb = SimDB()
    yarn_type = 'MAG'
    data = YMBCutData( source = YMBSource( yarn_type = yarn_type ), cf_limit = 0.5 )
    ymbmicro = YMBMicro( data = data )
    ymbmicro.configure_traits()
