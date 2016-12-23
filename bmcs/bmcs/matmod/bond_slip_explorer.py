'''
Created on 9. 7. 2014

@author: Vancikv

This module contains tree nodes defined for specific classes
as oposed to the two generic node types defined in mxn_tree_view
that are sufficient for displaying the tree structure but
don't provide different classes with specific node behaviour.
'''

import os
import pickle

from traits.etsconfig.api import ETSConfig
from traitsui.api import \
    TreeNode
from traitsui.menu import \
    Menu

from bond_slip_model import \
    Material, LoadingScenario, BondSlipModel

from view.window import BMCSWindow
from view.window.bmcs_tree_view_handler import \
    plot_self, new_material, del_material


if ETSConfig.toolkit == 'wx':
    from traitsui.wx.tree_editor import \
        NewAction, DeleteAction, CopyAction, PasteAction
if ETSConfig.toolkit == 'qt4':
    from traitsui.qt4.tree_editor import \
        NewAction, DeleteAction, CopyAction, PasteAction
else:
    raise ImportError, "tree actions for %s toolkit not availabe" % \
        ETSConfig.toolkit


# =========================================================================
# Special TreeNode classes
# =========================================================================


material_node = TreeNode(node_for=[Material],
                         auto_open=False,
                         children='tree_node_list',
                         label='node_name',
                         view='tree_view',
                         menu=Menu(del_material),
                         )

loading_scenario_node = TreeNode(node_for=[LoadingScenario],
                                 auto_open=True,
                                 children='tree_node_list',
                                 label='node_name',
                                 view='tree_view',
                                 menu=Menu(CopyAction, plot_self),
                                 )

bond_slip_model_node = TreeNode(node_for=[BondSlipModel],
                                auto_open=True,
                                children='tree_node_list',
                                label='node_name',
                                view='tree_view',
                                menu=Menu(plot_self, NewAction),
                                )

# =========================================================================
# List of all custom nodes
# =========================================================================

custom_node_list = [material_node, loading_scenario_node,
                    bond_slip_model_node]

model = BondSlipModel()
w = BMCSWindow(root=model)
w.configure_traits()