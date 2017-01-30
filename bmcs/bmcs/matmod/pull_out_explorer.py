

import os
import pickle

from traits.etsconfig.api import ETSConfig
from traitsui.api import \
    TreeNode
from traitsui.menu import \
    Menu

    
from pull_out_simulation import\
 Material, LoadingScenario, Geometry, PullOutSimulation

from view.window import BMCSWindow
from view.window.bmcs_tree_view_handler import \
    plot_self, new_material, del_material
    
from mats_bondslip import MATSEvalFatigue
from fets1d52ulrhfatigue import FETS1D52ULRHFatigue
from tloop import TLoop
from tstepper import TStepper
from ibvpy.api import BCDof



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

geometry_node = TreeNode(node_for=[Geometry],
                                 auto_open=True,
                                 children='tree_node_list',
                                 label='node_name',
                                 view='tree_view',
                                 menu=Menu(CopyAction, plot_self),
                                 )

pull_out_simulation_node = TreeNode(node_for=[PullOutSimulation],
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
                    geometry_node, pull_out_simulation_node]


ts = TStepper()
n_dofs = ts.domain.n_dofs
loading_scenario = LoadingScenario()


ts.bc_list = [BCDof(var='u', dof=0, value=0.0), BCDof(
        var='f', dof=n_dofs - 1, time_function=loading_scenario.time_func)]
tl = TLoop(ts=ts)

#loading_scenario = LoadingScenario()
geometry = Geometry()

model = PullOutSimulation(
        mats_eval=ts.mats_eval, fets_eval=ts.fets_eval
        ,time_stepper=ts, time_loop = tl  ,geometry = geometry,loading_scenario = loading_scenario)
       
w = BMCSWindow(root=model)
w.configure_traits()
