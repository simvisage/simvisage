'''
Created on 28. 3. 2014

The module defines the framework for tree visualization of the model
in a ModelView window. The components of the model classes should
inherit from the BMCSTreeNode and supply the attributes

 - node_name - a label to appear in the tree editor
 - tree_node_list  - subnodes of the current nodes

 Further, the plot method can be defined to plot the current node
 in the plot window of the model view.

@author: rch
'''

from matplotlib.figure import \
    Figure
from traits.api import \
    HasStrictTraits, Instance, Button, Event, Range, cached_property 
from traits.etsconfig.api import ETSConfig
from traitsui.api import \
    TreeEditor, TreeNode, View, Item, Group, \
    HSplit, HGroup
from traitsui.menu import \
    Menu, MenuBar, Separator

from bmcs_tree_view_handler import \
    BMCSTreeViewHandler, plot_self, menu_save, \
    menu_open, menu_exit
from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor
from view.ui.bmcs_tree_node import \
    BMCSTreeNode, BMCSLeafNode


if ETSConfig.toolkit == 'wx':
    from traitsui.wx.tree_editor import \
        DeleteAction
if ETSConfig.toolkit == 'qt4':
    from traitsui.qt4.tree_editor import \
        DeleteAction
else:
    raise ImportError, "tree actions for %s toolkit not availabe" % \
        ETSConfig.toolkit

tree_node = TreeNode(node_for=[BMCSTreeNode],
                     auto_open=False,
                     children='tree_node_list',
                     label='node_name',
                     view='tree_view',
                     menu=Menu(plot_self, DeleteAction),
                     )

leaf_node = TreeNode(node_for=[BMCSLeafNode],
                     auto_open=True,
                     children='',
                     label='node_name',
                     view='tree_view',
                     menu=Menu(plot_self)
                     )

tree_editor = TreeEditor(
    nodes=[tree_node, leaf_node],
    selected='selected_node',
    orientation='vertical'
)


class BMCSWindow(HasStrictTraits):

    '''View object for a cross section state.
    '''
    root = Instance(BMCSTreeNode)

    selected_node = Instance(HasStrictTraits)

    figure = Instance(Figure)

    def _figure_default(self):
        figure = Figure()
        return figure

    data_changed = Event

    replot = Button

    def _replot_fired(self):
        self.figure.clear()
        self.selected_node.plot(self.figure)
        self.data_changed = True

    clear = Button()

    def _clear_fired(self):
        self.figure.clear()
        self.data_changed = True
        
    #time = self.root.time

    view = View(HSplit(Group(Item('root',
                                  editor=tree_editor,
                                  resizable=True,
                                  show_label=False,
                                  width=400,
                                  height=400),
                             ),
                       Group(HGroup(Item('replot', show_label=False),
                                    Item('clear', show_label=False)
                                    ),
                             Item('figure', editor=MPLFigureEditor(),
                                  resizable=True, show_label=False,
                                  springy=True),#Item('self.root.time', label='t/T_max'),
                             label='plot sheet',
                             dock='tab',
                             )
                       ),
                id='bmcstreeview_id',
                width=0.9,
                height=0.5,
                title='BMCS',
                resizable=True,
                handler=BMCSTreeViewHandler(),
                menubar=MenuBar(Menu(menu_exit, Separator(),
                                     menu_save, menu_open,
                                     name='File'))
                )

if __name__ == '__main__':

    tr = BMCSTreeNode(node_name='root',
                      tree_node_list=[BMCSTreeNode(node_name='subnode 1'),
                                      BMCSTreeNode(node_name='subnode 2'),
                                      ])

    tv = BMCSWindow(root=tr)
    tv.configure_traits()
