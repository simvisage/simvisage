
from traits.api import \
    HasTraits, \
    Instance,  \
    List

from traitsui.api import \
    View, Item, VSplit, \
    TableEditor, ObjectColumn

from i_bcond import IBCond


# The definition of the demo TableEditor:
bcond_list_editor = TableEditor(
    columns=[ObjectColumn(label='Type', name='var'),
             ObjectColumn(label='Value', name='value'),
             ObjectColumn(label='DOF', name='dof')
             ],
    editable=False,
    selected='object.selected_bcond',
)


class BCondMngr(HasTraits):

    bcond_list = List(IBCond)

    selected_bcond = Instance(IBCond)

    def setup(self, sctx):
        '''
        '''
        for bc in self.bcond_list:
            bc.setup(sctx)

    def apply_essential(self, K):
        '''Register the boundary condition in the equation system.
        '''
        for bcond in self.bcond_list:
            bcond.apply_essential(K)

    def apply(self, step_flag, sctx, K, R, t_n, t_n1):

        for bcond in self.bcond_list:
            bcond.apply(step_flag, sctx, K, R, t_n, t_n1)

    traits_view = View(VSplit(Item('bcond_list', style='custom', editor=bcond_list_editor,
                                   show_label=False),
                              Item('selected_bcond@', show_label=False)),
                       resizable=True,
                       kind='subpanel',
                       )
