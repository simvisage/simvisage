'''
Created on Aug 23, 2012

@author: rch
'''

from etsproxy.traits.api import \
    HasStrictTraits, List, Float, Property, \
    WeakRef, cached_property, Str, \
    Instance, Button, Event

from util.traits.editors.mpl_figure_editor import  \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from matplotlib.ticker import AutoMinorLocator

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, ModelView, VGroup, HGroup, RangeEditor, InstanceEditor

import numpy as np

from math import exp, log

from mathkit.mfn import MFnLineArray

class CLBase(HasStrictTraits):
    '''Base class for Effective Crack Bridge Laws.'''

    cs = WeakRef

    def __init__(self, *args, **kw):
        super(HasStrictTraits, self).__init__(*args, **kw)
        self.on_trait_change(self._notify_cs, '+input')

    def _notify_cs(self):
        if self.cs:
            self.cs.set_modified()

    def set_cparams(self, *args):
        for name, value in zip(self.cnames, args):
            setattr(self, name, value)
        self._notify_cs()

    arr = Property()
    def _get_arr(self):
        return self.eps_arr, self.sig_arr

    mfn = Property()
    def _get_mfn(self):
        return MFnLineArray(xdata = self.eps_arr,
                            ydata = self.sig_arr)

    mfn_vct = Property()
    def _get_mfn_vct(self):
        return np.vectorize(self.mfn.get_value, otypes = [np.float])

    def plot(self, ax, **kw):
        ax.plot(*self.arr, **kw)
#        ax.autoscale(tight = True)
#        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#        ax.tick_params(direction = 'out', length = 2, width = 4)
        ax.ticklabel_format(axis = 'x', style = 'sci')

    def default_traits_view(self):

        input_traits = self.traits(input = lambda x: x != None)

        citems = [Item(name) for name in input_traits ]
        return View(*citems,
                    buttons = ['OK', 'Cancel']
                    )

class ConstitutiveLawModelView(ModelView):

    model = Instance(CLBase)

    data_changed = Event

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor = 'white')
        return figure

    replot = Button()
    def _replot_fired(self):
        ax = self.figure.add_subplot(1, 1, 1)
        self.model.plot(ax)
        self.data_changed = True

    clear = Button()
    def _clear_fired(self):
        self.figure.clear()
        self.data_changed = True

    traits_view = View(HSplit(
                Group(
                      Item('model', style = 'custom', show_label = False, resizable = True),
                      scrollable = True,
                      ),
                Group(HGroup(
                             Item('replot', show_label = False),
                             Item('clear', show_label = False),
                      ),
                      Item('figure', editor = MPLFigureEditor(),
                           resizable = True, show_label = False),
                      id = 'simexdb.plot_sheet',
                      label = 'plot sheet',
                      dock = 'tab',
                      ),
                       ),
                width = 0.5,
                height = 0.4,
                resizable = True,
                buttons = ['OK', 'Cancel'])

