'''
Created on Sep 4, 2012

@author: rch
'''
from etsproxy.traits.api import \
    HasTraits, Int, Instance, Property, cached_property, DelegatesTo, \
    Event, Button

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, VGroup, HGroup, RangeEditor, InstanceEditor

from ecb_law_calib import \
    ECBLCalib

import numpy as np

class ECBLMNDiagram(HasTraits):

    # calibrator supplying the effective material law
    calib = Instance(ECBLCalib)
    def _calib_default(self):
        return ECBLCalib(notify_change=self.set_modified)
    def _calib_changed(self):
        self.calib.notify_change = self.set_modified

    modified = Event
    def set_modified(self):
        print 'MN:set_modifeid'
        self.modified = True

    # cross section
    cs = DelegatesTo('calib')

    calibrated_ecb_law = Property(depends_on='modified')
    @cached_property
    def _get_calibrated_ecb_law(self):
        print 'NEW CALIBRATION'
        return self.calib.calibrated_ecb_law

    eps_cu = Property()
    def _get_eps_cu(self):
        return -self.cs.cc_law.eps_c_u

    eps_tu = Property()
    def _get_eps_tu(self):
        return self.calibrated_ecb_law.eps_tex_u

    n_eps = Int(5, auto_set=False, enter_set=True)
    eps_range = Property(depends_on='n_eps')
    @cached_property
    def _get_eps_range(self):
        eps_c_space = np.linspace(self.eps_cu, 0, self.n_eps)
        eps_t_space = np.linspace(0, self.eps_tu, self.n_eps)

        eps_ccu = 0.8 * self.eps_cu

        #eps_cc = self.eps_cu * np.ones_like(eps_c_space)
        eps_cc = np.linspace(eps_ccu, self.eps_cu, self.n_eps)
        eps_ct = self.eps_cu * np.ones_like(eps_t_space)
        eps_tc = self.eps_tu * np.ones_like(eps_c_space)
        eps_tt = self.eps_tu * np.ones_like(eps_t_space)

        eps1 = np.vstack([eps_c_space, eps_cc])
        eps2 = np.vstack([eps_t_space, eps_ct])
        eps3 = np.vstack([eps_tc, eps_c_space])
        eps4 = np.vstack([eps_tt, eps_t_space])

        return np.hstack([eps1, eps2, eps3, eps4])

    n_eps_range = Property(depends_on='n_eps')
    @cached_property
    def _get_n_eps_range(self):
        return self.eps_range.shape[1]
    #===========================================================================
    # MN Diagram
    #===========================================================================

    def _get_MN_fn(self, eps_lo, eps_up):
        self.cs.set(eps_lo=eps_lo,
                    eps_up=eps_up)
        return (self.cs.M, self.cs.N)

    MN_vct = Property(depends_on='modified')
    def _get_MN_vct(self):
        return np.vectorize(self._get_MN_fn)

    MN_arr = Property(depends_on='modified')
    @cached_property
    def _get_MN_arr(self):
        return self.MN_vct(self.eps_range[0, :], self.eps_range[1, :])

    #===========================================================================
    # f_eps Diagram
    #===========================================================================

    current_eps_idx = Int(0) # , auto_set = False, enter_set = True)

    def _current_eps_idx_changed(self):
        self._clear_fired()
        self._replot_fired()

    current_eps = Property(depends_on='current_eps_idx')
    @cached_property
    def _get_current_eps(self):
        return self.eps_range[(0, 1), self.current_eps_idx]

    current_MN = Property(depends_on='current_eps_idx')
    @cached_property
    def _get_current_MN(self):
        return self._get_MN_fn(*self.current_eps)

    #===========================================================================
    # Plotting
    #===========================================================================

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor='white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    data_changed = Event

    clear = Button
    def _clear_fired(self):
        self.figure.clear()
        self.data_changed = True

    replot = Button
    def _replot_fired(self):

        ax = self.figure.add_subplot(2, 2, 1)

        ax.plot(-self.eps_range, [0, 0.06], color='black')

        ax.plot(-self.current_eps, [0, 0.06], lw=3, color='red')

        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax = self.figure.add_subplot(2, 2, 2)

        ax.plot(self.MN_arr[0], -self.MN_arr[1], lw=2, color='blue')

        ax.plot(self.current_MN[0], -self.current_MN[1], 'g.', markersize=20.0, color='red')

        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.grid(b=None, which='major')

        self.cs.set(eps_lo=self.current_eps[0],
                    eps_up=self.current_eps[1])

        ax = self.figure.add_subplot(2, 2, 3)

        self.cs.plot_eps(ax)

        ax = self.figure.add_subplot(2, 2, 4)

        self.cs.plot_sig(ax)

        self.data_changed = True

    view = View(HSplit(Group(
                HGroup(
                Group(Item('n_eps', springy=True),
                      label='Discretization',
                      springy=True
                      ),
                springy=True,
                ),
                HGroup(
                Group(VGroup(
                      Item('cs', label='Cross section', show_label=False, springy=True,
                           editor=InstanceEditor(kind='live'),
                           ),
                      Item('calib', label='Calibration', show_label=False, springy=True,
                           editor=InstanceEditor(kind='live'),
                           ),
                      springy=True,
                      ),
                      label='Cross sectoin',
                      springy=True
                      ),
                springy=True,
                ),
                scrollable=True,
                             ),
                Group(HGroup(
                             Item('replot', show_label=False),
                             Item('clear', show_label=False),
                      ),
                      Item('current_eps_idx', editor=RangeEditor(low=0,
                                   high_name='n_eps_range',
                                   format='(%s)',
                                   mode='slider',
                                   auto_set=False,
                                   enter_set=False,
                                   ),
                           show_label=False,
                           ),
                      Item('figure', editor=MPLFigureEditor(),
                           resizable=True, show_label=False),
                      id='simexdb.plot_sheet',
                      label='plot sheet',
                      dock='tab',
                      ),
                       ),
                width=1.0,
                height=0.8,
                resizable=True,
                buttons=['OK', 'Cancel'])

if __name__ == '__main__':
    c = ECBLCalib(
                  Mu=3.49,
                  width=0.20,
                  n_rovings=23,
                  ecb_law_type='fbm',
                  cc_law_type='quadratic'             #eps_tu 0.0137279096658                              
                  )

    mn = ECBLMNDiagram(calib=c,
                       n_eps=30,
                      )

    mn.configure_traits()

