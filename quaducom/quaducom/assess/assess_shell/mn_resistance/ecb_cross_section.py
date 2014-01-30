'''
Created on Sep 4, 2012

@todo: introduce the dock feature for the views
@todo: classify the state changes and provide examples.


@author: rch
'''
from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, \
    Trait, Event, on_trait_change, Instance, Button, Callable, \
    DelegatesTo, Constant, List

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, VGroup, HGroup

from ecb_cross_section_state import \
    ECBCrossSectionState

from ecb_cross_section_geo import \
    ECBCrossSectionComponent

import numpy as np

class ECBCrossSection(ECBCrossSectionState):
    '''Cross section characteristics needed for tensile specimens
    '''
    components = List(ECBCrossSectionComponent)
    '''Components of the cross section including the matrix and reinforcement.
    '''

    components_with_state = Property(depends_on='components')
    '''Components linked to the strain state of the cross section
    '''
    @cached_property
    def _get_components_with_state(self):
        for c in self.components:
            c.state = self
        return self.components

    unit_conversion_factor = Constant(1000.0)

    '''Convert the MN to kN
    '''

    #===========================================================================
    # State management
    #===========================================================================
    changed = Event
    '''Notifier of a changed in some component of a cross section
    '''

    @on_trait_change('+eps_input')
    def _notify_eps_change(self):
        self.changed = True
        for c in self.components:
            c.eps_changed = True

    #===========================================================================
    # Cross-sectional stress resultants
    #===========================================================================

    N = Property(depends_on='changed')
    '''Get the resulting normal force.
    '''
    @cached_property
    def _get_N(self):
        return np.sum([c.N for c in self.components_with_state])

    M = Property(depends_on='changed')
    '''Get the resulting moment.
    '''
    @cached_property
    def _get_M(self):
        M = np.sum([c.M for c in self.components_with_state])
        return M - self.N * self.height / 2.

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor='white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    data_changed = Event

    replot = Button
    def _replot_fired(self):

        self.figure.clear()
        ax = self.figure.add_subplot(2, 2, 1)
        self.plot_eps(ax)

        ax = self.figure.add_subplot(2, 2, 2)
        self.plot_sig(ax)

        ax = self.figure.add_subplot(2, 2, 3)
        self.cc_law.plot(ax)

        ax = self.figure.add_subplot(2, 2, 4)
        self.ecb_law.plot(ax)

        self.data_changed = True

    def plot_eps(self, ax):
        #ax = self.figure.gca()

        d = self.thickness
        # eps ti
        ax.plot([-self.eps_lo, -self.eps_up], [0, self.thickness], color='black')
        ax.hlines(self.zz_ti_arr, [0], -self.eps_ti_arr, lw=4, color='red')

        # eps cj
        ec = np.hstack([self.eps_cj_arr] + [0, 0])
        zz = np.hstack([self.zz_cj_arr] + [0, self.thickness ])
        ax.fill(-ec, zz, color='blue')

        # reinforcement layers
        eps_range = np.array([max(0.0, self.eps_lo),
                              min(0.0, self.eps_up)], dtype='float')
        z_ti_arr = np.ones_like(eps_range)[:, None] * self.z_ti_arr[None, :]
        ax.plot(-eps_range, z_ti_arr, 'k--', color='black')

        # neutral axis
        ax.plot(-eps_range, [d, d], 'k--', color='green', lw=2)

        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

    def plot_sig(self, ax):

        d = self.thickness
        # f ti
        ax.hlines(self.zz_ti_arr, [0], -self.f_ti_arr, lw=4, color='red')

        # f cj
        f_c = np.hstack([self.f_cj_arr] + [0, 0])
        zz = np.hstack([self.zz_cj_arr] + [0, self.thickness ])
        ax.fill(-f_c, zz, color='blue')

        f_range = np.array([np.max(self.f_ti_arr), np.min(f_c)], dtype='float_')
        # neutral axis
        ax.plot(-f_range, [d, d], 'k--', color='green', lw=2)

        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

    view = View(HSplit(Group(
                HGroup(
                Group(Item('thickness', springy=True),
                      Item('width'),
                      Item('n_layers'),
                      Item('n_rovings'),
                      Item('A_roving'),
                      label='Geometry',
                      springy=True
                      ),
                Group(Item('eps_up', label='Upper strain', springy=True),
                      Item('eps_lo', label='Lower strain'),
                      label='Strain',
                      springy=True
                      ),
                springy=True,
                ),
                HGroup(
                Group(VGroup(
                      Item('cc_law_type', show_label=False, springy=True),
                      Item('cc_law', label='Edit', show_label=False, springy=True),
                      Item('show_cc_law', label='Show', show_label=False, springy=True),
                      springy=True
                      ),
                      Item('f_ck', label='Compressive strength'),
                      Item('n_cj', label='Discretization'),
                      label='Concrete',
                      springy=True
                      ),
                Group(VGroup(
                      Item('ecb_law_type', show_label=False, springy=True),
                      Item('ecb_law', label='Edit', show_label=False, springy=True),
                      Item('show_ecb_law', label='Show', show_label=False, springy=True),
                      springy=True,
                      ),
                      label='Reinforcement',
                      springy=True
                      ),
                springy=True,
                ),
                Group(Item('s_tex_z', label='vertical spacing', style='readonly'),
                      label='Layout',
                      ),
                Group(
                HGroup(Item('M', springy=True, style='readonly'),
                       Item('N', springy=True, style='readonly'),
                       ),
                       label='Stress resultants'
                       ),
                scrollable=True,
                             ),
                Group(Item('replot', show_label=False),
                      Item('figure', editor=MPLFigureEditor(),
                           resizable=True, show_label=False),
                      id='simexdb.plot_sheet',
                      label='plot sheet',
                      dock='tab',
                      ),
                       ),
                width=0.8,
                height=0.7,
                resizable=True,
                buttons=['OK', 'Cancel'])

if __name__ == '__main__':
    ecs = ECBCrossSection(# 7d: f_ck,cube = 62 MPa; f_ck,cyl = 62/1.2=52
                           # 9d: f_ck,cube = 66.8 MPa; f_ck,cyl = 55,7
                           f_ck=55.7,

                           ecb_law_type='fbm',
                           cc_law_type='quadratic'
                           )

    print 'initial'
    ecs.eps_up = -0.0033
    ecs.eps_lo = 0.014

    print 'z_cj'
    print ecs.z_cj_arr
    print ecs.eps_cj_arr
    print 'x', ecs.x

    eps_tex_u = ecs.convert_eps_lo_2_tex_u(0.014)
    print 'eps_tex_u', eps_tex_u
    eps_lo = ecs.convert_eps_tex_u_2_lo(eps_tex_u)
    print 'eps_lo', eps_lo

    print 'M', ecs.M
    print 'N', ecs.N

    ecs.configure_traits()
