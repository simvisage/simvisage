'''
Created on Sep 4, 2012

@todo: introduce the dock feature for the views
@todo: classify the state changes and provide examples.


@author: rch
'''
from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, \
    Trait, Event, on_trait_change, Instance, Button, Callable, \
    DelegatesTo, Constant

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, VGroup, HGroup

from ecb_law import \
    ECBLBase, ECBLLinear, ECBLFBM, ECBLCubic, ECBLBilinear

from constitutive_law import \
    ConstitutiveLawModelView

from ecb_reinf_component import \
    ECBReinfComponent, \
    ECB_COMPONENT_CHANGE, \
    ECB_COMPONENT_AND_EPS_CHANGE

import numpy as np

class ECBReinfTexUniform(ECBReinfComponent):
    '''Cross section characteristics needed for tensile specimens
    '''

    height = DelegatesTo('matrix_cs')
    '''height of reinforced cross section
    '''

    n_layers = Int(12, auto_set=False, enter_set=True, geo_input=True)
    '''total number of reinforcement layers [-]
    '''

    n_rovings = Int(23, auto_set=False, enter_set=True, geo_input=True)
    '''number of rovings in 0-direction of one composite layer of the
    bending test [-]:
    '''

    A_roving = Float(0.461, auto_set=False, enter_set=True, geo_input=True)
    '''cross section of one roving [mm**2]'''

    def convert_eps_tex_u_2_lo(self, eps_tex_u):
        '''Convert the strain in the lowest reinforcement layer at failure
        to the strain at the bottom of the cross section'''
        eps_up = self.state.eps_up
        return eps_up + (eps_tex_u - eps_up) / self.z_ti_arr[0] * self.height

    def convert_eps_lo_2_tex_u(self, eps_lo):
        '''Convert the strain at the bottom of the cross section to the strain
        in the lowest reinforcement layer at failure'''
        eps_up = self.state.eps_up
        return (eps_up + (eps_lo - eps_up) / self.height * self.z_ti_arr[0])

    '''Convert the MN to kN
    '''

    #===========================================================================
    # material properties 
    #===========================================================================

    sig_tex_u = Float(1216., auto_set=False, enter_set=True,
                      tt_input=True)
    '''Ultimate textile stress measured in the tensile test [MPa]
    '''

    #===========================================================================
    # Distribution of reinforcement
    #===========================================================================

    s_tex_z = Property(depends_on=ECB_COMPONENT_AND_EPS_CHANGE)
    '''spacing between the layers [m]'''
    @cached_property
    def _get_s_tex_z(self):
        return self.height / (self.n_layers + 1)


    z_ti_arr = Property(depends_on=ECB_COMPONENT_AND_EPS_CHANGE)
    '''property: distance of each reinforcement layer from the top [m]:
    '''
    @cached_property
    def _get_z_ti_arr(self):
        return np.array([ self.height - (i + 1) * self.s_tex_z
                         for i in range(self.n_layers) ],
                      dtype=float)

    zz_ti_arr = Property
    '''property: distance of reinforcement layers from the bottom
    '''
    def _get_zz_ti_arr(self):
        return self.height - self.z_ti_arr

    #===========================================================================
    # Discretization conform to the tex layers
    #===========================================================================

    eps_i_arr = Property(depends_on=ECB_COMPONENT_AND_EPS_CHANGE)
    '''Strain at the level of the i-th reinforcement layer
    '''
    @cached_property
    def _get_eps_i_arr(self):
        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_t'
        # ------------------------------------------------------------------------                
        height = self.height
        eps_lo = self.state.eps_lo
        eps_up = self.state.eps_up
        # strain at the height of each reinforcement layer [-]:
        #
        return eps_up + (eps_lo - eps_up) * self.z_ti_arr / height

    eps_ti_arr = Property(depends_on=ECB_COMPONENT_AND_EPS_CHANGE)
    '''Tension strain at the level of the i-th layer of the fabrics
    '''
    @cached_property
    def _get_eps_ti_arr(self):
        return (np.fabs(self.eps_i_arr) + self.eps_i_arr) / 2.0

    eps_ci_arr = Property(depends_on=ECB_COMPONENT_AND_EPS_CHANGE)
    '''Compression strain at the level of the i-th layer.
    '''
    @cached_property
    def _get_eps_ci_arr(self):
        return (-np.fabs(self.eps_i_arr) + self.eps_i_arr) / 2.0

    #===========================================================================
    # Effective crack bridge law
    #===========================================================================
    ecb_law_type = Trait('fbm', dict(fbm=ECBLFBM,
                                  cubic=ECBLCubic,
                                  linear=ECBLLinear,
                                  bilinear=ECBLBilinear),
                      tt_input=True)
    '''Selector of the effective crack bridge law type
    ['fbm', 'cubic', 'linear', 'bilinear']'''

    ecb_law = Property(Instance(ECBLBase), depends_on='+tt_input')
    '''Effective crack bridge law corresponding to ecb_law_type'''
    @cached_property
    def _get_ecb_law(self):
        return self.ecb_law_type_(sig_tex_u=self.sig_tex_u, cs=self)

    show_ecb_law = Button
    '''Button launching a separate view of the effective crack bridge law.
    '''
    def _show_ecb_law_fired(self):
        ecb_law_mw = ConstitutiveLawModelView(model=self.ecb_law)
        ecb_law_mw.edit_traits(kind='live')
        return

    tt_modified = Event

    sig_ti_arr = Property(depends_on=ECB_COMPONENT_AND_EPS_CHANGE)
    '''Stresses at the i-th fabric layer.
    '''
    @cached_property
    def _get_sig_ti_arr(self):
        return self.ecb_law.mfn_vct(self.eps_ti_arr)

    f_ti_arr = Property(depends_on=ECB_COMPONENT_AND_EPS_CHANGE)
    '''force at the height of each reinforcement layer [kN]:
    '''
    @cached_property
    def _get_f_ti_arr(self):
        sig_ti_arr = self.sig_ti_arr
        n_rovings = self.n_rovings
        A_roving = self.A_roving
        return sig_ti_arr * n_rovings * A_roving / self.unit_conversion_factor

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

        d = self.height
        # eps ti
        ax.plot([-self.eps_lo, -self.eps_up], [0, self.height], color='black')
        ax.hlines(self.zz_ti_arr, [0], -self.eps_ti_arr, lw=4, color='red')

        # eps cj
        ec = np.hstack([self.eps_cj_arr] + [0, 0])
        zz = np.hstack([self.zz_cj_arr] + [0, self.height ])
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

        d = self.height
        # f ti
        ax.hlines(self.zz_ti_arr, [0], -self.f_ti_arr, lw=4, color='red')

        # f cj
        f_c = np.hstack([self.f_cj_arr] + [0, 0])
        zz = np.hstack([self.zz_cj_arr] + [0, self.height ])
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
                Group(Item('height', springy=True),
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
    ecs = ECBReinfTexUniform(# 7d: f_ck,cube = 62 MPa; f_ck,cyl = 62/1.2=52
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
