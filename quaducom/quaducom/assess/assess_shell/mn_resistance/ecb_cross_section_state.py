'''
Created on Sep 4, 2012

@todo: introduce the dock feature for the views
@todo: classify the state changes and provide examples.


@author: rch
'''
from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, \
    Trait, Event, on_trait_change, Instance, Button, Callable

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, VGroup, HGroup

from ecb_cross_section_geo import ECBCrossSectionGeo

import numpy as np

class ECBCrossSectionState(HasStrictTraits):

    cs_geo = Instance(ECBCrossSectionGeo)

    #---------------------------------------------------------------
    # Cross section characteristics needed for tensile specimens 
    #---------------------------------------------------------------

    # thickness of reinforced cross section
    #
    thickness = Property
    def _get_thickness(self):
        return self.cs_geo.thickness

    # total number of reinforcement layers [-]
    # 
    n_layers = Property
    def _get_n_layers(self):
        return self.cs_geo.n_layers

    #---------------------------------------------------------------
    # Cross section characteristics needed for bending specimens 
    #---------------------------------------------------------------

    # width of the cross section [m]
    #
    width = Property
    def _get_width(self):
        return self.cs_geo.width

    # number of rovings in 0-direction of one composite 
    # layer of the bending test [-]:
    #
    n_rovings = Property
    def _get_n_rovings(self):
        return self.cs_geo.n_rovings

    # cross section of one roving [mm**2]:
    #
    A_roving = Property
    def _get_A_roving(self):
        return self.cs_geo.A_roving

    #===========================================================================
    # material properties 
    #===========================================================================

    f_ck = Property
    def _get_f_ck(self):
        return self.cs_geo.f_ck

    eps_c_u = Property
    def _get_c_u(self):
        return self.cs_geo.c_u

    # ultimate textile stress measured in the tensile test [MPa]
    #
    sig_tex_u = Property
    def _get_sig_tex_u(self):
        return self.cs_geo.sig_tex_u

    notify_change = Callable

    #===========================================================================
    # State management
    #===========================================================================
    modified = Event
    @on_trait_change('+eps_input,cs_geo.modified')
    def set_modified(self):
        self.modified = True
        if self.notify_change:
            self.notify_change()

    #===========================================================================
    # Strain state
    #===========================================================================

    eps_up = Float(-0.0033, auto_set = False, enter_set = True, eps_input = True)
    eps_lo = Float(0.0140, auto_set = False, enter_set = True, eps_input = True)

    #===========================================================================
    # Discretization conform to the tex layers
    #===========================================================================

    eps_i_arr = Property(depends_on = 'modified')
    @cached_property
    def _get_eps_i_arr(self):
        '''CALIBRATION: derive the unknown constitutive law of the layer
        (effective crack bridge law)
        '''
        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_t'
        # ------------------------------------------------------------------------                
        thickness = self.thickness

        # strain at the height of each reinforcement layer [-]:
        #
        return self.eps_up + (self.eps_lo - self.eps_up) * self.z_ti_arr / thickness

    x = Property(depends_on = 'modified')
    @cached_property
    def _get_x(self):
        # heights of the compressive zone:
        #
        if self.eps_up == self.eps_lo:
            return (abs(self.eps_up) / (abs(self.eps_up - self.eps_lo * 1e-9)) *
                     self.thickness)
        else:
            return (abs(self.eps_up) / (abs(self.eps_up - self.eps_lo)) *
                     self.thickness)

    eps_ti_arr = Property(depends_on = 'modified')
    @cached_property
    def _get_eps_ti_arr(self):
        return (np.fabs(self.eps_i_arr) + self.eps_i_arr) / 2.0

    def convert_eps_tex_u_2_lo(self, eps_tex_u):
        eps_up = self.eps_up
        return eps_up + (eps_tex_u - eps_up) / self.z_ti_arr[0] * self.thickness

    def convert_eps_lo_2_tex_u(self, eps_lo):
        eps_up = self.eps_up
        return (eps_up + (eps_lo - eps_up) / self.thickness * self.z_ti_arr[0])

    eps_ci_arr = Property(depends_on = 'modified')
    @cached_property
    def _get_eps_ci_arr(self):
        return (-np.fabs(self.eps_i_arr) + self.eps_i_arr) / 2.0

    eps_cj_arr = Property(depends_on = 'modified')
    @cached_property
    def _get_eps_cj_arr(self):
        '''get compressive strain at each integration layer of the compressive zone [-]:
        for 'stress_case' flexion
        '''
        # for calibration us measured compressive strain
        # @todo: use mapped traits instead
        #
        eps_j_arr = (self.eps_up + (self.eps_lo - self.eps_up) * self.z_cj_arr /
                     self.thickness)
        return (-np.fabs(eps_j_arr) + eps_j_arr) / 2.0

    z_ti_arr = Property
    def _get_z_ti_arr(self):
        return self.cs_geo.z_ti_arr

    zz_ti_arr = Property
    def _get_zz_ti_arr(self):
        return self.cs_geo.zz_ti_arr

    n_cj = Property
    def _get_n_cj(self):
        return self.cs_geo.n_cj

    z_cj_arr = Property(depends_on = 'modified')
    @cached_property
    def _get_z_cj_arr(self):
        '''Get the discretizaton of the  compressive zone
        '''
        if self.eps_up <= 0:
            zx = min(self.thickness, self.x)
            return np.linspace(0, zx, self.n_cj)
        elif self.eps_lo <= 0:
            return np.linspace(self.x, self.thickness, self.n_cj)
        else:
            return np.array([0], dtype = 'f')


    # distance of reinforcement layers from the bottom 
    #
    zz_cj_arr = Property(depends_on = 'modified')
    @cached_property
    def _get_zz_cj_arr(self):
        return self.thickness - self.z_cj_arr

    #===========================================================================
    # Compressive concrete constitutive law
    #===========================================================================
    cc_law = Property
    def _get_cc_law(self):
        '''Construct the compressive concrete law'''
        return self.cs_geo.cc_law

    sig_cj_arr = Property(depends_on = 'modified')
    @cached_property
    def _get_sig_cj_arr(self):
        return -self.cc_law.mfn_vct(-self.eps_cj_arr)

    f_cj_arr = Property(depends_on = 'modified')
    @cached_property
    def _get_f_cj_arr(self):
        return self.width * self.sig_cj_arr
    #===========================================================================
    # Effective crack bridge law
    #===========================================================================
    ecb_law = Property
    def _get_ecb_law(self):
        return self.cs_geo.ecb_law

    sig_ti_arr = Property(depends_on = 'modified')
    @cached_property
    def _get_sig_ti_arr(self):
        '''force at the height of each reinforcement layer [kN]:
        '''
        return self.ecb_law.mfn_vct(self.eps_ti_arr)

    #===========================================================================
    # Layer conform discretization of the tensile zone
    #===========================================================================
    f_ti_arr = Property(depends_on = 'modified')
    @cached_property
    def _get_f_ti_arr(self):
        '''force at the height of each reinforcement layer [kN]:
        '''
        sig_ti_arr = self.sig_ti_arr
        n_rovings = self.n_rovings
        A_roving = self.A_roving
        return sig_ti_arr * n_rovings * A_roving / 1000.

    #===========================================================================
    # Stress resultants
    #===========================================================================

    N = Property(depends_on = 'modified')
    @cached_property
    def _get_N(self):
        '''Get the resulting normal force.
        '''
        N_tk = sum(self.f_ti_arr)
        N_ck = np.trapz(self.sig_cj_arr * self.width, self.z_cj_arr) * 1000.0

        N_internal = N_ck + N_tk
        return N_internal

    M = Property(depends_on = 'modified')
    @cached_property
    def _get_M(self):
        M_tk = np.dot(self.f_ti_arr, self.z_ti_arr)
        M_ck = np.trapz(self.sig_cj_arr * self.width * self.z_cj_arr, self.z_cj_arr) * 1000.0

        M_internal_ = M_tk + M_ck

        # moment evaluated with respect to the center line
        #
        M_internal = M_internal_ - self.N * self.thickness / 2.

        # return the internal stress resultants

        return M_internal

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor = 'white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    data_changed = Event

    replot = Button
    def _replot_fired(self):
        self.plot(self.figure)
        self.data_changed = True

    def plot(self, figure, eps_range = None, f_range = None):

        d = self.thickness - self.x

        figure.clear()
        ax = figure.add_subplot(2, 2, 1)
        self.plot_eps(ax, eps_range)

        ax = figure.add_subplot(2, 2, 2)
        self.plot_sig(ax, f_range)

        ax = figure.add_subplot(2, 2, 3)
        self.cc_law.plot(ax)

        ax = figure.add_subplot(2, 2, 4)
        self.ecb_law.plot(ax)

    def plot_eps(self, ax, eps_range = None):
        #ax = self.figure.gca()

        d = self.thickness
        # eps ti
        ax.plot([-self.eps_lo, -self.eps_up], [0, self.thickness], color = 'black')
        ax.hlines(self.zz_ti_arr, [0], -self.eps_ti_arr, lw = 4, color = 'red')

        # eps cj
        ec = np.hstack([self.eps_cj_arr] + [0, 0])
        zz = np.hstack([self.zz_cj_arr] + [0, self.thickness ])
        ax.fill(-ec, zz, color = 'blue')

        # reinforcement layers
        if eps_range == None:
            eps_range = np.array([-max(0.0, self.eps_lo),
                                  - min(0.0, self.eps_up)], dtype = 'float')
        else:
            ax.set_xlim(eps_range)

        eps_range = np.asarray(eps_range)

        z_ti_arr = np.ones_like(eps_range)[:, None] * self.z_ti_arr[None, :]
        ax.plot(eps_range, z_ti_arr, 'k--', color = 'black')

        # neutral axis
        x = self.x
        ax.plot(eps_range, [d - x, d - x], 'k--', color = 'green', lw = 2)

        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

    def plot_sig(self, ax, f_range = None):

        d = self.thickness
        # f ti
        ax.hlines(self.zz_ti_arr, [0], -self.f_ti_arr, lw = 4, color = 'red')

        # f cj
        f_c = np.hstack([self.f_cj_arr] + [0, 0])
        zz = np.hstack([self.zz_cj_arr] + [0, self.thickness ])

        ax.fill(-f_c, zz, color = 'blue')

        # reinforcement layers
        if f_range == None:
            f_range = np.array([-np.max(self.f_ti_arr), -np.min(f_c)], dtype = 'float_')
        else:
            f_orig_range = np.array([-np.max(self.f_ti_arr), -np.min(f_c)], dtype = 'float_')
            print 'orig f_range', f_orig_range
            f_range = np.asarray(f_range)
            print 'supplied f_range', f_range

        ax.set_xlim(f_range)

        z_ti_arr = np.ones_like(f_range)[:, None] * self.z_ti_arr[None, :]
        ax.plot(f_range, z_ti_arr, 'k--', color = 'black')

        # neutral axis
        x = self.x
        ax.plot(f_range, [d - x, d - x], 'k--', color = 'green', lw = 2)

        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

    view = View(HSplit(Group(
                HGroup(
                Group(Item('cs_geo', springy = True),
                      ),
                Group(Item('eps_up', label = 'Upper strain', springy = True),
                      Item('eps_lo', label = 'Lower strain'),
                      label = 'Strain',
                      springy = True
                      ),
                springy = True,
                ),
                Group(
                HGroup(Item('M', springy = True, style = 'readonly'),
                       Item('N', springy = True, style = 'readonly'),
                       ),
                       label = 'Stress resultants'
                       ),
                scrollable = True,
                             ),
                Group(Item('replot', show_label = False),
                      Item('figure', editor = MPLFigureEditor(),
                           resizable = True, show_label = False),
                      id = 'simexdb.plot_sheet',
                      label = 'plot sheet',
                      dock = 'tab',
                      ),
                       ),
                width = 0.8,
                height = 0.7,
                resizable = True,
                buttons = ['OK', 'Cancel'])

if __name__ == '__main__':
    cs_geo = ECBCrossSectionGeo(# 7d: f_ck,cube = 62 Mecs.    csPa; f_ck,cyl = 62/1.2=52
                                 # 9d: f_ck,cube = 66.8 MPa; f_ck,cyl = 55,7
                                   cc_law_params = dict(bilinear = dict(eps_c_u = -0.0033,
                                                                        f_ck = 55.7,
                                                                        E_c = 29e+3),
                                                        quadratic = dict(eps_c_u = -0.0033,
                                                                        f_ck = 55.7,
                                                                        )
                                                        ),
                                 ecb_law_type = 'linear',
                                 cc_law_type = 'quadratic'
                                 )

    ecs = ECBCrossSectionState(cs_geo = cs_geo)

    print 'initial'
    ecs.eps_up = -0.0033
    ecs.eps_lo = 0.014

    print 'z_cj'
    print ecs.z_cj_arr
    print ecs.eps_cj_arr
    print 'x', ecs.x

    print 'z_ti'
    print ecs.z_ti_arr
    print ecs.eps_ti_arr

    print 'changed n_rovings'
    ecs.cs_geo.n_layers = 6

    print 'z_ti'
    print ecs.z_ti_arr
    print ecs.eps_ti_arr

    eps_tex_u = ecs.convert_eps_lo_2_tex_u(0.014)
    print 'eps_tex_u', eps_tex_u
    eps_lo = ecs.convert_eps_tex_u_2_lo(eps_tex_u)
    print 'eps_lo', eps_lo

    print 'M', ecs.M
    print 'N', ecs.N

    ecs.configure_traits(kind = 'live')

