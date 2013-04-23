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
# Created on Dec 14, 2010 by: kelidas

from etsproxy.traits.api import \
    HasTraits, Int, Array, Str, implements, Range, Property, cached_property, \
     Float, Instance, Any, Interface, Event, on_trait_change, Button, Bool

from etsproxy.traits.ui.api import \
    View, Item, Group, VGroup, HGroup, HSplit, VSplit, Tabbed

from math import pi, e

from numpy import \
    sign, linspace, array, cos, sqrt, argmax, hstack, max, zeros_like, argwhere, loadtxt

from spirrid.i_rf import \
    IRF

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from spirrid.rf import RF

from ymb_pdistrib import YMB_RV

from matplotlib.figure import Figure

from ymb_data import YMBSource, YMBCutData, YMBSegmentData, YMBSlider, IYMBData

from ymb_pdistrib import YMBDistrib

from spirrid import SPIRRID, RV

from stats.pdistrib.pdistrib import PDistrib, IPDistrib

from yarn_symmetrical import DoublePulloutSym

class YarnPullOut(HasTraits):
    '''Idealization of the double sided pullout using the SPIRRID
    statistical integration tool.
    '''
    rf = Instance(DoublePulloutSym)
    def _rf_default(self):
        return DoublePulloutSym(tau_fr=2.6, l=0.0, d=25.5e-3, E_mod=72.0e3,
                                 theta=0.0, xi=0.0179, phi=1., L=30.0,
                                 free_fiber_end=True)
#        return DoublePulloutSym( tau_fr = 2.5, l = 0.01, d = 25.5e-3, E_mod = 70.0e3,
#                                 theta = 0.01, xi = 0.0179, phi = 1., n_f = 1723 )

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor='white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    pdf_theta_on = Bool(True)
    pdf_l_on = Bool(True)
    pdf_phi_on = Bool(True)
    pdf_xi_on = Bool(True)

    pdf_xi = Instance(IPDistrib)
    def _pdf_xi_default(self):
        pd = PDistrib(distr_choice='weibull_min', n_segments=30)
        pd.distr_type.set(shape=4.54, scale=0.017)
        return pd

    n_f = Float(1, auto_set=False, enter_set=True,
                 desc='Number of filaments in the yarn')

    pdf_theta = Instance(IPDistrib)

    pdf_l = Instance(IPDistrib)

    pdf_phi = Instance(IPDistrib)

    run = Button
    def _run_fired(self):
        self._redraw()

    clear = Button
    def _clear_fired(self):
        axes = self.figure.axes[0]
        axes.clear()
        self.data_changed = True

    w_max = Float(1.0, enter_set=True, auto_set=False)

    n_w_pts = Int(100, enter_set=True, auto_set=False)

    n_G_ipts = Int(30, enter_set=True, auto_set=False)

    e_arr = Property(Array, depends_on='w_max, n_w_pts')
    def _get_e_arr(self):
        return linspace(0.00, self.w_max, self.n_w_pts)

    lab = Str(' ', enter_set=True, auto_set=False)

    data_changed = Event(True)
    def _redraw(self):

        s = SPIRRID(q=self.rf,
                    sampling_type='LHS',
                    e_arr=self.e_arr,
                    n_int=self.n_G_ipts,
                    theta_vars=dict(tau_fr=2.6,
                                      l=0.0,
                                      d=25.5e-3,
                                      E_mod=72.0e3,
                                      theta=0.0,
                                      xi=0.0179,
                                      phi=1.,
                                      L=30.0
                                 ),
                    # codegen_type='weave'
                )
        # construct the random variables

        if self.pdf_xi_on:
            s.theta_vars['xi'] = RV('weibull_min', shape=4.54, scale=0.017)  # RV( pd = self.pdf_xi, name = 'xi', n_int = self.n_G_ipts )

        print self.pdf_theta.interp_ppf([0.01, 0.02])
        print YMB_RV('theta', distr=self.pdf_theta, n_int=self.n_G_ipts).distr.interp_ppf([0.01, 0.02])
        if self.pdf_theta_on:
            s.theta_vars['theta'] = YMB_RV('theta', distr=self.pdf_theta, n_int=self.n_G_ipts)

        if self.pdf_l_on:
            s.theta_vars['l'] = YMB_RV('l', distr=self.pdf_l, n_int=self.n_G_ipts)

        if self.pdf_phi_on:
            s.theta_vars['phi'] = YMB_RV('phi', distr=self.pdf_phi, n_int=self.n_G_ipts)

        # print 'checking unity', s.mu_q_arr()

        mu = s.mu_q_arr
        axes = self.figure.axes[0]
        # TODO:
        axes.plot(self.e_arr, mu * self.n_f,
                   linewidth=2, label=self.lab)

        axes.set_xlabel('crack opening w[mm]')
        axes.set_ylabel('force P[N]')
        axes.legend(loc='best')

        self.data_changed = True

    view = View(HSplit(
                    Group(Item('rf@', show_label=False),
                           label='Response function'
                           ),
                Tabbed(
                    Group(
                    Item('pdf_theta_on', show_label=False),
                    Item('pdf_theta@', show_label=False),
                    label='Slack',
                    ),
                    Group(
                    Item('pdf_l_on', show_label=False),
                    Item('pdf_l@', show_label=False),
                    label='Contact free length',
                    ),
                    Group(
                    Item('pdf_phi_on', show_label=False),
                    Item('pdf_phi@', show_label=False),
                    label='Contact fraction',
                    ),
                    Group(
                    Item('pdf_xi_on', show_label=False),
                    Item('pdf_xi@', show_label=False),
                    label='Strength',
                    ),
                    label='yarn data',
                    scrollable=True,
                    id='ymb.pullout.dist',
                    dock='tab',
                ),
                    Group(HGroup(Item('run', show_label=False, springy=True),
                                   Item('clear', show_label=False, springy=True)),
                           HGroup(Item('w_max', show_label=True, springy=True,
                                         tooltip='maximum crack-opening displacement'),
                                   Item('n_w_pts', show_label=True, springy=True,
                                         tooltip='number of points for crack-opening'),
                                   Item('n_G_ipts', show_label=True, springy=True,
                                         tooltip='number of integration points for the random variables'),
                                Item('lab', show_label=True, springy=True,
                                         tooltip='label of pull-out curve'),
                                 ),
                           Item('figure', style='custom',
                                  editor=MPLFigureEditor(),
                                  show_label=False),
                    label='Pull-out response',
                    id='ymb.pullout.figure',
                    dock='tab',
                ),
                id='ymb.pullout.split',
                dock='tab',
                ),
                id='ymb.pullout',
                resizable=True,
                scrollable=True,
                dock='tab',
                width=0.8,
                height=0.4
        )

class YMBPullOut(YarnPullOut):

    data = Instance(IYMBData)

    slider_max = Property()
    def _get_slider_max(self):
        return self.data.n_cuts - 1

    @on_trait_change('data.source_change')
    def _set_n_f(self):
        self.n_f = self.data.n_continuous_filaments_aver  # n_continuous_filaments

    pdf_theta = Property(Instance(IPDistrib), depends_on='data, data.input_change')
    @cached_property
    def _get_pdf_theta(self):
        middle_cut = int(self.slider_max / 2)
        slider = YMBSlider(var_enum='slack', data=self.data,
                            cut_slider_on=True,
                            cut_slider=middle_cut)
        return YMBDistrib(slider=slider)

    pdf_l = Property(Instance(IPDistrib), depends_on='data, data.input_change')
    @cached_property
    def _get_pdf_l(self):
        middle_cut = int(self.slider_max / 2)
        slider = YMBSlider(var_enum='bond free length', data=self.data,
                            cut_slider_on=True,
                            cut_slider=middle_cut)
        return YMBDistrib(slider=slider)

    pdf_phi = Property(Instance(IPDistrib), depends_on='data, data.input_change')
    @cached_property
    def _get_pdf_phi(self):
        slider = YMBSlider(var_enum='contact fraction', data=self.data,
                            cut_slider_on=False)
        return YMBDistrib(slider=slider)

if __name__ == '__main__':

    yarn_type = 'MAG'
    source = YMBSource(yarn_type=yarn_type)

    data = YMBSegmentData(source=source, cf_limit=1.)

    po = YMBPullOut(data=data)
    po.configure_traits()

