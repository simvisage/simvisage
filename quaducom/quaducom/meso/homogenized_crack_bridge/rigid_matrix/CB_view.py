'''
Created on 03.07.2013

@author: rostar
'''
import numpy as np
from scipy.interpolate import interp1d
import os
from etsproxy.traits.api import HasTraits, Property, Array, \
     cached_property, Float, Int, Instance, Event, List
from etsproxy.traits.ui.api import Item, View, Group, HSplit, VGroup, Tabbed
from etsproxy.traits.ui.menu import OKButton, CancelButton
from matplotlib.figure import Figure
from quaducom.micro.resp_func.CB_clamped_rand_xi import CBClampedRandXi
from spirrid.spirrid import SPIRRID
from spirrid.rv import RV
from util.traits.editors.mpl_figure_editor import MPLFigureEditor
from etsproxy.traits.ui.api import ModelView
from etsproxy.util.home_directory import get_home_directory
from scipy.special import gamma as gamma_func
from scipy.optimize import newton
from matplotlib import pyplot as plt


class Model(HasTraits):

    test_xdata = Array
    test_ydata = Array
    sV0 = Float(auto_set=False, enter_set=True, params=True)
    m = Float(auto_set=False, enter_set=True, params=True)
    w_min = Float(auto_set=False, enter_set=True, params=True)
    w_max = Float(auto_set=False, enter_set=True, params=True)
    w_pts = Int(auto_set=False, enter_set=True, params=True)
    n_int = Int(auto_set=False, enter_set=True, params=True)
    w2_min = Float(auto_set=False, enter_set=True, params=True)
    w2_max = Float(auto_set=False, enter_set=True, params=True)
    w2_pts = Int(auto_set=False, enter_set=True, params=True)
    #tau_scale = Float(auto_set=False, enter_set=True, params=True)
    tau_loc = Float(auto_set=False, enter_set=True, params=True)
    #tau_shape = Float(auto_set=False, enter_set=True, params=True)
    Ef = Float(auto_set=False, enter_set=True, params=True)
    lm = Float(auto_set=False, enter_set=True, params=True)
    V_f = Float(.01, params=True)
    r = Float(3.5e-3, params=True)
    CS = Float(8.0, auto_set=False, enter_set=True, params=True)

    tau_scale = Property(Float, depends_on='test_xdata,test_ydata,Ef,r,V_f,CS')
    @cached_property
    def _get_tau_scale(self):
        w_hat = 0.075
        sigmaf_hat = self.interpolate_experiment(w_hat)
        mu_sqrt_tau = sigmaf_hat / np.sqrt(2. * self.Ef * w_hat / self.r)
        mu_tau = 1.3 * self.r *3.0* (1.-self.V_f) / (2. * self.V_f * self.CS)
        def scale_res(scale):
            res = scale - ((mu_sqrt_tau) / (gamma_func(mu_tau/scale + 0.5)/gamma_func(mu_tau/scale)))**2
            return res
        scale = newton(scale_res, 0.52)
        return scale
    
    tau_shape = Property(Float, depends_on='test_xdata,test_ydata,Ef,r,V_f,CS')
    @cached_property
    def _get_tau_shape(self):
        mu_tau = 1.3 * self.r * 3.0 * (1.-self.V_f) / (2. * self.V_f * self.CS)
        return mu_tau/self.tau_scale

    w = Property(Array)
    def _get_w(self):
        return np.linspace(self.w_min, self.w_max, self.w_pts)

    w2 = Property(Array)
    def _get_w2(self):
        return np.linspace(self.w2_min, self.w2_max, self.w2_pts)

    interpolate_experiment = Property(depends_on='test_xdata, test_ydata')
    @cached_property
    def _get_interpolate_experiment(self):
        return interp1d(self.test_xdata, self.test_ydata,
                        bounds_error=False, fill_value=0.0)

    model_rand = Property(Array)
    def _get_model_rand(self):
        cb = CBClampedRandXi(pullout=False)
        spirrid = SPIRRID(q=cb, sampling_type='LHS')
        sV0 = self.sV0
        tau_scale = self.tau_scale
        V_f = 1.0
        r = 3.5e-3
        m = self.m
        tau = RV('gamma', shape=self.tau_shape, scale=tau_scale, loc=self.tau_loc)
        n_int = self.n_int
        w = self.w
        lm = 1000.
        spirrid.eps_vars = dict(w=w)
        spirrid.theta_vars = dict(tau=tau, E_f=self.Ef, V_f=V_f, r=r, m=m, sV0=sV0, lm=lm)
        spirrid.n_int = n_int
        sigma_c = spirrid.mu_q_arr / self.r ** 2
        return sigma_c

    model_extrapolate = Property(Array)
    def _get_model_extrapolate(self):
        cb = CBClampedRandXi(pullout=False)
        spirrid = SPIRRID(q=cb, sampling_type='LHS')
        sV0 = self.sV0
        V_f = 1.0
        r = 3.5e-3
        m = self.m
        tau = RV('gamma', shape=self.tau_shape, scale=self.tau_scale, loc=self.tau_loc)
        n_int = self.n_int
        w = self.w2
        lm = self.lm
        spirrid.eps_vars = dict(w=w)
        spirrid.theta_vars = dict(tau=tau, E_f=self.Ef, V_f=V_f, r=r, m=m, sV0=sV0, lm=lm)
        spirrid.n_int = n_int
        sigma_c = spirrid.mu_q_arr / self.r ** 2
        return sigma_c

class CBView(ModelView):

    def __init__(self, **kw):
        super(CBView, self).__init__(**kw)
        self.on_trait_change(self.refresh, 'model.+params')
        self.refresh()

    model = Instance(Model)

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor='white')
        return figure

    figure2 = Instance(Figure)
    def _figure2_default(self):
        figure = Figure(facecolor='white')
        return figure

    data_changed = Event

    def plot(self, fig, fig2):
        figure = fig
        figure.clear()
        axes = figure.gca()
        # plot PDF
        axes.plot(self.model.w, self.model.model_rand, lw=2.0, color='blue', \
                  label='model')
        axes.plot(self.model.w, self.model.interpolate_experiment(self.model.w), lw=1.0, color='black', \
                  label='experiment')
        axes.legend()

#         figure2 = fig2
#         figure2.clear()
#         axes = figure2.gca()
#         # plot PDF
#         axes.plot(self.model.w2, self.model.model_extrapolate, lw=2.0, color='red', \
#                   label='model')
#         axes.legend()

    def refresh(self):
        self.plot(self.figure, self.figure2)
        self.data_changed = True

    traits_view = View(HSplit(VGroup(Group(Item('model.tau_scale'),
                                           Item('model.tau_shape'),
                                           Item('model.tau_loc'),
                                           Item('model.m'),
                                           Item('model.sV0'),
                                           Item('model.Ef'),
                                           Item('model.w_min'),
                                           Item('model.w_max'),
                                           Item('model.w_pts'),
                                           Item('model.n_int'),
                                           Item('model.w2_min'),
                                           Item('model.w2_max'),
                                           Item('model.w2_pts'),
                                           Item('model.lm'),
                                           Item('model.CS'),
                                           ),
                                      id='pdistrib.distr_type.pltctrls',
                                      label='Distribution parameters',
                                      scrollable=True,
                                      ),
                                Tabbed(Group(Item('figure',
                                            editor=MPLFigureEditor(),
                                            show_label=False,
                                            resizable=True),
                                            scrollable=True,
                                            label='Plot',
                                            ),
                                        label='Plot',
                                        id='pdistrib.figure.params',
                                        dock='tab',
                                       ),
                              Tabbed(Group(Item('figure2',
                                            editor=MPLFigureEditor(),
                                            show_label=False,
                                            resizable=True),
                                            scrollable=True,
                                            label='Plot',
                                            ),
                                        label='Plot',
                                        id='pdistrib.figure2',
                                        dock='tab',
                                       ),
                                dock='tab',
                                id='pdistrib.figure.view'
                                ),
                                id='pdistrib.view',
                                dock='tab',
                                title='Statistical distribution',
                                buttons=[OKButton, CancelButton],
                                scrollable=True,
                                resizable=True,
                                width=600, height=400
                        )

if __name__ == '__main__':

    model = Model(w_min=0.0, w_max=3.0, w_pts=200,
                  w2_min=0.0, w2_max=.5, w2_pts=200,
                  sV0=7.7e-3, m=9.0, tau_loc=0.0, Ef=200e3,
                  lm=20., n_int=100, CS=13.0)

    home_dir = get_home_directory()
    path = [home_dir, 'git',  # the path of the data file
            'rostar',
            'scratch',
            'diss_figs',
            'CB1.txt']
    filepath = os.path.join(*path)
    file1 = open(filepath, 'r')
    cb = np.loadtxt(file1, delimiter=';')
    model.test_xdata = -cb[:, 2] / 4. - cb[:, 3] / 4. - cb[:, 4] / 2.
    model.test_ydata = cb[:, 1] / (11. * 0.445) * 1000

    cb = CBView(model=model)
    cb.refresh()
    cb.configure_traits()
