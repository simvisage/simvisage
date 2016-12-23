'''
Created on 12.01.2016

@author: Yingxiong
'''
from math import *

from matplotlib.figure import Figure
from scipy.interpolate import interp1d
from traits.api import \
    HasTraits, Property, Instance, cached_property, Str, Button, Enum, \
    Range, on_trait_change, Array, List, Float
from traitsui.api import \
    View, Item, Group, VGroup, HSplit, TreeEditor, TreeNode

from fets1d52ulrhfatigue import FETS1D52ULRHFatigue
from ibvpy.api import BCDof, TStepper, TLoop
from mats_bondslip import MATSEvalFatigue
import numpy as np
from util.traits.editors import MPLFigureEditor
from view.ui import BMCSTreeNode, BMCSLeafNode


class Material(BMCSLeafNode):

    node_name = Str('material parameters')
    E_b = Float(2000,
                label="E_b ",
                desc="Bond Stiffness",
                enter_set=True,
                auto_set=False)

    gamma = Float(0,
                  label="Gamma ",
                  desc="Kinematic hardening modulus",
                  enter_set=True,
                  auto_set=False)

    K = Float(0,
              label="K ",
              desc="Isotropic harening",
              enter_set=True,
              auto_set=False)

    S = Float(1.0,
              label="S ",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    r = Float(1,
              label="r ",
              desc="Damage cumulation parameter",
              enter_set=True,
              auto_set=False)

    tau_pi_bar = Float(5,
                       label="Tau_pi_bar ",
                       desc="Reversibility limit",
                       enter_set=True,
                       auto_set=False)

    view = View(VGroup(Group(Item('E_b'),
                             Item('tau_pi_bar'), show_border=True, label='Bond Stiffness and reversibility limit'),
                       Group(Item('gamma'),
                             Item('K'), show_border=True, label='Hardening parameters'),
                       Group(Item('S'),
                             Item('r'), show_border=True, label='Damage cumulation parameters')))


class Geometry(BMCSLeafNode):
    L_x = Range(1, 100, value=10)
    A_m = Float(100 * 8 - 9 * 1.85, desc='matrix area [mm2]')
    A_f = Float(9 * 1.85, desc='reinforcement area [mm2]')
    P_b = Float(9 * np.sqrt(np.pi * 4 * 1.85),
                desc='perimeter of the bond interface [mm]')


class LoadingScenario(BMCSLeafNode):

    name = Str('<unknown>')
    number_of_cycles = Float(1.0)
    maximum_slip = Float(0.2)
    number_of_increments = Float(10)
    loading_type = Enum("Monotonic", "Cyclic")
    amplitude_type = Enum("Increased_Amplitude", "Constant_Amplitude")
    loading_range = Enum("Non_symmetric", "Symmetric")

    d_t = Float(0.005)
    t_max = Float(1.)
    k_max = Float(100)
    tolerance = Float(1e-4)

    d_array = Property(
        depends_on=' maximum_slip , number_of_cycles , loading_type , loading_range , amplitude_type ')

    @cached_property
    def _get_d_array(self):

        if self.loading_type == "Monotonic":
            self.number_of_cycles = 1
        d_levels = np.linspace(0, self.maximum_slip, self.number_of_cycles * 2)
        d_levels[0] = 0

        if self.amplitude_type == "Increased_Amplitude" and self.loading_range == "Symmetric":
            d_levels.reshape(-1, 2)[:, 0] *= -1
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], self.number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.amplitude_type == "Increased_Amplitude" and self.loading_range == "Non_symmetric":
            d_levels.reshape(-1, 2)[:, 0] *= 0
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], self.number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.amplitude_type == "Constant_Amplitude" and self.loading_range == "Symmetric":
            d_levels.reshape(-1, 2)[:, 0] = -self.maximum_slip
            d_levels[0] = 0
            d_levels.reshape(-1, 2)[:, 1] = self.maximum_slip
            d_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(d_history[i], d_history[i + 1], self.number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

        if self.amplitude_type == "Constant_Amplitude" and self.loading_range == "Non_symmetric":
            d_levels.reshape(-1, 2)[:, 0] *= 0
            d_levels.reshape(-1, 2)[:, 1] = self.maximum_slip
            s_history = d_levels.flatten()
            d_arr = np.hstack([np.linspace(s_history[i], s_history[i + 1], self.number_of_increments)
                               for i in range(len(d_levels) - 1)])

            return d_arr

    time_func = Property(depends_on='maximum_slip, t_max , d_array ')

    @cached_property
    def _get_time_func(self):
        t_arr = np.linspace(0, self.t_max, len(self.d_array))
        return interp1d(t_arr, self.d_array)

    figure = Instance(Figure)

    def _figure_default(self):
        figure = Figure()
        return figure

    update = Button()

    def _update_fired(self):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        x = np.arange(0, self.t_max, self.d_t)
        ax.plot(x, self.time_func(x))
        ax.set_xlabel('time')
        ax.set_ylabel('displacement')
        self.figure.canvas.draw()

    view = View(VGroup(Group(Item('loading_type')),
                       Group(Item('maximum_slip'),
                             Item('number_of_increments')),
                       Group(Item('number_of_cycles'),
                             Item('amplitude_type'),
                             Item('loading_range'), show_border=True, label='Cyclic load inputs'),
                       Group(Item('d_t'),
                             Item('t_max'),
                             Item('k_max'), show_border=True, label='Solver Settings')),
                Group(Item('update', label='Plot Loading senario')),
                Item('figure', editor=MPLFigureEditor(),
                     dock='horizontal', show_label=False))


class BondSlipModel(BMCSTreeNode):

    node_name = Str('pull out simulation')

    tree_node_list = List([])

    def _tree_node_list_default(self):
        print 'NODE', self.material
        return [self.material]

    material = Instance(Material)

    def _material_default(self):
        return Material()

    geometry = Instance(Geometry())

    mats_eval = Instance(MATSEvalFatigue)

    loading_scenario = Instance(LoadingScenario)

    fets_eval = Instance(FETS1D52ULRHFatigue)

    time_stepper = Instance(TStepper)

    time_loop = Instance(TLoop)

    t_record = Array
    U_record = Array
    F_record = Array
    sf_record = Array
    eps_record = List
    sig_record = List
    w_record = List

    figure = Instance(Figure)

    def _figure_default(self):
        figure = Figure()
        return figure

    plot = Button()

    def _plot_fired(self):
        # assign the material parameters
        self.mats_eval.E_b = self.material.E_b
        self.mats_eval.gamma = self.material.gamma
        self.mats_eval.S = self.material.S
        self.mats_eval.tau_pi_bar = self.material.tau_pi_bar
        self.mats_eval.r = self.material.r
        self.mats_eval.K = self.material.K

        # assign the geometry parameters
        self.fets_eval.A_m = self.geometry.A_m
        self.fets_eval.P_b = self.geometry.P_b
        self.fets_eval.A_f = self.geometry.A_f
        self.time_stepper.L_x = self.geometry.L_x

        # assign the parameters for solver and loading_scenario
        self.time_loop.t_max = self.loading_scenario.t_max
        self.time_loop.d_t = self.loading_scenario.d_t
        self.time_loop.k_max = self.loading_scenario.k_max
        self.time_loop.tolerance = self.loading_scenario.tolerance

        self.loading_scenario.maximum_slip = self.loading_scenario[
            0].maximum_slip
        self.loading_scenario.number_of_increments = self.loading_scenario[
            0].number_of_increments

        # assign the bc
        self.time_stepper.bc_list[1].value = 1
        self.time_stepper.bc_list[
            1].time_function = self.loading_scenario.time_func

        self.draw()
        self.time = 1.00
#         self.figure.canvas.draw()

    ax2 = Property()

    @cached_property
    def _get_ax1(self):

        return self.figure.add_subplot(231)

    ax1 = Property()

    @cached_property
    def _get_ax2(self):

        return self.figure.add_subplot(232)

    ax3 = Property()

    @cached_property
    def _get_ax3(self):
        return self.figure.add_subplot(234)

    ax4 = Property()

    @cached_property
    def _get_ax4(self):
        return self.figure.add_subplot(235)

    ax5 = Property()

    @cached_property
    def _get_ax5(self):
        return self.figure.add_subplot(233)

    ax6 = Property()

    @cached_property
    def _get_ax6(self):
        return self.figure.add_subplot(236)

    def draw(self):

        self.U_record, self.F_record, self.sf_record, self.t_record, self.eps_record, self.sig_record, self.w_record = self.time_loop.eval()
        n_dof = 2 * self.time_stepper.domain.n_active_elems + 1

        self.ax1.cla()
        l_po, = self.ax1.plot(self.U_record[:, n_dof], self.F_record[:, n_dof])
        marker_po, = self.ax1.plot(
            self.U_record[-1, n_dof], self.F_record[-1, n_dof], 'ro')
        self.ax1.set_title('pull-out force-displacement curve')

        self.ax2.cla()
        X = np.linspace(
            0, self.time_stepper.L_x, self.time_stepper.n_e_x + 1)
        X_ip = np.repeat(X, 2)[1:-1]
        l_w, = self.ax2.plot(X_ip, self.w_record[-1].flatten())
        self.ax2.set_title('Damage')

        self.ax3.cla()
        X = np.linspace(
            0, self.time_stepper.L_x, self.time_stepper.n_e_x + 1)
        X_ip = np.repeat(X, 2)[1:-1]

        l_sf, = self.ax3.plot(X_ip, self.sf_record[-1, :])
        self.ax3.set_title('shear flow in the bond interface')

        self.ax4.cla()
        U = np.reshape(self.U_record[-1, :], (-1, 2)).T
        l_u0, = self.ax4.plot(X, U[0])
        l_u1, = self.ax4.plot(X, U[1])
        l_us, = self.ax4.plot(X, U[1] - U[0])
        self.ax4.set_title('displacement and slip')

        self.ax5.cla()
        l_eps0, = self.ax5.plot(X_ip, self.eps_record[-1][:, :, 0].flatten())
        l_eps1, = self.ax5.plot(X_ip, self.eps_record[-1][:, :, 2].flatten())
        self.ax5.set_title('strain')

        self.ax6.cla()
        l_sig0, = self.ax6.plot(X_ip, self.sig_record[-1][:, :, 0].flatten())
        l_sig1, = self.ax6.plot(X_ip, self.sig_record[-1][:, :, 2].flatten())
        self.ax6.set_title('stress')

        self.ax2.set_ylim(0, 1)
        self.ax3.set_ylim(np.amin(self.sf_record), np.amax(self.sf_record))
        self.ax4.set_ylim(np.amin(self.U_record), np.amax(self.U_record))
        self.ax6.set_ylim(np.amin(self.sig_record), np.amax(self.sig_record))

        self.figure.canvas.draw()

    time = Range(0.00, 1.00, value=1.00)

    @on_trait_change('time')
    def draw_t(self):
        idx = (np.abs(self.time * max(self.t_record) - self.t_record)).argmin()
        n_dof = 2 * self.time_stepper.domain.n_active_elems + 1

        self.ax1.cla()
        l_po, = self.ax1.plot(self.U_record[:, n_dof], self.F_record[:, n_dof])
        marker_po, = self.ax1.plot(
            self.U_record[idx, n_dof], self.F_record[idx, n_dof], 'ro')
        self.ax1.set_title('pull-out force-displacement curve')

        self.ax2.cla()
        X = np.linspace(
            0, self.time_stepper.L_x, self.time_stepper.n_e_x + 1)
        X_ip = np.repeat(X, 2)[1:-1]
        l_w, = self.ax2.plot(X_ip, self.w_record[idx].flatten())
        self.ax2.set_title('Damage')

        self.ax3.cla()
        X = np.linspace(
            0, self.time_stepper.L_x, self.time_stepper.n_e_x + 1)
        X_ip = np.repeat(X, 2)[1:-1]
        l_sf, = self.ax3.plot(X_ip, self.sf_record[idx, :])
        self.ax3.set_title('shear flow in the bond interface')

        self.ax4.cla()
        U = np.reshape(self.U_record[idx, :], (-1, 2)).T
        l_u0, = self.ax4.plot(X, U[0])
        l_u1, = self.ax4.plot(X, U[1])
        l_us, = self.ax4.plot(X, U[1] - U[0])
        self.ax4.set_title('displacement and slip')

        self.ax5.cla()
        l_eps0, = self.ax5.plot(X_ip, self.eps_record[idx][:, :, 0].flatten())
        l_eps1, = self.ax5.plot(X_ip, self.eps_record[idx][:, :, 2].flatten())
        self.ax5.set_title('strain')

        self.ax6.cla()
        l_sig0, = self.ax6.plot(X_ip, self.sig_record[idx][:, :, 0].flatten())
        l_sig1, = self.ax6.plot(X_ip, self.sig_record[idx][:, :, 2].flatten())
        self.ax6.set_title('stress')

        self.ax2.set_ylim(0, 1)
        self.ax3.set_ylim(np.amin(self.sf_record), np.amax(self.sf_record))
        self.ax4.set_ylim(np.amin(self.U_record), np.amax(self.U_record))
        self.ax6.set_ylim(np.amin(self.sig_record), np.amax(self.sig_record))

        self.figure.canvas.draw()

    view = View(
        Group(
            Item('time', label='t/T_max'),
        ),
        dock='vertical',
        resizable=True,
        height=0.9, width=1.0
    )

if __name__ == '__main__':

    ts = TStepper()
    n_dofs = ts.domain.n_dofs
    loading_scenario = LoadingScenario()

    ts.bc_list = [BCDof(var='u', dof=0, value=0.0), BCDof(
        var='u', dof=n_dofs - 1, time_function=loading_scenario.time_func)]
    tl = TLoop(ts=ts)

    loading_scenario = LoadingScenario()

    window = BondSlipModel(
        mats_eval=ts.mats_eval,
        fets_eval=ts.fets_eval,
        time_stepper=ts,
        time_loop=tl, loading_scenario=loading_scenario)
#     window.draw()
#
    window.configure_traits()
