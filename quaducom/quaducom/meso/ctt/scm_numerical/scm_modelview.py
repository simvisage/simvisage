

from util.traits.editors.mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
from enthought.traits.ui.api import \
    View, Item, VGroup, HGroup, ModelView, HSplit, VSplit
from enthought.traits.api import \
    Instance, Enum, Bool, on_trait_change, Int, Event, Array, Tuple, List, \
    Float, HasTraits, Float, Property, Button
from enthought.traits.ui.menu import OKButton, CancelButton
from numpy import array, insert, mean, sort

from scm_model import SCMModel


#--------------------------------------------------------------------------
# MODEL_VIEW
#--------------------------------------------------------------------------

class SCMModelView (ModelView):

    model = Instance(SCMModel)
    def _model_default(self):
        return SCMModel()

    screen3 = Enum('no_of_cracks', 'cbs_sum', 'strain_profile', 'simulated_sigma_eps',
                   'random_matrix_strenght', modified = True)
    screen2 = Enum('strain_profile', 'analytical', 'simulated_sigma_eps',
                   'random_matrix_strenght', modified = True)
    screen1 = Enum('simulated_sigma_eps', 'analytical', 'strain_profile',
                   'random_matrix_strenght', modified = True)
    screen4 = Enum('random_matrix_strenght', 'strain_profile', 'simulated_sigma_eps',
                   'analytical', modified = True)


    launch = Button('launch')
    max_stress = Float(5000.)
    points = Int(20)
    clear_1 = Bool(True)
    clear_2 = Bool(True)
    clear_3 = Bool(True)
    clear_4 = Bool(True)


    figure1 = Instance(Figure)
    def _figure1_default(self):
        figure1 = Figure(facecolor = 'white')
        figure1.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure1

    figure2 = Instance(Figure)
    def _figure2_default(self):
        figure2 = Figure(facecolor = 'white')
        figure2.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure2

    figure3 = Instance(Figure)
    def _figure3_default(self):
        figure3 = Figure(facecolor = 'white')
        figure3.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure3

    figure4 = Instance(Figure)
    def _figure4_default(self):
        figure4 = Figure(facecolor = 'white')
        figure4.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure4

    data_changed = Event

    def _launch_fired(self):
        self.model.launch(self.max_stress, self.points)

    @on_trait_change('model.+ctrl_param, model.+modified, +modified')
    def refresh(self):

        scm = self.model

        figure1 = self.figure1
        if self.clear_1 == True:
            figure1.clear()
        axes1 = figure1.gca()

        figure2 = self.figure2
        if self.clear_2 == True:
            figure2.clear()
        axes2 = figure2.gca()

        figure3 = self.figure3
        if self.clear_3 == True:
            figure3.clear()
        axes3 = figure3.gca()

        figure4 = self.figure4
        if self.clear_4 == True:
            figure4.clear()
        axes4 = figure4.gca()
        scm.load_step()

        # ----------------------------------------------
        # STRAIN PROFILES

        def strain_profile(axes):

            axes.set_title('Dehnungsprofile', weight = 'bold')
            axes.set_xlabel('Position [mm]', weight = 'semibold')
            axes.set_ylabel('Dehnung [-]', weight = 'semibold')
            axes.plot(scm.x_arr, scm.eps_m_x, lw = 2, color = 'blue', label = 'Matrix')
            axes.plot(scm.x_arr, scm.eps_r_x, lw = 2, color = 'red', label = 'Bewehrung')
            axes.set_ylim(0, 0.001)
            #axes.legend( loc = 'upper right' )

        # ----------------------------------------------

        # ----------------------------------------------
        # SIMULATED STOCHASTIC CRACKING

        def simulated_sigma_eps(axes):
            #scm._get_cbs_responses()
            if len(scm.cracks) == 0:
                cs = 10e9
            else:
                cr = sort(scm.cracks)
                cs = mean((cr - insert(cr[:-1], [0], 0))[1:])

            eps_c, sigma_c = scm.sig_eps

            axes.set_title('Spannungs-Dehnungs Diagramm', weight = 'bold')
            axes.set_xlabel('Dehnung $\epsilon$', weight = 'semibold')
            axes.set_ylabel('Spannung $\sigma$ [N/mm$^2$]', weight = 'semibold')
            axes.plot(eps_c, sigma_c, lw = 2, color = 'red', label = 'simulation, CS = %.1f mm' % cs)
#            axes.plot( eps_c, array( eps_c ) * scm.cb.Kc / scm.cb.Ac, lw = 1, color = 'red',
#                       label = 'uncracked composite' )
            #axes.plot( eps_c, array( eps_c ) * scm.cb.Kr / scm.cb.Ac, lw = 1, color = 'black',
                       #label = 'reinforcement only' )
#            axes.plot( scm.cbs_sum_strain , scm.P_list , label = 'cbs sum',
#                       lw = 2, color = 'green' )
            axes.set_axis_bgcolor(color = 'white')
            axes.ticklabel_format(scilimits = (-3., 4.))
            axes.grid(color = 'gray', linestyle = '--', linewidth = 0.5, alpha = 0.7)
            axes.legend(['integ_s', 'cbs_sum'] , loc = 'lower right')


        # ----------------------------------------------
        # RANDOM MATRIX STRESS PROFILE

        def random_matrix_strenght(axes):
            axes.set_title('Zufallsfeld der Rissfestigkeit $\sigma_c$', weight = 'bold')
            axes.set_xlabel('Position [mm]', weight = 'semibold')
            axes.set_ylabel('Spannung $\sigma_{c,u}$ [N/mm$^2$]', weight = 'semibold')
            axes.plot(scm.x_arr, scm.random_field, color = 'black', label = 'matrix strength')
            axes.plot(scm.x_arr, scm.sigma_m_x, color = 'red', label = 'matrix stress')
            axes.set_ylim(0, max(scm.random_field) * 1.2)
            axes.legend( loc = 'upper right' )


        # ----------------------------------------------
        # No. OF CRACKS

        def no_of_cracks(axes):
            scm.get_no_of_cracks()
            axes.set_title('Rissanzahl', weight = 'bold')
            axes.set_xlabel('Spannung $\sigma_c$ [MPa]', weight = 'semibold')
            axes.set_ylabel('Anzahl an Rissen', weight = 'semibold')

            axes.plot(scm.no_of_cracks[0], scm.no_of_cracks[1],
                       lw = 2, color = 'black', label = 'No. of cracks')
            axes.plot([8.0, 8.0], [0.0, scm.no_of_cracks[1][-1]], ls = 'dashed', lw = 2, color = 'black')



        ##############################################
        '''Axel's addition'''
        ##############################################
        # ----------------------------------------------
        # cbs_sum

#        def cbs_sum( axes ):
#            scm._get_cbs_responses()
#
#            axes.set_title( 'l-d diagram', weight = 'bold' )
#            axes.set_xlabel( 'composite strain', weight = 'semibold' )
#            axes.set_ylabel( 'stress', weight = 'semibold' )
#            axes.plot( scm.cbs_sum_strain , scm.P_list , label = 'cbs sum',
#                       lw = 2, color = 'red' )
#            axes.legend( loc = 'best' )




            #
        # ----------------------------------------------

        eval(self.screen1 + '(axes1)')
        eval(self.screen2 + '(axes2)')
        eval(self.screen3 + '(axes3)')
        eval(self.screen4 + '(axes4)')

        self.data_changed = True

    traits_view = View(
                   HSplit(
                       VGroup(
                             Item('model@', show_label = False, resizable = True),
                             label = 'Material parameters',
                             id = 'scm.viewmodel.model',
                             dock = 'tab',
                             ),
                    VSplit(
                       VGroup(HSplit(VGroup(
                                       VGroup(Item('figure1',
                                                    editor = MPLFigureEditor(),
                                                    resizable = True, show_label = False),
                                             HGroup(Item('screen1'),
                                                    Item('clear_1', label = 'clear before relaunch'),),
                                                    id = 'scm.viewmodel.screen1',
                                                    dock = 'tab',
                                             ),
                                        VGroup(Item('figure2',
                                                      editor = MPLFigureEditor(),
                                                      resizable = True, show_label = False),
                                                HGroup(Item('screen2'),
                                                    Item('clear_2', label = 'clear before relaunch'),
                                                    id = 'scm.viewmodel.screen2',
                                                    dock = 'tab',),
                                                ),
                                        label = 'CDF and strain profile',
                                       id = 'scm.viewmodel.plot_window',
                                       dock = 'tab',
                                            ),
                                    VGroup(
                                       VGroup(Item('figure3',
                                                    editor = MPLFigureEditor(),
                                                    resizable = True, show_label = False),
                                             HGroup(Item('screen3'),
                                                    Item('clear_3', label = 'clear before relaunch'),),
                                             ),
                                        VGroup(
                                               Item('figure4',
                                                      editor = MPLFigureEditor(),
                                                      resizable = True, show_label = False),
                                                HGroup(Item('screen4'),
                                                    Item('clear_4', label = 'clear before relaunch'),),
                                                ),
                                        label = 'SCM_diagram',
                                       id = 'scm.diagram_plot_window',
                                       dock = 'tab',)
                                    )
                             ),
                       HGroup(
                               VGroup(Item('model.applied_force', label = 'current force'),
                                      HGroup(Item('launch', label = 'launch computation'),
                                              Item('max_stress', label = 'max stress'),
                                              Item('points', label = 'discretization points')),
                                      scrollable = False),
                                label = 'plot parameters',
                               id = 'scm.viewmodel.plot_params',
                               group_theme = 'blue',
                             ),
                             id = 'scm.viewmodel.right',
                            ),
                        id = 'scm.viewmodel.splitter',
                    ),
                    title = 'Stochastic Cracking Model',
                    id = 'scm.viewmodel',
                    dock = 'tab',
                    kind = 'live',
                    resizable = True,
                    buttons = [OKButton],
                    height = 0.8, width = 0.8
                    )

if __name__ == '__main__':
    s = SCMModelView(model = SCMModel())
    s.refresh()
    s.configure_traits()
