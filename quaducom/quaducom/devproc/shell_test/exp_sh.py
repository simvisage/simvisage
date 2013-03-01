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
# Created on Feb 15, 2010 by: rch

from etsproxy.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements, DelegatesTo

from etsproxy.util.home_directory import \
    get_home_directory

from etsproxy.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor, EnumEditor, Handler, FileEditor, VSplit, Group, \
    HGroup, Spring

## overload the 'get_label' method from 'Item' to display units in the label
from util.traits.ui.item import \
    Item

from etsproxy.traits.ui.table_column import \
    ObjectColumn

from etsproxy.traits.ui.menu import \
    OKButton, CancelButton

from etsproxy.traits.ui.tabular_adapter \
    import TabularAdapter

from util.traits.editors.mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure

import os

import csv

from numpy import array, fabs, where, copy, ones, argsort

from numpy import \
    loadtxt, argmax, polyfit, poly1d, frompyfunc, dot

from etsproxy.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter

#-- Tabular Adapter Definition -------------------------------------------------

from string import replace
from os.path import exists

#-----------------------------------------------------------------------------------
# ExDesignReader
#-----------------------------------------------------------------------------------
from etsproxy.traits.ui.file_dialog  \
    import open_file, FileInfo, TextInfo, ImageInfo

from etsproxy.traits.ui.api \
    import View, Item, TabularEditor, VGroup, HGroup

from etsproxy.traits.ui.tabular_adapter \
    import TabularAdapter

from matresdev.db.exdb.ex_type import ExType
from matresdev.db.exdb.i_ex_type import IExType

data_file_editor = FileEditor(filter = ['*.DAT'])

from mathkit.array.smoothing import smooth

from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp

from matresdev.db.matdb.trc.fabric_layout \
    import FabricLayOut

from matresdev.db.matdb.trc.concrete_mixture \
    import ConcreteMixture

from matresdev.db.matdb.trc.composite_cross_section import \
    CompositeCrossSection, plain_concrete

class ExpSH(ExType):
    '''Experiment: Shell Test
    '''
#    label = Str('slab test')
    implements(IExType)

    #--------------------------------------------------------------------
    # register a change of the traits with metadata 'input'
    #--------------------------------------------------------------------

    input_change = Event
    @on_trait_change('+input, ccs.input_change, +ironing_param')
    def _set_input_change(self):
        self.input_change = True

    #--------------------------------------------------------------------------------
    # specify inputs:
    #--------------------------------------------------------------------------------

    edge_length = Float(1.25, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)
    thickness = Float(0.03, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)

    # age of the concrete at the time of testing
    age = Int(28, unit = 'd', input = True, table_field = True,
                             auto_set = False, enter_set = True)
    loading_rate = Float(0.50, unit = 'mm/min', input = True, table_field = True,
                            auto_set = False, enter_set = True)

    #--------------------------------------------------------------------------
    # composite cross section
    #--------------------------------------------------------------------------

    ccs = Instance(CompositeCrossSection)
    def _ccs_default(self):
        '''default settings correspond to 
        setup '7u_MAG-07-03_PZ-0708-1'
        '''
        fabric_layout_key = 'MAG-07-03'
#        fabric_layout_key = '2D-02-06a'
        concrete_mixture_key = 'PZ-0708-1'
#        concrete_mixture_key = 'FIL-10-09'
#        orientation_fn_key = 'all0'                                           
        orientation_fn_key = '90_0'
        n_layers = 10
        s_tex_z = 0.030 / (n_layers + 1)
        ccs = CompositeCrossSection (
                    fabric_layup_list = [
                            plain_concrete(s_tex_z * 0.5),
                            FabricLayUp (
                                   n_layers = n_layers,
                                   orientation_fn_key = orientation_fn_key,
                                   s_tex_z = s_tex_z,
                                   fabric_layout_key = fabric_layout_key
                                   ),
                            plain_concrete(s_tex_z * 0.5)
                                        ],
                    concrete_mixture_key = concrete_mixture_key
                    )
        return ccs

    #--------------------------------------------------------------------------
    # Get properties of the composite 
    #--------------------------------------------------------------------------

    # E-modulus of the composite at the time of testing 
    E_c = Property(Float, unit = 'MPa', depends_on = 'input_change', table_field = True)
    def _get_E_c(self):
        return self.ccs.get_E_c_time(self.age)

    # E-modulus of the composite after 28 days
    E_c28 = DelegatesTo('ccs', listenable = False)

    # reinforcement ration of the composite 
    rho_c = DelegatesTo('ccs', listenable = False)

    #--------------------------------------------------------------------------------
    # define processing
    #--------------------------------------------------------------------------------

    # put this into the ironing procedure processor
    #
    jump_rtol = Float(0.2,
                      auto_set = False, enter_set = True,
                      ironing_param = True)


    data_array_ironed = Property(Array(float),
                                  depends_on = 'data_array, +ironing_param, +axis_selection')
    @cached_property
    def _get_data_array_ironed(self):
        '''remove the jumps in the displacement curves 
        due to resetting the displacement gauges. 
        '''
        print '*** curve ironing activated ***'

        # each column from the data array corresponds to a measured parameter 
        # e.g. displacement at a given point as function of time u = f(t))
        #
        data_array_ironed = copy(self.data_array)

        for idx in range(self.data_array.shape[1]):

            # use ironing method only for columns of the displacement gauges.
            #
            if self.names_and_units[0][ idx ] != 'Kraft' and \
                self.names_and_units[0][ idx ] != 'Bezugskanal' and \
                self.names_and_units[0][ idx ] != 'D_1o'and \
                self.names_and_units[0][ idx ] != 'D_2o'and \
                self.names_and_units[0][ idx ] != 'D_3o'and \
                self.names_and_units[0][ idx ] != 'D_4o'and \
                self.names_and_units[0][ idx ] != 'D_1u'and \
                self.names_and_units[0][ idx ] != 'D_2u'and \
                self.names_and_units[0][ idx ] != 'D_3u'and \
                self.names_and_units[0][ idx ] != 'D_4u'and \
                self.names_and_units[0][ idx ] != 'W30_vo'and \
                self.names_and_units[0][ idx ] != 'W30_hi':

                # 1d-array corresponding to column in data_array
                data_arr = copy(data_array_ironed[:, idx])

                # get the difference between each point and its successor
                jump_arr = data_arr[1:] - data_arr[0:-1]

                # get the range of the measured data 
                data_arr_range = max(data_arr) - min(data_arr)

                # determine the relevant criteria for a jump
                # based on the data range and the specified tolerances:
                jump_crit = self.jump_rtol * data_arr_range

                # get the indexes in 'data_column' after which a 
                # jump exceeds the defined tolerance criteria
                jump_idx = where(fabs(jump_arr) > jump_crit)[0]

                print 'number of jumps removed in data_arr_ironed for', self.names_and_units[0][ idx ], ': ', jump_idx.shape[0]
                # glue the curve at each jump together
                for jidx in jump_idx:
                    # get the offsets at each jump of the curve
                    shift = data_arr[jidx + 1] - data_arr[jidx]
                    # shift all succeeding values by the calculated offset
                    data_arr[jidx + 1:] -= shift

                data_array_ironed[:, idx] = data_arr[:]

        return data_array_ironed


    @on_trait_change('+ironing_param')
    def process_source_data(self):
        '''read in the measured data from file and assign
        attributes after array processing.        
        
        NOTE: if center displacement gauge ('WA_M') is missing the measured 
        displacement of the cylinder ('Weg') is used instead.
        A minor mistake is made depending on how much time passes
        before the cylinder has contact with the slab.
        '''
        print '*** process source data ***'

        self._read_data_array()
        # curve ironing:
        #
        self.processed_data_array = self.data_array_ironed

        # set attributes:
        #
        self._set_array_attribs()



    #--------------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------------

    plot_templates = {'force / time' : '_plot_force_time',
                      'force / deflections center' : '_plot_force_deflection_center',
                      'force / smoothed deflections center' : '_plot_smoothed_force_deflection_center',
                      'force / deflections' : '_plot_force_deflections',
                      'force / eps_t' : '_plot_force_eps_t',
                      'force / eps_c' : '_plot_force_eps_c',
                      'time / deflections' : '_plot_time_deflections'
                     }
    default_plot_template = 'force / time'

    def _plot_force_time(self, axes):

        xkey = 'time [s]'
        ykey = 'force [kN]'
        xdata = self.Bezugskanal
        ydata = self.Kraft

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata
                       # color = c, linewidth = w, linestyle = s 
                       )

    def _plot_force_deflection_center(self, axes):

        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        xdata = -self.WD4
        ydata = self.Kraft

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata
                       # color = c, linewidth = w, linestyle = s 
                       )

    n_fit_window_fraction = Float(0.1)

    def _plot_smoothed_force_deflection_center(self, axes):

        xkey = 'deflection [mm]'
        ykey = 'force [kN]'

        # get the index of the maximum stress
        max_force_idx = argmax(self.Kraft)
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:max_force_idx + 1]
        w_asc = -self.WD4[:max_force_idx + 1]

        f_max = f_asc[-1]
        w_max = w_asc[-1]

        n_points = int(self.n_fit_window_fraction * len(w_asc))
        f_smooth = smooth(f_asc, n_points, 'flat')
        w_smooth = smooth(w_asc, n_points, 'flat')

        axes.plot(w_smooth, f_smooth, color = 'blue', linewidth = 2)

        secant_stiffness_w10 = (f_smooth[10] - f_smooth[0]) / (w_smooth[10] - w_smooth[0])
        w0_lin = array([0.0, w_smooth[10] ], dtype = 'float_')
        f0_lin = array([0.0, w_smooth[10] * secant_stiffness_w10 ], dtype = 'float_')

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))


    def _plot_force_deflections(self, axes):

        xkey = 'displacement [mm]'
        ykey = 'force [kN]'

        # get the index of the maximum stress
        max_force_idx = argmax(self.Kraft)
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:max_force_idx + 1]

        WD1 = -self.WD1[:max_force_idx + 1]
        WD2 = -self.WD2[:max_force_idx + 1]
        WD3 = -self.WD3[:max_force_idx + 1]
        WD4 = -self.WD4[:max_force_idx + 1]
        WD5 = -self.WD5[:max_force_idx + 1]
        WD6 = -self.WD6[:max_force_idx + 1]
        WD7 = -self.WD7[:max_force_idx + 1]

        axes.plot(WD1, f_asc, color = 'blue', linewidth = 1 )
        axes.plot(WD2, f_asc, color = 'green', linewidth = 1 )
        axes.plot(WD3, f_asc, color = 'red', linewidth = 1 )
        axes.plot(WD4, f_asc, color = 'black', linewidth = 3 )
        axes.plot(WD5, f_asc, color = 'red', linewidth = 1 )
        axes.plot(WD6, f_asc, color = 'green', linewidth = 1 )
        axes.plot(WD7, f_asc, color = 'blue', linewidth = 1 )

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))


    def _plot_time_deflections(self, axes):

        xkey = 'time [s]'
        ykey = 'deflections [mm]'

        time = self.Bezugskanal

        WD1 = -self.WD1
        WD2 = -self.WD2
        WD3 = -self.WD3
        WD4 = -self.WD4
        WD5 = -self.WD5
        WD6 = -self.WD6
        WD7 = -self.WD7

        axes.plot(time, WD1, color = 'blue', linewidth = 1 )
        axes.plot(time, WD2, color = 'green', linewidth = 1 )
        axes.plot(time, WD3, color = 'red', linewidth = 1 )
        axes.plot(time, WD4, color = 'black', linewidth = 3 )
        axes.plot(time, WD5, color = 'red', linewidth = 1 )
        axes.plot(time, WD6, color = 'green', linewidth = 1 )
        axes.plot(time, WD7, color = 'blue', linewidth = 1 )

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))


    def _plot_force_eps_c(self, axes):

        xkey = 'force [kN]'
        ykey = 'strain [mm/m]'

        # get the index of the maximum stress
        max_force_idx = argmax(self.Kraft)
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:max_force_idx + 1]
        
        D_1u = -self.D_1u[:max_force_idx + 1]
        D_2u = -self.D_2u[:max_force_idx + 1]
        D_3u = -self.D_3u[:max_force_idx + 1]
        D_4u = -self.D_4u[:max_force_idx + 1]
        
        D_1o = -self.D_1o[:max_force_idx + 1]
        D_2o = -self.D_2o[:max_force_idx + 1]
        D_3o = -self.D_3o[:max_force_idx + 1]
        D_4o = -self.D_4o[:max_force_idx + 1]

        axes.plot(f_asc, D_1u, color = 'blue', linewidth = 1 )
        axes.plot(f_asc, D_2u, color = 'blue', linewidth = 1 )
        axes.plot(f_asc, D_3u, color = 'blue', linewidth = 1 )
        axes.plot(f_asc, D_4u, color = 'blue', linewidth = 1 )

        axes.plot(f_asc, D_1o, color = 'red', linewidth = 1 )
        axes.plot(f_asc, D_2o, color = 'red', linewidth = 1 )
        axes.plot(f_asc, D_3o, color = 'red', linewidth = 1 )
        axes.plot(f_asc, D_4o, color = 'red', linewidth = 1 )

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))


    def _plot_force_eps_t(self, axes):

        xkey = 'strain [mm/m]'
        ykey = 'force [kN]'

        # get the index of the maximum stress
        max_force_idx = argmax(self.Kraft)
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:max_force_idx + 1]

        eps_vo = -self.W30_vo[:max_force_idx + 1] / 300. 
        eps_hi = -self.W30_hi[:max_force_idx + 1] / 300. 

        axes.plot(eps_vo, f_asc, color = 'blue', linewidth = 1 )
        axes.plot(eps_hi, f_asc, color = 'green', linewidth = 1 )

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))


    #--------------------------------------------------------------------------------
    # view
    #--------------------------------------------------------------------------------

    traits_view = View(VGroup(
                         Group(
                              Item('jump_rtol', format_str = "%.4f"),
                              label = 'curve_ironing'
                              ),
                         Group(
                              Item('thickness', format_str = "%.3f"),
                              Item('edge_length', format_str = "%.3f"),
                              label = 'geometry'
                              ),
                         Group(
                              Item('loading_rate'),
                              Item('age'),
                              label = 'loading rate and age'
                              ),
                         Group(
                              Item('E_c', show_label = True, style = 'readonly', format_str = "%.0f"),
                              Item('ccs@', show_label = False),
                              label = 'composite cross section'
                              )
                         ),
                        scrollable = True,
                        resizable = True,
                        height = 0.8,
                        width = 0.6
                        )

if __name__ == '__main__':

    from matresdev.db.exdb.ex_run_table import ExRunClassExt
    ex = ExRunClassExt(klass = ExpSH)
    ex.configure_traits()

