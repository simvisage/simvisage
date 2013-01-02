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
    Dict, Button, Bool, Enum, Event, implements, DelegatesTo, \
    Callable

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

import numpy as np

from numpy import \
    loadtxt, argmax, polyfit, poly1d, frompyfunc, dot

from etsproxy.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter

from mathkit.mfn import MFnLineArray
from mathkit.mfn.mfn_line.mfn_matplotlib_editor import \
    MFnMatplotlibEditor

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

from mathkit.array.smoothing import smooth

from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp

from matresdev.db.matdb.trc.fabric_layout \
    import FabricLayOut

from matresdev.db.matdb.trc.concrete_mixture \
    import ConcreteMixture

from matresdev.db.matdb.trc.composite_cross_section import \
    CompositeCrossSection, plain_concrete

from matresdev.db.exdb.loadtxt_bending import loadtxt_bending

from matresdev.db.simdb import \
    SimDB

# Access to the toplevel directory of the database
#
simdb = SimDB()

#class ExpBendingTestThreePoint(ExType):
class ExpBT3PT(ExType):
    '''Experiment: Bending Test Three Point
    '''
#    label = Str('three point bending test')

    implements(IExType)

    file_ext = 'raw'

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

    length = Float(0.42, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)
    width = Float(0.1, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)
    thickness = Float(0.02, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)

    # age of the concrete at the time of testing
    age = Int(43, unit = 'd', input = True, table_field = True,
                             auto_set = False, enter_set = True)
    loading_rate = Float(4.0, unit = 'mm/min', input = True, table_field = True,
                            auto_set = False, enter_set = True)

    #--------------------------------------------------------------------------
    # composite cross section
    #--------------------------------------------------------------------------

    ccs = Instance(CompositeCrossSection)
    def _ccs_default(self):
        '''default settings
        '''
#        fabric_layout_key = 'MAG-07-03'
#        fabric_layout_key = '2D-02-06a'
#        fabric_layout_key = '2D-05-11'
        fabric_layout_key = '2D-09-12'
#        concrete_mixture_key = 'PZ-0708-1'
#        concrete_mixture_key = 'FIL-10-09'
        concrete_mixture_key = 'barrelshell'
#        orientation_fn_key = 'all0'
        orientation_fn_key = 'all90'                                           
#        orientation_fn_key = '90_0'
        n_layers = 8
        s_tex_z = 0.020 / (n_layers + 1)
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

    # flag distinguishes weather data from a displacement gauge is available
    # stored in a separate ASC-file with a corresponding file name
    # 
    flag_dat_file = 'true'

    def _read_data_array(self):
        ''' Read the experiment data.
        '''
        if exists(self.data_file):

            print 'READ FILE'
            # change the file name dat with asc  
            file_split = self.data_file.split('.')

            file_name = file_split[0] + '.csv'
            if not os.path.exists(file_name):

                file_name = file_split[0] + '.raw'
                if not os.path.exists(file_name):
                    raise IOException, 'file %s does not exist' % file_name

            _data_array = loadtxt_bending(file_name)

            print 'file_name', file_name

            file_name = file_split[0] + '.ASC'
            if not os.path.exists(file_name):
                print 'NOTE: no data from displacement gauge is available (no .ASC file); use machine displacement as value for center deflection'
                self.flag_dat_file = 'false'
            else: 
                print 'NOTE: data from displacement gauge for center deflection is available (stored in .ASC-file)!'

            self.data_array = _data_array


    names_and_units = Property
    @cached_property
    def _get_names_and_units(self):
        ''' Set the names and units of the measured data.
        '''
        names = ['w_raw', 'eps_c_raw', 'F_raw']
        units = ['mm', '1*E-3', 'N']
        return names, units

    F = Property(Array('float_'), depends_on = 'input_change')
    @cached_property
    def _get_F(self):
        # convert data in .raw-file from 'N' to 'kN' and change sign to positive value
        #
        return -0.001 * self.F_raw

    eps_c = Property(Array('float_'), depends_on = 'input_change')
    @cached_property
    def _get_eps_c(self):
        # convert the promile values to dimensionless
        #
        return 0.001 * self.eps_c_raw

    elastomer_law = Property(depends_on = 'input_change')
    @cached_property
    def _get_elastomer_law(self):

        elastomer_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE', 'elastomer_f-w.raw')
        _data_array_elastomer = loadtxt_bending(elastomer_path)

        # force [kN]:
        #
        F_elastomer = -0.001 * _data_array_elastomer[:, 2].flatten()

        # displacement [mm]:
        #
        w_elastomer = -1.0 * _data_array_elastomer[:, 0].flatten()

        mfn_displacement_elastomer = MFnLineArray(xdata = F_elastomer, ydata = w_elastomer)
        return frompyfunc(mfn_displacement_elastomer.get_value, 1, 1)

    w_wo_elast = Property(depends_on = 'input_change')
    @cached_property
    def _get_w_wo_elast(self):
        # use the machine displacement for the center displacement:
        # subtract the deformation of the elastomer cushion between the cylinder
        # and change sign in positive values for vertical displacement [mm]
        #
        w = -1.0 * self.w_raw
        return w - self.elastomer_law(self.F)

    M_kNm = Property(Array('float_'), depends_on = 'input_change')
    @cached_property
    def _get_M_kNm(self):
        return self.F * self.length / 4.0

    Fw_asc = Property(depends_on = 'input_change')
    @cached_property
    def _get_Fw_asc(self):
        # if data from displacement gauge is available (stored in .ASC file use this data for center displacement
        #
        # change the file name dat with asc  
        file_split = self.data_file.split('.')
        file_name = file_split[0] + '.ASC'
        return np.loadtxt(file_name,
                            delimiter = ';',
                            usecols = [1,2])
    
    F_asc = Property(depends_on = 'input_change')
    @cached_property
    def _get_F_asc(self):
        # force [kN] stored in the .ASC-file:
        #
        F_asc = self.Fw_asc[:, 0].flatten()
        # change sign to positive value  
        F_asc *= -1.
        return F_asc

    w_asc = Property(depends_on = 'input_change')
    @cached_property
    def _get_w_asc(self):
        # displacement [mm] stored in the .ASC-file:
        #
        w_asc = self.Fw_asc[:, 1].flatten()
        # change sign to positive value for deflection 
        w_asc *= -1.
        return w_asc


    def process_source_data(self):
        '''read in the measured data from file and assign
        attributes after array processing.
        '''
        super(ExpBT3PT, self).process_source_data()

        # access the derived arrays to initiate their processing
        self.M_kNm
        self.w_wo_elast
        # only if separate ASC.-file with force-displacement data from displacement gauge is available
        if self.flag_dat_file == 'true':
            print 'XXX', self.flag_dat_file
            self.F_asc
            self.w_asc

    #--------------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------------

    plot_templates = {'force / deflection (machine displacement)' : '_plot_force_deflection_machine_displ',
                      'force / deflection (displacement gauge)' : '_plot_force_deflection_gauge_displ',
                      'smoothed force / deflection' : '_plot_smoothed_force_deflection',
                      'moment / eps_c' : '_plot_moment_eps_c',
                      'smoothed moment / eps_c' : '_plot_smoothed_moment_eps_c'
                     }

    default_plot_template = 'force / deflection (displacement gauge)'

    def _plot_force_deflection_machine_displ(self, axes):
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'

        xdata = self.w_wo_elast
        ydata = self.F

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata
                       # color = c, linewidth = w, linestyle = s 
                       )

    def _plot_force_deflection_gauge_displ(self, axes):
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'

        xdata = self.w_asc
        ydata = self.F_asc

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata
                       # color = c, linewidth = w, linestyle = s 
                       )


    def _plot_moment_eps_c(self, axes):
        xkey = 'compressive strain [1*E-3]'
        ykey = 'moment [kNm]'
        # NOTE: processed data returns positive values for force and displacement
        #
        xdata = self.eps_c
        ydata = self.M_kNm

        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata
                       # color = c, linewidth = w, linestyle = s 
                       )

    n_fit_window_fraction = Float(0.1)

    def _plot_smoothed_force_deflection(self, axes):
#        if self.flag_dat_file == 'false':
#            w = self.w_wo_elast
#        else:
#            w = self.w_dat

        w = self.w_wo_elast

        # get the index of the maximum stress
        max_force_idx = argmax(self.F)
        # get only the ascending branch of the response curve
        f_asc = self.F[:max_force_idx + 1]
        w_asc = w[:max_force_idx + 1]

        n_points = int(self.n_fit_window_fraction * len(w_asc))
        f_smooth = smooth(f_asc, n_points, 'flat')
        w_smooth = smooth(w_asc, n_points, 'flat')

        axes.plot(w_smooth, f_smooth, color = 'blue', linewidth = 2)

#        secant_stiffness_w10 = ( f_smooth[10] - f_smooth[0] ) / ( w_smooth[10] - w_smooth[0] )
#        w0_lin = array( [0.0, w_smooth[10] ], dtype = 'float_' )
#        f0_lin = array( [0.0, w_smooth[10] * secant_stiffness_w10 ], dtype = 'float_' )

        #axes.plot( w0_lin, f0_lin, color = 'black' )

    M_eps_c_smoothed = Property(depends_on = 'input_change')
    @cached_property
    def _get_M_eps_c_smoothed(self):
        # get the index of the maximum stress
        max_idx = argmax(self.M_kNm)
        # get only the ascending branch of the response curve
        m_asc = self.M_kNm[:max_idx + 1]
        eps_c_asc = self.eps_c[:max_idx + 1]

        n_points = int(self.n_fit_window_fraction * len(eps_c_asc))
        m_smoothed = smooth(m_asc, n_points, 'flat')
        eps_c_smoothed = smooth(eps_c_asc, n_points, 'flat')
        return m_smoothed, eps_c_smoothed

    eps_c_smoothed = Property
    def _get_eps_c_smoothed(self):
        return self.M_eps_c_smoothed[1]

    M_smoothed = Property
    def _get_M_smoothed(self):
        return self.M_eps_c_smoothed[0]

    def _plot_smoothed_moment_eps_c(self, axes):
        axes.plot(self.eps_c_smoothed, self.M_smoothed, color = 'blue', linewidth = 2)

    #--------------------------------------------------------------------------------
    # view
    #--------------------------------------------------------------------------------

    traits_view = View(VGroup(
                         Group(
                              Item('length', format_str = "%.3f"),
                              Item('width', format_str = "%.3f"),
                              Item('thickness', format_str = "%.3f"),
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
    ex = ExRunClassExt(klass = ExpBT3PT)
    ex.configure_traits()
