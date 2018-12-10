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

from numpy import \
    loadtxt, argmax, polyfit, poly1d, frompyfunc, dot

import numpy as np


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
class ExpEM(ExType):
    '''Experiment: Elastic Modulus Test
    '''
#    label = Str('three point bending test')

    implements(IExType)

    file_ext = 'TRA'

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

    edge_length = Float(0.06, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)
    height = Float(0.12, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)
    gauge_length = Float(0.10, unit = 'm', input = True, table_field = True,
                           auto_set = False, enter_set = True)

    # age of the concrete at the time of testing
    age = Int(39, unit = 'd', input = True, table_field = True,
                             auto_set = False, enter_set = True)
    loading_rate = Float(0.6, unit = 'MPa/s', input = True, table_field = True,
                            auto_set = False, enter_set = True)

    #--------------------------------------------------------------------------
    # composite cross section
    #--------------------------------------------------------------------------

    ccs = Instance(CompositeCrossSection)
    def _ccs_default(self):
        '''default settings correspond to 
        setup '7u_MAG-07-03_PZ-0708-1'
        '''
#        fabric_layout_key = 'MAG-07-03'
#        fabric_layout_key = '2D-02-06a'
        fabric_layout_key = '2D-05-11'
#        concrete_mixture_key = 'PZ-0708-1'
        concrete_mixture_key = 'FIL-10-09'
        orientation_fn_key = 'all0'
#        orientation_fn_key = 'all90'                                           
#        orientation_fn_key = '90_0'
        n_layers = 12
        s_tex_z = 0.060 / (n_layers + 1)
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

    def _read_data_array(self):
        ''' Read the experiment data. 
        '''
        print('READ FILE')
        _data_array = np.loadtxt( self.data_file, delimiter = ';', skiprows = 2, usecols = (4,5,10,17,20) )
        self.data_array = _data_array

    names_and_units = Property
    @cached_property
    def _get_names_and_units(self):
        ''' Set the names and units of the measured data.
        '''
        names = ['delta_u1', 'delta_u2', 'w', 'F', 'time']
        units = ['mm', 'mm', 'mm', 'kN', 's']
        return names, units

#0#"Arbeit";
#1#"Dehn. abs";
#2#"Dehn. abs (2. Kanal)";
#3#"Dehnung";

#4#"Dehnung (1. Kanal)";
#5#"Dehnung (2. Kanal)";

#6#"Dehnung nominell";
#7#"DeltaE";
#8#"E1 korr";
#9#"E2 korr";

#10#"Kolbenweg";

#11#"Kolbenweg abs.";
#12#"Lastrahmen";
#13#"LE-Kanal";
#14#"PrXXXfzeit";
#15#"Querdehnung";
#16#"S korr";
#17#"Standardkraft";
#18#"Vorlaufzeit";
#19#"Weg";
#20#"Zeit";
#21#"Zyklus"

#
#"Nm";
#"mm";
#"mm";
#"mm";
#"mm";
#"mm";
#"mm";
#"%";
#"mm";
#" ";
#"mm";
#"mm";
#" ";
#"mm";
#"s";
#"mm";
#" ";
#"N";
#"s";
#"mm";
#"s";
#" "

    #--------------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------------

    plot_templates = {'force / displacement' : '_plot_force_displacement',
                      'stress / strain'      : '_plot_stress_strain',
                      'stress / time'      : '_plot_stress_time',
                     }

    default_plot_template = 'force / displacement'

    def _plot_force_displacement(self, axes):
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        xdata = self.w 
        ydata = self.F / 1000. # convert from [N] to [kN]
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata
                       # color = c, linewidth = w, linestyle = s 
                       )

    def _plot_stress_strain(self, axes):
        sig = ( self.F / 1000000.) / (self.edge_length**2) # convert from [N] to [MN]
        eps1 = ( self.delta_u1 / 1000.) / self.gauge_length
        eps2 = ( self.delta_u2 / 1000.) / self.gauge_length
        eps_m = (eps1 + eps2) / 2.
        axes.plot(eps_m, sig, color = 'blue', linewidth = 2)

    def _plot_stress_time(self, axes):
        sig = ( self.F / 1000000.) / (self.edge_length**2) # convert from [N] to [MN]
        axes.plot( self.time, sig, color = 'blue', linewidth = 2)

    def _plot_displ_time(self, axes):
        axes.plot( self.time, self.displ, color = 'blue', linewidth = 2)


    #--------------------------------------------------------------------------------
    # view
    #--------------------------------------------------------------------------------

    traits_view = View(VGroup(
                         Group(
                              Item('edge_length', format_str = "%.3f"),
                              Item('height', format_str = "%.3f"),
                              Item('gauge_length', format_str = "%.3f"),
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
    ex = ExRunClassExt(klass = ExpEM)
    ex.configure_traits()
