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

    length = Float(0.46, unit='m', input=True, table_field=True,
                           auto_set=False, enter_set=True)
    width = Float(0.1, unit='m', input=True, table_field=True,
                           auto_set=False, enter_set=True)
    thickness = Float(0.02, unit='m', input=True, table_field=True,
                           auto_set=False, enter_set=True)

    # age of the concrete at the time of testing
    age = Int(33, unit='d', input=True, table_field=True,
                             auto_set=False, enter_set=True)
    loading_rate = Float(4.0, unit='mm/min', input=True, table_field=True,
                            auto_set=False, enter_set=True)

    #--------------------------------------------------------------------------
    # composite cross section
    #--------------------------------------------------------------------------

    ccs = Instance(CompositeCrossSection)
    def _ccs_default(self):
        '''default settings
        '''
#        fabric_layout_key = 'MAG-07-03'
#        fabric_layout_key = '2D-02-06a'
        fabric_layout_key = '2D-05-11'
#        fabric_layout_key = '2D-09-12'
#        concrete_mixture_key = 'PZ-0708-1'
#        concrete_mixture_key = 'FIL-10-09'
        concrete_mixture_key = 'barrelshell'
        orientation_fn_key = 'all0'
#        orientation_fn_key = 'all90'                                           
#        orientation_fn_key = '90_0'
        n_layers = 6
        s_tex_z = 0.020 / (n_layers + 1)
        ccs = CompositeCrossSection (
                    fabric_layup_list=[
                            plain_concrete(s_tex_z * 0.5),
                            FabricLayUp (
                                   n_layers=n_layers,
                                   orientation_fn_key=orientation_fn_key,
                                   s_tex_z=s_tex_z,
                                   fabric_layout_key=fabric_layout_key
                                   ),
                            plain_concrete(s_tex_z * 0.5)
                                        ],
                    concrete_mixture_key=concrete_mixture_key
                    )
        return ccs

    #--------------------------------------------------------------------------
    # Get properties of the composite 
    #--------------------------------------------------------------------------

    # E-modulus of the composite at the time of testing 
    E_c = Property(Float, unit='MPa', depends_on='input_change', table_field=True)
    def _get_E_c(self):
        return self.ccs.get_E_c_time(self.age)

    # E-modulus of the composite after 28 days
    E_c28 = DelegatesTo('ccs', listenable=False)

    # reinforcement ration of the composite 
    rho_c = DelegatesTo('ccs', listenable=False)

    #--------------------------------------------------------------------------------
    # define processing
    #--------------------------------------------------------------------------------

    # flag distinguishes weather data from a displacement gauge is available
    # stored in a separate ASC-file with a corresponding file name
    # 
    flag_ASC_file = Bool(False)

    def _read_data_array(self):
        ''' Read the experiment data.
        '''
        if exists(self.data_file):

            print 'READ FILE'
            file_split = self.data_file.split('.')

            # first check if a '.csv' file exists. If yes use the 
            # data stored in the '.csv'-file and ignore 
            # the data in the '.raw' file!
            #
            file_name = file_split[0] + '.csv'
            if not os.path.exists(file_name):
                file_name = file_split[0] + '.raw'
                if not os.path.exists(file_name):
                    raise IOError, 'file %s does not exist' % file_name

            print 'file_name', file_name
            _data_array = loadtxt_bending(file_name)
            self.data_array = _data_array

            # check if a '.ASC'-file exists. If yes append this information 
            # to the data array.
            #
            file_name = file_split[0] + '.ASC'
            if not os.path.exists(file_name):
                print 'NOTE: no data from displacement gauge is available (no .ASC file)'
                self.flag_ASC_file = False
            else:
                print 'NOTE: additional data from displacement gauge for center deflection is available (.ASC-file loaded)!'
                self.flag_ASC_file = True
                # add data array read in from .ASC-file; the values are assigned by '_set_array_attribs' based on the 
                # read in values in 'names_and_units' read in from the corresponding .DAT-file
                #
                self.data_array_ASC = loadtxt(file_name,
                                              delimiter = ';')
 
    names_and_units = Property(depends_on = 'data_file')
    @cached_property
    def _get_names_and_units(self):
        '''names and units corresponding to the returned '_data_array' by 'loadtxt_bending'
        '''
        names = ['w_raw', 'eps_c_raw', 'F_raw']
        units = ['mm', '1*E-3', 'N']
        print 'names, units from .raw-file', names, units
        return names, units
    
    names_and_units_ASC = Property(depends_on = 'data_file')
    @cached_property
    def _get_names_and_units_ASC(self):
        ''' Extract the names and units of the measured data.
        The order of the names in the .DAT-file corresponds 
        to the order of the .ASC-file.   
        '''
        file_split = self.data_file.split('.')
        file_name = file_split[0] + '.DAT'
        data_file = open( file_name, 'r' )
        lines = data_file.read().split()
        names = []
        units = []
        for i in range( len( lines ) ):
            if lines[i] == '#BEGINCHANNELHEADER':
                name = lines[i + 1].split( ',' )[1]
                unit = lines[i + 3].split( ',' )[1]
                names.append( name )
                units.append( unit )

        print 'names, units extracted from .DAT-file', names, units
        return names, units

    factor_list_ASC = Property(depends_on = 'data_file')
    def _get_factor_list_ASC(self):
        return self.names_and_units_ASC[0]
    
    def _set_array_attribs( self ):
        '''Set the measured data as named attributes defining slices into 
        the processed data array.
        '''
        for i, factor in enumerate(self.factor_list):
            self.add_trait(factor, Array(value = self.processed_data_array[:, i], transient = True))

        if self.flag_ASC_file:
            for i, factor in enumerate(self.factor_list_ASC):
                self.add_trait(factor, Array(value = self.data_array_ASC[:, i], transient = True))
                


    elastomer_law = Property(depends_on='input_change')
    @cached_property
    def _get_elastomer_law(self):

        elastomer_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE', 'elastomer_f-w.raw')
        _data_array_elastomer = loadtxt_bending(elastomer_path)

        # force [kN]:
        # NOTE: after conversion 'F_elastomer' is a positive value
        #
        F_elastomer = -0.001 * _data_array_elastomer[:, 2].flatten()

        # displacement [mm]:
        # NOTE: after conversion 'w_elastomer' is a positive value
        #
        w_elastomer = -1.0 * _data_array_elastomer[:, 0].flatten()

        mfn_displacement_elastomer = MFnLineArray(xdata=F_elastomer, ydata=w_elastomer)
        return frompyfunc(mfn_displacement_elastomer.get_value, 1, 1)

    w_wo_elast = Property(depends_on='input_change')
    @cached_property
    def _get_w_wo_elast(self):
        # use the machine displacement for the center displacement:
        # subtract the deformation of the elastomer cushion between the cylinder
        # and change sign in positive values for vertical displacement [mm]
        #
        return self.w_raw - self.elastomer_law( self.F_raw )

    M_ASC = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_M_ASC(self):
        return self.F_ASC * self.length / 4.0
    
    M_raw = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_M_raw(self):
        return self.F_raw * self.length / 4.0

#    # get only the ascending branch of the response curve
#    #
#    max_force_idx = Property(Int)
#    def _get_max_force_idx(self):
#        '''get the index of the maximum force'''
#        return argmax(-self.Kraft)
#
#    f_asc = Property(Array)
#    def _get_f_asc(self):
#        '''get only the ascending branch of the response curve'''
#        return -self.Kraft[:self.max_force_idx + 1]
        

    def process_source_data(self):
        '''read in the measured data from file and assign
        attributes after array processing.
        '''
        super(ExpBT3PT, self).process_source_data()
        
        #---------------------------------------------
        # process data from .raw file (machine data)
        #---------------------------------------------
        
        # convert machine force [N] to [kN] and return only positive values
        #
        self.F_raw *= -0.001 

        # convert machine displacement [mm] to positive values
        # and remove offset
        #
        self.w_raw *= - 1.0 
        self.w_raw -= self.w_raw[0] 

        # convert [permille] to [-] and return only positive values
        #
        self.eps_c_raw *= - 0.001 
        
        # access the derived arrays to initiate their processing
        #
        self.w_wo_elast
        self.M_raw

        #---------------------------------------------
        # process data from .ASC file (displacement gauge)
        #---------------------------------------------

        # only if separate ASC.-file with force-displacement data from displacement gauge is available
        #
        if self.flag_ASC_file == True:

            self.F_ASC = - 1.0 * self.Kraft 

            # remove offset and change sign to return positive displacement values
            #
            if hasattr(self, "WA50"):
                self.WA50 *= -1
                self.WA50 -= self.WA50[0]
                WA50_avg = np.average( self.WA50 )
            
            if hasattr(self, "W10_u"):
                self.W10_u *= -1
                self.W10_u -= self.W10_u[0]
                W10_u_avg = np.average( self.W10_u )

            # check which displacement gauge has been used depending on weather two names are listed in .DAT file or only one
            # and assign values to 'w_ASC'
            #
            if hasattr(self, "W10_u") and hasattr(self, "WA50"):
                if W10_u_avg > WA50_avg: 
                    self.w_ASC = self.W10_u
                    print 'self.W10_u assigned to self.w_ASC'
                else:
                    self.w_ASC = self.WA50
                    print 'self.WA50 assigned to self.w_ASC'
            elif hasattr(self, "W10_u"):
                self.w_ASC = self.W10_u
                print 'self.W10_u assigned to self.w_ASC'
            elif hasattr(self, "WA50"):
                self.w_ASC = self.WA50
                print 'self.WA50 assigned to self.w_ASC'
               
            # convert strain from [permille] to [-], 
            # switch to positive values for compressive strains 
            # and remove offset
            #
            self.eps_c_ASC = - 0.001 * self.DMS_l 
            self.eps_c_ASC -= self.eps_c_ASC[0] 

            # access the derived arrays to initiate their processing
            #
            self.M_ASC



    #--------------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------------

    plot_templates = {'force / machine displacement (without w_elast)' : '_plot_force_machine_displacement_wo_elast',
                      'force / machine displacement (incl. w_elast)'   : '_plot_force_machine_displacement',
                      'force / gauge displacement'                     : '_plot_force_gauge_displacement',
                      
                      'smoothed force / gauge displacement'            : '_plot_smoothed_force_gauge_displacement',
                      'smoothed force / machine displacement'          : '_plot_smoothed_force_machine_displacement_wo_elast',

                      'moment / eps_c (ASC)'                           : '_plot_moment_eps_c_ASC',
                      'moment / eps_c (raw)'                           : '_plot_moment_eps_c_raw',

                      'smoothed moment / eps_c (ASC)'                  : '_plot_smoothed_moment_eps_c_ASC',
                      'smoothed moment / eps_c (raw)'                  : '_plot_smoothed_moment_eps_c_raw'
                     }

    default_plot_template = 'force / deflection (displacement gauge)'

    def _plot_force_machine_displacement_wo_elast(self, axes):

        xdata = self.w_wo_elast
        ydata = self.F_raw

        axes.plot( xdata, ydata )
#        xkey = 'deflection [mm]'
#        ykey = 'force [kN]'
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))


    def _plot_force_machine_displacement(self, axes):
        xdata = self.w_raw
        ydata = self.F_raw
        axes.plot( xdata, ydata )
#        xkey = 'deflection [mm]'
#        ykey = 'force [kN]'
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))


    def _plot_force_gauge_displacement(self, axes):
        xdata = self.w_ASC
        ydata = self.F_ASC
        axes.plot( xdata, ydata )
#        xkey = 'deflection [mm]'
#        ykey = 'force [kN]'
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))


    def _plot_smoothed_force_gauge_displacement(self, axes):

        # get the index of the maximum stress
        #
        max_force_idx = argmax(self.F_ASC)

        # get only the ascending branch of the response curve
        #
        F_asc = self.F_ASC[:max_force_idx + 1]
        w_asc = self.w_ASC[:max_force_idx + 1]

        n_points = int(self.n_fit_window_fraction * len(w_asc))
        F_smooth = smooth(F_asc, n_points, 'flat')
        w_smooth = smooth(w_asc, n_points, 'flat')

        axes.plot(w_smooth, F_smooth, color='blue', linewidth=2)


    def _plot_smoothed_force_machine_displacement_wo_elast(self, axes):

        # get the index of the maximum stress
        #
        max_force_idx = argmax(self.F_raw)

        # get only the ascending branch of the response curve
        #
        F_asc = self.F_raw[:max_force_idx + 1]
        w_asc = self.w_wo_elast[:max_force_idx + 1]

        n_points = int(self.n_fit_window_fraction * len(w_asc))
        F_smooth = smooth(F_asc, n_points, 'flat')
        w_smooth = smooth(w_asc, n_points, 'flat')

        axes.plot(w_smooth, F_smooth, color='blue', linewidth=2)

#        secant_stiffness_w10 = ( f_smooth[10] - f_smooth[0] ) / ( w_smooth[10] - w_smooth[0] )
#        w0_lin = array( [0.0, w_smooth[10] ], dtype = 'float_' )
#        f0_lin = array( [0.0, w_smooth[10] * secant_stiffness_w10 ], dtype = 'float_' )

        #axes.plot( w0_lin, f0_lin, color = 'black' )


    def _plot_moment_eps_c_ASC(self, axes):
        xkey = 'compressive strain [1*E-3]'
        ykey = 'moment [kNm]'
        xdata = self.eps_c_ASC
        ydata = self.M_ASC
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata)


    def _plot_moment_eps_c_raw(self, axes):
        xkey = 'compressive strain [1*E-3]'
        ykey = 'moment [kNm]'
        xdata = self.eps_c_raw
        ydata = self.M_raw
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata)


    n_fit_window_fraction = Float(0.1)

    smoothed_M_eps_c_ASC = Property(depends_on='input_change')
    @cached_property
    def _get_smoothed_M_eps_c_ASC(self):
        # get the index of the maximum stress
        max_idx = argmax(self.M_ASC)
        # get only the ascending branch of the response curve
        m_asc = self.M_ASC[:max_idx + 1]
        eps_c_asc = self.eps_c_ASC[:max_idx + 1]
        n_points = int(self.n_fit_window_fraction * len(eps_c_asc))
        m_smoothed = smooth(m_asc, n_points, 'flat')
        eps_c_smoothed = smooth(eps_c_asc, n_points, 'flat')
        return m_smoothed, eps_c_smoothed

    smoothed_eps_c_ASC = Property
    def _get_smoothed_eps_c_ASC(self):
        return self.smoothed_M_eps_c_ASC[1]

    smoothed_M_ASC = Property
    def _get_smoothed_M_ASC(self):
        return self.smoothed_M_eps_c_ASC[0]

    def _plot_smoothed_moment_eps_c_ASC(self, axes):
        axes.plot(self.smoothed_eps_c_ASC, self.smoothed_M_ASC, color='blue', linewidth=2)

    smoothed_M_eps_c_raw = Property(depends_on='input_change')
    @cached_property
    def _get_smoothed_M_eps_c_raw(self):
        # get the index of the maximum stress
        max_idx = argmax(self.M_raw)
        # get only the ascending branch of the response curve
        m_asc = self.M_raw[:max_idx + 1]
        eps_c_asc = self.eps_c_raw[:max_idx + 1]
        n_points = int(self.n_fit_window_fraction * len(eps_c_asc))
        m_smoothed = smooth(m_asc, n_points, 'flat')
        eps_c_smoothed = smooth(eps_c_asc, n_points, 'flat')
        return m_smoothed, eps_c_smoothed

    smoothed_eps_c_raw = Property
    def _get_smoothed_eps_c_raw(self):
        return self.smoothed_M_eps_c_raw[1]

    smoothed_M_raw = Property
    def _get_smoothed_M_raw(self):
        return self.smoothed_M_eps_c_raw[0]

    def _plot_smoothed_moment_eps_c_raw(self, axes):
        axes.plot(self.smoothed_eps_c_raw, self.smoothed_M_raw, color='blue', linewidth=2)

    #--------------------------------------------------------------------------------
    # view
    #--------------------------------------------------------------------------------

    traits_view = View(VGroup(
                         Group(
                              Item('length', format_str="%.3f"),
                              Item('width', format_str="%.3f"),
                              Item('thickness', format_str="%.3f"),
                              label='geometry'
                              ),
                         Group(
                              Item('loading_rate'),
                              Item('age'),
                              label='loading rate and age'
                              ),
                         Group(
                              Item('E_c', show_label=True, style='readonly', format_str="%.0f"),
                              Item('ccs@', show_label=False),
                              label='composite cross section'
                              )
                         ),
                        scrollable=True,
                        resizable=True,
                        height=0.8,
                        width=0.6
                        )

if __name__ == '__main__':

    from matresdev.db.exdb.ex_run_table import ExRunClassExt
    ex = ExRunClassExt(klass=ExpBT3PT)
    ex.configure_traits()
