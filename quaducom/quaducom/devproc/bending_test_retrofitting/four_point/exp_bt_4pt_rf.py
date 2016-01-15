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

# # overload the 'get_label' method from 'Item' to display units in the label
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

from numpy import array, fabs, where, copy, ones, argsort, \
    loadtxt, argmax, polyfit, poly1d, frompyfunc, dot, unique, around

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

from matresdev.db.exdb.loadtxt_novalue import loadtxt_novalue

# Access to the toplevel directory of the database
#
simdb = SimDB()


class ExpBT4PTRF(ExType):
    '''Experiment: Bending Test Four Point with RetroFitting beam
    '''
#    label = Str('four point bending test')

    implements(IExType)

    file_ext = 'DAT'

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

    # effective length of the bending test specimen
    # (does not include the part at each side of the specimens that leaps over the support lines)
    #
    length = Float(3.30, unit='m', input=True, table_field=True,
                           auto_set=False, enter_set=True)
    length_loadintroduction = Float(0.70, unit='m', input=True, table_field=True,
                           auto_set=False, enter_set=True)
    width = Float(0.50, unit='m', input=True, table_field=True,
                           auto_set=False, enter_set=True)
    thickness = Float(0.20, unit='m', input=True, table_field=True,
                           auto_set=False, enter_set=True)

    # age of the concrete at the time of testing
    age = Int(28, unit='d', input=True, table_field=True,
                             auto_set=False, enter_set=True)
    loading_rate = Float(1.0, unit='mm/min', input=True, table_field=True,
                            auto_set=False, enter_set=True)
    gauge_length_horizontal = Float(0.40, unit='m', input=True, table_field=True,
                            auto_set=False, enter_set=True)

    #--------------------------------------------------------------------------
    # composite cross section
    #--------------------------------------------------------------------------

    ccs = Instance(CompositeCrossSection)
    def _ccs_default(self):
        '''default settings'
        '''
        # SFB 532 - demonstrator textil and concrete:
        fabric_layout_key = '2D-05-11'
        concrete_mixture_key = 'FIL-10-09'
        orientation_fn_key = 'all0'
        n_layers = 2
        s_tex_z = 0.015 / (n_layers + 1)
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

    # put this into the ironing procedure processor
    #
    jump_rtol = Float(0.9,
                      auto_set=False, enter_set=True,
                      ironing_param=True)

    data_array_ironed = Property(Array(float),
                                  depends_on='data_array, +ironing_param, +axis_selection')
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
#            print 'self.names_and_units[0]',self.names_and_units[0]
#            print 'self.names_and_units',self.names_and_units
            if self.names_and_units[0][ idx ] in {'WA_M1', 'WA_M2', 'WA_L', 'WA_R'}:

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
                print 'force', unique(around(-self.data_array[jump_idx, 1], 2))
                # glue the curve at each jump together
                for jidx in jump_idx:
                    # get the offsets at each jump of the curve
                    shift = data_arr[jidx + 1] - data_arr[jidx]
                    # shift all succeeding values by the calculated offset
                    data_arr[jidx + 1:] -= shift

                data_array_ironed[:, idx] = data_arr[:]

        return data_array_ironed

    def process_source_data(self):
        '''read in the measured data from file and assign
        attributes after array processing.
        If necessary modify the assigned data, i.e. change
        the sign or specify an offset for the specific test setup.
        '''
        print '*** process source data ***'

        super(ExpBT4PTRF, self).process_source_data()

        self._read_data_array()

        # curve ironing:
        #
        self.processed_data_array = self.data_array_ironed

        # set attributes:
        #
        self._set_array_attribs()

        # PEEKEL-measuring software:

        # convert units and change signs
        self.Kraft -= self.Kraft[0]
        self.Kraft *= -1

        # (reset displacement gauges by their initial values and change sign
        # in order to return a positive value for a displacement)

        # vertical displacements [mm]:
        self.WA_M1 -= self.WA_M1[0]
        self.WA_M1 *= -1
        self.WA_M2 -= self.WA_M2[0]
        self.WA_M2 *= -1
        self.WA_L -= self.WA_L[0]
        self.WA_L *= -1
        self.WA_R -= self.WA_R[0]
        self.WA_R *= -1

        # horizontal displacements at the bottom side of the bending specimen [mm]
        self.WA_HR -= self.WA_HR[0]
        self.WA_HR *= -1
        self.WA_HM -= self.WA_HM[0]
        self.WA_HM *= -1
        self.WA_HL -= self.WA_HL[0]
        self.WA_HL *= -1

        # optional additional displacement gauges for crack width measuring of shear cracks
        self.Schub1 -= self.Schub1[0]
        self.Schub1 *= -1
        self.Schub1 -= self.Schub1[0]
        self.Schub1 *= -1

        # DMS for steel bars on vertical shear reinforcement
        self.VL2 -= self.VL2[0]
        self.VL2 *= -1
        self.VL3 -= self.VL3[0]
        self.VL3 *= -1
        self.VR2 -= self.VR2[0]
        self.VR2 *= -1
        self.VR3 -= self.VR3[0]
        self.VR3 *= -1

        # compressive strain at the upper side of the bending specimen at midspan [mm]
        # change unite from [nm/m], i.e. [10^(-6)*m / m], to [mm/m=permille]
        self.CM -= self.CM[0]
        self.CM /= 1000.
        self.CML -= self.CML[0]
        self.CML /= 1000.
        self.CMR -= self.CMR[0]
        self.CMR /= 1000.

        # DMS for longitudinal steel bars
        self.SL -= self.SL[0]
        self.SL /= 1000.
        self.SML -= self.SML[0]
        self.SML /= 1000.
        self.SR -= self.SR[0]
        self.SR /= 1000.
        self.SMR -= self.SMR[0]
        self.SMR /= 1000.
        self.SM1 -= self.SM1[0]
        self.SM1 /= 1000.
        self.SM2 -= self.SM2[0]
        self.SM2 /= 1000.

        # set attributes of displacement (original data before ironing):
        #
        WA_M1_orig = np.copy(self.WA_M1)
        self.add_trait("WA_M1_orig", Array(value=WA_M1_orig, transient=True))
        WA_M2_orig = np.copy(self.WA_M2)
        self.add_trait("WA_M2_orig", Array(value=WA_M2_orig, transient=True))
        WA_R_orig = np.copy(self.WA_R)
        self.add_trait("WA_R_orig", Array(value=WA_R_orig, transient=True))
        WA_L_orig = np.copy(self.WA_L)
        self.add_trait("WA_L_orig", Array(value=WA_L_orig, transient=True))

    K_bending_elast_c = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_K_bending_elast_c(self):
        '''calculate the analytical bending stiffness of the beam (4 point bending)
        relation between center deflection and 2 * load/2 in the thirdpoints (sum up to F)
        '''
        # @todo: update formula to be calculated based on 'with',total 'length' and 'lenght_loadintroduction'
        t = self.thickness
        w = self.width
        L = self.length
        L_load = self.length_loadintroduction

        # coposite E-modulus
        #
        E_c = self.E_c

        # moment of inertia
        #
        I_yy = t ** 3 * w / 12.

        delta_11 = (L ** 3) / 56.348 / E_c / I_yy

        # [MN/m]=[kN/mm] bending stiffness with respect to a force applied at center of the beam
        #
        K_bending_elast_c = 1 / delta_11
#         print 'K_bending_elast_c', K_bending_elast_c

        print 'K_bending_elast_c', K_bending_elast_c
        return K_bending_elast_c

    #--------------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------------

    plot_templates = {'force / deflection (center)'            : '_plot_force_deflection_center',
                      'force / deflection (center) - original' : '_plot_force_deflection_center_orig',
                      'force / deflection (thirdpoints)'       : '_plot_force_deflection_thirdpoints',
                      'strain (top/bottom) / force'            : '_plot_strain_top_bottom_force',
                      'displacement (ironed/original - center)': '_plot_ironed_orig_force_deflection_center',
                      'displacement (ironed/original - left)'  : '_plot_ironed_orig_force_deflection_left',
                      'displacement (ironed/original - right)' : '_plot_ironed_orig_force_deflection_right',
                     }

    default_plot_template = 'force / deflection (center)'

    # get only the ascending branch of the response curve
    #
    max_force_idx = Property(Int)
    def _get_max_force_idx(self):
        '''get the index of the maximum force'''
        # NOTE: processed data returns positive values for force and displacement
        return argmax(self.Kraft)

    def _plot_force_deflection_center(self, axes, offset_w=0., color='black', linewidth=1., label=None):
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:self.max_force_idx + 1]
        w_asc_M1 = self.WA_M1[:self.max_force_idx + 1]
        w_asc_M2 = self.WA_M2[:self.max_force_idx + 1]
        w_asc_Mavg = (w_asc_M1 + w_asc_M2) / 2.

        # add curves
        #
        axes.plot(w_asc_Mavg, f_asc, linewidth=linewidth, label=label, color=color)

        # add axes labels
        #
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))


    def _plot_force_deflection_center_orig(self, axes, offset_w=0., color='black', linewidth=1., label=None):
        '''plot the original data before jumps has been processed out
        '''
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:self.max_force_idx + 1]
        w_asc_M1_orig = self.WA_M1_orig[:self.max_force_idx + 1]
        w_asc_M2_orig = self.WA_M2_orig[:self.max_force_idx + 1]
        w_asc_Mavg_orig = (w_asc_M1_orig + w_asc_M2_orig) / 2.

        # add curves
        #
        axes.plot(w_asc_Mavg_orig, f_asc, linewidth=linewidth, label=label, color=color)

        # add axes labels
        #
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))


    def _plot_force_deflection_thirdpoints(self, axes):
        '''deflection at the third points (under the loading points)
        '''
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:self.max_force_idx + 1]
        # displacement left
        w_l_asc = self.WA_L[:self.max_force_idx + 1]
        # displacement right
        w_r_asc = self.WA_R[:self.max_force_idx + 1]

        axes.plot(w_l_asc, f_asc, color='green', linewidth=1)
        axes.plot(w_r_asc, f_asc, color='green', linewidth=1)


    def _plot_strain_top_bottom_force(self, axes, color='black', linewidth=1., label=None):
        '''plot compressive strains at top and
        average tensile strains based on horizontal displacement gauges
        '''
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:self.max_force_idx + 1]

        # compressive strains (top) [permile]
        eps_CM = self.CM [:self.max_force_idx + 1]
        eps_CMR = self.CMR [:self.max_force_idx + 1]
        eps_CML = self.CML [:self.max_force_idx + 1]

        # tensile strain (bottom) [permile];
        eps_HL = self.WA_HL[:self.max_force_idx + 1] / self.gauge_length_horizontal
        eps_HM = self.WA_HM[:self.max_force_idx + 1] / self.gauge_length_horizontal
        eps_HR = self.WA_HR[:self.max_force_idx + 1] / self.gauge_length_horizontal

        # add curves
        #
        axes.plot(eps_CM, f_asc, linewidth=linewidth, label='compression: eps_CM', color='grey')
        axes.plot(eps_CMR, f_asc, linewidth=linewidth, label='compression: eps_CMR', color='grey')
        axes.plot(eps_CML, f_asc, linewidth=linewidth, label='compression: eps_CML', color='grey')
        axes.plot(eps_HL, f_asc, linewidth=linewidth, label='tension: eps_HL', color='k')
        axes.plot(eps_HM, f_asc, linewidth=linewidth, label='tension: eps_HM', color='k')
        axes.plot(eps_HR, f_asc, linewidth=linewidth, label='tension: eps_HR', color='k')

        # add axes labels
        #
        xkey = 'strain [1*e-3]'
        ykey = 'force [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))


    def _plot_ironed_orig_force_deflection_center(self, axes):
        '''plot original displacement (center) as measured by the displacement gauge
        and compare with curve after data has been processed by ironing procedure
        '''
        # get only the ascending branch of the response curve
        F_asc = self.Kraft[:self.max_force_idx + 1]
        w_ironed_asc = self.WA_M1[:self.max_force_idx + 1]
        w_orig_asc = self.WA_M1_orig[:self.max_force_idx + 1]

        # add curves
        #
        axes.plot(w_ironed_asc, F_asc, color='blue', linewidth=1.5)
        axes.plot(w_orig_asc, F_asc, color='grey', linewidth=1.5)

        # add axes labels
        #
        xkey = 'deflection (original data / ironed data) [mm]'
        ykey = 'force [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))


    def _plot_ironed_orig_force_deflection_left(self, axes):
        '''plot original displacement (left) as measured by the displacement gauge
        and compare with curve after data has been processed by ironing procedure
        '''
        w_ironed = self.WA_L
        w_orig = self.WA_L_orig
        F = self.Kraft
        axes.plot(w_ironed, F)
        axes.plot(w_orig, F)
        xkey = 'deflection (original data / ironed data) [mm]'
        ykey = 'force [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))


    def _plot_ironed_orig_force_deflection_right(self, axes):
        '''plot original displacement (left) as measured by the displacement gauge
        and compare with curve after data has been processed by ironing procedure
        '''
        w_ironed = self.WA_R
        w_orig = self.WA_R_orig
        F = self.Kraft
        axes.plot(w_ironed, F)
        axes.plot(w_orig, F)
        xkey = 'deflection (original data / ironed data) [mm]'
        ykey = 'force [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))

    #--------------------------------------------------------------------------------
    # view
    #--------------------------------------------------------------------------------

    traits_view = View(VGroup(
                         Group(
                              Item('length', format_str="%.3f"),
                              Item('length_loadintroduction', format_str="%.3f"),
                              Item('width', format_str="%.3f"),
                              Item('thickness', format_str="%.3f"),
                              label='geometry'
                              ),
                         Group(
                              Item('loading_rate'),
                              Item('gauge_length_horizontal'),
                              Item('age'),
                              label='loading rate and age'
                              ),
                         Group(
                              Item('jump_rtol', format_str="%.4f"),
                              label='curve_ironing'
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
    ex = ExRunClassExt(klass=ExpBT4PTRF)
    ex.configure_traits()
