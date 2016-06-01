#-------------------------------------------------------------------------
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

import csv
from numpy import array, fabs, where, copy, ones, argsort, \
    loadtxt, argmax, polyfit, poly1d, frompyfunc, dot, unique, around
from traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements, DelegatesTo, \
    Callable
from traitsui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor, EnumEditor, Handler, FileEditor, VSplit, Group, \
    HGroup, Spring

from mathkit.array.smoothing import smooth
from matresdev.db.exdb.ex_type import ExType
from matresdev.db.exdb.i_ex_type import IExType
from matresdev.db.exdb.loadtxt_bending import loadtxt_bending
from matresdev.db.exdb.loadtxt_novalue import loadtxt_novalue
from matresdev.db.matdb.trc.composite_cross_section import \
    CompositeCrossSection, plain_concrete
from matresdev.db.matdb.trc.concrete_mixture \
    import ConcreteMixture
from matresdev.db.matdb.trc.fabric_layout \
    import FabricLayOut
from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp
import numpy as np


class ExpBT4PT(ExType):
    '''Experiment: Bending Test Four Point
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

    #-------------------------------------------------------------------------
    # specify inputs:
    #-------------------------------------------------------------------------

    # effective length of the bending test specimen
    # (does not include the 5cm part at each side of the specimens that leaps over the support lines)
    #
    length = Float(1.75, unit='m', input=True, table_field=True,
                   auto_set=False, enter_set=True)
    width = Float(0.2, unit='m', input=True, table_field=True,
                  auto_set=False, enter_set=True)
    thickness = Float(0.055, unit='m', input=True, table_field=True,
                      auto_set=False, enter_set=True)

    # age of the concrete at the time of testing
    age = Int(28, unit='d', input=True, table_field=True,
              auto_set=False, enter_set=True)
    loading_rate = Float(1.0, unit='mm/min', input=True, table_field=True,
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
        n_layers = 1
#         s_tex_z = 0.060 / (n_layers + 1)
        ccs = CompositeCrossSection(
            fabric_layup_list=[
                #                             plain_concrete(s_tex_z * 0.5),
                FabricLayUp(
                    n_layers=n_layers,
                    orientation_fn_key=orientation_fn_key,
                    #                                    s_tex_z=s_tex_z,
                    s_tex_z=0.030,  # [m]
                    fabric_layout_key=fabric_layout_key
                ),
                plain_concrete(0.025),  # [m]
                #                             plain_concrete(s_tex_z * 0.5)
            ],
            concrete_mixture_key=concrete_mixture_key
        )
        return ccs

    #--------------------------------------------------------------------------
    # Get properties of the composite
    #--------------------------------------------------------------------------

    # E-modulus of the composite at the time of testing
    E_c = Property(
        Float, unit='MPa', depends_on='input_change', table_field=True)

    def _get_E_c(self):
        return self.ccs.get_E_c_time(self.age)

    # E-modulus of the composite after 28 days
    E_c28 = DelegatesTo('ccs', listenable=False)

    # reinforcement ration of the composite
    rho_c = DelegatesTo('ccs', listenable=False)

    #-------------------------------------------------------------------------
    # define processing
    #-------------------------------------------------------------------------

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
            if self.names_and_units[0][idx] not in {'Kraft', 'Bezugskanal', 'DMS_o',
                                                    # for PEEKEL software check
                                                    # the following names:
                                                    'Datum/Uhrzeit', 'Kraft', 'Weg', 'DMS'}:

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

#                 print 'number of jumps removed in data_arr_ironed for', self.names_and_units[0][ idx ], ': ', jump_idx.shape[0]
#                 print 'force', unique(around(-self.data_array[jump_idx, 1], 2))
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

        super(ExpBT4PT, self).process_source_data()

        self._read_data_array()

        # curve ironing:
        #
        self.processed_data_array = self.data_array_ironed

        # set attributes:
        #
        self._set_array_attribs()

        # DIADAGO-measuring software:
        #
        if hasattr(self, "Kraft") and hasattr(self, "DB_mi") \
                and hasattr(self, "DB_re") and hasattr(self, "DB_li") \
                and hasattr(self, "W10_u") and hasattr(self, "DMS_o"):

            # convert units and change signs
            #
            self.Kraft -= self.Kraft[0]
            self.Kraft *= -1
            # vertical displacement at midspan [mm]:
            # (reset displacement gauge by its initial value and change sign
            # in order to return a positive value for a displacement)
            self.DB_mi -= self.DB_mi[0]
            self.DB_mi *= -1
            # vertical displacements at one third of the span (displacement under loading point) [mm]:
            # (left)
            self.DB_li -= self.DB_li[0]
            self.DB_li *= -1
            # (right)
            self.DB_re -= self.DB_re[0]
            self.DB_re *= -1
            # horizontal displacements at the bottom side of the bending specimen [mm]
            # (measuring length l_0 = 0.45 m)
            self.W10_u -= self.W10_u[0]
            self.W10_u *= -1
            # compressive strain at the upper side of the bending specimen at
            # midspan [mm]
            self.DMS_o -= self.DMS_o[0]
            # change unite from [nm/m], i.e. [10^(-6)*m / m], to [mm]
            self.DMS_o /= 1000.

            # set attributes of displacement (original data before ironing):
            #
            DB_mi_orig = self.data_array[:, 2]
            DB_mi_orig -= DB_mi_orig[0]
            DB_mi_orig *= -1
            self.add_trait(
                "DB_mi_orig", Array(value=DB_mi_orig, transient=True))

            DB_li_orig = self.data_array[:, 4]
            DB_li_orig -= DB_li_orig[0]
            DB_li_orig *= -1
            self.add_trait(
                "DB_li_orig", Array(value=DB_li_orig, transient=True))

            DB_re_orig = self.data_array[:, 5]
            DB_re_orig -= DB_re_orig[0]
            DB_re_orig *= -1
            self.add_trait(
                "DB_re_orig", Array(value=DB_re_orig, transient=True))

        # PEEKEL-measuring software:
        #
        if hasattr(self, "Kraft") and hasattr(self, "Weg") \
                and hasattr(self, "WA_Links") and hasattr(self, "WA_Mitte_Links") \
                and hasattr(self, "WA_Mitte") and hasattr(self, "WA_Mitte_Rechts") \
                and hasattr(self, "WA_Rechts") and hasattr(self, "WA_unten") and hasattr(self, "DMS"):

            # convert units and change signs
            #
            self.Kraft -= self.Kraft[0]
            self.Kraft *= -1
            # vertical displacements [mm]:
            # (reset displacement gauges by their initial values and change sign
            # in order to return a positive value for a displacement)
            self.WA_Links -= self.WA_Links[0]
            self.WA_Links *= -1
            self.DB_li = copy(self.WA_Links)  # set alias
            self.WA_Mitte -= self.WA_Mitte[0]
            self.WA_Mitte *= -1
            self.DB_mi = copy(self.WA_Mitte)  # set alias
            self.WA_Rechts -= self.WA_Rechts[0]
            self.WA_Rechts *= -1
            self.DB_re = copy(self.WA_Rechts)  # set alias
            self.WA_Mitte_Links -= self.WA_Mitte_Links[0]
            self.WA_Mitte_Links *= -1
            self.WA_Mitte_Rechts -= self.WA_Mitte_Rechts[0]
            self.WA_Mitte_Rechts *= -1

            # horizontal displacements at the bottom side of the bending specimen [mm]
            # (measuring length l_0 = 0.30 m)
            self.WA_unten -= self.WA_unten[0]
            self.WA_unten *= -1
            # compressive strain at the upper side of the bending specimen at
            # midspan [mm]
            self.DMS -= self.DMS[0]
            # change unite from [nm/m], i.e. [10^(-6)*m / m], to [mm]
            self.DMS /= 1000.
            self.DMS_o = copy(self.DMS)  # set alias

            # set attributes of displacement (original data before ironing):
            #
            DB_mi_orig = np.copy(self.data_array[:, 5])
            DB_mi_orig -= DB_mi_orig[0]
            DB_mi_orig *= -1
            self.add_trait(
                "DB_mi_orig", Array(value=DB_mi_orig, transient=True))

            DB_li_orig = np.copy(self.data_array[:, 3])
            DB_li_orig -= DB_li_orig[0]
            DB_li_orig *= -1
            self.add_trait(
                "DB_li_orig", Array(value=DB_li_orig, transient=True))

            DB_re_orig = np.copy(self.data_array[:, 7])
            DB_re_orig -= DB_re_orig[0]
            DB_re_orig *= -1

            self.add_trait(
                "DB_re_orig", Array(value=DB_re_orig, transient=True))

    K_bending_elast_c = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_K_bending_elast_c(self):
        '''calculate the analytical bending stiffness of the beam (4 point bending)
        relation between center deflection and 2 * load/2 in the thirdpoints (sum up to F)
        '''
        t = self.thickness
        w = self.width
        L = self.length

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

    K_bending_elast_thirdpoints = Property(
        Array('float_'), depends_on='input_change')

    @cached_property
    def _get_K_bending_elast_thirdpoints(self):
        '''calculate the analytical bending stiffness of the beam (4 point bending)
        relation between thirdpoint deflection and 2 * load/2 in the thirdpoints (sum up to F)
        '''
        t = self.thickness
        w = self.width
        L = self.length

        # coposite E-modulus
        #
        E_c = self.E_c

        # moment of inertia
        #
        I_yy = t ** 3 * w / 12.

        delta_11 = (L ** 3) * 5 / 324 / E_c / I_yy

        # [MN/m]=[kN/mm] bending stiffness with respect to a force applied at center of the beam
        #
        K_bending_elast_thirdpoints = 1 / delta_11
#         print 'K_bending_elast', K_bending_elast

        print 'K_bending_elast_thirdpoints', K_bending_elast_thirdpoints
        return K_bending_elast_thirdpoints

    #-------------------------------------------------------------------------
    # plot templates
    #-------------------------------------------------------------------------

    plot_templates = {'force / deflection (center)': '_plot_force_deflection_center',
                      'force / deflection (center) - original': '_plot_force_deflection_center_orig',
                      'smoothed force / deflection (center)': '_plot_smoothed_force_deflection_center',
                      'force / deflection (thirdpoints)': '_plot_force_deflection_thirdpoints',
                      'strain (top/bottom) / force': '_plot_strain_top_bottom_force',
                      'displacement (ironed/original - center)': '_plot_ironed_orig_force_deflection_center',
                      'displacement (ironed/original - left)': '_plot_ironed_orig_force_deflection_left',
                      'displacement (ironed/original - right)': '_plot_ironed_orig_force_deflection_right',
                      'displacement (center, thirdpoints) (ironed, average)': '_plot_ironed_force_deflection_avg'
                      }

    default_plot_template = 'force / deflection (center)'

    # get only the ascending branch of the response curve
    #
    max_force_idx = Property(Int)

    def _get_max_force_idx(self):
        '''get the index of the maximum force'''
        # NOTE: processed data returns positive values for force and
        # displacement
        return argmax(self.Kraft)

    def _plot_force_deflection_center(self, axes, offset_w=0., color='black', linewidth=1., label=None):
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:self.max_force_idx + 1]
        w_asc = self.DB_mi[:self.max_force_idx + 1]

        # add curves
        #
        axes.plot(w_asc, f_asc, linewidth=linewidth, label=label, color=color)

        # add axes labels
        #
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))
        # draw linear stiffness for 2 mm range
        f_max = f_asc[-1]
        K_c = self.K_bending_elast_c
        w_linear = np.array([0., f_max / K_c])
        F_linear = np.array([0., f_max])
        axes.plot(w_linear, F_linear, linestyle='--',
                  color='black', linewidth=linewidth)

    def _plot_force_deflection_center_orig(self, axes, offset_w=0., color='black', linewidth=1., label=None):
        '''plot the original data before jumps has been processed out
        '''
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:self.max_force_idx + 1]
        w_asc = self.DB_mi_orig[:self.max_force_idx + 1]

        # add curves
        #
        axes.plot(w_asc, f_asc, linewidth=linewidth, label=label, color=color)

        # add axes labels
        #
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))
        # draw linear stiffness for 2 mm range
        f_max = f_asc[-1]
        K_c = self.K_bending_elast_c
        w_linear = np.array([0., f_max / K_c])
        F_linear = np.array([0., f_max])
        axes.plot(w_linear, F_linear, linestyle='--',
                  color='black', linewidth=linewidth)

    n_fit_window_fraction = Float(0.1)

    def _plot_smoothed_force_deflection_center(self, axes):
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:self.max_force_idx + 1]
        w_asc = self.DB_mi[:self.max_force_idx + 1]

        # add axes labels
        #
        n_points = int(self.n_fit_window_fraction * len(w_asc))
        f_smooth = smooth(f_asc, n_points, 'flat')
        w_smooth = smooth(w_asc, n_points, 'flat')

        # add curves
        #
        axes.plot(w_smooth, f_smooth, color='blue', linewidth=2)

    def _plot_force_deflection_thirdpoints(self, axes):
        '''deflection at the third points (under the loading points)
        '''
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:self.max_force_idx + 1]
        # displacement left
        w_l_asc = self.DB_li[:self.max_force_idx + 1]
        # displacement rigth
        w_r_asc = self.DB_re[:self.max_force_idx + 1]

#        # average
#        w_lr_asc = (w_l_asc + w_r_asc) / 2

#        axes.plot( w_lr_asc, f_asc, color = 'green', linewidth = 2 )
        axes.plot(w_l_asc, f_asc, color='green', linewidth=1)
        axes.plot(w_r_asc, f_asc, color='green', linewidth=1)
        w_linear = 2 * np.array([0., 1.])
        F_linear = 2 * np.array([0., self.K_bending_elast_thirdpoints])
        axes.plot(w_linear, F_linear, linestyle='--')

    def _plot_strain_top_bottom_force(self, axes):
        '''deflection at the third points (under the loading points)
        '''
        # get only the ascending branch of the response curve
        f_asc = self.Kraft[:self.max_force_idx + 1]
        # compressive strain (top) [permile]
        eps_c = self.DMS_o[:self.max_force_idx + 1]
        # tensile strain (bottom) [permile];

        # NOTE: measuring length is not specified in 'exp_bt'-setup!
        # @todo: add this as configurable trait in the 'exp_bt' setup!
        if hasattr(self, "W10_u"):
            # NOTE: only valid for constant measuring length l_0 = 0.45m
            eps_t = self.W10_u[:self.max_force_idx + 1] / 0.45

        if hasattr(self, "WA_unten"):
            # NOTE: only valid for constant measuring length l_0 = 0.30m
            eps_t = self.WA_unten[:self.max_force_idx + 1] / 0.30

        # add curves
        #
        axes.plot(eps_c, f_asc, color='blue', linewidth=1)
        axes.plot(eps_t, f_asc, color='red', linewidth=1)

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
        w_ironed_asc = self.DB_mi[:self.max_force_idx + 1]
        w_orig_asc = self.DB_mi_orig[:self.max_force_idx + 1]

        # add curves
        #
        axes.plot(w_ironed_asc, F_asc, color='blue', linewidth=1.5)
        axes.plot(w_orig_asc, F_asc, color='grey', linewidth=1.5)

        # add axes labels
        #
        xkey = 'deflection (original data / ironed data) [mm]'
        ykey = 'force [kN]'
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))

#        fw_arr = np.hstack([F_asc[:, None], w_ironed_asc[:, None]])
#        print 'fw_arr.shape', fw_arr.shape
#        np.savetxt('BT-4PT-12c-6cm-TU-SH4-V1_f-w_interpolated.csv', fw_arr, delimiter='    ')

    def _plot_ironed_orig_force_deflection_left(self, axes):
        '''plot original displacement (left) as measured by the displacement gauge
        and compare with curve after data has been processed by ironing procedure
        '''
        w_ironed = self.DB_li
        w_orig = self.DB_li_orig
        F = self.Kraft
        axes.plot(w_ironed, F)
        axes.plot(w_orig, F)
        xkey = 'deflection (original data / ironed data) [mm]'
        ykey = 'force [kN]'
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))

    def _plot_ironed_orig_force_deflection_right(self, axes):
        '''plot original displacement (left) as measured by the displacement gauge
        and compare with curve after data has been processed by ironing procedure
        '''
        w_ironed = self.DB_re
        w_orig = self.DB_re_orig
        F = self.Kraft
        axes.plot(w_ironed, F)
        axes.plot(w_orig, F)
        xkey = 'deflection (original data / ironed data) [mm]'
        ykey = 'force [kN]'
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))

    def _plot_ironed_force_deflection_avg(self, axes):
        '''plot center displacement and the average of the thirdpoints displacement measured by the displacement gauge
        after data has been processed by ironing procedure (remove resetting jumps of the displacement gauges)
        '''
        # get only the ascending branch of the response curve
        F_asc = self.Kraft[:self.max_force_idx + 1]
        w_c_ironed_asc = self.DB_mi[:self.max_force_idx + 1]
        w_l_ironed_asc = self.DB_li[:self.max_force_idx + 1]
        w_r_ironed_asc = self.DB_re[:self.max_force_idx + 1]
        w_lr_avg = (w_l_ironed_asc + w_r_ironed_asc) / 2.
        axes.plot(w_lr_avg, F_asc, color='black', linewidth=1.5, linestyle='-')
        axes.plot(
            w_c_ironed_asc, F_asc, color='black', linewidth=1.5, linestyle='-')
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))

    #-------------------------------------------------------------------------
    # view
    #-------------------------------------------------------------------------

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
    ex = ExRunClassExt(klass=ExpBT4PT)
    ex.configure_traits()
