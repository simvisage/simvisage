# -------------------------------------------------------------------------
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
#

from traits.api import \
    Int, Str, \
    Array, Property, cached_property

from traitsui.api \
    import View, Item, HSplit, Group, VSplit

import numpy as np
from scipy.interpolate import interp1d

from aramis_cdt import AramisInfo, AramisFieldData, AramisCDT

from quaducom.devproc.tensile_test.dog_bone.exp_tt_db import ExpTTDB

from matresdev.db.exdb.ex_run_table import ExRunClassExt

DEPENDENCY_STR = 'data_file, start_t, delta_t, aramis_resolution_key'


def f_interp1d(x, xp, yp):
    F_t_interp = interp1d(xp, yp, bounds_error=False,
                          fill_value=0)
    return F_t_interp(x)


class ExpATTDB(ExpTTDB):

    '''Tensile test with aramis data
    '''

    # =========================================================================
    # 2D-ARAMIS PROCESSING
    # =========================================================================

    eps = Property(Array('float_'), output=True,
                   depends_on='input_change')
    '''Strain calculated using the displacement gauges
    placed along the edges of the specimen.
    '''
    @cached_property
    def _get_eps(self):

        W10_li = np.copy(self.W10_li)
        W10_re = np.copy(self.W10_re)

        # get the minimum value of the displacement gauges
        # used to reset the displacement gauges if they do not start at
        # zero
        min_W10_li = np.min(W10_li[:10])
        min_W10_re = np.min(W10_re[:10])

        # reset displacement gauges
        #
        W10_li -= min_W10_li
        W10_re -= min_W10_re

        # measured strains
        eps_li = W10_li / (self.gauge_length * 1000.)  # [mm/mm]
        eps_re = W10_re / (self.gauge_length * 1000.)

        # NOTE: if only 2 displacement gauges are used instead of 3
        # (only 'WA_re' for front and 'WA_li' for back)
        # below the average is performed
        # as = 0.5*( 0.5*(W10_re + W10_li) + W10_vo)
        #
        if np.average(eps_re) < 0.0001:
            print "displacement gauge 'WA_re' has not been used. Use value of \
                'WA_li' instead"
            eps_re = eps_li
        if np.average(eps_li) < 0.0001:
            print "displacement gauge 'WA_li' has not been used. Use value of\
                 'WA_re' instead"
            eps_li = eps_re

        # average strains
        #
        eps_m = (eps_li + eps_re) / 2.

        return eps_m

    aramis_resolution_key = Str('Xf15s3-Yf15s3')
    '''Specification of the resolution of the measured aramis field
    '''

    n_steps_aramis = Property(Int)
    '''
    '''

    def _get_n_steps_aramis(self):
        return self.aramis_info.number_of_steps

    t_aramis = Property(Array('float'),
                        depends_on=DEPENDENCY_STR)

    @cached_property
    def _get_t_aramis(self):
        step_times = self.aramis_field_data.step_times + \
            self.aramis_start_offset
        print 'self.aramis_start_offset', self.aramis_start_offset
        print '-----------------------------t_max_ARAMIS', step_times[-1]
        return step_times

    t_aramis_asc = Property(Array('float'),
                            depends_on=DEPENDENCY_STR)

    @cached_property
    def _get_t_aramis_asc(self):
        t_max = self.time_asc[-1]
        print '-----------------------------t_max_ASCII', t_max
        # @todo: make this using the first occurrence
        # of the condition and cut the array using slice
        return self.t_aramis[:self.maxF_idx_aramis]

    n_steps_asc = Property

    @cached_property
    def _get_n_steps_asc(self):
        'number of time steps in aramis after limiting the steps with t_max'
        # x = self.n_steps_aramis - len(self.t_aramis_cut)
        # if x == 0:
        # x = 1
        # print 'x ', x
        # return self.aramis_info.number_of_steps - x
        print '-----------------------------n_steps_cut2', \
            len(self.t_aramis_asc)
        return len(self.t_aramis_asc)

    aramis_info = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_aramis_info(self):
        af = self.get_cached_aramis_file(self.aramis_resolution_key)
        if af is None:
            return None
        return AramisInfo(data_dir=af)

    aramis_field_data = Property(depends_on='data_file,aramis_resolution_key')
    '''Field data including strains and displacements.
    '''
    @cached_property
    def _get_aramis_field_data(self):
        t_fail = self.time_asc[-1]
        ad = AramisFieldData(aramis_info=self.aramis_info)
        current_step = (np.abs(ad.step_times - t_fail).argmin())
        # print 'ad.step_times - t_fail', ad.step_times - t_fail
        ad.current_step = current_step
        return ad

    aramis_cdt = Property(depends_on='data_file,aramis_resolution_key')
    '''Field data including strains and displacements.
    '''
    @cached_property
    def _get_aramis_cdt(self):
        ad = self.aramis_field_data
        crack_detection_step = ad.current_step
        print '---------------------CRACK DETECTION STEP', crack_detection_step
        return AramisCDT(aramis_info=self.aramis_info,
                         crack_detection_step=crack_detection_step,
                         aramis_data=ad)

    number_of_cracks_aramis = Property(
        depends_on='data_file,aramis_resolution_key')

    def _get_number_of_cracks_aramis(self):
        self.aramis_cdt.aramis_data.current_step = self.maxF_idx_aramis
        n_slices = 10
        x = []
        for i in range(n_slices):
            ad = self.aramis_cdt.aramis_data
            ad.top_j = np.floor(ad.j_max / 10 * i).astype(int)
            ad.bottom_j = np.floor(ad.j_max / 10 * (i + 1)).astype(int)
            cf = self.aramis_cdt.crack_filter_avg
            x.append(np.sum(cf))

#         cf = self.aramis_cdt.crack_filter
#         print cf.shape
#         x = np.sum(cf, axis=0)

        mu = np.mean(x)
        std = np.std(x, ddof=1)
        mn = np.min(x)
        mx = np.max(x)
        print '---------------', mu, std, mn, mx
        return mu, std, mn, mx

    F_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_F_t_aramis(self):
        '''force interpolated to the time steps of aramis'''
        return f_interp1d(self.t_aramis, self.Bezugskanal, self.Kraft)

    maxF_idx_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_maxF_idx_aramis(self):
        return np.argmax(self.F_t_aramis)

    F_t_aramis_asc = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_F_t_aramis_asc(self):
        '''force interpolated to the time steps of aramis'''
        return self.F_t_aramis[:self.maxF_idx_aramis]

    eps_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_eps_aramis(self):
        '''strain from aramis'''
        return self.aramis_cdt.control_strain_t

    eps_aramis_asc = Property(depends_on='data_file,aramis_resolution_key')

    Finit_aramis = Property(depends_on='data_file,aramis_resolution_key')

    def _get_Finit_aramis(self):
        return self.F_t_aramis[self.aramis_cdt.init_step_avg_lst]

    eps_init_aramis = Property(depends_on='data_file,aramis_resolution_key')

    def _get_eps_init_aramis(self):
        return self.eps_aramis[self.aramis_cdt.init_step_avg_lst]

    @cached_property
    def _get_eps_aramis_asc(self):
        '''strain from aramis'''
        return self.eps_aramis[:self.maxF_idx_aramis]

    def process_aramisCDT_data(self):
        self.aramis_cdt._run_t_fired()
        self.aramis_cdt._run_back_fired()

    def _plot_aramis_F_t_asc(self, ax, **kwds):
        ax.plot(self.t_aramis_asc, self.F_t_aramis_asc, **kwds)

    def _plot_aramis_t_strain_asc(self, ax, **kwds):
        ax.plot(self.t_aramis_asc, self.eps_aramis_asc, **kwds)

    def _plot_aramis_F_strain_asc(self, ax, **kwds):
        ax.plot(self.eps_aramis_asc, self.F_t_aramis_asc, **kwds)

    def _plot_aramis_Finit_strain(self, ax, **kwds):
        ax.plot(self.eps_init_aramis, self.Finit_aramis, **kwds)

    # ---------------------------------
    # view
    # ---------------------------------

    traits_view = View(VSplit(
        HSplit(Group(
            Item('width', format_str="%.3f"),
            Item('length', format_str="%.3f"),
            Item('aramis_start_offset', format_str="%.0f"),
            springy=True,
            label='geometry',
            id='matresdev.db.exdb.ex_composite_bending_tensile_test.geometry',
            dock='tab',
        ),
            Group(
            Item('loading_rate_N'),
            Item('loading_rate_F'),
            Item('age'),
            Item('production_date'),
            Item('testing_date'),
            springy=True,
            label='loading rate and age',
            id='matresdev.db.exdb.ex_composite_bending_tensile_test.loading',
            dock='tab',),
            id='matresdev.db.exdb.ex_composite_bending_tensile_test.xxx',
            dock='tab',
        ),
        Group(
            Item('ccs@', resizable=True, show_label=False),
            label='composite cross section'
        ),
        #     label = 'input variables',
        #     id = 'matresdev.db.exdb.ex_composite_tensile_test.vgroup.inputs',
        #     dock = 'tab',
        #     scrollable = True,
        #                               ),
        Group(
            Item('E_c', style='readonly', show_label=True, format_str="%.0f"),
            Item('N_max', style='readonly',
                 emphasized=True, format_str="%.2f"),
            Item('F_max', style='readonly',
                 emphasized=True, format_str="%.2f"),
            Item('w_max', style='readonly',
                 emphasized=True, format_str="%.2f"),
            Item('u_max', style='readonly',
                 emphasized=True, format_str="%.2f"),
            Item('M_max', style='readonly',
                 emphasized=True, format_str="%.3f"),
            Item('MN_max', style='readonly',
                 emphasized=True, format_str="%.3f"),
            Item('MF_max', style='readonly',
                 emphasized=True, format_str="%.3f"),
            Item('w_pred', style='readonly',
                 emphasized=True, format_str="%.3f"),
            label='output characteristics',
            id='matresdev.db.exdb.ex_composite_bending_tensile_test.\
                vgroup.outputs',
            dock='tab',
            scrollable=True,
        ),
        scrollable=True,
        id='matresdev.db.exdb.ex_composite_bending_tensile_test.vgroup',
        dock='tab',
    ),
        id='matresdev.db.exdb.ex_composite_bending_tensile_test',
        dock='tab',
        scrollable=True,
        resizable=True,
        height=0.8,
        width=0.5,
    )


ExpATTDB.db = ExRunClassExt(klass=ExpATTDB)

# --------------------------------------------------------------

if __name__ == '__main__':

    from matresdev.db.simdb import SimDB
    from matresdev.db.exdb import ExRunView
    simdb = SimDB()
    import os

    ex_path = os.path.join(simdb.exdata_dir,
                           'tensile_tests', 'buttstrap_clamping',
                           '2013-12-01_TTb-4c-2cm-0-TU_Aramis2d_RR',
                           'TTb-4c-2cm-0-TU-V1.DAT')

    doe_reader = ExRunView(data_file=ex_path)
    doe_reader.configure_traits()

    ExpATTDB.db.configure_traits()
    # to see all experiments in one picture
