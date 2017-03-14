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
# Created on Feb 15, 2010 by: rch, ascholzen, laura-si

# @todo - construct the class for fabric layout calculating the
#         cs-area of the reinforcement.
#       - instead of processed array - construct the array traits accessible
#         with the name of the measured channels
#       - reread the pickle file without processing the data (take care to reestablish
#         the link from the ex_type to the ex_run
#       - define the exdb_browser showing the inputs and outputs in a survey
#       - define the ExTreatment class with cumulative evaluation of the response values.
#
#

from traits.api import \
    Int, Float, \
    on_trait_change, Instance, \
    Array, Property, cached_property, \
    Bool, Event, implements, \
    DelegatesTo, Str
from traitsui.api \
    import View, Item, HSplit, Group, VSplit
from aramis_cdt import AramisInfo, AramisFieldData, AramisCDT
from matresdev.db.exdb.ex_run_table import ExRunClassExt
from matresdev.db.exdb.ex_type import ExType
from matresdev.db.exdb.i_ex_type import IExType

from quaducom.devproc.bending_test.four_point.exp_bt_4pt import ExpBT4PT

from matresdev.db.matdb.trc.composite_cross_section \
    import CompositeCrossSection, plain_concrete
from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp
import numpy as np


# class ExpBTTDB(ExType):
class ExpBT4PTAramis2d(ExpBT4PT):

    '''Experiment: Bending Test with Aramis2d measurements
    '''

    #--------------------------------------------------------------------------
    # specify inputs:
    #--------------------------------------------------------------------------

    # start of 2D-Aramis
    start_time_aramis = Float(0, unit='m', input=True, table_field=True,
                              auto_set=False, enter_set=True)

    integ_radius = Int(3, unit='m', input=True, table_field=True,
                              auto_set=False, enter_set=True)

    integ_radius_crack = Int(5, unit='m', input=True, table_field=True,
                              auto_set=False, enter_set=True)

    ddd_ux_avg_threshold = Float(-0.5e-3, unit='m', input=True, table_field=True,
                              auto_set=False, enter_set=True)

    ddd_ux_threshold = Float(-0.5e-3, unit='m', input=True, table_field=True,
                              auto_set=False, enter_set=True)

    dd_ux_avg_threshold = Float(-0.5e-3, unit='m', input=True, table_field=True,
                              auto_set=False, enter_set=True)

    dd_ux_threshold = Float(-0.5e-3, unit='m', input=True, table_field=True,
                              auto_set=False, enter_set=True)

    scale_data_factor = Float(1.0, unit='', input=True, table_field=True,
                              auto_set=False, enter_set=True)

    crack_detaction_step = Int(-1, unit='', input=True, table_field=True,
                              auto_set=False, enter_set=True)

    crack_detaction_current_step = Bool(True, unit='', input=True, table_field=True,
                              auto_set=False, enter_set=True)

    #--------------------------------------------------------------------------
    # define processing
    #--------------------------------------------------------------------------

    def process_source_data(self):
        '''read in the measured data from file and assign
        attributes after array processing.
        If necessary modify the assigned data, e.i. change
        the sign or specify an offset for the specific test setup.
        '''
        super(ExpBT4PTAramis2d, self).process_source_data()

        self.F = self.Kraft
        self.w = self.Weg
        self.w -= self.w[0]
        self.w *= -1
        self.t = self.Bezugskanal

    #--------------------------------------------------------------------------
    # Get the maximum force index to cut off the descending part of the curves
    #--------------------------------------------------------------------------

    max_F_idx = Property(Int, depends_on='input_change')

    @cached_property
    def _get_max_F_idx(self):
        return np.argmax(self.F)

    #--------------------------------------------------------------------------
    # Get only the ascending branch of the response curve
    #--------------------------------------------------------------------------

    F_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_F_asc(self):
        return self.F[:self.max_F_idx + 1]

    t_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_t_asc(self):
        return self.t[:self.max_F_idx + 1]

    w_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_w_asc(self):
        return self.w[:self.max_F_idx + 1]


    #--------------------------------------------------------------------------
    # Get the moment arrays
    #--------------------------------------------------------------------------

    # distance between the support and the load introduction (F/2) of the four point bending test
    length_0 = 0.70  # [m] for 'BT_cyc'

    M_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_M_asc(self):
        'resulting moment'
        return self.F_asc / 2. * self.length_0

    #--------------------------------------------------------------------------
    # Get maximum values of the variables
    #--------------------------------------------------------------------------

    F_max = Property(Float, depends_on='input_change',
                     output=True, table_field=True, unit='kN')

    @cached_property
    def _get_F_max(self):
        return self.F_asc[-1]

    t_max = Property(Float, depends_on='input_change',
                     output=True, table_field=True, unit='kNm')

    @cached_property
    def _get_t_max(self):
        return self.t_asc[-1]

    w_max = Property(Float, depends_on='input_change',
                     output=True, table_field=True, unit='mm')

    @cached_property
    def _get_w_max(self):
        return self.w_asc[-1]

    M_max = Property(Float, depends_on='input_change',
                     output=True, table_field=True, unit='kNm')

    @cached_property
    def _get_M_max(self):
        return self.M_asc[-1]

    #--------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------

    plot_templates = {'w(t)': '_plot_w_t',
                      'F(w)': '_plot_F_w',
                      }

    default_plot_template = 'F(w)'

    def _plot_w_t(self, axes):
        '''Displacement versus time
        '''
        print 'self.w_asc', self.w_asc
        axes.plot(self.t_asc, self.w_asc, color='black')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('w [mm]')
        axes.set_ylim(0, 10)

    def _plot_F_w(self, axes):
        '''Bending force_cut versus displacement_cut
        '''
        axes.plot(self.w_asc, self.F_asc, color='blue')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('w [mm]')
        axes.set_ylabel('F [kN]')
        axes.set_ylim(0, 6)
        axes.set_xlim(0, 10)

    #=========================================================================
    # 2D-ARAMIS PROCESSING
    #=========================================================================

    aramis_resolution_key = Str('Xf15s3-Yf15s3')
    '''Specification of the resolution of the measured aramis field
    '''

    n_steps_aramis = Property(Int)

    def _get_n_steps_aramis(self):
        return self.aramis_info.number_of_steps

    t_aramis = Property(
        Array('float'), depends_on='data_file, start_t, delta_t, aramis_resolution_key')

    @cached_property
    def _get_t_aramis(self):
        step_times = self.aramis_field_data.step_times
        # step_times = self.aramis_field_data.step_times + self.start_time_aramis
        print 'self.start_time_aramis', self.start_time_aramis
        print '-----------------------------t_max_ARAMIS', step_times[-1]
        t_max = self.t[-1]
        print '-----------------------------t_max_ASCII', t_max
        # @todo: make this using the first occurence of the condition and cut the array using slice
        return step_times[np.where(step_times < t_max)]

    t_aramis_cut = Property(
        Array('float'), depends_on='data_file, start_t, delta_t, aramis_resolution_key')

    @cached_property
    def _get_t_aramis_cut(self):
        # print 't_aramis', self.t_aramis
        # print 't_aramis_cut', self.t_aramis[:]
        return self.t_aramis[:-1]

    n_steps = Property

    @cached_property
    def _get_n_steps(self):
        'number of time steps in aramis after limiting the steps with t_max'
        # x = self.n_steps_aramis - len(self.t_aramis_cut)
        # if x == 0:
        # x = 1
        # print 'x ', x
        # return self.aramis_info.number_of_steps - x
        print '-----------------------------n_steps_cut2', len(self.t_aramis_cut)
        return len(self.t_aramis_cut)

    aramis_info = Property(depends_on='data_file,aramis_resolution_key,input_change')

    @cached_property
    def _get_aramis_info(self):
        af = self.get_cached_aramis_file(self.aramis_resolution_key)
        if af == None:
            return None
        return AramisInfo(data_dir=af)

    aramis_field_data = Property(depends_on='data_file,aramis_resolution_key,scale_data_factor,input_change')
    '''Field data including strains and displacements.
    '''
    @cached_property
    def _get_aramis_field_data(self):
        print '_get_aramis_field_data'
        t_fail = self.t_asc[-1]
        ad = AramisFieldData(aramis_info=self.aramis_info,
                             integ_radius=self.integ_radius,
                             integ_radius_crack=self.integ_radius_crack,
                             transform_data=True,
                             scale_data_factor=self.scale_data_factor)
        current_step = (np.abs(ad.step_times - t_fail).argmin())
        # print 'ad.step_times - t_fail', ad.step_times - t_fail
        ad.current_step = current_step
        return ad

    aramis_cdt = Property(depends_on='data_file,aramis_resolution_key,aramis_field_data,input_change')
    '''Field data including strains and displacements.
    '''
    @cached_property
    def _get_aramis_cdt(self):
        print '_get_aramis_cdt'
        ad = self.aramis_field_data
        crack_detection_step = self.crack_detaction_step
        # if crack-detections-step is set to current step:
        if self.crack_detaction_current_step:
            crack_detection_step = ad.current_step
        print '---------------------CRACK DETECTION STEP', crack_detection_step
        return AramisCDT(aramis_info=self.aramis_info,
                         crack_detection_step=crack_detection_step,
                         aramis_data=ad,
                         ddd_ux_avg_threshold=self.ddd_ux_avg_threshold,
                         ddd_ux_threshold=self.ddd_ux_threshold,
                         dd_ux_avg_threshold=self.dd_ux_avg_threshold,
                         dd_ux_threshold=self.dd_ux_threshold
                         )

    F_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_F_t_aramis(self):
        'bending force interpolated to the time steps of aramis'
        return np.interp(self.t_aramis_cut, self.t, self.F)

    w_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_w_t_aramis(self):
        'displacement (with eliminated predeformation) interpolated to the time steps of aramis'
        return np.interp(self.t_aramis_cut, self.t_cut_asc, self.w)

    M_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_M_t_aramis(self):
        'resulting moment interpolated to the time steps of aramis'
        return self.F_t_aramis * self.length_0 / 2

    #-------------------------------------------------------------------------
    # crack bridge strain
    #-------------------------------------------------------------------------
    crack_filter_avg = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_crack_filter_avg(self):
        ''' method to get number and position of cracks
        '''
        ai = self.aramis_info
        if ai == None:
            return None
        field_data = self.aramis_field_data
        cdt = self.aramis_cdt
        # print field_data.x_arr_0[0, cdt.crack_filter_avg]
        return field_data.x_arr_0[0, cdt.crack_filter_avg]

    crack_bridge_strain_all = Property(
        depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_crack_bridge_strain_all(self):
        '''method to get crack bridge strain for the cracks determined by crack_filter_avg
        '''
        ai = self.aramis_info
        if ai == None:
            return None
        field_data = self.aramis_field_data

        # get the indices oft the cracks
        b = field_data.x_arr_0[0]
        c = self.crack_filter_avg
        crack_idx_list = []
        for c in c:
            crack_idx = np.where(b == c)
            crack_idx_list.append(crack_idx[0])

        crack_idx_arr1 = np.array(crack_idx_list)
        crack_idx_arr2 = np.rollaxis(crack_idx_arr1, 1, 0)
        crack_idx_arr = crack_idx_arr2[0]
        print 'crack_idx_arr', crack_idx_arr

        # get the indices of the middle between two cracks
        crack_mid_idx_arr = (crack_idx_arr[0:-1] + crack_idx_arr[1:]) / 2
        print 'crack_mid_idx_arr', crack_mid_idx_arr
        i_max = len(crack_mid_idx_arr)
        # print 'i_max', i_max

        if i_max <= 1:
            return None

        # ux = field_data.ux_arr
        # x_0 = field_data.x_arr_0

        # get crack bridge strain
        for step, t in enumerate(self.t_aramis_cut):
            field_data.current_step = step
            # eps_range = 1
            eps_list = []
            for i in range(0, i_max - 1, 1):
                # ux1 = np.mean(ux[:, crack_mid_idx_arr[i] - eps_range : crack_mid_idx_arr[i] + eps_range ], axis=1)
                # ux2 = np.mean(ux[:, crack_mid_idx_arr[i + 1] - eps_range: crack_mid_idx_arr[i + 1] + eps_range ], axis=1)
                # x_und1 = np.mean(x_und [:, crack_mid_idx_arr[i] - eps_range: crack_mid_idx_arr[i] + eps_range ], axis=1)
                # x_und2 = np.mean(x_und [:, crack_mid_idx_arr[i + 1] - eps_range: crack_mid_idx_arr[i + 1] + eps_range ], axis=1)

                # ux1 = ux[:, 0]
                # ux2 = ux[:, -1]
                # x_0_1 = x_0[:, 0]
                # x_0_2 = x_0[:, -1]

                # print 'ux1', ux1
                # eps_crack_bridge_i = (ux2 - ux1) / (x_0_2 - x_0_1)

                eps_crack_bridge_i = np.mean(
                    field_data.d_ux[:, crack_mid_idx_arr[i]:crack_mid_idx_arr[i + 1]], axis=1)
                eps_list.append(eps_crack_bridge_i)
            # print 'crack_bridge_strain_all', np.array(eps_list, dtype='f')
        return np.array(eps_list, dtype='f')

    idx_failure_crack = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_idx_failure_crack(self):
        '''method to get the index of the failure crack and the index of the corresponding
           middle to the neighbor cracks
        '''

        ai = self.aramis_info
        if ai == None:
            return None
        field_data = self.aramis_field_data

        # identify failure crack
        a = self.crack_bridge_strain_all  # only strain in the last step
        # print 'a', a
        if a == None:
            return None

        max_strain = a.max()
        # print 'max_strain', max_strain
        max_idx = np.where(a == max_strain)[0]
        # print 'max_idx', max_idx

        # get the indices oft all cracks
        b = field_data.x_arr_0[0]
        c = self.crack_filter_avg
        crack_idx_list = []
        for c in c:
            crack_idx = np.where(b == c)
            crack_idx_list.append(crack_idx[0])

        crack_idx_arr1 = np.array(crack_idx_list)
        crack_idx_arr2 = np.rollaxis(crack_idx_arr1, 1, 0)
        crack_idx_arr = crack_idx_arr2[0]
        # print 'crack_idx_arr', crack_idx_arr2[0]

        # get the indices of the middle between two cracks
        crack_mid_idx_arr = (crack_idx_arr[0:-1] + crack_idx_arr[1:]) / 2
        # print 'crack_mid_idx_arr', crack_mid_idx_arr

        idx_failure_crack = crack_idx_arr[max_idx + 1]
        print 'idx_failure_crack', idx_failure_crack
        print 'shape_ad.x_arr_0', np.shape(field_data.x_arr_0)
        print 'number_failure_crack', max_idx + 2
        idx_border1 = crack_mid_idx_arr[max_idx]
        print 'idx_border1 ', idx_border1
        idx_border2 = crack_mid_idx_arr[max_idx + 1]
        print 'idx_border2 ', idx_border2

        return idx_failure_crack, idx_border1, idx_border2

    #-------------------------------------------------------------------------
    # get max tensile strain in first reinforcement layer
    #-------------------------------------------------------------------------

    h_re1_threshold = Float(15., auto_set=False, enter_set=True)
    '''Threshold for position of first reinforcement layer (1 layers); c_nom = 15mm.
    '''

#     h_re1_6_threshold = Float(2.86, auto_set=False, enter_set=True)
#     '''Threshold for position of first reinforcement layer (6 layers).
#     '''
#
#     h_re1_4_threshold = Float(4.0, auto_set=False, enter_set=True)
#     '''Threshold for position of first reinforcement layer (4 layers).
#     '''

    meas_field = Property(depends_on='data_file,aramis_resolution_key')
    '''Get length and height of measuring field (calculated from the center of the facets)
    '''
    @cached_property
    def _get_meas_field(self):
        ai = self.aramis_info
        if ai == None:
            return None

        field_data = self.aramis_field_data
        field_data.current_step = self.n_steps

        x = field_data.x_arr_0
        y = field_data.y_arr_0

        l_mf1 = x[0, -1] - x[0, 0]
        l_mf2 = x[-1, -1] - x[-1, 0]

        h_mf1 = y[0, 0] - y[-1, 0]
        h_mf2 = y[0, -1] - y[-1, -1]

        l_mf = (l_mf1 + l_mf2) / 2
        h_mf = (h_mf1 + h_mf2) / 2
        # print 'meas_field', l_mf, h_mf
        return l_mf, h_mf

    h_dis = Property(depends_on='data_file,aramis_resolution_key')
    '''Get the distance between specimen edge and first / last facet node in y-direction.
    '''
    @cached_property
    def _get_h_dis(self):
        ai = self.aramis_info
        if ai == None:
            return None
        h_dis = (self.thickness * 1000. - self.meas_field[1]) / 2
        # print 'h_dis', h_dis
        return h_dis

    pos_fa = Property(depends_on='data_file,aramis_resolution_key')
    '''Get position of first and last facet-node in y-direction
    '''
    @cached_property
    def _get_pos_fa(self):
        ai = self.aramis_info
        if ai == None:
            return None

        pos_no_f = self.thickness * 1000. - self.h_dis
        # position of last node in y-direction
        pos_no_l = self.h_dis

        return pos_no_f, pos_no_l

    eps1_t_aramis = Property(depends_on='data_file,aramis_resolution_key')
    '''Tensile strain in first reinforcement layer
    '''
    @cached_property
    def _get_eps1_t_aramis(self):
        ai = self.aramis_info
        if ai == None:
            return None

        field_data = self.aramis_field_data

        # get number of facet nodes in y-direction
        n_fa = field_data.d_ux.shape[0]
        # print 'n_fa', n_fa

        # distance between facet nodes in y direction
        dis_fa = self.meas_field[1] / (n_fa - 1)
        # print 'dis_fa', dis_fa

        # get distance form top edge of mask to first reinforcement layer
        # for 55mm bending specimens
        pos_re1_6 = self.h_re1_threshold - self.h_dis
#         # get distance form top edge of mask to first reinforcement layer
#         pos_re1_6 = self.h_re1_6_threshold - self.h_dis
#         # print 'pos_re1_6', pos_re1_6  # -> negative value
#         pos_re1_4 = self.h_re1_4_threshold - self.h_dis
#         # print 'pos_re1_4', pos_re1_4

        # indices of strain next to position of reinforcement layer
        if pos_re1_6 < 0:
            idx_6a = 0
        else:
            idx_6a = pos_re1_6 / dis_fa
            idx_6a = round(idx_6a)
        idx_6b = idx_6a + 1

        a = self.crack_bridge_strain_all
        eps_t_list = []

        # get eps array
        for step, t in enumerate(self.t_aramis_cut):
            field_data.current_step = step
            if a == None:
                mid_idx = field_data.d_ux.shape[1] / 2
                eps_range = 3
                eps = np.mean(
                    field_data.d_ux[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
                # print 'strain in the middle of the measuring field'
            else:
                idx_border1 = self.idx_failure_crack[1]
                idx_border2 = self.idx_failure_crack[2]
                ux = field_data.ux_arr
                x_0 = field_data.x_arr_0
                eps_range = 1
                ux1 = np.mean(
                    ux[:, idx_border1 - eps_range: idx_border1 + eps_range], axis=1)
                ux2 = np.mean(
                    ux[:, idx_border2 - eps_range: idx_border2 + eps_range], axis=1)
                x_0_1 = np.mean(
                    x_0[:, idx_border1 - eps_range: idx_border1 + eps_range], axis=1)
                x_0_2 = np.mean(
                    x_0[:, idx_border2 - eps_range: idx_border2 + eps_range], axis=1)

                eps = (ux2 - ux1) / (x_0_2 - x_0_1)

                # eps = np.mean(ad.d_ux[:, idx_border1:idx_border2], axis=1)
                # print 'strain in the failure crack'

            # if 6layers
            x1 = pos_re1_6 - idx_6a * dis_fa
            # print 'x1', x1
            x_re1 = (dis_fa - x1) * (eps[idx_6a] - eps[idx_6b]) / dis_fa
            # print 'x_re1', x_re1
            eps_re1 = eps[idx_6b] + x_re1
            # print 'eps_re1', eps_re1
            eps_t_list.append(eps_re1)
            # print 'eps_t_list', np.array(eps_t_list, dtype='f')
        return np.array(eps_t_list, dtype='f')

    #-------------------------------------------------------------------------
    # get max and min strain of cross section
    #-------------------------------------------------------------------------

    eps_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_eps_t_aramis(self):

        ai = self.aramis_info
        if ai == None:
            return None

        field_data = self.aramis_field_data

        a = self.crack_bridge_strain_all
        n_fa = field_data.d_ux.shape[0]
        h = np.linspace(self.pos_fa[0], self.pos_fa[1], num=n_fa)
        eps_t_list = []
        eps_c_list = []

        # get eps
        for step, t in enumerate(self.t_aramis_cut):
            field_data.current_step = step
            if a == None:
                mid_idx = field_data.d_ux.shape[1] / 2
                eps_range = 3
                eps = np.mean(
                    field_data.d_ux[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
            else:
                idx_border1 = self.idx_failure_crack[1]
                idx_border2 = self.idx_failure_crack[2]
                ux = field_data.ux_arr
                x_0 = field_data.x_arr_0
                eps_range = 1
                ux1 = np.mean(
                    ux[:, idx_border1 - eps_range: idx_border1 + eps_range], axis=1)
                ux2 = np.mean(
                    ux[:, idx_border2 - eps_range: idx_border2 + eps_range], axis=1)
                x_0_1 = np.mean(
                    x_0[:, idx_border1 - eps_range: idx_border1 + eps_range], axis=1)
                x_0_2 = np.mean(
                    x_0[:, idx_border2 - eps_range: idx_border2 + eps_range], axis=1)

                eps = (ux2 - ux1) / (x_0_2 - x_0_1)

                # eps = np.mean(ad.d_ux[:, idx_border1:idx_border2], axis=1)

            # extrapolate eps for the specimen edges
            # 55mm bending specimen
            x = ((55 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
#             x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
            eps_ed_up = x + eps[-1]
            eps_ed_lo = eps[0] - x
            eps_to1 = np.append(eps, eps_ed_lo)
            eps_to2 = np.append(eps_ed_up, eps_to1)

            # eps_t_list.append(np.max(eps_to2))
            # eps_c_list.append(np.min(eps_to2))
            eps_t_list.append(eps_to2[0])
            eps_c_list.append(eps_to2[-1])

            # print 'eps_t_list', np.array(eps_t_list, dtype='f')
        return np.array(eps_t_list, dtype='f'), np.array(eps_c_list, dtype='f')



    #-------------------------------------------------------------------------
    # get strain(N) and strain(M)
    #-------------------------------------------------------------------------

    F_beg_threshold = Float(0.13, auto_set=False, enter_set=True)
    '''Threshold to find the beginning of raising F
    '''
    F_beg_idx = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_F_beg_idx(self):
        # get the index where the raising of F begins
        ai = self.aramis_info
        if ai == None:
            return None

        # print 'F', self.F_t_aramis
        # print 'N', self.N_t_aramis
        # print 't', self.t_aramis_cut

        F_idx = np.where(self.F_t_aramis > self.F_beg_threshold)[0]
        # print 'F_idx', F_idx

        if len(F_idx) == 0:
            F_beg_idx = []

        elif max(self.N_t_aramis) <= 2:
            F_beg_idx = 0
        else:
            F_beg_idx = F_idx[0]

        print 'F_beg_idx', F_beg_idx
        return F_beg_idx

    t_F_arr = Property(depends_on='data_file,aramis_resolution_key')
    '''Get the time where F rises
    '''
    @cached_property
    def _get_t_F_arr(self):
        ai = self.aramis_info
        if ai == None:
            return None

        if self.F_beg_idx == 0:
            t_F = self.t_aramis_cut
        elif self.F_beg_idx == []:
            t_F = []
        else:
            t_F = self.t_aramis_cut[self.F_beg_idx:]

        # print 't_F_arr', t_F
        return t_F


    eps_M = Property(depends_on='data_file,aramis_resolution_key')
    '''Get the strain corresponding to M
    '''
    @cached_property
    def _get_eps_M(self):
        ai = self.aramis_info
        if ai == None:
            return None

        field_data = self.aramis_field_data

        a = self.crack_bridge_strain_all
        n_fa = field_data.d_ux.shape[0]
        h = np.linspace(self.pos_fa[0], self.pos_fa[1], num=n_fa)
        t_F = self.t_F_arr
        F_beg_idx = self.F_beg_idx
        # print 'F_beg_idx', F_beg_idx
        eps_min_M_list = []
        eps_max_M_list = []
        max_step = self.n_steps
        # print 'max_step', max_step

        if t_F != []:
            for step in range(F_beg_idx, max_step, 1):
                field_data.current_step = step
                if a == None:
                    mid_idx = field_data.d_ux.shape[1] / 2
                    eps_range = 3
                    eps = np.mean(
                        field_data.d_ux[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)

                else:
                    idx_border1 = self.idx_failure_crack[1]
                    idx_border2 = self.idx_failure_crack[2]
                    eps = np.mean(
                        field_data.d_ux[:, idx_border1:idx_border2], axis=1)
                x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
                eps_ed_up = x + eps[-1]
                eps_ed_lo = eps[0] - x
                eps_to1 = np.append(eps, eps_ed_lo)
                eps_to2 = np.append(eps_ed_up, eps_to1)

                eps_max_M_list.append(np.max(eps_to2))
                eps_min_M_list.append(np.min(eps_to2))

            eps_max_M = np.array(eps_max_M_list, dtype='f')
            # print 'len_eps_max_M', len(eps_max_M)
            eps_min_M = np.array(eps_min_M_list, dtype='f')
            # print 'len_eps_min_M', len(eps_min_M)

            # print 'eps_max_M_list', eps_max_M
            print 'eps_max_M[-1]', eps_max_M[-1] * 1000
            print 'eps_max_M[0]', eps_max_M[0] * 1000

            # print 'eps_min_M_list', eps_min_M
            print 'eps_min_M[-1]', eps_min_M[-1] * 1000
            print 'eps_min_M[0]', eps_min_M[0] * 1000

            return eps_max_M, eps_min_M

    F_t_F = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_F_t_F(self):
        # get F in the area of t_F

        if self.F_beg_idx == []:
            # print 'F_tF', []
            return []
        else:
            F_t_F = self.F_t_aramis[self.F_beg_idx:]
            # print 'self.F_t_aramis[self.F_beg_idx:]', F_t_F
            # print 'len_F_t_F', len(F_t_F)
            print 'F_t_F[-1]', F_t_F[-1]
            return F_t_F

    M_t_F = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_M_t_F(self):
        # get M in the area of t_F

        if self.F_beg_idx == []:
            # print 'F_tF', []
            return []
        else:
            M_t_F = self.M_t_aramis[self.F_beg_idx:]
            # print 'self.F_t_aramis[self.F_beg_idx:]', F_t_F
            # print 'len_F_t_F', len(F_t_F)
            print 'M_t_F[-1]', M_t_F[-1]
            return M_t_F

    #-------------------------------------------------------------------------
    # get height of compression zone
    #-------------------------------------------------------------------------

    x_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_x_t_aramis(self):

        ai = self.aramis_info
        if ai == None:
            return None

        field_data = self.aramis_field_data
        a = self.crack_bridge_strain_all

        n_fa = field_data.d_ux.shape[0]
        dis_fa = self.meas_field[1] / (n_fa - 1)

        h = np.linspace(self.pos_fa[0], self.pos_fa[1], num=n_fa)
        h_1 = np.append(h, 0)
            # 55mm specimens
        h_2 = np.append(self.thickness * 1000., h_1)

        x_list = []

        # get eps
        for step, t in enumerate(self.t_aramis_cut):
            field_data.current_step = step
            if a == None:
                mid_idx = field_data.d_ux.shape[1] / 2
                eps_range = 3
                eps = np.mean(
                    field_data.d_ux[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
            else:
                idx_border1 = self.idx_failure_crack[1]
                idx_border2 = self.idx_failure_crack[2]
                eps = np.mean(
                    field_data.d_ux[:, idx_border1:idx_border2], axis=1)

            # extrapolate eps for the specimen edges
            # 55mm specimens
            x = ((self.thickness * 1000. - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
            eps_ed_up = x + eps[-1]
            eps_ed_lo = eps[0] - x
            eps_to1 = np.append(eps, eps_ed_lo)
            eps_to2 = np.append(eps_ed_up, eps_to1)

            # print 'eps_to2', eps_to2

            idx_neg = np.array(np.where(eps_to2 < 0)[0])
            # print 'idx_neg', idx_neg
            d = idx_neg.shape[0]
            # print 'shape', d

            if d == 0:
                x = 0

            elif eps[-1] > 0:
                x = 0

            elif idx_neg[0] == 0:
                x = self.thickness * 1000.

            else:
                idx_1 = idx_neg[0]
                # print 'idx_1', idx_1
                idx_2 = idx_1 - 1
                # print 'idx_2', idx_2

                x_a = (
                    dis_fa * abs(eps_to2[idx_1])) / (eps_to2[idx_2] - eps_to2[idx_1])
                # print 'x_a', x_a
                x = h_2[idx_1] + x_a
                # print 'x', x

            x_list.append(x)
            # print 'x', np.array(x_list, dtype='f')

        return np.array(x_list, dtype='f')

    #---------------------------------
    # view
    #---------------------------------

    traits_view = View(VSplit(
        HSplit(Group(
            Item('width', format_str="%.3f"),
            Item('length', format_str="%.3f"),
            Item('start_time_aramis', format_str="%.0f"),
            springy=True,
            label='geometry',
            id='matresdev.db.exdb.ex_composite_bending_tensile_test.geometry',
            dock='tab',
        ),
            Group(
            Item('loading_rate'),
            Item('age'),
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
        #                               label = 'input variables',
        #                               id = 'matresdev.db.exdb.ex_composite_tensile_test.vgroup.inputs',
        #                               dock = 'tab',
        #                               scrollable = True,
        #                               ),
        Group(
            Item('E_c', style='readonly', show_label=True, format_str="%.0f"),
            Item('F_max', style='readonly',
                 emphasized=True, format_str="%.2f"),
            Item('w_max', style='readonly',
                 emphasized=True, format_str="%.2f"),
            Item('M_max', style='readonly',
                 emphasized=True, format_str="%.3f"),
            label='output characteristics',
            id='matresdev.db.exdb.ex_composite_bending_tensile_test.vgroup.outputs',
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


ExpBT4PTAramis2d.db = ExRunClassExt(klass=ExpBT4PTAramis2d)
# ExpBTTDB.db = ExRunClassExt(klass=ExpBTTDB)

#--------------------------------------------------------------

if __name__ == '__main__':

    from matresdev.db.simdb import SimDB
    from matresdev.db.exdb import ExRunView
    simdb = SimDB()
    import os

    test_files = ['BT-1C-55mm-0-3300EP-V2_S3P2(11)-Aramis2d.DAT']

    test_file_path = os.path.join(simdb.exdata_dir,
                                  'bending_tests', 'four_point',
                                  '2015-09-02_BT-1C-55mm-0-3300SBR_cyc-Aramis2d',
                                  )

    doe_reader = ExRunView(data_file=test_file_path)
    doe_reader.configure_traits()

#     ExpBT4PTAramis2d.db.configure_traits()
#    to see all experiments in one picture
