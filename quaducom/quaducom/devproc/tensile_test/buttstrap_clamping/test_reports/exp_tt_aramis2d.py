'''
Script evaluating the aramis data
for tensile tests
Created on Jan 4, 2015 by scho
'''

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

from quaducom.devproc.tensile_test.dog_bone import ExpTTDB


from matresdev.db.matdb.trc.composite_cross_section \
    import CompositeCrossSection, plain_concrete
from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp
import numpy as np


class ExpTTAramis2d(ExpTTDB):
    '''Experiment: Tensile Test with Aramis2d measurements
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
        super(ExpTTAramis2d, self).process_source_data()

        self.F = self.Kraft
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
        print('self.start_time_aramis', self.start_time_aramis)
        print('-----------------------------t_max_ARAMIS', step_times[-1])
        t_max = self.t[-1]
        print('-----------------------------t_max_ASCII', t_max)
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
        print('-----------------------------n_steps_cut2', len(self.t_aramis_cut))
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
        print('_get_aramis_field_data')
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
        print('_get_aramis_cdt')
        ad = self.aramis_field_data
        crack_detection_step = self.crack_detaction_step
        # if crack-detections-step is set to current step:
        if self.crack_detaction_current_step:
            crack_detection_step = ad.current_step
        print('---------------------CRACK DETECTION STEP', crack_detection_step)
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
        print('crack_idx_arr', crack_idx_arr)

        # get the indices of the middle between two cracks
        crack_mid_idx_arr = (crack_idx_arr[0:-1] + crack_idx_arr[1:]) / 2
        print('crack_mid_idx_arr', crack_mid_idx_arr)
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
        print('idx_failure_crack', idx_failure_crack)
        print('shape_ad.x_arr_0', np.shape(field_data.x_arr_0))
        print('number_failure_crack', max_idx + 2)
        idx_border1 = crack_mid_idx_arr[max_idx]
        print('idx_border1 ', idx_border1)
        idx_border2 = crack_mid_idx_arr[max_idx + 1]
        print('idx_border2 ', idx_border2)

        return idx_failure_crack, idx_border1, idx_border2

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
        Group(
            Item('E_c', style='readonly', show_label=True, format_str="%.0f"),
            Item('F_max', style='readonly',
                 emphasized=True, format_str="%.2f"),
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


ExpTTAramis2d.db = ExRunClassExt(klass=ExpTTAramis2d)
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
