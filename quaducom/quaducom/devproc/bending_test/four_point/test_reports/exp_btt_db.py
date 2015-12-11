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
from matresdev.db.matdb.trc.composite_cross_section \
    import CompositeCrossSection, plain_concrete
from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp
import numpy as np


class ExpBTTDB(ExType):

    '''Experiment: Bending Tensile Test Dog Bone
    '''
#    label = Str('dog bone tensile test')

    implements(IExType)

    # ------------------------------------------------------- -------------
    # register a change of the traits with metadata 'input'
    # --------------------------------------------------------------------

    input_change = Event

    @on_trait_change('+input, ccs.input_change')
    def _set_input_change(self):
        print '*** raising input change in CTT'
        self.input_change = True

    #--------------------------------------------------------------------------
    # specify inputs:
    #--------------------------------------------------------------------------

    width = Float(0.100, unit='m', input=True, table_field=True,
                  auto_set=False, enter_set=True)
    length = Float(0.35, unit='m', input=True, table_field=True,
                   auto_set=False, enter_set=True)
    # start of 2D-Aramis
    start_time_aramis = Float(5, unit='m', input=True, table_field=True,
                              auto_set=False, enter_set=True)
    # age of the concrete at the time of testing
    age = Int(28, unit='d', input=True, table_field=True,
              auto_set=False, enter_set=True)
    loading_rate_N = Float(2.0, unit='mm/min', input=True, table_field=True,
                           auto_set=False, enter_set=True)
    loading_rate_F = Float(4.0, unit='kN/10min', input=True, table_field=True,
                           auto_set=False, enter_set=True)
    # date of the production of the specimen
    prodcution_date = Float(27.0, unit='DD.MM', input=True, table_field=True,
                            auto_set=False, enter_set=True)
    # date of the testing of the specimen
    testing_date = Float(27.0, unit='DD.MM', input=True, table_field=True,
                         auto_set=False, enter_set=True)

    #--------------------------------------------------------------------------
    # composite cross section
    #--------------------------------------------------------------------------

    ccs = Instance(CompositeCrossSection)

    def _ccs_default(self):
        '''default settings correspond to
        setup '9u_MAG-07-03_PZ-0708-1'
        '''
        print 'ccs default used'
        fabric_layout_key = '2D-05-11'
        concrete_mixture_key = 'barrelshell'
        orientation_fn_key = 'all0'
        n_layers = 6
        thickness = 0.02

        s_tex_z = thickness / (n_layers + 1)
        ccs = CompositeCrossSection(
            fabric_layup_list=[
                plain_concrete(s_tex_z * 0.5),
                FabricLayUp(
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
    # Indicate whether the test is suitable and prepared for
    # calibration.
    #--------------------------------------------------------------------------
    ready_for_calibration = Property(Bool)

    def _get_ready_for_calibration(self):
        # return False by default
        # the subclasses shall overload this
        # and define the rules
        return self.ccs.is_regular

    #--------------------------------------------------------------------------
    # Get properties of the composite
    #--------------------------------------------------------------------------

    # E-modulus of the composite at the time of testing
    E_c = Property(
        Float, unit='MPa', depends_on='input_change', table_field=True)

    @cached_property
    def _get_E_c(self):
        return self.ccs.get_E_c_time(self.age)

    # E-modulus of the concrete at the time of testing
    E_m = Property(
        Float, unit='MPa', depends_on='input_change', table_field=True)

    @cached_property
    def _get_E_m(self):
        return self.ccs.get_E_m_time(self.age)

    # cross-sectional-area of the composite
    A_c = Property(Float, unit='m^2', depends_on='input_change')

    @cached_property
    def _get_A_c(self):
        return self.width * self.ccs.thickness

    # total cross-sectional-area of the textile reinforcement
    A_tex = Property(Float, unit='mm^2', depends_on='input_change')

    @cached_property
    def _get_A_tex(self):
        return self.ccs.a_tex * self.width

    # E-modulus of the composite after 28 days
    E_c28 = DelegatesTo('ccs', listenable=False)

    # reinforcement ration of the composite
    rho_c = DelegatesTo('ccs', listenable=False)

    #--------------------------------------------------------------------------
    # define processing
    #--------------------------------------------------------------------------

    def process_source_data(self):
        '''read in the measured data from file and assign
        attributes after array processing.
        If necessary modify the assigned data, e.i. change
        the sign or specify an offset for the specific test setup.
        '''
        super(ExpBTTDB, self).process_source_data()

        # ad 0.265 kN for the weight of the test set up (26,5 kg)
        self.F = self.Kraft + 0.265
        self.F = self.Weg
        self.t = self.Bezugskanal

        if hasattr(self, 'WA_Mitte'):
            self.w = self.WA_Mitte

        print 'wwww',
        print self.w

        self.w -= self.w[0]
        self.w *= -1

    #--------------------------------------------------------------------------
    # Get the maximum force index to cut off the descending part of the curves
    #--------------------------------------------------------------------------

    max_F_idx = Property(Int, depends_on='input_change')

    @cached_property
    def _get_max_F_idx(self):
        return np.argmax(self.F)

    max_N_idx = Property(Int, depends_on='input_change')

    @cached_property
    def _get_max_N_idx(self):
        return np.argmax(self.N)

    F_max1 = Property(Float, depends_on='input_change',
                      output=True, table_field=True, unit='kN')

    @cached_property
    def _get_F_max1(self):
        '''Only needed for the distinction between pure tensile test and bending test
        '''
        return max(self.F)

    #--------------------------------------------------------------------------
    # Get only the ascending branch of the response curve
    #--------------------------------------------------------------------------

    F_max_threshold = Float(0.15, auto_set=False, enter_set=True)
    '''Threshold to distinguish between pure tensile test and bending test.
    '''

    N_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_N_asc(self):
        if self.F_max1 >= self.F_max_threshold:
            # maximum load is determined by the transverse force
            return self.N[:self.max_F_idx + 1]
        elif self.F_max1 <= self.F_max_threshold:
            # maximum load is determined by the normal force
            return self.N[:self.max_N_idx + 1]

    F_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_F_asc(self):
        if self.F_max1 >= self.F_max_threshold:
            # maximum load is determined by the transverse force
            return self.F[:self.max_F_idx + 1]
        elif self.F_max1 <= self.F_max_threshold:
            # maximum load is determined by the normal force
            return self.F[:self.max_N_idx + 1]

    t_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_t_asc(self):
        if self.F_max1 >= self.F_max_threshold:
            # maximum load is determined by the transverse force
            return self.t[:self.max_F_idx + 1]
        elif self.F_max1 <= self.F_max_threshold:
            # maximum load is determined by the normal force
            return self.t[:self.max_N_idx + 1]

    w_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_w_asc(self):
        if self.F_max1 >= self.F_max_threshold:
            # maximum load is determined by the transverse force
            return self.w[:self.max_F_idx + 1]
        elif self.F_max1 <= self.F_max_threshold:
            # maximum load is determined by the normal force
            return self.w[:self.max_N_idx + 1]

    u_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_u_asc(self):
        if self.F_max1 >= self.F_max_threshold:
            # maximum load is determined by the transverse force
            return self.u[:self.max_F_idx + 1]
        elif self.F_max1 <= self.F_max_threshold:
            # maximum load is determined by the normal force
            return self.u[:self.max_N_idx + 1]

    #--------------------------------------------------------------------------
    # Cut the curves at the end when the gradient of displacement is to big
    #--------------------------------------------------------------------------

    gra_max_threshold = Float(0.12, auto_set=False, enter_set=True)
    '''Threshold to limit the gradient of displacement-curve at the end.
    '''

    w_cut_idx = Property(Int, depends_on='input_change')

    @cached_property
    def _get_w_cut_idx(self):

        w_asc = self.w_asc
        t_asc = self.t_asc

        # Calculate the deltas with beginning at indize 840 that is equivalent
        # to t = 420 sec. The calculation of deltas start at t = 420 sec because there are
        # sometimes big gradients in the beginning or during the test due to the controlling
        # of the bending force by hand
        delta_w_arr = w_asc[841:] - w_asc[840:-1]
        delta_t_arr = t_asc[841:] - t_asc[840:-1]

        # Calculate the gradient for every index
        gra_arr = delta_w_arr[:] / delta_t_arr[:]

        # Examine the indices where the gradient is bigger than the threshold
        # gradient
        gra_idx_arr = np.where(gra_arr > self.gra_max_threshold)[0]

#        print '*** gradient is bigger than the threshold gradient at the following indices and time: ***'
#        print 'gra_idx_arr', gra_idx_arr
#        print 'gra_rr', gra_arr[gra_idx_arr]
#        print 'w_asc[gra_idx_arr]', w_asc[gra_idx_arr]

        if len(gra_idx_arr) > 0:
            return gra_idx_arr[0] + 840
        else:
            return len(self.w_asc)

    N_cut_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_N_cut_asc(self):
        ''' Method to cut the end of the displacement curve if the gradient is
        to big due to the kind of failure in the testing
        '''
        return self.N[:self.w_cut_idx]

    F_cut_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_F_cut_asc(self):
        ''' Method to cut the end of the displacement curve if the gradient is
        to big due to the kind of failure in the testing
        '''
        return self.F[:self.w_cut_idx]

    t_cut_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_t_cut_asc(self):
        ''' Method to cut the end of the displacement curve if the gradient is
        to big due to the kind of failure in the testing
        '''
        return self.t[:self.w_cut_idx]

    w_cut_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_w_cut_asc(self):
        ''' Method to cut the end of the displacement curve if the gradient is
        to big due to the kind of failure in the testing
        '''
        return self.w_asc[:self.w_cut_idx]

    u_cut_asc = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_u_cut_asc(self):
        ''' Method to cut the end of the displacement curve if the gradient is
        to big due to the kind of failure in the testing
        '''
        return self.u_asc[:self.w_cut_idx]

    #--------------------------------------------------------------------------
    # Method to eliminate the negative or positive deformation at the beginning
    # of the experiment (while F ~ 0.1 kN (constant) and N is increasing)
    # due to predeformation of the specimen
    #--------------------------------------------------------------------------

    # get the minimum value of deformation that is equivalent to the
    # predeformation
    w_pred = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_w_pred(self):
        w_pred = min(self.w_cut_asc)
        return w_pred

    F_lim_threshold = Float(0.142, auto_set=False, enter_set=True)
    '''Threshold for the lower limit of force. For all force values smaller
        than this limit, equate the deformation value with zero
    '''

    w_el_pred = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_w_el_pred(self):
        w = self.w_cut_asc
        F = self.F_cut_asc

        # get the indices for the force smaller than force threshold
        idx_F_lim = np.where(F <= self.F_lim_threshold)[0]

        if len(idx_F_lim) > 0:
            idx = idx_F_lim[-1]

            w_lim = w[idx_F_lim[-1]]

            w_0 = w[0:idx]
            w_0 = np.zeros_like(w_0)
            w_1 = w[idx:]
            w_1 = w_1 - w_lim
            w_el_pred = np.append(w_0, w_1)

            return w_el_pred

        else:
            return w

    #--------------------------------------------------------------------------
    # Get the moment arrays
    #--------------------------------------------------------------------------

    M = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_M(self):
        'resulting moment'
        return self.F_cut_asc * self.length / 4 - self.N_cut_asc * self.w_el_pred / 1000

    MN = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_MN(self):
        'moment due to normal force an eccentricity'
        return self.N_cut_asc * self.w_el_pred / 1000 * -1

    MF = Property(Array('float_'), depends_on='input_change')

    @cached_property
    def _get_MF(self):
        'moment due to bending force'
        return self.F_cut_asc * self.length / 4

    #--------------------------------------------------------------------------
    # Get maximum values of the variables
    #--------------------------------------------------------------------------

    N_max = Property(Float, depends_on='input_change',
                     output=True, table_field=True, unit='kN')

    @cached_property
    def _get_N_max(self):
        return self.N_cut_asc[-1]

    F_max = Property(Float, depends_on='input_change',
                     output=True, table_field=True, unit='kN')

    @cached_property
    def _get_F_max(self):
        return self.F_cut_asc[-1]

    t_max = Property(Float, depends_on='input_change',
                     output=True, table_field=True, unit='kNm')

    @cached_property
    def _get_t_max(self):
        return self.t_cut_asc[-1]

    w_max = Property(Float, depends_on='input_change',
                     output=True, table_field=True, unit='mm')

    @cached_property
    def _get_w_max(self):
        return self.w_el_pred[-1]

    u_max = Property(Float, depends_on='input_change',
                     output=True, table_field=True, unit='mm')

    @cached_property
    def _get_u_max(self):
        return self.u_cut_asc[-1]

    M_max = Property(Float, depends_on='input_change',
                     output=True, table_field=True, unit='kNm')

    @cached_property
    def _get_M_max(self):
        return self.M[-1]

    MN_max = Property(Float, depends_on='input_change',
                      output=True, table_field=True, unit='kNm')

    @cached_property
    def _get_MN_max(self):
        return self.MN[-1]

    MF_max = Property(Float, depends_on='input_change',
                      output=True, table_field=True, unit='kNm')

    @cached_property
    def _get_MF_max(self):
        return self.MF[-1]

    #--------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------

    plot_templates = {'N(t), F(t)': '_plot_N_F_t',
                      'N_cut(t), F_cut(t)': '_plot_N_F_t_cut',
                      'N(u)_cut': '_plot_N_u_cut',
                      'w(t)': '_plot_w_t',
                      'w_cut(t)': '_plot_w_t_cut',
                      'w_el_pred(t)': '_plot_w_el_pred_t',
                      'u(t)': '_plot_u_t',
                      'u_cut(t)': '_plot_u_t_cut',
                      'F(w)_cut': '_plot_F_w_cut',
                      'F(w_el_pred)': '_plot_F_w_el_pred',
                      'M(t), M_II(t), M_0(t)': '_plot_M_MN_MF_t',
                      'N(M), N(M_II), N(M_0)': '_plot_N_M_MN_MF',
                      'N-M': '_plot_N_M',
                      }

    default_plot_template = 'N(t), F(t)'

    def _plot_N_F_t(self, axes):
        '''Normal force and bending force versus time
        '''
        axes.plot(self.t_asc, self.N_asc, color='blue', label='N')
        axes.plot(self.t_asc, self.F_asc, color='red', label='F')
        axes.xaxis.tick_bottom()
        axes.set_ylim(0, 50)
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('N / F[kN]')
        # ax2.set_ylim(0, 3)
        # ax2.set_ylabel('F [kN]')
        axes.legend(loc=2)

    def _plot_N_F_t_cut(self, axes):
        '''Normal force_cut and bending force_cut versus time_cut
        '''
        axes.plot(self.t_cut_asc, self.N_cut_asc, color='darkblue', label='N')
        axes.grid()
        axes.set_ylim(0, 55)
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('N [kN]')
        # ax2 = axes.twinx()
        axes.plot(self.t_cut_asc, self.F_cut_asc, color='darkred', label='F')
        # axes.set_ylim(0, 5.5)
        axes.set_ylabel('F [kN]')
        # axes.legend(loc=9)
        axes.legend(loc=2)

    def _plot_N_u_cut(self, axes):
        '''Normal force_cut versus displacement_u_cut
        '''
        # ax1 = axes
        axes.plot(self.u_cut_asc, self.N_cut_asc, color='blue', label='N')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('u [mm]')
        axes.set_ylabel('N [kN]')
        axes.set_ylim(0, 50)
        axes.legend(loc=2)

    def _plot_w_t(self, axes):
        '''Displacement versus time
        '''
        print 'self.w_asc', self.w_asc
        axes.plot(self.t_asc, self.w_asc, color='black')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('w [mm]')
        axes.set_ylim(0, 10)

    def _plot_w_t_cut(self, axes):
        '''Displacement_cut versus time_cut
        '''
        axes.plot(self.t_cut_asc, self.w_cut_asc, color='black')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('w [mm]')
        axes.set_ylim(0, 10)

    def _plot_w_el_pred_t(self, axes):
        '''Displacement with eliminated predeformation versus time_cut
        '''
        axes.plot(self.t_cut_asc, self.w_el_pred, color='black')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('w [mm]')
        axes.set_ylim(0, 10)

    def _plot_u_t(self, axes):
        '''Displacement versus time
        '''
        axes.plot(self.t_asc, self.u_asc, color='black')
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('u [mm]')
        axes.set_ylim(0, 30)

    def _plot_u_t_cut(self, axes):
        '''Displacement_cut versus time_cut
        '''
        axes.plot(self.t_cut_asc, self.u_cut_asc, color='black')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('u [mm]')
        axes.set_ylim(0, 30)

    def _plot_F_w_cut(self, axes):
        '''Bending force_cut versus displacement_cut
        '''
        axes.plot(self.w_cut_asc, self.F_cut_asc, color='blue')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('w [mm]')
        axes.set_ylabel('F [kN]')
        axes.set_ylim(0, 6)
        axes.set_xlim(0, 10)

    def _plot_F_w_el_pred(self, axes):
        '''Normal force_cut versus displacement with eliminated predeformation
        '''
        axes.plot(self.w_el_pred, self.F_cut_asc, color='blue')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('w [mm]')
        axes.set_ylabel('F [kN]')
        axes.set_ylim(0, 6)
        axes.set_xlim(0, 10)

    def _plot_M_MN_MF_t(self, axes, *args, **kw):
        '''Moment versus time
        '''
        axes.plot(self.t_cut_asc, self.N_cut_asc, color='darkblue', label='N')
        axes.grid()
        axes.set_ylim(0, 55)
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('N [kN]')
        axes.legend(loc=2)
        ax2 = axes.twinx()
        ax2.plot(self.t_cut_asc, self.M, color='green', label='M')
        ax2.plot(self.t_cut_asc, self.MF, color='darkred', label='M_0')
        ax2.plot(self.t_cut_asc, self.MN, color='aqua', label='M_II')
        ax2.set_ylabel(' M [kNm]')
        ax2.set_ylim(-0.33, 0.58)
        ax2.legend(loc=1)
        ax2.legend(ncol=3)

    def _plot_N_M_MN_MF(self, axes):
        '''Normal force versus moment
        '''
        axes.plot(self.MF, self.N_cut_asc, color='red', label='M_0')
        axes.plot(self.MN, self.N_cut_asc, color='blue', label='M_II')
        axes.plot(self.M, self.N_cut_asc, color='green', label='M')
        # axes.xaxis.tick_top()
        axes.set_xlabel('M / M_II / M_0 [kNm]')
        axes.xaxis.set_label_position('top')
        axes.set_ylabel('N [kN]')
        axes.set_ylim(24, 0)
        axes.set_xlim(-0.25, 0.45)
        axes.legend(loc=2)

    # def _plot_N_M(self, axes):
        '''M-M-interaction diagramm
        '''
        '''axes.plot(self.M, self.N_cut_asc, color='green', label='M')
        x = [0, 0.35]
        y = [44, 0]
        axes.plot(x, y, color='grey', linestyle='--', label='M-N-Interaction')
        axes.grid()
        axes.xaxis.tick_top()
        axes.xaxis.set_label_position('top')
        axes.set_ylim(55 , -1)
        axes.set_xlim(-0.01 , 0.4)
        axes.set_xlabel('M [kNm]')
        axes.set_ylabel('N [kN]')
        axes.legend(loc=4)'''

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
        step_times = self.aramis_field_data.step_times + 5
        # step_times = self.aramis_field_data.step_times + self.start_time_aramis
        print 'self.start_time_aramis', self.start_time_aramis
        print '-----------------------------t_max_ARAMIS', step_times[-1]
        t_max = self.t[self.w_cut_idx]
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

    aramis_info = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_aramis_info(self):
        af = self.get_cached_aramis_file(self.aramis_resolution_key)
        if af == None:
            return None
        return AramisInfo(data_dir=af)

    aramis_field_data = Property(depends_on='data_file,aramis_resolution_key')
    '''Field data including strains and displacements.
    '''
    @cached_property
    def _get_aramis_field_data(self):
        t_fail = self.t_cut_asc[-1]
        ad = AramisFieldData(aramis_info=self.aramis_info,
                             integ_radius=3)
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
                         aramis_data=ad,
                         ddd_ux_avg_threshold=-0.5e-3,
                         ddd_ux_threshold=-0.5e-3)

    N_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_N_t_aramis(self):
        'normal force interpolated to the time steps of aramis'
        # print 'np.interp(self.t_aramis, self.t, self.N)',
        # np.interp(self.t_aramis, self.t, self.N)
        return np.interp(self.t_aramis_cut, self.t, self.N)

    F_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_F_t_aramis(self):
        'bending force interpolated to the time steps of aramis'
        return np.interp(self.t_aramis_cut, self.t, self.F)

    w_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_w_t_aramis(self):
        'displacement (with eliminated predeformation) interpolated to the time steps of aramis'
        return np.interp(self.t_aramis_cut, self.t_cut_asc, self.w_el_pred)

    M_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_M_t_aramis(self):
        'resulting moment interpolated to the time steps of aramis'
        return self.F_t_aramis * self.length / 4 - self.N_t_aramis * self.w_t_aramis / 1000

    MN_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_MN_t_aramis(self):
        'moment due to normal force interpolated to the time steps of aramis'
        return self.N_t_aramis * self.w_t_aramis / 1000 * -1

    MF_t_aramis = Property(depends_on='data_file,aramis_resolution_key')

    @cached_property
    def _get_MF_t_aramis(self):
        'moment due to bending force interpolated to the time steps of aramis'
        return self.F_t_aramis * self.length / 4

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

    h_re1_6_threshold = Float(2.86, auto_set=False, enter_set=True)
    '''Threshold for position of first reinforcement layer (6 layers).
    '''

    h_re1_4_threshold = Float(4.0, auto_set=False, enter_set=True)
    '''Threshold for position of first reinforcement layer (4 layers).
    '''

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

        h_dis = (20 - self.meas_field[1]) / 2
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

        # position of first node in y-direction
        pos_no_f = 20.0 - self.h_dis
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
        pos_re1_6 = self.h_re1_6_threshold - self.h_dis
        # print 'pos_re1_6', pos_re1_6  # -> negative value
        pos_re1_4 = self.h_re1_4_threshold - self.h_dis
        # print 'pos_re1_4', pos_re1_4

        # indices of strain next to position of reinforcement layer
        if pos_re1_6 < 0:
            idx_6a = 0
        else:
            idx_6a = pos_re1_6 / dis_fa
            idx_6a = round(idx_6a)
        idx_6b = idx_6a + 1

        if pos_re1_4 < 0:
            idx_4a = 0
        else:
            idx_4a = pos_re1_4 / dis_fa
            idx_4a = round(idx_4a)
        idx_4b = idx_4a + 1

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

        # if 4 layers
        # x1 = pos_re1_4 - idx_4a * dis_fa
        # print 'x1', x1
        # x_re1 = (dis_fa - x1) * (eps[idx_4a] - eps[idx_4b]) / dis_fa
        # print 'x_re1', x_re1
        # eps_re1 = eps[idx_4b] + x_re1
        # print 'eps_re1', eps_re1
        # eps_t_list.append(eps_re1)
        # print 'eps_t_list', np.array(eps_t_list, dtype='f')
        # return np.array(eps_t_list, dtype='f')

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
            x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
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

    # TO-DO:
    #--------------------------------------------------------------------------
    # Cut the strain  at the end when the gradient of displacement is to big
    #--------------------------------------------------------------------------

    # gra_max_eps_threshold = Float(0.0002, auto_set=False, enter_set=True)
    '''Threshold to limit the gradient of displacement-curve at the end.
    '''

    '''eps_cut_idx = Property(Int, depends_on='input_change')
    @cached_property
    def _get_eps_cut_idx(self):

        eps_max_arr = self.eps_t_aramis[0]
        print 'eps_max_arr', eps_max_arr
        t_arr = self.t_aramis_cut

        delta_eps_arr = eps_max_arr[1:] - eps_max_arr[:-1]
        print 'delta_eps_arr', delta_eps_arr
        delta_t_arr = t_arr[1:] - t_arr[:-1]
        print 'delta_t_arr', delta_t_arr

        # Calculate the gradient for every index
        gra_arr = delta_eps_arr[:] / delta_t_arr[:]
        print 'gra_rr', gra_arr

        # Examine the indices where the gradient is bigger than the threshold gradient
        gra_idx_arr = np.where(gra_arr > self.gra_max_eps_threshold)[0]

        print '*** gradient is bigger than the threshold gradient at the following indices and time: ***'
        print 'gra_idx_arr', gra_idx_arr

        if len(gra_idx_arr) > 0:
            return gra_idx_arr[0]
        else:
            return len(self.eps_t_aramis[0])'''

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

    t_N_arr = Property(depends_on='data_file,aramis_resolution_key')
    '''Get the time where N rises
    '''
    @cached_property
    def _get_t_N_arr(self):
        ai = self.aramis_info
        if ai == None:
            return None

        if self.F_beg_idx == 0:
            t_N = []
        elif self.F_beg_idx == []:
            t_N = self.t_aramis_cut
        else:
            N_end_idx = self.F_beg_idx
            # print 'N_end_idx', N_end_idx
            t_N = self.t_aramis_cut[0: N_end_idx]

        # print 't_N_arr', t_N
        return t_N

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

    eps_N = Property(depends_on='data_file,aramis_resolution_key')
    '''get the strain corresponding to N
    '''
    @cached_property
    def _get_eps_N(self):
        ai = self.aramis_info
        if ai == None:
            return None

        field_data = self.aramis_field_data

        a = self.crack_bridge_strain_all
        n_fa = field_data.d_ux.shape[0]
        h = np.linspace(self.pos_fa[0], self.pos_fa[1], num=n_fa)
        t_N = self.t_N_arr
        eps_N_list = []

        if len(t_N) == 0:
            N_end_idx = 0
        elif self.F_beg_idx == []:
            N_end_idx = len(self.N_t_aramis)
        else:
            N_end_idx = self.F_beg_idx
        print 'N_end_idx', N_end_idx

        if t_N != []:
            for step in range(0, N_end_idx, 1):
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
                    # print 'eps', eps

                    # extrapolate eps for the specimen edges
                    x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
                    eps_ed_up = x + eps[-1]
                    eps_ed_lo = eps[0] - x
                    eps_to1 = np.append(eps, eps_ed_lo)
                    eps_to2 = np.append(eps_ed_up, eps_to1)

                    # print 'eps_to2', eps_to2

                    eps_N_list.append(np.mean(eps_to2))
                    eps_N = np.array(eps_N_list, dtype='f')

            # print 'eps_N_list', eps_N
            print 'eps_N[-1]', eps_N[-1] * 1000
            return eps_N

    N_t_N = Property(depends_on='data_file,aramis_resolution_key')
    '''Get N in the range of t_N
    '''
    @cached_property
    def _get_N_t_N(self):

        if self.t_N_arr == None or len(self.t_N_arr) == 0:
            # print 'self.N_t_aramis[0:t_N_idx]', []
            return []
        elif self.F_beg_idx == []:
            N_t_N = self.N_t_aramis
            return N_t_N
        else:
            N_end_idx = self.F_beg_idx
            # print 'N_end_idx' , N_end_idx
            N_t_N = self.N_t_aramis[0:N_end_idx]
            # print 'self.N_t_aramis[0:t_N_idx]', N_t_N
            print 'N_t_N[-1]', N_t_N[-1]
            return N_t_N

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
        h_2 = np.append(20, h_1)

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
            x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
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
                x = 20

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

    #-------------------------------------------------------------------------
    # get curvature of the specimen
    #-------------------------------------------------------------------------

    # TO-DO: check if correct
    '''cu_t_aramis = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_cu_t_aramis(self):

        ai = self.aramis_info
        if ai == None:
            return None

        a = self.crack_bridge_strain_all
        eps = self.eps1_t_aramis[0]
        z_arr = 20 - self.x_t_aramis
        cu_list = []
        # print 'self.x_t_aramis', self.x_t_aramis
        # print 'z_arr', z_arr
        # print 'self.eps1_t_aramis[0]', self.eps1_t_aramis[0]

        if a == None:
            return None

        else:
            for i in z_arr:
                if i == 0:
                    cu = 0
                else:
                    cu = eps[i] / z_arr[i]
                cu_list.append(cu)

        return np.array(cu_list, dtype='f')'''

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
            Item('loading_rate_N'),
            Item('loading_rate_F'),
            Item('age'),
            Item('prodcution_date'),
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
        #                               label = 'input variables',
        #                               id = 'matresdev.db.exdb.ex_composite_tensile_test.vgroup.inputs',
        #                               dock = 'tab',
        #                               scrollable = True,
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


ExpBTTDB.db = ExRunClassExt(klass=ExpBTTDB)

#--------------------------------------------------------------

if __name__ == '__main__':

    from matresdev.db.simdb import SimDB
    from matresdev.db.exdb import ExRunView
    simdb = SimDB()
    import os

    ex_path = os.path.join(simdb.exdata_dir,
                           'bending_tensile_test',
                           '2014-06-12_BTT-6c-2cm-0-TU_MxN2',
                           'BTT-6c-2cm-TU-0-V03_MxN2.DAT')

    test_file = os.path.join(simdb.exdata_dir,
                             'bending_tensile_test',
                             '2014-06-12_BTT-4c-2cm-0-TU_MxN2',
                             'BTT-4c-2cm-TU-0-V02_MxN2.DAT')

    doe_reader = ExRunView(data_file=ex_path)
    doe_reader.configure_traits()

    # ExpBTTDB.db.configure_traits()
    # to see all experiments in one picture
