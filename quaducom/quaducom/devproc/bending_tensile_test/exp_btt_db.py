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
# Created on Feb 15, 2010 by: rch, ascholzen

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

from etsproxy.traits.api import \
    Int, Float, \
    on_trait_change, Instance, \
    Array, Property, cached_property, \
    Bool, Event, implements, \
    DelegatesTo, Str

import numpy as np

from etsproxy.traits.ui.api \
    import View, Item, HSplit, Group, VSplit

from matresdev.db.exdb.ex_type import ExType
from matresdev.db.exdb.i_ex_type import IExType

from aramis_cdt import AramisInfo, AramisData, AramisBSA, AramisUI

from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp

from matresdev.db.matdb.trc.composite_cross_section \
    import CompositeCrossSection, plain_concrete

from matresdev.db.exdb.ex_run_table import ExRunClassExt

class ExpBTTDB(ExType):
    '''Experiment: Bending Tensile Test Dog Bone
    '''
#    label = Str('dog bone tensile test')

    implements(IExType)

    #--------------------------------------------------------------------
    # register a change of the traits with metadata 'input'
    #--------------------------------------------------------------------

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
    E_c = Property(Float, unit='MPa', depends_on='input_change', table_field=True)
    @cached_property
    def _get_E_c(self):
        return self.ccs.get_E_c_time(self.age)

    # E-modulus of the concrete at the time of testing
    E_m = Property(Float, unit='MPa', depends_on='input_change', table_field=True)
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
        self.N = self.F160kN + 0.265
        self.F = self.Zquer
        self.t = self.Bezugskanal

        if hasattr(self, 'W10_Ho'):
            self.w = self.W10_Ho
        elif hasattr(self, 'W10_li'):
            self.w = self.W10_li
        else:
            self.w = np.zeros_like(self.F)

        if hasattr(self, 'Weg'):
            self.u = self.Weg
        else:
            self.u = np.zeros_like(self.F)

        # u is the machine control displacement corresponding the F160kN

    # get the first value of the displacement gauges
    # used to reset the displacement gauges if they do not start at zero
        self.w -= self.w[0]
        self.w *= -1

        self.u -= self.u[0]

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
    '''Threshold to distniguish between pure tensile test and bending test.
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

        # Examine the indices where the gradient is bigger than the threshold gradient
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

    # get the minimum value of deformation that is equivalent to the predeformation
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

    plot_templates = {'N(t), F(t)' : '_plot_N_F_t',
                      'N_cut(t), F_cut(t)' : '_plot_N_F_t_cut',
                      'N(u)_cut': '_plot_N_u_cut',
                      'w(t)' : '_plot_w_t',
                      'w_cut(t)' : '_plot_w_t_cut',
                      'w_el_pred(t)' : '_plot_w_el_pred_t',
                      'u(t)' : '_plot_u_t',
                      'u_cut(t)' : '_plot_u_t_cut',
                      'F(w)_cut' : '_plot_F_w_cut',
                      'F(w_el_pred)' : '_plot_F_w_el_pred',
                      'M(t), M_II(t), M_0(t)' : '_plot_M_MN_MF_t',
                      'N(M), N(M_II), N(M_0)' : '_plot_N_M_MN_MF',
                                             }

    default_plot_template = 'N(t), F(t)'

    def _plot_N_F_t(self, axes):
        '''Normal force and bending force versus time
        '''
        # ax1 = axes
        axes.plot(self.t_asc, self.N_asc, color='blue', label='N')
        axes.xaxis.tick_bottom()
        axes.set_ylim(0, 50)
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('N / F[kN]')
        # ax2 = axes.twinx()
        axes.plot(self.t_asc, self.F_asc, color='red', label='F')
        # ax2.set_ylim(0, 3)
        # ax2.set_ylabel('F [kN]')
        axes.legend(loc=2)

    def _plot_N_F_t_cut(self, axes):
        '''Normal force_cut and bending force_cut versus time_cut
        '''
        # ax1 = axes
        axes.plot(self.t_cut_asc, self.N_cut_asc, color='blue', label='N')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('N [kN]')
        axes.set_ylim(0, 50)
        # ax2 = axes.twinx()
        axes.plot(self.t_cut_asc, self.F_cut_asc, color='red', label='F')
        # ax2.set_ylim(0, 3)
        # ax2.set_ylabel('F [kN]')
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
        axes.plot(self.t_cut_asc, self.M, color='green', label='M')
        axes.plot(self.t_cut_asc, self.MN, color='blue', label='M_II')
        axes.plot(self.t_cut_asc, self.MF, color='red', label='M_0')
        axes.xaxis.tick_bottom()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('M / M_II / M_0 [kNm]')
        axes.set_ylim(-0.25, 0.45)
        axes.legend(loc=2)

    def _plot_N_M_MN_MF(self, axes):
        '''Normal force versus moment
        '''
        axes.plot(self.MF, self.N_cut_asc, color='red', label='M_0')
        axes.plot(self.MN, self.N_cut_asc, color='blue', label='M_II')
        axes.plot(self.M, self.N_cut_asc, color='green', label='M')
        axes.xaxis.tick_top()
        axes.set_xlabel('M / M_II / M_0 [kNm]')
        axes.xaxis.set_label_position('top')
        axes.set_ylabel('N [kN]')
        axes.set_ylim(24, 0)
        axes.set_xlim(-0.25, 0.45)
        axes.legend(loc=2)

    #===========================================================================
    # 2D-ARAMIS PROCESSING
    #===========================================================================

    aramis_resolution_key = Str('Xf15s3-Yf15s3')
    '''Specification of the resolution of the measured aramis field
    '''

    start_t_aramis = Float(5.0)
    '''Start time of aramis measurement.
    '''

    delta_t_aramis = Float(5.0)
    '''Delta between aramis snapshots.
    '''

    n_steps_aramis = Property(Int)
    def _get_n_steps_aramis(self):
        return self.aramis_info.number_of_steps

    t_aramis = Property(Array('float'), depends_on='data_file, start_t, delta_t, aramis_resolution_key')
    @cached_property
    def _get_t_aramis(self):
        start_t = self.start_t_aramis
        delta_t = self.delta_t_aramis
        n_steps = self.n_steps_aramis
        t_max = self.t[self.w_cut_idx]
        # print 'n_steps', n_steps
        # print 't_max', t_max
        t_aramis_full_range = np.linspace(start_t, n_steps * delta_t, n_steps)
        # print 't-aramis_full_range', t_aramis_full_range
        # print 't_aramis_full_range[t_aramis_full_range < t_max]', t_aramis_full_range[t_aramis_full_range < t_max]
        return t_aramis_full_range[t_aramis_full_range < t_max]

    # TO-DO: def some method to cut the time at a variable position, depending on strains at the end
    # cut the last values from t_aramis,because of very high strains at the end due to failure
    t_aramis_cut = Property(Array('float'), depends_on='data_file, start_t, delta_t, aramis_resolution_key')
    @cached_property
    def _get_t_aramis_cut(self):
        return self.t_aramis[ :-2]

    n_steps = Property
    @cached_property
    def _get_n_steps(self):
        'number of time steps in aramis after limiting the steps with t_max'
        # print 'n_steps', len(self.t_aramis)
        return len(self.t_aramis_cut)

    aramis_info = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_aramis_info(self):
        af = self.get_cached_aramis_file(self.aramis_resolution_key)
        if af == None:
            return None
        return AramisInfo(data_dir=af)

    N_t_aramis = Property
    @cached_property
    def _get_N_t_aramis(self):
        'normal force interpolated to the time steps of aramis'
        # print 'np.interp(self.t_aramis, self.t, self.N)', np.interp(self.t_aramis, self.t, self.N)
        return np.interp(self.t_aramis_cut, self.t, self.N)

    F_t_aramis = Property
    @cached_property
    def _get_F_t_aramis(self):
        'bending force interpolated to the time steps of aramis'
        return np.interp(self.t_aramis_cut, self.t, self.F)

    w_t_aramis = Property
    @cached_property
    def _get_w_t_aramis(self):
        'displacement (with eliminated predeformation) interpolated to the time steps of aramis'
        return np.interp(self.t_aramis_cut, self.t_cut_asc, self.w_el_pred)

    M_t_aramis = Property
    @cached_property
    def _get_M_t_aramis(self):
        'resulting moment interpolated to the time steps of aramis'
        return self.F_t_aramis * self.length / 4 - self.N_t_aramis * self.w_t_aramis / 1000

    MN_t_aramis = Property
    @cached_property
    def _get_MN_t_aramis(self):
        'moment due to normal force interpolated to the time steps of aramis'
        return self.N_t_aramis * self.w_t_aramis / 1000 * -1

    MF_t_aramis = Property
    @cached_property
    def _get_MF_t_aramis(self):
        'moment due to bending force interpolated to the time steps of aramis'
        return self.F_t_aramis * self.length / 4

    #--------------------------------------------------------------------------------
    # crack bridge strain
    #--------------------------------------------------------------------------------

    crack_filter_avg = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_crack_filter_avg(self):
        ''' method to get number and position of cracks
        '''

        ai = self.aramis_info
        if ai == None:
            return None
        ad = AramisData(aramis_info=self.aramis_info,
                        evaluated_step_idx=self.n_steps)
        absa = AramisBSA(aramis_info=self.aramis_info,
                         crack_detect_idx=self.n_steps,
                         aramis_data=ad,
                         integ_radius=10)

        print 'crack detect step', absa.crack_detect_idx
        print 'testing n_cracks', absa.number_of_cracks_avg
        # print 'ad.x_arr_undeformed', ad.x_arr_undeformed [0]
        print 'ad.x_arr_undeformed[0, absa.crack_filter_avg]', ad.x_arr_undeformed[0, absa.crack_filter_avg]
        # print 'ad.length_x_undeformed', ad.length_x_undeformed
        # print 'ad.length_y_undeformed', ad.length_y_undeformed
        # print absa.crack_filter_avg
        return ad.x_arr_undeformed[0, absa.crack_filter_avg]


    crack_bridge_strain_all = Property  # (depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_crack_bridge_strain_all(self):
        '''method to get crack bridge strain for the cracks determined by crack_filter_avg
        '''

        ai = self.aramis_info
        if ai == None:
            return None
        ad = AramisData(aramis_info=self.aramis_info)
        absa = AramisBSA(aramis_info=self.aramis_info,
                         aramis_data=ad,
                         integ_radius=10)

        # get the indices oft the cracks
        b = ad.x_arr_undeformed [0]
        c = self.crack_filter_avg
        crack_idx_list = []
        for c in c:
            crack_idx = np.where(b == c)
            crack_idx_list.append(crack_idx[0])

        crack_idx_arr1 = np.array(crack_idx_list)
        # print 'crack_idx_arr1', crack_idx_arr1
        crack_idx_arr2 = np.rollaxis(crack_idx_arr1, 1, 0)
        # print 'crack_idx_arr2', crack_idx_arr2
        crack_idx_arr = crack_idx_arr2[0]
        # print 'crack_idx_arr', crack_idx_arr

        # get the indices of the middle between two cracks
        crack_mid_idx_arr = (crack_idx_arr[0:-1] + crack_idx_arr[1:]) / 2
        # print 'crack_mid_idx_arr', crack_mid_idx_arr
        i_max = len(crack_mid_idx_arr)
        # print 'i_max', i_max

        if i_max <= 1:
            # print ' less than three cracks'
            return None

        # get crack bridge strain
        for step, t in enumerate(self.t_aramis_cut):
            ad.evaluated_step_idx = step
            eps_list = []
            for i in range(0, i_max - 1, 1):
                eps_crack_bridge_i = np.mean(absa.d_ux_arr[:, crack_mid_idx_arr[i]:crack_mid_idx_arr[i + 1]], axis=1)
                eps_list.append(eps_crack_bridge_i)
            # print 'np.array', np.array(eps_list, dtype='f')
        return np.array(eps_list, dtype='f')


    idx_failure_crack = Property  # (depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_idx_failure_crack(self):
        '''method to get the index of the failure crack and the index of the corresponding
           middle to the neighbor cracks
        '''

        ai = self.aramis_info
        if ai == None:
            return None
        ad = AramisData(aramis_info=self.aramis_info)

        # identify failure crack
        a = self.crack_bridge_strain_all  # only strain in the last step

        if a == None:
            return None

        max_strain = a.max()
        # print 'max_strain', max_strain
        max_idx = np.where(a == max_strain)[0]
        # print 'max_idx', max_idx

        # get the indices oft all cracks
        b = ad.x_arr_undeformed [0]
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
        # print 'idx_failure_crack', idx_failure_crack
        idx_border1 = crack_mid_idx_arr[max_idx]
        # print 'idx_border1 ', idx_border1
        idx_border2 = crack_mid_idx_arr[max_idx + 1]
        # print 'idx_border2 ', idx_border2

        return idx_failure_crack, idx_border1, idx_border2


    #--------------------------------------------------------------------------------
    # get tensile strain in first reinforcement layer
    #--------------------------------------------------------------------------------

    h_re1_6_threshold = Float(2.86, auto_set=False, enter_set=True)
    '''Threshold for position of first reinforcement layer (6 layers).
    '''
    h_re1_4_threshold = Float(4.0, auto_set=False, enter_set=True)
    '''Threshold for position of first reinforcement layer (4 layers).
    '''

    meas_field = Property
    @cached_property
    def _get_meas_field(self):
        '''method to get length and height of measuring field (calculated from the center of the facets)
        '''
        ai = self.aramis_info
        if ai == None:
            return None

        AUI = AramisUI(aramis_info=self.aramis_info)
        AUI.aramis_data.evaluated_step_idx = self.n_steps - 2

        x = AUI.aramis_data.x_arr_undeformed
        y = AUI.aramis_data.y_arr_undeformed

        l_mf1 = x[0, -1] - x[0, 0]
        l_mf2 = x[-1, -1] - x[-1, 0]

        h_mf1 = y[0, 0] - y[-1, 0]
        h_mf2 = y[0, -1] - y[-1, -1]

        l_mf = (l_mf1 + l_mf2) / 2
        h_mf = (h_mf1 + h_mf2) / 2
        # print 'meas_field', l_mf, h_mf
        return l_mf, h_mf

    h_dis = Property
    @cached_property
    def _get_h_dis(self):
        ''' method to get distance between specimen edge and first / last facet node in y-direction
        '''
        ai = self.aramis_info
        if ai == None:
            return None

        h_dis = (20 - self.meas_field[1]) / 2
        # print 'h_dis', h_dis
        return h_dis

    pos_fa = Property
    @cached_property
    def _get_pos_fa(self):
        # method to get position of first and last facet-node in y-direction

        ai = self.aramis_info
        if ai == None:
            return None

        # position of first node in y-direction
        pos_no_f = 20.0 - self.h_dis
        # position of last node in y-direction
        pos_no_l = self.h_dis

        return pos_no_f, pos_no_l

    eps1_t_aramis = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_eps1_t_aramis(self):
        '''method to get tensile strain in first reinforcement layer
        '''
        ai = self.aramis_info
        if ai == None:
            return None

        ad = AramisData(aramis_info=self.aramis_info)
        absa = AramisBSA(aramis_info=self.aramis_info,
                         aramis_data=ad,
                         integ_radius=10)

        # get number of facet nodes in y-direction
        n_fa = absa.d_ux_arr.shape[0]
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
            idx_6a = round (idx_6a)
        idx_6b = idx_6a + 1

        if pos_re1_4 < 0:
            idx_4a = 0
        else:
            idx_4a = pos_re1_4 / dis_fa
            idx_4a = round (idx_4a)
        idx_4b = idx_4a + 1

        a = self.crack_bridge_strain_all
        eps_t_list = []

        # get eps array
        for step, t in enumerate(self.t_aramis_cut):
            ad.evaluated_step_idx = step
            if a == None:
                mid_idx = absa.d_ux_arr.shape[1] / 2
                eps_range = 3
                eps = np.mean(absa.d_ux_arr[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
                # print 'strain in the middle of the measuring field'
            else:
                idx_border1 = self.idx_failure_crack[1]
                idx_border2 = self.idx_failure_crack[2]
                eps = np.mean(absa.d_ux_arr[:, idx_border1:idx_border2], axis=1)
                # print 'strain in the failure crack'

            # if 6layers
            '''x1 = pos_re1_6 - idx_6a * dis_fa
            print 'x1', x1
            x_re1 = (dis_fa - x1) * (eps[idx_6a] - eps[idx_6b]) / dis_fa
            print 'x_re1', x_re1
            eps_re1 = eps[idx_6b] + x_re1
            print 'eps_re1', eps_re1
            eps_t_list.append(eps_re1)
            print 'eps_t_list', np.array(eps_t_list, dtype='f')
        return np.array(eps_t_list, dtype='f')'''

            # if 4 layers
            x1 = pos_re1_4 - idx_4a * dis_fa
            # print 'x1', x1
            x_re1 = (dis_fa - x1) * (eps[idx_4a] - eps[idx_4b]) / dis_fa
            # print 'x_re1', x_re1
            eps_re1 = eps[idx_4b] + x_re1
            # print 'eps_re1', eps_re1
            eps_t_list.append(eps_re1)
            # print 'eps_t_list', np.array(eps_t_list, dtype='f')
        return np.array(eps_t_list, dtype='f')

    #--------------------------------------------------------------------------------
    # get max and min strain of cross section
    #--------------------------------------------------------------------------------

    eps_t_aramis = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_eps_t_aramis(self):

        ai = self.aramis_info
        if ai == None:
            return None

        ad = AramisData(aramis_info=self.aramis_info)
        absa = AramisBSA(aramis_info=self.aramis_info,
                         aramis_data=ad,
                         integ_radius=10)

        a = self.crack_bridge_strain_all
        n_fa = absa.d_ux_arr.shape[0]
        h = np.linspace(self.pos_fa[0], self.pos_fa[1], num=n_fa)
        eps_t_list = []
        eps_c_list = []

        # get eps
        for step, t in enumerate(self.t_aramis_cut):
            ad.evaluated_step_idx = step
            if a == None:
                mid_idx = absa.d_ux_arr.shape[1] / 2
                eps_range = 3
                eps = np.mean(absa.d_ux_arr[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
            else:
                idx_border1 = self.idx_failure_crack[1]
                idx_border2 = self.idx_failure_crack[2]
                eps = np.mean(absa.d_ux_arr[:, idx_border1:idx_border2], axis=1)

            # extrapolate eps for the specimen edges
            x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
            eps_ed_up = x + eps[-1]
            eps_ed_lo = eps[0] - x
            eps_to1 = np.append(eps, eps_ed_lo)
            eps_to2 = np.append(eps_ed_up, eps_to1)

            eps_t_list.append(np.max(eps_to2))
            eps_c_list.append(np.min(eps_to2))

            # print 'np.array', np.array(eps_t_list, dtype='f')
        return np.array(eps_t_list, dtype='f'), np.array(eps_c_list, dtype='f')

    #--------------------------------------------------------------------------------
    # get strain(N) and strain(M)
    #--------------------------------------------------------------------------------

    F_beg_threshold = Float(0.11, auto_set=False, enter_set=True)
    '''Threshold to find the beginning of raising F
    '''

    t_F_beg = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_t_F_beg(self):
        # get the time of the the beginning of raising F

        ai = self.aramis_info
        if ai == None:
            return None

        F = self.F_t_aramis
        # print 'F', F
        N = self.N_t_aramis
        # print 'N', N
        t = self.t_aramis_cut
        # print 't', t

        F_idx = np.where(F > self.F_beg_threshold)[0]
        # print 'F_idx', F_idx

        if len(F_idx) == 0:
            # print 't_F_beg', []
            return []
        elif max(N) <= 2:
            # print 't_F_beg', t[0]
            return t[0]
        else:
            # print 't_F_beg', t[F_idx[0]]
            return t[F_idx[0]]

    t_N_arr = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_t_N_arr(self):
        # get the time where N rises

        ai = self.aramis_info
        if ai == None:
            return None

        t = self.t_aramis_cut

        if self.t_F_beg == 0:
            t_N = []
        elif self.t_F_beg == []:
            t_N = t
        else:
            t_N_end = self.t_F_beg - 5
            t_N = t [t <= t_N_end]

        # print 't_N_arr', t_N
        return t_N

    t_F_arr = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_t_F_arr(self):
        # get the time where F rises

        ai = self.aramis_info
        if ai == None:
            return None

        t = self.t_aramis_cut

        if self.t_F_beg == 0:
            t_F = t
        elif self.t_F_beg == []:
            t_F = []
        else:
            t_F = t[t >= self.t_F_beg]

        # print 't_F_arr', t_F
        return t_F


    eps_N = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_eps_N(self):
        # get the strain corresponding to N

        ai = self.aramis_info
        if ai == None:
            return None

        ad = AramisData(aramis_info=self.aramis_info)
        absa = AramisBSA(aramis_info=self.aramis_info,
                         aramis_data=ad,
                         integ_radius=10)

        a = self.crack_bridge_strain_all

        t = self.t_aramis_cut
        t_N = self.t_N_arr

        if len(self.t_N_arr) == 0:
            N_end_idx = []
        else:
            N_end_idx = np.where(t == self.t_N_arr[-1]) [0]

        eps_N_list = []

        if t_N != []:
            for step, t in enumerate(t_N):
                ad.evaluated_step_idx = step
                if a == None:
                    mid_idx = absa.d_ux_arr.shape[1] / 2
                    eps_range = 3
                    eps = np.mean(absa.d_ux_arr[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)

                else:
                    idx_border1 = self.idx_failure_crack[1]
                    idx_border2 = self.idx_failure_crack[2]
                    eps = np.mean(absa.d_ux_arr[:, idx_border1:idx_border2], axis=1)

                eps_N_list.append(np.mean(eps))

                # print 'eps_N_list', np.array(eps_N_list, dtype='f')
            return np.array(eps_N_list, dtype='f')


    N_t_N = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_N_t_N(self):
        # get N in the area of t_N

        t = self.t_aramis_cut

        if len(self.t_N_arr) == 0:
            # print 'self.N_t_aramis[0:t_N_idx]', []
            return []
        else:
            N_end_idx = np.where(t == self.t_N_arr[-1]) [0]
            # print 'N_end_idx' , N_end_idx
            # print 'self.N_t_aramis[0:t_N_idx+1]', self.N_t_aramis[0:N_end_idx + 1]
            return self.N_t_aramis[0:N_end_idx + 1]

    eps_M = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_eps_M(self):
        # get the strain corresponding to M

        ai = self.aramis_info
        if ai == None:
            return None

        ad = AramisData(aramis_info=self.aramis_info)
        absa = AramisBSA(aramis_info=self.aramis_info,
                         aramis_data=ad,
                         integ_radius=10)

        a = self.crack_bridge_strain_all

        t = self.t_aramis_cut
        t_F = self.t_F_arr

        if len(self.t_F_arr) == 0:
            F_beg_idx = []
        else:
            F_beg_idx = np.where(t == self.t_F_arr[0]) [0]

        eps_M_max_list = []
        eps_M_min_list = []

        step_1 = F_beg_idx
        max_step = self.n_steps

        if t_F != []:
            for step in range(step_1, max_step, 1):
                ad.evaluated_step_idx = step
                if a == None:
                    mid_idx = absa.d_ux_arr.shape[1] / 2
                    eps_range = 3
                    eps = np.mean(absa.d_ux_arr[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)

                else:
                    idx_border1 = self.idx_failure_crack[1]
                    idx_border2 = self.idx_failure_crack[2]
                    eps = np.mean(absa.d_ux_arr[:, idx_border1:idx_border2], axis=1)

                eps_M_max_list.append(np.max(eps))
                eps_M_min_list.append(np.min(eps))
                # print 'eps_M_max_list', np.array(eps_M_max_list, dtype='f')
                # print 'eps_M_min_list', np.array(eps_M_min_list, dtype='f')

            return np.array(eps_M_max_list, dtype='f'), np.array(eps_M_min_list, dtype='f')

    F_t_F = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_F_t_F(self):
        # get F in the area of t_F

        t = self.t_aramis_cut

        if self.t_F_arr == []:
            # print 'self.F_t_aramis[F_beg_idx:]', []
            return []
        else:
            F_beg_idx = np.where(t == self.t_F_arr[0]) [0]
            # print 'F_beg_idx', F_beg_idx
            # print 'self.F_t_aramis[F_beg_idx:]', self.F_t_aramis[F_beg_idx:]
            return self.F_t_aramis[F_beg_idx:]

    #--------------------------------------------------------------------------------
    # get height of compression zone
    #--------------------------------------------------------------------------------

    x_t_aramis = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_x_t_aramis(self):

        ai = self.aramis_info
        if ai == None:
            return None

        ad = AramisData(aramis_info=self.aramis_info)
        absa = AramisBSA(aramis_info=self.aramis_info,
                         aramis_data=ad,
                         integ_radius=10)

        a = self.crack_bridge_strain_all

        n_fa = absa.d_ux_arr.shape[0]
        dis_fa = self.meas_field[1] / (n_fa - 1)

        h = np.linspace(self.pos_fa[0], self.pos_fa[1], num=n_fa)
        h_1 = np.append(h, 0)
        h_2 = np.append(20, h_1)

        x_list = []

        # get eps
        for step, t in enumerate(self.t_aramis_cut):
            ad.evaluated_step_idx = step
            if a == None:
                mid_idx = absa.d_ux_arr.shape[1] / 2
                eps_range = 3
                eps = np.mean(absa.d_ux_arr[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
            else:
                idx_border1 = self.idx_failure_crack[1]
                idx_border2 = self.idx_failure_crack[2]
                eps = np.mean(absa.d_ux_arr[:, idx_border1:idx_border2], axis=1)

            # extrapolate eps for the specimen edges
            x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
            eps_ed_up = x + eps[-1]
            eps_ed_lo = eps[0] - x
            eps_to1 = np.append(eps, eps_ed_lo)
            eps_to2 = np.append(eps_ed_up, eps_to1)

            print 'eps_to2', eps_to2

            idx_neg = np.array(np.where(eps_to2 < 0) [0])
            print 'idx_neg', idx_neg
            d = idx_neg.shape [0]
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

                x_a = (dis_fa * abs(eps_to2[idx_1])) / (eps_to2[idx_2] - eps[idx_1])
                # print 'x_a', x_a
                x = h_2[idx_1] + x_a
                # print 'x', x

            x_list.append(x)
            print 'x', np.array(x_list, dtype='f')

        return np.array(x_list, dtype='f')

    #--------------------------------------------------------------------------------
    # get curvature of the specimen
    #--------------------------------------------------------------------------------

    cu_t_aramis = Property(depends_on='data_file,aramis_resolution_key')
    @cached_property
    def _get_cu_t_aramis(self):

        ai = self.aramis_info
        if ai == None:
            return None
        z = 20 - self.x_t_aramis

        cu = self.eps1_t_aramis[0] / z
        return cu

    #---------------------------------
    # view
    #---------------------------------

    traits_view = View(VSplit(
                         HSplit(Group(
                                  Item('width', format_str="%.3f"),
                                  Item('length', format_str="%.3f"),
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
                               Item('E_c', visible_when='derived_data_available',
                                                style='readonly', show_label=True , format_str="%.0f"),
                               Item('N_max', visible_when='derived_data_available',
                                                style='readonly', emphasized=True , format_str="%.2f"),
                               Item('F_max', visible_when='derived_data_available',
                                                style='readonly', emphasized=True , format_str="%.2f"),
                               Item('w_max', visible_when='derived_data_available',
                                                style='readonly', emphasized=True , format_str="%.2f"),
                               Item('u_max', visible_when='derived_data_available',
                                                style='readonly', emphasized=True , format_str="%.2f"),
                               Item('M_max', visible_when='derived_data_available',
                                                style='readonly', emphasized=True , format_str="%.3f"),
                               Item('MN_max', visible_when='derived_data_available',
                                                style='readonly', emphasized=True , format_str="%.3f"),
                               Item('MF_max', visible_when='derived_data_available',
                                                style='readonly', emphasized=True , format_str="%.3f"),
                               Item('w_pred', visible_when='derived_data_available',
                                                style='readonly', emphasized=True , format_str="%.3f"),
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

#    ExpTTDB.add_class_trait('production_date', Date(input=True, table_field=True,))
#    for inst in ExpTTDB.db.inst_list:
#        print inst.key
#        print inst.add_trait('production_date', Date('14/9/2011', input=True, table_field=True,))
#        print inst.production_date
#        inst.save()

    import pylab as p


    from matresdev.db.simdb import SimDB
    from matresdev.db.exdb import ExRunView
    simdb = SimDB()
    import os

    ex_path = os.path.join(simdb.exdata_dir,
                           'bending_tensile_test',
                           '2014-06-12_BTT-4c-2cm-0-TU_MxN2',
                           'BTT-4c-2cm-TU-0-V02_MxN2.DAT')

    test_file = os.path.join(simdb.exdata_dir,
                           'bending_tensile_test',
                           '2014-06-12_BTT-4c-2cm-0-TU_MxN2',
                           'BTT-4c-2cm-TU-0-V02_MxN2.DAT')

    doe_reader = ExRunView(data_file=ex_path)
    doe_reader.configure_traits()

    # ExpBTTDB.db.configure_traits()
    # to see all experiments in one picture
