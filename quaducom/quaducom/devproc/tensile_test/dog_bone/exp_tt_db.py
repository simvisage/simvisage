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
# Created on Feb 15, 2010 by: rch, ascholzen

# @todo - construct the class for fabric layout calculating the
#         cs-area of the reinforcement.
#       - instead of processed array - construct the array traits accessible
#         with the name of the measured channels
#       - reread the pickle file without processing the data (take care to
#         reestablish
#         the link from the ex_type to the ex_run
#       - define the exdb_browser showing the inputs and outputs in a survey
#       - define the ExTreatment class with cumulative evaluation of the
#         response values.
#
#

from traits.api import \
    Int, Float, \
    on_trait_change, Instance, \
    Array, Property, cached_property, \
    Bool, Event, implements, \
    DelegatesTo

from numpy import \
    array, where, argmax

import numpy as np

from traitsui.api \
    import View, Item, HSplit, Group, VSplit

from matresdev.db.exdb.ex_type import ExType
from matresdev.db.exdb.i_ex_type import IExType

from mathkit.array.smoothing import smooth

from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp

from matresdev.db.matdb.trc.fabric_layout \
    import FabricLayOut

from matresdev.db.matdb.trc.concrete_mixture \
    import ConcreteMixture

from matresdev.db.matdb.trc.composite_cross_section \
    import CompositeCrossSection, plain_concrete

from matresdev.db.exdb.ex_run_table import ExRunClassExt


class ExpTTDB(ExType):

    '''Experiment: Tensile Test with dog bond shape
    '''
#    label = Str('dog bone tensile test')

    implements(IExType)

    # --------------------------------------------------------------------
    # register a change of the traits with metadata 'input'
    # --------------------------------------------------------------------

    input_change = Event

    @on_trait_change('+input, ccs.input_change')
    def _set_input_change(self):
        print '*** raising input change in CTT'
        self.input_change = True

    # -------------------------------------------------------------------------
    # specify inputs:
    # -------------------------------------------------------------------------

    width = Float(0.100, unit='m', input=True, table_field=True,
                  auto_set=False, enter_set=True)
    '''Width of the specimen cross section.
    '''
    gauge_length = Float(0.250, unit='m', input=True, table_field=True,
                         auto_set=False, enter_set=True)
    '''Gauge length for measuring the elongation of the specimen.
    '''

    age = Int(29, unit='d', input=True, table_field=True,
              auto_set=False, enter_set=True)
    '''Age of the concrete at the time of testing.
    '''

    loading_rate = Float(2.0, unit='mm/min', input=True, table_field=True,
                         auto_set=False, enter_set=True)
    '''Applied loading rate.
    '''

    # --------------------------------------------------------------------------
    # composite cross section
    # --------------------------------------------------------------------------

    ccs = Instance(CompositeCrossSection)
    '''Link to the composite cross section object containing
    defining the layout and lay-up of the reinforcement fabrics.
    '''

    def _ccs_default(self):
        '''default settings correspond to
        setup '9u_MAG-07-03_PZ-0708-1'
        '''
        print 'ccs default used'
#        fabric_layout_key = 'MAG-07-03'
#        fabric_layout_key = '2D-02-06a'
#        fabric_layout_key2 = 'C-Grid-C50'
#        fabric_layout_key1 = '2D-14-10'
#        fabric_layout_key = '2D-14-10'
#        fabric_layout_key = '2D-18-10'
#        fabric_layout_key = '2D-04-11'
        fabric_layout_key = '2D-05-11'
#        fabric_layout_key = 'Grid-600'
#        fabric_layout_key = '2D-15-10'
#        concrete_mixture_key = 'PZ-0708-1'
        concrete_mixture_key = 'barrelshell'
#        concrete_mixture_key = 'sto-100'
#        concrete_mixture_key = 'FIL-10-09'
        orientation_fn_key = 'all0'
#        orientation_fn_key = 'all90'
#        orientation_fn_key = '90_0'
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

    # --------------------------------------------------------------------------
    # Indicate whether the test is suitable and prepared for
    # calibration.
    # --------------------------------------------------------------------------
    ready_for_calibration = Property(Bool)
    '''Indicator if the test can be used for the calibration
    of the smeared composite-microplane-damage model.
    The composite cross section can be only used if the cross section has
    a regular layout.
    '''

    def _get_ready_for_calibration(self):
        # return False by default
        # the subclasses shall overload this
        # and define the rules
        return self.ccs.is_regular

    # --------------------------------------------------------------------------
    # Get properties of the composite
    # --------------------------------------------------------------------------

    E_c = Property(
        Float, unit='MPa', depends_on='input_change', table_field=True)
    '''E-modulus of the composite at the time of testing.
    '''

    @cached_property
    def _get_E_c(self):
        return self.ccs.get_E_c_time(self.age)

    E_m = Property(
        Float, unit='MPa', depends_on='input_change', table_field=True)
    '''E-modulus of the concrete matrix at the time of testing.
    '''

    @cached_property
    def _get_E_m(self):
        return self.ccs.get_E_m_time(self.age)

    A_c = Property(Float, unit='m^2', depends_on='input_change')
    '''Cross-sectional-area of the composite.
    '''

    @cached_property
    def _get_A_c(self):
        return self.width * self.ccs.thickness

    A_tex = Property(Float, unit='mm^2', depends_on='input_change')
    '''Total cross-sectional-area of the textile reinforcement.
    '''

    @cached_property
    def _get_A_tex(self):
        return self.ccs.a_tex * self.width

    E_c28 = DelegatesTo('ccs', listenable=False)
    '''E-modulus of the composite after 28 days.
    '''

    rho_c = DelegatesTo('ccs', listenable=False)
    '''Reinforcement ratio of the composite.
    '''

    # -------------------------------------------------------------------------
    # define processing
    # -------------------------------------------------------------------------

    def process_source_data(self):
        '''Extend the default data processing with
        the handling of gauges placed at the front and back side of
        the specimen. If necessary modify the assigned data, i. e. change
        the sign or specify an offset for the specific test setup.
        @todo: make the distinction in the subclasses.
        '''
        super(ExpTTDB, self).process_source_data()

        # NOTE: the small tensile tests (INSTRON) with width = 0.10 m have
        # only 3 displacement gauges
        #
        if hasattr(self, "W10_re") and hasattr(self, "W10_li") \
                and hasattr(self, "W10_vo"):
            self.W10_re -= self.W10_re[0]
            self.W10_re *= -1
            self.W10_li -= self.W10_li[0]
            self.W10_li *= -1
            self.W10_vo -= self.W10_vo[0]
            self.W10_vo *= -1

        if hasattr(self, "W10_vli"):
            print 'change_varname'
            self.WA_VL = self.W10_vli
        if hasattr(self, "W10_vre"):
            self.WA_VR = self.W10_vre
        if hasattr(self, "W10_hli"):
            self.WA_HL = self.W10_hli
        if hasattr(self, "W20_hre"):
            self.WA_HR = self.W20_hre

        # NOTE: the large tensile tests (PSB1000) with width = 0.14 m have
        # 4 displacement gauges
        #
        if hasattr(self, "WA_VL") and hasattr(self, "WA_VR") \
                and hasattr(self, "WA_HL") and hasattr(self, "WA_HR"):
            self.WA_VL -= self.WA_VL[0]
            self.WA_VL *= -1
            self.WA_VR -= self.WA_VR[0]
            self.WA_VR *= -1
            self.WA_HL -= self.WA_HL[0]
            self.WA_HL *= -1
            self.WA_HR -= self.WA_HR[0]
            self.WA_HR *= -1

    eps = Property(Array('float_'), output=True,
                   depends_on='input_change')
    '''Strains calculated from the input gauge displacements.
    '''

    @cached_property
    def _get_eps(self):

        if hasattr(self, "W10_re") and hasattr(self, "W10_li") and \
                hasattr(self, "W10_vo"):
            W10_li = np.copy(self.W10_li)
            W10_re = np.copy(self.W10_re)
            W10_vo = np.copy(self.W10_vo)

            # get the minimum value of the displacement gauges
            # used to reset the displacement gauges if they do not start at
            # zero
            min_W10_li = np.min(W10_li[:10])
            min_W10_re = np.min(W10_re[:10])
            min_W10_vo = np.min(W10_vo[:10])

            # reset displacement gauges
            #
            W10_li -= min_W10_li
            W10_re -= min_W10_re
            W10_vo -= min_W10_vo

            # measured strains
            eps_li = W10_li / (self.gauge_length * 1000.)  # [mm/mm]
            eps_re = W10_re / (self.gauge_length * 1000.)
            eps_vo = W10_vo / (self.gauge_length * 1000.)

            # NOTE: if only 2 displacement gauges are used instead of 3 (only 'WA_re' for front and 'WA_li' for back)
            # below the average is performed as = 0.5*( 0.5*(W10_re + W10_li) + W10_vo)
            #
            if np.average(eps_re) < 0.0001:
                print "displacement gauge 'WA_re' has not been used. Use value of 'WA_li' instead"
                eps_re = eps_li
            if np.average(eps_li) < 0.0001:
                print "displacement gauge 'WA_li' has not been used. Use value of 'WA_re' instead"
                eps_li = eps_re
            if np.average(eps_vo) < 0.0001:
                print "displacement gauge 'WA_vo' has not been used. Use average value of 'WA_li' and 'WA_re' instead"
                eps_vo = (eps_li + eps_re) / 2.

            # average strains
            #
            eps_m = ((eps_li + eps_re) / 2. + eps_vo) / 2.

        if hasattr(self, "WA_VL") and hasattr(self, "WA_VR") and hasattr(self, "WA_HL") and hasattr(self, "WA_HR"):
            WA_VL = np.copy(self.WA_VL)
            WA_VR = np.copy(self.WA_VR)
            WA_HL = np.copy(self.WA_HL)
            WA_HR = np.copy(self.WA_HR)

            # get the minimum value of the displacement gauges
            # used to reset the displacement gauges if they do not start at
            # zero
            min_WA_VL = np.min(WA_VL[:10])
            min_WA_VR = np.min(WA_VR[:10])
            min_WA_HL = np.min(WA_HL[:10])
            min_WA_HR = np.min(WA_HR[:10])

            # reset displacement gauges
            #
            WA_VL -= min_WA_VL
            WA_VR -= min_WA_VR
            WA_HL -= min_WA_HL
            WA_HR -= min_WA_HR

            # measured strains
            #
            eps_V = (self.WA_VL + self.WA_VR) / 2. / \
                (self.gauge_length * 1000.)  # [mm/mm]
            eps_H = (self.WA_HL + self.WA_HR) / \
                2. / (self.gauge_length * 1000.)

            # average strains
            #
            eps_m = (eps_V + eps_H) / 2.

        return eps_m

    sig_c = Property(Array('float_'), output=True,
                     depends_on='input_change')
    '''Composite stress obtained as force related to the total area.
    '''
    @cached_property
    def _get_sig_c(self):
        print 'CALCULATING COMPOSITE STRESS'
        # measured force:
        force = self.Kraft  # [kN]
        # cross sectional area of the concrete [m^2]:
        A_c = self.A_c
        # calculated stress:
        sig_c = (force / 1000.) / A_c  # [MPa]
        return sig_c

    sig_tex = Property(Array('float_'),
                       output=True, depends_on='input_change')
    '''Stress in the textile fabrics cross section.
    '''
    @cached_property
    def _get_sig_tex(self):
        # measured force:
        force = self.Kraft  # [kN]
        # cross sectional area of the reinforcement:
        A_tex = self.A_tex
        # calculated stresses:
        sig_tex = (force * 1000.) / A_tex  # [MPa]
        return sig_tex

    max_stress_idx = Property(Int, depends_on='input_change')
    '''Maximum stress index to cut off the descending part of the curves.
    '''
    @cached_property
    def _get_max_stress_idx(self):
        return argmax(self.sig_c)

    eps_asc = Property(Array('float_'), depends_on='input_change')
    '''Strain corresponding to the ascending branch of the measured stress curve.
    '''
    @cached_property
    def _get_eps_asc(self):
        return self.eps[:self.max_stress_idx + 1]

    sig_c_asc = Property(Array('float_'), depends_on='input_change')
    '''Ascending branch of the composite stress.
    '''
    @cached_property
    def _get_sig_c_asc(self):
        return self.sig_c[:self.max_stress_idx + 1]

    sig_tex_asc = Property(Array('float_'), depends_on='input_change')
    '''Ascending branch of the textile stress
    '''
    @cached_property
    def _get_sig_tex_asc(self):
        return self.sig_tex[:self.max_stress_idx + 1]

    F_asc = Property(Array('float_'), depends_on='input_change')
    '''Ascending branch of the force
    '''
    @cached_property
    def _get_F_asc(self):
        return self.Kraft[:self.max_stress_idx + 1]

    time_asc = Property(Array('float_'), depends_on='input_change')
    '''Time corresponding to the ascending branch of the stress curve.
    '''
    @cached_property
    def _get_time_asc(self):
        return self.Bezugskanal[:self.max_stress_idx + 1]

    jump_rtol = Float(0.0001, ironing_param=True)
    '''Relative tolerance for ironing the jumps in the force response.
    '''

    F_w_ironed = Property(Array('float_'), depends_on='input_change')
    ''''Force curve with removes jumps in the stress strain curve
    due to sliding in the buttstrap clamping. The small unloading/loading
    branches are removed from the data and the smoothed curve is
    determined based on the remaining data.
    In order to smoothen out the effect of the jump the pieces of the
    F-w-curve that contain jumps in the force (=unloading/reloading) path
    are removed from the data in a range up to the double range of the jump,
    so that the smoothing does not change the shape of the F-w- curve.
    '''
    @cached_property
    def _get_F_w_ironed(self):
        F_asc = self.F_asc
        eps_asc = self.eps_asc
        time_asc = self.time_asc

        # get the differences of the force values between two adjacent
        # time steps
        #
        jump_arr = F_asc[1:] - F_asc[0:-1]

        # determine the criteria for a jump
        # based on the data range and the specified tolerances:
        #
        jump_rtol = self.jump_rtol
        jump_crit = jump_rtol * F_asc[-1]

        # get the indices of the measurement data at which a
        # force jump exceeds (last step before the jump) the defined tolerance
        # criteria
        # i.e. negative jump that exceeds the defined tolerance magnitude
        #
        jump_idx_arr = np.where(jump_arr < -jump_crit)[0] - 1

        # index of the measurement data where the force reaches
        # the same magnitude before the sudden value drop due to the jump
        #
        jump_idx2_arr = np.zeros_like(jump_idx_arr)

        # amount of indices between the sudden value drop of the force and
        # the reloading to the same load level; delta value indicate
        # the strain range that will be removed in order to smoothen out
        # the influence of the jump in the force curve
        #
        delta_jump_idx_arr = np.zeros_like(jump_idx_arr)

        # search at which index the force reaches its old value before
        # the jump again check that this value is index wise and time
        # wise reached after the jump occurred
        #
        for n_idx, jump_idx in enumerate(jump_idx_arr):
            delta_F = F_asc - F_asc[jump_idx]
            delta_eps = eps_asc - eps_asc[jump_idx]
            delta_t = time_asc - time_asc[jump_idx]
            bool_arr_F = delta_F > 0.
            bool_arr_eps = delta_eps > 0.
            bool_arr_t = delta_t > 0.
            bool_arr = bool_arr_F * bool_arr_eps * bool_arr_t
            jump_idx2 = np.where(bool_arr)[0][1]
            delta_jump_idx = jump_idx2 - jump_idx
            jump_idx2_arr[n_idx] = jump_idx2
            delta_jump_idx_arr[n_idx] = delta_jump_idx

        # remove jumps from the jump index when a succeeding jump still
        # lays within the influence range of an earlier jump
        # this can happen when jumps occur within the re-mounting
        # branch of the force
        #
        remove_idx = []
        for i in range(jump_idx2_arr.shape[0] - 1):
            if np.any(jump_idx2_arr[:i + 1] > jump_idx2_arr[i + 1]):
                remove_idx += [i + 1]

        jump_idx_arr = np.delete(jump_idx_arr, remove_idx)
        jump_idx2_arr = np.delete(jump_idx2_arr, remove_idx)
        delta_jump_idx_arr = np.delete(delta_jump_idx_arr, remove_idx)

        # specify the factor by with the index delta range of a jump
        # (i.e. displacement range of the jump)
        # is multiplied, i.e. up to which index the values of the
        # F-w- curve are removed
        #
        jump_smooth_fact = 2

        # remove the values of the curve within the jump and the neighboring
        # region
        #
        F_asc_ironed_list = []
        eps_asc_ironed_list = []

        jump_idx_arr_ = np.hstack(
            [np.array([0.]), jump_idx_arr, np.array([self.max_stress_idx])])
        delta_jump_idx_arr_ = np.hstack(
            [np.array([0]), delta_jump_idx_arr, np.array([0])])

        for i in range(jump_idx_arr_.shape[0] - 1):
            F_asc_ironed_list += \
                [F_asc[jump_idx_arr_[i] +
                       jump_smooth_fact *
                       delta_jump_idx_arr_[i]: jump_idx_arr_[i + 1]]]
            eps_asc_ironed_list += \
                [eps_asc[jump_idx_arr_[i] +
                         jump_smooth_fact *
                         delta_jump_idx_arr_[i]: jump_idx_arr_[i + 1]]]

        # remove the values of the curve within the jump
        # and the neighboring region
        #
        F_asc_ironed = np.hstack(F_asc_ironed_list)
        eps_asc_ironed = np.hstack(eps_asc_ironed_list)
        F_asc_ironed *= 0.001  # convert units from [kN] to [MN]
        return F_asc_ironed, eps_asc_ironed

    sig_c_ironed = Property(Float, depends_on='input_change', unit='MPa')
#                                      ,output=False, table_field=False,
    '''Stress in composite - ironed
    '''
    @cached_property
    def _get_sig_c_ironed(self):
        F_asc_ironed = self.F_w_ironed[0]
        sig_c_asc_ironed = F_asc_ironed / self.A_c
#        sig_c_asc_smoothed_ironed = smooth(sig_c_asc_ironed, self.n_points, 'flat')
        return sig_c_asc_ironed

    sig_c_interpolated = Property(Float, depends_on='input_change')
#                                      ,output=False, table_field=False, unit='MPa')
    '''Stress in composite - ironed
    '''
    @cached_property
    def _get_sig_c_interpolated(self):
        sig_c_ironed = np.copy(self.sig_c_ironed)
        sig_c_interpolated = np.hstack([0., sig_c_ironed])
        return sig_c_interpolated

    sig_c_interpolated_smoothed = Property(Float, depends_on='input_change')

    @cached_property
    def _get_sig_c_interpolated_smoothed(self):
        return smooth(self.sig_c_interpolated, 15, 'flat')

    eps_c_interpolated = Property(Float, depends_on='input_change')
#                                      ,output=False, table_field=False, unit='MPa')

    @cached_property
    def _get_eps_c_interpolated(self):
        eps_c_interpolated = np.copy(self.eps_ironed)
        K_I = self.E_c  # depending of the testing age
        offset_eps_c = self.sig_c_ironed[0] / K_I
        # remove initial strain value and shift starting point of measured
        # curve into the analytical stiffness
        eps_c_interpolated -= eps_c_interpolated[0]
        eps_c_interpolated += offset_eps_c
        eps_c_interpolated = np.hstack([0., eps_c_interpolated])
        return eps_c_interpolated

    eps_c_interpolated_smoothed = Property(Float, depends_on='input_change')

    @cached_property
    def _get_eps_c_interpolated_smoothed(self):
        return smooth(self.eps_c_interpolated, 15, 'flat')

    sig_tex_ironed = Property(Float, depends_on='input_change')
#                                      ,output=False, table_field=False, unit='MPa')

    @cached_property
    def _get_sig_tex_ironed(self):
        F_asc_ironed = self.F_w_ironed[0]
        sig_tex_asc_ironed = F_asc_ironed / (self.A_tex / 1000000.)
        return sig_tex_asc_ironed

    eps_ironed = Property(Float, depends_on='input_change')
#                                    ,output=False, table_field=False, unit='-')

    @cached_property
    def _get_eps_ironed(self):
        eps_asc_ironed = self.F_w_ironed[1]
#        eps_asc_smoothed_ironed = smooth(eps_asc_ironed, self.n_points, 'flat')
        return eps_asc_ironed

    sig_tex_interpolated = Property(Float, depends_on='input_change')
#                                      ,output=False, table_field=False, unit='MPa')

    @cached_property
    def _get_sig_tex_interpolated(self):
        sig_tex_ironed = np.copy(self.sig_tex_ironed)
        sig_tex_interpolated = np.hstack([0., sig_tex_ironed])
        return sig_tex_interpolated

    eps_tex_interpolated = Property(Float, depends_on='input_change')
#                                      ,output=False, table_field=False, unit='MPa')

    @cached_property
    def _get_eps_tex_interpolated(self):
        eps_tex_interpolated = np.copy(self.eps_ironed)
        K_I = self.E_c / self.rho_c
        offset_eps_tex = self.sig_tex_ironed[0] / K_I
        # remove initial strain value and shift starting point of measured
        # curve into the analytical stiffness
        eps_tex_interpolated -= eps_tex_interpolated[0]
        eps_tex_interpolated += offset_eps_tex
        eps_tex_interpolated = np.hstack([0., eps_tex_interpolated])
        return eps_tex_interpolated

    eps_max = Property(Float, depends_on='input_change',
                       output=True, table_field=True, unit='MPa')
    '''Strain at maximum stress
    '''
    @cached_property
    def _get_eps_max(self):
        return np.max(self.eps_asc)

    sig_c_max = Property(Float, depends_on='input_change',
                         output=True, table_field=True, unit='MPa')
    '''Maximum composite stress.
    '''
    @cached_property
    def _get_sig_c_max(self):
        return self.sig_c_asc[-1]

    sig_tex_max = Property(Float, depends_on='input_change',
                           output=True, table_field=True, unit='-')
    '''Maximum textile stress.
    '''
    @cached_property
    def _get_sig_tex_max(self):
        return self.sig_tex_asc[-1]

    # -------------------------------------------------------------------------
    # Smoothing parameters
    # -------------------------------------------------------------------------
    n_smooth_window_fraction = Float(0.1)
    '''Processing parameter for smoothing.
    '''

    n_points = Property(Int)
    '''Fit with n-th-order polynomial
    '''

    def _get_n_points(self):
        n_points = int(self.n_smooth_window_fraction * len(self.eps))
        return n_points

    eps_smooth = Property(Array('float_'), output=True,
                          depends_on='input_change')
    '''Smoothied strains.
    '''
    @cached_property
    def _get_eps_smooth(self):
        return smooth(self.eps_asc, self.n_points, 'flat')

    cut_sig_c_eps_smoothed_with_E_c_linear = Bool(False)
    '''Processing parameter for smoothing - if true, initial
    linear branch of the smoothed response consists
    of the calculated linear stiffness and not of the measured
    response that can be distorted by geometric imperfections.
    '''

    sig_c_smooth = Property(Array('float_'), output=True,
                            depends_on='input_change')
    '''Smoothed composite stress.
    '''
    @cached_property
    def _get_sig_c_smooth(self):
        sig_c_smooth = smooth(self.sig_c_asc, self.n_points, 'flat')
        if self.cut_sig_c_eps_smoothed_with_E_c_linear:
            sig_lin = self.E_c * self.eps_smooth
            cut_sig = where(sig_c_smooth > sig_lin)
            sig_c_smooth[cut_sig] = sig_lin[cut_sig]
        return sig_c_smooth

    sig_tex_smooth = Property(Array('float_'), output=True,
                              depends_on='input_change')
    '''Smoothed stress in textile fabrics.
    '''
    @cached_property
    def _get_sig_tex_smooth(self):
        return smooth(self.sig_tex_asc, self.n_points, 'flat')

    # -------------------------------------------------------------------------
    # plot templates
    # -------------------------------------------------------------------------

    plot_templates = {'force / gauge displacement': '_plot_force_displacement',
                      'force / gauge displacement (ascending)': '_plot_force_displacement_asc',
                      'composite stress / strain': '_plot_sigc_eps',
                      'ironed composite stress / strain': '_plot_sigc_eps_ironed',
                      'interpolated composite stress / strain': '_plot_sigc_eps_interpolated',
                      'interpolated composite stress / strain (smoothed)': '_plot_sigc_eps_interpolated_smoothed',
                      'ironed textile stress / strain': '_plot_sigtex_eps_ironed',
                      'interpolated textile stress / strain': '_plot_sigtex_eps_interpolated',
                      'smoothed composite stress / strain': '_plot_sigc_eps_smoothed',
                      'textile stress / strain': '_plot_sigtex_eps',
                      'smoothed textile stress / strain': '_plot_sigtex_eps_smoothed',
                      }

    default_plot_template = 'force / gauge displacement'

    def _plot_force_displacement(self, axes):
        '''plot force-displacement diagram
        '''
        if hasattr(self, "W10_re") and hasattr(self, "W10_li") and hasattr(self, "W10_vo"):
            #
            axes.plot(self.W10_re, self.Kraft)
            axes.plot(self.W10_li, self.Kraft)
            axes.plot(self.W10_vo, self.Kraft)
#            axes.set_xlabel('%s' % ('displacement [mm]',))
#            axes.set_ylabel('%s' % ('force [kN]',))
        if hasattr(self, "WA_VL") and hasattr(self, "WA_VR") and hasattr(self, "WA_HL") and hasattr(self, "WA_HR"):
            #
            axes.plot(self.WA_VL, self.Kraft)
            axes.plot(self.WA_VR, self.Kraft)
            axes.plot(self.WA_HL, self.Kraft)
            axes.plot(self.WA_HR, self.Kraft)
#            axes.set_xlabel('%s' % ('displacement [mm]',))
#            axes.set_ylabel('%s' % ('force [kN]',))

    def _plot_force_displacement_asc(self, axes):
        '''plot force-displacement diagram (only the ascending branch)
        '''
        if hasattr(self, "W10_re") and hasattr(self, "W10_li") and hasattr(self, "W10_vo"):
            #
            axes.plot(self.W10_re[:self.max_stress_idx + 1], self.F_asc)
            axes.plot(self.W10_li[:self.max_stress_idx + 1], self.F_asc)
            axes.plot(self.W10_vo[:self.max_stress_idx + 1], self.F_asc)
#            axes.set_xlabel('%s' % ('displacement [mm]',))
#            axes.set_ylabel('%s' % ('force [kN]',))
        if hasattr(self, "WA_VL") and hasattr(self, "WA_VR") and hasattr(self, "WA_HL") and hasattr(self, "WA_HR"):
            #
            axes.plot(self.WA_VL[:self.max_stress_idx + 1], self.F_asc)
            axes.plot(self.WA_VR[:self.max_stress_idx + 1], self.F_asc)
            axes.plot(self.WA_HL[:self.max_stress_idx + 1], self.F_asc)
            axes.plot(self.WA_HR[:self.max_stress_idx + 1], self.F_asc)
#            axes.set_xlabel('%s' % ('displacement [mm]',))
#            axes.set_ylabel('%s' % ('force [kN]',))

    def _plot_sigc_eps(self, axes, color='black', linewidth=1., linestyle='-'):
        '''plot composite stress-strain diagram
        '''
        axes.plot(self.eps_asc, self.sig_c_asc,
                  color=color, linewidth=linewidth, linestyle=linestyle)
        axes.set_xlabel('strain [-]')
        axes.set_ylabel('composite stress [MPa]')

    def _plot_sigc_eps_smoothed(self, axes, color='black', linewidth=1., linestyle='-'):
        '''plot smoothed composite stress-strain diagram
        '''
        axes.plot(
            self.eps_smooth, self.sig_c_smooth, color='blue', linewidth=2)
        axes.set_xlabel('strain [-]')
        axes.set_ylabel('composite stress [MPa]')
        # original curve
        #
        axes.plot(self.eps_asc, self.sig_c_asc,
                  color=color, linewidth=linewidth, linestyle=linestyle)
# secant  stiffness
#        #
#        sig_lin = array([0, self.sig_c_max], dtype='float_')
#        eps_lin = array([0, self.sig_c_max / self.E_c ], dtype='float_')
#        axes.plot(eps_lin, sig_lin, color='red')

    def _plot_sigc_eps_ironed(self, axes):
        '''plot smoothed composite stress-strain diagram without unwanted unloading/reloading paths
        due to sliding in the buttstrap clamping
        '''
        axes.plot(self.eps_ironed, self.sig_c_ironed, color='green')
        axes.set_xlabel('strain [-]')
        axes.set_ylabel('composite stress [MPa]')
        # plot the stiffness of the composite (K_I) - uncracked state)
        #
        K_I = self.E_c  # depending of the testing age
#        K_I = self.E_c28
        eps_lin = array([0, self.sig_c_max / K_I], dtype='float_')
        sig_lin = array([0, self.sig_c_max], dtype='float_')
        axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')
        # plot the stiffness of the garn (K_IIb - cracked state)
        #
        E_tex = self.ccs.E_tex
        K_III = E_tex * self.A_tex / (self.A_c * 1000000.)
        eps_lin = array([0, self.eps_max], dtype='float_')
        sig_lin = array([0, self.eps_max * K_III], dtype='float_')
        axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')
        # original curve
        #
#        axes.plot(self.eps_asc, self.sig_c_asc, color='black')

    def _plot_sigc_eps_interpolated(self, axes):
        '''plot ironed composite stress-strain diagram starting at the origin,
        i.e. shift the strain by the offset resulting from the
        initial strain and the analytic composite stiffness
        '''
        axes.plot(
            self.eps_c_interpolated, self.sig_c_interpolated, color='green')
        axes.set_xlabel('strain [-]')
        axes.set_ylabel('composite stress [MPa]')
        # plot the stiffness of the composite (K_I) - uncracked state)
        #
        K_I = self.E_c  # depending of the testing age
#        K_I = self.E_c28
        eps_lin = array([0, self.sig_c_max / K_I], dtype='float_')
        sig_lin = array([0, self.sig_c_max], dtype='float_')
        axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')
        # plot the stiffness of the garn (K_IIb - cracked state)
        #
        E_tex = self.ccs.E_tex
        K_III = E_tex * self.rho_c
        eps_lin = array([0, self.eps_max], dtype='float_')
        sig_lin = array([0, self.eps_max * K_III], dtype='float_')
        axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')

    def _plot_sigc_eps_interpolated_smoothed(self, axes):
        '''plot ironed composite stress-strain diagram starting at the origin,
        i.e. shift the strain by the offset resulting from the
        initial strain and the analytic composite stiffness
        '''
        axes.plot(
            self.eps_c_interpolated, self.sig_c_interpolated, color='black')
        axes.plot(self.eps_c_interpolated_smoothed,
                  self.sig_c_interpolated_smoothed, color='green')
        axes.set_xlabel('strain [-]')
        axes.set_ylabel('composite stress [MPa]')
        # plot the stiffness of the composite (K_I) - uncracked state)
        #
        K_I = self.E_c  # depending of the testing age
#        K_I = self.E_c28
        eps_lin = array([0, self.sig_c_max / K_I], dtype='float_')
        sig_lin = array([0, self.sig_c_max], dtype='float_')
        axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')
        # plot the stiffness of the garn (K_IIb - cracked state)
        #
        E_tex = self.ccs.E_tex
        K_III = E_tex * self.rho_c
        eps_lin = array([0, self.eps_max], dtype='float_')
        sig_lin = array([0, self.eps_max * K_III], dtype='float_')
        axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')

    def _plot_sigtex_eps_ironed(self, axes):
        '''plot smoothed (textile) stress-strain diagram without unwanted unloading/reloading paths
        due to sliding in the buttstrap clamping
        '''
        axes.plot(self.eps_ironed, self.sig_tex_ironed, color='green')
        axes.set_xlabel('strain [-]')
        axes.set_ylabel('textile stress [MPa]')
        # yarn stiffness
        #
        E_tex = self.ccs.E_tex
        eps_lin = array([0, self.eps_max], dtype='float_')
        sig_lin = array([0, self.eps_max * E_tex], dtype='float_')
        axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')

    def _plot_sigtex_eps_interpolated(self, axes):
        '''plot ironed textile stress-strain diagram starting at the origin,
        i.e. shift the strain by the offset resulting from the
        initial strain and the analytic yarn stiffness
        '''
        axes.plot(
            self.eps_tex_interpolated, self.sig_tex_interpolated, color='green')
        axes.set_xlabel('strain [-]')
        axes.set_ylabel('textile stress [MPa]')
        # plot the stiffness of the composite (K_I) - uncracked state)
        #
        K_I = self.E_c / self.rho_c
        eps_lin = array([0, self.sig_tex_max / K_I], dtype='float_')
        sig_lin = array([0, self.sig_tex_max], dtype='float_')
        axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')
        # plot the stiffness of the garn (K_IIb - cracked state)
        #
        E_tex = self.ccs.E_tex
        K_III = E_tex
        eps_lin = array([0, self.eps_max], dtype='float_')
        sig_lin = array([0, self.eps_max * K_III], dtype='float_')
        axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')

    def _plot_sigtex_eps(self, axes, color='blue', linewidth=1., linestyle='-'):
        axes.plot(self.eps_asc, self.sig_tex_asc,
                  color=color, linewidth=linewidth, linestyle=linestyle)
        axes.set_xlabel('strain [-]')
        axes.set_ylabel('textile stress [MPa]')
        # original curve
        #
#        axes.plot(self.eps_asc, self.sig_tex_asc)
        # plot the textile secant stiffness at fracture state
        #
        eps_lin = array([0, self.eps_smooth[-1]], dtype='float_')
        sig_lin = self.ccs.E_tex_arr[1] * eps_lin
        axes.plot(eps_lin, sig_lin)

    def _plot_sigtex_eps_smoothed(self, axes, color='blue', linewidth=1., linestyle='-'):
        axes.plot(self.eps_smooth, self.sig_tex_smooth,
                  color=color, linewidth=linewidth, linestyle=linestyle)
        axes.set_xlabel('strain [-]')
        axes.set_ylabel('textile stress [MPa]')
        # original curve
        #
#        axes.plot(self.eps_asc, self.sig_tex_asc)
        # plot the textile secant stiffness at fracture state
        #
        eps_lin = array([0, self.eps_smooth[-1]], dtype='float_')
        sig_lin = self.ccs.E_tex * eps_lin
        axes.plot(eps_lin, sig_lin)

    # scaleable plotting methods
    #
    def _plot_tex_stress_strain_asc(self, axes, color='blue', linewidth=1.0, linestyle='-', label=None, f=None, xscale=1., k_rho=1.0, plot_analytical_stiffness_I=True, plot_analytical_stiffness_II=True, plot_analytical_stiffness=False, interpolated=True):
        '''plot the textile stress-strain curve; plot styles are configurable; analytical stiffness values are displayed if desired;
        '''
        #---------------
        # plot stiffness KI and KII if option is set to True
        #---------------
        K_I = self.E_c / self.rho_c * k_rho
        print 'K_I = E_c (simdb)', self.E_c
        rho_new = self.rho_c / k_rho
        K_I = self.E_c / rho_new
        print 'K_I = E_c (new)', K_I
        E_tex = self.ccs.E_tex
        K_IIb = E_tex
        print 'K_IIb = E_tex ', K_IIb
        #---------------
        if plot_analytical_stiffness == True:
            plot_analytical_stiffness_I = True
            plot_analytical_stiffness_II = True
        #---------------
        if plot_analytical_stiffness_I == True:
            print 'plot analytical stiffness (K_I)'
            # plot the stiffness of the composite (K_I) - uncracked state)
            eps_lin = array(
                [0, self.sig_tex_max / K_I], dtype='float_') * xscale
            sig_lin = array([0, self.sig_tex_max], dtype='float_')
            axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')
        #---------------
        if plot_analytical_stiffness_II == True:
            print 'plot analytical stiffness (K_IIb)'
            # plot the stiffness of the garn (K_IIb - cracked state)
            eps_lin = array([0, self.eps_max], dtype='float_') * xscale
            sig_lin = array([0, self.eps_max * K_IIb], dtype='float_')
            axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')
        #---------------
        # plot stress-strain curves
        #---------------
        if interpolated == True:
            # use ironed date (without initial offset)
            # scale by scale-factor scale_factor = 1000. for setting strain
            # unite to "permile"
            eps_asc_scaled = self.eps_c_interpolated * xscale
            sig_tex_interpolated = k_rho * self.sig_c_interpolated / self.rho_c
        else:
            # use ironed date (still contains initial offset)
            # scale by scale-factor scale_factor = 1000. for setting strain
            # unite to "permile"
            eps_asc_scaled = self.eps_asc * xscale
            sig_tex_interpolated = k_rho * self.sig_c_asc / self.rho_c
        axes.plot(eps_asc_scaled, sig_tex_interpolated, color=color,
                  linewidth=linewidth, linestyle=linestyle, label=label)

    def _plot_comp_stress_strain_asc(self, axes, color='blue', linewidth=1.0, linestyle='-', label=None, f=None, xscale=1., k_rho=1.0, plot_analytical_stiffness_I=True, plot_analytical_stiffness_II=True, plot_analytical_stiffness=False, interpolated=True):
        '''plot the composite stress-strain curve; plot styles are configurable; analytical stiffness values are displayed if desired;
        '''
        #---------------
        # plot stiffness KI and KII if option is set to True
        #---------------
        K_I = self.E_c  # depending of the testing age
        print 'K_I = E_c (simdb)', self.E_c
        print 'E_tex (simdb)', self.ccs.E_tex
        print 'E_m (simdb)', self.E_m
        print 'age (simdb)', self.age
        rho_new = self.ccs.rho_c / k_rho
        K_I = (1 - rho_new) * self.E_m + rho_new * self.ccs.E_tex
        print 'K_I (new)', K_I
        print 'rho_c (new)', self.ccs.rho_c / k_rho
#        K_I = self.E_c28
        E_tex = self.ccs.E_tex
        K_IIb = E_tex * self.rho_c
        print 'K_IIb = E_tex * self.rho_c', K_IIb
        K_IIb = E_tex * rho_new
        print 'K_IIb = E_tex * rho_new', K_IIb
        # ---------------
        if plot_analytical_stiffness == True:
            plot_analytical_stiffness_I = True
            plot_analytical_stiffness_II = True
        # ---------------
        if plot_analytical_stiffness_I == True:
            print 'plot analytical stiffness (K_I)'
            eps_lin = array([0, self.sig_c_max / K_I], dtype='float_') * xscale
            sig_lin = array([0, self.sig_c_max], dtype='float_')
            axes.plot(
                eps_lin, sig_lin, color='grey', linestyle='--', linewidth=linewidth)
        # ---------------
        if plot_analytical_stiffness_II == True:
            print 'plot analytical stiffness (K_IIb)'
            # plot the stiffness of the garn (K_IIb - cracked state)
            eps_lin = array([0, self.eps_max], dtype='float_') * xscale
            sig_lin = array([0, self.eps_max * K_IIb], dtype='float_')
            axes.plot(
                eps_lin, sig_lin, color='grey', linestyle='--', linewidth=linewidth)
        # ---------------
        # plot stress-strain curves
        # ---------------
        if interpolated == True:
            # use ironed date (without initial offset)
            # scale by scale-factor scale_factor = 1000. for setting strain
            # unite to "permile"
            eps_asc_scaled = self.eps_c_interpolated * xscale
            axes.plot(eps_asc_scaled, self.sig_c_interpolated, color=color,
                      linewidth=linewidth, linestyle=linestyle, label=label)

            xdata = eps_asc_scaled
            ydata = self.sig_c_interpolated
#            TT_arr = np.hstack([xdata[:, None], ydata[:, None]])
#            print 'TT_arr ', TT_arr.shape, ' save to file "TT_arr.csv"'
#            np.savetxt('TT_arr.csv', TT_arr, delimiter=';')

#            idx_sig_c_mean = np.where(self.sig_c_interpolated >= 16.607)[0][0]
#            print 'sig_c_mean = ', self.sig_c_interpolated[idx_sig_c_mean]
#            eps_c_mean = self.eps_c_interpolated[idx_sig_c_mean]
#            print 'eps_c_mean = ', eps_c_mean

#            sig_c_interpolated_smoothed = self.sig_c_interpolated_smoothed
#            eps_c_interpolated_smoothed = self.eps_c_interpolated_smoothed * xscale
#            axes.plot(eps_c_interpolated_smoothed, sig_c_interpolated_smoothed, color='black', linewidth=linewidth, linestyle=linestyle, label=label)
        else:
            # use ironed date (still contains initial offset)
            # scale by scale-factor scale_factor = 1000. for setting strain
            # unite to "permile"
            eps_asc_scaled = self.eps_asc * xscale
            axes.plot(eps_asc_scaled, self.sig_c_asc, color=color,
                      linewidth=linewidth, linestyle=linestyle, label=label)
# eps_asc_scaled = self.eps_ironed * xscale  # scale by scale-factor scale_factor = 1000. for setting strain unite to "permile"
#            axes.plot(eps_asc_scaled, self.sig_c_ironed, color=color, linewidth=linewidth, linestyle=linestyle, label=label)

    # ---------------------------------
    # view
    # ---------------------------------

    traits_view = View(VSplit(
        HSplit(Group(
            Item('width', format_str="%.3f"),
            Item('gauge_length', format_str="%.3f"),
            springy=True,
            label='geometry',
            id='matresdev.db.exdb.ex_composite_tensile_test.geometry',
            dock='tab',
        ),
            Group(
            Item('loading_rate'),
            Item('age'),
            springy=True,
            label='loading rate and age',
            id='matresdev.db.exdb.ex_composite_tensile_test.loading',
            dock='tab',),
            id='matresdev.db.exdb.ex_composite_tensile_test.xxx',
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
                 style='readonly', show_label=True, format_str="%.0f"),
            Item('sig_c_max', visible_when='derived_data_available',
                 style='readonly', emphasized=True, format_str="%.2f"),
            Item('sig_tex_max', visible_when='derived_data_available',
                 style='readonly', emphasized=True, format_str="%.2f"),
            Item('eps_max', visible_when='derived_data_available',
                 style='readonly', emphasized=True, format_str="%.4f"),
            label='output characteristics',
            id='matresdev.db.exdb.ex_composite_tensile_test.vgroup.outputs',
            dock='tab',
            scrollable=True,
        ),
        scrollable=True,
        id='matresdev.db.exdb.ex_composite_tensile_test.vgroup',
        dock='tab',
    ),
        id='matresdev.db.exdb.ex_composite_tensile_test',
        dock='tab',
        scrollable=True,
        resizable=True,
        height=0.8,
        width=0.5,
    )

ExpTTDB.db = ExRunClassExt(klass=ExpTTDB)

if __name__ == '__main__':

    #    ExpTTDB.add_class_trait('production_date', Date(input=True, table_field=True,))
    #    for inst in ExpTTDB.db.inst_list:
    #        print inst.key
    #        print inst.add_trait('production_date', Date('14/9/2011', input=True, table_field=True,))
    #        print inst.production_date
    #        inst.save()

    print ExpTTDB.db['TT-12c-6cm-0-TU-WMTC'].sig_c_interpolated_smoothed

    # ExpTTDB.db.configure_traits()
