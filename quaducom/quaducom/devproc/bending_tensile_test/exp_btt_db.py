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
    DelegatesTo, Date

from numpy import \
    array, where, argmax

import numpy as np

from etsproxy.traits.ui.api \
    import View, Item, HSplit, Group, VSplit

from matresdev.db.exdb.ex_type import ExType
from matresdev.db.exdb.i_ex_type import IExType

from mathkit.array.smoothing import smooth
from mathkit.mfn import MFnLineArray

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

    #--------------------------------------------------------------------------------
    # specify inputs:
    #--------------------------------------------------------------------------------

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

    #--------------------------------------------------------------------------------
    # define processing
    #--------------------------------------------------------------------------------

    def process_source_data(self):
        '''read in the measured data from file and assign
        attributes after array processing.
        If necessary modify the assigned data, e.i. change
        the sign or specify an offset for the specific test setup.
        '''
        super(ExpBTTDB, self).process_source_data()

        self.N = self.F160kN
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

    #-------------------------------------------------------------------------------
    # Get the maximum force index to cut off the descending part of the curves
    #-------------------------------------------------------------------------------

    max_F_idx = Property(Int, depends_on='input_change')
    @cached_property
    def _get_max_F_idx(self):
        return argmax(self.F)

    max_N_idx = Property(Int, depends_on='input_change')
    @cached_property
    def _get_max_N_idx(self):
        return argmax(self.N)

    F_max1 = Property(Float, depends_on='input_change',
                            output=True, table_field=True, unit='kN')
    '''Only needed for the distinction between pure tensile test and the bending test
    '''
    @cached_property
    def _get_F_max1(self):
        return max(self.F)

    #-------------------------------------------------------------------------------
    # Get only the ascending branch of the response curve
    #-------------------------------------------------------------------------------

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

    #-------------------------------------------------------------------------------
    # Cut the curves at the end when the gradient of displacement is to big
    #-------------------------------------------------------------------------------

    gra_max_threshold = Float(0.12, auto_set=False, enter_set=True)
    '''Threshold to limit the gradient of the displacement curve at the end.
    '''

    w_cut_idx = Property(Int, depends_on='input_change')
    @cached_property
    def _get_w_cut_idx(self):

        w_asc = self.w_asc
        t_asc = self.t_asc

        # Calculate the deltas with beginning at indize 840 that is equivalent to t = 420 sec,
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

    w_cut_asc = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_w_cut_asc(self):
        ''' Method to cut the end of the displacement curve if the gradient is to big
        due to the kind of failure in the testing
        '''
        return self.w_asc[:self.w_cut_idx]

    u_cut_asc = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_u_cut_asc(self):
        ''' Method to cut the end of the displacement curve if the gradient is to big
        due to the kind of failure in the testing
        '''
        return self.u_asc[:self.w_cut_idx]

    F_cut_asc = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_F_cut_asc(self):
        ''' Method to cut the end of the displacement curve if the gradient is to big
        due to the kind of failure in the testing
        '''
        return self.F[:self.w_cut_idx]

    N_cut_asc = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_N_cut_asc(self):
        ''' Method to cut the end of the displacement curve if the gradient is to big
        due to the kind of failure in the testing
        '''
        return self.N[:self.w_cut_idx]

    t_cut_asc = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_t_cut_asc(self):
        ''' Method to cut the end of the displacement curve if the gradient is to big
        due to the kind of failure in the testing
        '''
        return self.t[:self.w_cut_idx]

    #-------------------------------------------------------------------------------
    # Method to eliminate the negative or positive deformation at the beginning of the experiment
    # due to predeformation of the specimen
    #-------------------------------------------------------------------------------

    # get the minimum value of deformation that is equivalent to the predeformation
    w_pred = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_w_pred(self):
        w_pred = min(self.w_cut_asc)
        return w_pred

    F_lim_threshold = Float(0.142, auto_set=False, enter_set=True)
    '''Threshold for the lower limit of force. For all force values smaller than this limit,
        equate the deformation value with zero
    '''

    w_el_pred = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_w_el_pred(self):
        w = self.w_cut_asc
        F = self.F_cut_asc

        # get the indices for the force smaller than force threshold
        idx_F_lim = where(F <= self.F_lim_threshold)[0]

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

    #-------------------------------------------------------------------------------
    # Get the moment arrays
    #-------------------------------------------------------------------------------

    M = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_M(self):
        return self.F_cut_asc * self.length / 4 - self.N_cut_asc * self.w_el_pred / 1000

    MN = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_MN(self):
        return self.N_cut_asc * self.w_el_pred / 1000 * -1

    MF = Property(Array('float_'), depends_on='input_change')
    @cached_property
    def _get_MF(self):
        return self.F_cut_asc * self.length / 4

    #-------------------------------------------------------------------------------
    # Get maximum values of the variables
    #-------------------------------------------------------------------------------

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

    t_max = Property(Float, depends_on='input_change',
                            output=True, table_field=True, unit='kNm')
    @cached_property
    def _get_t_max(self):
        return self.t_cut_asc[-1]

    #--------------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------------

    plot_templates = {'N(t), F(t)' : '_plot_N_F_t',
                      'N_cut(t), F_cut(t)' : '_plot_N_F_t_cut',
                      'M(t), MN(t), MF(t)' : '_plot_M_MN_MF_t',
                      'N(M), N(MN), N(MF)' : '_plot_N_M_MN_MF',
                      'w(t)' : '_plot_w_t',
                      'w_cut(t)' : '_plot_w_t_cut',
                      'u(t)' : '_plot_u_t',
                      'u_cut(t)' : '_plot_u_t_cut',
                      'F(w)_cut' : '_plot_F_w_cut',
                      'F(w)_el_pred' : '_plot_F_w_el_pred',
                                             }

    default_plot_template = 'N(t), F(t)'

    def _plot_N_F_t(self, axes):
        '''Normal force versus time
        '''
        # ax = twiny()
        axes.plot(self.t_asc, self.N_asc, color='blue')
        axes.plot(self.t_asc, self.F_asc, color='red')

        axes.set_xlabel('t [sec]')
        axes.set_ylabel('N / F [kN]')
        # ax.set_ylabel('F [kN]')
        axes.set_ylim(0, 50)
        # ax.set_ylim(0, 5)

    def _plot_N_F_t_cut(self, axes):
        '''Normal force_cut versus time_cut
        '''
        axes.plot(self.t_cut_asc, self.N_cut_asc, color='blue')
        axes.plot(self.t_cut_asc, self.F_cut_asc, color='red')

        axes.set_xlabel('t [sec]')
        axes.set_ylabel('N / F [kN]')
        axes.set_ylim(0, 50)

    def _plot_M_MN_MF_t(self, axes):
        '''Moment versus time
        '''
        axes.plot(self.t_cut_asc, self.M, color='green')
        axes.plot(self.t_cut_asc, self.MN, color='blue')
        axes.plot(self.t_cut_asc, self.MF, color='red')

        axes.set_xlabel('t [sec]')
        axes.set_ylabel('M / MN / MF [kNm]')
        axes.set_ylim(-0.3, 0.45)

    def _plot_N_M_MN_MF(self, axes):
        '''Normal force versus moment
        '''
        axes.plot(self.M, self.N_cut_asc, color='green')
        axes.plot(self.MN, self.N_cut_asc, color='blue')
        axes.plot(self.MF, self.N_cut_asc, color='red')

        axes.set_xlabel('M / MN / MF [kNm]')
        axes.set_ylabel('N [kN]')
        axes.set_ylim(50, 0)
        axes.set_xlim(-0.3, 0.45)

    def _plot_w_t(self, axes):
        '''displacement versus time
        '''
        axes.plot(self.t_asc, self.w_asc, color='black')

        axes.set_xlabel('t [sec]')
        axes.set_ylabel('w [mm]')
        axes.set_ylim(0, 10)

    def _plot_w_t_cut(self, axes):
        '''displacement_cut versus time_cut
        '''
        axes.plot(self.t_cut_asc, self.w_cut_asc, color='black')

        axes.set_xlabel('t [sec]')
        axes.set_ylabel('w [mm]')
        axes.set_ylim(0, 10)

    def _plot_u_t(self, axes):
        '''displacement versus time
        '''
        axes.plot(self.t_asc, self.u_asc, color='black')

        axes.set_xlabel('t [sec]')
        axes.set_ylabel('u [mm]')
        axes.set_ylim(0, 30)

    def _plot_u_t_cut(self, axes):
        '''displacement_cut versus time_cut
        '''
        axes.plot(self.t_cut_asc, self.u_cut_asc, color='black')

        axes.set_xlabel('t [sec]')
        axes.set_ylabel('u [mm]')
        axes.set_ylim(0, 30)

    def _plot_F_w_cut(self, axes):
        '''Normal force_cut versus displacement_cut
        '''
        axes.plot(self.w_cut_asc, self.F_cut_asc, color='blue')

        axes.set_xlabel('w [mm]')
        axes.set_ylabel('F [kN]')
        axes.set_ylim(0, 7)
        axes.set_xlim(-2, 10)

    def _plot_F_w_el_pred(self, axes):
        '''Normal force_cut versus displacement with eliminated predeformation
        '''
        axes.plot(self.w_el_pred, self.F_cut_asc, color='blue')

        axes.set_xlabel('w [mm]')
        axes.set_ylabel('F [kN]')
        axes.set_ylim(0, 7)
        axes.set_xlim(-2, 10)

    #---------------------------------
    # view
    #---------------------------------

    traits_view = View(VSplit(
                         HSplit(Group(
                                  Item('width'       , format_str="%.3f"),
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

    test_file = os.path.join(simdb.exdata_dir,
                           'bending_tensile_test',
                           '2014-06-12_BTT-6c-2cm-0-TU_MxN2',
                           'BTT-6c-2cm-TU-0-V03_MxN2.DAT')

#     e1 = ExpBTTDB(data_file=test_file)
#     e1.process_source_data()
#     e1
#     print 'F_max1', e1.F_max1
#
#     p.plot(e1.w_cut_asc, e1.F_cut_asc, linewidth=3)
#     p.plot(e1.w_asc, e1.F_asc)
#     p.show()

    doe_reader = ExRunView(data_file=test_file)
    doe_reader.configure_traits()