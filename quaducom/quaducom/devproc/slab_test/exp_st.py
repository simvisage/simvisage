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
    Dict, Button, Bool, Enum, Event, implements, DelegatesTo

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

from matresdev.db.simdb import SimDB
simdb = SimDB()

import csv

from numpy import array, fabs, where, copy, ones, argsort

from numpy import \
    loadtxt, argmax, polyfit, poly1d, frompyfunc, dot

import numpy as np

from etsproxy.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter

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

data_file_editor = FileEditor(filter=['*.DAT'])

from mathkit.array.smoothing import smooth

from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp

from matresdev.db.matdb.trc.fabric_layout \
    import FabricLayOut

from matresdev.db.matdb.trc.concrete_mixture \
    import ConcreteMixture

from matresdev.db.matdb.trc.composite_cross_section import \
    CompositeCrossSection, plain_concrete

class ExpST(ExType):
    '''Experiment: Slab Test
    '''
#    label = Str('slab test')
    implements(IExType)

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

    edge_length = Float(1.25, unit='m', input=True, table_field=True,
                           auto_set=False, enter_set=True)
    thickness = Float(0.06, unit='m', input=True, table_field=True,
                           auto_set=False, enter_set=True)

    # age of the concrete at the time of testing
    age = Int(28, unit='d', input=True, table_field=True,
                             auto_set=False, enter_set=True)
    loading_rate = Float(2.0, unit='mm/min', input=True, table_field=True,
                            auto_set=False, enter_set=True)

    #--------------------------------------------------------------------------
    # composite cross section
    #--------------------------------------------------------------------------

    ccs = Instance(CompositeCrossSection)
    def _ccs_default(self):
        '''default settings correspond to
        setup '7u_MAG-07-03_PZ-0708-1'
        '''
        fabric_layout_key = '2D-05-11'
#        fabric_layout_key = 'MAG-07-03'
#        fabric_layout_key = '2D-02-06a'
#        concrete_mixture_key = 'PZ-0708-1'
#        concrete_mixture_key = 'FIL-10-09'
#        concrete_mixture_key = 'barrelshell'
        concrete_mixture_key = 'PZ-0708-1'
#        orientation_fn_key = 'all0'
        orientation_fn_key = '90_0'
        n_layers = 2
        s_tex_z = 0.015 / (n_layers + 1)
        ccs = CompositeCrossSection (
                    fabric_layup_list=[
                            plain_concrete(0.035),
#                            plain_concrete(s_tex_z * 0.5),
                            FabricLayUp (
                                   n_layers=n_layers,
                                   orientation_fn_key=orientation_fn_key,
                                   s_tex_z=s_tex_z,
                                   fabric_layout_key=fabric_layout_key
                                   ),
#                            plain_concrete(s_tex_z * 0.5)
                            plain_concrete(0.5)
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
    jump_rtol = Float(0.1,
                      auto_set=False, enter_set=True,
                      ironing_param=True)


    data_array_ironed = Property(Array(float),
                                  depends_on='data_array, +ironing_param, +axis_selection')
    @cached_property
    def _get_data_array_ironed(self):
        '''remove the jumps in the displacement curves
        due to resetting the displacement gauges.
        '''
        print('*** curve ironing activated ***')

        # each column from the data array corresponds to a measured parameter
        # e.g. displacement at a given point as function of time u = f(t))
        #
        data_array_ironed = copy(self.data_array)

        for idx in range(self.data_array.shape[1]):

            # use ironing method only for columns of the displacement gauges.
            #
            if self.names_and_units[0][ idx ] != 'Kraft' and \
                self.names_and_units[0][ idx ] != 'Bezugskanal' and \
                self.names_and_units[0][ idx ] != 'Weg':

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

                print('number of jumps removed in data_arr_ironed for', self.names_and_units[0][ idx ], ': ', jump_idx.shape[0])
                # glue the curve at each jump together
                for jidx in jump_idx:
                    # get the offsets at each jump of the curve
                    shift = data_arr[jidx + 1] - data_arr[jidx]
                    # shift all succeeding values by the calculated offset
                    data_arr[jidx + 1:] -= shift

                data_array_ironed[:, idx] = data_arr[:]

        return data_array_ironed


    @on_trait_change('+ironing_param')
    def process_source_data(self):
        '''read in the measured data from file and assign
        attributes after array processing.

        NOTE: if center displacement gauge ('WA_M') is missing the measured
        displacement of the cylinder ('Weg') is used instead.
        A minor mistake is made depending on how much time passes
        before the cylinder has contact with the slab.
        '''
        print('*** process source data ***')

        self._read_data_array()
        # curve ironing:
        #
        self.processed_data_array = self.data_array_ironed

        # set attributes:
        #
        self._set_array_attribs()

        if 'WA_M' not in self.factor_list:
            print('*** NOTE: Displacement gauge at center ("WA_M") missing. Cylinder displacement ("Weg") is used instead! ***')
            self.WA_M = self.Weg


    #--------------------------------------------------------------------------------
    # plot templates
    #--------------------------------------------------------------------------------

    plot_templates = {'force / deflection (all)' : '_plot_force_deflection',
#                      'force / average deflection (c; ce; e) interpolated'    : '_plot_force_deflection_avg_interpolated',
#                      'force / deflection (center)'               : '_plot_force_center_deflection',
#                      'smoothed force / deflection (center)'      : '_plot_force_center_deflection_smoothed',
#                      'force / deflection (edges)'                : '_plot_force_edge_deflection',
#                      'force / deflection (center-edges)'         : '_plot_force_center_edge_deflection',
#                      'force / average deflection (edges)'        : '_plot_force_edge_deflection_avg',
#                      'force / average deflection (center-edges)' : '_plot_force_center_edge_deflection_avg',
#                      'force / average deflection (c; ce; e)'     : '_plot_force_deflection_avg',
#                      'force / deflection (center) interpolated'  : '_plot_force_center_deflection_interpolated',
#                      'force / deflection (corner)'               : '_plot_force_corner_deflection',
#                      'edge_deflection_avg (front/back) / edge_deflection_avg (left/right)' : '_plot_edge_deflection_edge_deflection_avg',

                     }

    default_plot_template = 'force / deflection (center)'


    n_fit_window_fraction = Float(0.1)

    # get only the ascending branch of the response curve
    #
    max_force_idx = Property(Int)
    def _get_max_force_idx(self):
        '''get the index of the maximum force'''
        return argmax(-self.Kraft)

    f_asc = Property(Array)
    def _get_f_asc(self):
        '''get only the ascending branch of the response curve'''
        return -self.Kraft[:self.max_force_idx + 1]

    w_asc = Property(Array)
    def _get_w_asc(self):
        '''get only the ascending branch of the response curve'''
        return -self.WA_M[:self.max_force_idx + 1]

    n_points = Property(Int)
    def _get_n_points(self):
        return int(self.n_fit_window_fraction * len(self.w_asc))

    f_smooth = Property()
    def _get_f_smooth(self):
        return smooth(self.f_asc, self.n_points, 'flat')

    w_smooth = Property()
    def _get_w_smooth(self):
        return smooth(self.w_asc, self.n_points, 'flat')

    def _plot_force_deflection(self, axes, offset_w=0., color='blue', linewidth=1.5, linestyle='-'):
        '''plot the F-w-diagramm for all displacement gauges including maschine displacement
        '''
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'

        max_force_idx = self.max_force_idx
#        max_force_idx = -2
        f_asc = -self.Kraft[:max_force_idx + 1]

        print('self.factor_list', self.factor_list)
        header_string = ''
        for i in self.factor_list[2:]:
            header_string = header_string + i + '; '
            T_arr = -getattr(self, i)
            T_asc = T_arr[:max_force_idx + 1]
            axes.plot(T_asc, f_asc, label=i)

        axes.legend()

        img_dir = os.path.join(simdb.exdata_dir, 'img_dir')
        # check if directory exist otherwise create
        #
        if os.path.isdir(img_dir) == False:
            os.makedirs(img_dir)
        test_series_dir = os.path.join(img_dir, 'slab_test_astark')
        # check if directory exist otherwise create
        #
        if os.path.isdir(test_series_dir) == False:
            os.makedirs(test_series_dir)
        out_file = os.path.join(test_series_dir, self.key)

        np.savetxt(out_file, self.processed_data_array, delimiter=';')

        # workaround for old numpy.savetxt version:
        f = open(out_file, 'r')
        temp = f.read()
        f.close()

        f = open(out_file, 'w')
        f.write(header_string)

        f.write(temp)
        f.close()



    def _plot_force_center_deflection(self, axes, offset_w=0., color='blue', linewidth=1.5, linestyle='-'):
        '''plot the F-w-diagramm for the center (c) deflection
        '''
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        xdata = -self.WA_M
        ydata = -self.Kraft

        xdata += offset_w
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata, color=color, linewidth=linewidth, linestyle=linestyle)

    def _plot_force_corner_deflection(self, axes):
        '''plot the F-w-diagramm for the corner deflection (at the center of one of the supports)
        '''
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        xdata = -self.WA_Eck
        ydata = -self.Kraft

#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))
        axes.plot(xdata, ydata
                       # color = c, linewidth = w, linestyle = s
                       )

    def _plot_force_center_deflection_smoothed(self, axes):
        '''plot the F-w-diagramm for the center (c) deflection (smoothed curves)
        '''
        axes.plot(self.w_smooth, self.f_smooth, color='blue', linewidth=1)
#        secant_stiffness_w10 = (f_smooth[10] - f_smooth[0]) / (w_smooth[10] - w_smooth[0])
#        w0_lin = array([0.0, w_smooth[10] ], dtype = 'float_')
#        f0_lin = array([0.0, w_smooth[10] * secant_stiffness_w10 ], dtype = 'float_')
        # axes.plot( w0_lin, f0_lin, color = 'black' )

    def _plot_force_center_deflection_interpolated(self, axes, linestyle='-', linewidth=1.5, plot_elastic_stiffness=True):
        '''plot the F-w-diagramm for the center (c) deflection (interpolated initial stiffness)
        '''
        # get the index of the maximum stress
        max_force_idx = argmax(-self.Kraft)
        # get only the ascending branch of the response curve
        f_asc = -self.Kraft[:max_force_idx + 1]
        w_m = -self.WA_M[:max_force_idx + 1]

#        w_m -= 0.17
#        axes.plot(w_m, f_asc, color = 'blue', linewidth = 1)

        # move the starting point of the center deflection curve to the point where the force starts
        # (remove offset in measured displacement where there is still no force measured)
        #
        idx_0 = np.where(f_asc > 0.0)[0][0]

        f_asc_cut = np.hstack([f_asc[ idx_0: ]])
        w_m_cut = np.hstack([w_m[ idx_0: ]]) - w_m[ idx_0 ]
#        print 'f_asc_cut.shape', f_asc_cut.shape
#        fw_arr = np.hstack([f_asc_cut[:, None], w_m_cut[:, None]])
#        print 'fw_arr.shape', fw_arr.shape
#        np.savetxt('ST-6c-2cm-TU_bs2_f-w_asc.csv', fw_arr, delimiter=';')
        axes.plot(w_m_cut, f_asc_cut, color='k', linewidth=linewidth, linestyle=linestyle)

        # composite E-modulus
        #
        E_c = self.E_c
        print('E_c', E_c)

        if self.thickness == 0.02 and self.edge_length == 0.80 and plot_elastic_stiffness == True:
            K_linear = E_c / 24900. * 1.056  # [MN/m]=[kN/mm] bending stiffness with respect to center force
            max_f = f_asc_cut[-1]
            w_linear = np.array([0., max_f / K_linear])
            f_linear = np.array([0., max_f])
            axes.plot(w_linear, f_linear, linewidth=linewidth, linestyle='--', color='k')
        if self.thickness == 0.03 and self.edge_length == 1.25 and plot_elastic_stiffness == True:
            K_linear = E_c / 24900. * 1.267  # [MN/m]=[kN/mm] bending stiffness with respect to center force
            max_f = f_asc_cut[-1]
            w_linear = np.array([0., max_f / K_linear])
            f_linear = np.array([0., max_f])
            axes.plot(w_linear, f_linear, linewidth=linewidth, linestyle='--', color='k')
        if self.thickness == 0.06 and self.edge_length == 1.25 and plot_elastic_stiffness == True:
            K_linear = E_c / 24900. * 10.14  # [MN/m]=[kN/mm] bending stiffness with respect to center force
            max_f = f_asc_cut[-1]
            w_linear = np.array([0., max_f / K_linear])
            f_linear = np.array([0., max_f])
            axes.plot(w_linear, f_linear, linewidth=linewidth, linestyle='--', color='k')

    def _plot_force_edge_deflection(self, axes):
        '''plot the F-w-diagramm for the edge (e) deflections
        '''
        max_force_idx = self.max_force_idx
        f_asc = self.f_asc
        w_v_asc = -self.WA_V[:max_force_idx + 1]
        w_h_asc = -self.WA_H[:max_force_idx + 1]
        w_l_asc = -self.WA_L[:max_force_idx + 1]
        w_r_asc = -self.WA_R[:max_force_idx + 1]
        axes.plot(w_v_asc, f_asc, color='blue', linewidth=1)
        axes.plot(w_h_asc, f_asc, color='blue', linewidth=1)
        axes.plot(w_l_asc, f_asc, color='green', linewidth=1)
        axes.plot(w_r_asc, f_asc, color='green', linewidth=1)

    def _plot_force_edge_deflection_avg(self, axes):
        '''plot the average F-w-diagramm for the edge (e) deflections
        '''
        max_force_idx = self.max_force_idx
        f_asc = self.f_asc
        w_v_asc = -self.WA_V[:max_force_idx + 1]
        w_h_asc = -self.WA_H[:max_force_idx + 1]
        w_l_asc = -self.WA_L[:max_force_idx + 1]
        w_r_asc = -self.WA_R[:max_force_idx + 1]
        # get the average displacement values of the corresponding displacement gauges
        w_vh_asc = (w_v_asc + w_h_asc) / 2
        w_lr_asc = (w_l_asc + w_r_asc) / 2
        axes.plot(w_vh_asc, f_asc, color='blue', linewidth=1, label='w_vh')
        axes.plot(w_lr_asc, f_asc, color='blue', linewidth=1, label='w_lr')
        axes.legend()

    def _plot_edge_deflection_edge_deflection_avg(self, axes):
        '''plot the average edge (e) deflections for the front/back and left/right
        '''
        max_force_idx = self.max_force_idx
        w_v_asc = -self.WA_V[:max_force_idx + 1]
        w_h_asc = -self.WA_H[:max_force_idx + 1]
        w_l_asc = -self.WA_L[:max_force_idx + 1]
        w_r_asc = -self.WA_R[:max_force_idx + 1]
        # get the average displacement values of the corresponding displacement gauges
        w_vh_asc = (w_v_asc + w_h_asc) / 2
        w_lr_asc = (w_l_asc + w_r_asc) / 2
        axes.plot(w_vh_asc, w_lr_asc, color='blue', linewidth=1, label='w_vh / w_lr')
        axes.plot(np.array([0, w_vh_asc[-1]]), np.array([0, w_vh_asc[-1]]), color='k', linewidth=1, linestyle='-')
        axes.legend()

    def _plot_force_center_edge_deflection(self, axes):
        '''plot the F-w-diagramm for the center-edge (ce) deflections
        '''
        max_force_idx = argmax(-self.Kraft)
        # get only the ascending branch of the response curve
        f_asc = -self.Kraft[:max_force_idx + 1]
        w_ml_asc = -self.WA_ML[:max_force_idx + 1]
        w_mr_asc = -self.WA_MR[:max_force_idx + 1]
        axes.plot(w_ml_asc, f_asc, color='blue', linewidth=1, label='w_ml')
        axes.plot(w_mr_asc, f_asc, color='blue', linewidth=1, label='w_mr')
        axes.legend()

    def _plot_force_center_edge_deflection_avg(self, axes):
        '''plot the average F-w-diagramm for the center-edge (ce) deflections
        '''
        # get the index of the maximum stress
        max_force_idx = argmax(-self.Kraft)
        # get only the ascending branch of the response curve
        f_asc = -self.Kraft[:max_force_idx + 1]
        w_ml_asc = -self.WA_ML[:max_force_idx + 1]
        w_mr_asc = -self.WA_MR[:max_force_idx + 1]
        # get the average displacement values of the corresponding displacement gauges
        w_mlmr_asc = (w_ml_asc + w_mr_asc) / 2
        axes.plot(w_mlmr_asc, f_asc, color='blue', linewidth=1)

    def _plot_force_deflection_avg(self, axes):
        '''plot the average F-w-diagramms for the center(c), center-edge (ce) and edge(vh) and edge (lr) deflections
        '''
        # get the index of the maximum stress
        max_force_idx = argmax(-self.Kraft)
        # get only the ascending branch of the response curve
        f_asc = -self.Kraft[:max_force_idx + 1]

        w_m = -self.WA_M[:max_force_idx + 1]
        axes.plot(w_m, f_asc, color='blue', linewidth=1)

        # ## center-edge deflection (ce)
        w_ml_asc = -self.WA_ML[:max_force_idx + 1]
        w_mr_asc = -self.WA_MR[:max_force_idx + 1]
        # get the average displacement values of the corresponding displacement gauges
        w_mlmr_asc = (w_ml_asc + w_mr_asc) / 2
        axes.plot(w_mlmr_asc, f_asc, color='red', linewidth=1)

        # ## edge deflections (e)
        w_v_asc = -self.WA_V[:max_force_idx + 1]
        w_h_asc = -self.WA_H[:max_force_idx + 1]
        w_l_asc = -self.WA_L[:max_force_idx + 1]
        w_r_asc = -self.WA_R[:max_force_idx + 1]
        # get the average displacement values of the corresponding displacement gauges
        w_vh_asc = (w_v_asc + w_h_asc) / 2
        w_lr_asc = (w_l_asc + w_r_asc) / 2
        axes.plot(w_vh_asc, f_asc, color='green', linewidth=1, label='w_vh')
        axes.plot(w_lr_asc, f_asc, color='blue', linewidth=1, label='w_lr')

#        # set axis-labels
#        xkey = 'deflection [mm]'
#        ykey = 'force [kN]'
#        axes.set_xlabel('%s' % (xkey,))
#        axes.set_ylabel('%s' % (ykey,))

    def _plot_force_deflection_avg_interpolated(self, axes, linewidth=1):
        '''plot the average F-w-diagrams for the center(c), center-edge (ce) and edge(vh) and edge (lr) deflections
        NOTE: center deflection curve is cut at its starting point in order to remove offset in the dispacement meassurement
        '''
        # get the index of the maximum stress
        max_force_idx = argmax(-self.Kraft)
        # get only the ascending branch of the response curve
        f_asc = -self.Kraft[:max_force_idx + 1]
        w_m = -self.WA_M[:max_force_idx + 1]

#        w_m -= 0.17
#        axes.plot(w_m, f_asc, color = 'blue', linewidth = 1)

        # move the starting point of the center deflection curve to the point where the force starts
        # (remove offset in measured displacement where there is still no force measured)
        #
        idx_0 = np.where(f_asc > 0.05)[0][0]

        f_asc_cut = f_asc[ idx_0: ]
        w_m_cut = w_m[ idx_0: ] - w_m[ idx_0 ]
        axes.plot(w_m_cut, f_asc_cut, color='k', linewidth=linewidth)

        # plot machine displacement (hydraulic cylinder)
        #
#        Weg_asc = -self.Weg[ :max_force_idx + 1 ]
#        axes.plot(Weg_asc, f_asc, color='k', linewidth=1.5)

        # ## center-edge deflection (ce)
        w_ml_asc = -self.WA_ML[:max_force_idx + 1]
        w_mr_asc = -self.WA_MR[:max_force_idx + 1]
        # get the average displacement values of the corresponding displacement gauges
        w_mlmr_asc = (w_ml_asc + w_mr_asc) / 2
#        axes.plot(w_mlmr_asc, f_asc, color='red', linewidth=1)
        axes.plot(w_ml_asc, f_asc, color='k', linewidth=linewidth)
        axes.plot(w_mr_asc, f_asc, color='k', linewidth=linewidth)

        # ## edge deflections (e)
        w_v_asc = -self.WA_V[:max_force_idx + 1]
        w_h_asc = -self.WA_H[:max_force_idx + 1]
        w_l_asc = -self.WA_L[:max_force_idx + 1]
        w_r_asc = -self.WA_R[:max_force_idx + 1]
        axes.plot(w_v_asc, f_asc, color='grey', linewidth=linewidth)
        axes.plot(w_h_asc, f_asc, color='grey', linewidth=linewidth)
        axes.plot(w_l_asc, f_asc, color='k', linewidth=linewidth)
        axes.plot(w_r_asc, f_asc, color='k', linewidth=linewidth)

        # get the average displacement values of the corresponding displacement gauges
#        w_vh_asc = (w_v_asc + w_h_asc) / 2
#        w_lr_asc = (w_l_asc + w_r_asc) / 2
#        axes.plot(w_vh_asc, f_asc, color='green', linewidth=1, label='w_vh')
#        axes.plot(w_lr_asc, f_asc, color='blue', linewidth=1, label='w_lr')

        # save 'F-w-arr_m-mlmr-vh-lr' in directory "/simdb/simdata/exp_st"
#        simdata_dir = os.path.join(simdb.simdata_dir, 'exp_st')
#        if os.path.isdir(simdata_dir) == False:
#            os.makedirs(simdata_dir)
#        filename = os.path.join(simdata_dir, 'F-w-arr_m-mlmr-vh-lr_' + self.key + '.csv')
#        Fw_m_mlmr_vh_lr_arr = np.hstack([f_asc[:, None], w_m[:, None] - w_m[ idx_0 ], w_mlmr_asc[:, None], w_vh_asc[:, None], w_lr_asc[:, None]])
#        print 'Fw_m_mlmr_vh_lr_arr'
#        np.savetxt(filename, Fw_m_mlmr_vh_lr_arr, delimiter=';')
#        print 'F-w-curves for center, middle, edges saved to file %s' % (filename)

    #--------------------------------------------------------------------------------
    # view
    #--------------------------------------------------------------------------------

    traits_view = View(VGroup(
                         Group(
                              Item('jump_rtol', format_str="%.4f"),
                              label='curve_ironing'
                              ),
                         Group(
                              Item('thickness', format_str="%.3f"),
                              Item('edge_length', format_str="%.3f"),
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
    ex = ExRunClassExt(klass=ExpST)
    ex.configure_traits()

