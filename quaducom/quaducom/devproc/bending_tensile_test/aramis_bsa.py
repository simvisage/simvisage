#-------------------------------------------------------------------------------
#
# Copyright (c) 2012
# IMB, RWTH Aachen University,
# ISM, Brno University of Technology
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in the Spirrid top directory "licence.txt" and may be
# redistributed only under the conditions described in the aforementioned
# license.
#
# Thanks for using Simvisage open source!
#
#-------------------------------------------------------------------------------

from etsproxy.traits.api import \
    Int, Float, \
    Property, cached_property, Array

import numpy as np
import os

import platform
import time
if platform.system() == 'Linux':
    sysclock = time.time
elif platform.system() == 'Windows':
    sysclock = time.clock

from aramis_cdt import AramisInfo, AramisData, AramisCDT

def get_d(u_arr, integ_radius):
    ir = integ_radius
    du_arr = np.zeros_like(u_arr)
    du_arr[:, ir:-ir] = (u_arr[:, 2 * ir:] - u_arr[:, :-2 * ir])
    return du_arr


class AramisBSA(AramisCDT):
    '''Crack Detection Tool for detection of cracks, etc. from Aramis data.
    '''

    d_ux_arr2 = Property(Array, depends_on='aramis_info_changed, aramis_data.+params_changed')
    '''The first derivative of displacement in x-direction
    '''
    @cached_property
    def _get_d_ux_arr2(self):
        x_und = self.aramis_data.x_arr_undeformed
        # print x_und
        # print np.ndim (x_und)
        # print np.shape (x_und)
        ux = self.aramis_data.ux_arr
        # print ux
        # print np.ndim (ux)
        # print np.shape (ux)
        du_arr = np.zeros_like(x_und)
        ir = self.integ_radius
        du_arr[:, ir:-ir] = (ux[:, 2 * ir:] - ux[:, :-2 * ir]) / (x_und[:, 2 * ir:] - x_und[:, :-2 * ir])
        # print 'ux'
        # print ux[:, :]
        # print 'ux[:, 2 * ir:]'
        # print ux[:, 2 * ir:]
        # print 'ux[:, :-2 * ir]'
        # print ux[:, :-2 * ir]
        # print 'x_und'
        # print x_und[:, :]
        # print 'x_und[:, 2 * ir:]'
        # print x_und[:, 2 * ir:]
        # print 'x_und[:, :-2 * ir]'
        # print x_und[:, :-2 * ir]
        # print 'du_arr[:, ir:-ir]'
        # print du_arr[:, ir:-ir]
        return du_arr

    #--------------------------------------------------------------------------------
    # get max tensile strain in the reinforcement layer
    #--------------------------------------------------------------------------------

    h_top_threshold = Float(1.5, auto_set=False, enter_set=True)
    '''Threshold for distance between specimen edge and aramis-mask on top side.
    '''
    h_bot_threshold = Float(1.0, auto_set=False, enter_set=True)
    '''Threshold for distance between specimen edge and aramis-mask on bottom side.
    '''
    h_re_6_threshold = Float(2.86, auto_set=False, enter_set=True)
    '''Threshold for position of first reinforcement layer (6 layers).
    '''
    h_re_4_threshold = Float(4.0, auto_set=False, enter_set=True)
    '''Threshold for position of first reinforcement layer (4 layers).
    '''

    d_ux_max_ten = Property(Array, depends_on='aramis_info_changed, aramis_data.+params_changed')
    @cached_property
    def _get_d_ux_max_ten(self):

        # get number of values in y-direction
        num_ele = self.d_ux_arr2.shape[0]

        # size of facet in y direction
        h_mask = 20.0 - self.h_top_threshold - self.h_bot_threshold
        size_fa_y = h_mask / num_ele

        # get indices of strain before strain of reinforcement layer
        pos_re_6 = self.h_re1_6_threshold - self.h_top_threshold
        pos_re_4 = self.h_re1_4_threshold - self.h_top_threshold

        idx_6 = pos_re_6 / size_fa_y
        idx_6 = round (idx_6)
        idx_4 = pos_re_4 / size_fa_y
        idx_4 = round (idx_4)

        # position of first node in y-direction
        pos_no_f = 20. - self.h_top_threshold - (size_fa_y / 2)
        # position of last node in y-direction
        pos_no_l = (size_fa_y / 2) + self.h_bot_threshold

        # get all values for x- axis and y-axis
        mid_idx = self.d_ux_arr2.shape[1] / 2
        x = np.mean(self.d_ux_arr2[:, mid_idx - 0:mid_idx + 1], axis=1)
        y = np.linspace(pos_no_f, pos_no_l, num=25)

        if '6layers':
            x1 = x[idx_6]
            y1 = y[idx_6]
            y_re = 20. - self.h_re_6_threshold
            x_re = (x1 * y_re) / y1
            return x_re
        else:
            x1 = x[idx_4]
            y1 = y[idx_4]
            y_re = 20. - self.h_re_4_threshold
            x_re = (x1 * y_re) / y1
            return x_re

    #--------------------------------------------------------------------------------
    # get heigth of compression zone
    #--------------------------------------------------------------------------------

    x_cz = Property(Array, depends_on='aramis_info_changed, aramis_data.+params_changed')

    @cached_property
    def _get_x_cz(self):

        # get indices of the strain before and after zero on x-axis
        x_idx = np.where (self.x < 0)
        x1_idx = x_idx[0]
        x2_idx = x1_idx - 1

        x1 = self.x[x1_idx]
        x2 = self.x[x2_idx]
        y2 = self.y[x2_idx]

        x_cz = (y2 * x1) / (x2 + x1)

        print 'x_cz'
        print x_cz
        return x_cz

if __name__ == '__main__':
    from os.path import expanduser
    home = expanduser("~")

    data_dir = os.path.join(home, '.simdb_cache', 'exdata/bending_tests',
                             'three_point', '2013-07-09_BT-6c-2cm-0-TU_bs4-Aramis3d',
                             'aramis', 'BT-6c-V4-bs4-Xf19s15-Yf19s15')

    data_dir = os.path.join(home, '.simdb_cache', 'exdata',
                             'bending_tensile_test', '2014-06-12_BTT-4c-2cm-0-TU_MxN2',
                             'aramis', 'A2d_BTT-4c-2cm-TU-0-V02_MxN2-Xf15s3-Yf15s3')

    AI = AramisInfo(data_dir=data_dir)
    # AI = AramisInfo()
    AD = AramisData(aramis_info=AI,
                    evaluated_step_idx=50)
    AC = AramisBSA(aramis_info=AI,
                   aramis_data=AD,
                   integ_radius=10)

    for step in range(0, 85, 5):  # [225]:
        AD.evaluated_step_idx = step
        mid_idx = AC.d_ux_arr2.shape[1] / 2
        x = np.mean(AC.d_ux_arr2[:, mid_idx - 0:mid_idx + 1], axis=1)
        # print x
        # x = AC.d_ux_arr2[:, mid_idx - 4:mid_idx + 4]

        # y = AD.y_arr_undeformed[:, mid_idx ]

        c = AC.d_ux_arr2.shape[0]
        b = 20. / (c + 1)
        a = 20. - b
        y = np.linspace(a, b, num=c)

        import matplotlib.pyplot as plt
        plt.plot(x, y)

#     plt.plot([-0.004, 0.014], [0, 0], color='black')
#     plt.plot([0, 0], [-10, 10], color='black')
    plt.show()
    # AC.configure_traits()
