'''
Created on Jan 28, 2015

'''

from exp_att_db import ExpATTDB, f_interp1d
from matresdev.db.simdb import SimDB
simdb = SimDB()
import os
import numpy as np
import pylab as p
params = {'legend.fontsize': 10,
          #'legend.linewidth': 2
          }
p.rcParams.update(params)

test_files = ['TTb-4c-2cm-0-TU-V%d.DAT' % i for i in range(1, 3)]
test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2013-12-01_TTb-4c-2cm-0-TU_Aramis2d_RR'
                              )

# test_files = ['TTb-6c-2cm-0-TU-V%d.DAT' % i for i in range(1, 2)]
# test_file_path = os.path.join(simdb.exdata_dir,
#                               'tensile_tests', 'buttstrap_clamping',
#                               '2013-12-02_TTb-6c-2cm-0-TU_Aramis2d_RR'
#                               )
#
res_key = 'Xf15s1-Yf15s4'

# test_files = ['TTb-4c-2cm-0-TU-V%d_bs4.DAT' % i for i in range(1, 2)]
#
# test_file_path = os.path.join(simdb.exdata_dir,
#                               'tensile_tests', 'buttstrap_clamping',
#                               '2013-07-09_TTb-4c-2cm-0-TU_bs4-Aramis3d'
#                               )
#
# res_key = 'Xf19s1-Yf19s4'

e_list = [ExpATTDB(data_file=os.path.join(test_file_path, test_file),
                   aramis_resolution_key=res_key,
                   delta_t_aramis=5)
          for test_file in test_files]


def plot_all():

    fig = p.figure()
    fig.subplots_adjust(
        left=0.03, right=0.97, bottom=0.04, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e in enumerate(e_list):
        e.process_source_data()
        print 'step_times', e.aramis_field_data.step_times
        # print 'crack filter', e.crack_filter_avg
        e.process_source_data()
        e.aramis_cdt.ddd_ux_avg_threshold = -5e-4
        e.aramis_cdt.crack_detection_step = 116

        axes = p.subplot(231)

        axes.plot(e.t_aramis_cut, e.F_t_aramis, color='darkblue', label='F')
        axes.plot(e.time_asc, e.F_asc, color='gray', label='F2')

        x = e.aramis_cdt.aramis_data.step_times + e.aramis_start_time
        stress = e.aramis_cdt.aramis_data.ad_channels_arr[:, 1]
        axes.plot(x, stress, color='red')
        axes.grid()
        #axes.set_ylim(0, 55)
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('F [kN]')
        axes.legend(loc=2)

        axes = p.subplot(232)

        e.aramis_cdt.run_t = True
        idx = np.argwhere(e.aramis_cdt.force == np.max(e.aramis_cdt.force))
        axes.plot(e.aramis_cdt.aramis_data.step_times[:idx + 1],
                  e.aramis_cdt.control_strain_t[:idx + 1], color='lightblue', label='aramis eps')
        axes.plot(e.time_asc, -e.eps_asc, color='darkblue', label='eps')
        axes.grid()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('eps [-]')
        axes.legend(loc=2)
        axes.set_title(test_files)

        axes = p.subplot(233)
        axes.plot(-e.eps_asc, e.F_asc, color='darkblue', label='eps')
        axes.plot(e.aramis_cdt.control_strain_t[:idx + 1], e.aramis_cdt.force[:idx + 1],
                  color='lightblue', label='aramis eps')
        axes.plot(e.aramis_cdt.control_strain_t[:idx + 1], f_interp1d(e.aramis_cdt.aramis_data.step_times[:idx + 1] + e.aramis_start_time, e.time_asc, e.F_asc),
                  color='red', label='aramis eps')
        axes.grid()
        axes.set_xlabel('eps [-]')
        axes.set_ylabel('F [kN]')
        axes.legend(loc=2)

        print '-----------------------------n_steps', e.aramis_info.number_of_steps
        print '-----------------------------n_steps_cut1', len(e.t_aramis)

        axes = p.subplot(234)

        axes = p.subplot(235)

        axes = p.subplot(236)


if __name__ == '__main__':
    plot_all()
    p.show()
