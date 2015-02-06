'''
Created on Jan 28, 2015

'''

from exp_att_db import ExpATTDB, f_interp1d
from matresdev.db.simdb import SimDB
simdb = SimDB()
from matresdev.db.exdb import ExRun

import os
import numpy as np
import pylab as p
params = {'legend.fontsize': 10,
          # 'legend.linewidth': 2
          }
p.rcParams.update(params)

test_files = ['TTb-6c-2cm-0-TU-V%d.DAT' % i for i in range(1, 2)]
test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2013-12-02_TTb-6c-2cm-0-TU_Aramis2d_RR'
                              )
e_list_6c = [ExRun(data_file=os.path.join(test_file_path, test_file),
                   delta_t_aramis=5)
             for test_file in test_files]

res_key = 'Xf15s1-Yf15s4'


def plot_all():

    fig = p.figure()
    fig.subplots_adjust(
        left=0.03, right=0.97, bottom=0.04, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_list_6c):
        e = e_run.ex_type
        e.aramis_resolution_key = res_key
        # print 'crack filter', e.crack_filter_avg
        e.aramis_cdt.ddd_ux_avg_threshold = -5e-4
        e.aramis_cdt.crack_detection_step = 116
        e.process_aramisCDT_data()

        axes = p.subplot(231)

        e._plot_aramis_F_t_asc(axes)
        axes.plot(e.time_asc, e.F_asc, color='gray', label='F2')

        axes.grid()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('F [kN]')
        axes.legend(loc=2)

        axes = p.subplot(232)

        e._plot_aramis_t_strain_asc(axes)
        axes.plot(e.time_asc, -e.eps_asc, color='darkblue', label='eps')
        axes.grid()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('eps [-]')
        axes.legend(loc=2)
        axes.set_title(test_files)

        axes = p.subplot(233)
        axes.plot(-e.eps_asc, e.F_asc, color='darkblue', label='eps')
        e._plot_aramis_F_strain_asc(axes)
        axes.grid()
        axes.set_xlabel('eps [-]')
        axes.set_ylabel('F [kN]')
        axes.legend(loc=2)

        print '-----------------------------n_steps', e.aramis_info.number_of_steps
        print '-----------------------------n_steps_cut1', len(e.t_aramis)

        axes = p.subplot(234)

        e._plot_force_displacement(axes)

        axes = p.subplot(235)

        axes = p.subplot(236)


if __name__ == '__main__':
    plot_all()
    p.show()
