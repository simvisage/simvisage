'''
Created on Jan 28, 2015

'''

from exp_att_db import ExpATTDB
from matresdev.db.simdb import SimDB
simdb = SimDB()
import os
import numpy as np
import pylab as p

test_files = ['TTb-4c-2cm-0-TU-V%d.DAT' % i for i in range(1, 6)]

test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2013-12-01_TTb-4c-2cm-0-TU_Aramis2d_RR'
                              )

res_key = 'Xf15s1-Yf15s4'

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
        axes.grid()
        axes.set_ylim(0, 55)
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('F [kN]')
        axes.legend(loc=2)

        axes = p.subplot(232)

        axes.plot(e.time_asc, -e.eps_asc, color='darkblue', label='eps')
        axes.grid()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('eps [-]')
        axes.legend(loc=2)
        axes.set_title(test_files)

        axes = p.subplot(233)

        print '-----------------------------n_steps', e.aramis_info.number_of_steps
        print '-----------------------------n_steps_cut1', len(e.t_aramis)

        axes = p.subplot(234)

        axes = p.subplot(235)

        axes = p.subplot(236)


if __name__ == '__main__':
    fig = p.figure()
    plot_all()
    p.show()
