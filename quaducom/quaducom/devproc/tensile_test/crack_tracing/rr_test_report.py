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
e_list_6c = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

test_files = ['TTb-4c-2cm-0-TU-V%d.DAT' % i for i in range(1, 2)]
test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2013-12-01_TTb-4c-2cm-0-TU_Aramis2d_RR'
                              )
e_list_4c = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

res_key = 'Xf15s1-Yf15s4'

# (integ_radius, dux_avg, ddd_ux_avg)
param_dict = {'TTb-4c-2cm-0-TU-V1.DAT': (22, -4e-4, 0.003),
              'TTb-4c-2cm-0-TU-V2.DAT': (20, -4e-4, 0.0025),
              'TTb-4c-2cm-0-TU-V3.DAT': (20, -4e-4, 0.0035),
              'TTb-4c-2cm-0-TU-V4.DAT': (22, -0.00028, 0.003),
              'TTb-4c-2cm-0-TU-V5.DAT': (20, -4e-4, 0.0035),
              'TTb-6c-2cm-0-TU-V1.DAT': (20, -2e-4, 0.002),
              'TTb-6c-2cm-0-TU-V2.DAT': (22, -2.5e-4, 0.0015),
              'TTb-6c-2cm-0-TU-V3.DAT': (20, -2e-4, 0.0025),
              'TTb-6c-2cm-0-TU-V4.DAT': (20, -2.5e-4, 0.0015),
              'TTb-6c-2cm-0-TU-V5.DAT': (20, -2.5e-4, 0.0015)
              }

e_list = e_list_4c


def plot_all():

    fig = p.figure()
    fig.subplots_adjust(
        left=0.03, right=0.97, bottom=0.04, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_list):
        e = e_run.ex_type
        e.aramis_resolution_key = res_key
        # print 'crack filter', e.crack_filter_avg
        params = param_dict.get(os.path.basename(e.data_file))
        e.aramis_cdt.aramis_data.integ_radius = params[0]
        e.aramis_cdt.ddd_ux_avg_threshold = params[1]
        e.aramis_cdt.ddd_ux_threshold = params[1]
        e.aramis_cdt.d_ux_avg_threshold = params[2]
        e.number_of_cracks_aramis

        axes = p.subplot(231)

        e._plot_aramis_F_t_asc(axes)
        axes.plot(e.time_asc, e.F_asc, color='gray', label='F2')

        axes.grid()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('F [kN]')
        axes.legend(loc=2)

        axes = p.subplot(232)

        e._plot_aramis_t_strain_asc(axes)
        axes.plot(e.time_asc, e.eps_asc, color='darkblue', label='eps')
        axes.grid()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('eps [-]')
        axes.legend(loc=2)
        axes.set_title(test_files)

        axes = p.subplot(233)
        axes.plot(e.eps_asc, e.F_asc, color='darkblue', label='eps')
        e._plot_aramis_Finit_strain(axes, color='g', marker='o', ls='None')
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
