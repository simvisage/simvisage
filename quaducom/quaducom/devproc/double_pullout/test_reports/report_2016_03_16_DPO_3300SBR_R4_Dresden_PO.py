'''
Created on Mar 17, 2016

'''

import os
from matresdev.db.exdb import ExRun
from matresdev.db.simdb.simdb import simdb
import numpy as np
import pylab as p


params = {'legend.fontsize': 10,
          # 'legend.linewidth': 2
          }
p.rcParams.update(params)

test_files = [
    'B3-Ia-EV3-B-50-1.DAT',
    'B3-Ia-EV3-B-50-7.DAT',
    'B3-Ia-EV3-B-100-2.DAT',
    'B3-Ia-EV3-B-100-8.DAT',
    'B3-Ia-EV3-B-150-3.DAT',
    'B3-Ia-EV3-B-150-9.DAT',
    'B3-Ia-EV3-B-200-10.DAT',
    'B3-Ia-EV3-B-200-4.DAT',
    'B3-Ia-EV3-B-250-11.DAT',
    'B3-Ia-EV3-B-250-5.DAT',
    'B3-Ia-EV3-B-300-6.DAT',
    'B3-Ia-EV3-B-350-12.DAT',
]

test_file_path = os.path.join(simdb.exdata_dir,
                              'double_pullout_tests',
                              '2016-03-16_DPO-15mm-0-3300SBR_R4',
                              'Dresden_PO')

e_array = np.array([ExRun(data_file=os.path.join(test_file_path, test_file))
                    for test_file in test_files]).reshape(-1, 2)

color_list = [
    'r',  'g',  'b',  'k', 'mediumturquoise', 'firebrick', 'darkblue',
]

linestyle_list = [
    '-',    '-',    '-',    '-',    '-',    '-',
]


n_roving_array = np.array([
    9, 9,
    9, 9,
    9, 9,
    9, 9,
    9, 9,
    9, 9,
    9, 9,
]).reshape(-1, 2)

l_v_array = np.array([
    57, 52,
    100, 100,
    150, 150,
    200, 200,
    250, 250,
    300, 300,
    306, 368,
]).reshape(-1, 2)

gauge_dist = 65  # cm
left_gauge_name = 'IWA1'
right_gauge_name = 'IWA2'


def plot_all():

    fig = p.figure(facecolor='white', figsize=(12, 9))
    fig.subplots_adjust(
        left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, (e_treatment, e_treatment_n_r) in enumerate(zip(e_array, n_roving_array)):
        for e_run, n_r in zip(e_treatment, e_treatment_n_r):
            e = e_run.ex_type
            axes = p.subplot(111)
            axes.plot(e.w, e.Kraft / n_r, linewidth=1.5,
                      linestyle=linestyle_list[idx], color=color_list[idx])
            #label = test_files[idx].split('.')[0]
            #  e._plot_force_displacement_asc(axes)
    axes.grid()
    axes.set_xlabel('$\Delta$ w [mm]')
    axes.set_ylabel('force per roving [kN]')

    axes.legend(loc=2)
    axes.axis([0., 15, 0., 2.5])

if __name__ == '__main__':
    plot_all()
    p.show()
