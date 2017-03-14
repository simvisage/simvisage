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
    'DPO-10cm-0-3300SBR-V1_R4.DAT',
    'DPO-10cm-0-3300SBR-V2_R4.DAT',
    'DPO-10cm-0-3300SBR-V3_R4.DAT',
    'DPO-20cm-0-3300SBR-V1_R4.DAT',
    'DPO-20cm-0-3300SBR-V2_R4.DAT',
    'DPO-20cm-0-3300SBR-V3_R4.DAT',
    'DPO-30cm-0-3300SBR-V1_R4.DAT',
    'DPO-30cm-0-3300SBR-V2_R4.DAT',
    'DPO-30cm-0-3300SBR-V3_R4.DAT',
    'DPO-40cm-0-3300SBR-V1_R4.DAT',
    'DPO-40cm-0-3300SBR-V2_R4.DAT',
    'DPO-40cm-0-3300SBR-V3_R4.DAT',
    'DPO-50cm-0-3300SBR-V1_R4.DAT',
    'DPO-50cm-0-3300SBR-V2_R4.DAT',
    'DPO-50cm-0-3300SBR-V3_R4.DAT',
    'DPO-60cm-0-3300SBR-V1_R4.DAT',
    'DPO-60cm-0-3300SBR-V2_R4.DAT',
    'DPO-60cm-0-3300SBR-V3_R4.DAT',
    'DPO-70cm-0-3300SBR-V1_R4.DAT',
    'DPO-70cm-0-3300SBR-V2_R4.DAT',
    'DPO-70cm-0-3300SBR-V3_R4.DAT',
]

test_file_path = os.path.join(simdb.exdata_dir,
                              'double_pullout_tests',
                              '2016-03-16_DPO-15mm-0-3300SBR_R4',
                              'raw_data')

e_array = np.array([ExRun(data_file=os.path.join(test_file_path, test_file))
                    for test_file in test_files]).reshape(-1, 3)

color_list = [
    'r',  'g',  'b',  'k', 'mediumturquoise', 'firebrick', 'darkblue',
]

linestyle_list = [
    '-',    '-',    '-',    '-',    '-',    '-',  '-',
]


n_roving_array = np.array([
    9, 9, 9,
    8, 8, 8,
    8, 9, 9,
    8, 9, 9,
    8, 8, 8,
    9, 9, 9,
    9, 8, 8,
]).reshape(-1, 3)

l_v_array = np.array([
    50, 50, 50,
    100, 100, 100,
    150, 150, 150,
    200, 200, 200,
    250, 250, 250,
    300, 300, 300,
    350, 350, 350,
]).reshape(-1, 3)

gauge_dist = 120
left_gauge_name = 'W10_li'
right_gauge_name = 'W10_re'
data_name = 'DPO Aachen'


def plot_all():

    fig = p.figure(facecolor='white', figsize=(12, 9))
    fig.subplots_adjust(
        left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, (e_treatment, e_treatment_n_r) in enumerate(zip(e_array,
                                                             n_roving_array)):
        for e_run, n_r in zip(e_treatment, e_treatment_n_r):
            e = e_run.ex_type
            axes = p.subplot(111)
            axes.plot(e.w, e.Kraft / n_r, linewidth=1.5,
                      linestyle=linestyle_list[idx], color=color_list[idx])
            #label = test_files[idx].split('.')[0]
            #  e._plot_force_displacement_asc(axes)
    axes.grid()
    axes.set_xlabel('$\Delta$ w [mm]')
    axes.set_ylabel('Kraft je Roving [kN]')

    axes.legend(loc=2)
    axes.axis([0., 15, 0., 2.5])

if __name__ == '__main__':
    plot_all()
    p.show()
