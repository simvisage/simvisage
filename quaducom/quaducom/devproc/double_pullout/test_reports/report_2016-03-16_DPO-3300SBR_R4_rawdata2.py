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

e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
          for test_file in test_files]

color_list = [
    'r', 'r', 'r', 'g', 'g', 'g', 'b', 'b', 'b',
    'k', 'k', 'k', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'firebrick',
    'firebrick', 'firebrick', 'darkblue', 'darkblue', 'darkblue'
]

linestyle_list = [
    '-',    '-',    '-',    '-',    '-',    '-',    '-',    '-',
    '-',    '-',    '-',    '-',    '-',    '-',    '-',    '-',
    '-',    '-',    '-',    '-',    '-',
]


n_roving_list = [
    9, 9, 9,
    8, 8, 8,
    8, 9, 9,
    8, 9, 9,
    8, 8, 8,
    9, 9, 9,
    9, 8, 8,
]

l_v_list = [
    50, 50, 50,
    100, 100, 100,
    150, 150, 150,
    200, 200, 200,
    250, 250, 250,
    300, 300, 300,
    350, 350, 350,
]

label_list = [
    0, 0, 1,
    0, 0, 1,
    0, 1, 0,
    0, 0, 1,
    0, 0, 1,
    0, 0, 1,
    0, 0, 1,
]

def plot_all():

    fig = p.figure(facecolor='white', figsize=(12, 9))
    fig.subplots_adjust(
        left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, (e_run, n_r) in enumerate(zip(e_list, n_roving_list)):
        e = e_run.ex_type
        axes = p.subplot(111)
        w_lim_idx_array = np.where(e.w > 10)[0]
        if len(w_lim_idx_array) == 0:
            w_lim_idx = len(e.w)
        else:
            w_lim_idx = w_lim_idx_array[0]
        if label_list[idx] == 1:
            axes.plot(e.w[:w_lim_idx], e.Kraft[:w_lim_idx] / n_r, linewidth=1.5,
                  linestyle=linestyle_list[idx], color=color_list[idx], label = 'length: %s mm' % (l_v_list[idx]))
        axes.grid()
        axes.set_xlabel('$\Delta$ w [mm]')
        axes.set_ylabel('Force per roving [kN]')
    handles, labels = axes.get_legend_handles_labels()
    axes.legend(handles[::-1], labels[::-1],loc=2)
    axes.axis([0., 12, 0., 2.5])

if __name__ == '__main__':
    plot_all()
    p.show()
