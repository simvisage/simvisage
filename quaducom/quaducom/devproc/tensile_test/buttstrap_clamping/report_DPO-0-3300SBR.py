'''
Created on Jan 28, 2015

'''

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

test_files = [
    'DPO-30cm-0-3300SBR-V1.DAT',
    'DPO-30cm-0-3300SBR-V2.DAT',
    'DPO-40cm-0-3300SBR-V2.DAT',
    'DPO-40cm-0-3300SBR-V3.DAT',
    'DPO-50cm-0-3300SBR-V1.DAT',
    'DPO-50cm-0-3300SBR-V2.DAT',
    'DPO-60cm-0-3300SBR-V1.DAT',
    'DPO-60cm-0-3300SBR-V2.DAT',
]

test_file_path = os.path.join(simdb.exdata_dir,
                              'double_pullout',
                              '2015-07-15_DPO-30cm-0-3300SBR')

e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
          for test_file in test_files]

color_list = [
    'r',
    'r',
    'g',
    'g',
    'b',
    'b',
    'k',
    'k',
]

linestyle_list = [
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
]


def plot_all():

    fig = p.figure(facecolor='white', figsize=(12, 9))
    fig.subplots_adjust(
        left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_list):
        e = e_run.ex_type

        axes = p.subplot(111)

        label = test_files[idx].split('.')[0]
        e._plot_force_displacement_asc(axes)
        axes.grid()
        axes.set_xlabel('$\Delta$ w [mm]')
        axes.set_ylabel('F [kN]')

    axes.legend(loc=2)


if __name__ == '__main__':
    plot_all()
    p.show()
