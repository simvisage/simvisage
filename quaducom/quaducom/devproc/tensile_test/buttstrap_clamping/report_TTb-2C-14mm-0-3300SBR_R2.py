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
#               'TTb-2C-14mm-0-3300SBR-V1_R2.DAT',
#               'TTb-2C-14mm-0-3300SBR-V1b_R2.DAT',
#              'TTb-2C-14mm-0-3300SBR-V2_R2.DAT',
              'TTb-2C-14mm-0-3300SBR-V3_R2.DAT',
#               'TTb-2C-14mm-0-3300SBR-V4_R2.DAT',
              'TTb-2C-14mm-0-3300SBR-V5_R2.DAT',
              ]
test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-07-10_TTb-2C-14mm-0-3300SBR_R2'
                              )
e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

def plot_all():

    fig = p.figure()
    fig.subplots_adjust(
        left=0.03, right=0.97, bottom=0.04, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_list):
        e = e_run.ex_type

        axes = p.subplot(121)
        e._plot_sigc_eps(axes)
        axes.grid()
        axes.set_xlabel('eps [-]')
        axes.set_ylabel('sig_comp [MPa]')
        axes.legend(loc=2)
        axes.set_title(test_files)

        axes = p.subplot(122)
        e._plot_sigtex_eps(axes)
        axes.grid()
        axes.set_xlabel('eps [-]')
        axes.set_ylabel('sig_tex [MPa]')
        axes.legend(loc=2)


if __name__ == '__main__':
    plot_all()
    p.show()
