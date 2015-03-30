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

test_files = ['TTb1-2C-9mm-0-800SBR-V1.DAT',
              'TTb1-2C-9mm-0-800SBR-V2.DAT',
              'TTb1-2C-9mm-0-800SBR-V3.DAT',
              'TTb2-2C-9mm-0-800SBR-V1.DAT',
              ]
test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-03-23_TTb-2C-9mm-0-800SBR_Ring1'
                              )
e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

def plot_all():

    fig = p.figure()
    fig.subplots_adjust(
        left=0.03, right=0.97, bottom=0.04, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_list):
        e = e_run.ex_type

        axes = p.subplot(221)

        e._plot_force_displacement_asc(axes)

        axes.grid()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('F [kN]')
        axes.legend(loc=2)

        axes = p.subplot(222)

        e._plot_sigc_eps(axes)
        axes.grid()
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('eps [-]')
        axes.legend(loc=2)
        axes.set_title(test_files)

        axes = p.subplot(223)
        e._plot_sigtex_eps(axes)
        axes.grid()
        axes.set_xlabel('eps [-]')
        axes.set_ylabel('sig_tex [MPa]')
        axes.legend(loc=2)

        axes = p.subplot(224)


if __name__ == '__main__':
    plot_all()
    p.show()
