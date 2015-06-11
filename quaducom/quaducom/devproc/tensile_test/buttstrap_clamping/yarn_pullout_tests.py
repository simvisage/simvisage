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
#              'YPO-CAR800SBR-V2.DAT',
              'YPO-CAR800SBR-V3.DAT',
#              'YPO-CAR800SBR-V4.DAT',
              'YPO-CAR800SBR-V5.DAT'
              ]

test_file_path = os.path.join(simdb.exdata_dir,
                              'yarn_pullout',
                              '2015-04-13_YPO-CAR800SBR_R1'
                              )

e_list = [ExRun(data_file=os.path.join(test_file_path, test_file)) for test_file in test_files]

fig = p.figure(facecolor='white')
fig.subplots_adjust(left=0.09, right=0.97, bottom=0.14, top=0.96, wspace=0.25, hspace=0.2)
axes = p.subplot(111)

n = 1
for e_run in e_list:
    e = e_run.ex_type
    e_name = e_run.ex_type.key
    axes = p.subplot(111)
    F = e.Kraft * 1000  # force [N]
    w_r = (e.W10_li + e.W10_re) / 2.  # average crack opening
    axes.plot(w_r, F, linewidth=2, label=e_name)
#    if n == 1:
#        color = 'red'
#    if n == 2:
#        color = 'grey'
#    axes.plot(e.W10_li, F, linewidth=2, label=e_name, color=color)
#    axes.plot(e.W10_re, F, linewidth=2, label=e_name, color=color)
#    n += 1
    axes.grid()
    axes.set_xlabel('w_r [mm]')
    axes.set_ylabel('F [N]')
    axes.legend(loc=1)
    axes.axis([-0.5, 10, 0., 250])

p.show()
