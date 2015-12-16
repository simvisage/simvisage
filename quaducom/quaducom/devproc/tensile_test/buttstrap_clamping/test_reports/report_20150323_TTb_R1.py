'''
Created on Jan 28, 2015

'''

from matresdev.db.simdb import SimDB
simdb = SimDB()
from matresdev.db.exdb import ExRun

import time

import os
import numpy as np
import pylab as p
params = {'legend.fontsize': 10,
          # 'legend.linewidth': 2
          }
p.rcParams.update(params)

test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-03-23_TTb_R1_all'
                              )

test_files = [
               'Aachen_TTb3-2C-9mm-0-800SBR-V2_R1S2.csv',
               'Aachen_TTb3-2C-9mm-0-800SBR-V3_R1S2.csv',
#               #
#               'Dresden_1C-DK1.csv',
#               'Dresden_1C-DK2.csv',
#               'Dresden_1C-DK3.csv',
#               'Dresden_1C-DK4.csv',
#               'Dresden_1C-DK5.csv',
#               #
#               'Leipzig_D1.csv',
#               'Leipzig_D2.csv',
#               'Leipzig_D3.csv',
#               'Leipzig_D4.csv',
#               'Leipzig_D5.csv',
               #
#               'Leipzig_G539.csv',  # A6
#               'Leipzig_G540.csv',  # A7
               'Leipzig_G546.csv',  # B6
               'Leipzig_G548.csv',  # B7
#               'Leipzig_G549.csv',  # C7
#               'Leipzig_G550.csv',  # C5
              ]


e_list = [os.path.join(test_file_path, test_file) for test_file in test_files]

fig = p.figure(facecolor='white', figsize=(18, 12))

A_tex = 0.446 * 24  # mm^2

sig_tex_max = 0.
n = 0
for e_run in e_list:
    if test_files[n].split('.')[0][0] == 'A':
        color = 'b'
    else:
        color = 'r'
    data_arr = np.loadtxt(e_run, delimiter=",", skiprows=1)
    F = data_arr[:, 0]
    idx_max = np.argmax(F)
    eps_asc = data_arr[:idx_max, 2]
    F_asc = data_arr[:idx_max, 0]
    sig_tex_asc = F_asc * 1000. / A_tex
    p.plot(eps_asc, sig_tex_asc, label=test_files[n].split('.')[0], linewidth=2, color=color)
    sig_tex_max = max(sig_tex_max, sig_tex_asc[-1])
    n += 1

E_tex = 245000.
eps_E_arr = np.array([0., sig_tex_max / E_tex]) * 1000.
sig_E_arr = np.array([0., sig_tex_max])
p.plot(eps_E_arr, sig_E_arr, color='grey', linestyle='--', linewidth=2, label='E_tex = 245 GPa')
E_tex = 200000.
eps_E_arr = np.array([0., sig_tex_max / E_tex]) * 1000. + 3
sig_E_arr = np.array([0., sig_tex_max])
p.plot(eps_E_arr, sig_E_arr, color='k', linestyle='--', linewidth=2, label='E_tex = 200 GPa')
p.xlabel('eps [1E-3]')
p.ylabel('sig_tex [MPa]')
p.legend(loc=2)
p.show()
