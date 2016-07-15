'''
Created on Jan 28, 2015

'''

from matresdev.db.simdb.simdb import simdb
from matresdev.db.exdb import ExRun

import os
import numpy as np
import pylab as p
params = {'legend.fontsize': 10,
          # 'legend.linewidth': 2
          }
p.rcParams.update(params)

test_files = [
               'TTb-2C-14mm-0-3300EP-V1_stextz3mm.DAT',
               'TTb-2C-14mm-0-3300EP-V2_stextz3mm.DAT',
               'TTb-2C-14mm-0-3300EP-V1_stextz4mm.DAT',
               'TTb-2C-14mm-0-3300EP-V2_stextz4mm.DAT',
               'TTb-2C-14mm-0-3300EP-V1_stextz5mm.DAT',
               'TTb-2C-14mm-0-3300EP-V2_stextz5mm.DAT',
              ]

test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-07-14_TTb-2C-14mm-0-3300EP_stexz'
                              )

e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

# test_files = ['TTb-2C-14mm-0-3300EP-V1_cyc-Aramis2d.DAT']
#
# test_file_path = os.path.join(simdb.exdata_dir,
#                              'tensile_tests', 'buttstrap_clamping',
#                              '2015-08-11_TTb-2c-14mm-0-3300EP_cyc-Aramis2d'
#                              )
#
# e_list = e_list + [os.path.join(test_file_path, test_file)]

def plot_all():

    fig = p.figure(facecolor='white', figsize=(12, 9))
    fig.subplots_adjust(
        left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_list):
        e = e_run.ex_type

        axes = p.subplot(111)

        e._plot_sigtex_eps(axes, label=test_files[idx], plot_analytical_stiffness_II=False)
        axes.grid()
        axes.set_xlabel('eps [-]')
        axes.set_ylabel('sig_tex [MPa]')
#        axes.set_xlim([0., 0.012])
#        axes.set_ylim([0., 2500])

    # material stiffness carbon (E=245GPa)
    xarr = np.array([0., 0.010])
    yarr = np.array([0., 2450.])
    axes.plot(xarr, yarr, linestyle='--', color='grey', linewidth=1.5, label='E_tex = 245 GPa')

    axes.legend(loc=2)


def get_sig_tex_max():

    sig_tex_max_list = []
    for idx, e_run in enumerate(e_list):
        e = e_run.ex_type
        sig_tex_max = e.sig_tex_asc[-1]
        sig_tex_max_list = sig_tex_max_list + [sig_tex_max]
    sig_tex_max_arr = np.array(sig_tex_max_list)
    print 'sig_tex_max_list', sig_tex_max_list
    print 'sig_tex_max_arr', sig_tex_max_arr
    sig_tex_max_average = np.average(sig_tex_max_arr)
    print 'sig_tex_max_average', sig_tex_max_average


if __name__ == '__main__':
    plot_all()
    get_sig_tex_max()

    p.show()
