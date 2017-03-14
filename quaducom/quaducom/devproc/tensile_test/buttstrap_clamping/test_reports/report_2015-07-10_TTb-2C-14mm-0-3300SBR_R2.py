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
               'TTb-2C-14mm-0-3300SBR-V1_R2.DAT',
               'TTb-2C-14mm-0-3300SBR-V1b_R2.DAT',
              'TTb-2C-14mm-0-3300SBR-V2_R2.DAT',
              'TTb-2C-14mm-0-3300SBR-V3_R2.DAT',
               'TTb-2C-14mm-0-3300SBR-V4_R2.DAT',
              'TTb-2C-14mm-0-3300SBR-V5_R2.DAT',
              ]
test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-07-10_TTb-2C-14mm-0-3300SBR_R2'
                              )
e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

color_list = [
              'r',
              'r',
              'g',
              'b',
              'k',
              'magenta',
              ]
linestyle_list = [
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

        e._plot_sigtex_eps(axes, color=color_list[idx], linestyle=linestyle_list[idx], label=test_files[idx], plot_analytical_stiffness_II=False)
        axes.grid()
        axes.set_xlabel('eps [-]')
        axes.set_ylabel('sig_tex [MPa]')
        axes.set_xlim([0., 0.012])
        axes.set_ylim([0., 2500])

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
#    get_sig_tex_max()

#    #--------
#    # V1-V1b: glue reloading path to the end of first loading path
#    #--------
#    # stresses of first loading path
#    sig_tex_asc_V1 = e_list[0].ex_type.sig_tex_asc
#    eps_asc_V1 = e_list[0].ex_type.eps_asc
#    # stresses of second loading path
#    sig_tex_asc_V1b = e_list[1].ex_type.sig_tex_asc
#    eps_asc_V1b = e_list[0].ex_type.eps_asc
#    # index where reloading path reaches max load level from first loading path
#    sig_tex_max_V1 = np.max(sig_tex_asc_V1)
#    idx_rl = np.where(sig_tex_asc_V1b >= sig_tex_max_V1)[0]
#    print 'idx_rl', idx_rl
#    # remaining values of reloading path with greater values then max of first loading path
#    sig_tex_rl = sig_tex_asc_V1b[idx_rl]
#    eps_asc_rl = eps_asc_V1b[idx_rl]
#    # shift strains to last strain value of 1st loading path
#    eps_asc_rl -= eps_asc_rl[0]
#    eps_asc_rl += eps_asc_V1[-1]
#    # glue reloading path to the end of first loading path
#    sig_tex_asc_V1_glued = np.hstack([sig_tex_asc_V1, sig_tex_rl])
#    eps_asc_V1_glued = np.hstack([eps_asc_V1, eps_asc_rl])
#    p.plot(eps_asc_V1_glued, sig_tex_asc_V1_glued)

#    #--------
#    # V2-V2b: glue reloading path to the end of first loading path
#    #--------
#    # stresses of first loading path
#    sig_tex_V2 = e_list[2].ex_type.sig_tex
#    eps_V2 = e_list[2].ex_type.eps
#    # index where reloading path reaches max load level from first loading path
#    idx_rl = np.where(0 > sig_tex_V2[:-10] - sig_tex_V2[10:])[0]
#    print 'idx_rl', idx_rl
#    # remaining values of reloading path with greater values then max of first loading path
#    sig_tex_rl = sig_tex_V2[idx_rl[0]:]
#    eps_asc_rl = eps_asc_V2[idx_rl[0]:]
#    # shift strains to last strain value of 1st loading path
#    eps_asc_rl -= eps_asc_rl[0]
#    eps_asc_rl += eps_V1[-1]
#    # glue reloading path to the end of first loading path
#    sig_tex_asc_V1_glued = np.hstack([sig_tex_asc_V1, sig_tex_rl])
#    eps_asc_V1_glued = np.hstack([eps_asc_V1, eps_asc_rl])
#    p.plot(eps_asc_V1_glued, sig_tex_asc_V1_glued)


    p.show()
