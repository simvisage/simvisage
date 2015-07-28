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
        e._plot_force_displacement_asc(axes, color=color_list[idx], linestyle=linestyle_list[idx], label=label)
        axes.grid()
        axes.set_xlabel('$\Delta$ w [mm]')
        axes.set_ylabel('F [kN]')
#        axes.set_xlim([0., 0.012])
#        axes.set_ylim([0., 2500])

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
