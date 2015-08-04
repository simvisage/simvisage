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
               'TTb-2C-14mm-0-3300SBR-V2_cyc-Aramis2d.DAT',
               'TTb-2C-14mm-0-3300SBR-V3_cyc-Aramis2d.DAT',
               'TTb-2C-14mm-0-3300SBR-V5_cyc-Aramis2d.DAT',
              ]

test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-08-03_TTb-2C-14mm-0-3300SBR_cyc-Aramis2d'
                              )

e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

color_list = [
              'r',
              'g',
              'b',
              ]
linestyle_list = [
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

#        sig_tex_asc = e.sig_tex_asc
#        WA1_hinten = np.copy(e.WA1_hinten)
#        WA2_links = np.copy(e.WA2_links)
#        WA3_rechts = np.copy(e.WA3_rechts)
#        min_WA1_hinten = np.min(WA1_hinten[:10])
#        min_WA2_links = np.min(WA2_links[:10])
#        min_WA3_rechts = np.min(WA3_rechts[:10])
#        WA1_hinten -= min_WA1_hinten
#        WA2_links -= min_WA2_links
#        WA3_rechts -= min_WA3_rechts
#        eps_hi = WA1_hinten / (e.gauge_length * 1000.)  # [mm/mm]
#        eps_li = WA2_links / (e.gauge_length * 1000.)  # [mm/mm]
#        eps_re = WA3_rechts / (e.gauge_length * 1000.)
#        if idx == 2:
#            eps_m = (eps_hi + eps_li + eps_re) / 2.
#        else:
#            eps_m = (eps_li + eps_re) / 2.
#        eps_asc = eps_m[:e.max_stress_idx + 1]
#        p.plot(eps_asc, sig_tex_asc, label=test_files[idx])

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
