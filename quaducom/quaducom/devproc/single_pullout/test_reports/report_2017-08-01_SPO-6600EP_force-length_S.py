'''
Created on October 10, 2017
by jbielak
'''

import os
import scipy
from matresdev.db.exdb import ExRun
from matresdev.db.simdb.simdb import simdb
import numpy as np
import pylab as p


params = {'legend.fontsize': 20,
          'font.size': 20,
          'font.family': 'Times New Roman',
          'font.style': 'normal'}
p.rcParams.update(params)

test_files = [
    'SPOS19-1.DAT',
    'SPOS19-2.DAT',
    'SPOS38-1.DAT',
    'SPOS38-2.DAT',

]

test_file_path = os.path.join(simdb.exdata_dir,
                              'single_pullout_tests',
                              '2017-08-03_SPO-1C-40mm-90-3,62EP')

e_array = np.array([ExRun(data_file=os.path.join(test_file_path, test_file))
                    for test_file in test_files]).reshape(-1, 1)

n_roving_array = np.array([1])
n_roving_array = np.repeat(n_roving_array, 4).reshape(-1, 1)

l_v_array = np.array([20.5, 22, 38.5, 39.5, ]).reshape(-1, 1)

color_list = ['r', 'r', 'k', 'k', ]
marker_list = ['^', '^', 'o', 'o', ]


def plot_all():

    fig = p.figure(
        facecolor='white', figsize=(30 / 2.54, 20 / 2.54), dpi=100)
    fig.suptitle(
        'Q95/95-CCE-38 // C3-HF2-165-4 90 Grad // 40mm Probendicke', fontsize=20)
    fig.subplots_adjust(
        left=0.1, right=0.96, bottom=0.1, top=0.93, wspace=0.25, hspace=0.2)

    for idx, (e_treatment, e_treatment_n_r, e_treatment_l_v) in enumerate(zip(e_array,
                                                                              n_roving_array, l_v_array)):
        for e_run, n_r, l_v in zip(e_treatment, e_treatment_n_r, e_treatment_l_v):
            e = e_run.ex_type
            axes = p.subplot(111)
            axes.plot(l_v, e.Kraft.max() / n_r,
                      marker=marker_list[idx], markersize=8, color=color_list[idx],
                      #label = test_files[idx].split('.')[0]
                      )
            print(e.Kraft.max())
# Print Ultimate tensile stress (average value) = 3221 N/mm2
    F_r = 3221 * 3.62 / 1000
    f_r = [F_r, F_r]
    lv_r = [0, 200]
    axes.plot(lv_r, f_r, color='k', linestyle='--', linewidth=2,
              label='Bruchkraft Faserstrang (Mittelwert)')

    axes.grid()
    axes.set_xlabel('Verankerungslaenge [mm]')
    axes.set_ylabel('Kraft je Faserstrang [kN]')

    axes.tick_params(axis='both', which='major', labelsize=16)
    first_legend = axes.legend(
        loc=2, markerscale=1., fontsize=20, numpoints=1,)
    newlegend = p.gca().add_artist(first_legend)

    markercircle, =  axes.plot(
        [], '^', color='r', markersize=12, label='Ohne Querroving')
    markertriangle, = axes.plot(
        [], 'o', color='k', markersize=12, label='1 Querroving')
#    markerdiamond, = axes.plot(
#        [], 'D', color='teal', markersize=12, label='2 Querrovings')
    handles = [markercircle, markertriangle, ]
    labels = [h.get_label() for h in handles]
    axes.legend(
        handles=handles, markerscale=1., fontsize=20, numpoints=1, loc=4)

    major_xticks = np.arange(0, 51, 10)
    major_yticks = np.arange(0, 16, 1)
    axes.axis([0., 50, 0., 15])
    axes.set_xticks(major_xticks)
    axes.set_yticks(major_yticks)

if __name__ == '__main__':
    plot_all()
    p.show()
