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
          u'font.size': 20,
          u'font.family': 'Times New Roman',
          u'font.style': 'normal'}
p.rcParams.update(params)

test_files = [
    'SPOK19A1.DAT',
    'SPOK19A2.DAT',
    'SPOK19A3.DAT',
    'SPOK19A4.DAT',
    'SPOK19A6.DAT',
    'SPOK19B1.DAT',
    'SPOK19B2.DAT',
    'SPOK19B4.DAT',
    'SPOK19B5.DAT',
    'SPOK38-1.DAT',
    'SPOK38-2.DAT',
    'SPOK38-3.DAT',
    'SPOK38-4.DAT',
    'SPOK38-5.DAT',
    'SPOK38-6.DAT',
    'SPOK76-1.DAT',
    'SPOK76-2.DAT',
    'SPOK76-3.DAT',
    'SPOK1521.DAT',
    'SPOK1522.DAT',
    'SPOK1523.DAT',
]

test_file_path = os.path.join(simdb.exdata_dir,
                              'single_pullout_tests',
                              '2017-08-01_SPO-1C-40mm-0-3,62EP')

e_array = np.array([ExRun(data_file=os.path.join(test_file_path, test_file))
                    for test_file in test_files]).reshape(-1, 1)

n_roving_array = np.array([1])
n_roving_array = np.repeat(n_roving_array, 21).reshape(-1, 1)

l_v_array = np.array([19, 18, 19, 20, 19.5, 19.5, 20,
                      19.5, 22, 37, 39, 37, 38, 37, 38.5, 78, 77, 77, 152, 153, 152, ]).reshape(-1, 1)

color_list = ['k', 'k', 'k', 'k', 'k', 'r', 'r', 'r', 'r', 'g', 'g', 'g', 'g',
              'g', 'g', 'darkblue', 'darkblue', 'darkblue', 'purple', 'purple', 'purple']
marker_list = ['o', 'o', 'o', 'o', 'o', '^', '^', '^', '^', 'o', 'o', 'o', 'o',
               'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']


def plot_all():

    fig = p.figure(
        facecolor='white', figsize=(30 / 2.54, 20 / 2.54), dpi=100)
    fig.suptitle(
        'Q95/95-CCE-38 // C3-HF2-165-4 // 40mm Probendicke', fontsize=20)
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

    markertriangle, = axes.plot(
        [], 'o', color='k', markersize=12, label='Mit Querroving')
    markercircle, =  axes.plot(
        [], '^', color='r', markersize=12, label='Ohne Querroving')
    handles = [markertriangle, markercircle]
    labels = [h.get_label() for h in handles]
    axes.legend(
        handles=handles, markerscale=1., fontsize=20, numpoints=1, loc=4)

    major_xticks = np.arange(0, 201, 20)
    major_yticks = np.arange(0, 16, 1)
    axes.axis([0., 200, 0., 15])
    axes.set_xticks(major_xticks)
    axes.set_yticks(major_yticks)

if __name__ == '__main__':
    plot_all()
    p.show()
