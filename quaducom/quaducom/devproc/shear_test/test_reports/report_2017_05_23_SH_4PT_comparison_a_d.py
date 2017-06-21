'''
Created on Jun 20, 2017

'''

import os

from matresdev.db.exdb import ExRun
from matresdev.db.simdb.simdb import simdb
import numpy as np
import pylab as p


params = {'legend.fontsize': 16,
          # 'legend.linewidth': 2
          }
p.rcParams.update(params)

data_name = '4-Punkt-Schubversuch'

test_files = [
    'SH4-9-V1.DAT',
    'SH4-9-V2.DAT',
    'SH4-9-V3.DAT',
    'SH475V1.DAT',
    'SH475V2.DAT',
    'SH475V3.DAT',
]

test_file_path = os.path.join(simdb.exdata_dir,
                              'shear_tests', 'four_point', '2017-05-23_SH-4PT-1G-3cm')

e_array = np.array([ExRun(data_file=os.path.join(test_file_path, test_file))
                    for test_file in test_files]).reshape(-1, 3)

# n_roving_array = np.array([
#     5, 5, 5,
#     5, 5, 5,
#     5, 5, 5,
#     5, 5, 5,
# ]).reshape(-1, 3)

marker_array = np.array([
    'o', '^', 's',
    's', 'o', 'o',
]).reshape(-1, 3)

color_array = np.array([
    'r', 'r', 'r',
    'r', 'r', 'r',
]).reshape(-1, 3)

ad_array = np.array([
    4.90, 5.00, 4.97,
    5.39, 5.56, 5.92,
]).reshape(-1, 3)


def plot_all():

    fig = p.figure(facecolor='white', figsize=(12, 9))
    fig.subplots_adjust(
        left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, (e_treatment, e_treatment_ad, e_treatment_marker, e_treatment_color) in enumerate(zip(e_array, ad_array, marker_array, color_array)):
        for e_run, ad, marker, color in zip(e_treatment, e_treatment_ad, e_treatment_marker, e_treatment_color):
            e = e_run.ex_type
            axes = p.subplot(111)
            axes.plot(ad, 0.5*e.Kraft.max(),
                      marker= marker, markersize=12, color= color,
                      #label = test_files[e_run].split('.')[0]
                      )

    axes.grid(b=True, which='major', color='gray', linestyle='-', linewidth = .5,)
    axes.set_xlabel('a/d [-]', fontsize = 16)
    axes.set_ylabel('Querkraft [kN]', fontsize = 16)
    axes.tick_params(axis='both', which='major', labelsize=16)
    
    
    # Create custom legend

    markertriangle, = axes.plot([],'^', color='k', markersize=12, label='Querkraftversagen')
    markercircle, =  axes.plot([],'o', color='k', markersize=12, label='Biegeversagen')
    markersquare, = axes.plot([],"s", color='k', markersize=12, label='uneindeutiges Versagen')
    handles = [markertriangle, markercircle, markersquare]
    labels = [h.get_label() for h in handles] 
    axes.legend(handles=handles, labels=labels, numpoints=1,)
    axes.axis([0., 8, 0., 10])

if __name__ == '__main__':
    plot_all()
    p.show()
