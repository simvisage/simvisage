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

data_name = '3-Punkt-Schubversuch'

test_files = [
    'SH3-60V1.DAT',
    'SH3-60V2.DAT',
    'SH3-60V3.DAT',
    'SH3-75V1.DAT',
    'SH3-75V2.DAT',
    'SH3-75V3.DAT',
    'SH3-9-V1.DAT',
    'SH3-9-V2.DAT',
    'SH3-9-V3.DAT',
    'SH3105V1.DAT',
    'SH3105V2.DAT',
    'SH3105V3.DAT',
]

test_file_path = os.path.join(simdb.exdata_dir,
                              'shear_tests', 'three_point', '2017-05-04_SH-3PT-1G-3cm')

e_array = np.array([ExRun(data_file=os.path.join(test_file_path, test_file))
                    for test_file in test_files]).reshape(-1, 3)

# n_roving_array = np.array([
#     5, 5, 5,
#     5, 5, 5,
#     5, 5, 5,
#     5, 5, 5,
# ]).reshape(-1, 3)

marker_array = np.array([
    's', '^', '^',
    '^', '^', '^',
    's', '^', 'o',
    '^', 's', '^',
]).reshape(-1, 3)

color_array = np.array([
    'b', 'b', 'b',
    'b', 'b', 'b',
    'b', 'b', 'b',
    'b', 'b', 'b',
]).reshape(-1, 3)

ad_array = np.array([
    3.53, 3.55, 3.58,
    4.38, 4.72, 4.63,
    5.12, 5.47, 5.45,
    6.91, 6.75, 6.33,
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
