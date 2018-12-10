'''
Created on 01.12.2016

@author: janb
'''
# TO DO: Specify the cross sectional area of one roving within the
# fabric_layup, delete the 1.84 mm^2 used in the k_roh - calculation

import os

from matresdev.db.exdb import ExRun
from matresdev.db.simdb.simdb import simdb
import numpy as np
import pylab as p


# specify font options for plots
params = {'legend.fontsize': 12,
          #         'legend.linewidth': 2,
          'font.size': 12,
          'font.family': 'Times New Roman',
          'font.style': 'normal'}
p.rcParams.update(params)

test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-06-19_TTb-1C-30mm-0-6600EP_Posata')

test_files = ['DK1-0-1.DAT',
              'DK1-0-2.DAT',
              'DK1-0-3.DAT',
              'DK1-0-4.DAT',
              'DK1-0-5.DAT', ]

label_list = ['hinten Rechts', 'hinten Links', 'Vorne']
e_array = [ExRun(data_file=os.path.join(test_file_path, test_file))
           for test_file in test_files]


def plot_all():

    fig = p.figure(
        facecolor='white', figsize=(8.8 / 2.54, 8.8 / 2.54), dpi=100)
    fig.subplots_adjust(
        left=0.15, right=0.96, bottom=0.15, top=0.96, wspace=0.25, hspace=0.2)

    e = e_array[2].ex_type

    axes = p.subplot(111)

    axes.plot(e.W10_re[:e.max_stress_idx + 1], e.F_asc,
              color='b', linewidth=1.5, linestyle='-', label=label_list[0])
    axes.plot(e.W10_li[:e.max_stress_idx + 1], e.F_asc,
              color='r', linewidth=1.5, linestyle='-', label=label_list[1])
    axes.plot(e.W10_vo[:e.max_stress_idx + 1], e.F_asc,
              color='g', linewidth=1.5, linestyle='-', label=label_list[2])
    axes.plot(((e.W10_re[:e.max_stress_idx + 1] + e.W10_li[:e.max_stress_idx + 1]) / 2 + e.W10_vo[
              :e.max_stress_idx + 1]) / 2, e.F_asc, color='k', linewidth=2.5, linestyle='-', label='Mittelwert')
    axes.grid(
        b=True, which='major', color='gray', linestyle='-', linewidth=.5,)
    axes.set_xlabel('Verformung [mm]')
    axes.set_ylabel('Kraft [kN]')
    axes.axis([0., 6, 0., 35])
    axes.legend(loc=4, handletextpad=0.1)
    # --------------------------------
    # save figure
    # --------------------------------
    save_fig_to_file = False

    if save_fig_to_file:
        # create a report-subfolder with the name of the script (without file extension '.py')
        # and save it in the test-type-subfolders with the name and path as
        # ex_type
        test_series_name = os.path.basename(__file__)[:-3]
        subfolder_list = __file__.split(os.path.sep)
        devproc_idx = np.where(np.array(subfolder_list) == 'devproc')[0]
        subfolder_path = subfolder_list[
            devproc_idx + 1:-2] + [test_series_name]
        test_series_dir = os.path.join(simdb.report_dir)
        for subfolder_name in subfolder_path:
            test_series_dir = os.path.join(test_series_dir, subfolder_name)

        if not os.path.exists(test_series_dir):
            os.makedirs(test_series_dir)
        filename = os.path.join(test_series_dir, '.eps')
        p.savefig(filename)
        print('figure saved to file %s' % (filename))

if __name__ == '__main__':
    plot_all()
    p.show()
