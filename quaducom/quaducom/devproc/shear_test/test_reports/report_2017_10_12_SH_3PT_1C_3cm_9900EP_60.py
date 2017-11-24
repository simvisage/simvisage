'''
Created on 19.06.2017

@author: janb
'''


import os

from matresdev.db.exdb import ExRun
from matresdev.db.simdb.simdb import simdb
import numpy as np
import pylab as p


# specify font options for plots
params = {'legend.fontsize': 12,
          #         'legend.linewidth': 2,
          u'font.size': 15,
          u'font.family': 'serif',
          u'font.style': 'normal'}
p.rcParams.update(params)
# print p.rcParams.keys()


test_file_path = os.path.join(simdb.exdata_dir,
                              'shear_tests', 'three_point', '2017-10-12_SH-3PT-1C-3cm')

test_files = ['SH3-6-V1.DAT',
              'SH3-6-V2.DAT',
              'SH3-6-V3.DAT', ]

color_list = ['red', 'blue', 'green', ]
linestyle_list = ['-', '-', '-', ]
plot_asc_list = [0, 0, 0, ]
e_array = [ExRun(data_file=os.path.join(test_file_path, test_file))
           for test_file in test_files]


def plot_all():

    fig = p.figure(facecolor='white', figsize=(6, 6))
#     fig.subplots_adjust(
# left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_array):

        e = e_run.ex_type

        axes = p.subplot(111)
        if plot_asc_list[idx]:
            e._plot_shearforce_deflection_asc(axes, color=color_list[idx], linewidth=1.5,
                                              label=test_files[idx].split('.')[0],)
        else:
            e._plot_shearforce_deflection(
                axes, color=color_list[idx], linewidth=1.5, label=test_files[idx].split('.')[0])
        print e.Kraft.max()
        axes.set_ylabel('Querkraft [KN]')
        axes.set_xlabel('Mitteldurchbiegung [mm]')
        xlim = 10
        ylim = 20.
        axes.axis([0., xlim, 0., ylim])

    axes.grid(
        b=True, which='major', color='gray', linestyle='-', linewidth=.5,)
    axes.legend(loc=2)
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
        filename = os.path.join(test_series_dir, '90.pdf')
        p.savefig(filename)
        print 'figure saved to file %s' % (filename)

if __name__ == '__main__':
    plot_all()
    p.show()
