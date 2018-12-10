'''
Created on Jan 28, 2015

'''
from matresdev.db.simdb.simdb import simdb
from matresdev.db.exdb import ExRun
from matresdev.db.simdb.simdb import simdb
import numpy as np
import pylab as p
import os
import quaducom.devproc.bending_test_retrofitting.four_point.exp_bt_4pt_rf

# specify font options for plots
params = {'legend.fontsize': 12,
          'ps.fonttype': 42,
          'font.size':15,
          'font.family':'sans-serif',
          'font.style':'normal',
          'font.serif': 'Arial'}

p.rcParams.update(params)
# print p.rcParams.keys()

#--------------------
# select plot
# NOTE: 'do'-key is used for file name of saved image
#--------------------
do = 'F-w'  # force-displacement
# do = 'crack-opening-shear-cracks' 

#--------------------
# bending test three point (with retrofitting / cyclic loading)
#--------------------
test_file_path = [
    os.path.join(simdb.exdata_dir,
                 'bending_tests_retrofitting', 'three_point',
                 '2016-01-21_BT-3PT-1s-20cm-d8mm-RF1_ref-mon-Aramis2d'),
    os.path.join(simdb.exdata_dir,
                 'bending_tests_retrofitting', 'three_point',
                 '2016-01-22_BT-3PT-1s-20cm-d8mm-RF3_2C-mon-Aramis2d'),
    os.path.join(simdb.exdata_dir,
                 'bending_tests_retrofitting', 'three_point',
                 '2016-01-25_BT-3PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d'),
    os.path.join(simdb.exdata_dir,
                 'bending_tests_retrofitting', 'three_point',
                 '2016-01-28_BT-3PT-1s-20cm-d8mm-RF2_2C-mon-Aramis2d')
]

test_files = [
    'BT-3PT-1s-20cm-d8mm-RF1_ref-mon-Aramis2d.DAT',
    'BT-3PT-1s-20cm-d8mm-RF3_2C-mon-Aramis2d.DAT',
    'BT-3PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d.DAT',
    'BT-3PT-1s-20cm-d8mm-RF2_2C-mon-Aramis2d.DAT',
]

#--------------------
# format plot
#--------------------
color_list = ['grey', 'r', 'k', 'g']
linestyle_list = ['-', ':', '-', '--']
plot_orig_list = [1, 1, 1, 1]
# cutoff long label names at the end for cleaner legend display
label_cutoff = [-9, -9, -9, -9]
#--------------------

e_list = [ExRun(data_file=os.path.join(test_file_path[i], test_files[i]))
          for i in range(len(test_files))]


def plot_all():

    fig = p.figure(facecolor='white', figsize=(8, 6))

    for idx, e_run in enumerate(e_list):

        e = e_run.ex_type

        axes = p.subplot(111)

        if do == 'F-w':
            e._plot_force_deflection(axes,
                                     linewidth=1.5,
                                     linestyle=linestyle_list[idx],
                                     color=color_list[idx],
                                     label=e_list[idx].ex_type.key[0:label_cutoff[idx]])
            axes.set_xlabel('$w$ [mm]')
            axes.set_ylabel('$F$ [kN]')
            xlim = 15.
            ylim = 250
            axes.axis([0, xlim, 0., ylim])

    axes.grid()
    axes.legend(loc=4)

    # --------------------------------
    # save figure
    # --------------------------------
    save_fig_to_file = True

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
        filename = os.path.join(test_series_dir, do + '.eps')
        p.savefig(filename)
        print('figure saved to file %s' % (filename))

if __name__ == '__main__':
    plot_all()
    p.show()
