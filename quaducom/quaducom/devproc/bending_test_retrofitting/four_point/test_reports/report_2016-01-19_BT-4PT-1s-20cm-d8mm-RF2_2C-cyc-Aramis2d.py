'''
Created on Jan 28, 2015

'''
import os

from matresdev.db.exdb import ExRun
from matresdev.db.simdb import SimDB
import numpy as np
import pylab as p
simdb = SimDB()

import quaducom.devproc.bending_test_retrofitting.four_point.exp_bt_4pt_rf

# specify font options for plots
params = {'legend.fontsize': 12,
          u'font.size': 15,
          u'font.family': 'serif',
          u'font.style': 'normal'}
p.rcParams.update(params)

#--------------------
# select plot
# NOTE: 'do'-key is used for file name of saved image
#--------------------
do = 'F-w-center'  # force-displacement
# do = 'strains-top-bottom'

#--------------------
# four point bending test (test with retrofitting / monotonic loading)
#--------------------
test_file_path = [os.path.join(simdb.exdata_dir,
                              'bending_tests_retrofitting', 'four_point',
                              '2016-01-19_BT-4PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d')]

test_files = ['BT-4PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d.DAT']

color_list = ['blue']
linestyle_list = ['-']
plot_orig_list = [1]
# cutoff long label names at the end for cleaner legend display
label_cutoff = [-9]
#--------------------


#--------------------
# script fills up empty values (last available value of lines above) if:
# 1)if measuring error ocurred
# 2)if measuring frequency was different for displacement gauges and strain gauges
# script merges the first two columns to one
#--------------------

# data_file = os.path.join(test_file_path[0], test_files[0])
#
# file_split = data_file.split('.')
# file_ = open(file_split[0] + '.csv', 'r')
#
# data_file_str = file_.readlines()
#
# data_file_str[0] = data_file_str[0].replace("Datum;Zeit;", "Datum/Uhrzeit;")
# data_file_str[1] = data_file_str[1].replace(";;kN", ";kN")
# data_file_str = [data_file_str[i].replace(".2016;", ".2016 ") for i in range(len(data_file_str))]
#
# data_file_str_3_split = data_file_str[3].split(';')
# # print 'data_file_str_3_split', data_file_str_3_split
# for i in range(3, len(data_file_str)):
#     data_file_str_i_split = data_file_str[i].split(';')
# #     print 'data_file_str_i_split', data_file_str_i_split
# #     print 'len(data_file_str_i_split)', len(data_file_str_i_split)
#     for j in range(len(data_file_str_i_split)):
#        if data_file_str_i_split[j].strip() == '':
# #            print 'replace', j
#            # replace value in current line with value of previous line
#            data_file_str_i_split[j] = data_file_str[i - 1].split(';')[j]
# #            print 'data_file_str_i_split[j]', data_file_str_i_split[j]
#            data_file_str[i] = ''
#            for j in range(len(data_file_str_i_split) - 1):
#                data_file_str[i] += data_file_str_i_split[j] + ';'
#            data_file_str[i] += data_file_str_i_split[-1]
#        if len(data_file_str_3_split) > len(data_file_str_i_split):
# #            print 'len(data_file_str_i_split)', len(data_file_str_i_split)
#            # replace last value in current line with value of previous line
#            data_file_str[i] += data_file_str[i - 1].split(';')[j + 1]
# file_002 = open(file_split[0] + '_processed.csv', 'w')
# file_002 = file_002.writelines(data_file_str[:])
# test_files = ['BT-4PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d_processed.DAT']

e_list = [ExRun(data_file=os.path.join(test_file_path[i], test_files[i]))
          for i in range(len(test_files))]

def plot_all():

    fig = p.figure(facecolor='white', figsize=(8, 6))

    for idx, e_run in enumerate(e_list):

        e = e_run.ex_type

        axes = p.subplot(111)

        if do == 'F-w-center':
            e._plot_force_deflection_center_orig(axes, linewidth=1.5, color=color_list[
                                                 idx], label=e_list[idx].ex_type.key[0:label_cutoff[idx]])
            axes.set_xlabel('$w$ [mm]')
            axes.set_ylabel('$F$ [kN]')
            xlim = 140
            ylim = 100.
            axes.axis([-1., xlim, 0., ylim])

        if do == 'strains-top-bottom':
            e._plot_strain_top_bottom_force(axes, linewidth=1.5, color=color_list[
                                                idx], label=e_list[idx].ex_type.key[0:label_cutoff[idx]])

            axes.set_xlabel('strain $\epsilon$ [1E-3]')
            axes.set_ylabel('vertical load $F$ [kN]')
            axes.axis([-5., 60, 0., 40])

    axes.grid()
    axes.legend(loc=4)

    # --------------------------------
    # save figure
    # --------------------------------
    save_fig_to_file = True

    if save_fig_to_file:
        # create a report-subfolder with the name of the script (without file extension '.py')
        # and save it in the test-type-subfolders with the name and path as ex_type
        test_series_name = os.path.basename(__file__)[:-3]
        subfolder_list = __file__.split(os.path.sep)
        devproc_idx = np.where(np.array(subfolder_list) == 'devproc')[0]
        subfolder_path = subfolder_list[devproc_idx + 1:-2] + [test_series_name]
        test_series_dir = os.path.join(simdb.report_dir)
        for subfolder_name in subfolder_path:
            test_series_dir = os.path.join(test_series_dir, subfolder_name)

        if not os.path.exists(test_series_dir):
            os.makedirs(test_series_dir)
        filename = os.path.join(test_series_dir, do + '.png')
        p.savefig(filename)
        print 'figure saved to file %s' % (filename)

if __name__ == '__main__':
    plot_all()
    p.show()
