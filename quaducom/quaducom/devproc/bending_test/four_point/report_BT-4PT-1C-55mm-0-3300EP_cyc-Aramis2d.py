'''
Created on Jan 28, 2015

'''
from matresdev.db.simdb import SimDB
simdb = SimDB()
from matresdev.db.exdb import ExRun

import os
import numpy as np
import pylab as p

# specify font options for plots
params = {'legend.fontsize': 12,
#         'legend.linewidth': 2,
          u'font.size':15,
          u'font.family':'serif',
          u'font.style':'normal'}
p.rcParams.update(params)
# print p.rcParams.keys()


test_file_path = os.path.join(simdb.exdata_dir,
                             'bending_tests', 'four_point',
                             '2015-09-02_BT-1C-55mm-0-3300SBR_cyc-Aramis2d')

#--------------------
# BT-4PT-1C-55cm-3300EP
#--------------------
test_files = ['BT-1C-55mm-0-3300EP-V2_S3P2(11)-Aramis2d.DAT',
              'BT-1C-55mm-0-3300EP-V2_S4P2(13)_cyc-Aramis2d.DAT']

n_rov_list = [9, 9]
color_list = ['grey', 'k']
linestyle_list = ['-', '-']
plot_orig_list = [0, 1]
label_cutoff = [-9, -9]  # cutoff long label names at the end for cleaner legend display
xlim = 160
ylim = 8.

#--------------------
# BT-4PT-1C-55cm-3300SBR
#--------------------
# test_files = ['BT-1C-55mm-0-3300SBR-V2_S4P2(14)-Aramis2d.DAT',
#               'BT-1C-55mm-0-3300SBR-V3_S2P1(12)-cyc-Aramis2d.DAT']
#
# n_rov_list = [15, 15]
# color_list = ['grey', 'k']
# linestyle_list = ['-', '-']
# plot_orig_list = [0, 1]
# label_cutoff = [-9, -9]  # cutoff long label names at the end for cleaner legend display
# xlim = 160
# ylim = 8.

#--------------------


e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

def plot_all():

    fig = p.figure(facecolor='white', figsize=(8, 6))

    for idx, e_run in enumerate(e_list):

        e = e_run.ex_type

        axes = p.subplot(111)

        if plot_orig_list[idx]:
            e._plot_force_deflection_center_orig(axes, linewidth=1.5, color=color_list[idx], label=e_list[idx].ex_type.key[0:label_cutoff[idx]])
        else:
            e._plot_force_deflection_center(axes, linewidth=1.5, color=color_list[idx], label=e_list[idx].ex_type.key[0:label_cutoff[idx]])

        axes.set_xlabel('$w$ [mm]')
        axes.set_ylabel('$F$ [kN]')
        axes.axis([0., xlim, 0., ylim])

    axes.grid()
    axes.legend(loc=4)

    # --------------------------------
    # save figure
    # --------------------------------
    save_fig_to_file = True
    test_series_name = 'BT_cyc-Aramis2d'
    if save_fig_to_file:
        img_dir = os.path.join(simdb.exdata_dir, 'img_dir')
        # check if directory exist otherwise create
        #
        if os.path.isdir(img_dir) == False:
            os.makedirs(img_dir)
        test_series_dir = os.path.join(img_dir, test_series_name)
        # check if directory exist otherwise create
        #
        if os.path.isdir(test_series_dir) == False:
            os.makedirs(test_series_dir)
        filename = os.path.join(test_series_dir, 'F-w.png')
        p.savefig(filename, format='png')
        print 'figure saved to file %s' % (filename)

if __name__ == '__main__':
    plot_all()
    p.show()
