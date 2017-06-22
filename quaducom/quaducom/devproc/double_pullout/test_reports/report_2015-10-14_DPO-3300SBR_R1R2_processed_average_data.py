'''
Created on Jul 30, 2015

@author: Yingxiong
'''
from itertools import cycle
from matplotlib.widgets import CheckButtons
from os import listdir
import os
import os

import matplotlib.pyplot as p
from matresdev.db.simdb import SimDB
import numpy as np
simdb = SimDB()

colors = ['b', 'g', 'r', 'c', 'm',  # 'y', 'k',
          'hotpink', 'darkblue', 'mediumaquamarine', 'mediumturquoise', 'tan', 'lightcoral',
          'coral', 'indianred', 'palegreen', 'crimson', 'orange', 'burlywood', 'deepskyblue',
          'firebrick', 'sienna', 'sandybrown']

colorcycler = cycle(colors)

# specify color dictionary in order to use same color for the same test in all figures
#
color_dict = {  # R2:
                '30-V2g_R2_f': next(colorcycler), '30-V4_R2_f': next(colorcycler), '30-V1_R2': next(colorcycler),
                '30-V3_R2': next(colorcycler), '30-V5g_R2': next(colorcycler), '40-V4_R2': next(colorcycler),
                '40-V3_R2_f': next(colorcycler), '40-V1_R2': next(colorcycler), '40-V2_R2_f': next(colorcycler),
                '40-V5_R2_f': next(colorcycler), '50-V1_R2': next(colorcycler), '50-V6g_R2_f': next(colorcycler),
                '50-V2_R2': next(colorcycler), '50-V4_R2_f': next(colorcycler), '50-V3_R2_f': next(colorcycler),
                '50-V5_R2_f': next(colorcycler), '60-V1_R2': next(colorcycler), '60-V2_R2': next(colorcycler),
                '60-V3_R2': next(colorcycler), '70-V2_R2_f': next(colorcycler), '70-V3g_R2_f': next(colorcycler),
                '70-V1_R2_f': next(colorcycler), '80-V2_R2_f': next(colorcycler), '80-V1_R2_f': next(colorcycler),
                # R1:
                '60-V1_R1': next(colorcycler), '60-V2_R1': next(colorcycler), '50-V1_R1': next(colorcycler),
                '50-V2_R1': next(colorcycler), '40-V2_R1': next(colorcycler), '30-V1_R1': next(colorcycler),
                '30-V2_R1': next(colorcycler)
}

fig, ax = p.subplots()

d_set = {}
labels = []
# 'processed_averaged_data_R1+R2')
fpath = os.path.join(
    simdb.exdata_dir, 'double_pullout', '2015-10-14_DPO-15mm-0-3300SBR_R2', 'processed_averaged_data_R1+R2')

do = 'R2_f_ALL'
# do = 'R2_f_TREND-30-40-50-70-80'
# do = 'R2_f_TREND-30-40-50-80'
# do = 'R2_f_GLUE-30'
# do = 'R2_f_GLUE-50'
# do = 'R2_f_GLUE-70'
# do = 'R2_RATE-30'
# do = 'R2_RATE-40'
# do = 'R2_RATE-50'
# do = 'R2_s_ALL'
# do = 'R2_s_TREND-50-40-30'
# do = 'R2_s_TREND-60-50-40-30'
# do = 'R1_all'

if do == 'R2_f_ALL':
    # ALL 14 fast tests
    # DPO-15mm-0-3300SBR_R2 (all fast ones; including glue + 70cm)
    #
    plot_set = {
        '30-V2g_R2_f',  # glued
        '30-V4_R2_f',
        '40-V3_R2_f',
        '40-V2_R2_f',
        '40-V5_R2_f',
        '50-V6g_R2_f',  # glued
        '50-V4_R2_f',
        '50-V3_R2_f',
        '50-V5_R2_f',
        '70-V2_R2_f',  # predamage (?)
        '70-V3g_R2_f',
        '70-V1_R2_f',
        '80-V2_R2_f',
        '80-V1_R2_f',
    }

if do == 'R2_s_ALL':
    # ALL 8 tests
    # DPO-15mm-0-3300SBR_R2 (all slow ones)
    #
    plot_set = {
        '30-V1_R2',
        '30-V3_R2',
        '40-V4_R2',
        '40-V1_R2',
        '50-V1_R2',
        '50-V2_R2',
        '60-V1_R2',
        '60-V2_R2',
        '60-V3_R2',
    }

if do == 'R2_s_TREND-50-40-30':
    # DPO-15mm-0-3300SBR_R2 (TREND for slow ones; no DPO-60cm)
    #
    plot_set = {
        '30-V1_R2',
        '40-V4_R2',
        '50-V1_R2',
    }

if do == 'R2_s_TREND-60-50-40-30':
    # DPO-15mm-0-3300SBR_R2 (TREND for slow ones; with DPO-60cm)
    #
    plot_set = {
        '30-V1_R2',
        '40-V4_R2',
        '50-V1_R2',
        '60-V1_R2',
    }

if do == 'R2_f_TREND-30-40-50-80':
    # GENERELL TREND -f
    # DPO-15mm-0-3300SBR_R2 (pressed+fast; no 70cm)
    #
    plot_set = {
        '30-V4_R2_f',
        '40-V5_R2_f',
        '50-V5_R2_f',
        '80-V1_R2_f',
    }

if do == 'R2_f_TREND-30-40-50-70-80':
    # GENERELL TREND -f
    # DPO-15mm-0-3300SBR_R2 (pressed+fast; with 70cm)
    #
    plot_set = {
        '30-V4_R2_f',
        '40-V5_R2_f',
        '50-V5_R2_f',
        '70-V2_R2_f',
        '80-V1_R2_f',
    }

if do == 'R2_f_GLUE-30':
    # INFLUENCE GLUE
    # DPO-15mm-0-3300SBR_R2 (compare glued/pressed + only fast ones;DPO-30)
    #
    plot_set = {
        '30-V2g_R2_f',  # glued
        '30-V4_R2_f',
    }

if do == 'R2_f_GLUE-50':
    # INFLUENCE GLUE
    # DPO-15mm-0-3300SBR_R2 (compare glued/pressed + only fast ones;DPO-50)
    #
    plot_set = {
        '50-V6g_R2_f',  # glued
        '50-V4_R2_f',
        '50-V3_R2_f',
        '50-V5_R2_f',
    }

if do == 'R2_f_GLUE-70':
    # INFLUENCE GLUE
    # DPO-15mm-0-3300SBR_R2 (compare glued/pressed + only fast ones; only DPO-70)
    #
    plot_set = {
        '70-V2_R2_f',  # predamage (?)
        '70-V3g_R2_f',
        '70-V1_R2_f',
    }

if do == 'R2_RATE-30':
    # RATE DEPENDENCY for DPO-30 only
    # DPO-15mm-0-3300SBR_R2 (fast ones compared with slow ones; no glued; only 30cm)
    #
    plot_set = {
        '30-V4_R2_f',
        '30-V1_R2',
        '30-V3_R2',
    }

if do == 'R2_RATE-40':
    # RATE DEPENDENCY for DPO-40 only
    # DPO-15mm-0-3300SBR_R2 (fast ones compared with slow ones; no glued; only 40cm)
    #
    plot_set = {
        '40-V4_R2',
        '40-V1_R2',
        '40-V2_R2_f',
        '40-V3_R2_f',
        '40-V5_R2_f',
    }

if do == 'R2_RATE-50':
    # RATE DEPENDENCY for DPO-50 only
    # DPO-15mm-0-3300SBR_R2 (fast ones compared with slow ones; no glued; only 50cm)
    #
    plot_set = {
        '50-V1_R2',
        '50-V2_R2',
        '50-V4_R2_f',
        '50-V3_R2_f',
        '50-V5_R2_f',
    }


if do == 'R2_all':
    # DPO-15mm-0-3300SBR_R2
    #
    plot_set = {
        '30-V2g_R2_f',
        '30-V4_R2_f',

        '30-V1_R2',
        '30-V3_R2',
        '30-V5g_R2',

        '40-V4_R2',
        '40-V3_R2_f',
        '40-V1_R2',
        '40-V2_R2_f',
        '40-V5_R2_f',

        '50-V1_R2',
        '50-V6g_R2_f',
        '50-V2_R2',
        '50-V4_R2_f',
        '50-V3_R2_f',
        '50-V5_R2_f',

        '60-V1_R2',
        '60-V2_R2',
        '60-V3_R2',

        '70-V2_R2_f',
        '70-V3g_R2_f',
        '70-V1_R2_f',

        '80-V2_R2_f',
        '80-V1_R2_f',
    }

if do == 'R1_all':
    # DPO-8mm-0-3300SBR_R1
    #
    plot_set = {'30-V1_R1',
                '30-V2_R1',
                '40-V2_R1',
                '50-V1_R1',
                '50-V2_R1',
                '60-V1_R1',
                '60-V2_R1',
                }

for fname in sorted(listdir(fpath), reverse=True):

    print 'file', fpath, fname
    data = np.loadtxt(os.path.join(fpath, fname), delimiter=';')
    flabel = fname.replace('-0-3300SBR', '')
    flabel = flabel.replace('DPO-', '')
    flabel = flabel.replace('.txt', '')
    flabel = flabel.replace('cm', '')
    labels.append(flabel)
#     print 'flabel', flabel

    if flabel in plot_set:
        #         print 'flabel', flabel
        if 'g' in flabel:
            marker = 'x'
        else:
            marker = None

        if '_f' in flabel:
            ls = '--'
        else:
            ls = '-'

        color = next(colorcycler)
#         color = str(float(flabel[:2]) * 0.01 + 0.2)
        d_set[flabel], = ax.plot(data[0][data[0] < 25], data[1][data[0] < 25],
                                 linestyle=ls, lw=1.5, marker=marker, markevery=0.03, color=color_dict[flabel], label=flabel)

p.xlabel('crack opening [mm]')
p.ylabel('force [kN]')
p.xlim((-2, 25))
p.ylim((0, 35))
p.legend(loc='best', prop={'size': 11})

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


p.show()
