'''
Created on 22.10.2015

@author: Yingxiong
'''
from itertools import cycle
from matplotlib.widgets import CheckButtons
from os import listdir
import os
import os

import matplotlib.pyplot as plt
from matresdev.db.simdb.simdb import simdb
import numpy as np

lines = ["-", "-.", ":"]
linecycler = cycle(lines)

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorcycler = cycle(colors)

fig, ax = plt.subplots()
plt.xlim((-2, 25))
plt.ylim((0, 35))

rax = {}

# 1rd round
rax['30-V1_R1'] = plt.axes([0.05, 0.1, 0.06, 0.03])
rax['30-V2_R1'] = plt.axes([0.12, 0.1, 0.06, 0.03])
rax['40-V2_R1'] = plt.axes([0.19, 0.1, 0.08, 0.03])
rax['50-V1_R1'] = plt.axes([0.05, 0.14, 0.10, 0.03])
rax['50-V2_R1'] = plt.axes([0.16, 0.14, 0.1, 0.03])
rax['60-V1_R1'] = plt.axes([0.05, 0.18, 0.12, 0.03])
rax['60-V2_R1'] = plt.axes([0.18, 0.18, 0.12, 0.03])

# 2rd round - plate 1
rax['60-V1_R2'] = plt.axes([0.05, 0.23, 0.12, 0.03])
rax['60-V2_R2'] = plt.axes([0.17, 0.23, 0.12, 0.03])
rax['60-V3_R2'] = plt.axes([0.17, 0.26, 0.12, 0.03])
rax['30-V4_R2_f'] = plt.axes([0.05, 0.26, 0.06, 0.03])
rax['30-V5g_R2'] = plt.axes([0.11, 0.26, 0.06, 0.03])
rax['30-V3_R2'] = plt.axes([0.05, 0.29, 0.06, 0.03])
rax['30-V2g_R2_f'] = plt.axes([0.05, 0.32, 0.06, 0.03])
rax['30-V1_R2'] = plt.axes([0.05, 0.35, 0.06, 0.03])
rax['40-V1_R2'] = plt.axes([0.11, 0.35, 0.08, 0.03])
rax['40-V2_R2_f'] = plt.axes([0.11, 0.32, 0.08, 0.03])
rax['40-V3_R2_f'] = plt.axes([0.11, 0.29, 0.08, 0.03])
rax['50-V3_R2_f'] = plt.axes([0.19, 0.29, 0.10, 0.03])
rax['50-V2_R2'] = plt.axes([0.19, 0.32, 0.10, 0.03])
rax['50-V1_R2'] = plt.axes([0.19, 0.35, 0.10, 0.03])

# 2rd round - plate 2
rax['40-V5_R2_f'] = plt.axes([0.05, 0.43, 0.08, 0.03])
rax['40-V4_R2'] = plt.axes([0.05, 0.46, 0.08, 0.03])
rax['50-V6g_R2_f'] = plt.axes([0.05, 0.49, 0.10, 0.03])
rax['50-V5_R2_f'] = plt.axes([0.05, 0.52, 0.10, 0.03])
rax['50-V4_R2_f'] = plt.axes([0.05, 0.55, 0.10, 0.03])
rax['80-V2_R2_f'] = plt.axes([0.13, 0.43, 0.16, 0.03])
rax['80-V1_R2_f'] = plt.axes([0.13, 0.46, 0.16, 0.03])
rax['70-V3g_R2_f'] = plt.axes([0.15, 0.49, 0.14, 0.03])
rax['70-V2_R2_f'] = plt.axes([0.15, 0.52, 0.14, 0.03])
rax['70-V1_R2_f'] = plt.axes([0.15, 0.55, 0.14, 0.03])

# 3rd round
rax['50-V1g_R3_f'] = plt.axes([0.05, 0.61, 0.10, 0.03])
rax['50-V3_R3_f'] = plt.axes([0.15, 0.61, 0.10, 0.03])

rax['60-V2_R3_f'] = plt.axes([0.05, 0.64, 0.12, 0.03])
rax['60-V3_R3_f'] = plt.axes([0.05, 0.67, 0.12, 0.03])
rax['60-V1g_R3_f'] = plt.axes([0.05, 0.70, 0.12, 0.03])
rax['40-V2_R3_f'] = plt.axes([0.17, 0.64, 0.08, 0.03])
rax['40-V3_R3_f'] = plt.axes([0.17, 0.67, 0.08, 0.03])
rax['40-V1g_R3_f'] = plt.axes([0.17, 0.70, 0.08, 0.03])

rax['70-V2_R3_f'] = plt.axes([0.05, 0.73, 0.14, 0.03])
rax['70-V1g_R3_f'] = plt.axes([0.05, 0.79, 0.14, 0.03])
rax['30-V2_R3_f'] = plt.axes([0.19, 0.73, 0.06, 0.03])
rax['30-V3_R3_f'] = plt.axes([0.19, 0.76, 0.06, 0.03])
rax['30-V1g_R3_f'] = plt.axes([0.19, 0.79, 0.06, 0.03])

rax['80-V2_R3_f'] = plt.axes([0.05, 0.82, 0.16, 0.03])
rax['80-V1g_R3_f'] = plt.axes([0.05, 0.85, 0.16, 0.03])
rax['20-V2_R3_f'] = plt.axes([0.21, 0.82, 0.04, 0.03])
rax['20-V1_R3_f'] = plt.axes([0.21, 0.85, 0.04, 0.03])

d_set = {}
labels = []

# 'processed_averaged_data_R1+R2+R3')
fpath = os.path.join(
    simdb.exdata_dir, 'double_pullout_tests', '2015-12-14_DPO-15mm-0-3300SBR_R3', 'processed_averaged_data_R1+R2+R3')

for fname in listdir(fpath):

    data = np.loadtxt(os.path.join(fpath, fname), delimiter=';')
    flabel = fname.replace('-0-3300SBR', '')
    flabel = flabel.replace('DPO-', '')
    flabel = flabel.replace('.txt', '')
    flabel = flabel.replace('cm', '')
    flabel = flabel.replace('.asc', '')
    labels.append(flabel)
    print('flabel', flabel)

    if 'g' in flabel:
        marker = 'x'
    else:
        marker = None

    if 'R1' in flabel:
        alpha = 1.0
        color = str(float(flabel[:2]) * 0.01 + 0.2)
    else:
        alpha = 1.0
        color = next(colorcycler)

    if '_f' in flabel:
        ls = '--'
    else:
        ls = '-'

    d_set[flabel], = ax.plot(data[0][data[0] < 25], data[1][data[0] < 25],
                             linestyle=ls, lw=1.5, color=color, marker=marker, markevery=0.03, label=flabel, alpha=alpha, visible=False)
    rax[flabel].plot([0, 2], [0, 0], linestyle=ls, lw=1.5, color=color,
                     marker=marker, markevery=0.1, label=flabel, alpha=alpha)
    rax[flabel].set_ylim((-0.2, 1))
    rax[flabel].set_xlim((-2, 4))

print('d_set.keys()', list(d_set.keys()))
# plt.legend(loc=2, ncol=2)
plt.subplots_adjust(left=0.35)

checkers = []
for key in d_set:
    checkers.append(CheckButtons(rax[key], [key], [False]))


def func(label):
    d_set[label].set_visible(not d_set[label].get_visible())
    plt.draw()

for checker in checkers:
    checker.on_clicked(func)

plt.show()
