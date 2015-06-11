'''
Created on Jan 28, 2015

'''

from matresdev.db.simdb import SimDB
simdb = SimDB()
from matresdev.db.exdb import ExRun

import os
import numpy as np
import pylab as p
params = {'legend.fontsize': 10,
          # 'legend.linewidth': 2
          }
p.rcParams.update(params)

test_files = [('TTb1-2C-9mm-0-800SBR-V1_R1S1.DAT', 'clamp 1.0', 'series 1'),
              ('TTb1-2C-9mm-0-800SBR-V2_R1S1.DAT', 'clamp 1.0', 'series 1'),
              ('TTb1-2C-9mm-0-800SBR-V3_R1S1.DAT', 'clamp 1.0', 'series 1'),
              ('TTb2-2C-9mm-0-800SBR-V4_R1S1.DAT', 'clamp 2.0', 'series 1'),
              ('TTb3-2C-9mm-0-800SBR-V5_R1S1.DAT', 'clamp 3.0', 'series 1'),
              ('TTb3-2C-9mm-0-800SBR-V1_R1S2.DAT', 'clamp 3.0', 'series 2'),
              ('TTb3-2C-9mm-0-800SBR-V2_R1S2.DAT', 'clamp 3.0', 'series 2'),
              ('TTb3-2C-9mm-0-800SBR-V3_R1S2.DAT', 'clamp 3.0', 'series 2'),
              ('TTb1-2C-9mm-0-800SBR-V4_R1S2.DAT', 'clamp 1.0', 'series 2'),
              ('TTb1-2C-9mm-0-800SBR-V5_R1S2.DAT', 'clamp 1.0', 'series 2'),
              ]


clamp_colors = {'clamp 1.0': 'grey',
                'clamp 2.0': 'blue',
                'clamp 3.0': 'blue'}

series_colors = {'series 1': 'grey',
                 'series 2': 'red'}

series_styles = {'series 1': '--',
                 'series 2': '-'}

test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-03-23_TTb-2C-9mm-0-800SBR_Ring1'
                              )
e_list = [(ExRun(data_file=os.path.join(test_file_path, test_file)),
           clamp, series)
          for test_file, clamp, series in test_files]

sig_tex_max_arr = np.array([e_run.ex_type.sig_tex_max
                            for e_run, clamp, series in e_list], dtype='float_')

tf_array = np.array(test_files)
series_1_imap = np.where(tf_array[:, 2] == 'series 1')
series_2_imap = np.where(tf_array[:, 2] == 'series 2')
clamp_1_imap = np.where(tf_array[:, 1] == 'clamp 1.0')

clamp_23_map = np.logical_or(
    tf_array[:, 1] == 'clamp 2.0', tf_array[:, 1] == 'clamp 3.0-alpha')
clamp_23_imap = np.where(clamp_23_map)
print clamp_23_imap

print 'series 1', np.average(sig_tex_max_arr[series_1_imap]), np.std(sig_tex_max_arr[series_1_imap])
print 'series 2', np.average(sig_tex_max_arr[series_2_imap]), np.std(sig_tex_max_arr[series_2_imap])
print 'clamp - Jesse', np.average(sig_tex_max_arr[clamp_1_imap]), np.std(sig_tex_max_arr[clamp_1_imap])
print 'clamp - Scholzen', np.average(sig_tex_max_arr[clamp_23_imap]), np.std(sig_tex_max_arr[clamp_23_imap])


def plot_all():

    fig = p.figure(facecolor='white')
    fig.subplots_adjust(
        left=0.09, right=0.97, bottom=0.14, top=0.96, wspace=0.25, hspace=0.2)

    axes = p.subplot(111)

#    axes.plot([0.0, 0.01], [0.0, 1800], '--', color='gray')
#    axes.plot([0.0, 0.01], [0.0, 2450], '--', color='gray')

    for e_run, clamp, series in e_list:
        e = e_run.ex_type
        color_clamps = clamp_colors[clamp]
        style_series = series_styles[series]
        axes = p.subplot(111)

#        axes.plot(e.eps_asc, e.sig_tex_asc, style_series,
#                  color=color_clamps, linewidth=2)
        axes.plot(e.eps_smooth, e.sig_tex_smooth, style_series,
                  color=color_clamps, linewidth=2)

        axes.grid()
        axes.set_xlabel('eps [-]')
        axes.set_ylabel('sig_tex [MPa]')
        axes.legend(loc=2)

if __name__ == '__main__':
    plot_all()
    p.show()
