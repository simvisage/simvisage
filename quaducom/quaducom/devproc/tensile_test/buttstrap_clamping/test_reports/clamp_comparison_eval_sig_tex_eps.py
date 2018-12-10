'''
Created on Jan 28, 2015

'''

from matresdev.db.simdb import SimDB
simdb = SimDB()
from matresdev.db.exdb import ExRun
import string
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
    tf_array[:, 1] == 'clamp 2.0', tf_array[:, 1] == 'clamp 3.0')
clamp_23_imap = np.where(clamp_23_map)
print(clamp_23_imap)

print('series 1', np.average(sig_tex_max_arr[series_1_imap]), np.std(sig_tex_max_arr[series_1_imap]))
print('series 2', np.average(sig_tex_max_arr[series_2_imap]), np.std(sig_tex_max_arr[series_2_imap]))
print('clamp - Jesse', np.average(sig_tex_max_arr[clamp_1_imap]), np.std(sig_tex_max_arr[clamp_1_imap]))
print('clamp - Scholzen', np.average(sig_tex_max_arr[clamp_23_imap]), np.std(sig_tex_max_arr[clamp_23_imap]))

header_template = '''===============================================================
IMB Tensile test report - %s
        
thickness: %g [m]    \twidth: %g [m]
gauge length: %g [m] \tA_tex: %g [mm2]
reinforcement ratio: %g [-]

textile fabrics: %s
number of reinforcement layers: %g
===============================================================
\teps [-]\t\t;\tsig_tex [MPa]
        '''


def eval_sig_tex_eps(smooth=False):

    for e_run, clamp, series in e_list:
        e = e_run.ex_type
        data_file = e_run.data_file
        dir_name = os.path.dirname(data_file)
        file_name = os.path.basename(data_file)
        base_name, ext = os.path.splitext(file_name)

        if smooth:
            out_arr = np.vstack([e.eps_smooth, e.sig_tex_smooth]).T
            out_file = os.path.join(
                dir_name, base_name + '-eps-sigtex_smooth.csv')
        else:
            out_arr = np.vstack([e.eps_asc, e.sig_tex_asc]).T
            out_file = os.path.join(
                dir_name, base_name + '-eps-sigtex_raw.csv')

        print('writing output to', out_file)

        mid_layer = e.ccs.fabric_layup_list[1]
        fabric_type = mid_layer.fabric_layout_key
        n_fabric_layers = mid_layer.n_layers

        header_string = header_template % (base_name, e.ccs.thickness, e.width,
                                           e.gauge_length, e.A_tex,
                                           e.rho_c, fabric_type,
                                           n_fabric_layers)
        np.savetxt(out_file,
                   out_arr, delimiter=';',
                   header=header_string)

if __name__ == '__main__':
    eval_sig_tex_eps(smooth=True)
