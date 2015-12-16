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

test_files = [
#              'TTb-2C-14mm-0-3300SBR-V1_R2.DAT',
              'TTb-2C-14mm-0-3300SBR-V1b_R2.DAT',
              'TTb-2C-14mm-0-3300SBR-V2_R2.DAT',
              'TTb-2C-14mm-0-3300SBR-V3_R2.DAT',
              'TTb-2C-14mm-0-3300SBR-V4_R2.DAT',
              'TTb-2C-14mm-0-3300SBR-V5_R2.DAT',
              ]

test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-07-10_TTb-2C-14mm-0-3300SBR_R2')

e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
          for test_file in test_files]

thickness_list = [
#                  0.0142,  # [m]
                  0.0142,  # [m]
                  0.0145,
                  0.0145,
                  0.0146,
                  0.0144,
                  ]

width_list = [
#               0.098,  # [m]]
               0.098,  # [m]]
               0.100,
               0.0996,
               0.0988,
               0.0999,
             ]

n_rov_list = [
#               16,  # [m]]
               16,  # [m]]
               16,
               16,
               16,
               16,
             ]

color_list = [
#              'r',
              'r',
              'g',
              'b',
              'k',
              'magenta',
              ]

linestyle_list = [
#                  '-',
                  '-',
                  '-',
                  '-',
                  '-',
                  '-',
                  ]

def plot_eps_sigtex():
    fig = p.figure(facecolor='white', figsize=(12, 9))
    axes = p.subplot(111)
    for idx, e_run in enumerate(e_list):
        e = e_run.ex_type
        e._plot_sigtex_eps(axes, color=color_list[idx], linestyle=linestyle_list[idx], label=test_files[idx], plot_analytical_stiffness_II=False)
        axes.grid()
        axes.set_xlabel('eps [-]')
        axes.set_ylabel('sig_tex [MPa]')
        axes.legend(loc=2)
        axes.set_xlim([0., 0.012])
        axes.set_ylim([0., 2500])

#    p.savefig('eps-sigtex.pdf', format='pdf')
    p.show()

# sig_tex_max_arr = np.array([e_run.ex_type.sig_tex_max
#                            for e_run in e_list], dtype='float_')
# print 'sigtex (average) = %g' % (np.average(sig_tex_max_arr))
# print 'sigtex (standard deviation) = %g' % (np.std(sig_tex_max_arr))

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

    # plot
    fig = p.figure(facecolor='white', figsize=(12, 9))
    axes = p.subplot(111)

    sigtex_max_list = []
    for idx, e_run in enumerate(e_list):
        e = e_run.ex_type
        data_file = e_run.data_file
        dir_name = os.path.dirname(data_file)
        file_name = os.path.basename(data_file)
        base_name, ext = os.path.splitext(file_name)

        # use measured values for number of rovings and thickness of the cross section
        #
        A_rov = 1.84  # [mm^2]
        n_rov = n_rov_list[idx]  # [-]
        A_tex = n_rov * A_rov  # [mm^2]
        thickness = thickness_list[idx]  # [m]
        width = width_list[idx]  # [m]
        A_c = thickness * width
        rho_c = A_tex / A_c / 1000000.
        print 'A_rov', A_rov
        print 'n_rov', n_rov
        print 'A_tex', A_tex
        print 'thickness', thickness
        print 'width', width
        print 'A_c', A_c
        print 'rho_c', rho_c

        mid_layer = e.ccs.fabric_layup_list[1]
        fabric_type = mid_layer.fabric_layout_key
        flo_ref = mid_layer.fabric_layout_ref
        n_fabric_layers = mid_layer.n_layers

        a_tex_db = e.ccs.a_tex
        print 'a_tex_db', a_tex_db

        width_db = e.width
        print 'width_db', width_db

        A_tex_db = a_tex_db * width_db
        print 'A_tex_db', A_tex_db

        scale_factor = A_tex_db / A_tex

        if smooth:
            out_arr = np.vstack([e.eps_smooth, scale_factor * e.sig_tex_smooth]).T
            out_file = os.path.join(
                dir_name, base_name + '-eps-sigtex_smooth.csv')
        else:
            out_arr = np.vstack([e.eps_asc, scale_factor * e.sig_tex_asc]).T
            out_file = os.path.join(
                dir_name, base_name + '-eps-sigtex_raw.csv')

        print 'writing output to', out_file

#        header_string = header_template % (base_name, e.ccs.thickness, e.width,
#                                           e.gauge_length, e.A_tex,
#                                           e.rho_c, fabric_type,
#                                           n_fabric_layers)

        header_string = header_template % (base_name, thickness, width,
                                           e.gauge_length, A_tex,
                                           rho_c, fabric_type,
                                           n_fabric_layers)

        # get textile stress based on measured roving number
        #
        sigtex_max_list = sigtex_max_list + [e_run.ex_type.F_asc[-1] * 1000. / A_tex ]

        np.savetxt(out_file,
                   out_arr, delimiter=';')  # ,
#                   header=header_string)

        # workaround for old numpy.savetxt version:
        f = open(out_file, 'r')
        temp = f.read()
        f.close()

        f = open(out_file, 'w')
        f.write(header_string)

        f.write(temp)
        f.close()

        # plot
        e._plot_sigtex_eps(axes, color=color_list[idx], linestyle=linestyle_list[idx], label=test_files[idx], plot_analytical_stiffness_II=False)
        axes.grid()
        axes.set_xlabel('eps [-]')
        axes.set_ylabel('sig_tex [MPa]')
        axes.legend(loc=2)
        axes.set_xlim([0., 0.012])
        axes.set_ylim([0., 2500])

    print 'sigtex_max_list: ', sigtex_max_list
    print 'sigtex (average) = %g' % (np.average(np.array(sigtex_max_list)))
    print 'sigtex (standard deviation) = %g' % (np.std(np.array(sigtex_max_list)))

    p.savefig('eps-sigtex.pdf', format='pdf')
    p.show()

sig_tex_max_arr_db = np.array([e_run.ex_type.sig_tex_max
                            for e_run in e_list], dtype='float_')
print 'sigtex_db (average) = %g' % (np.average(sig_tex_max_arr_db))
print 'sigtex_db (standard deviation) = %g' % (np.std(sig_tex_max_arr_db))

if __name__ == '__main__':
#    plot_eps_sigtex()
    eval_sig_tex_eps(smooth=False)
