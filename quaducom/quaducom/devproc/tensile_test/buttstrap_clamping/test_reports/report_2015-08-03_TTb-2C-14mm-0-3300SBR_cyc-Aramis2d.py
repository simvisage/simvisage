'''
Created on 14.12.2015

@author: alexander
'''
from matresdev.db.simdb.simdb import simdb
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

#--------------------
# TTb-2C-14mm-3300SBR
#--------------------
# do = 'sigtex-eps'  # gauge displacement
do = 'F-w'  # machine displacement

test_file_path = os.path.join(simdb.exdata_dir,
                             'tensile_tests', 'buttstrap_clamping',
                             '2015-08-03_TTb-2C-14mm-0-3300SBR_cyc-Aramis2d')

test_files = [ 'TTb-2C-14mm-0-3300SBR-V1_cyc-Aramis2d.DAT']
#               'TTb-2C-14mm-0-3300SBR-V3_Aramis2d.DAT']
#                 'TTb-2C-14mm-0-3300SBR-V2_Aramis2d.DAT']
#                 'TTb-2C-14mm-0-3300SBR-V5_Aramis2d.DAT']


A_rov = 1.84  # 3300 tex
n_rov_list = [14, 18]
color_list = ['k', 'grey']
linestyle_list = ['-', '-']
plot_asc_list = [0, 1]
label_cutoff = [-9, -9]  # cutoff long label names at the end for cleaner legend display

k_rho_list = [1.080,  # = 27.82 mm^2 / 25.76 mm^2 = (2layers * 144.9 mm^2/m * 0.096 m) / (14 rovings * 1,84 mm^2)
              0.998]  # = 33.04 mm^2 / 33.12 mm^2 = (2layers * 144.9 mm^2/m * 0.114 m) / (18 rovings * 1,84 mm^2)

e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

# compare 'A_tex' as calculated in simdb based on 'a_tex_0 [mm^2/m]' and the specimen 'width [m]'
# and the real number of rovings as counted in the test:
#for n, n_rov in enumerate(n_rov_list):
#    k_rho_list[n] = e_list[n].ex_type.A_tex / (n_rov * A_rov)
#print '\n'
#print 'k_rho_list (caclulated based on simdb data and n_rovings specified in script: \n', k_rho_list
#print '\n'

def plot_all():

    fig = p.figure(facecolor='white', figsize=(8, 6))

    for idx, e_run in enumerate(e_list):

        e = e_run.ex_type

        axes = p.subplot(111)

        if do == 'sigtex-eps':
            if plot_asc_list[idx]:
                e._plot_tex_stress_strain_asc(axes, xscale=1000., k_rho=k_rho_list[idx], color=color_list[idx], linewidth=1.5,
                                              plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=e_list[idx].ex_type.key[0:label_cutoff[idx]])
            else:
                e._plot_tex_stress_strain(axes, xscale=1000., k_rho=k_rho_list[idx], color=color_list[idx], linewidth=1.5, plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=e_list[idx].ex_type.key[0:label_cutoff[idx]])
            axes.set_ylabel('Textilspannung [MPa]')
            axes.set_xlabel('Dehnung [1E+3]')
            xlim = 10
            ylim = 1800.

        if do == 'F-w':
            e._plot_force_displacement_machine(axes, color=color_list[idx], linewidth=1.5, label=e_list[idx].ex_type.key[0:label_cutoff[idx]])
            axes.set_xlabel('Maschinenweg [mm]')
            axes.set_ylabel('Kraft [kN]')
            xlim = 20
            ylim = 30.

        axes.axis([0., xlim, 0., ylim])

#     # material stiffness carbon (E=245GPa)
#     xarr = np.array([0., 0.010])
#     yarr = np.array([0., 2450.])
#     axes.plot(xarr, yarr, linestyle='--', color='grey', linewidth=1.5, label='E_tex = 245 GPa')

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
        filename = os.path.join(test_series_dir, do + '.eps')
        p.savefig(filename, dpi=300)
        print 'figure saved to file %s' % (filename)

if __name__ == '__main__':
    plot_all()
    p.show()
