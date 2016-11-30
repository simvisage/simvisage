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

data_name = 'R4 Aachen'

#--------------------
# TTb-2C-9mm-0-3300SBR-DK3_A5_R4
#--------------------
do = 'sigtex-eps'  # stress / strain
#do = 'F-w'  # Force / gauge displacement 

test_file_path = os.path.join(simdb.exdata_dir,
                             'tensile_tests', 'buttstrap_clamping',
                             '2016-03-14_TTb-2C-9mm-0-3300SBR_R4')

test_files = ['TTb-2C-9mm-0-3300SBR-DK3_A5_R4.DAT', 
              'TTb-2C-9mm-0-3300SBR-DK3_A6_R4.DAT', 
              'TTb-2C-9mm-0-3300SBR-DK3_B1_R4.DAT', 
              'TTb-2C-9mm-0-3300SBR-DK3_B2_R4.DAT',
              'TTb-2C-9mm-0-3300SBR-DK3_C3_R4.DAT',
              'TTb-2C-9mm-0-3300SBR-DK3_C4_R4.DAT',]

#n_rov_list = [10, 10]
color_list = ['k', 'grey','red','blue',"green","teal"]
linestyle_list = ['-', '-', '-', '-', '-', '-']
plot_asc_list = [1, 1, 1, 1, 1, 1]
label_cutoff = [-3, -3, -3, -3, -3, -3]  # cutoff long label names at the end for cleaner legend display
k_rho_list = [1, 1, 1, 1, 1, 1]  # modification factor to ajust to actual textile area in the given width of the specimen
n_roving_array = [18, 18, 18, 18, 18, 18]
e_array = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

def plot_all():

    fig = p.figure(facecolor='white', figsize=(8, 6))
#     fig.subplots_adjust(
#         left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_array):

        e = e_run.ex_type

        axes = p.subplot(111)

        if do == 'sigtex-eps':
            if plot_asc_list[idx]:
                e._plot_tex_stress_strain_asc(axes, xscale=1000., k_rho=k_rho_list[idx], color=color_list[idx], linewidth=1.5,
                                              plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=e_array[idx].ex_type.key[0:label_cutoff[idx]])
            else:
                e._plot_tex_stress_strain(axes, xscale=1000., k_rho=k_rho_list[idx], color=color_list[idx], linewidth=1.5, plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=e_array[idx].ex_type.key[0:label_cutoff[idx]])
            axes.set_ylabel('Textile Stress [MPa]')
            axes.set_xlabel('Strain [1E+3]')
            xlim = 15
            ylim = 1500.

        if do == 'F-w':
            e._plot_force_displacement_asc(axes, color=color_list[idx], linewidth=1.5, label=e_array[idx].ex_type.key[0:label_cutoff[idx]])
            axes.set_xlabel('Weg [mm]')
            axes.set_ylabel('Kraft [kN]')
            xlim = 5
            ylim = 50.

        axes.axis([0., xlim, 0., ylim])

    axes.grid()
    axes.legend(loc=4)

    # --------------------------------
    # save figure
    # --------------------------------
    save_fig_to_file = False

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
        p.savefig(filename)
        print 'figure saved to file %s' % (filename)

if __name__ == '__main__':
    plot_all()
    p.show()

# displacement gauge dropped off --> precess data
#--------------------
# TTb-2C-3cm-3300EP
#--------------------
# idx_1 = np.where(sig_tex_asc > 1091.)[0][0]
# idx_2 = np.where(sig_tex_asc > 977.)[0][0]
# idx_cut = np.where(sig_tex_asc > 1091.)[0]
# K_cracked = (sig_tex_asc[idx_2] - sig_tex_asc[idx_1]) / (eps_asc[idx_2] - eps_asc[idx_1])
# sig_tex_asc[idx_cut] = K_cracked * eps_asc[idx_cut]
# idx_cut = np.where(sig_tex_asc < 1761.)[0]
# sig_tex_asc = sig_tex_asc[idx_cut]
# eps_asc = eps_asc[idx_cut]
# print 'idx_1', idx_1
# print 'idx_2', idx_2
# print 'K_cracked', K_cracked
# print 'sig_tex_asc[idx_1]', sig_tex_asc[idx_1]

