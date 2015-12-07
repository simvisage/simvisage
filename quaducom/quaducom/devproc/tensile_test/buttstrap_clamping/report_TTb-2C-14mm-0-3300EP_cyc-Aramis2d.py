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


#--------------------
# TTb-2C-3cm-3300EP
#--------------------
# test_files = ['TTb-2C-3cm-0-3300EP-V4_Aramis2d.DAT',
#             'TTb-2C-3cm-0-3300EP-V5_cyc-Aramis2d.DAT']
#
# test_file_path = os.path.join(simdb.exdata_dir,
#                              'tensile_tests', 'buttstrap_clamping',
#                              '2015-08-18_TTb-2C-3cm-0-3300EP_cyc-Aramis-2d')
#
# n_rov_list = [9, 9]
# color_list = ['g', 'b']
# linestyle_list = ['-', '-']
# plot_asc_list = [1, 1]
# xlim = 0.02
# ylim = 2000.

#--------------------
# TTb-2C-14mm-3300EP
#--------------------
# test_files = ['TTb-2C-14mm-0-3300EP-V2_cyc-Aramis2d.DAT',
#               'TTb-2C-14mm-0-3300EP-V1_Aramis2d.DAT']
#
# test_file_path = os.path.join(simdb.exdata_dir,
#                              'tensile_tests', 'buttstrap_clamping',
#                              '2015-08-11_TTb-2C-14mm-0-3300EP_cyc-Aramis2d')
#
# n_rov_list = [9, 9]
# color_list = ['k', 'grey']
# linestyle_list = ['-', '-']
# plot_asc_list = [1, 1]
# label_cutoff = [-9, -9]  # cutoff long label names at the end for cleaner legend display
# xlim = 0.025
# ylim = 3500.

#--------------------
# TTb-1C-3cm-3300EP
#--------------------
test_files = [  # 'TTb-1C-3cm-0-3300EP-V6_cyc-Aramis2d-sideview-notched.DAT',
              'TTb-1C-3cm-0-3300EP-V4_cyc-Aramis2d-sideview.DAT',
#                'TTb-1C-3cm-0-3300EP-V1_Aramis2d.DAT',
               'TTb-1C-3cm-0-3300EP-V2_Aramis2d.DAT']

test_file_path = os.path.join(simdb.exdata_dir,
                             'tensile_tests', 'buttstrap_clamping',
                             '2015-08-10_TTb-1C-3cm-0-3300EP_cyc-Aramis2d')
print 'test_file_path', test_file_path

n_rov_list = [5, 5]
color_list = ['k', 'grey']
linestyle_list = ['-', '-']
plot_asc_list = [1, 1]
label_cutoff = [-18, -9]  # cutoff long label names at the end for cleaner legend display
xlim = 20
ylim = 3500.
k_rho_list = [1.187, 1.187]  # = 9,1 mm^2 / 10,8 mm^2 = (5 rovings * 1,84 mm^2) / (90 mm^2/m * 0.12 m)

#--------------------
# TTb-2C-14mm-3300SBR
#--------------------
# test_files = [ 'TTb-2C-14mm-0-3300SBR-V1_cyc-Aramis2d.DAT',
#                'TTb-2C-14mm-0-3300SBR-V3_Aramis2d.DAT']
#
# test_file_path = os.path.join(simdb.exdata_dir,
#                              'tensile_tests', 'buttstrap_clamping',
#                              '2015-08-03_TTb-2C-14mm-0-3300SBR_cyc-Aramis2d')
#
# n_rov_list = [14, 18]
# color_list = ['k', 'grey']
# linestyle_list = ['-', '-']
# plot_asc_list = [0, 1]
# label_cutoff = [-9, -9]  # cutoff long label names at the end for cleaner legend display
# xlim = 10
# ylim = 1800.
# k_rho_list = [1.002, 1.08]  # = 9,1 mm^2 / 10,8 mm^2 = (5 rovings * 1,84 mm^2) / (90 mm^2/m * 0.12 m)

#--------------------

# test_files = [ 'TTb-2C-14mm-0-3300SBR-V1_cyc-Aramis2d.DAT']
#
# test_file_path = os.path.join(simdb.exdata_dir,
#                              'tensile_tests', 'buttstrap_clamping',
#                              '2015-08-03_TTb-2C-14mm-0-3300SBR_cyc-Aramis2d')
#
# color_list = ['k']
# linestyle_list = ['-']
# label_cutoff = [-9]  # cutoff long label names at the end for cleaner legend display
# xlim = 20.
# ylim = 30.


#--------------------
# TTb-2C-14mm-800SBR
#--------------------
# test_files = [ 'TTb-2C-1cm-0-800SBR-V4_cyc-Aramis2d.DAT',
#                'TTb-2C-1cm-0-800SBR-V1_Aramis2d.DAT']
#
# test_file_path = os.path.join(simdb.exdata_dir,
#                              'tensile_tests', 'buttstrap_clamping',
#                              '2015-04-20_TTb-2C-1cm-0-800SBR_cyc-Aramis2d')
#
# n_rov_list = [9, 9]
# color_list = ['k', 'grey']
# linestyle_list = ['-', '-']
# plot_asc_list = [0, 1]
# xlim = 0.015
# ylim = 2000.
#--------------------




e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

def plot_all():

    fig = p.figure(facecolor='white', figsize=(8, 6))
#     fig.subplots_adjust(
#         left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_list):

        e = e_run.ex_type

        axes = p.subplot(111)

        if plot_asc_list[idx]:
            e._plot_tex_stress_strain_asc(axes, xscale=1000., k_rho=k_rho_list[idx], color=color_list[idx], linewidth=1.5, plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=e_list[idx].ex_type.key[0:label_cutoff[idx]])
        else:
            e._plot_tex_stress_strain(axes, xscale=1000., k_rho=k_rho_list[idx], color=color_list[idx], linewidth=1.5, plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=e_list[idx].ex_type.key[0:label_cutoff[idx]])
        axes.set_ylabel('Textilspannung [MPa]')
        axes.set_xlabel('Dehnung [1E+3]')

#         e._plot_force_displacement_machine(axes, color=color_list[idx], linewidth=1.5, label=e_list[idx].ex_type.key[0:label_cutoff[idx]])
#         axes.set_xlabel('Maschinenweg [mm]')
#         axes.set_ylabel('Kraft [kN]')

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
    test_series_name = 'TTb_cyc-Aramis2d'
    if save_fig_to_file:
        print 'XXX'
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
        filename = os.path.join(test_series_dir, 'sigtex-epsu.png')
        p.savefig(filename, format='png')
        print 'figure saved to file %s' % (filename)

if __name__ == '__main__':
    plot_all()
    p.show()


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




