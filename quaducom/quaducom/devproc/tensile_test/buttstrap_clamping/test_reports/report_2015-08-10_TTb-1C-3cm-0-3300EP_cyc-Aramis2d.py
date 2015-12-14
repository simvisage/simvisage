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
# TTb-1C-3cm-3300EP
#--------------------
do = 'sigtex-eps'  # gauge displacement
# do = 'F-w' # machine displacement

test_file_path = os.path.join(simdb.exdata_dir,
                             'tensile_tests', 'buttstrap_clamping',
                             '2015-08-10_TTb-1C-3cm-0-3300EP_cyc-Aramis2d')

test_files = [  # 'TTb-1C-3cm-0-3300EP-V6_cyc-Aramis2d-sideview-notched.DAT',
              'TTb-1C-3cm-0-3300EP-V4_cyc-Aramis2d-sideview.DAT',
#                'TTb-1C-3cm-0-3300EP-V1_Aramis2d.DAT',
               'TTb-1C-3cm-0-3300EP-V2_Aramis2d.DAT']

n_rov_list = [5, 5]
color_list = ['k', 'grey']
linestyle_list = ['-', '-']
plot_asc_list = [1, 1]
label_cutoff = [-18, -9]  # cutoff long label names at the end for cleaner legend display
k_rho_list = [1.187, 1.187]  # = 9,1 mm^2 / 10,8 mm^2 = (5 rovings * 1,84 mm^2) / (90 mm^2/m * 0.12 m)


e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

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
            xlim = 20
            ylim = 3500.

        if do == 'F-w':
            e._plot_force_displacement_machine(axes, color=color_list[idx], linewidth=1.5, label=e_list[idx].ex_type.key[0:label_cutoff[idx]])
            axes.set_xlabel('Maschinenweg [mm]')
            axes.set_ylabel('Kraft [kN]')
            xlim = 12
            ylim = 35.

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
        img_dir = os.path.join(simdb.exdata_dir, 'img_dir')
        test_series_dir = os.path.join(simdb.report_dir, test_series_name)
        if not os.path.exists(test_series_dir):
            os.makedirs(test_series_dir)
        filename = os.path.join(test_series_dir, 'sigtex-epsu.png')
        p.savefig(filename, format='png')
        print 'figure saved to file %s' % (filename)

if __name__ == '__main__':
    plot_all()
    p.show()
