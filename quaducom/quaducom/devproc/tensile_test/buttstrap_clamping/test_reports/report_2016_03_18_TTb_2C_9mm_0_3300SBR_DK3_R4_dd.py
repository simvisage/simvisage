'''
Created on 01.12.2016

@author: janb
'''
#TO DO: Specify the cross sectional area of one roving within the fabric_layup, delete the 1.84 mm^2 used in the k_roh - calculation

from matresdev.db.simdb.simdb import simdb
from matresdev.db.exdb import ExRun

import os
import numpy as np
import pylab as p

# specify font options for plots
params = {'legend.fontsize': 12,
#         'legend.linewidth': 2,
          'font.size':15,
          'font.family':'serif',
          'font.style':'normal'}
p.rcParams.update(params)
# print p.rcParams.keys()

data_name = 'R4 Dresden'

test_files = ['B3-Ia-DK3-B-3.DAT', 
              'B3-Ia-DK3-B-4.DAT', 
              'B3-Ia-DK3-B-7.DAT', 
              'B3-Ia-DK3-C-5.DAT',
              'B3-Ia-DK3-C-6.DAT',
              ]

test_file_path = os.path.join(simdb.exdata_dir,
                             'tensile_tests', 'buttstrap_clamping',
                             '2016-03-18_TTb-2C-9mm-0-3300SBR_R4_Dresden')


#n_rov_list = [10, 10]
color_list = ['k', 'grey','red','blue',"green",]
linestyle_list = ['-', '-', '-', '-', '-',]
plot_asc_list = [1, 1, 1, 1, 1,]
label_cutoff = [-6, -3, -3, -3, -3,]  # cutoff long label names at the end for cleaner legend display
e_array = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]
tension_max = np.zeros((len(test_files),))
def plot_all():

    fig = p.figure(facecolor='white', figsize=(8, 6))
#     fig.subplots_adjust(
#         left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_array):

        e = e_run.ex_type

        axes = p.subplot(111)
        a_roving_0 = e.ccs.fabric_layup_list[1].a_roving_0      
        if plot_asc_list[idx]:
            e._plot_tex_stress_strain_asc(axes, xscale=1000., k_rho=e.A_tex/(e.n_rovings*a_roving_0), color=color_list[idx], linewidth=1.5,
                                              plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=e_array[idx].ex_type.key[0:])
#
        else:
            e._plot_tex_stress_strain(axes, xscale=1000., k_rho=e.A_tex/(e.n_rovings*a_roving_0), color=color_list[idx], linewidth=1.5, plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=e_array[idx].ex_type.key[0:label_cutoff[idx]])
        axes.set_ylabel('Textile Stress [MPa]')
        axes.set_xlabel('Strain [1E+3]')
        xlim = 15
        ylim = 1500.
        tension_max[idx] = 1000.*max(e.Kraft)/(e.n_rovings*a_roving_0)
        axes.axis([0., xlim, 0., ylim])

    axes.grid()
    axes.legend(loc=4)
    print('maximum tension =', tension_max)
    print('average maximum tension =', np.average(tension_max))
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
        filename = os.path.join(test_series_dir, '.eps')
        p.savefig(filename)
        print('figure saved to file %s' % (filename))

if __name__ == '__main__':
    plot_all()
    p.show()

