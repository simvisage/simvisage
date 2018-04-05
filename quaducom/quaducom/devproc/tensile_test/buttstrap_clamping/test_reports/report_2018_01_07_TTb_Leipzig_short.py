'''
Created on 04.04.2018

@author: janb
'''
#TO DO: 

from matresdev.db.simdb.simdb import simdb
from matresdev.db.exdb import ExRun

import os
from numpy import \
    array
import numpy as np
import pylab as p

# specify font options for plots
params = {'legend.fontsize': 12,
#         'legend.linewidth': 2,
          u'font.size':12,
          u'font.family':'Times New Roman',
          u'font.style':'normal'}
p.rcParams.update(params)
# print p.rcParams.keys()

data_name = 'Leipzig_short'

test_file_path = os.path.join(simdb.exdata_dir,
                             'tensile_tests', 'buttstrap_clamping',
                             '2018-01-07_TTb_Leipzig', 'raw')

test_files = ['DKBC11.DAT', 
              'DKBC12.DAT', 
              'DKBC13.DAT', 
              'DKBC14.DAT',
              'DKBC15.DAT',
              'DKBC16.DAT',]

#n_rov_list = [10, 10]
color_list = ['blue', 'deepskyblue','green','palegreen','red','orange']
linestyle_list = ['-', '-', '-', '-', '-','-',]
plot_asc_list = [1, 1, 1, 1, 1, 1,]

#a_roving_list = [3.62, 3.62, 3.62, 3.62, 3.62]
label_cutoff = [1, 1, 1, 1, 1, 1]  # cutoff long label names at the end for cleaner legend display
e_array = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]
tension_max = np.zeros((len(test_files),))
def plot_all():

    fig = p.figure(facecolor='white', figsize=(8.8/2.54, 8.8/2.54), dpi=100)
    fig.subplots_adjust(left=0.19, right=0.96, bottom=0.15, top=0.93, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_array):

        e = e_run.ex_type

        axes = p.subplot(111)
        a_roving_0 = e.ccs.fabric_layup_list[1].a_roving_0
               
        if plot_asc_list[idx]:
            e._plot_tex_stress_strain_asc(axes, xscale=1000.,  k_rho=e.A_tex/(e.n_rovings*a_roving_0), color=color_list[idx], linewidth=1.5,
                                              plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=None)
        else:
            e._plot_tex_stress_strain(axes, xscale=1000.,  k_rho=e.A_tex/(e.n_rovings*a_roving_0), color=color_list[idx], linewidth=1.5, plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=e_array[idx].ex_type.key[0:label_cutoff[idx]])
        axes.set_ylabel('Textile Stress [MPa]')
        axes.set_xlabel('Strain [1E+3]')
        xlim = 12
        ylim = 2500.      
        tension_max[idx] = 1000.*max(e.Kraft)/(e.n_rovings*a_roving_0)
        axes.axis([0., xlim, 0., ylim],)

# Print Analytical Roving stiffness. Ultimate tensile stress (average value) = 3221 N/mm2
#     E_tex = e.ccs.E_tex
#     K_IIb = E_tex
#     eps_ru = [0, 3221/K_IIb*1000]
#     sig_ru = [0, 3221]
#     axes.plot(eps_ru, sig_ru, color='k', linestyle='--', label='Faserstrang-Zug')    
# 
#     axes.grid(b=True, which='major', color='gray', linestyle='-', linewidth = .5,)
#     axes.legend(loc=4, handletextpad=0.09)
#     print 'maximum tension =', tension_max
#     print 'average maximum tension =', np.average(tension_max)
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
        filename = os.path.join(test_series_dir, 'fig1.pdf')
        p.savefig(filename)
        print 'figure saved to file %s' % (filename)

if __name__ == '__main__':
    plot_all()
    p.show()

