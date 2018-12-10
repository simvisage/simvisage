'''
Created on 19.06.2017

@author: janb
'''


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


test_file_path = os.path.join(simdb.exdata_dir,
                             'shear_tests', 'four_point', '2017-05-23_SH-4PT-1G-3cm')

test_files = ['SH475V1.DAT', 
              'SH475V2.DAT', 
              'SH475V3.DAT',]

color_list = ['red','blue','green',]
linestyle_list = ['-', '-', '-',]
plot_asc_list = [0, 0, 0,]
e_array = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]
def plot_all():

    fig = p.figure(facecolor='white', figsize=(6, 6))
#     fig.subplots_adjust(
#         left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_array):

        e = e_run.ex_type

        axes = p.subplot(111)
        if plot_asc_list[idx]:
            e._plot_shearforce_deflection_asc(axes, color=color_list[idx], linewidth=1.5,
                                              label= test_files[idx].split('.')[0],)
        else:
            e._plot_shearforce_deflection(axes, color=color_list[idx], linewidth=1.5, label= test_files[idx].split('.')[0] )
        axes.set_ylabel('Querkraft [KN]')
        axes.set_xlabel('Mitteldurchbiegung [mm]')
        xlim = 20
        ylim = 10.      
        axes.axis([0., xlim, 0., ylim])

    axes.grid(b=True, which='major', color='gray', linestyle='-', linewidth = .5,)
    axes.legend(loc=1)
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
        filename = os.path.join(test_series_dir, '75.pdf')
        p.savefig(filename)
        print('figure saved to file %s' % (filename))

if __name__ == '__main__':
    plot_all()
    p.show()

