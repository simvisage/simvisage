'''
Created on Nov 30, 2016
by janb
'''

import os
import scipy
from traits.api import HasTraits

from matresdev.db.exdb import ExRun
from matresdev.db.simdb.simdb import simdb
import numpy as np
import pylab as p
import report_2016_03_18_TTb_2C_9mm_0_3300SBR_DK3_R4_dd as dd
import report_2016_03_18_TTb_2C_9mm_0_3300SBR_DK3_R4_MFPA as mfpa
import report_2016_03_15_TTb_2C_9mm_0_3300SBR_DK3_R4_ac as ac
print scipy.__version__

# specify font options for plots
params = {'legend.fontsize': 30,
          'ps.fonttype': 42,
          u'font.size':30,
          u'font.family':'serif',
          u'font.style':'normal',
          u'font.serif': 'Times New Roman'}

p.rcParams.update(params)


class PlotBase(HasTraits):

    def figure(self):
        fig = p.figure(facecolor='white', figsize=(12, 9))
        fig.subplots_adjust(
            left=0.12, right=0.97, bottom=0.09, top=0.96, wspace=0.25, hspace=0.2)
        axes = p.subplot(111)
        return axes

class PlotSE(PlotBase):

    def plot(self, dataset, axes, linestyle=None, color=None, label=None):

     for idx, (e_run) in enumerate(dataset.e_array):
            if idx == 1:
                lb = label
            else:
                lb = None
            e = e_run.ex_type
            a_roving_0 = e.ccs.fabric_layup_list[1].a_roving_0
            e._plot_tex_stress_strain_asc(axes, xscale=1000., k_rho=e.A_tex/(e.n_rovings*a_roving_0), linestyle=linestyle, color=color,
                                           linewidth=1.5, plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=lb)
        
        
    def decorate(self, axes):
        axes.grid()
        axes.set_xlabel('Strain [1E+3]')
        axes.set_ylabel('Textile Stress [MPa]')
        axes.legend(loc=2)
        axes.axis([0., 15, 0., 1500])

if __name__ == '__main__':

    pw = PlotSE()

    ax = pw.figure()
    pw.plot(dd, ax, linestyle='dashed', color='red',  label='R4 Dresden')
    pw.plot(mfpa, ax, linestyle='dashdot', color='blue',  label='R4 MFPA Leipzig')
    pw.plot(ac, ax, linestyle='-', color='black', label= 'R4 Aachen')
    pw.decorate(ax)
    
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
        filename = os.path.join(test_series_dir, '1'+ '.eps')
        p.savefig(filename)
        print 'figure saved to file %s' % (filename)

    p.show()



