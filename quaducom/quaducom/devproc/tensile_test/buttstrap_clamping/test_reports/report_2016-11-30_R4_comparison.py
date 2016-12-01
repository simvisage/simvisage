'''
Created on Nov 30, 2016

'''

import os
import scipy
from traits.api import HasTraits

from matresdev.db.exdb import ExRun
from matresdev.db.simdb.simdb import simdb
import numpy as np
import pylab as p
#import report_2016_03_16_DPO_3300SBR_R4_Dresden as dd
#import report_2016_03_16_DPO_3300SBR_R4_Dresden_PO as lz
import report_2016_03_15_TTb_2C_9mm_0_3300SBR_DK3_R4_ac as ac
print scipy.__version__


class PlotBase(HasTraits):

    def figure(self):
        fig = p.figure(facecolor='white', figsize=(12, 9))
        fig.subplots_adjust(
            left=0.09, right=0.97, bottom=0.09, top=0.96, wspace=0.25, hspace=0.2)
        axes = p.subplot(111)
        return axes

class PlotSE(PlotBase):

    def plot(self, dataset, axes, linestyle=None, color=None, label=None):

       # gauge_dist = dataset.gauge_dist
       # e_list = np.array(dataset.e_array).reshape(3, -1)

        for idx, (e_run, n_r) in enumerate(zip(dataset.e_array, dataset.n_roving_array)):
            if idx == 1:
                lb = label
            else:
                lb = None
            e = e_run.ex_type
            e._plot_tex_stress_strain_asc(axes, xscale=1000., k_rho=dataset.k_rho_list[idx], color=color, linewidth=1.5,
                                           plot_analytical_stiffness_II=False, plot_analytical_stiffness_I=False, label=lb)
        
        
    def decorate(self, axes):
        axes.grid()
        axes.set_xlabel('Textile Stress [MPa]')
        axes.set_ylabel('Strain [1E+3]')
        axes.legend(loc=2)
        axes.axis([0., 15, 0., 1500])



if __name__ == '__main__':

    pw = PlotSE()
#    pw = PlotRotation()
#    pw = PlotFW()
    ax = pw.figure()
    #pw.plot(dd, ax, color='red', markerstyle='v', label='DPO Dresden')
    #pw.plot(ac, ax, linestyle='dashed', color='black', markerstyle='o', linewidth=2, label='DPO Aachen')
    pw.plot(ac, ax, linestyle='dashed', color='black', label= 'R4 Aachen')
    pw.decorate(ax)

    p.show()
