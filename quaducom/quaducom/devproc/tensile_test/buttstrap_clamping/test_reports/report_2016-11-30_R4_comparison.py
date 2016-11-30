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
import report_2016_03_15_TTb_2C_9mm_0_3300SBR_DK3_R4_all as ac
print scipy.__version__


class PlotBase(HasTraits):

    def figure(self):
        fig = p.figure(facecolor='white', figsize=(12, 9))
        fig.subplots_adjust(
            left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)
        axes = p.subplot(111)
        return axes

class Plot_sig_eps(PlotBase):

    def plot(self, dataset, axes, markerstyle=None, linestyle=None, color=None):

       # gauge_dist = dataset.gauge_dist
       # e_list = np.array(dataset.e_array).reshape(3, -1)

        for idx, (e_run, n_r) in enumerate(zip(dataset.e_array, dataset.n_roving_array)):
            if linestyle == None:
                linestyle = dataset.linestyle_list[idx]
            if color == None:
                color = dataset.color_list[idx]
            e = e_run.ex_type
            max_F_idx = np.argmax(e.Kraft)
            axes.plot(e.w[:max_F_idx],
                      e.Kraft[:max_F_idx] / n_r,
                      linewidth=1.5,
                      linestyle=linestyle, color=color)

    def decorate(self, axes):
        axes.grid()
        axes.set_xlabel('Textile Stress [MPa]')
        axes.set_ylabel('Strain [1E+3]')
        axes.legend(loc=2)
        axes.axis([0., 15, 0., 1500])



if __name__ == '__main__':

    pw = Plot_sig_eps()
#    pw = PlotRotation()
#    pw = PlotFW()
    ax = pw.figure()
    #pw.plot(dd, ax, color='red', markerstyle='v', label='DPO Dresden')
    #pw.plot(dd_po, ax, color='blue', markerstyle='D', label='PO Dresden')
    #pw.plot(ac, ax, linestyle='dashed', color='black', markerstyle='o', linewidth=2, label='DPO Aachen')
    pw.plot(ac, ax, linestyle='dashed', color='black', markerstyle='o', linewidth=2)
    pw.decorate(ax)

    p.show()
