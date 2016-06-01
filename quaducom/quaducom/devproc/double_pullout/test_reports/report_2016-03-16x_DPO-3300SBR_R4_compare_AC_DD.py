'''
Created on Mar 17, 2016

'''

import os
import scipy
from traits.api import HasTraits

from matresdev.db.exdb import ExRun
from matresdev.db.simdb.simdb import simdb
import numpy as np
import pylab as p
import report_2016_03_16_DPO_3300SBR_R4_Dresden as dd
import report_2016_03_16_DPO_3300SBR_R4_Dresden_PO as dd_po
import report_2016_03_16_DPO_3300SBR_R4_rawdata as ac
print scipy.__version__


class PlotBase(HasTraits):

    def figure(self):
        fig = p.figure(facecolor='white', figsize=(12, 9))
        fig.subplots_adjust(
            left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)
        axes = p.subplot(111)
        return axes


class PlotRotation(PlotBase):

    def plot(self, dataset, axes, markerstyle=None, linestyle=None, color=None, linewidth=1.5):
        gauge_dist = dataset.gauge_dist
        for idx, (e_treatment, e_treatment_n_r) in enumerate(zip(dataset.e_array,
                                                                 dataset.n_roving_array)):
            for e_run, n_r in zip(e_treatment, e_treatment_n_r):
                e = e_run.ex_type
                max_F_idx = np.argmax(e.Kraft)
                if linestyle == None:
                    linestyle = dataset.linestyle_list[idx]
                if color == None:
                    color = dataset.color_list[idx]
                w_left = getattr(e, dataset.left_gauge_name)
                w_right = getattr(e, dataset.right_gauge_name)
                axes.plot(e.Kraft[:max_F_idx] / n_r,
                          (np.fabs(w_right - w_left)[:max_F_idx]) / gauge_dist,
                          linewidth=linewidth,
                          linestyle=linestyle, color=color)

    def decorate(self, axes):
        axes.grid()
        axes.set_ylabel('$\Delta$ w [mm] / L [mm]')
        axes.set_xlabel('force per roving [kN]')
        axes.legend(loc=2)
        axes.axis([0., 2.5, 0., 0.025])


class PlotFW(PlotBase):

    def plot(self, dataset, axes, markerstyle=None, linestyle=None, color=None, linewidth=1.5):

        gauge_dist = dataset.gauge_dist
        for idx, (e_treatment, e_treatment_n_r) in enumerate(zip(dataset.e_array,
                                                                 dataset.n_roving_array)):
            for e_run, n_r in zip(e_treatment, e_treatment_n_r):
                if linestyle == None:
                    linestyle = dataset.linestyle_list[idx]
                if color == None:
                    color = dataset.color_list[idx]
                e = e_run.ex_type
                max_F_idx = np.argmax(e.Kraft)
                axes.plot(e.w[:max_F_idx],
                          e.Kraft[:max_F_idx] / n_r,
                          linewidth=linewidth,
                          linestyle=linestyle, color=color)

    def decorate(self, axes):
        axes.grid()
        axes.set_xlabel('$\Delta$ w [mm]')
        axes.set_ylabel('force per roving [kN]')
        axes.legend(loc=2)
        axes.axis([0., 8, 0., 2.5])


class PlotFWAvg(PlotBase):

    def plot(self, dataset, axes, markerstyle=None, linestyle=None, color=None):

        gauge_dist = dataset.gauge_dist
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
        axes.set_xlabel('$\Delta$ w [mm]')
        axes.set_ylabel('force per roving [kN]')
        axes.legend(loc=2)
        axes.axis([0., 8, 0., 2.5])

class PlotFL(PlotBase):

    def plot(self, dataset, axes, markerstyle=None, linestyle=None, color=None, linewidth=1.5):

#        gauge_dist = dataset.gauge_dist
        for idx, (e_treatment, e_treatment_n_r, e_treatment_l_v) in enumerate(zip(dataset.e_array,
                                                                 dataset.n_roving_array, dataset.l_v_array)):
            for e_run, n_r, l_v in zip(e_treatment, e_treatment_n_r, e_treatment_l_v):
                if linestyle == None:
                    linestyle = dataset.linestyle_list[idx]
                #if color == None:
                    #color = dataset.color_list[idx]
                e = e_run.ex_type
#                max_F_idx = np.argmax(e.Kraft)
                axes.plot(l_v, e.Kraft.max() / n_r,
                 marker=markerstyle, markersize=8, color = color)
#               axes.plot(e.w[:max_F_idx],
#                        e.Kraft[:max_F_idx] / n_r,
#                         linewidth=linewidth,
#                         linestyle=linestyle, color=color)

    def decorate(self, axes):
        axes.grid()
        axes.set_ylabel('Length [mm]')
        axes.set_xlabel('force per roving [kN]')
        axes.legend(loc=2)
        axes.axis([0., 400, 0., 2.5])          

if __name__ == '__main__':

    pw = PlotFL()
#    pw = PlotRotation()
#    pw = PlotFW()
    ax = pw.figure()
    pw.plot(dd, ax, color='black', markerstyle='v')
    pw.plot(dd_po, ax, color='blue', markerstyle='D')
    pw.plot(ac, ax, linestyle='dashed', color='red', markerstyle='o', linewidth=2)
    pw.decorate(ax)

    p.show()
