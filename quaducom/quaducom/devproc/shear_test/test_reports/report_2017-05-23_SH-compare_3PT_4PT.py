'''
Created on Jun 21, 2017

'''

import os
import scipy
from traits.api import HasTraits

from matresdev.db.exdb import ExRun
from matresdev.db.simdb.simdb import simdb
import numpy as np
import pylab as p
from . import report_2017_05_04_SH_3PT_comparison_a_d as threePT
from . import report_2017_05_23_SH_4PT_comparison_a_d as fourPT

class PlotBase(HasTraits):

    def figure(self):
        fig = p.figure(facecolor='white', figsize=(12, 9))
        fig.subplots_adjust(
            left=0.07, right=0.97, bottom=0.08, top=0.96, wspace=0.25, hspace=0.2)
        axes = p.subplot(111)
        return axes

class PlotQW(PlotBase):

    def plot(self, dataset, axes, markerstyle=None, linestyle=None, color=None, label=None):
    
#        gauge_dist = dataset.gauge_dist
        for idx, (e_treatment, e_treatment_ad, e_treatment_marker, e_treatment_color) in enumerate(zip(dataset.e_array, dataset.ad_array, 
                                                                                                       dataset.marker_array, dataset.color_array)):
            if idx == 1:
                    axes.plot(ad, 0.5*e.Kraft.max(), marker=marker, markersize=12, color= color, label = label)
            for e_run, ad, marker, color in zip(e_treatment, e_treatment_ad, e_treatment_marker, e_treatment_color):
                if color == None:
                        color = dataset.color_array[idx]
                e = e_run.ex_type
                axes = p.subplot(111)
                axes.plot(ad, 0.5*e.Kraft.max(), marker= marker, markersize=12, color= color)
                

    def decorate(self, axes):
        axes.axis([0., 8, 0., 10]) 
        axes.grid(b=True, which='major', color='gray', linestyle='-', linewidth = .5,)
        axes.set_xlabel('a/d [-]', fontsize = 16)
        axes.set_ylabel('Querkraft [kN]', fontsize = 16)
        axes.tick_params(axis='both', which='major', labelsize=16)
        first_legend = axes.legend(loc=1, markerscale=1., fontsize = 16, numpoints=1,)
        newlegend = p.gca().add_artist(first_legend)
                                       
        markertriangle, = axes.plot([],'^', color='k', markersize=12, label='Querkraftversagen')
        markercircle, =  axes.plot([],'o', color='k', markersize=12, label='Biegeversagen')
        markersquare, = axes.plot([],"s", color='k', markersize=12, label='uneindeutiges Versagen')
        handles = [markertriangle, markercircle, markersquare]
        labels = [h.get_label() for h in handles] 
        axes.legend(handles = handles, markerscale=1., fontsize = 16, numpoints=1, loc=2)    


                 

if __name__ == '__main__':

    pw = PlotQW()

    ax = pw.figure()
    pw.plot(threePT, ax, color='b', label='3-Punkt-Schubversuch')
    pw.plot(fourPT, ax, color='r',  label='4-Punkt-Schubversuch')
    pw.decorate(ax)
    
    p.show()
