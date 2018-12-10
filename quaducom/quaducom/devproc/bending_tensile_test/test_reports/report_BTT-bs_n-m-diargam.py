'''
Created on Aug 26, 2013

@author: alexander
'''
import numpy as np

import os

from quaducom.devproc.format_plot import format_plot
from matplotlib.font_manager import FontProperties

# from matresdev.db.exdb.ex_run import ExRun
import pylab as p

# from matresdev.db.exdb.ex_run_view import \
#    ExRunView
#
from matresdev.db.simdb import \
    SimDB

from matresdev.db.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

simdb = SimDB()

if __name__ == '__main__':

    fig = p.figure(facecolor='white')
    fig.set_size_inches(8, 6)
    fig.set_size_inches(6, 6)

    ax = fig.add_subplot(111)

    # shift xlabes to top axis
    #
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_smart_bounds(True)
    ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('top')
    ax.yaxis.set_ticks_position('left')
    ax.axhline(linewidth=1., color='k')
    ax.axvline(linewidth=1., color='k')

    fontsize = 18
    from matplotlib.font_manager import FontProperties
    font = FontProperties()
    font.set_family('serif')
    font.set_style('normal')
    font.set_size(fontsize)

    m = np.array([0, 0.38])
    n = np.array([46.62, 1.26])
    ax.plot(m, n, color='grey', linestyle='--', linewidth=1.5)

#    p.annotate('$M$', xy=(0, 46.62), xytext=(0.1, 10), textcoords='offset points', fontsize=12)
    xoffset = 0.002
    yoffset = 2.6
    # N_max
    p.plot([0.00], [46.62], 'kD')
    p.text(0.00 + xoffset, 46.62 + yoffset, '$N_\mathrm{max}$', fontproperties=font)
    p.plot([0.00], [43.08], 'kx')
    p.plot([0.00], [50.15], 'kx')
    # I
    p.plot([0.14], [33.26], 'kD')
    p.text(0.14 + xoffset, 33.26 + yoffset, '$I$', fontproperties=font)
    p.plot([0.168], [33.25], 'kx')
    p.plot([0.074], [33.26], 'kx')
    p.plot([0.171], [33.26], 'kx')
    # II
    p.plot([0.24], [22.26], 'kD')
    p.text(0.24 + xoffset, 22.26 + yoffset, '$II$', fontproperties=font)
    p.plot([0.233], [22.25], 'kx')
    p.plot([0.227], [22.26], 'kx')
    p.plot([0.263], [22.26], 'kx')
    # III
    p.plot([0.32], [11.26], 'kD')
    p.text(0.32 + xoffset, 11.26 + yoffset, '$III$', fontproperties=font)
    p.plot([0.364], [11.26], 'kx')
    p.plot([0.344], [11.26], 'kx')
    p.plot([0.256], [11.26], 'kx')
    # M_max
    p.text(0.38 + xoffset, 1.26 + yoffset, '$M_\mathrm{max}$', fontproperties=font)
    p.plot([0.385], [1.25], 'kx')
    p.plot([0.373], [1.26], 'kx')
    p.plot([0.395], [1.26], 'kx')
    p.plot([0.38], [1.26], 'kD', label='mean values')

    p.legend(loc=4, frameon=True, numpoints=1, markerscale=1.1)

    # ranges and labels
    #
#    xlim = 0.4
#    ylim = 50

    p.axis([-0.01, 0.4, 50.8, 0])
    p.xticks([0, 0.1, 0.2, 0.3, 0.4])
    p.yticks([0, 10, 20, 30, 40, 50])

    xlabel = 'bending moment $M$ [kNm]'
    ylabel = 'normal force $N$ [kN]'
    format_plot(p, fontsize=fontsize, xformat="%.2f", xlabel=xlabel, ylabel=ylabel)  # , xlim=xlim, ylim=ylim)

    ax = p.gca()
    ax.xaxis.grid(True, which='major')
    ax.yaxis.grid(True, which='major')


    save_fig_to_file = True
    test_series_name = 'BTT-6c_NXM2'
    if save_fig_to_file:
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
        filename = os.path.join(test_series_dir, 'n-m-diagr.pdf')
        p.savefig(filename, format='pdf')
        filename = os.path.join(test_series_dir, 'n-m-diagr.png')
        p.savefig(filename, format='png')
        print('figure saved to file %s' % (filename))

    p.show()
