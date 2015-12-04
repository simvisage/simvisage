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

    path_a = os.path.join(simdb.simdata_dir, 'SimST', 'a.csv')
    path_b = os.path.join(simdb.simdata_dir, 'SimST', 'b.csv')
    path_c = os.path.join(simdb.simdata_dir, 'SimST', 'c.csv')
    path_d = os.path.join(simdb.simdata_dir, 'SimST', 'd.csv')
    path_e = os.path.join(simdb.simdata_dir, 'SimST', 'e.csv')
    path_f = os.path.join(simdb.simdata_dir, 'SimST', 'f.csv')
    path_g = os.path.join(simdb.simdata_dir, 'SimST', 'g.csv')
    path_h = os.path.join(simdb.simdata_dir, 'SimST', 'h.csv')

    path_list = [path_a, path_b, path_c, path_d, path_e, path_f, path_g, path_h]

#     data_arr = np.loadtxt(path_a, delimiter=";")
    for path_i in path_list:
        print 'path_i', path_i[-5]
        data_arr = np.loadtxt(path_i, delimiter=";")
        f = data_arr[:, 0]
        w = data_arr[:, 1]
        ax.plot(w, f, label=path_i[-5])

#     ax.plot(f, w, color='grey', linestyle='--', linewidth=1.5)
        p.legend(loc=4, frameon=True, numpoints=1, markerscale=1.1)

    xlabel = 'deflection $w$ [m]'
    ylabel = 'force $F$ [kN]'

    save_fig_to_file = True
    test_series_name = 'slabtest_pstudy'
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
        filename = os.path.join(test_series_dir, 'pstudy_slabtests_solids.pdf')
        p.savefig(filename, format='pdf')
        filename = os.path.join(test_series_dir, 'pstudy_slabtests_solids.png')
        p.savefig(filename, format='png')
        print 'figure saved to file %s' % (filename)

    p.show()
