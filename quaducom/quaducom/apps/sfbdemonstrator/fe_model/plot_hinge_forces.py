'''
Created on Aug 26, 2013

@author: alexander
'''
import numpy as np

from os.path import join

from matresdev.db.simdb import \
    SimDB

import os

simdb = SimDB()

from quaducom.devproc.tensile_test.dog_bone.test_reports import format_plot
# from matplotlib.font_manager import FontProperties
# font = FontProperties()

if __name__ == '__main__':

    import pylab as p
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2])

    fig = p.figure(facecolor='white')
    fig.set_size_inches(15, 5)
    #
    filename = 'exc_V_ip_equal_100cm_7m_hf_lc_w_asym'
    N = 14
    width = 0.15
    step = 1.0  # m
    #
    filename = 'exc_V_ip_equal_25cm_7m_hf_lc_w_asym'
    N = 56
    width = 0.15
    step = 0.25  # m
    #
    path = filename + '.csv'
    #
    data_file = open(path)
    headings = data_file.readline()
    data_arr = np.loadtxt(path, skiprows=1, delimiter=";")
    #
    X_sym = data_arr[:N, 0]
    Y_sym = data_arr[:N, 1]
    Z_sym = data_arr[:N, 2]
    N_ip_sym = data_arr[:N, 3] * 1000.
    V_ip_sym = data_arr[:N, 4] * 1000.
    V_op_sym = data_arr[:N, 5] * 1000.
    print('X_sym.shape', X_sym.shape)
    #
    X_r01 = data_arr[N:, 0]
    Y_r01 = data_arr[N:, 1]
    Z_r01 = data_arr[N:, 2]
    N_ip_r01 = data_arr[N:, 3] * 1000.
    V_ip_r01 = data_arr[N:, 4] * 1000.
    V_op_r01 = data_arr[N:, 5] * 1000.
    print('Y_r01.shape', Y_r01.shape)
    #

    do = 'N_ip'
#    do = 'V_op'

    fontsize = 20
    color = 'grey'

    if do == 'N_ip':
        ax1 = p.subplot(gs[1])
        ax1.set_title('Symmetrieachse')
        ind = np.arange(N)
        ax1.plot(ind * step + step / 2. - width / 2., np.zeros_like(X_sym), linestyle='--', color=color)
        ax1.bar(ind * step + step / 2. - width / 2., N_ip_sym, width, color='grey', edgecolor='none', linewidth=0)
        format_plot(p, fontsize=fontsize, xlabel='$X\,\mathrm{[m]}$', ylabel='$N_\mathrm{ip}\,\mathrm{[kN]}$')
        #
        ax2 = p.subplot(gs[0], sharey=ax1)
        ax2.set_title('Schnitt A-A')
        ind = np.arange(N / 2)
        ax2.plot(ind * step + step / 2. - width / 2., np.zeros_like(Y_r01), linestyle='--', color=color)
        ax2.bar(ind * step + step / 2. - width / 2., N_ip_r01, width, color='grey', edgecolor='none', linewidth=0)
        format_plot(p, fontsize=fontsize, xlabel='$Y\,\mathrm{[m]}$', ylabel='$N_\mathrm{ip}\,\mathrm{[kN]}$')

    if do == 'V_op':
        ax2 = p.subplot(gs[0])
        ax2.set_title('Schnitt A-A')
        ind = np.arange(N / 2)
        ax2.plot(ind * step + step / 2. - width / 2., np.zeros_like(Y_r01), linestyle='--', color=color)
        ax2.bar(ind * step + step / 2. - width / 2., V_op_r01, width, color='grey', edgecolor='none', linewidth=0)
        format_plot(p, fontsize=fontsize, yformat="%.1f", xlabel='$Y\,\mathrm{[m]}$', ylabel='$V_\mathrm{op}\,\mathrm{[kN]}$')
        #
        ax1 = p.subplot(gs[1], sharey=ax2)
        ax1.set_title('Symmetrieachse')
        ind = np.arange(N)
        ax1.plot(ind * step + step / 2. - width / 2., np.zeros_like(X_sym), linestyle='--', color=color)
        ax1.bar(ind * step + step / 2. - width / 2., V_op_sym, width, color=color, edgecolor='none', linewidth=0)
        format_plot(p, fontsize=fontsize, yformat="%.1f", xlabel='$X\,\mathrm{[m]}$', ylabel='$V_\mathrm{op}\,\mathrm{[kN]}$')

    p.ylim([-20., 20.])
#    fig.tight_layout()
    simdata_dir = os.path.join(simdb.simdata_dir, 'show_results')
    p.savefig(os.path.join(simdata_dir, filename + '_' + do), format='pdf', dpi=600.)
    print('png saved to file ' + filename)
    p.show()


#    ax1 = fig.add_subplot(1, 2, 1)
#    ax1.set_xlabel('$X\,\mathrm{[m]}$')
#    ax1.set_ylabel('$N_\mathrm{ip}\,\mathrm{[kN]}$')
