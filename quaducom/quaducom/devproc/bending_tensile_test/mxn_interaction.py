'''
Created on Jul 1, 2014

Script evaluating the mxn_simple_script interaction between
tensile and bending load
'''

from .exp_btt_db import ExpBTTDB
from matresdev.db.simdb import SimDB
from aramis_cdt import AramisInfo, AramisUI, AramisFieldData, AramisCDT, AramisUI
# AramisData, AramisBSA,
simdb = SimDB()
import os
import numpy as np

import pylab as p

test_files = [  # 'BTT-6c-2cm-TU-0-V04_MxN2.DAT',
               'BTT-6c-2cm-TU-0-V02_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V03_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V04_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V05_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V06_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V07_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V08_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V09_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V10_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V11_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V12_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V13_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V14_MxN2.DAT',
              # 'BTT-6c-2cm-TU-0-V15_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V01_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V02_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V03_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V04_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V05_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V06_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V07_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V08_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V09_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V10_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V11_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V12_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V13_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V14_MxN2.DAT',
              # 'BTT-4c-2cm-TU-0-V15_MxN2.DAT',
              # 'BT-6c-2cm-0-TU-V4_bs4.DAT'
              ]

test_file_path = os.path.join(simdb.exdata_dir,
                             'bending_tensile_test',
                             '2014-06-12_BTT-6c-2cm-0-TU_MxN2',
                             # '2013-07-09_BT-6c-2cm-0-TU_bs4-Aramis3d'
                             )

res_key = 'Xf15s3-Yf15s3'
res_key = 'Xf19s1-Yf5s4'
res_key = 'Xf19s15-Yf19s15'

e_list = [ExpBTTDB(data_file=os.path.join(test_file_path, test_file),
                   delta_t_aramis=5)
       for test_file in test_files]

def plot_all():

    fig = p.figure()
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.04, top=0.96, wspace=0.25, hspace=0.2)
    axes = p.subplot(231)

    for e in e_list:
        e.process_source_data()
        # print 'step_times', e.aramis_field_data.step_times
        # print 'crack filter', e.crack_filter_avg
        e.process_source_data()
        e.aramis_cdt.ddd_ux_avg_threshold = -5e-4
        e.aramis_cdt.crack_detection_step = 116

        # print 'step_times', e.aramis_field_data.step_times
        # print 'crack filter', e.crack_filter_avg
        axes.plot(e.t_aramis_cut, e.N_t_aramis, color='darkblue', label='N')
        axes.grid()
        axes.set_ylim(0, 55)
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('N [kN]')
        axes.legend(loc=2)
        ax2 = axes.twinx()
        ax2.plot(e.t_aramis_cut, e.F_t_aramis, color='darkred', label='F')
        # ax2.set_xlim(0, 650)
        ax2.set_ylim(0, 5.5)
        ax2.set_ylabel('F [kN]')
        ax2.legend(loc=9)

    axes = p.subplot(232)
    for e in e_list:
        axes.plot(e.t_aramis_cut, e.N_t_aramis, color='darkblue', label='N')
        axes.grid()
        axes.set_ylim(0, 55)
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('N [kN]')
        axes.legend(loc=2)
        ax2 = axes.twinx()
        ax2.plot(e.t_aramis_cut, e.M_t_aramis, color='green', label='M')
        ax2.plot(e.t_aramis_cut, e.MF_t_aramis, color='darkred', label='M0')
        ax2.plot(e.t_aramis_cut, e.MN_t_aramis, color='aqua', label='MII')
        ax2.set_ylim(-0.33, 0.58)
        # ax2.set_xlim(0, 650)
        ax2.set_ylabel('M [kNm]')
        ax2.legend(ncol=3)
        axes.set_title(test_files)

    axes = p.subplot(233)
    for e in e_list:
        axes.plot(e.t_aramis_cut, e.eps_t_aramis[0] * 1000, color='grey', label='cbs_max')
        axes.plot(e.t_aramis_cut, e.eps1_t_aramis * 1000, color='green', label='cbs_1re')
        axes.plot(e.t_aramis_cut, e.eps_t_aramis[1] * 1000, color='black', label='cbs_min')
        axes.grid()
        axes.set_ylim(-6, 30)
        axes.set_xlabel('t [sec]')
        axes.set_ylabel('crack bridge strain [1*E-3]')
        axes.legend(loc=2)
        ax2 = axes.twinx()
        ax2.plot(e.t_aramis_cut, e.w_t_aramis, color='darkred', label='w')
        ax2.set_ylim(0, 10)
        # ax2.set_xlim(0, 650)
        ax2.set_ylabel('w [mm]')
        ax2.legend(loc=9)

        print('-----------------------------n_steps', e.aramis_info.number_of_steps)
        print('-----------------------------n_steps_cut1', len(e.t_aramis))
        print('max tension strain', (e.eps_t_aramis[0] * 1000)[-1])
        print('min compression strain', (e.eps_t_aramis[1] * 1000)[-1])
        print('max tension strain in first reinforcement layer', max(e.eps1_t_aramis * 1000))

    axes = p.subplot(234)
    for e in e_list:
        axes.plot(e.M_t_aramis, e.N_t_aramis, color='green', label='M')
        x = [0, 0.35]
        y = [44, 0]
        axes.plot(x, y, color='grey', linestyle='--', label='M-N-Interaction')
        axes.grid()
        axes.xaxis.tick_top()
        axes.xaxis.set_label_position('top')
        axes.set_ylim(55 , -1)
        axes.set_xlim(-0.01 , 0.4)
        axes.set_xlabel('M [kNm]')
        axes.set_ylabel('N [kN]')
        axes.legend(loc=4)


    axes = p.subplot(235)
    for e in e_list:
        ad = e.aramis_field_data

        max_step = e.n_steps
        a = e.crack_bridge_strain_all
        F_max = np.max(e.F_t_aramis)
        eps_list = []

        for step in range(0, max_step, 1):
            ad.current_step = step

            eps = np.mean(ad.d_ux[:, :], axis=1)
            eps_list.append(np.mean(eps))

        eps_mean = np.array(eps_list, dtype='f')
        sig_c = e.N_t_aramis / (e.A_c * 1000)
        sig_tex = e.N_t_aramis * 1000 / (e.A_tex)

        # print'eps_mean', eps_mean
        # print 'shape_eps_mean', np.shape(eps_mean)
        # print 'N_t_aramis', e.N_t_aramis
        # print 'N_t_aramis', np.shape(e.N_t_aramis)
        # print 'sig_c', sig_c
        # print 'e.A_tex', e.A_tex
        # print  'sig_tex', sig_tex

        print('tensile_text_max tension strain', eps_mean[-1] * 1000)

        # p.plot(eps_mean * 1000, sig_c, label='sig_c')

        if F_max > 1.2:
            ax2 = axes.twinx()
            ax2.plot(e.eps_M[0] * 1000, e.M_t_F, color='grey', label='cbs_max (M)')
            ax2.plot(e.eps_M[1] * 1000, e.M_t_F, color='black', label='cbs_min (M)')
            # ax2.plot(e.eps_M[1] * 1000, e.M_t_F, color='green', label=r'$\Delta$ cbs_min (M)')
            ax2.set_ylim(0, 0.45)
            ax2.set_ylabel('M [kNm]')
            ax2.legend(loc=9)

            if e.N_t_N == []:
                p.plot(1, 1, color='black')
                axes.grid()
                axes.set_ylim(0, 55)
                axes.set_xlabel('crack bridge strain [1*E-3]')
            else:
                axes.plot(e.eps_N * 1000, e.N_t_N, color='darkblue', label='cbs (N)')
                axes.grid()
                axes.set_ylim(0, 55)
                axes.set_xlabel('crack bridge strain [1*E-3]')
                axes.set_ylabel('N [kN]')
                axes.legend(loc=2)

        else:
            axes.plot(eps_mean * 1000, sig_tex, label='sig_tex')

            E_tex = e.ccs.E_tex
            K_III = E_tex
            eps_max = max(eps_mean)
            # print 'K_III', K_III
            eps_lin = np.array([0, eps_max * 1000], dtype='float_')
            sig_lin = np.array([0, eps_max * K_III], dtype='float_')
            axes.plot(eps_lin, sig_lin, color='grey', linestyle='--')
            axes.set_xlim(-1.5, 8)
            axes.set_ylim(0, 1600)
            axes.set_xlabel('strain [1*E-3]')
            axes.set_ylabel('sigma tex [N/mm]')
            axes.legend(loc=2)


    axes = p.subplot(236)
    for e in e_list:
        ad = e.aramis_field_data

        max_step = e.n_steps
        a = e.crack_bridge_strain_all
        n_fa = ad.d_ux.shape[0]
        h = np.linspace(e.pos_fa[0], e.pos_fa[1], num=n_fa)
        # print 'h', h

        # print 'crack filter', e.crack_filter_avg

        for step in range(0, max_step, 10):

            ad.current_step = step

            if a == None:
                mid_idx = ad.d_ux.shape[1] / 2
                eps_range = 3
                eps = np.mean(ad.d_ux[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
            else:
                ux = ad.ux_arr
                x_0 = ad.x_arr_0
                idx_border1 = e.idx_failure_crack[1]
                idx_border2 = e.idx_failure_crack[2]
                eps_range = 2
                ux1 = np.mean(ux[:, idx_border1 - eps_range: idx_border1 + eps_range ], axis=1)
                ux2 = np.mean(ux[:, idx_border2 - eps_range: idx_border2 + eps_range ], axis=1)
                x_0_1 = np.mean(x_0[:, idx_border1 - eps_range: idx_border1 + eps_range ], axis=1)
                x_0_2 = np.mean(x_0[:, idx_border2 - eps_range: idx_border2 + eps_range ], axis=1)

                # ux1 = ux[:, 0]
                # ux2 = ux[:, -1]
                # x_0_1 = x_0[:, 0]
                # x_0_2 = x_0[:, -1]
                # print 'ux1', ux1
                # eps = np.mean(ac.d_ux_arr[:, idx_border1:idx_border2], axis=1)

                eps = (ux2 - ux1) / (x_0_2 - x_0_1)

                # print 'eps', eps
                p.title('crack bridge strain in the failure crack')

            x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
            eps_ed_up = x + eps[-1]
            eps_ed_lo = eps[0] - x
            eps_to1 = np.append(eps, eps_ed_lo)
            eps_to2 = np.append(eps_ed_up, eps_to1)
            # print 'eps_ed_lo', eps_ed_lo
            # print 'eps_ed_up', eps_ed_up
            # print 'x', x
            # print 'eps_to2', eps_to2

            h_1 = np.append(h, 0)
            h_2 = np.append(20, h_1)
            # print 'h_2', h_2

            eps_rev = eps_to2[::-1]

            step_time = e.t_aramis[step]
            axes.plot(eps_rev * 1000, h_2, label='%i' % step_time)
            axes.grid()
            axes.set_xlim(-6, 30)
            axes.set_ylim(0, 20)
            axes.set_xlabel('crack bridge strain [1*E-3]')
            axes.set_ylabel('h [mm]')
            # p.legend(bbox_to_anchor=(0.55, 0.02), borderaxespad=0., ncol=2, loc=3)
            p.legend(bbox_to_anchor=(0.99, 0.98), borderaxespad=0., ncol=2, loc=1)


    '''p.subplot(235)
    for e in e_list:
        p.plot(e.t_aramis_cut, e.N_t_aramis, color='blue', label='N')
        p.xlim(0, 500)
        p.ylim(0, 50)
        p.grid()
        p.xlabel('t [sec]')
        p.ylabel('N [kN]')
        p.legend(loc=2)
        p.twinx()
        p.plot(e.t_aramis_cut, e.M_t_aramis, color='black', label='M')
        p.plot(e.t_aramis_cut, e.MF_t_aramis, color='red', label='M0')
        p.plot(e.t_aramis_cut, e.MN_t_aramis, color='aqua', label='MII')
        p.ylim(-0.25, 0.5)
        # p.xlim(0, 950)
        p.ylabel('M [kNm]')
        p.legend(ncol=3)'''


#        AUI = AramisUI(aramis_info=e.aramis_info)
#        AUI.configure_traits()


if __name__ == '__main__':

    fig = p.figure()
    plot_all()
    p.show()
