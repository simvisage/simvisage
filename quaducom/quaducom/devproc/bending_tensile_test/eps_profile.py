'''
Created on Jul 1, 2014

Script evaluating the mxn interaction between
tensile and bending load
'''

if __name__ == '__main__':

    from exp_btt_db import ExpBTTDB
    from matresdev.db.simdb import SimDB
    from aramis_cdt import AramisInfo, AramisUI, AramisFieldData, AramisCDT, AramisUI
    simdb = SimDB()
    import os
    import numpy as np

    import pylab as p


    test_files = ['BTT-6c-2cm-TU-0-V09_MxN2.DAT',
                  ]

    test_file_path = os.path.join(simdb.exdata_dir,
                             'bending_tensile_test',
                             '2014-06-12_BTT-6c-2cm-0-TU_MxN2')

    e_list = [ExpBTTDB(data_file=os.path.join(test_file_path, test_file),
                    delta_t_aramis=5)
           for test_file in test_files]
    for e in e_list:
        e.process_source_data()
        # print 'w_idx_cut', e.t[e.w_cut_idx]
        # print 'F_max1', e.F_max1
        # print 'crack filter', e.crack_filter_avg

    p.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.2, hspace=0.2)

    p.subplot(221)
    for e in e_list:
        p.plot(e.t_aramis_cut, e.N_t_aramis, color='blue', label='N')
        p.ylim(0, 50)
        p.grid()
        p.xlabel('t [sec]')
        p.ylabel('N [kN]')
        p.legend(loc=2)
        p.twinx()
        p.plot(e.t_aramis_cut, e.F_t_aramis, color='red', label='F')
        p.ylabel('F [kN]')
        p.ylim(0, 6)
        p.title(test_files)
        # p.xlim(0, 950)
        p.legend()

    p.subplot(222)
    for e in e_list:
        p.plot(e.t_aramis_cut, e.eps_t_aramis[0] * 1000, color='grey', label='strain_max')
        p.plot(e.t_aramis_cut, e.eps1_t_aramis * 1000, color='green', label='strain_1re')
        p.plot(e.t_aramis_cut, e.eps_t_aramis[1] * 1000, color='black', label='strain_min')
        # p.xlim(0, 500)
        p.ylim(-5, 30)
        p.grid()
        p.xlabel('t [sec]')
        p.ylabel('strain [1*E-3]')
        p.legend(loc=2)
        p.twinx()
        p.plot(e.t_aramis_cut, e.w_t_aramis, color='darkred', label='w')
        p.ylim(0, 10)
        # p.xlim(0, 950)
        p.ylabel('w [mm]')
        p.legend(loc=1)

        print '-----------------------------n_steps', e.aramis_info.number_of_steps
        print '-----------------------------n_steps_cut1', len(e.t_aramis)
        print 'max tension strain', max(e.eps_t_aramis[0] * 1000)
        print 'min compression strain', min(e.eps_t_aramis[1] * 1000)
        print 'max tension strain in first reinforcement layer', max(e.eps1_t_aramis * 1000)


    p.subplot(223)
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
        p.legend(ncol=3)

    p.subplot(224)
    for e in e_list:
        aramis_file_path = e.get_cached_aramis_file('Xf15s3-Yf15s3')
        AI = AramisInfo(data_dir=aramis_file_path)
        ad = AramisFieldData(aramis_info=AI, integ_radius=3)

        max_step = e.n_steps
        a = e.crack_bridge_strain_all
        n_fa = ad.d_ux.shape[0]
        h = np.linspace(e.pos_fa[0], e.pos_fa[1], num=n_fa)
        # print 'h', h

        for step in range(0, max_step, 10):

            ad.current_step = step

            if a == None:
                mid_idx = ad.d_ux.shape[1] / 2
                eps_range = 3
                eps = np.mean(ad.d_ux[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
                # print 'eps', eps
                p.title('strain in the middle of the measuring field')

            else:
                ux = ad.ux_arr
                x_und = ad.x_arr_0
                idx_border1 = e.idx_failure_crack[1]
                idx_border2 = e.idx_failure_crack[2]
                eps_range = 1
                ux1 = np.mean(ux[:, idx_border1 - eps_range: idx_border1 + eps_range ], axis=1)
                ux2 = np.mean(ux[:, idx_border2 - eps_range: idx_border2 + eps_range ], axis=1)
                x_01 = np.mean(x_und [:, idx_border1 - eps_range: idx_border1 + eps_range ], axis=1)
                x_02 = np.mean(x_und [:, idx_border2 - eps_range: idx_border2 + eps_range ], axis=1)
                eps = (ux2 - ux1) / (x_02 - x_01)

                p.title('strain in the failure crack')

            x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
            # print 'x', x
            eps_ed_up = x + eps[-1]
            # print 'eps_ed_up', eps_ed_up
            eps_ed_lo = eps[0] - x
            # print 'eps_ed_lo', eps_ed_lo
            eps_to1 = np.append(eps, eps_ed_lo)
            eps_to2 = np.append(eps_ed_up, eps_to1)
            # print 'eps_to2', eps_to2

            h_1 = np.append(h, 0)
            h_2 = np.append(20, h_1)
            # print 'h_2', h_2

            eps_rev = eps_to2[::-1]

            step_time = e.t_aramis[step]
            p.plot(eps_rev * 1000, h_2, label='%i' % step_time)
            p.xlim(-5, 25)
            p.ylim(0, 20)
            p.xlabel('strain [1*E-3]')
            p.ylabel('h [mm]')
            # p.legend(bbox_to_anchor=(0.66, 0.02), borderaxespad=0., ncol=2, loc=3)
            p.legend(bbox_to_anchor=(0.99, 0.98), borderaxespad=0., ncol=2, loc=1)


        p.show()

        # AUI = AramisUI(aramis_info=e.aramis_info)
        # AUI.configure_traits()
