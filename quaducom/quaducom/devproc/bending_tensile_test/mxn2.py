'''
Created on Jul 1, 2014

Script evaluating the mxn interaction between
tensile and bending load
'''

if __name__ == '__main__':

    from exp_btt_db import ExpBTTDB
    from matresdev.db.simdb import SimDB
    from aramis_cdt import AramisInfo, AramisFieldData, AramisUI, AramisCDT
    simdb = SimDB()
    import os
    import numpy as np

    import pylab as p


    test_files = ['BTT-6c-2cm-TU-0-V06_MxN2.DAT']


    test_file_path = os.path.join(simdb.exdata_dir,
                             'bending_tensile_test',
                             '2014-06-12_BTT-6c-2cm-0-TU_MxN2')

    e_list = [ExpBTTDB(data_file=os.path.join(test_file_path, test_file),
                    delta_t_aramis=5)
           for test_file in test_files]
    for e in e_list:
        e.process_source_data()


    p.subplot(221)
    for e in e_list:
        if e.N_t_N == []:
            p.plot(1, 1, color='black')
        else:
            p.tick_params
            p.plot(e.eps_N * 1000, e.N_t_N, color='black', label='eps_N')
            p.ylim(0, 50)
            p.grid()
            p.xlabel('strain [1*E-3]')
            p.ylabel('N [kN]')
            p.legend(loc=1)
            p.title(test_files)

    # eps(N)
    p.subplot(222)
    for e in e_list:
        aramis_file_path = e.get_cached_aramis_file('Xf15s3-Yf15s3')
        AI = AramisInfo(data_dir=aramis_file_path)
        ad = AramisFieldData(aramis_info=AI,
                             integ_radius=3)
        absa = AramisCDT(aramis_info=AI,
                         aramis_data=ad)

        max_step = e.n_steps
        a = e.crack_bridge_strain_all
        n_fa = ad.d_ux.shape[0]
        h = np.linspace(e.pos_fa[0], e.pos_fa[1], num=n_fa)
        t_N = e.t_N_arr

        if len(t_N) == 0:
            N_end_idx = 0
        else:
            N_end_idx = np.shape(t_N) [0]

        for step in range(0, N_end_idx, 5):
            ad.current_step = step

            if a == None:
                mid_idx = ad.d_ux.shape[1] / 2
                eps_range = 3
                eps = np.mean(ad.d_ux[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
                p.title('eps(N) in the middle of the measuring field')
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

                eps = (ux2 - ux1) / (x_0_2 - x_0_1)
                # eps = np.mean(ad.d_ux[:, idx_border1:idx_border2], axis=1)
                p.title('eps(N) in the failure crack')

            x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
            eps_ed_up = x + eps[-1]
            eps_ed_lo = eps[0] - x
            eps_to1 = np.append(eps, eps_ed_lo)
            eps_to2 = np.append(eps_ed_up, eps_to1)

            h_1 = np.append(h, 0)
            h_2 = np.append(20, h_1)

            eps_rev = eps_to2[::-1]
            step_time = e.t_aramis[step]
            p.plot(eps_rev * 1000, h_2, label='%i' % step_time)
            p.xlim(-5, 25)
            p.ylim(0, 20)
            p.xlabel('strain [1*E-3]')
            p.ylabel('h [mm]')
            # p.legend(bbox_to_anchor=(0.66, 0.02), borderaxespad=0., ncol=2, loc=3)
            p.legend(bbox_to_anchor=(0.99, 0.98), borderaxespad=0., ncol=2, loc=1)


    p.subplot(223)
    for e in e_list:
        if e.F_t_F == []:
            p.plot(1, 1, color='black')
        else:
            p.tick_params
            p.plot(e.eps_M[0] * 1000, e.M_t_F, color='grey', label='eps_M_max')
            p.plot(e.eps_M[1] * 1000, e.M_t_F, color='black', label='eps_M_min')
            p.ylim(0, 0.4)
            p.grid()
            p.xlabel('strain [1*E-3]')
            p.ylabel('M [kNm]')
            p.legend(loc=1, ncol=2)
            p.title(test_files)

    #  eps(M)
    p.subplot(224)
    for e in e_list:
        aramis_file_path = e.get_cached_aramis_file('Xf15s3-Yf15s3')
        AI = AramisInfo(data_dir=aramis_file_path)
        ad = AramisFieldData(aramis_info=AI,
                             integ_radius=3)
        absa = AramisCDT(aramis_info=AI,
                         aramis_data=ad)

        max_step = e.n_steps
        a = e.crack_bridge_strain_all
        n_fa = ad.d_ux.shape[0]
        h = np.linspace(e.pos_fa[0], e.pos_fa[1], num=n_fa)
        t_N = e.t_N_arr
        F_beg_idx = e.F_beg_idx

        if F_beg_idx == []:
            F_beg_idx = max_step

        elif len(t_N) == 0:
            F_beg_idx = 0

        print 'F_beg_idx', F_beg_idx

        # else:
            # F_beg_idx_n = np.shape(t_N) [0]
            # print 'F_beg_idx_n', F_beg_idx_n


        for step in range(F_beg_idx , max_step, 5):
            ad.current_step = step

            if a == None:
                mid_idx = ad.d_ux.shape[1] / 2
                eps_range = 3
                eps = np.mean(ad.d_ux[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
                p.title('eps(M)in the middle of the measuring field')
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
                eps = np.mean(ad.d_ux[:, idx_border1:idx_border2], axis=1)
                p.title('eps(M) in the failure crack')

            x = ((20 - h[-1]) * (eps[0] - eps[-1])) / (h[0] - h[-1])
            eps_ed_up = x + eps[-1]
            eps_ed_lo = eps[0] - x
            eps_to1 = np.append(eps, eps_ed_lo)
            eps_to2 = np.append(eps_ed_up, eps_to1)

            h_1 = np.append(h, 0)
            h_2 = np.append(20, h_1)

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



    ''''p.subplot(221)
    for e in e_list:
        p.plot(e.cu_t_aramis, e.M_t_aramis, color='black')
        p.ylim(0, 0.4)
        p.xlabel('curvature [1/mm]')
        p.ylabel('M [kNm]')
        p.legend(loc=4)
    for e in e_list:
        p.plot(e.t_aramis_cut, e.x_t_aramis, color='black')
        p.ylim(-1, 21)
        p.xlabel('t [step]')
        p.ylabel('x [mm]')
        p.legend(loc=4)'''

        # AUI = AramisUI(aramis_info=e.aramis_info)
        # AUI.configure_traits()
