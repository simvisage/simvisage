'''
Created on Jul 1, 2014

Script evaluating the mxn interaction between
tensile and bending load
'''

if __name__ == '__main__':

    from exp_btt_db import ExpBTTDB
    from matresdev.db.simdb import SimDB
    from aramis_cdt import AramisInfo, AramisData, AramisBSA
    simdb = SimDB()
    import os
    import numpy as np

    import pylab as p


    test_files = ['BTT-6c-2cm-TU-0-V02_MxN2.DAT',
                  # 'BTT-4c-2cm-TU-0-V09_MxN2.DAT',
                  # 'BTT-4c-2cm-TU-0-V13_MxN2.DAT',
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

        # @todo:
        # 1) define the mapping between astage_idx(t).
        # 2) define a method evaluating strain profile within the object.

        AI = AramisInfo(data_dir=aramis_file_path)
        # AD = AramisData(aramis_info=AI, evaluated_step_idx=60)
        AD = AramisData(aramis_info=AI)

        ac = AramisBSA(aramis_info=AI,
                        aramis_data=AD,
                        integ_radius=10)

        max_step = e.n_steps
        a = e.crack_bridge_strain_all
        n_fa = ac.d_ux_arr.shape[0]
        h = np.linspace(e.pos_fa[0], e.pos_fa[1], num=n_fa)
        # print 'h', h

        for step in range(0, max_step, 10):

            AD.evaluated_step_idx = step

            if a == None:
                mid_idx = ac.d_ux_arr.shape[1] / 2
                eps_range = 3
                eps = np.mean(ac.d_ux_arr[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
                # print 'eps', eps
                p.title('strain in the middle of the measuring field')

            else:
                ux = AD.ux_arr
                x_und = AD.x_arr_undeformed
                idx_border1 = e.idx_failure_crack[1]
                idx_border2 = e.idx_failure_crack[2]
                eps_range = 1
                ux1 = np.mean(ux[:, idx_border1 - eps_range: idx_border1 + eps_range ], axis=1)
                ux2 = np.mean(ux[:, idx_border2 - eps_range: idx_border2 + eps_range ], axis=1)
                x_und1 = np.mean(x_und [:, idx_border1 - eps_range: idx_border1 + eps_range ], axis=1)
                x_und2 = np.mean(x_und [:, idx_border2 - eps_range: idx_border2 + eps_range ], axis=1)
                eps = (ux2 - ux1) / (x_und2 - x_und1)

                u_left, u_right = ux[0, (idx_border1, idx_border2) ]
                print 'u', u_left, u_right
                x_left, x_right = x_und[0, (idx_border1, idx_border2) ]
                print 'x', x_left, x_right
                eps_cb = np.fabs(u_right - u_left) / np.fabs(x_right - x_left)
                print 'eps', eps_cb
                # print 'x_und[0, idx_border1:idx_border2', x_und[0, idx_border1:idx_border2]
                # print 'ux1', ux1
                # print 'ux[:, idx_border2]', ux[:, idx_border2]
                # print 'x_und [:, idx_border1]', x_und [:, idx_border1]
                # print 'x_und [:, idx_border2]', x_und [:, idx_border2]
                # print 'eps', eps

                # eps = np.mean(ac.d_ux_arr[:, idx_border1:idx_border2], axis=1)
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
