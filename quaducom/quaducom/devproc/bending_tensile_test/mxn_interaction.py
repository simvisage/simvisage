'''
Created on Jul 1, 2014

Script evaluating the mxn interaction between
tensile and bending load
'''

if __name__ == '__main__':

    from exp_btt_db import ExpBTTDB
    from matresdev.db.simdb import SimDB
    from aramis_cdt import AramisInfo, AramisData, AramisBSA, AramisUI
    simdb = SimDB()
    import os
    import numpy as np

    import pylab as p


    test_files = ['BTT-4c-2cm-TU-0-V03_MxN2.DAT']

    test_file_path = os.path.join(simdb.exdata_dir,
                             'bending_tensile_test',
                             '2014-06-12_BTT-4c-2cm-0-TU_MxN2')

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
        p.plot(e.t_aramis_cut, e.eps_t_aramis[0] * 1000, color='grey', label='eps_tension')
        p.plot(e.t_aramis_cut, e.eps1_t_aramis * 1000, color='black', label='eps_tension_1re')
        p.plot(e.t_aramis_cut, e.eps_t_aramis[1] * 1000, color='green', label='eps_compression')
        p.xlim(0, 500)
        p.ylim(-4, 20)
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
        print 'max compression strain', max(e.eps_t_aramis[1] * 1000)
        print 'max tension strain in first reinforcement layer', max(e.eps1_t_aramis * 1000)


    p.subplot(223)
    for e in e_list:
        p.plot(e.t_aramis_cut, e.N_t_aramis, color='blue', label='N')
        p.xlim(0, 500)
        p.ylim(0, 50)
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
        p.legend()

    p.subplot(224)
    for e in e_list:
        aramis_file_path = e.get_cached_aramis_file('Xf15s3-Yf15s3')

        # @todo:
        # 1) define the mapping between astage_idx(t).
        # 2) define a method evaluating strain profile within the object.

        AI = AramisInfo(data_dir=aramis_file_path)
        AD = AramisData(aramis_info=AI, evaluated_step_idx=60)

        ac = AramisBSA(aramis_info=AI,
                        aramis_data=AD,
                        integ_radius=10)

        max_step = e.n_steps

        for step in range(0, max_step, 5):

            AD.evaluated_step_idx = step

            # y = AD.y_arr_undeformed[:, mid_idx]

            a = e.crack_bridge_strain_all
            n_fa = ac.d_ux_arr.shape[0]
            h = np.linspace(e.pos_fa[0], e.pos_fa[1], num=n_fa)

            if a == None:
                mid_idx = ac.d_ux_arr.shape[1] / 2
                eps_range = 3
                eps = np.mean(ac.d_ux_arr[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)

                step_time = e.t_aramis[step]
                p.plot(eps * 1000, h, label='%i' % step_time)
                p.xlim(-5, 25)
                p.ylim(0, 20)
                p.xlabel('strain [1*E-3]')
                p.ylabel('h [mm]')
                p.legend(bbox_to_anchor=(0.66, 0.02), borderaxespad=0., ncol=2, loc=3)
                p.title('strain in the middle of the measuring field')

            else:
                idx_border1 = e.idx_failure_crack[1]
                idx_border2 = e.idx_failure_crack[2]

                eps = np.mean(ac.d_ux_arr[:, idx_border1:idx_border2], axis=1)

                step_time = e.t_aramis[step]
                p.plot(eps * 1000, h, label='%i' % step_time)
                p.xlim(-5, 25)
                p.ylim(0, 20)
                p.xlabel('strain [1*E-3]')
                p.ylabel('h [mm]')
                p.legend(bbox_to_anchor=(0.66, 0.02), borderaxespad=0., ncol=2, loc=3)
                p.title('strain in the failure crack')

        p.show()

        # AUI = AramisUI(aramis_info=e.aramis_info)
        # AUI.configure_traits()
