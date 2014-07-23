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

    test_files = ['BTT-4c-2cm-TU-0-V03_MxN2.DAT']

    test_file_path = os.path.join(simdb.exdata_dir,
                             'bending_tensile_test',
                             '2014-06-12_BTT-4c-2cm-0-TU_MxN2')

    e_list = [ExpBTTDB(data_file=os.path.join(test_file_path, test_file),
                    delta_t_aramis=5)
           for test_file in test_files]
    for e in e_list:
        e.process_source_data()
        print 'w_idx_cut', e.t[e.w_cut_idx]
        print 'F_max1', e.F_max1

    p.subplot(221)
    for e in e_list:
        p.plot(e.t_aramis, e.F_t_aramis, color='blue', label='F')
        p.plot(e.t_aramis, e.N_t_aramis, color='red', label='N')
    p.legend()

    p.subplot(222)

    for e in e_list:
        p.plot(e.t_aramis, e.eps_t_aramis[0], color='blue', label='eps_tension')
        p.plot(e.t_aramis, e.eps_t_aramis[1], color='green', label='eps_compression')
        p.twinx()
        p.plot(e.t_aramis, e.w_t_aramis, color='black', label='w')
        p.legend()

    p.subplot(223)

    for e in e_list:
        p.plot(e.t_cut_asc, e.M, color='blue', label='eps_tension')
        p.twinx()
        p.plot(e.t_cut_asc, e.N_cut_asc, color='black', label='w')
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

        mid_idx = ac.d_ux_arr2.shape[1] / 2
        eps_range = 10
        eps = np.mean(ac.d_ux_arr2[:, mid_idx - eps_range:mid_idx + eps_range], axis=1)
        y = AD.y_arr_undeformed[:, mid_idx]
        p.plot(eps, y)
        p.show()
