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

    test_files = ['BTT-4c-2cm-TU-0-V05_MxN2.DAT']

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
        print 'crack filter', e.crack_filter_avg
        AUI = AramisUI(aramis_info=e.aramis_info)
        print e.n_steps
        AUI.aramis_cdt.integ_radius = 5
        AUI.aramis_data.evaluated_step_idx = e.n_steps
        AUI.aramis_cdt.crack_detect_idx = e.n_steps
        AUI.aramis_cdt.ddd_ux_avg_threshold = -1e-4
        AUI.aramis_cdt.ddd_ux_threshold = -1e-4

        x = AUI.aramis_data.x_arr_undeformed
        y = AUI.aramis_data.y_arr_undeformed
        print 'Measuring field length (bottom, top) =', x[0, -1] - x[0, 0], x[-1, -1] - x[-1, 0]
        print 'Measuring field height (left, right) =', y[0, 0] - y[-1, 0], y[0, -1] - y[-1, -1]
        xs = 3  # step in x-direction
        dist = x[:, 1:] - x[:, :-1]
        print 'Mean, std, min, max of facets center distance [mm]: ', np.mean(dist), np.std(dist), np.min(dist), np.max(dist)
        print 'Mean and std of pixel size [mm]', np.mean(dist / xs), np.std(dist / xs)
        AUI.configure_traits()
