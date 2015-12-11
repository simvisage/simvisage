'''
Created on Jul 1, 2014

Script evaluating the mxn interaction between
tensile and bending load
'''

if __name__ == '__main__':

    from exp_btt_db import ExpBTTDB
    from matresdev.db.simdb import SimDB
    from aramis_cdt import AramisUI
    simdb = SimDB()
    import os
    import numpy as np

    test_files = [
        'BT-1C-55mm-0-3300EP-V2_S3P2(11)-Aramis2d.DAT',
    ]

    test_file_path = os.path.join(simdb.exdata_dir,
                                  'bending_tests', 'four_point',
                                  '2015-09-02_BT-1C-55mm-0-3300SBR_cyc-Aramis2d',
                                  )

    #res_key = 'Xf15s13-Yf15s13'
    res_key = 'Xf15s1-Yf15s4'

    e_list = [ExpBTTDB(data_file=os.path.join(test_file_path, test_file),
                       aramis_resolution_key=res_key,
                       delta_t_aramis=5)
              for test_file in test_files]
    for e in e_list:
        e.process_source_data()
        print 'crack filter', e.crack_filter_avg
        print e.n_steps
        e.aramis_field_data.integ_radius = 3
#         e.aramis_field_data.current_step = e.n_steps
#         e.aramis_cdt.crack_detection_step = e.n_steps
        e.aramis_cdt.ddd_ux_avg_threshold = -0.5e-3
        e.aramis_cdt.ddd_ux_threshold = -0.5e-3

        x = e.aramis_field_data.x_arr_0
        y = e.aramis_field_data.y_arr_0
        print 'Measuring field length (bottom, top) =', x[0, -1] - x[0, 0], x[-1, -1] - x[-1, 0]
        print 'Measuring field height (left, right) =', y[0, 0] - y[-1, 0], y[0, -1] - y[-1, -1]
        xs = 3  # step in x-direction
        dist = x[:, 1:] - x[:, :-1]
        print 'Mean, std, min, max of facets center distance [mm]: ', np.mean(dist), np.std(dist), np.min(dist), np.max(dist)
        print 'Mean and std of pixel size [mm]', np.mean(dist / xs), np.std(dist / xs)
        aui = AramisUI(aramis_info=e.aramis_info,
                       aramis_data=e.aramis_field_data,
                       aramis_cdt=e.aramis_cdt)
        aui.configure_traits()
