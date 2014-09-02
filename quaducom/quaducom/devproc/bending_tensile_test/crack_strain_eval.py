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

    test_files = ['BTT-4c-2cm-TU-0-V02_MxN2.DAT']

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
        AUI.aramis_data.evaluated_step_idx = e.n_steps - 2
        AUI.configure_traits()
