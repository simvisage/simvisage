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

    import pylab as p

    test_files = ['BTT-4c-2cm-TU-0-V14_MxN2',
                  'BTT-4c-2cm-TU-0-V13_MxN2']

    test_file = os.path.join(simdb.exdata_dir,
                             'bending_tensile_test',
                             '2014-06-12_BTT-4c-2cm-0-TU_MxN2',
                             'BTT-4c-2cm-TU-0-V01_MxN2.DAT')

    e1 = ExpBTTDB(data_file=test_file)
    e1.process_source_data()
    print 'F_max1', e1.F_max1

    print '4' * 20
    p.plot(e1.t_aramis, e1.w_t_aramis, color='blue')
    p.plot(e1.t, e1.w, color='red')
    p.show()

    aramis_file_path = e1.get_cached_aramis_file('Xf15s3-Yf15s3')

    # @todo:
    # 1) define the mapping between astage_idx(t).
    # 2) define a method evaluating strain profile within the object.

    AI = AramisInfo(data_dir=aramis_file_path)
    AD = AramisData(aramis_info=AI, evaluated_step_idx=60)

    print AI.number_of_steps

    AC = AramisBSA(aramis_info=AI,
                   aramis_data=AD,
                   integ_radius=10)

    print AC.d_ux_arr2.shape
    mid_idx = AC.d_ux_arr2.shape[1] / 2
    # eps = np.mean(AC.d_ux_arr2[:, mid_idx - 0:mid_idx + 1], axis=1)
    eps = AC.d_ux_arr2[:, mid_idx]
    y = AD.y_arr_undeformed[:, mid_idx]
    p.plot(eps, y)
#     p.plot([-0.004, 0.014], [0, 0], color='black')
#     p.plot([0, 0], [-10, 10], color='black')
    p.show()
    # AC.configure_traits()
