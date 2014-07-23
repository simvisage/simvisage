'''
Created on Jul 1, 2014

Script evaluating the mxn interaction between
tensile and bending load
'''

if __name__ == '__main__':

    from exp_btt_db import ExpBTTDB
    from matresdev.db.simdb import SimDB
    simdb = SimDB()
    import os

    test_file = os.path.join(simdb.exdata_dir,
                           'bending_tensile_test',
                           '2014-06-12_BTT-6c-2cm-0-TU_MxN2',
                           'BTT-6c-2cm-TU-0-V14_MxN2.DAT')

    e1 = ExpBTTDB(data_file=test_file)
    e1.process_source_data()
    print 'F_max1', e1.F_max1


    e1.extended_data_cfg

    aramis_data = e1.cache_data('19s15')

    # call the processor

    from aramis_cdt.aramis_bsa import AramisBSA

    aramis_bsa = AramisBSA(aramis_data=aramis_data)

    aramis_bsa.d_ux_arr2

    # ExpBTTDB.db.configure_traits()
