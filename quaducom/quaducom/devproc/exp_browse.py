'''
Created on Mar 15, 2012

@author: rch
'''

import os.path

from matresdev.db.simdb import \
    SimDB

simdb = SimDB()


from matresdev.db.exdb.ex_run_view import ExRunView

# ex_path = os.path.join(simdb.exdata_dir, 'slab_tests', '2013-07-10_ST-6c-2cm-TU_bs2', 'ST-6c-2cm-TU_bs2.DAT')

ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-03-20_TT-12c-4cm-0-TU_SH3',
                         'TT-12c-4cm-TU-0-SH3-V2.DAT')

ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2',
                         'TT-12c-6cm-0-TU-SH2-V1.DAT')

# ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d',
#                        'TTb-6c-2cm-0-TU-V2_bs4.DAT')
#
# ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-02_BT-6c-2cm-0-TU_bs4',
#                                       'BT-6c-2cm-0-TU-V1_bs4.DAT')

ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-02_BT-6c-2cm-0-TU_bs4',
                                       'BT-6c-2cm-0-TU-V1_bs4.DAT')
ex_path = os.path.join(simdb.exdata_dir, 'bending_tensile_test', '2014-06-12_BTT-4c-2cm-0-TU_MxN2',
                                       'BTT-4c-2cm-TU-0-V01_MxN2.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE',
#                                       'BT-3PT-12c-6cm-0-Tu-V1.raw')

# ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU',
#                                       'BT-4PT-12c-6cm-SH4', 'BT-4PT-12c-6cm-SH4-V1.DAT')
#
# ex_path = os.path.join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU',
#                                       'ST-12c-6cm-u-TU.DAT')

print os.path.exists(ex_path)

doe_reader = ExRunView(data_file=ex_path)

# print 'model', doe_reader.model.ex_type.Kraft
# print 'model', doe_reader.model.ex_type.WA_VL
# print 'model', doe_reader.model.ex_type.WA_VR
# print 'model', doe_reader.model.ex_type.WA_HL
# print 'model', doe_reader.model.ex_type.WA_HR

doe_reader.configure_traits()
