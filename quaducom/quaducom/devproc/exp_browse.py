'''
Created on Mar 15, 2012

@author: rch
'''

import os.path
from matresdev.db.exdb.ex_run_view import ExRunView
from matresdev.db.simdb.simdb import simdb


# ex_path = os.path.join(simdb.exdata_dir, 'slab_tests', '2013-07-10_ST-6c-2cm-TU_bs2', 'ST-6c-2cm-TU_bs2.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'single_pullout_tests',  '2017-08-03_SPO-1C-40mm-90-3,62EP',
#                        'SPOS19-1.DAT')
ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2',
                       'TT-12c-6cm-0-TU-SH2-V1.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone',
#                       '2015-02-24_TT-12c-6cm-0-TU_QS-2Jahre', 'TT-12c-6cm-0-TU-SH3F-V1_QS-2Jahre.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'slab_tests',
#                        '2015-08-27_ST-2C-6cm-3300EP_SPP', 'B-C1-3-3.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2015-02-24_TT-12c-6cm-0-TU_QS-2Jahre', 'TT-12c-6cm-0-TU-SH3F-V1_QS-2Jahre.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d',
#                        'TTb-6c-2cm-0-TU-V2_bs4.DAT')
#
# ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-02_BT-6c-2cm-0-TU_bs4',
#                                       'BT-6c-2cm-0-TU-V1_bs4.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'four_point', '2015-09-02_BT-1C-55mm-0-3300SBR_cyc-Aramis2d',
# ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-02_BT-6c-2cm-0-TU_bs4',
#                       'BT-6c-2cm-0-TU-V1_bs4.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'bending_tensile_test', '2014-06-12_BTT-4c-2cm-0-TU_MxN2',
#                       'BTT-4c-2cm-TU-0-V01_MxN2.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE',
#                                       'BT-3PT-12c-6cm-0-Tu-V1.raw')

# ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU',
#                                       'BT-4PT-12c-6cm-SH4', 'BT-4PT-12c-6cm-SH4-V1.DAT')
#
# ex_path = os.path.join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU',
#                                       'ST-12c-6cm-u-TU.DAT')

# ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2014-04-30_TTb-6c-2cm-0-TU_NxM1',
#                                        'TTb-6c-2cm-0-TU-V16_NxM1.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests',
#                       'buttstrap_clamping',
#                       '2016-03-18_TTb-2C-9mm-0-3300SBR_R4_Dresden',
#                       'B3-Ia-DK3-A-1.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests',
#                        'buttstrap_clamping',
#                        '2016-03-18_TTb-2C-9mm-0-3300SBR_R4_Dresden',
#                        'B3-Ia-DK3-A-1.DAT')
# ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests',
#                       'buttstrap_clamping',
#                       '2016-03-22_TTb_TUWien',
#                       'P10_1.DAT')
#
print os.path.exists(ex_path)

doe_reader = ExRunView(data_file=ex_path)


# print 'model', doe_reader.model.ex_type.Kraft
# print 'model', doe_reader.model.ex_type.WA_VR
# print 'model', doe_reader.model.ex_type.WA_HL
# print 'model', doe_reader.model.ex_type.WA_HR

doe_reader.configure_traits()
