'''
Created on Mar 15, 2012

@author: rch
'''

import os.path

from matresdev.db.simdb import \
    SimDB

simdb = SimDB()

from matresdev.db.exdb.ex_run_view import ExRunView

#ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-12-07_TT-3g-1cm-a-FR_ibac', 'with_short_fibres', 'TT-3g-2cm-0-90-0-FR-wSF-ibac-V1.DAT')
#ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4', 'TT-12c-6cm-0-TU-SH4-V1.DAT')
ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V1.DAT')
#ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU',
#                                         'BT-4PT-12c-6cm-SH4', 'BT-4PT-12c-6cm-SH4-V1.DAT')


doe_reader = ExRunView(data_file = ex_path)

#print 'model', doe_reader.model.ex_type.Kraft
#print 'model', doe_reader.model.ex_type.WA_VL 
#print 'model', doe_reader.model.ex_type.WA_VR
#print 'model', doe_reader.model.ex_type.WA_HL
#print 'model', doe_reader.model.ex_type.WA_HR

doe_reader.configure_traits()
