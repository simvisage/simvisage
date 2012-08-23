'''
Created on Mar 15, 2012

@author: rch
'''

import os.path

from matresdev.db.simdb import \
    SimDB

simdb = SimDB()

from matresdev.db.exdb.ex_run_view import ExRunView

#ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'ZiE_2011-06-08_TT-12c-6cm-90-TU',
#                        'TT-12c-6cm-90-TU-V1.DAT')

#ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', '2011-06-10_TT-12c-6cm-90-TU_ZiE',
#                        'TT-12c-6cm-90-TU-V1.DAT')

ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU',
                                        'BT-4PT-12c-6cm-SH4', 'BT-4PT-12c-6cm-SH4-V1.DAT')
#ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU',
#                                       'BT-4PT-12c-6cm-SH4', 'BT-4PT-12c-6cm-SH4-V1.DAT')



print os.path.exists(ex_path)

doe_reader = ExRunView(data_file = ex_path)

#print 'model', doe_reader.model.ex_type.Kraft
#print 'model', doe_reader.model.ex_type.WA_VL 
#print 'model', doe_reader.model.ex_type.WA_VR
#print 'model', doe_reader.model.ex_type.WA_HL
#print 'model', doe_reader.model.ex_type.WA_HR

doe_reader.configure_traits()
