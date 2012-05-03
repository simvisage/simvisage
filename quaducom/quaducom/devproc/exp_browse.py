'''
Created on Mar 15, 2012

@author: rch
'''

import os.path

from matresdev.db.simdb import \
    SimDB

simdb = SimDB()

from matresdev.db.exdb.ex_run_view import ExRunView

ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-12-07_TT-3g-1cm-a-FR_ibac', 'with_short_fibres', 'TT-3g-2cm-0-90-0-FR-wSF-ibac-V1.DAT')
#ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU',
#                                         'BT-4PT-12c-6cm-SH4', 'BT-4PT-12c-6cm-SH4-V1.DAT')


doe_reader = ExRunView(data_file = ex_path)
doe_reader.configure_traits()
