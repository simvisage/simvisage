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

ex_path = os.path.join(simdb.exdata_dir, 'bending_tests', '2012-04-03_BT-4PT-12c-6cm-0-TU',
                                         'BT-4PT-12c-6cm-SH4', 'BT-4PT-12c-6cm-SH4-V1.DAT')


doe_reader = ExRunView(data_file = ex_path)
doe_reader.configure_traits()
