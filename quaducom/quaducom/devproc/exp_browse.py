'''
Created on Mar 15, 2012

@author: rch
'''

import os.path

from matresdev.db.simdb import \
    SimDB

simdb = SimDB()

from matresdev.db.exdb.ex_run_view import ExRunView

ex_path = os.path.join(simdb.exdata_dir, 'tensile_tests', 'ZiE_2011-06-08_TT-12c-6cm-90-TU',
                        'TT-12c-6cm-90-TU-V1.DAT')

#    ex_path = os.path.join( simdb.exdata_dir, 'plate_tests', 'PT-6a-ibac',
#                            'PTi-6a-woSF', 'PTi-6a-woSF-V1.DAT' )

#    ex_path = os.path.join( simdb.exdata_dir, 'plate_tests', 'PT-10a',
#                            'PT10-10a.DAT' )

doe_reader = ExRunView(data_file = ex_path)
doe_reader.configure_traits()
