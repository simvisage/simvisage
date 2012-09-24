'''
Created on Sep 11, 2012

The script demonstrates how to fetch a database object with
measured and postprocessed test results and how to 
access the data.

@author: rch
'''

from os.path import join
from matresdev.db import SimDB
from matresdev.db.exdb import ExRun

if __name__ == '__main__':    

    simdb = SimDB()

    # specify the path to the data file.
    test_file = join(simdb.exdata_dir,
                     'tensile_tests',
                     'dog_bone',
                     '2012-04-12_TT-12c-6cm-0-TU_SH4',
                     'TT-12c-6cm-0-TU-SH4-V1.DAT')

    from quaducom.devproc.tensile_test.dog_bone.exp_tt_db import ExpTTDB
    
    print ExpTTDB.db.instances.keys()
    print ExpTTDB.db['TT-12c-6cm-0-TU-SH1-V3']
    
    # construct the experiment
    tt = ExRun(data_file = test_file)

    # access the response values.
    print tt.ex_type.eps_asc
    print tt.ex_type.sig_c_asc
