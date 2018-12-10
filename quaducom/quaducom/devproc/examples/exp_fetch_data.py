'''
Created on Sep 11, 2012

The script demonstrates how to fetch a database object with
measured and postprocessed test results and how to 
access the data.

@author: rch
'''
#from traits.etsconfig.api import ETSConfig
#ETSConfig.toolkit = 'wx'
from matresdev.db import SimDB

if __name__ == '__main__':

    simdb = SimDB()

    # specify the path to the data file.

    from quaducom.devproc.tensile_test.dog_bone.exp_tt_db import ExpTTDB

    print(list(ExpTTDB.db.instances.keys()))
    print(ExpTTDB.db['TT-12c-6cm-0-TU-SH1-V3'])

    ExpTTDB.db['TT-12c-6cm-0-TU-SH1-V3'].configure_traits()
