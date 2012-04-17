'''
Created on Mar 6, 2012

@author: rch
'''

import fnmatch
import os.path

def replace_string_in_files(path,
                            replace_dict,
                            file_ext = '.pickle'
                            ):
    '''Walk through all files within the subdirectories of `path'
    and perform the replacement of strings specifified in `replace_dict\
    in all files ending with `file_ext'
    '''
    selected_files = []
    for path, dirs, files in os.walk(path):
        for f in fnmatch.filter(files, '*' + file_ext):
            selected_files.append(os.path.join(path, f))

    for file_name in selected_files:
        print 'replacement applied to', file_name
        fid = open(file_name, "r")
        text = fid.read()
        fid.close()
        # write the replacements back into the file.
        fid = open(file_name, "w")
        for old, new in replace_dict.items():
            text = text.replace(old, new)
        fid.write(text)
        fid.close()
        print 'finished'

def migrate_classes(migration_table):
    '''Walk through all files within the simdb storage'
    and perform the replacement of strings specifified in `replace_dict\
    in all files ending with pickle'
    '''
    from matresdev.db.simdb import SimDB
    simdb = SimDB()    

    replace_string_in_files(simdb.exdata_dir,
                            migration_table, '.pickle')
    replace_string_in_files(simdb.matdata_dir,
                            migration_table, '.pickle')
    replace_string_in_files(simdb.exdata_dir,
                            migration_table, 'ex_type.cls')
                    
if __name__ == '__main__':

    from matresdev.db.simdb import SimDB
    simdb = SimDB()

    migration_table = {'promod.matdb.trc' : 'matresdev.db.matdb.trc',
                       'promod.exdb.ex_composite_tensile_test' : 'quaducom.devproc.tt.dbtt.exp_dbtt',
                       'ExCompositeTensileTest' : 'ExpDogBoneTensileTest',
                       'promod.exdb.ex_bending_test' : 'quaducom.devproc.bt.p3.exp_bt_3pt',
                       'ExBendingTest' : 'ExpBendingTest3Pt',
                       'promod.exdb.ex_plate_test' : 'quaducom.devproc.st.exp_st',
                       'ExPlateTest' : 'ExpSlabTest'} 


    migrate_classes(migration_table)
