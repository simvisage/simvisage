from numpy import *
from promod.simdb import \
    SimDB

import os
import pickle
import string
from os.path import join

# Access to the top level directory of the database
#
simdb = SimDB()

data_dir = os.path.join( simdb.simdb_dir,
                         'simdata', 'input_data_NM_stb',
                         '1shell' )

file_name = os.path.join( data_dir, 'LC1.csv' )

print '*** read state data from file: %s ***' % ( file_name )

arr = loadtxt( file_name, delimiter = ";" )
state_data = arr
print state_data
print state_data.shape

