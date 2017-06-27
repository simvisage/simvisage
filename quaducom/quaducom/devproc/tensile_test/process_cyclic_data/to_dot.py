'''
Created on Aug 5, 2015

@author: Yingxiong
'''
from matresdev.db.simdb import SimDB
simdb = SimDB()
import os

test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-08-03_TTb-2C-14mm-0-3300SBR_cyc-Aramis2d'
                              )

# file_ = os.path.join(test_file_path, 'TTb-2C-14mm-0-3300SBR-V1_cyc-Aramis2d')
# file_ = os.path.join(test_file_path, 'TTb-2C-14mm-0-3300SBR-V6_cyc-Aramis2d')
file_ = os.path.join(test_file_path, 'TTb-2C-14mm-0-3300SBR-V4_cyc-Aramis2d')

with open(file_ + '.csv', 'r') as fin:
    with open(file_ + '-dot.csv', 'w') as fout:
        for line in fin:
            fout.write(line.replace(',', '.'))
