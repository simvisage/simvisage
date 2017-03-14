'''
Created on Aug 5, 2015

@author: Yingxiong
'''
import numpy as np
from matresdev.db.simdb import SimDB
simdb = SimDB()
import os

test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-08-03_TTb-2C-14mm-0-3300SBR_cyc-Aramis2d'
                              )

# file_ = os.path.join(test_file_path, 'TTb-2C-14mm-0-3300SBR-V1_cyc-Aramis2d-dot')
# file_ = os.path.join(test_file_path, 'TTb-2C-14mm-0-3300SBR-V6_cyc-Aramis2d-dot')
file_ = os.path.join(test_file_path, 'TTb-2C-14mm-0-3300SBR-V4_cyc-Aramis2d-dot')

# force
with open(file_ + '.csv', 'r') as fin:
    with open(file_[:-4] + '-force.csv', 'w') as fout:
        i = 0
        for line in fin:
            if i < 2:
                i += 1
            else:
                fout.write(line.split(';')[1] + '\n')

# displacement

# with open(file_ + '.csv', 'r') as fin:
#    with open(file_ + '-disp.csv', 'w') as fout:
#        i = 0
#        for line in fin:
#            if i < 2:
#                i += 1
#            else:
#                disp = -(
#                    float(line.split(';')[5]) + float(line.split(';')[6])) / 2.
#                fout.write(str(disp) + '\n')


# with open(file_ + '.csv', 'r') as fin:
#
# with open(file_ + '-disp.csv', 'w') as fout:
#     i = 0
#     for line in fin:
#         if i < 2:
#             i += 1
#         elif np.abs(float(line.split(';')[5]) - float(line.split(';')[6])) > np.abs((
#                 float(line.split(';')[5]) + float(line.split(';')[6])) / 2.) * 0.5:
#             print line.split(';')[5], line.split(';')[6]
