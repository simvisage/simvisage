'''
Created on Apr 25, 2012

@author: rch
'''

import numpy as np
import numpy.ma as ma
import enthought.mayavi.mlab as m
#
from matresdev.db.simdb import \
   SimDB
#
import os.path
#
## Access to the toplevel directory of the database
##
simdb = SimDB()

#file_name = os.path.join(simdb.exdata_dir, 'tensile_tests',
#                       'doge_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4',
#                       'ARAMIS', 'V1', 'TT-V1-Stufe-0-428.txt')


home_dir = os.path.expanduser('~')
print 'home_dir', home_dir

file_name = os.path.join(home_dir, 'uni', '6-Semester','Praktikum', 'Praktikum_Massivbau','ARAMIS','V1','TT-V1-Stufe-0-300.txt')
print 'file_name', file_name

#file_name = os.path.join('C:\\','Praktikum_Massivbau','ARAMIS','V1','TT-V1TT-V1-Stufe-0-300.txt')
file_name = 'C:\\Praktikum_Massivbau\ARAMIS\V1\TT-V1-Stufe-0-300.txt'
input_arr = np.loadtxt(file_name,
                    skiprows = 14,
                    usecols = [0, 1, 5, 6, 7, 8, 9, 10]
    )

print 'data', input_arr.shape

print input_arr[0, :]
select_idx = np.where(input_arr[:, 4] < -50.0)[0]
print 'select_idx', select_idx

data_arr = input_arr[ select_idx, :]
print 'data_arr', data_arr

x_idx = np.array(data_arr[:, 0], dtype = int)
y_idx = np.array(data_arr[:, 1], dtype = int)

print 'x_idx', x_idx

x_min, x_max = np.min(x_idx), np.max(x_idx)
y_min, y_max = np.min(y_idx), np.max(y_idx)

n_x = x_max - x_min
n_y = y_max - y_min

print 'n_x', n_x
print 'n_y', n_y

n_data = data_arr.shape[1] - 2

print 'n_data', n_data


daf = np.zeros((n_x + 1, n_y + 1, n_data), dtype = float)

select_idx_arr = (x_idx - x_min, y_idx - y_min, slice(None))
daf[select_idx_arr] = data_arr[:, 2:] 

#construct the mask for elements to be ignored
mask_idx_array = np.zeros((n_x + 1, n_y + 1), dtype = bool)
mask_idx_array[:, :] = True
mask_idx_array[(x_idx - x_min, y_idx - y_min)] = False

# generate the elem_node_map

elem_arr = np.array([daf[:-1, :-1], daf[1:, :-1], daf[1:, 1:], daf[ :-1, 1:]])
print elem_arr[:, 0, 0, :]
#
## print a surface - why is moved in the y direction?
#m.surf(daf[:, :, 0], daf[:, :, 1], daf[:, :, 3], mask = mask_idx_array)
m.surf(daf[:, :, 0], daf[:, :, 1], daf[:, :, 2], mask = mask_idx_array)

qargs = [ daf[:, :, i, np.newaxis] for i in range(0, 6)]

#m.quiver3d(*qargs)

s = m.points3d(daf[:, :, 0],             daf[:, :, 1],
              daf[:, :, 2],
              daf[:, :, 3])

m.show()

