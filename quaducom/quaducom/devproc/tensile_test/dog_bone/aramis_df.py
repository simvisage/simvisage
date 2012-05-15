'''
Created on Apr 25, 2012

@author: rch
'''

import numpy as np
import etsproxy.mayavi.mlab as m

from matresdev.db.simdb import \
    SimDB

import os.path

# Access to the toplevel directory of the database
#
simdb = SimDB()

file_name = os.path.join(simdb.exdata_dir, 'tensile_tests',
                       'dog_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4',
                       'ARAMIS', 'V1', 'TT-V1-Stufe-0-428.txt')

input_arr = np.loadtxt(file_name,
                    skiprows = 14,
                    # idx_x, idx_y, (x, y, z)_undef, (x, y, z)_def, (x, y, z)_displ 
                    usecols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    )

print 'data', input_arr.shape
print input_arr[0, :]

# @TODO: check in which txt-files are the coordinates of the origin (camera position) included? 
#select_idx = np.where(input_arr[:, 4] < -50.0)[0]
#print 'select_idx.shape', select_idx.shape
#data_arr = input_arr[ select_idx, :]

data_arr = input_arr

x_idx = np.array(data_arr[:, 0], dtype = int)
y_idx = np.array(data_arr[:, 1], dtype = int)

print 'x_idx', x_idx
print 'y_idx', y_idx

x_min, x_max = np.min(x_idx), np.max(x_idx)
print 'x_min, x_max', x_min, x_max

y_min, y_max = np.min(y_idx), np.max(y_idx)
print 'y_min, y_max', y_min, y_max

n_x = x_max - x_min
n_y = y_max - y_min

print 'n_x', n_x
print 'n_y', n_y

# n_data = 6 = (x,y,z)undeformed + (x,y,z)deformed
n_data = data_arr.shape[1] - 2

# displacement field

# generate a rectangular displacement field 
# (fill missing values of the grid with zeros)
daf = np.zeros((n_x + 1, n_y + 1, n_data), dtype = float)

select_idx_arr = (x_idx - x_min, y_idx - y_min, slice(None))
#select_idx_arr = (x_idx - x_min, y_idx - y_min, Ellipsis)
print 'select_idx_arr',select_idx_arr

daf[select_idx_arr] = data_arr[:, 2:] 
print 'daf', daf

# construct the mask for elements to be ignored
mask_idx_array = np.zeros((n_x + 1, n_y + 1), dtype = bool)
mask_idx_array[:, :] = True
mask_idx_array[(x_idx - x_min, y_idx - y_min)] = False

# generate the elem_node_map

#elem_arr = np.array([daf[:-1, :-1], daf[1:, :-1], daf[1:, 1:], daf[ :-1, 1:]])
#print elem_arr[:, 0, 0, :]

# print a surface - why is moved in the y direction?
#m.surf(daf[:, :, 0], daf[:, :, 1], daf[:, :, 3], mask = mask_idx_array)
#m.surf(daf[:, :, 0], daf[:, :, 1], daf[:, :, 2], mask = mask_idx_array)

#qargs = [ daf[:, :, i, np.newaxis] for i in range(0, 6)]

#m.quiver3d(*qargs)

m.figure( figure = "ARAMIS plot",
                 bgcolor = ( 1.0, 1.0, 1.0 ),
                 fgcolor = ( 0.0, 0.0, 0.0 ) )

Ux = daf[:, :, 8]
print 'Ux', Ux

m.points3d(daf[:, :, 0],
           daf[:, :, 1],
           daf[:, :, 2],
           Ux,
#           colormap = "YlOrBr",
           mode = "cube",
           scale_factor = 0.5 )

m.axes()
m.scalarbar( title = 'UX [m]', orientation = 'vertical' )

m.show()
