'''
Created on Apr 25, 2012

@author: rch
'''

import numpy as np
import numpy.ma as ma
import etsproxy.mayavi.mlab as m

from matresdev.db.simdb import \
    SimDB

import os.path

# Access to the toplevel directory of the database
#
simdb = SimDB()

file_name = os.path.join(simdb.exdata_dir, 'tensile_tests',
                       'doge_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4',
                       'ARAMIS', 'V1', 'TT-V1-Stufe-0-428.txt')

file_name = os.path.join(simdb.exdata_dir, 'tensile_tests',
                       'doge_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4',
                       'ARAMIS', 'V1_kurz', 'P1-s458-Xf15a1-Yf5a4-Ausschnitt-Stufe-0-430.txt')

input_arr = np.loadtxt(file_name,
                    skiprows = 14,
                    usecols = [0, 1, 2, 3, 4, 8, 9, 10]
    )

print 'data', input_arr.shape

print 'first row', input_arr[0, :]
# identify the points that do not belong to the specimen
# 
select_idx = np.where(input_arr[:, 4] < -50.0)
print 'select_idx', select_idx
data_arr = input_arr[ select_idx[0], :]
data_arr = input_arr

x_arr, y_arr, z_arr = data_arr[:, [2, 3, 4]].T

L_x = np.max(x_arr) - np.min(x_arr)
L_y = np.max(y_arr) - np.min(y_arr)

print 'L_x', L_x
print 'L_y', L_y

x_idx = np.array(data_arr[:, 0], dtype = int)
y_idx = np.array(data_arr[:, 1], dtype = int)

print 'x_idx', x_idx

# find the shape of the point grid
# by taking aramis indices stored in
# first two columns
x_min, x_max = np.min(x_idx), np.max(x_idx)
y_min, y_max = np.min(y_idx), np.max(y_idx)

n_x = x_max - x_min
n_y = y_max - y_min

print 'n_x', n_x
print 'n_y', n_y

# total number of data values associated with a single point
n_data = data_arr.shape[1] - 2

# define a three dimensional array
# first two indices are grid indices, third
# index selects one of the values associated with a point
daf = np.zeros((n_x + 1, n_y + 1, n_data), dtype = float)

# define an array of indexes performing the mapping
# from the data_arr into daf
select_idx_arr = (x_idx - x_min, y_idx - y_min, slice(None))
daf[select_idx_arr] = data_arr[:, 2:] 

ux_arr = daf[:, :, 3]
ux_mask = ux_arr == 0.0
ux_masked = ma.masked_array(ux_arr, mask = ux_mask)

ux_avg = ma.average(ux_masked, axis = 1)

x_idx_zeros, y_idx_zeros = np.where(ux_mask)

print 'ux_avg', ux_avg[x_idx_zeros].shape
print 'ux_arr', ux_arr[x_idx_zeros, y_idx_zeros].shape

ux_arr[x_idx_zeros, y_idx_zeros] = ux_avg[x_idx_zeros]

print 'x_idx', x_idx_zeros
print 'y_idx', y_idx_zeros

ux_min, ux_max = np.min(ux_arr), np.max(ux_arr)
dux = ux_max - ux_min
print '-----------------------------------'
print 'total displacement', dux
print 'middle crack width', dux / 17.0
print '-----------------------------------'

# construct the mask for elements to be ignored
mask_idx_array = np.zeros((n_x + 1, n_y + 1), dtype = bool)
mask_idx_array[:, :] = True
mask_idx_array[(x_idx - x_min, y_idx - y_min)] = False

# generate the elem_node_map

integ_range = 22
d_x = L_x / float(n_x)
eps1_x_arr = (ux_arr[1:, :] - ux_arr[:-1, :])
eps4_x_arr = (ux_arr[integ_range:, :] - ux_arr[:-integ_range, :])
eps_x_avg = np.average(eps1_x_arr, axis = 1) 
eps4_x_arr[ eps4_x_arr < 0.0 ] = 0.0
eps4_x_avg = np.average(eps4_x_arr, axis = 1)


import pylab as p

row_idx = (10, n_y - int(n_y / 2), n_y - 10)
plot_row_idx = 10
p.figure()
#p.plot(daf[:-1, plot_row_idx, 0], eps1_x_arr[:, row_idx], color = 'blue')
p.plot(daf[:-integ_range, plot_row_idx, 0], eps4_x_arr[:, row_idx], color = 'green')
p.plot(daf[:-integ_range, plot_row_idx, 0], eps4_x_avg, color = 'black')
#p.plot(daf[:-1, plot_row_idx, 0], eps_x_avg, color = 'red')

p.figure()
p.plot(daf[:, plot_row_idx, 0], ux_arr[:, :], color = 'green')
p.plot(daf[:, plot_row_idx, 0], ux_avg, color = 'red'
           )
p.show()

#qargs = [ daf[:, :, i, np.newaxis] for i in range(0, 6)]

#m.quiver3d(*qargs)

s = m.points3d(daf[:-integ_range, :, 0],
               daf[:-integ_range, :, 1],
               daf[:-integ_range, :, 2],
               eps4_x_arr)

## print a surface - why is moved in the y direction?
#m.surf(daf[:-integ_range, :, 0], daf[:-integ_range, :, 1], eps4_x_arr)#, mask = mask_idx_array)

m.show()

# @todo: 
# sensitivity with respect to the integ_range
# 
