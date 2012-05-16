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
import pylab as p


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

# find the shape of the point grid
# by taking aramis indices stored in
# first two columns
x_min, x_max = np.min(x_idx), np.max(x_idx)
y_min, y_max = np.min(y_idx), np.max(y_idx)

n_x = x_max - x_min
n_y = y_max - y_min

print 'n_x', n_x
print 'n_y', n_y

# construct the mask for elements to be ignored
mask_idx_array = np.zeros((n_x + 1, n_y + 1), dtype = bool)
mask_idx_array[:, :] = True
mask_idx_array[(x_idx - x_min, y_idx - y_min)] = False


# total number of data values associated with a single point
n_values = data_arr.shape[1] - 2

# define a three dimensional array
# first two indices are grid indices, third
# index selects one of the values associated with a point
daf = np.zeros((n_x + 1, n_y + 1, n_values), dtype = float)

# define an array of indexes performing the mapping
# from the data_arr into daf
select_idx_arr = (x_idx - x_min, y_idx - y_min, slice(None))
daf[select_idx_arr] = data_arr[:, 2:] 

ux_arr = daf[:, :, 3]
print 'ux_arr', ux_arr[:, 0]
print 'mask_idx_array', mask_idx_array
ux_masked = ma.masked_array(ux_arr, mask = mask_idx_array)

ux_avg = ma.average(ux_masked, axis = 1)

x_idx_zeros, y_idx_zeros = np.where(mask_idx_array)

print x_idx_zeros, y_idx_zeros

ux_arr[x_idx_zeros, y_idx_zeros] = ux_avg[x_idx_zeros]

ux_min, ux_max = np.min(ux_arr), np.max(ux_arr)
dux = ux_max - ux_min
print '-----------------------------------'
print 'total displacement', dux
print 'miDle crack width', dux / 17.0
print '-----------------------------------'

# generate the elem_node_map

integ_radius = 11
integ_range = 2 * integ_radius

d_x = L_x / float(n_x)

d4_ux_arr = np.zeros_like(ux_arr)
d4_ux_arr[integ_radius:-integ_radius, :] = (ux_arr[integ_range:, :] - 
                                            ux_arr[:-integ_range, :])

# cutoff the negative strains - noise
d4_ux_arr[ d4_ux_arr < 0.0 ] = 0.0
d4_ux_avg = np.average(d4_ux_arr, axis = 1)

dd44_ux_arr = np.zeros_like(ux_arr)
dd44_ux_arr[integ_radius:-integ_radius, :] = (d4_ux_arr[integ_range:, :] - 
                                                      d4_ux_arr[:-integ_range, :])
dd44_ux_avg = np.average(dd44_ux_arr, axis = 1)

# third difference field
ddd444_ux_arr = np.zeros_like(ux_arr)
ddd444_ux_arr[integ_radius:-integ_radius, :] = (dd44_ux_arr[integ_range:, :] - 
                                                        dd44_ux_arr[:-integ_range, :])
ddd444_ux_avg = np.average(ddd444_ux_arr, axis = 1)

#crack_idx = ((ddd444_ux_avg[1:] * ddd444_ux_avg[:-1] < 0.0) * 
#             (ddd444_ux_avg[1:] < ddd444_ux_avg[:-1])) 

crack_idx_arr = ((dd44_ux_arr[1:, :] * dd44_ux_arr[:-1, :] < 0.0) * 
                 ((ddd444_ux_arr[1:, :] + ddd444_ux_arr[:-1, :]) / 2.0 < -0.005)) 

crack_idx_avg = ((dd44_ux_avg[1:] * dd44_ux_avg[:-1] < 0.0) * 
             ((ddd444_ux_avg[1:] + ddd444_ux_avg[:-1]) / 2.0 < -0.005)) 

w_arr = d4_ux_arr[np.where(crack_idx_arr)]
w_avg = np.max(d4_ux_arr[np.where(crack_idx_avg)], axis = 1)

print 'crack_openings', w_avg

crack_avg = np.array(crack_idx_avg, dtype = float) * 0.2

row_idx = (10, n_y - int(n_y / 2), n_y - 10)
plot_row_idx = 10

x_arr = daf[:, plot_row_idx, 0]
x4_arr = daf[integ_radius:-integ_radius, plot_row_idx, 0]
x44_arr = daf[2 * integ_radius:-2 * integ_radius, plot_row_idx, 0]
x444_arr = daf[3 * integ_radius:-3 * integ_radius, plot_row_idx, 0]
xcrack = x444_arr[:-1] + d_x / 2.0

idx_arr = np.arange(ux_arr.shape[0], dtype = int)

p.figure()
p.subplot(2, 2, 1)
#p.plot(idx_arr, d4_ux_arr[:, row_idx], color = 'green')
p.plot(idx_arr, d4_ux_avg, color = 'black')
p.plot(idx_arr[:-1], crack_avg, color = 'magenta', linewidth = 4)

p.subplot(2, 2, 3)
#p.plot(idx_arr, dd44_ux_arr[:, row_idx], color = 'green')
p.plot(idx_arr, dd44_ux_avg, color = 'black')

#p.plot(idx_arr, ddd444_ux_arr[:, row_idx], color = 'red')
p.plot(idx_arr, ddd444_ux_avg, color = 'blue')
p.plot(idx_arr[:-1], crack_avg, color = 'magenta', linewidth = 4)

p.subplot(2, 2, 2)
p.plot(idx_arr, ux_arr[:, :], color = 'green')
p.plot(idx_arr, ux_avg, color = 'red')

#p.figure()
#p.hist(w_avg)

p.subplot(2, 2, 4)
p.hist(w_arr, bins = 40)

p.show()

#qargs = [ daf[:, :, i, np.newaxis] for i in range(0, 6)]

#m.quiver3d(*qargs)

#s = m.points3d(daf[:-integ_range, :, 0],
#               daf[:-integ_range, :, 1],
#               daf[:-integ_range, :, 2],
#               d4_ux_arr)

# print a surface - why is moved in the y direction?
m.surf(daf[:-integ_range, :, 0], daf[:-integ_range, :, 1], d4_ux_arr)#, mask = mask_idx_array)

m.show()

# @todo: 
# sensitivity with respect to the integ_range
# 
