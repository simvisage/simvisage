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
from matresdev.db.simdb import \
   SimDB
#
import os.path
#
## Access to the toplevel directory of the database
##
simdb = SimDB()

file_name = os.path.join(simdb.exdata_dir, 'tensile_tests',
                       'dog_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4',
                       'ARAMIS', 'V1', 'TT-V1-Stufe-0-428.txt')


home_dir = os.path.expanduser('~')
print 'home_dir', home_dir
file_name = os.path.join(home_dir, 'Trash', 'V3', 'TT-V3-Stufe-0-82.txt')

#file_name = os.path.join(home_dir, 'uni', '6-Semester', 'Praktikum', 'Praktikum_Massivbau', 'ARAMIS', 'V1', 'TT-V1-Stufe-0-300.txt')
print 'file_name', file_name

#file_name = os.path.join('C:\\', 'Praktikum_Massivbau', 'ARAMIS', 'V1', 'TT-V1TT-V1-Stufe-0-300.txt')
#file_name = 'C:\\Praktikum_Massivbau\ARAMIS\V1\TT-V1-Stufe-0-300.txt'
input_arr = np.loadtxt(file_name,
                    skiprows=14,
                    usecols=[0, 1, 2, 3, 4, 8, 9, 10]
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
X_arr = x_arr

L_x = np.max(x_arr) - np.min(x_arr)
L_y = np.max(y_arr) - np.min(y_arr)

print 'L_x', L_x
print 'L_y', L_y

x_idx = np.array(data_arr[:, 0], dtype=int)
y_idx = np.array(data_arr[:, 1], dtype=int)

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
mask_idx_array = np.zeros((n_x + 1, n_y + 1), dtype=bool)
mask_idx_array[:, :] = True
mask_idx_array[(x_idx - x_min, y_idx - y_min)] = False


# total number of data values associated with a single point
n_values = data_arr.shape[1] - 2

# define a three dimensional array
# first two indices are grid indices, third
# index selects one of the values associated with a point
daf = np.zeros((n_x + 1, n_y + 1, n_values), dtype=float)

# define an array of indexes performing the mapping
# from the data_arr into daf
select_idx_arr = (x_idx - x_min, y_idx - y_min, slice(None))
daf[select_idx_arr] = data_arr[:, 2:]



# Coordinate transformation:
#
#x_vec=np.average(np.array([daf[-1,:,0]-daf[0,:,0],[daf[-1,:,1]-daf[0,:,1],[daf[-1,:,2]-daf[0,:,2]]))
#y_vec=np.average(np.array([daf[:,-1,0]-daf[:,0,0],[daf[:,-1,1]-daf[:,0,1],[daf[0,-1,2]-daf[:,0,2]]))
daf_new = np.copy(daf)

# @todo: select a facette with none-zero coordinate values
# e.g. select a facette from the last middle axis and in a list of the first 10 or last 10 facettes
# 
x_vec_ = daf[-1, 10, :3] - daf[0, 10, :3]
#print'x_vec_', x_vec_

x_vec_normed = x_vec_ / np.math.sqrt(np.dot(x_vec_, x_vec_))
#print 'x_vec_normed', x_vec_normed

y_vec_ = daf[10, -1, :3] - daf[10, 0, :3]
y_vec_normed = y_vec_ / np.math.sqrt(np.dot(y_vec_, y_vec_))

z_vec_normed = np.cross(x_vec_normed, y_vec_normed)

x = np.array([1, 0, 0])
y = np.array([0, 1, 0])
z = np.array([0, 0, 1])

cos_xx_ = np.dot(x_vec_normed, x)
#angle_xx_ = np.arccos(cos_xx_)
#print 'angle_xx_',angle_xx_ 

cos_yx_ = np.dot(x_vec_normed, y)
#angle_yx_ = np.arccos(cos_yx_)
#print 'angle_yx_',angle_yx_ 

cos_zx_ = np.dot(x_vec_normed, z)
#angle_zx_ = np.arccos(cos_zx_)
#print 'angle_zx_',angle_zx_ 

cos_xy_ = np.dot(y_vec_normed, x)
#angle_xy_ = np.arccos(cos_xy_)
#print 'angle_xy_',angle_xy_ 

cos_yy_ = np.dot(y_vec_normed, y)
#angle_yy_ = np.arccos(cos_yy_)
#print 'angle_yy_',angle_yy_ 

cos_zy_ = np.dot(y_vec_normed, z)
#angle_zy_ = np.arccos(cos_zy_)
#print 'angle_zy_',angle_zy_ 

cos_xz_ = np.dot(z_vec_normed, x)
#angle_xz_ = np.arccos(cos_xz_)
#print 'angle_xz_',angle_xz_ 

cos_yz_ = np.dot(z_vec_normed, y)
#angle_yz_ = np.arccos(cos_yz_)
#print 'angle_yz_',angle_yz_ 

cos_zz_ = np.dot(z_vec_normed, z)
#angle_zz_ = np.arccos(cos_zz_) 
#print 'angle_zz_',angle_zz_

T_mtx = np.array([[cos_xx_, cos_yx_, cos_zx_],
               [cos_xy_, cos_yy_, cos_zy_],
               [cos_xz_, cos_yz_, cos_zz_]])

#daf_new[:,:,:3] = np.dot( T_mtx[np.newaxis, np.newaxis,:,:], daf[:,:,:3])

# rotatation:
#
daf_new[:, :, 0] = daf_new[:, :, 0] * cos_xx_ + daf_new[:, :, 1] * cos_yx_ + daf_new[:, :, 2] * cos_zx_
daf_new[:, :, 1] = daf_new[:, :, 0] * cos_xy_ + daf_new[:, :, 1] * cos_yy_ + daf_new[:, :, 2] * cos_zy_
daf_new[:, :, 2] = daf_new[:, :, 0] * cos_xz_ + daf_new[:, :, 1] * cos_yz_ + daf_new[:, :, 2] * cos_zz_

# translation:
#
daf_new[:, :, 0] = daf_new[:, :, 0] - daf[0, 0, 0]
daf_new[:, :, 1] = daf_new[:, :, 1] - daf[0, 0, 1]
daf_new[:, :, 2] = daf_new[:, :, 2] - daf[0, 0, 2]


















ux_arr = daf[:, :, 3]
ux_masked = ma.masked_array(ux_arr, mask=mask_idx_array)

ux_avg = ma.average(ux_masked, axis=1)

x_idx_zeros, y_idx_zeros = np.where(mask_idx_array)

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
d4_ux_avg = np.average(d4_ux_arr, axis=1)

dd44_ux_arr = np.zeros_like(ux_arr)
dd44_ux_arr[integ_radius:-integ_radius, :] = (d4_ux_arr[integ_range:, :] -
                                                      d4_ux_arr[:-integ_range, :])
dd44_ux_avg = np.average(dd44_ux_arr, axis=1)

# third difference field
ddd444_ux_arr = np.zeros_like(ux_arr)
ddd444_ux_arr[integ_radius:-integ_radius, :] = (dd44_ux_arr[integ_range:, :] -
                                                        dd44_ux_arr[:-integ_range, :])
ddd444_ux_avg = np.average(ddd444_ux_arr, axis=1)

#crack_idx = ((ddd444_ux_avg[1:] * ddd444_ux_avg[:-1] < 0.0) * 
#             (ddd444_ux_avg[1:] < ddd444_ux_avg[:-1])) 

crack_idx_arr = ((dd44_ux_arr[1:, :] * dd44_ux_arr[:-1, :] < 0.0) *
                 ((ddd444_ux_arr[1:, :] + ddd444_ux_arr[:-1, :]) / 2.0 < -0.005))

crack_idx_avg = ((dd44_ux_avg[1:] * dd44_ux_avg[:-1] < 0.0) *
             ((ddd444_ux_avg[1:] + ddd444_ux_avg[:-1]) / 2.0 < -0.005))

w_arr = d4_ux_arr[np.where(crack_idx_arr)]
w_avg = np.max(d4_ux_arr[np.where(crack_idx_avg)], axis=1)

print 'crack_openings', w_avg

crack_avg = np.array(crack_idx_avg, dtype=float) * 0.2

row_idx = (10, n_y - int(n_y / 2), n_y - 10)
plot_row_idx = 10

x_arr = daf[:, plot_row_idx, 0]
x4_arr = daf[integ_radius:-integ_radius, plot_row_idx, 0]
x44_arr = daf[2 * integ_radius:-2 * integ_radius, plot_row_idx, 0]
x444_arr = daf[3 * integ_radius:-3 * integ_radius, plot_row_idx, 0]
xcrack = x444_arr[:-1] + d_x / 2.0

idx_arr = np.arange(ux_arr.shape[0], dtype=int)

p.figure()
p.subplot(2, 2, 1)
#p.plot(idx_arr, d4_ux_arr[:, row_idx], color = 'green')
p.plot(idx_arr, d4_ux_avg, color='black')
p.plot(idx_arr[:-1], crack_avg, color='magenta', linewidth=4)

p.subplot(2, 2, 3)
#p.plot(idx_arr, dd44_ux_arr[:, row_idx], color = 'green')
p.plot(idx_arr, dd44_ux_avg, color='black')

#p.plot(idx_arr, ddd444_ux_arr[:, row_idx], color = 'red')
p.plot(idx_arr, ddd444_ux_avg, color='blue')
p.plot(idx_arr[:-1], crack_avg, color='magenta', linewidth=4)

p.subplot(2, 2, 2)
p.plot(idx_arr, ux_arr[:, :], color='green')
p.plot(idx_arr, ux_avg, color='red')

#p.figure()
#p.hist(w_avg)

p.subplot(2, 2, 4)
#p.hist(w_arr, bins=40)

p.show()

#qargs = [ daf[:, :, i, np.newaxis] for i in range(0, 6)]

#m.quiver3d(*qargs)

print 'X_arr', X_arr.shape
print 'y_arr', y_arr.shape
print 'z_arr', z_arr.shape
print 'd4_ux_arr', ux_arr.shape

#s = m.points3d(X_arr,
#               y_arr,
#               z_arr,
#               ux_arr)

## print a surface - why is moved in the y direction?
#m.surf(daf[:-integ_range, :, 0], daf[:-integ_range, :, 1], d4_ux_arr)#, mask = mask_idx_array)

m.show()

# @todo: 
# sensitivity with respect to the integ_range
# 
