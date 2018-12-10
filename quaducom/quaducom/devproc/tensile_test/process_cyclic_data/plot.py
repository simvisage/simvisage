'''
Created on Aug 5, 2015

@author: Yingxiong
'''
import numpy as np
import matplotlib.pyplot as p
from matresdev.db.simdb import SimDB
simdb = SimDB()
import os

def scaledown_data(data_arr, n_avg):
    '''scale down the number of rows in 'data_array' by the
    integer 'n_avg', i.e. if 'n_avg=2' reduce the number of rows in
    'data_arr' to half. The methods cuts of  up to 2*'n_avgs' rows from the end
    of the file in order to make sure that the sub arrays used for averaging have the same
    shape; in general that doesn't effect the data as the measuring is continued after rupture
    or high frequency measurements is used'''
    n_rows = data_arr.shape[0]
    print('n_rows', n_rows)
    n_steps = (n_rows - n_avg) / n_avg
    print('n_steps', n_steps)
    n_max = n_steps * n_avg
    print('n_max', n_max)
    avg_list = [data_arr[i:n_max:n_avg] for i in range(n_avg)]
    avg_arr = np.array(avg_list)
    data_arr_ = np.mean(avg_arr, 0)
    return data_arr_

def avg(a):
    n = 10
    return np.mean(a.reshape(-1, n), axis=1)

test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2015-08-03_TTb-2C-14mm-0-3300SBR_cyc-Aramis2d'
                              )

test_file = 'TTb-2C-14mm-0-3300SBR-V1_cyc-Aramis2d'
# test_file = 'TTb-2C-14mm-0-3300SBR-V6_cyc-Aramis2d'
# test_file = 'TTb-2C-14mm-0-3300SBR-V4_cyc-Aramis2d'
d_file = os.path.join(test_file_path, test_file + '-disp-binary')
f_file = os.path.join(test_file_path, test_file + '-force-binary')

disp = np.fromfile(d_file)
force = np.fromfile(f_file)

# gauge_length = 250mm; A_tex = 2layers * 7 Rov * 1.84 mm^2 / Rov = 25.76 mm^2
#
eps = disp / 250. * 1000.  # [permille]
sigtex = force / 25.76 * 1000  # [MPa]

fig = p.figure(facecolor='white', figsize=(16, 9))
fig.subplots_adjust(wspace=0.25, hspace=0.4)

# disp_ = scaledown_data(disp, 20)
# force_ = scaledown_data(force, 20)
# print 'force.shape[0]', force.shape[0]
# print 'force_.shape[0]', force_.shape[0]

######################
# n_start = 0
# n_stop = 1000000
# n_avg = 200  # reduce the measuring data by this factor (use average instead)
# disp_ = scaledown_data(disp[n_start:n_stop], n_avg)
# print 'np.mean(disp_)', np.mean(disp_)
# print 'np.average(disp_)', np.average(disp_)
# print 'np.max(disp_)', np.max(disp_)
# narr_ = np.arange(len(disp_))
# p.subplot(121)
# p.plot(narr_, disp_)
#
# # d_disp_ = disp_[1:] - disp_[:-1]
# # narr_d = np.arange(len(d_disp_))
# # p.subplot(122)
# # p.plot(narr_d, d_disp_)
# p.show()
######################

#----------------------------
# Sec-1 (monotone loading)
#----------------------------
n_start = 0
n_stop = 110000
n_avg = 200  # reduce the measuring data by this factor (use average instead)
#
eps_ = scaledown_data(eps[n_start:n_stop], n_avg)
sigtex_ = scaledown_data(sigtex[n_start:n_stop], n_avg)
p.subplot(121)
p.plot(eps_, sigtex_)
p.xlabel('strain')
p.ylabel('stress [MPa]')
#
# disp_ = scaledown_data(disp[n_start:n_stop], n_avg)
# narr_ = np.arange(len(disp_))
# p.subplot(122)
# p.plot(narr_, disp_)
# p.xlabel('n [-]')
# p.ylabel('displacement [mm]')

#----------------------------
# Sec-2 (cyclic loading)
#----------------------------
n_start = 110000
n_stop = 551500
n_avg = 50
#
p.subplot(121)
eps_ = scaledown_data(eps[n_start:n_stop], n_avg)
sigtex_ = scaledown_data(sigtex[n_start:n_stop], n_avg)
p.plot(eps_, sigtex_)
#
# disp_ = scaledown_data(disp[n_start:n_stop], n_avg)
# narr_ = np.arange(len(disp_)) + len(narr_)
# p.subplot(122)
# p.plot(narr_, disp_)

#----------------------------
# Sec-3 (cyclic loading)
#----------------------------
# n_start = 151500
# n_stop = 360000
# n_avg = 50
# #
# p.subplot(121)
# eps_ = scaledown_data(eps[n_start:n_stop], n_avg)
# sigtex_ = scaledown_data(sigtex[n_start:n_stop], n_avg)
# p.plot(eps_, sigtex_)
# #
# disp_ = scaledown_data(disp[n_start:n_stop], n_avg)
# narr_ = np.arange(len(disp_)) + len(narr_)
# p.subplot(122)
# p.plot(narr_, disp_)
#----------------------------
p.show()


## v1============================================================
# fe0 = force[0:134060]
# fe1 = force[200700:216400]
# fe2 = force[307400:323250]
# fe3 = force[391950:407700]
# fe4 = force[754250:770000]
# fe5 = force[3487829 + 540150:3487829 + 555900]
#
# de0 = disp[0:134060]
# de1 = disp[200700:216400]
# de2 = disp[307400:323250]
# de3 = disp[391950:407700]
# de4 = disp[754250:770000]
# de5 = disp[3487829 + 540150:3487829 + 555900]
#
# # A_tex = 14 Rov * 1.84 mm^2 / Rov = 25.76
# #
# p.plot(avg(de0) / 250., avg(fe0) / 25.76 * 1000., 'k')
# p.plot(avg(de1) / 250., avg(fe1) / 25.76 * 1000., 'k')
# p.plot(avg(de2) / 250., avg(fe2) / 25.76 * 1000., 'k')
# p.plot(avg(de3) / 250., avg(fe3) / 25.76 * 1000., 'k')
# p.plot(avg(de4) / 250., avg(fe4) / 25.76 * 1000., 'k')
# p.plot(avg(de5) / 250., avg(fe5) / 25.76 * 1000., 'k')
#
# p.xlabel('strain')
# p.ylabel('stress [MPa]')
#
#
# se0 = avg(de0) / 250.
# se1 = avg(de1) / 250.
# se2 = avg(de2) / 250.
# se3 = avg(de3) / 250.
# se4 = avg(de4) / 250.
# se5 = avg(de5) / 250.
#
#
# p.figure()
# p.plot(
#    [0, 1, 2, 3, 4, 5], [se0[-1], se1[-1], se2[-1], se3[-1], se4[-1], se5[-1]], marker='o', label='upper strain')
#
# p.plot(np.arange(6), [min(se0[-1300::]), min(se1), min(se2),
#                        min(se3), min(se4), min(se5)], marker='o', label='lower strain')
# p.ylim(0, 0.01)
# p.xlabel('10e')
# p.legend()


# v6 ======================================

# fe0 = force[0:133200]
# fe1 = force[176080:191840]
# fe2 = force[232020:247760]
# fe3 = force[321100:336950]
# fe4 = force[690670:706500]
# #
# de0 = disp[0:133200]
# de1 = disp[176080:191840]
# de2 = disp[232020:247760]
# de3 = disp[321100:336950]
# de4 = disp[690670:706500]
#
# p.plot(avg(de0) / 250., avg(fe0) / 25.76 * 1000., 'k')
# p.plot(avg(de1) / 250., avg(fe1) / 25.76 * 1000., 'k')
# p.plot(avg(de2) / 250., avg(fe2) / 25.76 * 1000., 'k')
# p.plot(avg(de3) / 250., avg(fe3) / 25.76 * 1000., 'k')
# p.plot(avg(de4) / 250., avg(fe4) / 25.76 * 1000., 'k')
#
# p.xlabel('strain')
# p.ylabel('stress [MPa]')
#
#
# se0 = avg(de0) / 250.
# se1 = avg(de1) / 250.
# se2 = avg(de2) / 250.
# se3 = avg(de3) / 250.
# se4 = avg(de4) / 250.
#
#
# p.figure()
# p.plot(
#     [0, 1, 2, 3, 4], [se0[-1], se1[-1], se2[-1], se3[-1], se4[-1]], marker='o', label='upper strain')
#
# p.plot(np.arange(5), [min(se0[-1300::]), min(se1), min(se2),
#                         min(se3), min(se4)], marker='o', label='lower strain')
# p.ylim(0, 0.01)
# p.xlabel('10e')
# p.legend()


# p.show()
