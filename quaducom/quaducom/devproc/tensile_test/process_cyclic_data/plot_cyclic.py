'''
Created on Aug 5, 2015

@author: Yingxiong
'''
import numpy as np
import matplotlib.pyplot as plt
from matresdev.db.simdb import SimDB
simdb = SimDB()
import os

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
f_file = os.path.join(test_file_path, test_file + '-force-binary')
d_file = os.path.join(test_file_path, test_file + '-disp-binary')

force = np.fromfile(f_file)
disp = np.fromfile(d_file)
print len(disp)
plt.plot(np.arange(800000), disp[-800000::])

# v1============================================================
fe0 = force[0:134060]
fe1 = force[200700:216400]
fe2 = force[307400:323250]
fe3 = force[391950:407700]
fe4 = force[754250:770000]
fe5 = force[3487829 + 540150:3487829 + 555900]

de0 = disp[0:134060]
de1 = disp[200700:216400]
de2 = disp[307400:323250]
de3 = disp[391950:407700]
de4 = disp[754250:770000]
de5 = disp[3487829 + 540150:3487829 + 555900]

# A_tex = 14 Rov * 1.84 mm^2 / Rov = 25.76
#
plt.plot(avg(de0) / 250., avg(fe0) / 25.76 * 1000., 'k')
plt.plot(avg(de1) / 250., avg(fe1) / 25.76 * 1000., 'k')
plt.plot(avg(de2) / 250., avg(fe2) / 25.76 * 1000., 'k')
plt.plot(avg(de3) / 250., avg(fe3) / 25.76 * 1000., 'k')
plt.plot(avg(de4) / 250., avg(fe4) / 25.76 * 1000., 'k')
plt.plot(avg(de5) / 250., avg(fe5) / 25.76 * 1000., 'k')

plt.xlabel('strain')
plt.ylabel('stress [MPa]')


se0 = avg(de0) / 250.
se1 = avg(de1) / 250.
se2 = avg(de2) / 250.
se3 = avg(de3) / 250.
se4 = avg(de4) / 250.
se5 = avg(de5) / 250.


plt.figure()
plt.plot(
    [0, 1, 2, 3, 4, 5], [se0[-1], se1[-1], se2[-1], se3[-1], se4[-1], se5[-1]], marker='o', label='upper strain')

plt.plot(np.arange(6), [min(se0[-1300::]), min(se1), min(se2),
                        min(se3), min(se4), min(se5)], marker='o', label='lower strain')
plt.ylim(0, 0.01)
plt.xlabel('10e')
plt.legend()


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
# plt.plot(avg(de0) / 250., avg(fe0) / 25.76 * 1000., 'k')
# plt.plot(avg(de1) / 250., avg(fe1) / 25.76 * 1000., 'k')
# plt.plot(avg(de2) / 250., avg(fe2) / 25.76 * 1000., 'k')
# plt.plot(avg(de3) / 250., avg(fe3) / 25.76 * 1000., 'k')
# plt.plot(avg(de4) / 250., avg(fe4) / 25.76 * 1000., 'k')
#
# plt.xlabel('strain')
# plt.ylabel('stress [MPa]')
#
#
# se0 = avg(de0) / 250.
# se1 = avg(de1) / 250.
# se2 = avg(de2) / 250.
# se3 = avg(de3) / 250.
# se4 = avg(de4) / 250.
#
#
# plt.figure()
# plt.plot(
#     [0, 1, 2, 3, 4], [se0[-1], se1[-1], se2[-1], se3[-1], se4[-1]], marker='o', label='upper strain')
#
# plt.plot(np.arange(5), [min(se0[-1300::]), min(se1), min(se2),
#                         min(se3), min(se4)], marker='o', label='lower strain')
# plt.ylim(0, 0.01)
# plt.xlabel('10e')
# plt.legend()


plt.show()
