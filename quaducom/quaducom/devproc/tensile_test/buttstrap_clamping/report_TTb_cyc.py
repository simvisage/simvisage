'''
Created on Jan 28, 2015

'''

from matresdev.db.simdb import SimDB
simdb = SimDB()
from matresdev.db.exdb import ExRun

import time

import os
import numpy as np
import pylab as p
params = {'legend.fontsize': 10,
          # 'legend.linewidth': 2
          }
p.rcParams.update(params)


def comma2dot(c):
    '''convert float with comma separator into float with dot separator'''
    return float((str(c)).replace(",", "."))


def time2sec(date):
    '''convert time format (hh:mm:ss) to seconds (s)'''
    d_list = str(date).split()
    t_list = d_list[1].split(':')
    t_sec = int(t_list[0]) * 60 * 60 + int(t_list[1]) * 60 + int(t_list[2])
    return t_sec


def scaledown_data(data_arr, n_avg):
    '''scale down the number of rows in 'data_array' by the
    integer 'n_avg', i.e. if 'n_avg=2' reduce the number of rows in
    'data_arr' to half. The methods cuts of  up to 2*'n_avgs' rows from the end
    of the file in order to make sure that the sub arrays used for averaging have the same
    shape; in general that doesn't effect the data as the measuring is continued after rupture
    or high frequency measurements is used'''
    n_rows = data_arr.shape[0]
    n_steps = (n_rows - n_avg) / n_avg
    n_max = n_steps * n_avg
    avg_list = [data_arr[i:n_max:n_avg, :] for i in range(n_avg)]
    avg_arr = np.array(avg_list)
    data_arr_ = np.mean(avg_arr, 0)
    return data_arr_

# test_file_path = os.path.join(simdb.exdata_dir,
#                               'tensile_tests',
#                               'buttstrap_clamping',
#                               '2015-04-20_TTb-2C-1cm-0-800SBR_cyc'
#                               )
#
# test_files = [
#              'TTb-2C-1cm-0-800SBR-V1_cyc-Aramis2d.csv',
#             ]
#
#
# test_file = os.path.join(test_file_path, test_files[0])

test_file = 'D:\\2015-04-20_TTb-2C-1cm-0-800SBR_cyc\\TTb-2C-1cm-0-800SBR-V4_cyc-Aramis2d-dot.csv'

downscale_data = True
plot_data = True

t0 = time.time()

# use 'genfromtxt' if file containes blank entries
# data_arr = np.loadtxt(test_file, delimiter=";", skiprows=2,
data_arr = np.genfromtxt(test_file, delimiter=";", skiprows=2,
                         converters={0: time2sec, 1: comma2dot, 2: comma2dot, 3: comma2dot,
                                     4: comma2dot, 5: comma2dot, 6: comma2dot})

print 'data_arr', data_arr.shape
# print 'data_arr[0]', data_arr[0]
dlist = [np.array([data_arr[i][0], data_arr[i][1], data_arr[i][2], data_arr[i][
                  3], data_arr[i][4], data_arr[i][5], data_arr[i][6]]) for i in range(data_arr.shape[0])]
data_arr = np.vstack(dlist)

t1 = time.time()
print 'time needed for loading original data array [min]: %g' % ((t1 - t0) / 60)

if downscale_data:
    n_avg_list = [2, 4, 10, 20, 40, 100, 200]
    for n_avg in n_avg_list:
        t0 = time.time()
        data_arr_ = scaledown_data(data_arr, n_avg)
        t1 = time.time()
        dt = (t1 - t0) / 60
        print 'time needed for down scaling of data array [min]: %g' % dt
        print 'data_arr_.shape', data_arr_.shape

        test_file_ = os.path.join(
            test_file_path, test_files[0].split('.')[0] + '_navg' + str(n_avg) + '.csv')
        print 'test_file_', test_file_

        header_string = '''Datum/Uhrzeit;Kraft;Clamping;Weg;WA1_vorne;WA2_links;WA3_rechts\r\n;kN;kN;mm;mm;mm;mm\r\n'''
        np.savetxt(test_file_, data_arr_, delimiter=";")
        # header=header_string) # use header_string option of new loadtxt
        # version if available

        #-------
        # workaround for old numpy.savetxt version:
        t0 = time.time()
        f = open(test_file_, 'r')
        temp = f.read()
        f.close()
        f = open(test_file_, 'w')
        f.write(header_string)
        f.write(temp)
        f.close()
        t1 = time.time()
        print 'time needed for adding header [min]: %g' % ((t1 - t0) / 60)
        #-------

if plot_data:
    e_list = [ExRun(data_file=os.path.join(test_file_path, test_file))
              for test_file in test_files]

    fig = p.figure(facecolor='white')
    fig.subplots_adjust(
        left=0.09, right=0.97, bottom=0.14, top=0.96, wspace=0.25, hspace=0.2)
    p.legend()

    for e_run in e_list:
        e = e_run.ex_type
        axes = p.subplot(121)
        p.plot(e.time_asc, e.F_asc)
        axes = p.subplot(122)
        p.plot(e.eps_asc, e.sig_c_asc)
        # for comparison of different compression levels (endings by "_navg200.csv" etc.)
    #    p.plot(e.eps_asc, e.sig_c_asc, label=e.key[-6:])
    #    print 'label=e.key[-6:]', e.key[-6:]
    #    p.legend()

    #    e._plot_comp_stress_strain_asc(axes)
    #    e._plot_tex_stress_strain_asc(axes)

    p.show()
