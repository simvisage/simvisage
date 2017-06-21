'''
Created on Jan 28, 2015

'''
import os

from matresdev.db.exdb import ExRun
from matresdev.db.simdb import SimDB
import numpy as np
import pylab as p
simdb = SimDB()

import quaducom.devproc.bending_test_retrofitting.four_point.exp_bt_4pt_rf

#--------------------
# script fills up empty values (last available value of lines above) if:
# 1) if measuring error occured and single values are not recorded (blank fields)
# 2) as the available measuring frequency for the displacement gauges and the strain gauges
#    varies the script need to fill up values;
# 3) the scipt also merges the first two columns ('date' and 'time' to one column in order to
#    provide the default format for 'ex_type' reading stanard PEEKEL/csv.-files
#--------------------

replace_missing_values = False
merge_data_and_time = True

#--------------------
# NOTE: the script may take a while (for approx. 189.000 lines several minutes)
# @todo: check efficiency of the scripting approach
#--------------------

#---------------------------------------
# test data file
#---------------------------------------

# 'RF1_ref-mon'
# -------------
data_file_path = os.path.join(simdb.exdata_dir,
                            'bending_tests_retrofitting',
                            'three_point',
                            '2016-01-21_BT-3PT-1s-20cm-d8mm-RF1_ref-mon-Aramis2d')
data_file_name = 'BT-3PT-1s-20cm-d8mm-RF1_ref-mon-Aramis2d.csv'

# 'RF3_2C-mon'
# -------------
data_file_path = os.path.join(simdb.exdata_dir,
                            'bending_tests_retrofitting',
                            'three_point',
                            '2016-01-22_BT-3PT-1s-20cm-d8mm-RF3_2C-mon-Aramis2d')
data_file_name = 'BT-3PT-1s-20cm-d8mm-RF3_2C-mon-Aramis2d.csv'

# 'RF2_2C-cyc'
# -------------
data_file_path = os.path.join(simdb.exdata_dir,
                            'bending_tests_retrofitting',
                            'three_point',
                            '2016-01-25_BT-3PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d')
data_file_name = 'BT-3PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d.csv'

# 'RF2_2C-mon'
# -------------
data_file_path = os.path.join(simdb.exdata_dir,
                            'bending_tests_retrofitting',
                            'three_point',
                            '2016-01-27_BT-3PT-1s-20cm-d8mm-RF2_2C-mon-Aramis2d')
data_file_name = 'BT-3PT-1s-20cm-d8mm-RF2_2C-mon-Aramis2d.csv'

# check for original data in subfolder 'original_data'
#
data_file = os.path.join(data_file_path, 'original_data', data_file_name)

# read data as list of strings using 'readlines()'
#
data_file_ = open(data_file, 'r')
data_file_str = data_file_.readlines()

#---------------------------------------
# merge first two columns
#---------------------------------------
if merge_data_and_time:
    # if two columns are used for 'date' and 'time' merge them to a single column 'date/time'
    # which is the standard format 'ex_type' expects for  PEEKEL/csv-files
    #
    data_file_str[0] = data_file_str[0].replace("Datum;Zeit;", "Datum/Uhrzeit;")
    data_file_str[1] = data_file_str[1].replace(";;kN", ";kN")
    data_file_str = [data_file_str[i].replace(".2016;", ".2016 ") for i in range(len(data_file_str))]

#---------------------------------------
# fill up blank fields with last available value
#---------------------------------------
if replace_missing_values:
    n_rows = len(data_file_str)
    for i in range(3, n_rows):
    #     if i % 1000 == 0:
    #         print '%i lines from total %i processed' % (i, n_rows)
        data_file_str_i_split = data_file_str[i].split(';')
        for j in range(len(data_file_str_i_split)):
            # check for missing values (= white spaces between semicolons or at end of line)
            if data_file_str_i_split[j].strip() == '':
                # if whitespace is found, replace missing value with in current line with
                # value of previous line
                data_file_str_i_split[j] = data_file_str[i - 1].split(';')[j]
                data_file_str[i] = ''
                for j in range(len(data_file_str_i_split)):
                    data_file_str[i] += data_file_str_i_split[j] + ';'
                # remove trailing semicolon
                data_file_str = data_file_str[:-1]

#---------------------------------------
# save processed data to file
#---------------------------------------
# using 'writelines()' (create file if it doesn't exist)
#
data_file_processed = os.path.join(data_file_path, data_file_name)
data_file_processed_ = open(data_file_processed, 'w+')
data_file_processed_.writelines(data_file_str[:])
