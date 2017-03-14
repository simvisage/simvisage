'''
Created on Jan 28, 2015

'''
import os

from matresdev.db.exdb import ExRun
from matresdev.db.simdb import SimDB
import numpy as np
import pylab as p
simdb = SimDB()

#--------------------
# script fills up empty values (last available value of lines above) if:
# 1) if measuring error occured and single values are not recorded (blank fields)
# 2) as the available measuring frequency for the displacement gauges and the strain gauges
#    varies the script need to fill up values;
# 3) the scipt also merges the first two columns ('date' and 'time' to one column in order to
#    provide the default format for 'ex_type' reading stanard PEEKEL/csv.-files
#--------------------

replace_missing_values = True
merge_data_and_time = True

#--------------------
# NOTE: the script may take a while (for approx. 189.000 lines several minutes)
# @todo: check efficiency of the scripting approach
#--------------------

#---------------------------------------
# test data file
#---------------------------------------
# read data as list of strings using 'readlines()'
#
data_file = os.path.join(simdb.exdata_dir,
                            'bending_tests_retrofitting',
                            'four_point',
                            '2016-01-19_BT-4PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d',
                            'original_data',
                            'BT-4PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d_ALL-CHANNELS_n0-n191560.csv')

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
data_file_processed = os.path.join(simdb.exdata_dir,
                                   'bending_tests_retrofitting',
                                   'four_point',
                                   '2016-01-19_BT-4PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d',
                                   'BT-4PT-1s-20cm-d8mm-RF2_2C-cyc-Aramis2d.csv')

data_file_processed_ = open(data_file_processed, 'w+')
data_file_processed_.writelines(data_file_str[:])
