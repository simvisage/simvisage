'''
Created on Jun 10, 2015

@author: Yingxiong
'''
import numpy as np
test_file = 'D:\\2015-04-20_TTb-2C-1cm-0-800SBR_cyc\\TTb-2C-1cm-0-800SBR-V5_cyc-Aramis2d-dot.csv'
data_arr = np.genfromtxt(test_file, delimiter=";", skiprows=2)
# data = np.memmap(test_file, dtype='float32')
print(data_arr[-1])
# with open(test_file, 'r+') as ins:
#     for line in ins:
#         print line
