'''
Created on May 31, 2012

@author: rostar
'''
import pickle
from matplotlib import pyplot as plt
import numpy as np

file = open('w5', 'r')
w5 = pickle.load(file)
file = open('w5_2', 'r')
w5_2 = pickle.load(file)
file = open('w5_3', 'r')
w5_3 = pickle.load(file)
file = open('w3_1', 'r')
w3_1 = pickle.load(file)
file = open('w3_2', 'r')
w3_2 = pickle.load(file)
file = open('w3_3', 'r')
w3_3 = pickle.load(file)
file = open('w1_1', 'r')
w1_1 = pickle.load(file)
file = open('w1_2', 'r')
w1_2 = pickle.load(file)
file = open('w1_3', 'r')
w1_3 = pickle.load(file)
#plt.hist(w5 + w5_2 + w5_3, bins = 21, normed = False, color = 'white', lw = 2)
#plt.hist(w3_1 + w3_2 + w3_2, bins = 25, normed = False, color = 'white', lw = 2)
plt.hist(w1_1 + w1_2 + w1_3, bins = 25, normed = False, color = 'white', lw = 2)
plt.xlim(0.04, 0.19)
plt.ylim(0.0, 100)
plt.show()