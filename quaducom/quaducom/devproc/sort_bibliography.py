'''
Created on Aug 26, 2013

@author: alexander
'''
import numpy as np

import os

from matresdev.db.exdb.ex_run_view import \
    ExRunView

from matresdev.db.simdb import \
    SimDB

from matresdev.db.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

simdb = SimDB()

from quaducom.devproc.format_plot import format_plot
from matplotlib.font_manager import FontProperties
font = FontProperties()

if __name__ == '__main__':

    import pylab as p

    img_dir = os.path.join(simdb.exdata_dir, 'img_dir')
    filename = os.path.join(img_dir, 'Literaturverzeichnis.tex')
    f = open(filename, 'r')
    bib_list = f.readlines()
    print bib_list
    print len(bib_list)

    sorted_bib_list = []
    sorted_bib_item = ''
    # remove lines with coments
    for name in bib_list:
        if len(name) < 5:
           print sorted_bib_item
           sorted_bib_list.append(sorted_bib_item + '\n')
           sorted_bib_item = ''
        else:
           sorted_bib_item += name

    bib_list_alpha = sorted(sorted_bib_list)
    print 'sorted_bib_list', bib_list_alpha

    bib_str = ''
    for i in bib_list_alpha:
        bib_str += i
    print 'len', len(bib_str)
    print 'XXX', bib_str

    img_dir = os.path.join(simdb.exdata_dir, 'img_dir')
    filename_out = os.path.join(img_dir, 'Literaturverzeichnis_sorted.tex')
    f_out = open(filename_out, 'w')
    f_out.write(bib_str)
