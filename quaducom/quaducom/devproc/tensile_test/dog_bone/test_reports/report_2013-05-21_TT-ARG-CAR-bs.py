'''
Created on Aug 26, 2013

@author: alexander
'''
from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt
from matplotlib.font_manager import FontProperties
from numpy import \
    sin, cos, c_, arange, hstack, array, max, frompyfunc, linspace
import os
from os.path import join
from pickle import dump, load
from time import time
from traits.api import \
    Array, Bool, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, \
    Callable, List, TraitDict, Any, Range, \
    Delegate, Event, on_trait_change, Button, \
    Interface, implements, Property, cached_property

from .format_plot import format_plot
from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray
from matresdev.db.exdb.ex_run_view import \
    ExRunView
from matresdev.db.matdb.trc.ccs_unit_cell import \
    CCSUnitCell, DamageFunctionEntry
from matresdev.db.simdb import \
    SimDB
from matresdev.db.simdb.simdb_class import \
    SimDBClass, SimDBClassExt
import numpy as np


simdb = SimDB()


font = FontProperties()
font.set_name('Script MT')
font.set_family('serif')
font.set_style('normal')
font.set_variant('normal')
font.set_weight('medium')

if __name__ == '__main__':

    from matresdev.db.exdb.ex_run import ExRun
    import pylab as p

    do = 'TT-ARG-CAR-bs'
    sig_flag = 'comp'
#     sig_flag = 'tex'

    save_fig_to_file = True

    #---------------------------
    # tensile test results (CAR-800-TU / AR-1200-TU / AR-1200-TR )
    #---------------------------
    #
    fig = p.figure(edgecolor='w', facecolor='w')
    fig.set_size_inches(8, 6)

    #-------------------------------
    # carbon, tissue, 800tex (bs2)
    #-------------------------------
    #
    path = join(simdb.exdata_dir, 'tensile_tests',
                'dog_bone', '2013-05-21-TT-6c-2cm-0-TU_bs2')
    tests = ['TT-6c-2cm-0-TU-V1_bs2.DAT']  # , 'TT-6c-2cm-0-TU-V2_bs2.DAT']
    for t in tests:
        ex_path = join(path, t)
        ex_run = ExRun(ex_path)

        if sig_flag == 'comp':
            ex_run.ex_type._plot_comp_stress_strain_asc(
                p, plot_analytical_stiffness=False, color='black', linewidth=1.3, linestyle='--', label='CAR-800-TU')
        if sig_flag == 'tex':
            ex_run.ex_type._plot_tex_stress_strain_asc(
                p, color='black', linewidth=1.3, linestyle='--', label='CAR-800-TU')

    #-------------------------------
    # AR-glas, tissue, 1200tex
    #-------------------------------
    #
    path = join(simdb.exdata_dir, 'tensile_tests',
                'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs')
    tests = ['TT-6g-2cm-0-V2.DAT']  # 'TT-6g-2cm-0-V1.DAT'
    for t in tests:
        ex_path = join(path, t)
        ex_run = ExRun(ex_path)
        if sig_flag == 'comp':
            ex_run.ex_type._plot_comp_stress_strain_asc(
                p, xscale=1000., color='black', linewidth=1.3, linestyle='-.', label='ARG-1200-TU')
        if sig_flag == 'tex':
            ex_run.ex_type._plot_tex_stress_strain_asc(
                p, xscale=1000., color='black', linewidth=1.3, linestyle='-.', label='2D-09-12')

    #-------------------------------
    # AR-glas, tricot, 1200tex
    #-------------------------------
    #
#        path = join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-03-05_TT-10g-3cm-0-TR')
#        tests = ['TT-10g-3cm-0-TR-V2.DAT']#, 'TT-10g-3cm-0-TR-V3.DAT']
    path = join(simdb.exdata_dir, 'tensile_tests',
                'dog_bone', '2010-03-17_TT-8g-3cm-0-TR')
    # , 'TT-8g-3cm-0-TR-V2.DAT', 'TT-8g-3cm-0-TR-V3.DAT']
    tests = ['TT-8g-3cm-0-TR-V2.DAT']
    for t in tests:
        ex_path = join(path, t)
        ex_run = ExRun(ex_path)
        if sig_flag == 'comp':
            ex_run.ex_type._plot_comp_stress_strain_asc(
                p, xscale=1000., color='black', linewidth=1.3, linestyle='-', label='ARG-1200-TR')
        if sig_flag == 'tex':
            ex_run.ex_type._plot_tex_stress_strain_asc(
                p, xscale=1000., color='black', linewidth=1.3, linestyle='-', label='2D-01-08')

    if sig_flag == 'tex':
        format_plot(
            p, xlabel='Dehnung [1E-3]', ylabel='Textilspannung [MPa]', xlim=14., ylim=1400.)
    if sig_flag == 'comp':
        format_plot(
            p, xlabel='Dehnung [1E-3]', ylabel='Kompositspannung [MPa]', xlim=14., ylim=25.)

    font.set_size('12')
    p.legend(prop=font)

    if save_fig_to_file:
        # create a report-subfolder with the name of the script (without file extension '.py')
        # and save it in the test-type-subfolders with the name and path as
        # ex_type
        test_series_name = os.path.basename(__file__)[:-3]
        subfolder_list = __file__.split(os.path.sep)
        devproc_idx = np.where(np.array(subfolder_list) == 'devproc')[0]
        subfolder_path = subfolder_list[
            devproc_idx + 1:-2] + [test_series_name]
        test_series_dir = os.path.join(simdb.report_dir)
        for subfolder_name in subfolder_path:
            test_series_dir = os.path.join(test_series_dir, subfolder_name)
        if not os.path.exists(test_series_dir):
            os.makedirs(test_series_dir)

        filename = os.path.join(
            test_series_dir, do + '_sig' + sig_flag + '-epsu.png')
        p.savefig(filename, format='png')
        print('figure saved to file %s' % (filename))

    p.show()
