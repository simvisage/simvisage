'''
Created on Aug 26, 2013

@author: alexander
'''
from etsproxy.traits.api import \
    Array, Bool, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, \
    Callable, List, TraitDict, Any, Range, \
    Delegate, Event, on_trait_change, Button, \
    Interface, implements, Property, cached_property

from ibvpy.api import \
    TStepper as TS, TLoop, TLine, \
    IBVModel, DOTSEval, \
    RTraceGraph, RTraceDomainListField, \
    BCDof, BCDofGroup, BCSlice, FERefinementGrid, FEDomain

from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U

from ibvpy.mesh.fe_grid import \
    FEGrid

from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
    MATS2D5MicroplaneDamage

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
    IPhiFn, PhiFnGeneralExtended, \
    PhiFnGeneral, PhiFnStrainHardening, PhiFnStrainHardeningLinear, \
    PhiFnGeneralExtendedExp

from mathkit.geo.geo_ndgrid import \
    GeoNDGrid

from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import \
    MFnNDGrid, GridPoint

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

from numpy import \
    sin, cos, c_, arange, hstack, array, max, frompyfunc, linspace

import numpy as np

from time import time
from os.path import join

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from simiter.sim_pstudy import \
    SimPStudy, SimOut, ISimModel

from matresdev.db.exdb.ex_run_view import \
    ExRunView

from matresdev.db.matdb.trc.ccs_unit_cell import \
    CCSUnitCell, DamageFunctionEntry

from matresdev.db.simdb import \
    SimDB

from matresdev.db.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

simdb = SimDB()

from pickle import dump, load

from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_name('Script MT')
font.set_family('serif')
font.set_style('normal')
font.set_size('large')
font.set_variant('normal')
font.set_weight('medium')

def format_plot(axes, xlim = None, ylim = None, xlabel = '', ylabel = ''):
    '''format 2d-plot black and with with times legends 
    '''
    #-------------------------------------------------------------------
    # configure the style of the font to be used for labels and ticks
    #-------------------------------------------------------------------
    #
    from matplotlib.font_manager import FontProperties
    font = FontProperties()
    font.set_name('Script MT')
    font.set_family('serif')
    font.set_style('normal')
#    font.set_size('small')
    font.set_size('large')
    font.set_variant('normal')
    font.set_weight('medium')
    
    if xlim != None and ylim != None:
        axes.axis([0, xlim, 0., ylim], fontproperties=font)

    # format ticks for plot
    #
    locs,labels = axes.xticks()
    axes.xticks(locs, map(lambda x: "%.0f" % x, locs), fontproperties=font)
    axes.xlabel(xlabel, fontproperties=font)

    locs,labels = axes.yticks()
    axes.yticks(locs, map(lambda x: "%.0f" % x, locs), fontproperties=font)
    axes.ylabel(ylabel, fontproperties=font)

if __name__ == '__main__':

    from matresdev.db.exdb.ex_run import ExRun
    import pylab as p

    do = 'show_test_results_ST'
#    do = 'show_test_results_SH'
#    do = 'show_test_results_TT-CAR'
#    do = 'show_test_results_TT'


    #---------------------------
    # shell test results
    #---------------------------
    #
    if do == 'show_test_results_ST':

        # carbon
        #
        path = join(simdb.exdata_dir, 'slab_tests', '2013-07-10_ST-6c-2cm-TU_bs2')
        tests = ['ST-6c-2cm-TU_bs2.DAT']
        
        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
            ex_run.ex_type._plot_force_center_deflection( p )

        format_plot(p, xlabel = 'strain [1E-3]', ylabel = 'applied force [kN]', xlim = 80., ylim = 20.)

#        p.legend(prop = font)
        
        p.show()

    #---------------------------
    # shell test results
    #---------------------------
    #
    if do == 'show_test_results_SH':

        # carbon
        #
        path = join(simdb.exdata_dir, 'shell_tests', '2013-02-27_SH-6c-2cm-TU_bs')
        tests = ['AS_BS_6C.DAT']
        
        # AR-glas
        #
#        path = join(simdb.exdata_dir, 'shell_tests', '2013-03-05_SH-6g-2cm-TU_bs')
#        tests = ['AS_BS_6G.DAT']#, 'TT-6c-2cm-0-TU-V2_bs2.DAT']

        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
#            ex_run.ex_type._plot_force_deflection_center( p, color = 'black', linewidth = 1., linestyle = '-', label = '2D-05-11' )
            ex_run.ex_type._plot_force_eps_t( p, linewidth = 1., linestyle = '-', label = '2D-05-11' )

        format_plot(p, xlabel = 'strain [1E-3]', ylabel = 'applied force [kN]', xlim = 10., ylim = 100.)

        p.legend(prop = font)
        
        p.show()


    #---------------------------
    # tensile test results (CAR-800tex-TU)
    #---------------------------
    #
    if do == 'show_test_results_TT-CAR':

        fig = p.figure()
        fig.set_size_inches(8, 6)

        # carbon, tissue, 800tex
        #
        path = join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-05-21-TT-6c-2cm-0-TU_bs2')
        tests = ['TT-6c-2cm-0-TU-V1_bs2.DAT']#, 'TT-6c-2cm-0-TU-V2_bs2.DAT']
        for t in tests:
            ex_path = join(path, t)
            ex_run = ExRun(ex_path)
#            ex_run.ex_type._plot_tex_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '--', label = '2D-05-11' )
#            ex_run.ex_type._plot_comp_stress_strain_asc( p, color = 'black', linewidth = 1., linestyle = '-', label = 'bs2', xscale = 1000. )
            ex_run.ex_type._plot_comp_stress_strain_asc( p, color = 'black', linewidth = 1., linestyle = '-', label = 'bs2', xscale = 1000. )

#        path = join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-06-12_TT-6c-2cm-0-TU_bs3')
#        tests = ['TT-6c-2cm-0-TU-V1_bs3.DAT', 'TT-6c-2cm-0-TU-V2_bs3.DAT', 'TT-6c-2cm-0-TU-V3_bs3.DAT']
#        for t in tests:
#            ex_path = join(path, t)
#            ex_run = ExRun(ex_path)
##            ex_run.ex_type._plot_tex_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '--', label = '2D-05-11' )
#            ex_run.ex_type._plot_comp_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '--', label = 'bs3', xscale = 1000. )
#
#        path = join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d')
#        tests = ['TTb-6c-2cm-0-TU-V1_bs4.DAT', 'TTb-6c-2cm-0-TU-V2_bs4.DAT', 'TTb-6c-2cm-0-TU-V3_bs4.DAT']
#        for t in tests:
#            ex_path = join(path, t)
#            ex_run = ExRun(ex_path)
##            ex_run.ex_type._plot_tex_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '--', label = '2D-05-11' )
#            ex_run.ex_type._plot_comp_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '--', label = 'bs4', xscale = 1000. )
#
#        path = join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-18_TTb-6c-2cm-0-TU_bs5')
#        tests = ['TTb-6c-2cm-0-TU-V1_bs5.DAT', 'TTb-6c-2cm-0-TU-V2_bs5.DAT', 'TTb-6c-2cm-0-TU-V3_bs5.DAT']
#        for t in tests:
#            ex_path = join(path, t)
#            ex_run = ExRun(ex_path)
##            ex_run.ex_type._plot_tex_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '--', label = '2D-05-11' )
#            ex_run.ex_type._plot_comp_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '--', label = 'bs5', xscale = 1000. )
 
        p.legend(prop = font)

        p.show()

    
    #---------------------------
    # tensile test results (CAR-800-TU / AR-1200-TU / AR-1200-TR )
    #---------------------------
    #
    if do == 'show_test_results_TT':

        fig = p.figure()
        fig.set_size_inches(8, 6)

        # carbon, tissue, 800tex
        #
#        path = join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-05-21-TT-6c-2cm-0-TU_bs2')
#        tests = ['TT-6c-2cm-0-TU-V1_bs2.DAT']#, 'TT-6c-2cm-0-TU-V2_bs2.DAT']
#        for t in tests:
#            ex_path = join(path, t)
#            ex_run = ExRun(ex_path)
##            ex_run.ex_type._plot_tex_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '--', label = '2D-05-11' )
##            ex_run.ex_type._plot_comp_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '--', label = 'CAR-800-TU', xscale = 1.  )
#            ex_run.ex_type._plot_tex_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '--', label = 'CAR-800-TU', xscale = 1.  )
#
#        ex_path = join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d',
#                        'TTb-6c-2cm-0-TU-V2_bs4.DAT')
#        ex_run = ExRun(ex_path)
#        ex_run.ex_type._plot_tex_stress_strain_asc( p )#, color = 'black', linewidth = 1.3, linestyle = '--', label = 'CAR-800-TU-bs4', xscale = 1000.  )

#        ex_path = join(simdb.exdata_dir, 'tensile_tests', 'dog_bone',
#                                         '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2F-V2.DAT')
#        ex_run = ExRun(ex_path)
#        ex_run.ex_type._plot_tex_stress_strain_asc( p )#, color = 'black', linewidth = 1.3, linestyle = '--', label = 'CAR-800-TU-bs4', xscale = 1000.  )

        ex_path = join(simdb.exdata_dir, 'tensile_tests', 'dog_bone',
                                         '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2F-V3.DAT')
        ex_run = ExRun(ex_path)
        ex_run.ex_type._plot_tex_stress_strain_asc( p )#, color = 'black', linewidth = 1.3, linestyle = '--', label = 'CAR-800-TU-bs4', xscale = 1000.  )
                                         
#        ex_path = join(simdb.exdata_dir, 'tensile_tests', 'dog_bone',
#                                         '2012-01-10_TT-12c-6cm-0-TU_SH1', 'TT-12c-6cm-0-TU-SH1F-V3.DAT')
#        ex_run = ExRun(ex_path)
#        ex_run.ex_type._plot_tex_stress_strain_asc( p )#, color = 'black', linewidth = 1.3, linestyle = '--', label = 'CAR-800-TU-bs4', xscale = 1000.  )

#        # AR-glas, tissue, 1200tex
#        #
#        path = join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs')
#        tests = ['TT-6g-2cm-0-V2.DAT']#'TT-6g-2cm-0-V1.DAT'
#        for t in tests:
#            ex_path = join(path, t)
#            ex_run = ExRun(ex_path)
##            ex_run.ex_type._plot_tex_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '-.', label = '2D-09-12'  )
#            ex_run.ex_type._plot_comp_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '-.', label = 'ARG-1200-TU', xscale = 1000.   )
#        
#        # AR-glas, tricot, 1200tex
#        #
##        path = join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-03-05_TT-10g-3cm-0-TR')
##        tests = ['TT-10g-3cm-0-TR-V2.DAT']#, 'TT-10g-3cm-0-TR-V3.DAT']
#        path = join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-03-17_TT-8g-3cm-0-TR')
#        tests = ['TT-8g-3cm-0-TR-V2.DAT']#, 'TT-8g-3cm-0-TR-V2.DAT', 'TT-8g-3cm-0-TR-V3.DAT']
#        for t in tests:
#            ex_path = join(path, t)
#            ex_run = ExRun(ex_path)
##            ex_run.ex_type._plot_tex_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '-', label = '2D-01-08'  )
#            ex_run.ex_type._plot_comp_stress_strain_asc( p, color = 'black', linewidth = 1.3, linestyle = '-', label = 'ARG-1200-TR', xscale = 1000.   )
#
##        format_plot(p, xlabel = 'Dehnung [1E-3]', ylabel = 'Textilspannung [MPa]', xlim = 14., ylim = 1400.)
#        format_plot(p, xlabel = 'Dehnung [1E-3]', ylabel = 'Kompositspannung [MPa]', xlim = 14., ylim = 25.)

        p.legend(prop = font)

        p.savefig('sig_c-eps', dpi=600.)
                    
        p.show()

