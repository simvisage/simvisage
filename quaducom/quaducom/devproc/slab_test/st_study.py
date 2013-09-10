'''
Created on Sep 9, 2013

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

from geo_st import \
    GeoST

from geo_lip import \
    GeoLIP
    
from geo_supprt import \
    GeoSUPPRT

from pickle import dump, load

from sim_st import SimST
from sim_st import SimSTDB

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

    test_series = 'ST-10g'
    test_series = 'ST-6c'
    test_series = 'ST-12c'
    
    #-----------------------------------------------
    # ST-10g: AG-glas slab tests (125 cm / 3 cm) with tricot-binding:
    #-----------------------------------------------
    if test_series == 'ST-10g':
     
        sim_model = SimSTDB( 
    
                            # calibration for: age = XXd; E_m = XX MPa; nu = XX   
                            #
                            ccs_unit_cell_key = 'FIL-10-09_2D-02-06a_0.00273_90_0',
                            calibration_test = 'TT11-10a-average',
                            thickness = 0.03,
                            length = 1.25,
                            age = 23,
                            #
                            elastmr_flag = False,
                            supprt_flag = False,
                            geo_st_flag = False,
                            #
                            shape_xy = 10,
                            shape_z = 3,
                            #
                            tstep = 0.05, 
                            tmax = 1.0, 
                            tolerance = 0.0005
                            )

    #-----------------------------------------------
    # ST-6c: carbon slab tests (80 cm / 2 cm): 
    #-----------------------------------------------
    if test_series == 'ST-10g':

        sim_model = SimSTDB(
                            
                            ccs_unit_cell_key = 'barrelshell_2D-05-11_0.00286_all0',
                            calibration_test = 'TTb-6c-2cm-0-TU-V1_bs5_a23d-nu02',
                            thickness = 0.02,
                            length = 0.80,
                            age = 23,
                            #
                            elastmr_flag = True,
                            supprt_flag = True,
                            geo_st_flag = True,
                            #
                            shape_xy = 28,
                            shape_R = 5,
                            shape_z = 3,
                            tstep = 0.05, 
                            tmax = 1.0, 
                            tolerance = 0.0005
                            )

    #-----------------------------------------
    # Var.3 ### ST-12c-6cm; L = 1,25m; t = 60 mm
    #-----------------------------------------
    if test_series == 'ST-10g':

        sim_model = SimSTDB( 
    
                            thickness = 0.06,
                            length = 1.25,
    
                            ccs_unit_cell_key = 'FIL-10-09_2D-05-11_0.00462_all0',
    
                            # calibration for: age = 23d; E_m = 28400 MPa; nu = 0.25   
                            #
                            calibration_test = 'TT-12c-6cm-0-TU-SH2F-V3_a23d-nu02',
                            # age of the slab at the time of testing
                            age = 23,
                            # NOTE: that the same phi-function is used independent of age. This assumes a 
                            # an afine/proportional damage evolution for different ages. 
                            #
                            elastmr_flag = True,
                            supprt_flag = True,
                            geo_st_flag = True,
    
                            shape_xy = 28,
                            shape_R = 5,
                            shape_z = 3,
                            #
                            tstep = 0.05, 
                            tmax = 1.0, 
                            tolerance = 0.0001
                            )

    #------------------------------
    # do
    #------------------------------
#    do = 'ui'
#    do = 'pstudy'
    do = 'validation'
#    do = 'show_last_results'

    #------------------------------
    # ui
    #------------------------------    
    if do == 'ui':
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource = sim_model)
        sim_model.tloop.eval()
        app.main()

    #------------------------------
    # validation
    #------------------------------
    if do == 'validation':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = 'pickle_files'
        png_path = 'png_files'

        p.figure(facecolor = 'white') # white background

        #--------------------        
        # simulation
        #--------------------        
        sim_model.tloop.eval()
 
        pickle_path = 'pickle_files'
        ccs_unit_cell_key = sim_model.ccs_unit_cell_key
        calibration_test = sim_model.calibration_test
        length = sim_model.length
        thickness = sim_model.thickness
        shape_xy = sim_model.shape_xy
        E_m = sim_model.E_m
        nu = sim_model.nu
        tolerance = sim_model.tolerance
        param_key = ccs_unit_cell_key + '_' + calibration_test + 'L=%g h=% sxy= %g Em=%g nu=%g tol=%g' %(length, thickness, shape_xy, E_m, nu, tolerance ) 

        # f-w-diagram_center
        #
        sim_model.f_w_diagram_center.refresh()
        file_name = 'f_w_diagram_c_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'w')
        dump(sim_model.f_w_diagram_center.trace, file)
        file.close()
        sim_model.f_w_diagram_center.trace.mpl_plot(p, color = 'red')

        # f-w-diagram_supprt
        #
        sim_model.f_w_diagram_supprt.refresh()
        file_name = 'f_w_diagram_supprt_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'w')
        dump(sim_model.f_w_diagram_supprt.trace, file)
        file.close()
        sim_model.f_w_diagram_supprt.trace.mpl_plot(p, color = 'blue')

#        # f-w-diagram_center-edge
#        #
#        sim_model.f_w_diagram_center.refresh()
#        file_name = 'f_w_diagram_ce_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'w')
#        dump(sim_model.f_w_diagram_center_edge.trace, file)
#        file.close()
#        sim_model.f_w_diagram_center_edge.trace.mpl_plot(p, color = 'red')
#
#        # f-w-diagram_edge
#        #
#        sim_model.f_w_diagram_center.refresh()
#        file_name = 'f_w_diagram_e_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'w')
#        dump(sim_model.f_w_diagram_edge.trace, file)
#        file.close()
#        sim_model.f_w_diagram_edge.trace.mpl_plot(p, color = 'red')

        #--------------------        
        # experiments
        #--------------------        

        # PT-12c-6cm-TU
#        path = join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU')
#        tests = [ 'ST-12c-6cm-u-TU.DAT' ]

        # PT-10a
#        path10 = join( simdb.exdata_dir, 'slab_tests', '2010-03-08_ST-10g-3cm-a-FR_TRC10', 'ST-10g-3cm-a-FR-TRC10.DAT' )
#        path11 = join( simdb.exdata_dir, 'slab_tests', '2010-03-09_ST-10g-3cm-a-FR_TRC11', 'ST-10g-3cm-a-FR-TRC11.DAT')
#        path12 = join( simdb.exdata_dir, 'slab_tests', '2010-03-10_ST-10g-3cm-a-FR_TRC12', 'ST-10g-3cm-a-FR-TRC12.DAT' )
#        tests = [path10, path11, path12]
#        for ex_path in tests:
#            ex_run = ExRun(ex_path)
#            ex_run.ex_type._plot_force_center_deflection_smoothed(p)

        # ST-6c-2cm-TU_bs2
        #
        ex_path = join( simdb.exdata_dir, 'slab_tests', '2013-07-10_ST-6c-2cm-TU_bs2', 'ST-6c-2cm-TU_bs2.DAT')
        ex_run = ExRun(ex_path)
        ex_run.ex_type._plot_force_center_deflection( p )

        # plot sim curve as time new roman within the predefined limits  
        #
#        format_plot(p, xlim = 34, ylim = 54, xlabel = 'displacement [mm]', ylabel = 'force [kN]')
        png_file_path = join(png_path, param_key + '.png')
        p.title( param_key )
        p.savefig( png_file_path, dpi=600. )
        p.show()

    #------------------------------
    # show last results
    #------------------------------

    if do == 'show_last_results':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        #------------------
        # simulation
        #------------------

        pickle_path = 'pickle_files'
#        param_key = 'SH2F-V3_nelems14-3-2_w_line_Ec28400'
#        param_key = 'SH2F-V3_nelems14-3-2_line_w_Em23d_a23d-nu02_ts0-0025'
#        param_key = 'SH2F-V3_nelems10-10-2_node_w_node_supprt'
        param_key = 'FIL-10-09_2D-02-06a_0.00273_90_0__TT11-10a-average'

        # f-w-diagram_supprt
        #
        file_name = 'f_w_diagram_supprt_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color = 'blue')

        # f-w-diagram_center
        #
        file_name = 'f_w_diagram_c_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color = 'red')

#        # f-w-diagram_center-edge
#        #
#        file_name = 'f_w_diagram_ce_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'r')
#        trace = load(file)
#        p.plot(trace.xdata, trace.ydata, color = 'red')
#
#        # f-w-diagram_edge
#        #
#        file_name = 'f_w_diagram_e_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'r')
#        trace = load(file)
#        p.plot(trace.xdata, trace.ydata, color = 'red')

        #------------------
        # experiments
        #------------------

#        path = join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU')
#        tests = [ 'ST-12c-6cm-u-TU.DAT' ]
##        path = join( simdb.exdata_dir, 'plate_tests', 'PT-10a' )
##        tests = [ 'PT10-10a.DAT', 'PT11-10a.DAT' , 'PT12-10a.DAT' ]
##        tests = [ 'PT10-10a.DAT' ]
#
#        for t in tests:
#            ex_path = join(path, t)
#            ex_run = ExRun(ex_path)
#            ex_run.ex_type._plot_force_deflection_avg_interpolated(p)
##            ex_run.ex_type._plot_force_deflection_avg( p )

        # ST-6c-2cm-TU_bs2
        #
        ex_path = join( simdb.exdata_dir, 'slab_tests', '2013-07-10_ST-6c-2cm-TU_bs2', 'ST-6c-2cm-TU_bs2.DAT')
        ex_run = ExRun(ex_path)
        ex_run.ex_type._plot_force_center_deflection( p )
#            ex_run.ex_type._plot_force_center_deflection_smoothed(p)

        # plot sim curve as time new roman within the predefined limits
        #
#        format_plot(p, xlim = 34, ylim = 54, xlabel = 'displacement [mm]', ylabel = 'force [kN]')

        p.show()


    if do == 'pstudy':
        sim_ps = SimPStudy(sim_model = sim_model)
        sim_ps.configure_traits()



