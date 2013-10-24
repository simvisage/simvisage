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

import os 

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

# from devproc.format_plot import format_plot

def format_plot(axes, xlim=None, ylim=None, xlabel='', ylabel=''):
    '''format 2d-plot black and with with times font 
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
    locs, labels = axes.xticks()
    axes.xticks(locs, map(lambda x: "%.0f" % x, locs), fontproperties=font)
    axes.xlabel(xlabel, fontproperties=font)

    locs, labels = axes.yticks()
    axes.yticks(locs, map(lambda x: "%.0f" % x, locs), fontproperties=font)
    axes.ylabel(ylabel, fontproperties=font)
    
#------------------------------------------------
# script for parameter study  
#------------------------------------------------

if __name__ == '__main__':

    #------------------------------
    # do
    #------------------------------
#    do = 'show_phi_fn'
#    do = 'ui'
    do = 'validation'
#    do = 'show_last_results'
    
#    test_series = 'ST-10g'
    test_series = 'ST-12c'
#    test_series = 'ST-6c'
    
    #-----------------------------------------------
    # ST-10g: AG-glas slab tests (125 cm / 3 cm) with tricot-binding:
    #-----------------------------------------------
    if test_series == 'ST-10g':
     
        sim_model = SimSTDB(
    
                            # calibration for: age = XXd; E_m = XX MPa; nu = XX   
                            #
                            ccs_unit_cell_key='FIL-10-09_2D-02-06a_0.00273_90_0',
                            calibration_test='TT11-10a-average',
#                            age = 23,
                            #
                            thickness=0.03,
                            length=1.25,
                            #
                            elstmr_flag=False,
                            supprt_flag=False,
                            geo_st_flag=False,
                            #
                            shape_xy=10,
                            shape_z=2,
                            #
                            tstep=0.05,
                            tmax=1.00,
                            tolerance=0.0005,
                            ord=np.inf
                            )

    #-----------------------------------------
    # ST-12c-6cm; L = 1,25m; t = 6 cm
    #-----------------------------------------
    if test_series == 'ST-12c':

        sim_model = SimSTDB(
    
                            thickness=0.06,
                            length=1.25,
                            radius_plate=0.095,  # D=8cm
    
                            ccs_unit_cell_key='FIL-10-09_2D-05-11_0.00462_all0',
    
                            # calibration for: age = 23d; E_m = 27975.8 MPa; nu = 0.20; nsteps = 100   
                            #
                            calibration_test='TT-12c-6cm-0-TU-V1_ZiE-S1_age23_Em27975.8_nu0.2_nsteps100',
                            
                            # calibration for: age = 23d; E_m = 27975 MPa; nu = 0.20; nsteps = 100   
                            #
#                             calibration_test='TT-12c-6cm-0-TU-SH2F-V3_a23d_nu02_s100',
#                            calibration_test = 'TT-12c-6cm-0-TU-SH2F-V3_a23d_nu02_s50',
#                            calibration_test = 'TT-12c-6cm-TU-SH1F-V1',
                            
                            n_mp=30,
                            
                            # age of the slab at the time of testing
                            age=23,
                            # NOTE: that the same phi-function is used independent of age. This assumes a 
                            # an afine/proportional damage evolution for different ages. 
                            #
                            elstmr_flag=False,
                            supprt_flag=True,
                            geo_st_flag=True,
                            #
                            # coarse mesh:
                            shape_xy=10,
                            shape_R=2,
                            shape_z=3,
                            shape_supprt_xy=2,
                            #
                            # fine mesh:
#                            shape_z = 2,
#                            shape_xy = 26,
#                            shape_R = 4,
#                            shape_supprt_xy = 4,
                            #
                            w_max=-0.040,
                            tstep=0.02,
#                            tstep = 1.00, 
                            tmax=1.00,
                            # 'NOTE: tloop.norm switched to "max(abs(x))"'
                            tolerance=0.0001,  # #[MN]0.0001#1e-6#1e-8#0.0005
                            ord=np.inf,  # "norm = max(abs(x_i))"
                            #
                            # 'factor_eps_fail' = 1.0 (default)
                            phi_fn_class=PhiFnGeneralExtended
                            )

    #-----------------------------------------------
    # ST-6c: carbon slab tests (80 cm / 2 cm): 
    #-----------------------------------------------
    if test_series == 'ST-6c':

        sim_model = SimSTDB(
                            
                            ccs_unit_cell_key='barrelshell_2D-05-11_0.00286_all0',
                            calibration_test='TTb-6c-2cm-0-TU-V1_bs5_nu02_s100',
                            age=28,
                            #
                            thickness=0.02,
                            length=0.80,
                            radius_plate=0.04,  # D=8cm
                            #
                            elstmr_flag=False,
                            supprt_flag=True,
                            geo_st_flag=True,
                            #
                            shape_xy=10,
                            shape_R=2,
                            shape_z=3,
                            #
                            w_max=-0.040,
                            tstep=0.01,
                            tmax=1.0,
                            # 'NOTE: tloop.norm switched to "max(abs(x))"'
                            tolerance=0.0001,  # #[MN]0.0001#1e-6#1e-8#0.0005
                            #
                            # 'factor_eps_fail' = 1.0 (default)
                            phi_fn_class=PhiFnGeneralExtended
                            )

    # print settings:
    #
    ccs_unit_cell_key = sim_model.ccs_unit_cell_key
    calibration_test = sim_model.calibration_test
    length = sim_model.length
    thickness = sim_model.thickness
    shape_xy = sim_model.shape_xy
    E_m = sim_model.E_m
    nu = sim_model.nu
    tolerance = sim_model.tolerance
    n_mp = sim_model.n_mp

    print '\n' 
    print '### calculation settings: ###'
    print 'ccs_unit_cell_key', ccs_unit_cell_key
    print 'calibration_test', calibration_test
    print 'length', length
    print 'thickness', thickness
    print 'shape_xy', shape_xy
    print 'E_m', E_m
    print 'nu', nu
    print 'tolerance', tolerance
    print 'n_mp', n_mp
    print '\n' 

#--------------------------------------------------------------
# do: ui / validation / show_last_result / pstudy
#--------------------------------------------------------------
     
    if do == 'show_phi_fn':
        import pylab as p
        p.figure(facecolor='white') 

        phi_fn = sim_model.phi_fn
#        phi_fn_ext = sim_model.phi_fn_ext
#        phi_fn_exp = sim_model.phi_fn_exp

#        damage_function = sim_model.damage_function
#        damage_function.plot(p, color = 'red', linewidth = 1)
#        print 'sim_model.damage_function', sim_model.damage_function
#        print 'self.ccs_unit_cell_ref.damage_function_list', [sim_model.ccs_unit_cell_ref.damage_function_list[i].calibration_test for i in range(len(sim_model.ccs_unit_cell_ref.damage_function_list))]
        
        phi_fn.mfn.plot(p, color='black', linewidth=3)
#        phi_fn_ext.mfn.plot(p, color = 'blue', linewidth = 2)
#        phi_fn_exp.mfn.plot(p, color = 'green', linewidth = 3)

        xmax = sim_model.damage_function.xdata[-1]
        print 'xmax', xmax
        x = linspace(0, 3 * xmax, 1000)
        phi_fn = frompyfunc(phi_fn, 1, 1)
        y = phi_fn(x)
        p.plot(x, y, color='grey', linewidth=2)

        p.show()

    #------------------------------
    # ui
    #------------------------------    
    if do == 'ui':
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource=sim_model)
        sim_model.tloop.eval()
        app.main()

    #------------------------------
    # validation
    #------------------------------
    if do == 'validation':
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource=sim_model)
        
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = join(simdb.simdata_dir, 'pickle_files')
        png_path = join(simdb.simdata_dir, 'png_files')

        if not os.path.exists(pickle_path):
            os.mkdir(pickle_path)
            os.mkdir(png_path)

        # pstudy: thickness
        #
        pst_list = [ 0.06 ]  # 0.0554, 
        
        # pstudy: n_mp
        #
#        pst_list = [ 30 ]
        
        # pstudy: calibration test
        #
#        pst_list = [ 'TT-12c-6cm-0-TU-SH2F-V3_a23d_nu02_s100' , 'TT-12c-6cm-TU-SH1F-V1' ]

        # pstudy: phi_fn
        #
#        pst_list = [ PhiFnGeneral, PhiFnGeneralExtended, PhiFnGeneralExtendedExp ] 

        
        for pst_param in pst_list:
            
            sim_model.thickness = pst_param
            
#            sim_model.n_mp = pst_param
#            sim_model.calibration_test = pst_param
#            sim_model.phi_fn_class = pst_param
            
            p.figure(facecolor='white')  # white background for diagram
    
            #--------------------        
            # simulation 
            #--------------------        
            sim_model.tloop.eval()
     
            # settings:
            #
            sim_model_name = sim_model.__class__.__name__
            ccs_unit_cell_key = sim_model.ccs_unit_cell_key
            calibration_test = sim_model.calibration_test
            length = sim_model.length
            thickness = sim_model.thickness
            print 'thickness', thickness
            radius_plate = sim_model.radius_plate
            shape_xy = sim_model.shape_xy
            shape_z = sim_model.shape_z
            shape_R = sim_model.shape_R
            E_m = sim_model.E_m
            nu = sim_model.nu
            tolerance = sim_model.tolerance
            phi_fn_class = sim_model.phi_fn_class.__name__
            supprt_flag = str(sim_model.supprt_flag)
            geo_st_flag = str(sim_model.geo_st_flag) 
            n_mp = sim_model.n_mp 
            tstep = sim_model.tstep
            w_max = sim_model.w_max

            # param_key 
            #
            param_key = sim_model_name + '_' + ccs_unit_cell_key + '_' + calibration_test + '_%s_L%gh%gR%g_sxy%gz%gR%g_s%sg%s_Em%g_nu%g_tol%g_w%g_ts%g_nmp%g' \
                        % (phi_fn_class, length, thickness, radius_plate, shape_xy, shape_z, shape_R, supprt_flag[0], geo_st_flag[0], E_m, nu, tolerance, w_max, tstep, n_mp) 
            print 'param_key = %s' % param_key
    
            # f-w-diagram_center
            #
            sim_model.f_w_diagram_center.refresh()
            file_name = 'f_w_diagram_c_' + param_key + '.pickle'
            pickle_file_path = join(pickle_path, file_name)
            file = open(pickle_file_path, 'w')
            dump(sim_model.f_w_diagram_center.trace, file)
            print 'pickle file saved to file: %s' % file_name
            file.close()
            sim_model.f_w_diagram_center.trace.mpl_plot(p, color='red')
    
            # f-w-diagram_supprt
            #
            sim_model.f_w_diagram_supprt.refresh()
            file_name = 'f_w_diagram_supprt_' + param_key + '.pickle'
            pickle_file_path = join(pickle_path, file_name)
            file = open(pickle_file_path, 'w')
            dump(sim_model.f_w_diagram_supprt.trace, file)
            print 'pickle file saved to file: %s' % file_name
            file.close()
            sim_model.f_w_diagram_supprt.trace.mpl_plot(p, color='blue')
    
            #--------------------        
            # experiments
            #--------------------        
    
            if test_series == 'ST-12c':
                # PT-12c-6cm-TU
                #
                ex_path = join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU', 'ST-12c-6cm-u-TU.DAT')
                ex_run = ExRun(ex_path)
                ex_run.ex_type._plot_force_center_deflection_interpolated(p)
    
            if test_series == 'ST-10g':
                # PT-10a
                #
                ex_path_TRC10 = join(simdb.exdata_dir, 'slab_tests', '2010-03-08_ST-10g-3cm-a-FR_TRC10', 'ST-10g-3cm-a-FR-TRC10.DAT')
                ex_path_TRC11 = join(simdb.exdata_dir, 'slab_tests', '2010-03-09_ST-10g-3cm-a-FR_TRC11', 'ST-10g-3cm-a-FR-TRC11.DAT')
                ex_path_TRC12 = join(simdb.exdata_dir, 'slab_tests', '2010-03-10_ST-10g-3cm-a-FR_TRC12', 'ST-10g-3cm-a-FR-TRC12.DAT')
                tests = [ex_path_TRC10, ex_path_TRC11, ex_path_TRC12]
                for ex_path in tests:
                    ex_run = ExRun(ex_path)
                    ex_run.ex_type._plot_force_center_deflection_interpolated(p)
    
            if test_series == 'ST-6c':
                # ST-6c-2cm-TU_bs2
                #
                ex_path = join(simdb.exdata_dir, 'slab_tests', '2013-07-10_ST-6c-2cm-TU_bs2', 'ST-6c-2cm-TU_bs2.DAT')
                ex_run = ExRun(ex_path)
                ex_run.ex_type._plot_force_center_deflection(p, offset_w=0.03)
    
            #----------------------------------------------------------------------
            # plot sim curve as time new roman within the predefined limits  
            #----------------------------------------------------------------------
            #
    #        format_plot(p, xlim = 34, ylim = 54, xlabel = 'displacement [mm]', ylabel = 'force [kN]')
            png_file_path = join(png_path, param_key + '.png')
            p.title(param_key, fontsize=8)
            p.savefig(png_file_path, dpi=600.)
            print 'png-file saved to file: %s' % png_file_path
            p.show()

        app.main()


    #------------------------------
    # show last results
    #------------------------------

    if do == 'show_last_results':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = join(simdb.simdata_dir, 'pickle_files')
        png_path = join(simdb.simdata_dir, 'png_files')

        # settings:
        #
        sim_model_name = sim_model.__class__.__name__
        ccs_unit_cell_key = sim_model.ccs_unit_cell_key
        calibration_test = sim_model.calibration_test
        length = sim_model.length
        thickness = sim_model.thickness
        radius_plate = sim_model.radius_plate
        shape_xy = sim_model.shape_xy
        shape_z = sim_model.shape_z
        shape_R = sim_model.shape_R
        E_m = sim_model.E_m
        nu = sim_model.nu
        tolerance = sim_model.tolerance
        phi_fn_class = sim_model.phi_fn_class.__name__
        supprt_flag = str(sim_model.supprt_flag)
        geo_st_flag = str(sim_model.geo_st_flag) 
        n_mp = sim_model.n_mp 
        tstep = sim_model.tstep
        w_max = sim_model.w_max

        # param_key 
        #
        param_key = sim_model_name + '_' + ccs_unit_cell_key + '_' + calibration_test + '_%s_L%gh%gR%g_sxy%gz%gR%g_s%sg%s_Em%g_nu%g_tol%g_w%g_ts%g_nmp%g' \
                    % (phi_fn_class, length, thickness, radius_plate, shape_xy, shape_z, shape_R, supprt_flag[0], geo_st_flag[0], E_m, nu, tolerance, w_max, tstep, n_mp) 
        print 'param_key = %s' % param_key
        
        #------------------
        # simulation
        #------------------

        # f-w-diagram_supprt
        #
        file_name = 'f_w_diagram_supprt_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color='blue')

#        # f-w-diagram_center
#        #
#        file_name = 'f_w_diagram_c_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'r')
#        trace = load(file)
#        p.plot(trace.xdata, trace.ydata, color='red')

        #--------------------        
        # experiments
        #--------------------        

        if test_series == 'ST-12c':
            # PT-12c-6cm-TU
            #
            ex_path = join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU', 'ST-12c-6cm-u-TU.DAT')
            ex_run = ExRun(ex_path)
            ex_run.ex_type._plot_force_center_deflection(p)

        if test_series == 'ST-10g':
            # PT-10a
            #
            ex_path_TRC10 = join(simdb.exdata_dir, 'slab_tests', '2010-03-08_ST-10g-3cm-a-FR_TRC10', 'ST-10g-3cm-a-FR-TRC10.DAT')
            ex_path_TRC11 = join(simdb.exdata_dir, 'slab_tests', '2010-03-09_ST-10g-3cm-a-FR_TRC11', 'ST-10g-3cm-a-FR-TRC11.DAT')
            ex_path_TRC12 = join(simdb.exdata_dir, 'slab_tests', '2010-03-10_ST-10g-3cm-a-FR_TRC12', 'ST-10g-3cm-a-FR-TRC12.DAT')
            tests = [ex_path_TRC10, ex_path_TRC11, ex_path_TRC12]
            for ex_path in tests:
                ex_run = ExRun(ex_path)
                ex_run.ex_type._plot_force_center_deflection(p)

        if test_series == 'ST-6c':
            # ST-6c-2cm-TU_bs2
            #
            ex_path = join(simdb.exdata_dir, 'slab_tests', '2013-07-10_ST-6c-2cm-TU_bs2', 'ST-6c-2cm-TU_bs2.DAT')
            ex_run = ExRun(ex_path)
            # specify offset as displacement measurement starts a little to late at a small force > 0 
            # shift experimental curve to simulation curve value for the displacement at this load level 
            # corresponding to linear elastic approximation of slab stiffness; 
            # 
            ex_run.ex_type._plot_force_center_deflection(p, offset_w=0.0385)  
        # plot sim curve as time new roman within the predefined limits
        #
#        format_plot(p, xlim = 34, ylim = 54, xlabel = 'displacement [mm]', ylabel = 'force [kN]')

        p.show()
