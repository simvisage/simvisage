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

from pickle import dump, load

from sim_bt_4pt import SimBT4PT
from sim_bt_4pt import SimBT4PTDB

# from devproc.format_plot import format_plot

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
    
    test_series = 'BT-12c'
    
    #-----------------------------------------
    # BT-4PT-12c-6cm; L_0 = 1.50m; t = 6 cm
    #-----------------------------------------
    if test_series == 'BT-12c':

        sim_model = SimBT4PTDB(
                               ccs_unit_cell_key='FIL-10-09_2D-05-11_0.00462_all0',
                               calibration_test='TT-12c-6cm-0-TU-SH2F-V3_a23d_nu02_s100',
                               age=26,
                               #
                               thickness=0.06,
                               length=1.50,
                               width=0.20,
                               #
                               elstmr_flag=False,
                               supprt_flag=False,
                               #
                               # fine mesh
#                               outer_zone_shape_x = 10,
#                               load_zone_shape_x = 2,
#                               mid_zone_shape_x = 8,
#                               shape_y = 2,
#                               shape_z = 3,
                               #
                               # coarse mesh
                               outer_zone_shape_x=6,
                               load_zone_shape_x=2,
                               mid_zone_shape_x=4,
                               shape_y=2,
                               shape_z=3,
                               #
                               w_max= -0.030,
                               tstep=0.04,
                               tmax=1.00,
                               tolerance=0.0001,
                               ord=np.inf,  # "norm = max(abs(x_i))"
                               #
                               # 'factor_eps_fail' = 1.0 (default)
                               phi_fn_class=PhiFnGeneralExtended
                               )

    # print settings:
    #
    ccs_unit_cell_key = sim_model.ccs_unit_cell_key
    calibration_test = sim_model.calibration_test
    length = sim_model.length
    width = sim_model.width
    thickness = sim_model.thickness
    outer_zone_shape_x = sim_model.outer_zone_shape_x
    load_zone_shape_x = sim_model.load_zone_shape_x
    mid_zone_shape_x = sim_model.mid_zone_shape_x
    shape_y = sim_model.shape_y
    shape_z = sim_model.shape_z
    E_m = sim_model.E_m
    nu = sim_model.nu
    tolerance = sim_model.tolerance
    n_mp = sim_model.n_mp

    print '\n' 
    print '### calculation settings: ###'
    print 'ccs_unit_cell_key', ccs_unit_cell_key
    print 'calibration_test', calibration_test
    print 'length', length
    print 'width', width
    print 'thickness', thickness
    print 'outer_zone_shape_x', outer_zone_shape_x 
    print 'load_zone_shape_x', load_zone_shape_x,
    print 'mid_zone_shape_x', mid_zone_shape_x
    print 'shape_y', shape_y
    print 'shape_z', shape_z
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

        # plot mfn-line function stored in phi_fn
        # 
        phi_fn = sim_model.phi_fn
        phi_fn.mfn.plot(p, color='black', linewidth=3)

        # plot phi_fn (extended plot range)
        #
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
#        sim_model.tloop.eval()
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

        # pstudy: n_mp
        #
        pst_list = [ 0.060, 0.0554, 0.065, 0.05 ]
        
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
#            sim_model.calibration_test = st_param
#            sim_model.phi_fn_class = st_param
            
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
            outer_zone_shape_x = sim_model.outer_zone_shape_x
            load_zone_shape_x = sim_model.load_zone_shape_x
            mid_zone_shape_x = sim_model.mid_zone_shape_x
            shape_y = sim_model.shape_y
            shape_z = sim_model.shape_z
            E_m = sim_model.E_m
            nu = sim_model.nu
            tolerance = sim_model.tolerance
            phi_fn_class = sim_model.phi_fn_class.__name__
            supprt_flag = str(sim_model.supprt_flag)
            elstmr_flag = str(sim_model.elstmr_flag) 
            n_mp = sim_model.n_mp 
            tstep = sim_model.tstep
            w_max = sim_model.w_max

            # param_key 
            #
            param_key = sim_model_name + '_' + ccs_unit_cell_key + '_' + calibration_test + '_%s_L%g_h%g_sxo%gl%gm%gy%gz%g_s%se%s_Em%g_nu%g_tol%g_w%g_ts%g_nmp%g' \
                        % (phi_fn_class, length, thickness, outer_zone_shape_x, load_zone_shape_x, mid_zone_shape_x, shape_y, shape_z, supprt_flag[0], elstmr_flag[0], E_m, nu, tolerance, w_max, tstep, n_mp) 
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
    
            # f-w-diagram_thirdpoint
            #
            sim_model.f_w_diagram_supprt.refresh()
            file_name = 'f_w_diagram_thirdpoint_' + param_key + '.pickle'
            pickle_file_path = join(pickle_path, file_name)
            file = open(pickle_file_path, 'w')
            dump(sim_model.f_w_diagram_supprt.trace, file)
            print 'pickle file saved to file: %s' % file_name
            file.close()
            sim_model.f_w_diagram_supprt.trace.mpl_plot(p, color='blue')
    
            #--------------------        
            # experiments
            #--------------------        
    
            if test_series == 'BT-12c':
                # continuous fabrication
                #
                path = join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU', 'BT-4PT-12c-6cm-SH4')
                tests = [ 'BT-4PT-12c-6cm-SH4-V1.DAT']  # , 'BT-4PT-12c-6cm-SH4-V2.DAT' ]
                for t in tests:
                    ex_path = join(path, t)
                    ex_run = ExRun(ex_path)
#                    ex_run.ex_type._plot_ironed_force_deflection_avg(p)
                    ex_run.ex_type._plot_force_deflection_center(p, offset_w= -0.017)
#                    ex_run.ex_type._plot_force_deflection_thirdpoints( p )

                # continuous fabrication with pause
                #
#                path = join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU', 'BT-4PT-12c-6cm-SH4FF')
#                tests = [ 'BT-4PT-12c-6cm-SH4FF-V1.DAT', 'BT-4PT-12c-6cm-SH4FF-V2.DAT' ]
#                for t in tests:
#                    ex_path = join(path, t)
#                    ex_run = ExRun(ex_path)
#                    ex_run.ex_type._plot_force_deflection_center(p)
    
            #----------------------------------------------------------------------
            # plot sim curve as time new roman within the predefined limits  
            #----------------------------------------------------------------------
            #
    #        format_plot(p, xlim = 34, ylim = 54, xlabel = 'displacement [mm]', ylabel = 'force [kN]')
            png_file_path = join(png_path, param_key + '.png')
            p.title(param_key, fontsize=8)
            p.savefig(png_file_path, dpi=1200.)
            print 'png-file saved to file: %s' % png_file_path
#            p.show()

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
        outer_zone_shape_x = sim_model.outer_zone_shape_x
        load_zone_shape_x = sim_model.load_zone_shape_x
        mid_zone_shape_x = sim_model.mid_zone_shape_x
        shape_y = sim_model.shape_y
        shape_z = sim_model.shape_z
        E_m = sim_model.E_m
        nu = sim_model.nu
        tolerance = sim_model.tolerance
        phi_fn_class = sim_model.phi_fn_class.__name__
        supprt_flag = str(sim_model.supprt_flag)
        elstmr_flag = str(sim_model.elstmr_flag) 
        n_mp = sim_model.n_mp 
        tstep = sim_model.tstep
        w_max = sim_model.w_max

        # param_key 
        #
        param_key = sim_model_name + '_' + ccs_unit_cell_key + '_' + calibration_test + '_%s_L%g_h%g_sxo%gl%gm%gy%gz%g_s%se%s_Em%g_nu%g_tol%g_w%g_ts%g_nmp%g' \
                    % (phi_fn_class, length, thickness, outer_zone_shape_x, load_zone_shape_x, mid_zone_shape_x, shape_y, shape_z, supprt_flag[0], elstmr_flag[0], E_m, nu, tolerance, w_max, tstep, n_mp) 
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

        # f-w-diagram_center
        #
        file_name = 'f_w_diagram_c_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color='red')

        # f-w-diagram_thirdpoint
        #
        file_name = 'f_w_diagram__thirdpoint_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color='red')

        #--------------------        
        # experiments
        #--------------------        

        if test_series == 'BT-12c':
            path = join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU', 'BT-4PT-12c-6cm-SH4')
            tests = [ 'BT-4PT-12c-6cm-SH4-V1.DAT' ]  # , 'BT-4PT-12c-6cm-SH4-V2.DAT' ]
            for t in tests:
                ex_path = join(path, t)
                ex_run = ExRun(ex_path)
                ex_run.ex_type._plot_force_deflection_center(p, offset_w= -0.017)
#                ex_run.ex_type._plot_ironed_force_deflection_avg(p)
#                ex_run.ex_type._plot_force_deflection_thirdpoints( p )
            path = join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU', 'BT-4PT-12c-6cm-SH4FF')
            tests = [ 'BT-4PT-12c-6cm-SH4FF-V1.DAT', 'BT-4PT-12c-6cm-SH4FF-V2.DAT' ]
            for t in tests:
                ex_path = join(path, t)
                ex_run = ExRun(ex_path)
                ex_run.ex_type._plot_force_deflection_center(p)


        # plot sim curve as time new roman within the predefined limits
        #
#        format_plot(p, xlim = 34, ylim = 54, xlabel = 'displacement [mm]', ylabel = 'force [kN]')

        p.show()

    if do == 'pstudy':
        sim_ps = SimPStudy(sim_model=sim_model)
        sim_ps.configure_traits()


