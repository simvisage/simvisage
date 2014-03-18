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

from sim_bt_3pt import SimBT3PT
from sim_bt_3pt import SimBT3PTDB

from quaducom.devproc.format_plot import format_plot

if __name__ == '__main__':

    #------------------------------
    # do
    #------------------------------
#    do = 'show_phi_fn'
#    do = 'ui'
    do = 'validation'
#     do = 'show_last_results'

#    test_series = 'BT-12c'
     test_series = 'BT-6c'

    #-----------------------------------------
    # BT-3PT-12c-6cm; L = 1.25m; L0 = 1.15 cm t = 6 cm
    #-----------------------------------------
    if test_series == 'BT-12c':

        sim_model = SimBT3PTDB(
                               ccs_unit_cell_key='FIL-10-09_2D-05-11_0.00462_all0',
#                               calibration_test = 'TT-12c-6cm-0-TU-SH2F-V3_a23d_nu02_s100',
#                               calibration_test='TT-12c-6cm-0-TU-V1_a9d_nu02_s100',

                               # ZiE-S1 (TT-age = 11d)
                               #
#                                calibration_test='TT-12c-6cm-0-TU-V2_ZiE-S1_age9_Em24524.1_nu0.2_nsteps100',
#                                 calibration_test='TT-12c-6cm-0-TU-V2_ZiE-S1_age9_Ec26080.6_nu0.2_nsteps100',

                               # ZiE-S2 (TT-age = 9d)
                               #
#                               calibration_test='TT-12c-6cm-0-TU-V2_ZiE-S2_age9_Em24524.1_nu0.2_nsteps100',
                               calibration_test='TT-12c-6cm-0-TU-V2_ZiE-S2_age9_Ec26080.6_nu0.2_nsteps100',

                               age=9,
                               #
#                                 thickness=0.06,
                                thickness=0.0558,
#                                 thickness=0.0554,

                               length=1.15,
                               width=0.20,
                               #
                               elstmr_flag=False,
                               supprt_flag=False,
                               #
                               shape_x=4,
                               mid_shape_x=1,
                               shape_y=1,
                               shape_z=3,
                               #
                               w_max= -0.020,
                               tstep=0.02,
                               tmax=1.00,
                               tolerance=0.0001,
                               #
                               # 'factor_eps_fail' = 1.0 (default)
                               phi_fn_class=PhiFnGeneralExtended
                               )

    #-----------------------------------------------
    # ST-6c: carbon slab tests (80 cm / 2 cm): 
    #-----------------------------------------------
    if test_series == 'BT-6c':

        sim_model = SimBT3PTDB(
                               # calibration performed with:
                               # 'n_mp' = 30; 
                               # 'E_c' = 28600; 
                               # 'nu' = 0.25; 
                               # 'compliance_version'
                               ccs_unit_cell_key='barrelshell_2D-05-11_0.00286_all0',

#                                calibration_test='TTb-6c-2cm-0-TU-V1_bs5_nu02_s100',

                               # Em = 19800. MPa
                               #
#                               calibration_test='TT-6c-2cm-0-TU-V3_bs1_age28_Em19800_nu0.2_nsteps100',
#                               calibration_test='TT-6c-2cm-0-TU-V1_bs2_age28_Em19800_nu0.2_nsteps100',
#                               calibration_test='TT-6c-2cm-0-TU-V1_bs3_age28_Em19800_nu0.2_nsteps100',
#                               calibration_test='TTb-6c-2cm-0-TU-V3_bs5_age28_Em19800_nu0.2_nsteps100',

                               # Ec = 22213.2 MPa
                               #
#                               calibration_test='TT-6c-2cm-0-TU-V3_bs1_age28_Ec22213.2_nu0.2_nsteps100',
#                               calibration_test='TT-6c-2cm-0-TU-V1_bs2_age28_Ec22213.2_nu0.2_nsteps100',
#                               calibration_test='TT-6c-2cm-0-TU-V1_bs3_age28_Ec22213.2_nu0.2_nsteps100',
                                calibration_test='TTb-6c-2cm-0-TU-V3_bs5_age28_Ec22213.2_nu0.2_nsteps100',

                               age=28,
                               #
#                                 thickness=0.020,
#                                 thickness=0.0186,
                               thickness=0.02,
                               length=0.46,
                               width=0.10,
                               #
                               elstmr_flag=False,
                               supprt_flag=False,
                               #
                               shape_x=4,
                               mid_shape_x=1,
                               shape_y=1,
                               shape_z=3,
                               #
                               w_max= -0.015,
                               tstep=0.02,
                               tmax=1.00,
                               tolerance=0.0001, # [MN]
                               ord=np.inf, # "norm = max(abs(x_i))"
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
    shape_x = sim_model.shape_x
    E = sim_model.specmn_mats.E
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
    print 'shape_x', shape_x
    print 'E_m', E_m
    print 'E_mats =', E
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
#         app = IBVPyApp(ibv_resource=sim_model)

        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = join(simdb.simdata_dir, 'pickle_files')
        png_path = join(simdb.simdata_dir, 'png_files')

        if not os.path.exists(pickle_path):
            os.mkdir(pickle_path)
            os.mkdir(png_path)

        pst_list = [1]
        # pstudy: thickness
        #
#         pst_list = [ 0.060 ]
#         pst_list = [ 0.020, 0.019, 0.018, 0.017, 0.021, 0.022, 0.023 ]
#         pst_list = [ 0.020, 0.0186, 0.0171 ]

        # pstudy: shape_z
        #
#        pst_list = [ 3, 6, 8 ]

        # pstudy: n_mp
        #
#        pst_list = [ 6, 12, 30 ]

        # pstudy: calibration test
        #
#         pst_list = [ 
#                     'TT-12c-6cm-0-TU-V2_ZiE-S1_age9_Ec26080.6_nu0.2_s100' ,
#                     'TT-12c-6cm-0-TU-V2_ZiE-S2_age9_Ec26080.6_nu0.2_s100' ,
#                    ]

        # pstudy: phi_fn
        #
#        st_list = [ PhiFnGeneral, PhiFnGeneralExtended, PhiFnGeneralExtendedExp ] 


        for pst_param in pst_list:

#            sim_model.thickness = pst_param
#            sim_model.shape_z = pst_param
#            sim_model.n_mp = pst_param
#            sim_model.calibration_test = pst_param
#            sim_model.phi_fn_class = pst_param

            p.figure(facecolor='white', figsize=(12, 9))  # white background for diagram

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
            shape_x = sim_model.shape_x
            mid_shape_x = sim_model.mid_shape_x
            shape_y = sim_model.shape_y
            shape_z = sim_model.shape_z
            E_m = sim_model.E_m
            # E-modulus used by 'specmn_mats'
            #
            E = sim_model.specmn_mats.E
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
            param_key = sim_model_name + '_' + ccs_unit_cell_key + '_' + calibration_test + '_%s_L%g_h%g_sx%gm%gy%gz%g_s%se%s_E%g_nu%g_tol%g_w%g_ts%g_nmp%g' \
                        % (phi_fn_class, length, thickness, shape_x, mid_shape_x, shape_y, shape_z, supprt_flag[0], elstmr_flag[0], E, nu, tolerance, w_max, tstep, n_mp)
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

            if test_series == 'BT-12c':
                ex_path_V1 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE',
                                       'BT-3PT-12c-6cm-0-Tu-V1.raw')
                ex_path_V2 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE',
                                       'BT-3PT-12c-6cm-0-Tu-V2.raw')
                ex_path_V3 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE',
                                       'BT-3PT-12c-6cm-0-Tu-V3.raw')
                ex_path_V4 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE',
                                       'BT-3PT-12c-6cm-0-Tu-V4.raw')
                tests = [ex_path_V1, ex_path_V2, ex_path_V3, ex_path_V4]
                for ex_path in tests:
                    ex_run = ExRun(ex_path)
                    # ex_run.ex_type._plot_force_machine_displacement_wo_elast_interpolated(p, color='black')
                    ex_run.ex_type._plot_force_machine_displacement_wo_elast_analytical_offset(p, color='black')

                # format plot font as time new roman within the predefined limits  
                #
                format_plot(p, xlim=25, ylim=14, xlabel='displacement [mm]', ylabel='force [kN]')


            if test_series == 'BT-6c':
                ex_path_V2 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-02_BT-6c-2cm-0-TU_bs4',
                                       'BT-6c-2cm-0-TU-V2_bs4.DAT')
                ex_path_V3 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-02_BT-6c-2cm-0-TU_bs4',
                                       'BT-6c-2cm-0-TU-V3_bs4.DAT')
#                 ex_path_V1 = join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-18_BT-6c-2cm-0-TU_bs5',
#                                        'BT-6c-2cm-0-TU-V1_bs5.DAT')
#                 ex_path_V2 = join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-18_BT-6c-2cm-0-TU_bs5',
#                                        'BT-6c-2cm-0-TU-V2_bs5.DAT')
#                 ex_path_V3 = join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-18_BT-6c-2cm-0-TU_bs5',
#                                        'BT-6c-2cm-0-TU-V3_bs5.DAT')
                tests = [ ex_path_V2, ex_path_V3 ]
                for ex_path in tests:
                    ex_run = ExRun(ex_path)
                    ex_run.ex_type._plot_force_gauge_displacement_with_analytical_offset(p, color='black')

                # plot sim curve as time new roman within the predefined limits
                #
                format_plot(p, xlim=20, ylim=3.5, xlabel='displacement [mm]', ylabel='force [kN]')

            #----------------------------------
            # show plot and save png-file
            #----------------------------------
            #
            png_file_path = join(png_path, param_key + '.png')
            p.title(param_key, fontsize=8)
            p.savefig(png_file_path, dpi=300.)
            print 'png-file saved to file: %s' % png_file_path
            p.show()

#         app.main()


    #------------------------------
    # show last results
    #------------------------------

    if do == 'show_last_results':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = join(simdb.simdata_dir, 'pickle_files')

        # settings:
        #
        sim_model_name = sim_model.__class__.__name__
        ccs_unit_cell_key = sim_model.ccs_unit_cell_key
        calibration_test = sim_model.calibration_test
        length = sim_model.length
        thickness = sim_model.thickness
        shape_x = sim_model.shape_x
        mid_shape_x = sim_model.mid_shape_x
        shape_y = sim_model.shape_y
        shape_z = sim_model.shape_z
        E_m = sim_model.E_m
        # E-modulus used by 'specmn_mats'
        #
        E = sim_model.specmn_mats.E
        nu = sim_model.nu
        tolerance = sim_model.tolerance
        phi_fn_class = sim_model.phi_fn_class.__name__
        supprt_flag = str(sim_model.supprt_flag)
        elstmr_flag = str(sim_model.elstmr_flag)
        n_mp = sim_model.n_mp
        tstep = sim_model.tstep
        w_max = sim_model.w_max

        # pst_list: thickness
        #
        pst_list = [1]
#        pst_list = [ 0.020, 0.019, 0.018, 0.017, 0.021, 0.022, 0.023 ]

#        for pst_param in pst_list:
        for pst_param in pst_list:

            sim_model.thickness = pst_param
            thickness = sim_model.thickness
            print 'thickness', thickness

            # param_key 
            #
            param_key = sim_model_name + '_' + ccs_unit_cell_key + '_' + calibration_test + '_%s_L%g_h%g_sx%gm%gy%gz%g_s%se%s_E%g_nu%g_tol%g_w%g_ts%g_nmp%g' \
                        % (phi_fn_class, length, thickness, shape_x, mid_shape_x, shape_y, shape_z, supprt_flag[0], elstmr_flag[0], E, nu, tolerance, w_max, tstep, n_mp)
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

#            # f-w-diagram_center
#            #
#            file_name = 'f_w_diagram_c_' + param_key + '.pickle'
#            pickle_file_path = join(pickle_path, file_name)
#            file = open(pickle_file_path, 'r')
#            trace = load(file)
#            p.plot(trace.xdata, trace.ydata, color='red')

        #--------------------        
        # experiments
        #--------------------        

        if test_series == 'BT-12c':
            ex_path_V1 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE',
                                   'BT-3PT-12c-6cm-0-Tu-V1.raw')
            ex_path_V2 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE',
                                   'BT-3PT-12c-6cm-0-Tu-V2.raw')
            ex_path_V3 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE',
                                   'BT-3PT-12c-6cm-0-Tu-V3.raw')
            ex_path_V4 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2011-06-10_BT-3PT-12c-6cm-0-TU_ZiE',
                                   'BT-3PT-12c-6cm-0-Tu-V4.raw')
            tests = [ex_path_V1, ex_path_V2, ex_path_V3, ex_path_V4]
            for ex_path in tests:
                ex_run = ExRun(ex_path)
#                 ex_run.ex_type._plot_force_center_deflection_interpolated(p, color='black')
                ex_run.ex_type._plot_force_machine_displacement_wo_elast_analytical_offset(p, color='black')

            # analytical bending stiffness of the beam
            #
            K_linear = 1 / 0.3318  # [MN/m] bending stiffness with respect to a force applied at center of the beam (for E_c(9d) = 26524)
#            K_linear = 1 / 0.3589  # [MN/m] bending stiffness with respect to a force applied at center of the beam (for E_m(9d) = 24524)
            w_linear = 2 * np.array([0., 1.])
            F_linear = 2 * np.array([0., K_linear])
            p.plot(w_linear, F_linear, linestyle='--')
    #        ex_run.ex_type._plot_analytical_bending_stiffness(p)

            # analytical cracking load of the beam
            #
            w_cr = 2 * np.array([0., 1.])
            F_cr = np.array([2.92, 2.92])  # corresponds to l = 1.15m and f_cfl = approx. 7 MPa
            p.plot(w_cr, F_cr, linestyle='--')

            # format plot font as time new roman within the predefined limits  
            #
            format_plot(p, xlim=25, ylim=14, xlabel='displacement [mm]', ylabel='force [kN]')


        if test_series == 'BT-6c':
            ex_path_bs4_V1 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-02_BT-6c-2cm-0-TU_bs4',
                                   'BT-6c-2cm-0-TU-V1_bs4.DAT')
            ex_path_bs4_V2 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-02_BT-6c-2cm-0-TU_bs4',
                                   'BT-6c-2cm-0-TU-V2_bs4.DAT')
            ex_path_bs4_V3 = os.path.join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-02_BT-6c-2cm-0-TU_bs4',
                                   'BT-6c-2cm-0-TU-V3_bs4.DAT')
            ex_path_bs5_V1 = join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-18_BT-6c-2cm-0-TU_bs5',
                                   'BT-6c-2cm-0-TU-V1_bs5.DAT')
            ex_path_bs5_V2 = join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-18_BT-6c-2cm-0-TU_bs5',
                                   'BT-6c-2cm-0-TU-V2_bs5.DAT')
            ex_path_bs5_V3 = join(simdb.exdata_dir, 'bending_tests', 'three_point', '2013-07-18_BT-6c-2cm-0-TU_bs5',
                                   'BT-6c-2cm-0-TU-V3_bs5.DAT')
            tests = [
                     ex_path_bs4_V2,
                     ex_path_bs4_V3,
#                     ex_path_bs5_V1,
#                     ex_path_bs5_V3,
                    ]
            for ex_path in tests:
                ex_run = ExRun(ex_path)
                ex_run.ex_type._plot_force_gauge_displacement_with_analytical_offset(p, linewidth=1.5)
#                ex_run.ex_type._plot_force_gauge_displacement(p, linewidth=1.5)
#                ex_run.ex_type._plot_force_gauge_displacement_interpolated(p, linewidth=2.)
#                ex_run.ex_type._plot_force_machine_displacement(p, linewidth=2.5)
#                ex_run.ex_type._plot_force_machine_displacement_wo_elast(p)
#                ex_run.ex_type._plot_force_machine_displacement_wo_elast_interpolated(p, linewidth=2.5)

            # analytical bending stiffness of the beam
            #
            K_linear = 1 / 1.369  # [MN/m] bending stiffness with respect to a force applied at center of the beam (E_c(28d) = 22213)
#             K_linear = 1 / 1.536  # [MN/m] bending stiffness with respect to a force applied at center of the beam (E_m(28d) = 19800)
            w_linear = 2 * np.array([0., 1.])
            F_linear = 2 * np.array([0., K_linear])
            p.plot(w_linear, F_linear, linestyle='--')
    #        ex_run.ex_type._plot_analytical_bending_stiffness(p)

            # analytical cracking load of the beam
            #
            w_cr = 2 * np.array([0., 1.])
            F_cr = np.array([0.48, 0.48])  # corresponds to l = 0.46m and f_cfl = approx. 8.4 MPa
            p.plot(w_cr, F_cr, linestyle='--')

            # plot sim curve as time new roman within the predefined limits
            #
            format_plot(p, xlim=20, ylim=3.5, xlabel='displacement [mm]', ylabel='force [kN]')

        p.show()

