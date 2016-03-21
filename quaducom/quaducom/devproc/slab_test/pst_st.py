'''
Created on Sep 9, 2013

@author: alexander
'''
from traits.api import \
    Array, Bool, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, \
    Callable, List, TraitDict, Any, Range, \
    Delegate, Event, on_trait_change, Button, \
    Interface, implements, Property, cached_property

from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
    IPhiFn, PhiFnGeneralExtended
from matresdev.db.simdb import \
    SimDB
import numpy as np
from quaducom.devproc.tensile_test.dog_bone.test_reports import format_plot
from sim_st import SimSTDB
import pickle

simdb = SimDB()


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
#    test_series = 'ST-12c'
    test_series = 'ST-6c'

    #-----------------------------------------------
    # ST-10g: AG-glas slab tests (125 cm / 3 cm) with tricot-binding:
    #-----------------------------------------------
    if test_series == 'ST-10g':

        sim_model = SimSTDB(

            # calibration for: age = XXd; E_m = XX MPa; nu = XX
            #
            ccs_unit_cell_key='FIL-10-09_2D-02-06a_0.00273_90_0',
            #                            calibration_test='TT11-10a-average',
            #                            calibration_test='TT-10g-3cm-a-TR-average_age28_Ec18709.5_nu0.2_nsteps200_smoothed',
            calibration_test='TT-10g-3cm-a-TR-average_age28_Ec28700_nu0.2_nsteps200_smoothed',
            #                            age = 23,
            #
            thickness=0.03,
            length=1.25,
            #
            radius_plate=0.095,  # D=8cm
            #
            elstmr_flag=False,
            supprt_flag=True,
            geo_st_flag=True,
            #                            supprt_flag=False,
            #                            geo_st_flag=False,
            #
            # coarse mesh:
            shape_xy=10,
            shape_R=2,
            shape_z=2,
            shape_supprt_xy=2,
            #
            tstep=0.05,
            tmax=1.0,
            #
            w_max=-0.055,
            # 'NOTE: tloop.norm switched to "max(abs(x))"'
            tolerance=0.0001,  # #[MN]0.0001#1e-6#1e-8#0.0005
            ord=np.inf,  # "norm = max(abs(x_i))"
        )

    #-----------------------------------------
    # ST-12c-6cm; L = 1,25m; t = 6 cm
    #-----------------------------------------
    if test_series == 'ST-12c':

        sim_model = SimSTDB(

            thickness=0.06,
            #                             thickness=0.0577,
            #                            thickness=0.0554,

            length=1.25,
            radius_plate=0.095,  # D=8cm

            ccs_unit_cell_key='FIL-10-09_2D-05-11_0.00462_all0',

            #                            calibration_test='TT-12c-6cm-0-TU-V1_ZiE-S1_age23_Em27975.8_nu0.2_nsteps100',
            #                            calibration_test='TT-12c-6cm-0-TU-SH2F-V3_a23d_nu02_s100',
            #                            calibration_test = 'TT-12c-6cm-0-TU-SH2F-V3_a23d_nu02_s50',
            #                            calibration_test = 'TT-12c-6cm-TU-SH1F-V1',

            # fresh-in-fresh fabrication
            #
            #                            calibration_test='TT-12c-6cm-0-TU-SH2F-V3_age23_Em27975.8_nu0.2_nsteps100',
            #                            calibration_test='TT-12c-6cm-0-TU-SH2F-V3_age23_Ec29494_nu0.2_nsteps100',

            # initial stiffness corresponds to tensile test; cut of at eps=0.007; based on smoothed experimental curve without oscillation
            #
            calibration_test='TT-12c-6cm-0-TU-SH2-V1_age26_Ec29100_nu0.2_nsteps100_maxeps0.007_smoothed',

            # number of microplanes
            #
            n_mp=30,

            # age of the slab at the time of testing
            age=26,
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
            shape_z=2,
            shape_supprt_xy=2,
            #
            # fine mesh:
            #                            shape_z=2,
            #                            shape_xy=26,
            #                            shape_R=4,
            #                            shape_supprt_xy=4,
            #
            w_max=-0.035,
            tstep=0.01,
            #                            tstep=1.00,
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

            # Em = 19800. MPa
            #
            #                            calibration_test='TT-6c-2cm-0-TU-V3_bs1_age28_Em19800_nu0.2_nsteps100',
            #                             calibration_test='TT-6c-2cm-0-TU-V1_bs2_age28_Em19800_nu0.2_nsteps100',
            #                             calibration_test='TT-6c-2cm-0-TU-V1_bs3_age28_Em19800_nu0.2_nsteps100',
            #                             calibration_test='TTb-6c-2cm-0-TU-V3_bs5_age28_Em19800_nu0.2_nsteps100',

            # Ec = 22213.2 MPa
            #
            #                             calibration_test='TT-6c-2cm-0-TU-V3_bs1_age28_Ec22213.2_nu0.2_nsteps100',
            #                             calibration_test='TT-6c-2cm-0-TU-V1_bs2_age28_Ec22213.2_nu0.2_nsteps100',
            #                             calibration_test='TT-6c-2cm-0-TU-V1_bs3_age28_Ec22213.2_nu0.2_nsteps100',
            #                             calibration_test='TTb-6c-2cm-0-TU-V3_bs5_age28_Ec22213.2_nu0.2_nsteps100',

            #                             calibration_test='TTb-6c-2cm-0-TU-V2_bs4_age28_Ec18709.5_nu0.2_nsteps200_smoothed',
            calibration_test='TTb-6c-2cm-0-TU-V1_bs5_age28_Ec18709.5_nu0.2_nsteps100_smoothed',
            #
            age=28,
            #
            thickness=0.02,
            #                            thickness_plain_concrete=0.00286,
            length=0.80,
            radius_plate=0.04,  # D=8cm
            #
            plain_concrete_flag=False,
            elstmr_flag=False,
            supprt_flag=True,
            geo_st_flag=True,
            #
#             shape_xy=6,  # coarse mesh
            shape_xy=11,  # medium mesh
#             shape_xy=20,  # fine mesh
            shape_R=2,
            shape_z=2,
            #
            w_max=-0.050,
            tstep=0.005,
            tmax=0.66,

            # quadratic elements (ngp=2x2x2) and w_max = 0.050m --> w_u=0.0289m
#             w_max=-0.050,
#             tmax=0.575,

            # cubic elements (ngp=4x4x4) and w_max = 0.070m --> w_u=0.0322m
#            w_max=-0.050,
#            tmax=0.6450, # no failure reached up to that load level!
            #
#            w_max=-0.070,
#             tmax=0.4600,
            #
            ord=np.inf,  # "norm = max(abs(x_i))"
            tolerance=0.00005,  # [MN]0.0001#1e-6#1e-8#0.0005
#             tolerance=0.00005,  # [MN]0.0001#1e-6#1e-8#0.0005
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
# print 'self.ccs_unit_cell_ref.damage_function_list',
# [sim_model.ccs_unit_cell_ref.damage_function_list[i].calibration_test
# for i in range(len(sim_model.ccs_unit_cell_ref.damage_function_list))]

        phi_fn.mfn.plot(p, color='black', linewidth=3)
#        phi_fn_ext.mfn.plot(p, color = 'blue', linewidth = 2)
#        phi_fn_exp.mfn.plot(p, color = 'green', linewidth = 3)

        xmax = sim_model.damage_function.xdata[-1]
        print 'xmax', xmax
        x = np.linspace(0, 3 * xmax, 1000)
        phi_fn = np.frompyfunc(phi_fn, 1, 1)
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
        import traits
        traits.has_traits.CHECK_INTERFACES = 0
        import os

        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource=sim_model)

        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = os.path.join(simdb.simdata_dir, 'pickle_files')
        png_path = os.path.join(simdb.simdata_dir, 'png_files')

        if not os.path.exists(pickle_path):
            os.mkdir(pickle_path)
            os.mkdir(png_path)

        # pstudy: None
        #
#        pst_list = [1]

        # pstudy: thickness
        #
#        pst_list = [ 0.0554 ]  # 0.06, ]
#        pst_list = [ 0.01714 ]  # 0.02 ]
        pst_list = [0.02]

        # pstudy: n_mp
        #
#        pst_list = [ 30 ]

        # pstudy: calibration test
        #
#        pst_list = [ 'TT-12c-6cm-0-TU-SH2F-V3_a23d_nu02_s100' , 'TT-12c-6cm-TU-SH1F-V1' ]

#        pst_list = [
#                       # Em = 19800. MPa
#                       #
#                       'TT-6c-2cm-0-TU-V3_bs1_age28_Em19800_nu0.2_nsteps100',
#                       'TT-6c-2cm-0-TU-V1_bs2_age28_Em19800_nu0.2_nsteps100',
#                       'TT-6c-2cm-0-TU-V1_bs3_age28_Em19800_nu0.2_nsteps100',
#                       'TTb-6c-2cm-0-TU-V3_bs5_age28_Em19800_nu0.2_nsteps100',
#
#                       # Ec = 22213.2 MPa
#                       #
#                       'TT-6c-2cm-0-TU-V3_bs1_age28_Ec22213.2_nu0.2_nsteps100',
#                       'TT-6c-2cm-0-TU-V1_bs2_age28_Ec22213.2_nu0.2_nsteps100',
#                       'TT-6c-2cm-0-TU-V1_bs3_age28_Ec22213.2_nu0.2_nsteps100',
#                       'TTb-6c-2cm-0-TU-V3_bs5_age28_Ec22213.2_nu0.2_nsteps100'

        # initial stiffness corresponds to tensile test; cut of at eps=0.007; based on smoothed experimental curve without oscillation
        #
#                       'TT-12c-6cm-0-TU-SH2-V1_age26_Ec29100_nu0.2_nsteps100_maxeps0.007_smoothed',
#                   ]

        # pstudy: phi_fn
        #
#        pst_list = [ PhiFnGeneral, PhiFnGeneralExtended, PhiFnGeneralExtendedExp ]

        for pst_param in pst_list:

            sim_model.thickness = pst_param
#            sim_model.n_mp = pst_param
#            sim_model.calibration_test = pst_param
#            sim_model.phi_fn_class = pst_param

            # white background for diagram
            p.figure(facecolor='white', figsize=(12, 9))

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
            E = sim_model.specmn_mats.E
            print 'E', E
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
            param_key = sim_model_name + '_' + ccs_unit_cell_key + '_' + calibration_test + '_%s_L%gh%gR%g_sxy%gz%gR%g_s%sg%s_E%g_nu%g_tol%g_w%g_ts%g_nmp%g' \
                % (phi_fn_class, length, thickness, radius_plate, shape_xy, shape_z, shape_R, supprt_flag[0], geo_st_flag[0], E, nu, tolerance, w_max, tstep, n_mp)
            print 'param_key = %s' % param_key

            # f-w-diagram_center
            #
            sim_model.f_w_diagram_center.refresh()
            file_name = 'f_w_diagram_c_' + param_key + '.pickle'
            pickle_file_path = os.path.join(pickle_path, file_name)
            file = open(pickle_file_path, 'w')
            pickle.dump(sim_model.f_w_diagram_center.trace, file)
            print 'pickle file saved to file: %s' % file_name
            file.close()
            sim_model.f_w_diagram_center.trace.mpl_plot(p, color='grey')

            # f-w-diagram_supprt
            #
            sim_model.f_w_diagram_supprt.refresh()
            file_name = 'f_w_diagram_supprt_' + param_key + '.pickle'
            pickle_file_path = os.path.join(pickle_path, file_name)
            file = open(pickle_file_path, 'w')
            pickle.dump(sim_model.f_w_diagram_supprt.trace, file)
            print 'pickle file saved to file: %s' % file_name
            file.close()
#            sim_model.f_w_diagram_supprt.trace.mpl_plot(p, color='grey')

            w_asc = sim_model.f_w_diagram_supprt.trace.xdata
            f_asc = sim_model.f_w_diagram_supprt.trace.ydata
            print 'f_asc.shape', f_asc.shape
            print 'w_asc.shape', w_asc.shape
            fw_arr = np.hstack([f_asc[:, None], w_asc[:, None]])
            print 'fw_arr.shape', fw_arr.shape
            np.savetxt(param_key + '.csv', fw_arr, delimiter=';')

            #--------------------
            # experiments
            #--------------------

            if test_series == 'ST-12c':
                # ST-12c-6cm-TU
                #
                ex_path = os.path.join(simdb.exdata_dir, 'slab_tests',
                               '2011-12-15_ST-12c-6cm-u-TU', 'ST-12c-6cm-u-TU.DAT')
                ex_run = ExRun(ex_path)
                ex_run.ex_type._plot_force_center_deflection_interpolated(p)
                format_plot(
                    p, xlim=35, ylim=100, xlabel='$w_\mathrm{Mitte}$ [mm]', ylabel='$F$ [kN]')

            if test_series == 'ST-10g':
                # ST-10a
                #
                ex_path_TRC10 = os.path.join(
                    simdb.exdata_dir, 'slab_tests', '2010-03-08_ST-10g-3cm-a-FR_TRC10', 'ST-10g-3cm-a-FR-TRC10.DAT')
                ex_path_TRC11 = os.path.join(
                    simdb.exdata_dir, 'slab_tests', '2010-03-09_ST-10g-3cm-a-FR_TRC11', 'ST-10g-3cm-a-FR-TRC11.DAT')
                ex_path_TRC12 = os.path.join(
                    simdb.exdata_dir, 'slab_tests', '2010-03-10_ST-10g-3cm-a-FR_TRC12', 'ST-10g-3cm-a-FR-TRC12.DAT')
                tests = [ex_path_TRC10, ex_path_TRC11, ex_path_TRC12]
                for ex_path in tests:
                    ex_run = ExRun(ex_path)
                    ex_run.ex_type._plot_force_center_deflection_interpolated(
                        p)

                format_plot(
                    p, xlim=35, ylim=70, xlabel='displacement [mm]', ylabel='force [kN]')

            if test_series == 'ST-6c':
                # ST-6c-2cm-TU_bs2
                #
                ex_path = os.path.join(simdb.exdata_dir, 'slab_tests',
                               '2013-07-10_ST-6c-2cm-TU_bs2', 'ST-6c-2cm-TU_bs2.DAT')
                ex_run = ExRun(ex_path)
                ex_run.ex_type._plot_force_center_deflection(
                    p, offset_w=0.0385)

                format_plot(
                    p, xlim=80, ylim=20, xlabel='displacement [mm]', ylabel='force [kN]')

            #------------------------------------------------------------------
            # plot sim curve as time new roman within the predefined limits
            #------------------------------------------------------------------
            #
            png_file_path = os.path.join(png_path, param_key + '.png')
            pdf_file_path = os.path.join(png_path, param_key + '.pdf')
#            p.title(param_key, fontsize=8)
            p.savefig(png_file_path, dpi=600.)
            p.savefig(pdf_file_path)
            print 'png-file saved to file: %s' % png_file_path
            print 'pdf-file saved to file: %s' % pdf_file_path
#            p.show()

        app.main()

    #------------------------------
    # show last results
    #------------------------------

    if do == 'show_last_results':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        # white background for diagram
        p.figure(facecolor='white', figsize=(12, 9))

        pickle_path = os.path.join(simdb.simdata_dir, 'pickle_files')
        png_path = os.path.join(simdb.simdata_dir, 'png_files')

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
        E = sim_model.specmn_mats.E
        E_m = sim_model.E_m
        nu = sim_model.nu
        tolerance = sim_model.tolerance
        phi_fn_class = sim_model.phi_fn_class.__name__
        supprt_flag = str(sim_model.supprt_flag)
        geo_st_flag = str(sim_model.geo_st_flag)
        n_mp = sim_model.n_mp
        tstep = sim_model.tstep
        w_max = sim_model.w_max

        pst_list = [  # Em = 19800. MPa
            'TT-6c-2cm-0-TU-V3_bs1_age28_Em19800_nu0.2_nsteps100',
            #                    'TT-6c-2cm-0-TU-V1_bs2_age28_Em19800_nu0.2_nsteps100',
            #                    'TT-6c-2cm-0-TU-V1_bs3_age28_Em19800_nu0.2_nsteps100',
            #                    'TTb-6c-2cm-0-TU-V3_bs5_age28_Em19800_nu0.2_nsteps100',

            # Ec = 22213.2 MPa
            #                    'TT-6c-2cm-0-TU-V3_bs1_age28_Ec22213.2_nu0.2_nsteps100',
            #                    'TT-6c-2cm-0-TU-V1_bs2_age28_Ec22213.2_nu0.2_nsteps100',
            #                     'TT-6c-2cm-0-TU-V1_bs3_age28_Ec22213.2_nu0.2_nsteps100',
            #                     'TTb-6c-2cm-0-TU-V3_bs5_age28_Ec22213.2_nu0.2_nsteps100'
        ]

        for pst_param in pst_list:

            calibration_test = pst_param

            # param_key
            #
            param_key = sim_model_name + '_' + ccs_unit_cell_key + '_' + calibration_test + '_%s_L%gh%gR%g_sxy%gz%gR%g_s%sg%s_E%g_nu%g_tol%g_w%g_ts%g_nmp%g' \
                % (phi_fn_class, length, thickness, radius_plate, shape_xy, shape_z, shape_R, supprt_flag[0], geo_st_flag[0], E, nu, tolerance, w_max, tstep, n_mp)
            param_key = 'SimSTDB_barrelshell_2D-05-11_0.00286_all0_TTb-6c-2cm-0-TU-V1_bs5_age28_Ec18709.5_nu0.2_nsteps100_smoothed_PhiFnGeneralExtended_L0.8h0.02R0.04_sxy10z2R2_sTgT_E18709.5_nu0.2_tol5e-05_w-0.04_ts0.05_nmp30'

            print 'param_key = %s' % param_key

            #------------------
            # simulation
            #------------------

            # f-w-diagram_supprt
            #
            file_name = 'f_w_diagram_supprt_' + param_key + '.pickle'

#            # h_eff = 0.0554 m
#            file_name = 'f_w_diagram_supprt_SimSTDB_FIL-10-09_2D-05-11_0.00462_all0_TT-12c-6cm-0-TU-SH2-V1_age26_Ec29100_nu0.2_nsteps100_maxeps0.007_smoothed_PhiFnGeneralExtended_L1.25h0.0554R0.095_sxy10z2R2_sTgT_E29100_nu0.2_tol0.0001_w-0.035_ts0.01_nmp30.pickle'
#            pickle_file_path = join(pickle_path, file_name)
#            file = open(pickle_file_path, 'r')
#            trace = load(file)
#            n_max = -28
#            p.plot(trace.xdata[:n_max], trace.ydata[:n_max], color='grey', linestyle='--', linewidth=2.)

            # h = 0.06 m
#            file_name = 'f_w_diagram_supprt_SimSTDB_FIL-10-09_2D-05-11_0.00462_all0_TT-12c-6cm-0-TU-SH2-V1_age26_Ec29100_nu0.2_nsteps100_maxeps0.007_smoothed_PhiFnGeneralExtended_L1.25h0.06R0.095_sxy10z2R2_sTgT_E29100_nu0.2_tol0.0001_w-0.035_ts0.01_nmp30.pickle'
            pickle_file_path = os.path.join(pickle_path, file_name)
            file = open(pickle_file_path, 'r')
            trace = pickle.load(file)
            p.plot(trace.xdata, trace.ydata, color='grey',
                   linestyle='-', linewidth=2.)

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
            ex_path = os.path.join(simdb.exdata_dir, 'slab_tests',
                           '2011-12-15_ST-12c-6cm-u-TU', 'ST-12c-6cm-u-TU.DAT')
            ex_run = ExRun(ex_path)
            ex_run.ex_type._plot_force_center_deflection_interpolated(
                p, linewidth=2., plot_elastic_stiffness=False)
            # PT-12c-6cm-TU
            format_plot(
                p, xlim=35, ylim=70, xlabel='$w_\mathrm{Mitte}$ [mm]', ylabel='$F$ [kN]')

        if test_series == 'ST-10g':
            # PT-10a
            #
            ex_path_TRC10 = os.path.join(
                simdb.exdata_dir, 'slab_tests', '2010-03-08_ST-10g-3cm-a-FR_TRC10', 'ST-10g-3cm-a-FR-TRC10.DAT')
            ex_path_TRC11 = os.path.join(
                simdb.exdata_dir, 'slab_tests', '2010-03-09_ST-10g-3cm-a-FR_TRC11', 'ST-10g-3cm-a-FR-TRC11.DAT')
            ex_path_TRC12 = os.path.join(
                simdb.exdata_dir, 'slab_tests', '2010-03-10_ST-10g-3cm-a-FR_TRC12', 'ST-10g-3cm-a-FR-TRC12.DAT')
            tests = [ex_path_TRC10, ex_path_TRC11, ex_path_TRC12]
            for ex_path in tests:
                ex_run = ExRun(ex_path)
                ex_run.ex_type._plot_force_center_deflection(p)
            # plot sim curve as time new roman within the predefined limits
            #
            format_plot(
                p, xlim=34, ylim=54, xlabel='Durchbiegung [mm]', ylabel='Kraft [kN]')

        if test_series == 'ST-6cX':
            # ST-6c-2cm-TU_bs2
            #
            ex_path = os.path.join(simdb.exdata_dir, 'slab_tests',
                           '2013-07-10_ST-6c-2cm-TU_bs2', 'ST-6c-2cm-TU_bs2.DAT')
            ex_run = ExRun(ex_path)
            # specify offset as displacement measurement starts a little to late at a small force > 0
            # shift experimental curve to simulation curve value for the displacement at this load level
            # corresponding to linear elastic approximation of slab stiffness;
            #
#            ex_run.ex_type._plot_force_center_deflection(p, offset_w=0.0385, color='black', linewidth=1.2)
            ex_run.ex_type._plot_force_center_deflection_interpolated(p)

            ex_path = os.path.join(simdb.exdata_dir, 'slab_tests',
                           '2013-07-11_ST-6c-2cm-TU_bs3-4-5-Aramis3d', 'ST-6c-2cm-TU_bs5.DAT')
            ex_run = ExRun(ex_path)
            # specify offset as displacement measurement starts a little to late at a small force > 0
            # shift experimental curve to simulation curve value for the displacement at this load level
            # corresponding to linear elastic approximation of slab stiffness;
            #
#            ex_run.ex_type._plot_force_center_deflection(p, offset_w=0.0385, color='black', linewidth=1.2)
            ex_run.ex_type._plot_force_center_deflection_interpolated(p)

            # plot sim curve as time new roman within the predefined limits
            #
            format_plot(
                p, xlim=80, ylim=20, xlabel='Durchbiegung [mm]', ylabel='Kraft [kN]')

#        p.plot(np.array([0, 1]), np.array([0, 1]), color='grey', linestyle='--', linewidth=1.5)
        p.show()
