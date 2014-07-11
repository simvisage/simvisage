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

    from matresdev.db.exdb.ex_run import ExRun
    import pylab as p

    #---------------------------
    # configure filename for outputs (save pdf-file and tex-file)
    #---------------------------
    #
    save_table_to_file = 'true'
    save_fig_to_file = 'true'

    #---------------------------
    # tensile test results ( AR-1200-TU )
    #---------------------------
    # AR-glas, tissue, 1200tex
    #
    test_series_name = 'TT-6g-2cm-0-TU_Serie-1-2-3'

    # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
    # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
    #
    k_rho = 1.098  # correction factor for ARG-1200-TU (a_tex*b = 54mm2/m*0.10m*6layers =32,4mm2; A_tex = 0.447mm2/Roving * 11 Rovings/layer * 6 layers=29,5mm2)

    # specify limits for the plot
    xlim = 14.
    ylim = 20.

    path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-V1.DAT')
    path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-V2.DAT')
    path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-V3.DAT')
    #
    path_V4 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-03-11_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-TU-V1.DAT')
    path_V5 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-03-11_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-TU-V2.DAT')
    path_V6 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-03-11_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-TU-V3.DAT')
    #
    path_V7 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V1_bs4.DAT')
    path_V8 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V2_bs4.DAT')
    path_V9 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V3_bs4.DAT')

    path_list = [
                 path_V1, path_V2, path_V3,
                 path_V4, path_V5, path_V6,
                 path_V7, path_V8, path_V9,
                 ]
    label_list = [
                  'Serie 1', None, None,
                  'Serie 2', None, None,
                  'Serie 3', None, None,
                 ]
    color_list = [
                  'g', 'g', 'g',
                  'b', 'b', 'b',
                  'r', 'r', 'r',
                  ]

    #---------------------------
    # tensile test results ( CAR-800-TU )
    #---------------------------
    # CAR-800-TU-1v1l, tissue, 800tex
    #
    test_series_name = 'TT-12c-6cm-0-TU_SH2'

    k_rho = 1.057  # correction factor for CAR-800-TU (a_tex*b = 53.9mm2/m*0.14m*12 =90,55mm2; A_tex = 0.446mm2/Roving * 16 Rovings/layer * 12 layers=85,63mm2)

    # specify limits for the plot
    #
#    xlim = 8.
#    ylim = 20.

    path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V1.DAT')
    path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V2.DAT')
    path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V3.DAT')

    path_list = [
                 path_V1, path_V2, path_V3,
                 ]
    label_list = [
                 'Serie 1', None, None,
                 ]
    color_list = [
                 'b', 'b', 'b',
                 ]
    #---------------------------
    # slab test results ( CAR-800-TU )
    #---------------------------
    # CAR-800-TU-1v1l, tissue, 800tex
    #
    test_series_name = 'ST-12c-6cm-0-TU_SH2'
    path_V1 = os.path.join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU', 'ST-12c-6cm-u-TU.DAT')
    path_list = [path_V1]
    label_list = [None]
    color_list = ['k']

#    #---------------------------
#    # bending test results ( CAR-800-TU )
#    #---------------------------
#    # CAR-800-TU-1v1l, tissue, 800tex
#    #
#    test_series_name = 'BT-4PT-12c-6cm-0-TU'
#    path_V1 = os.path.join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU', 'BT-4PT-12c-6cm-SH4', 'BT-4PT-12c-6cm-SH4-V1.DAT')
#    path_V2 = os.path.join(simdb.exdata_dir, 'bending_tests', 'four_point', '2012-04-03_BT-4PT-12c-6cm-0-TU', 'BT-4PT-12c-6cm-SH4', 'BT-4PT-12c-6cm-SH4-V2.DAT')
#    path_list = [path_V1, path_V2]
#    label_list = [None, None]
#    color_list = ['k', 'k']

    # --------------------------------
    # plot sig-eps-curves
    # --------------------------------
    fig = p.figure(facecolor='white')
#    fig.set_size_inches(8, 6)
    fig.set_size_inches(16, 12)
    fig.set_size_inches(12, 9)

    for i, ex_path in enumerate(path_list):
        ex_run = ExRun(ex_path)
#        ex_run.ex_type._plot_tex_stress_strain_asc(p, k_rho=k_rho, color=color_list, linewidth=2., linestyle='-', label=label_list[i], xscale=1000.)
#        ex_run.ex_type._plot_comp_stress_strain_asc(p, k_rho=k_rho, color=color_list[i], linewidth=2., linestyle='-', label=label_list[i], xscale=1000.)
        # ex_run.ex_type._plot_force_center_deflection_interpolated_avg(p)
        ex_run.ex_type._plot_force_deflection_avg_interpolated(p, linewidth=2.)
#        ex_run.ex_type._plot_force_deflection_center(p, linewidth=2.)

    format_plot(p, xlabel='Durchbiegung $w$ [mm]', ylabel='Kraft [kN]', xlim=35., ylim=60.)
#    format_plot(p, xlabel='Dehnung [1E-3]', ylabel='Kompositspannung [MPa]', xlim=xlim, ylim=ylim)

    axes = p.gca()
    axes.xaxis.grid(True, which='major')
#    axes.xaxis.grid(True, which='minor')
    axes.yaxis.grid(True, which='major')
#    axes.yaxis.grid(True, which='minor')

    p.legend(prop=font, loc=7)

    # --------------------------------
    # save figure
    # --------------------------------
    if save_fig_to_file:
        img_dir = os.path.join(simdb.exdata_dir, 'img_dir')
        # check if directory exist otherwise create
        #
        if os.path.isdir(img_dir) == False:
            os.makedirs(img_dir)
        test_series_dir = os.path.join(img_dir, test_series_name)
        # check if directory exist otherwise create
        #
        if os.path.isdir(test_series_dir) == False:
            os.makedirs(test_series_dir)
#        filename = os.path.join(test_series_dir, 'sigc-epsu.pdf')
        filename = os.path.join(test_series_dir, 'F-w.pdf')
        p.savefig(filename, format='pdf')
        filename = os.path.join(test_series_dir, 'F-w.png')
        p.savefig(filename, format='png',)
        print 'figure saved to file %s' % (filename)

    p.show()

    # --------------------------------
    # get tabular values
    # --------------------------------
    eps_u_list = []
    sig_c_list = []
    sig_tex_list = []
    for i, ex_path in enumerate(path_list):
        ex_run = ExRun(ex_path)
        # get maximum values for stress and strain for table summary
        #
        eps_u_list += [np.max(ex_run.ex_type.eps_c_interpolated) * 1000.]
        sig_c_list += [np.max(ex_run.ex_type.sig_c_interpolated)]
        # correct values for 'sig_tex' by real reinforcement ration determined by numer of rovings instead ob specimn width
        #
        sig_tex_list += [np.max(ex_run.ex_type.sig_tex_interpolated) * k_rho]

    eps_u_arr = np.hstack(eps_u_list)
    sig_c_arr = np.hstack(sig_c_list)
    sig_tex_arr = np.hstack(sig_tex_list)
    print 'eps_u_arr', eps_u_arr
    print 'sig_c_arr', sig_c_arr
    print 'sig_tex_arr', sig_tex_arr

    # --------------------------------
    # calculate average, standard deviation and coreficient of variation
    # --------------------------------
    sig_c_avg = np.average(sig_c_arr)
    sig_tex_avg = np.average(sig_tex_arr)
    eps_u_avg = np.average(eps_u_arr)
    sig_c_dev = np.std(sig_c_arr)
    sig_tex_dev = np.std(sig_tex_arr)
    eps_u_dev = np.std(eps_u_arr)
    sig_c_Vx = np.var(sig_c_arr)
    sig_tex_Vx = np.var(sig_tex_arr)
    eps_u_Vx = np.var(eps_u_arr)

    # --------------------------------
    # save tabular values as Latex formated table
    # --------------------------------
    if save_table_to_file:
        img_dir = os.path.join(simdb.exdata_dir, 'img_dir')
        # check if directory exist otherwise create
        #
        if os.path.isdir(img_dir) == False:
            os.makedirs(img_dir)
        test_series_dir = os.path.join(img_dir, test_series_name)
        # check if directory exist otherwise create
        #
        if os.path.isdir(test_series_dir) == False:
            os.makedirs(test_series_dir)
        filename = os.path.join(test_series_dir, 'sigc-sigtex-epsu.tex')
        f = open(filename, 'w')

        output_str = ""
        output_str = output_str + r"\begin{table} " + "\n"
        output_str = output_str + r"\begin{tabular}{c c c c} " + "\n"
        output_str = output_str + r"\toprule " + "\n"
        output_str = output_str + r"\bf{Vers.-Nr.} & \bf{$\sigma_\textrm{c}$} & \bf{$\sigma_\textrm{tex}$} & {$\varepsilon_\textrm{u}$} \\ " + "\n"
        output_str = output_str + r"\midrule " + "\n"
        #
        for i in range(len(sig_c_arr)):
            output_str = output_str + "{V" + str(i + 1) + "} & {%.1f} & {%.0f} & {%.1f}" % (sig_c_arr[i], sig_tex_arr[i], eps_u_arr[i])
            output_str = output_str + r" \\ " + "\n"
        output_str = output_str + r"\bottomrule " + "\n"
        output_str = output_str + "{$m$} & {%.1f} & {%.0f} & {%.1f} " % (sig_c_avg, sig_tex_avg, eps_u_avg)
        output_str = output_str + r" \\ " + "\n"
        output_str = output_str + "{$s$} & {%.1f} & {%.0f} & {%.1f} " % (sig_c_dev, sig_tex_dev, eps_u_dev)
        output_str = output_str + r" \\ " + "\n"
        output_str = output_str + "{$v$} & {%.1f} & {%.0f} & {%.1f} " % (sig_c_Vx, sig_tex_Vx, eps_u_Vx)
        output_str = output_str + r" \\ " + "\n"
        output_str = output_str + r"\bottomrule " + "\n"
        output_str = output_str + r"\end{tabular} " + "\n"
        output_str = output_str + r"\end{table} " + "\n"

        print 'output_str \n', output_str

        f.write(output_str)
        print 'table data saved to file %s' % (filename)

