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
#    p.rcParams['text.usetex'] = True
    p.rcParams['text.latex.unicode'] = True

    #---------------------------
    # configure filename for outputs (save pdf-file and tex-file)
    #---------------------------
    #
    test_series_name = 'TT-12c-6cm-0-TU_Serie-SH3-4cm-SH2-6cm'
    save_table_to_file = 'true'
    save_fig_to_file = 'true'

    # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
    # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
    #
    # correction factor for CAR-800-TU (a_tex*b = 53.9mm2/m*0.14m*12layers =90.5mm2; A_tex = 0.446mm2/Roving * 16 Rovings/layer * 12 layers=85,6mm2)
    k_rho = 1.057

    # plot composite or textile stress-strain curve
    #
    stress_flag = 'comp'
#    stress_flag = 'tex'

    #---------------------------
    # tensile test results ( CAR-800-TU )
    #---------------------------
    # 4 cm
    path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-03-20_TT-12c-4cm-0-TU_SH3', 'TT-12c-4cm-TU-0-SH3-V1.DAT')
    path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-03-20_TT-12c-4cm-0-TU_SH3', 'TT-12c-4cm-TU-0-SH3-V2.DAT')
    path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-03-20_TT-12c-4cm-0-TU_SH3', 'TT-12c-4cm-TU-0-SH3-V3.DAT')
    #
    # 6 cm
    path_V4 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V1.DAT')
    path_V5 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V2.DAT')
    path_V6 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V3.DAT')

    path_list = [
                 path_V1, path_V2, path_V3,
                 path_V4, path_V5, path_V6,
                 ]
    label_list = [
                  r"$h_{2}\,=\,6\,\rm{cm}\,(QS-2 Jahre;\,kontinuierlich)$", None, None,
                  r'$h_{1}\,=\,6\,\rm{cm}\,(QS-SH4,\,synchron)$', None, None,
                 ]
    color_list = [
                  'grey', 'grey', 'grey',
                  'k', 'k', 'k',
                  ]

    fig = p.figure(facecolor='white')
    fig.set_size_inches(8, 6)
    plot_method_str = '_plot_' + stress_flag + '_stress_strain_asc'

    for i, ex_path in enumerate(path_list):
        ex_run = ExRun(ex_path)
        plot_method = getattr(ex_run.ex_type, plot_method_str)
        plot_method(p, k_rho=k_rho, color=color_list[i], linewidth=1.2, linestyle='-', label=label_list[i], xscale=1000., plot_analytical_stiffness=True, interpolated=True)

    # set limits and labels for axis
    #
    if stress_flag == 'tex':
        format_plot(p, xlabel='strain [1E-3]', ylabel='textile stress [MPa]', xlim=8., ylim=1700.)
    if stress_flag == 'comp':
        format_plot(p, xlabel='strain [1E-3]', ylabel='composite stress [MPa]', xlim=8., ylim=25.)

    # plot grid lines
    #
    axes = p.gca()
    axes.xaxis.grid(True, which='major')
    axes.yaxis.grid(True, which='major')

    # plot legend
    #
    p.legend(prop=font, loc=4, frameon=False)  # (4: lower right) (7: center right)

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
        filename = os.path.join(test_series_dir, 'eps-sig' + stress_flag)
        p.savefig(filename + '.pdf', format='pdf')
        p.savefig(filename + '.png', format='png', dpi=600)
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

