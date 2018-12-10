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

from .format_plot import format_plot
from matplotlib.font_manager import FontProperties

from matresdev.db.exdb.ex_run import ExRun
import pylab as p

if __name__ == '__main__':

    #---------------------------
    # configure filename for outputs (save pdf-file and tex-file)
    #---------------------------
    #
    save_fig_to_file = True
    save_table_to_file = False
    stiffness_flag = False
    stiffness_flag_I = False
    stiffness_flag_II = True
    plot_sigtex_flag = True
    legend_flag = True
    legend_loc = 4  # upper right # set to "legend_loc = 4" for lower right
#    sig_flag = 'comp'
    sig_flag = 'tex'

    #---------------------------
#    do = 'ARG-2400-FR-1v1l'
#    do = 'ARG-2400-TR-1v1l'
#    do = 'ARG-1200-TR-2v1l'
    do = 'ARG-1200-TU-1v1l_TTb-2cm'
#    do = 'ARG-1200-TU-1v1l'
#    do = 'CAR-800-TR-2v1l'
#    do = 'CAR-800-TR-1v1l'
#    do = 'CAR-800-TU-1v1l_TTb-2cm'
#     do = 'CAR-800-TU-1v1l_TRC-pretests'
#    do = 'CAR-800-TU-1v1l_TT-2cm'
#     do = 'CAR-800-TU-1v1l_TT-6cm'
#    do = 'CAR-3300-TR-3v3l'
#    do = 'CAR-3300-TR-3v3l'
    #---------------------------
#    do = 'compare_CAR'
#    do = 'compare_ARG'
#    do = 'compare_ARG-CAR'

    #---------------------------
    # tensile test results (ARG-2400-FR-1v1l) / (TT-7g-3cm-0-FR_TRC07)
    #---------------------------
    if do == 'ARG-2400-FR-1v1l':
        interpolated_flag = False
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.0  # correction factor for ARG-2400-FR (a_tex*b = 107,2mm2/m*0.10m*7layers=75,0mm2; A_tex = 0.89mm2/Roving*12 Rovings/layer*7layers=75,1mm2)

        # specify limits for the plot
        xlim = 5.
        ylim = 8.
        #
        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2009-11-26_TT-7g-3cm-0-FR_TRC07', 'TT-7g-3cm-0-FR-TRC07-V1.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2009-11-26_TT-7g-3cm-0-FR_TRC07', 'TT-7g-3cm-0-FR-TRC07-V2.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2009-11-26_TT-7g-3cm-0-FR_TRC07', 'TT-7g-3cm-0-FR-TRC07-V3.DAT')
        #
#        path_V4 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2009-11-26_TT-3g-1cm-0-FR', 'TT-3g-1cm-0-FR-V1.DAT')
#        path_V5 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2009-11-26_TT-3g-1cm-0-FR', 'TT-3g-1cm-0-FR-V2.DAT')
#        path_V6 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2009-11-26_TT-3g-1cm-0-FR', 'TT-3g-1cm-0-FR-V3.DAT')

        path_list = [path_V1, path_V2, path_V3]

        label_list = ['Serie 1 (V1 bis V3)', None, None]

        color_list = ['k', 'k', 'k']

        linestyle_list = ['-', '-', '-']

    #---------------------------
    # tensile test results (ARG-2400-TR-1v1l) / (TT-2g-1cm-0-TR-2400)
    #---------------------------
    if do == 'ARG-2400-TR-1v1l':
        interpolated_flag = False
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.0  # correction factor for ARG-2400-TR (a_tex*b = 107,2mm2/m*0.10m*7layers=75,0mm2; A_tex = 0.89mm2/Roving*12 Rovings/layer*7layers=75,1mm2)

        # specify limits for the plot
        xlim = 8.
        ylim = 12.
        #
        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-11-02_TT-2g-1cm-0-TR_2400tex', 'TT-2g-1cm-0-TR-2400-V1.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-11-02_TT-2g-1cm-0-TR_2400tex', 'TT-2g-1cm-0-TR-2400-V2.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-11-02_TT-2g-1cm-0-TR_2400tex', 'TT-2g-1cm-0-TR-2400-V3.DAT')
#        path_V4 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-11-02_TT-2g-1cm-0-TR_2400tex', 'TT-2g-1cm-0-TR-2400-V4.DAT')
#        path_V5 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-11-02_TT-2g-1cm-0-TR_2400tex', 'TT-2g-1cm-0-TR-2400-V5.DAT')
#        path_V6 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-11-02_TT-2g-1cm-0-TR_2400tex', 'TT-2g-1cm-0-TR-2400-V6.DAT')
        #
        path_list = [path_V1, path_V2, path_V3]
        label_list = ['Serie 1 (V1 bis V3)', None, None]
        color_list = ['k', 'k', 'k']
        linestyle_list = ['-', '-', '-']

    #---------------------------
    # tensile test results (ARG-1200-TR-2v1l) / (TT-8g-3cm-0-TR)
    #---------------------------
    if do == 'ARG-1200-TR-2v1l':
        interpolated_flag = False
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.0  # correction factor for ARG-1200-TR-2v1l (a_tex*b = 71,7mm2/m*0.10m*8layers=57,4mm2; A_tex = 0.447mm2/Roving*16Rovings/layer*8layers=57,2mm2)

        # specify limits for the plot
        xlim = 10.
        ylim = 14.
        #
        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-03-17_TT-8g-3cm-0-TR', 'TT-8g-3cm-0-TR-V1.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-03-17_TT-8g-3cm-0-TR', 'TT-8g-3cm-0-TR-V2.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-03-17_TT-8g-3cm-0-TR', 'TT-8g-3cm-0-TR-V3.DAT')
        #
        path_list = [path_V1, path_V2, path_V3]
        label_list = ['Serie 1 (V1 bis V3)', None, None]
        color_list = ['none', 'k', 'k']
        linestyle_list = ['-', '-', '-']

    #---------------------------
    # tensile test results (AR-1200-TU-1v1l) / (TTb-6g-2cm-0-TU_bs4)
    #---------------------------
    if do == 'ARG-1200-TU-1v1l_TTb-2cm':
        interpolated_flag = True
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.006  # correction factor for ARG-1200-TU (a_tex*b = 54mm2/m*0.10m*6layers =32,4mm2; A_tex = 0.447mm2/Roving * 12 Rovings/layer * 6 layers=32,2mm2)

        # specify limits for the plot
        if sig_flag == 'comp':
            xlim = 12.
            ylim = 14.
        if sig_flag == 'tex':
            xlim = 12.
            ylim = 1500.

        # TTb_bs4
        path_V7 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V1_bs4.DAT')
        path_V8 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V2_bs4.DAT')
        path_V9 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V3_bs4.DAT')
        #
        path_list = [path_V7, path_V8, path_V9]
        label_list = ['Serie TTb', None, None]
        color_list = ['k', 'k', 'k']
        linestyle_list = ['-', '-', '-']

    #---------------------------
    # tensile test results (AR-1200-TU-1v1l) / (TT-6g-2cm-0-TU_bs)
    #---------------------------
    if do == 'ARG-1200-TU-1v1l':
        interpolated_flag = True
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.006  # correction factor for ARG-1200-TU (a_tex*b = 54mm2/m*0.10m*6layers =32,4mm2; A_tex = 0.447mm2/Roving * 12 Rovings/layer * 6 layers=32,2mm2)

        # specify limits for the plot
        xlim = 14.
        ylim = 18.

        # TT-S1
        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-V1.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-V2.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-V3.DAT')
        #
        # TT-S2
        path_V4 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-03-11_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-TU-V1.DAT')
        path_V5 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-03-11_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-TU-V2.DAT')
        path_V6 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-03-11_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-TU-V3.DAT')
        #
        path_list = [path_V1, path_V2, path_V3,
                     path_V4, path_V5, path_V6]
        label_list = ['Serie 1', None, None,
                      'Serie 2', None, None]
        color_list = ['grey', 'grey', 'grey',
                      'k', 'k', 'k']
        linestyle_list = ['-', '-', '-',
                          '-', '-', '-']

    #---------------------------
    # tensile test results (CAR-800-TU-1v1l)-(TTb-6c-2cm-0-TU_bs)
    #---------------------------
    if do == 'CAR-800-TU-1v1l_TTb-2cm':
        interpolated_flag = True
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.007  # correction factor for CAR-800-TU (a_tex*b = 53.9mm2/m*0.10m*6layers =32,34mm2; A_tex = 0.446mm2/Roving * 12 Rovings/layer * 6 layers=32,11mm2)

        # specify limits for the plot
        if sig_flag == 'comp':
            xlim = 9.
            ylim = 30.
        if sig_flag == 'tex':
            xlim = 9.
            ylim = 1500.
        #
        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d', 'TTb-6c-2cm-0-TU-V1_bs4.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d', 'TTb-6c-2cm-0-TU-V2_bs4.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d', 'TTb-6c-2cm-0-TU-V3_bs4.DAT')
        #
#        path_V4 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-18_TTb-6c-2cm-0-TU_bs5', 'TTb-6c-2cm-0-TU-V1_bs5.DAT')
#        path_V5 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-18_TTb-6c-2cm-0-TU_bs5', 'TTb-6c-2cm-0-TU-V2_bs5.DAT')
#        path_V6 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-18_TTb-6c-2cm-0-TU_bs5', 'TTb-6c-2cm-0-TU-V3_bs5.DAT')
        #
        path_list = [path_V1, path_V2, path_V3]
        label_list = ['Serie 1 (V1 bis V3)', None, None]
        color_list = ['k', 'r', 'b']
        linestyle_list = ['-', '-', '-']

    #---------------------------
    # tensile test results (CAR-800-TU-1v1l)-(TTb-2c-25cm-0-TU_TRC-pretests)
    #---------------------------
    if do == 'CAR-800-TU-1v1l_TRC-pretests':
        interpolated_flag = True
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.0  # correction factor for CAR-800-TU (a_tex*b = 53.9mm2/m*0.10m*6layers =32,34mm2; A_tex = 0.446mm2/Roving * 12 Rovings/layer * 6 layers=32,11mm2)

        # specify limits for the plot
        xlim = 9.
#        xlim = 18.
        if sig_flag == 'comp':
            ylim = 18.
        if sig_flag == 'tex':
            ylim = 1600.

        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-03-22_TTb_TRC-pretests', 'V12-1.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-03-22_TTb_TRC-pretests', 'V12-3.DAT')

#        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d', 'TTb-6c-2cm-0-TU-V1_bs4.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d', 'TTb-6c-2cm-0-TU-V2_bs4.DAT')
        path_V4 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d', 'TTb-6c-2cm-0-TU-V3_bs4.DAT')

        path_list = [path_V1, path_V2, path_V3, path_V4]
        label_list = [None, None, None, None]
        color_list = ['k', 'k', 'grey', 'grey']
        linestyle_list = ['-', '-', '-', '-']
        label_list = [r'fixe Endverankerung nach Bild 3.3', None, r'variable Endverankerung nach Bild 3.4', None]

    #---------------------------
    # tensile test results (CAR-800-TU-1v1l)-(TT-6c-2cm-0-TU_bs)
    #---------------------------
    if do == 'CAR-800-TU-1v1l_TT-2cm':
        interpolated_flag = True
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.007  # correction factor for CAR-800-TU (a_tex*b = 53.9mm2/m*0.10m*6layers =32,34mm2; A_tex = 0.446mm2/Roving * 12 Rovings/layer * 6 layers=32,11mm2)

        # specify limits for the plot
        xlim = 9.
        ylim = 30.
        #
        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-10-09_TT-6c-2cm-0-TU_bs', 'TT-6c-2cm-0-TU-V1.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-10-09_TT-6c-2cm-0-TU_bs', 'TT-6c-2cm-0-TU-V2.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-10-09_TT-6c-2cm-0-TU_bs', 'TT-6c-2cm-0-TU-V3.DAT')
        #
        path_V4 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-02-19_TT-6c-2cm-0-TU_bs', 'TT-6c-2cm-0-TU-V1.DAT')
        path_V5 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-02-19_TT-6c-2cm-0-TU_bs', 'TT-6c-2cm-0-TU-V2.DAT')
        path_V6 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2013-02-19_TT-6c-2cm-0-TU_bs', 'TT-6c-2cm-0-TU-V3.DAT')

        path_list = [path_V1, path_V2, path_V3,
                     path_V4, path_V5, path_V6]

        label_list = ['Serie 1 (V1 bis V3)', None, None,
                      'Serie 2 (V4 bis V6)', None, None]

        color_list = ['grey', 'grey', 'grey',
                      'k', 'k', 'k']

        linestyle_list = ['-', '-', '-',
                          '-', '-', '-']

    #---------------------------
    # tensile test results (CAR-3300-TR-3v3l) / (TT-6c-3cm-0-HT)
    #---------------------------
    if do == 'CAR-3300-TR-3v3l':
        interpolated_flag = False
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.0  # ignor correction factor for CAR-3300-TR-3v3l

        # specify limits for the plot
        xlim = 16.
        ylim = 9.
        #
        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-12-12_TT-6c-3cm-0-HT', 'TT-6c-3cm-0-HT-V1.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-12-12_TT-6c-3cm-0-HT', 'TT-6c-3cm-0-HT-V2.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-12-12_TT-6c-3cm-0-HT', 'TT-6c-3cm-0-HT-V3.DAT')
        #

        path_list = [path_V1, path_V2, path_V3]

        label_list = ['Serie 1 (V1 bis V3)', None, None]

        color_list = ['k', 'k', 'k']

        linestyle_list = ['-', '-', '-']

    #---------------------------
    # tensile test results (CAR-800-TR-2v1l) / (TT-12c-6cm-0-TR)
    #---------------------------
    if do == 'CAR-800-TR-2v1l':
        interpolated_flag = False
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        # correction factor for CAR-800-TR-2v1l (a_tex*b = 73.9mm2/m*0.14m*12layers =124mm2; A_tex = 0.446mm2/Roving * 22 Rovings/layer * 12 layers=117,7mm2)
        k_rho = 1.053

        # specify limits for the plot
        xlim = 18.
        ylim = 12.
        #
        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-03-22_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V1.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-03-22_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V2.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-03-22_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V3.DAT')
        #

        path_list = [path_V1, path_V2, path_V3]

        label_list = ['Serie 1 (V1 bis V3)', None, None]

        color_list = ['k', 'k', 'k']

        linestyle_list = ['-', '-', '-']

    #---------------------------
    # tensile test results (CAR-800-TR-1v1l) / (TT-12c-6cm-0-TR)
    #---------------------------
    if do == 'CAR-800-TR-1v1l':
        interpolated_flag = True
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        # correction factor for CAR-800-TR (a_tex*b = 53.9mm2/m*0.14m*12layers =91.6mm2; A_tex = 0.446mm2/Roving * 16 Rovings/layer * 12 layers=85,6mm2)
        k_rho = 1.057

        # specify limits for the plot
        xlim = 8.
        ylim = 12.
        #
        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-05-23_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V1.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-05-23_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V2.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-05-23_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V3.DAT')
        #

        path_list = [path_V1, path_V2, path_V3]

        label_list = ['Serie 1 (V1 bis V3)', None, None]

        color_list = ['k', 'k', 'k']

        linestyle_list = ['-', '-', '-']

    #---------------------------
    # tensile test results (CAR-800-TU-1v1l) / (TT-12c-6cm-0-TU-SH2)
    #---------------------------
    if do == 'CAR-800-TU-1v1l_TT-6cm':
        interpolated_flag = False
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        # correction factor for CAR-800-TU (a_tex*b = 53.9mm2/m*0.14m*12layers =90.5mm2; A_tex = 0.446mm2/Roving * 16 Rovings/layer * 12 layers=85,6mm2)
        k_rho = 1.057

        # specify limits for the plot
        xlim = 8.
        ylim = 18.
        #
        path_V1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V1.DAT')
        path_V2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V2.DAT')
        path_V3 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V3.DAT')
        #

        path_list = [path_V1, path_V2, path_V3]

        label_list = ['Serie 1 (V1 bis V3)', None, None]

        color_list = ['k', 'k', 'k']

        linestyle_list = ['-', '-', '-']

    #---------------------------
    # COMPARE (TT-12c-6cm-0-TR1-TR2-TU)
    #---------------------------
    if do == 'compare_CAR':
        save_table_to_file = False
        interpolated_flag = False
        stiffness_flag = False
        stiffness_flag_I = False
        stiffness_flag_II = True
        plot_sigtex_flag = False
        legend_flag = True
        legend_loc = 1
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.057  # correction factor

        # specify limits for the plot
        xlim = 18.
        if sig_flag == 'comp':
            ylim = 18.
        if sig_flag == 'tex':
            ylim = 1700.
        #
        # TU
        path_V1_TU = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V1.DAT')
        path_V2_TU = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V2.DAT')
        path_V3_TU = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-02-14_TT-12c-6cm-0-TU_SH2', 'TT-12c-6cm-0-TU-SH2-V3.DAT')

        # TR-1v1l
        path_V1_TR2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-05-23_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V1.DAT')
        path_V2_TR2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-05-23_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V2.DAT')
        path_V3_TR2 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-05-23_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V3.DAT')

        # TR-2v1l
        path_V1_TR1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-03-22_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V1.DAT')
        path_V2_TR1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-03-22_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V2.DAT')
        path_V3_TR1 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2011-03-22_TT-12c-6cm-0-TR', 'TT-12c-6cm-0-TR-V3.DAT')

        # HT
        path_V1_HT = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-12-12_TT-6c-3cm-0-HT', 'TT-6c-3cm-0-HT-V1.DAT')
        path_V2_HT = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-12-12_TT-6c-3cm-0-HT', 'TT-6c-3cm-0-HT-V2.DAT')
        path_V3_HT = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-12-12_TT-6c-3cm-0-HT', 'TT-6c-3cm-0-HT-V3.DAT')

        path_list = [path_V1_TU, path_V2_TU, path_V3_TU,
                     path_V1_TR2, path_V2_TR2, path_V3_TR2,
                     path_V1_TR1, path_V2_TR1, path_V3_TR1,
                     path_V1_HT, path_V2_HT, path_V3_HT
                    ]

        label_list = ['CAR-800-TU-1v1l', None, None,
                      'CAR-800-TR-1v1l', None, None,
                      'CAR-800-TR-2v1l', None, None,
                      'CAR-3300-TR-3v3l', None, None]

        color_list = ['k', 'k', 'k',
                      'grey', 'grey', 'grey',
                      'silver', 'silver', 'silver',
                      'k', 'k', 'k']

        linestyle_list = ['-', '-', '-',
                          '-', '-', '-',
                          '-', '-', '-',
                          ':', ':', ':']

    #---------------------------
    # COMPARE (ARG)
    #---------------------------
    if do == 'compare_ARG':
        save_table_to_file = False
        interpolated_flag = False
        stiffness_flag = False
        stiffness_flag_I = False
        stiffness_flag_II = True
        plot_sigtex_flag = False
        legend_flag = True
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.0  # ignor correction factor

        # specify limits for the plot
        xlim = 12.
        if sig_flag == 'comp':
            ylim = 14.
        if sig_flag == 'tex':
            ylim = 900.

        # FR
        path_V1_FR = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2009-11-26_TT-7g-3cm-0-FR_TRC07', 'TT-7g-3cm-0-FR-TRC07-V1.DAT')
        path_V2_FR = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2009-11-26_TT-7g-3cm-0-FR_TRC07', 'TT-7g-3cm-0-FR-TRC07-V2.DAT')
        path_V3_FR = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2009-11-26_TT-7g-3cm-0-FR_TRC07', 'TT-7g-3cm-0-FR-TRC07-V3.DAT')

        # TR24
        path_V1_TR24 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-11-02_TT-2g-1cm-0-TR_2400tex', 'TT-2g-1cm-0-TR-2400-V1.DAT')
        path_V2_TR24 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-11-02_TT-2g-1cm-0-TR_2400tex', 'TT-2g-1cm-0-TR-2400-V2.DAT')
        path_V3_TR24 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-11-02_TT-2g-1cm-0-TR_2400tex', 'TT-2g-1cm-0-TR-2400-V3.DAT')

        # TR12
        path_V1_TR12 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-03-17_TT-8g-3cm-0-TR', 'TT-8g-3cm-0-TR-V1.DAT')
        path_V2_TR12 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-03-17_TT-8g-3cm-0-TR', 'TT-8g-3cm-0-TR-V2.DAT')
        path_V3_TR12 = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2010-03-17_TT-8g-3cm-0-TR', 'TT-8g-3cm-0-TR-V3.DAT')

        # ARG-TU_TTb-2cm-bs4
        path_V1_TU = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V1_bs4.DAT')
        path_V2_TU = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V2_bs4.DAT')
        path_V3_TU = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V3_bs4.DAT')
#        # TU
#        path_V1_TU = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-V1.DAT')
#        path_V2_TU = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-V2.DAT')
#        path_V3_TU = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-12-10_TT-6g-2cm-0-TU_bs', 'TT-6g-2cm-0-V3.DAT')

        path_list = [path_V1_TU, path_V2_TU, path_V3_TU,
                     path_V2_TR12, path_V3_TR12,
                     path_V1_TR24, path_V2_TR24, path_V3_TR24,
                     path_V1_FR, path_V2_FR, path_V3_FR,
                    ]

        label_list = [
                      'ARG-1200-TU-1v1l', None, None,
                      'ARG-1200-TR-2v1l', None, None,
                      'ARG-2400-TR-1v1l', None, None,
                      'ARG-2400-FR-1v1l', None, None,
                      ]

        color_list = ['k', 'k', 'k',
                      'grey', 'grey',
                      'silver', 'silver', 'silver',
                      'k', 'k', 'k',
                      ]

        linestyle_list = ['-', '-', '-',
                          '-', '-',
                          '-', '-', '-',
                          ':', ':', ':',
                          ]

    #---------------------------
    # COMPARE (ARG-CAR-TU)
    #---------------------------
    if do == 'compare_ARG-CAR':
        save_table_to_file = False
        interpolated_flag = True
        stiffness_flag = False
        stiffness_flag_I = False
        stiffness_flag_II = True
        plot_sigtex_flag = False
        legend_flag = True
        # specify correction factor in order to calculate sigma_tex with the real reinforcement ratio
        # instead of the default value based on the specimen width multiplied with 'a_tex [mm2/m]'
        #
        k_rho = 1.0  # ignor correction factor as it is almost 1 for both 'ARG-1200-TU' and 'CAR-800-TU'

        # specify limits for the plot
        xlim = 12.
        if sig_flag == 'comp':
            ylim = 28.
        if sig_flag == 'tex':
            ylim = 1200.

        # CAR-TU_TTb-2cm-bs4
        path_V1_CAR = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d', 'TTb-6c-2cm-0-TU-V1_bs4.DAT')
        path_V2_CAR = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d', 'TTb-6c-2cm-0-TU-V2_bs4.DAT')
        path_V3_CAR = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6c-2cm-0-TU_bs4-Aramis3d', 'TTb-6c-2cm-0-TU-V3_bs4.DAT')

        # ARG-TU_TTb-2cm-bs4
        path_V1_ARG = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V1_bs4.DAT')
        path_V2_ARG = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V2_bs4.DAT')
        path_V3_ARG = os.path.join(simdb.exdata_dir, 'tensile_tests', 'buttstrap_clamping', '2013-07-09_TTb-6g-2cm-0-TU_bs4-Aramis3d', 'TTb-6g-2cm-0-TU-V3_bs4.DAT')
        #
        path_list = [path_V1_CAR, path_V2_CAR, path_V3_CAR,
                     path_V1_ARG, path_V2_ARG, path_V3_ARG]
        label_list = ['CAR-800-TU-1v1l', None, None,
                      'ARG-1200-TU-1v1l', None, None]
        color_list = ['k', 'k', 'k',
                      'grey', 'grey', 'grey']
        linestyle_list = ['-', '-', '-',
                          '-', '-', '-']

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
#        eps_u_list += [np.max(ex_run.ex_type.eps_c_interpolated) * 1000.]
        eps_u_list += [np.max(ex_run.ex_type.eps_asc) * 1000.]
#        sig_c_list += [np.max(ex_run.ex_type.sig_c_interpolated)]
        sig_c_list += [np.max(ex_run.ex_type.sig_c_asc)]
        # correct values for 'sig_tex' by real reinforcement ration determined by numer of rovings instead ob specimn width
        #
#        sig_tex_list += [np.max(ex_run.ex_type.sig_tex_interpolated) * k_rho]
        sig_tex_list += [np.max(ex_run.ex_type.sig_tex_asc) * k_rho]

    eps_u_arr = np.hstack(eps_u_list)
    sig_c_arr = np.hstack(sig_c_list)
    sig_tex_arr = np.hstack(sig_tex_list)
    print('eps_u_arr', eps_u_arr)
    print('sig_c_arr', sig_c_arr)
    print('sig_tex_arr', sig_tex_arr)

    # --------------------------------
    # calculate average, standard deviation and coreficient of variation
    # --------------------------------
    sig_c_avg = np.average(sig_c_arr)
    sig_tex_avg = np.average(sig_tex_arr)
    eps_u_avg = np.average(eps_u_arr)
    #
    sig_c_dev = np.std(sig_c_arr)
    sig_tex_dev = np.std(sig_tex_arr)
    eps_u_dev = np.std(eps_u_arr)
    #
    sig_c_Vx = sig_c_dev / sig_c_avg  # np.var(sig_c_arr)
    sig_tex_Vx = sig_tex_dev / sig_tex_avg
    eps_u_Vx = eps_u_dev / eps_u_avg
    # --------------------------------
    # plot sig-eps-curves
    # --------------------------------
    fig = p.figure(facecolor='white')
    fig.set_size_inches(8, 6)
#    fig.set_size_inches(12, 9)

    for i, ex_path in enumerate(path_list):
        if i <= 2:
            interpolated_flag = True
        else:
            interpolated_flag = False
        ex_run = ExRun(ex_path)
        if sig_flag == 'comp':
            ex_run.ex_type._plot_comp_stress_strain_asc(p, interpolated=interpolated_flag, k_rho=k_rho, color=color_list[i], linewidth=1.3, linestyle=linestyle_list[i], label=label_list[i], xscale=1000., plot_analytical_stiffness=stiffness_flag, plot_analytical_stiffness_I=stiffness_flag_I, plot_analytical_stiffness_II=stiffness_flag_II)
        if sig_flag == 'tex':
            ex_run.ex_type._plot_tex_stress_strain_asc(p, interpolated=True, k_rho=k_rho, color=color_list[i], linewidth=1.3, linestyle=linestyle_list[i], label=label_list[i], xscale=1000., plot_analytical_stiffness=stiffness_flag, plot_analytical_stiffness_I=stiffness_flag_I, plot_analytical_stiffness_II=stiffness_flag_II)

    # plot maximum average textile stress
    #
    if plot_sigtex_flag:
        xarr = np.array([0.4 * eps_u_avg, 1.2 * eps_u_avg])
        if sig_flag == 'comp':
            yarr = np.array([sig_c_avg, sig_c_avg])
        if sig_flag == 'tex':
            yarr = np.array([sig_tex_avg, sig_tex_avg])
        p.plot(xarr, yarr, linestyle='--', color='grey', linewidth=1.3)
        if sig_flag == 'comp':
            p.text(0.4 * eps_u_avg, sig_c_avg * 1.03, r'$\sigma_\mathrm{tex,u}\,=\,%0.0f\,\mathrm{MPa}$' % (sig_tex_avg), fontsize=16)  # , bbox={'facecolor':'white', 'edgecolor':'none'})
        if sig_flag == 'tex':
            p.text(0.4 * eps_u_avg, sig_tex_avg * 1.03, r'$\sigma_\mathrm{tex,u}\,=\,%0.0f\,\mathrm{MPa}$' % (sig_tex_avg), fontsize=16)  # , bbox={'facecolor':'white', 'edgecolor':'none'})
            # ## for figure compare pretest with barrelshell (TTb_2cm) only:
#            xarr = np.array([0.3 * eps_u_avg, 1.2 * eps_u_avg])
#            sig_tex_avg = 1077
#            yarr = np.array([sig_tex_avg, sig_tex_avg])
#            p.plot(xarr, yarr, linestyle='--', color='k', linewidth=1.)
#            p.text(0.3 * eps_u_avg, sig_tex_avg * 1.03, r'$\sigma_\mathrm{tex,u}\,=\,%0.0f\,\mathrm{MPa}$' % (sig_tex_avg), fontsize=16)  # , bbox={'facecolor':'white', 'edgecolor':'none'})
#            sig_tex_avg = 1483
#            yarr = np.array([sig_tex_avg, sig_tex_avg])
#            p.plot(xarr, yarr, linestyle='--', color='grey', linewidth=1.)
#            p.text(0.3 * eps_u_avg, sig_tex_avg * 1.03, r'$\sigma_\mathrm{tex,u}\,=\,%0.0f\,\mathrm{MPa}$' % (sig_tex_avg), color='grey', fontsize=16)  # , bbox={'facecolor':'white', 'edgecolor':'none'})
#            p.text(0.10 * eps_u_avg, sig_tex_avg * .04, r'$E_\mathrm{tex}$', color='k', fontsize=16)  # , bbox={'facecolor':'white', 'edgecolor':'none'})

    # format plot
    #
    if sig_flag == 'comp':
        format_plot(p, fontsize=15, xlabel='Dehnung $\epsilon$ [1E-3]', ylabel='Kompositspannung $\sigma_\mathrm{c}$ [MPa]', xlim=xlim, ylim=ylim)
    if sig_flag == 'tex':
        format_plot(p, fontsize=15, xlabel='Dehnung $\epsilon$ [1E-3]', ylabel='Textilspannung $\sigma_\mathrm{tex}$ [MPa]', xlim=xlim, ylim=ylim)
    axes = p.gca()
    axes.xaxis.grid(True, which='major')
    axes.yaxis.grid(True, which='major')

    # plot legend
    #
    if legend_flag:
        font = FontProperties()
        font.set_family('serif')
        font.set_style('normal')
        font.set_size(14)
        font.set_variant('normal')
#        font.set_weight('ultralight')
        #
        leg = p.legend(prop=font, loc=legend_loc)  # (loc=4: lower right) (7: center right)
        frame = leg.get_frame()
        frame.set_facecolor('white')
        frame.set_edgecolor('white')
        # set the linewidth of each legend object
        for legobj in leg.legendHandles:
            legobj.set_linewidth(1.50)

    # --------------------------------
    # save figure
    # --------------------------------
    if save_fig_to_file:
        # create a report-subfolder with the name of the script (without file extension '.py')
        # and save it in the test-type-subfolders with the name and path as ex_type
        test_series_name = os.path.basename(__file__)[:-3]
        subfolder_list = __file__.split(os.path.sep)
        devproc_idx = np.where(np.array(subfolder_list) == 'devproc')[0]
        subfolder_path = subfolder_list[devproc_idx + 1:-2] + [test_series_name]
        test_series_dir = os.path.join(simdb.report_dir)
        for subfolder_name in subfolder_path:
            test_series_dir = os.path.join(test_series_dir, subfolder_name)
        if not os.path.exists(test_series_dir):
            os.makedirs(test_series_dir)

        filename = os.path.join(test_series_dir, do + '_sig' + sig_flag + '-epsu.png')
        p.savefig(filename, format='png')
        print('figure saved to file %s' % (filename))

    p.show()

    # --------------------------------
    # save tabular values as Latex formated table
    # --------------------------------
    if save_table_to_file:
        test_series_name = os.path.basename(__file__)[:-3]
        subfolder_list = __file__.split(os.path.sep)
        devproc_idx = np.where(np.array(subfolder_list) == 'devproc')[0]
        subfolder_path = subfolder_list[devproc_idx + 1:-2] + [test_series_name]
        test_series_dir = os.path.join(simdb.report_dir)
        for subfolder_name in subfolder_path:
            test_series_dir = os.path.join(test_series_dir, subfolder_name)
        if not os.path.exists(test_series_dir):
            os.makedirs(test_series_dir)

        filename = os.path.join(test_series_dir, do + '_sig' + sig_flag + '-epsu.tex')
        f = open(filename, 'w')

        #---------------------
        # HORIZONTAL TABLE
        #---------------------
        print('sig_c_arr', sig_c_arr)
        if len(sig_c_arr) == 6:
            print('len(sig_c_arr):', len(sig_c_arr))
            output_str = ""
            output_str = output_str + r"\begin{tabularx}{\textwidth}{p{2,8cm} C{1cm} C{1cm} C{1cm} C{1cm} C{1cm} C{1cm} C{1cm} C{1cm} C{1cm}}" + "\n"
            output_str = output_str + r"\toprule" + "\n"
            output_str = output_str + r"{Versuchs Nr.:} & {V1} & {V2} & {V3} & {V4} & {V5} & {V6} & {($m$)} & {($s$)} & {($v$)} \\ " + "\n"
            output_str = output_str + r"\midrule" + "\n"

            output_str = output_str + r"{$\sigma_\textrm{c,u}\,[\textrm{MPa}]$:} & {%.1f} & {%.1f} & {%.1f} & {%.1f} & {%.1f} & {%.1f} & {(%.1f)} & {(%.1f)} & {(%.2f)}  \\" % (sig_c_arr[0], sig_c_arr[1], sig_c_arr[2], sig_c_arr[3], sig_c_arr[4], sig_c_arr[5], sig_c_avg, sig_c_dev, sig_c_Vx)
            output_str = output_str + "\n"
            output_str = output_str + r"{$\sigma_\textrm{tex,u}\,[\textrm{MPa}]$:} & {%.0f} & {%.0f} & {%.0f} & {%.0f} & {%.0f} & {%.0f} & {(%.0f)}  & {(%.0f)} & {(%.2f)} \\" % (sig_tex_arr[0], sig_tex_arr[1], sig_tex_arr[2], sig_tex_arr[3], sig_tex_arr[4], sig_tex_arr[5], sig_tex_avg, sig_tex_dev, sig_tex_Vx)
            output_str = output_str + "\n"
            output_str = output_str + r"{$\varepsilon_\textrm{u}\,[\textrm{\textperthousand}]$:}  & {%.1f} & {%.1f} & {%.1f} & {%.1f} & {%.1f}  & {%.1f} & {(%.1f)} & {(%.1f)} & {(%.2f)} \\" % (eps_u_arr[0], eps_u_arr[1], eps_u_arr[2], eps_u_arr[3], eps_u_arr[4], eps_u_arr[5], eps_u_avg, eps_u_dev, eps_u_Vx)
            output_str = output_str + "\n"
            output_str = output_str + r"\bottomrule" + "\n"
            output_str = output_str + r"\end{tabularx}" + "\n"

        if len(sig_c_arr) == 3:
            print('len(sig_c_arr):', len(sig_c_arr))
            output_str = ""
            output_str = output_str + r"\begin{tabularx}{\textwidth}{p{6cm} C{1.2cm} C{1.2cm} C{1.2cm} C{1.2cm} C{1.2cm} C{1.2cm}}" + "\n"
            output_str = output_str + r"\toprule" + "\n"
            output_str = output_str + r"{Versuchs Nr.:} & {V1} & {V2} & {V3} & {($m$)} & {($s$)} & {($v$)} \\ " + "\n"
            output_str = output_str + r"\midrule" + "\n"

            output_str = output_str + r"{Kompositspannung $\sigma_\textrm{c,u}\,[\textrm{MPa}]$:} & {%.1f} & {%.1f} & {%.1f} & {(%.1f)} & {(%.1f)} & {(%.2f)}  \\" % (sig_c_arr[0], sig_c_arr[1], sig_c_arr[2], sig_c_avg, sig_c_dev, sig_c_Vx)
            output_str = output_str + "\n"
            output_str = output_str + r"{Textilspannung $\sigma_\textrm{tex,u}\,[\textrm{MPa}]$:} & {%.0f} & {%.0f} & {%.0f} & {(%.0f)}  & {(%.0f)} & {(%.2f)} \\" % (sig_tex_arr[0], sig_tex_arr[1], sig_tex_arr[2], sig_tex_avg, sig_tex_dev, sig_tex_Vx)
            output_str = output_str + "\n"
            output_str = output_str + r"{Bruchdehnung $\varepsilon_\textrm{u}\,[\textrm{\textperthousand}]$:}  & {%.1f} & {%.1f} & {%.1f} & {(%.1f)} & {(%.1f)} & {(%.2f)} \\" % (eps_u_arr[0], eps_u_arr[1], eps_u_arr[2], eps_u_avg, eps_u_dev, eps_u_Vx)
            output_str = output_str + "\n"
            output_str = output_str + r"\bottomrule" + "\n"
            output_str = output_str + r"\end{tabularx}" + "\n"

        #---------------------
        # VERTICAL TABLE
        #---------------------
#        output_str = ""
#        output_str = output_str + r"\begin{tabular}{c c c c} " + "\n"
#        output_str = output_str + r"\toprule " + "\n"
#        output_str = output_str + r"\bf{Nr.} & \bf{$\sigma_\textrm{c}$} & \bf{$\sigma_\textrm{tex}$} & {$\varepsilon_\textrm{u}$} \\ " + "\n"
#        output_str = output_str + r"\midrule " + "\n"
#        #
#        for i in range(len(sig_c_arr)):
#            output_str = output_str + "{V" + str(i + 1) + "} & {%.1f} & {%.0f} & {%.1f}" % (sig_c_arr[i], sig_tex_arr[i], eps_u_arr[i])
#            output_str = output_str + r" \\ " + "\n"
#        output_str = output_str + r"\bottomrule " + "\n"
#        output_str = output_str + "{$m$} & {%.1f} & {%.0f} & {%.1f} " % (sig_c_avg, sig_tex_avg, eps_u_avg)
#        output_str = output_str + r" \\ " + "\n"
#        output_str = output_str + "{$s$} & {%.1f} & {%.0f} & {%.1f} " % (sig_c_dev, sig_tex_dev, eps_u_dev)
#        output_str = output_str + r" \\ " + "\n"
#        output_str = output_str + "{$v$} & {%.2f} & {%.2f} & {%.2f} " % (sig_c_Vx, sig_tex_Vx, eps_u_Vx)
#        output_str = output_str + r" \\ " + "\n"
#        output_str = output_str + r"\bottomrule " + "\n"
#        output_str = output_str + r"\end{tabular} " + "\n"

            print('output_str \n', output_str)

        f.write(output_str)
        print('table data saved to file %s' % (filename))


