
import numpy as np

if __name__ == '__main__':


    from matresdev.db.simdb import \
        SimDB

    import os

    from .lcc_table import LCCTableULS, LC

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    #---------------------------------------------
    # define filter:
    #---------------------------------------------

    # remove only the lowest point = connection shell/column
    # as this is a singularity of the FE-shell-model
    #                       cut_z_fraction = 0.05, # corresponds to 50cm x 50cm
    #                       cut_z_fraction = 0.10, # corresponds to 75cm x 75cm
    cut_z_fraction = 0.15,  # corresponds to 100cm x 100cm

    def remove_midpoints(lcc_table, arr):
        # remove the center points of the shells
        Zp = np.fabs(lcc_table.geo_data_orig['Z'])
        max_Z = max(Zp)
        min_Z = min(Zp)
        h = max_Z - min_Z
        z_active_idx = np.where(Zp >= (min_Z + cut_z_fraction * h))[0]
        return arr[ z_active_idx ]

    #---------------------------------------------
    # define data directory for evaluation of two shells with
    # new geometry and new loading cases (plus waterfilling)
    #---------------------------------------------

    data_dir = os.path.join(simdb.simdb_dir,
                            'simdata', 'input_data_mushroof_stb',
                            'ZiE_state_data_2shells_delta_h_865mm_2011-08-16')
#                            'state_data_2shells_delta_h_865mm_fine_mesh')

    #------------------------
    # define loading cases:
    #------------------------
    #
    #---------------------------------------------------------
    # "staendige und voruebergehende Bemessungssitauation":
    #---------------------------------------------------------
    lc_list = [
                 # LC1:
                 LC(name='g', category='dead-load', file_name='LC1.csv'
                    ),
                 # LC2:
                 LC(name='s_sym', category='imposed-load', file_name='LC2.csv',
                    exclusive_to=['s_asym', 'WF'],
                     psi_0=0.5, psi_1=0.2, psi_2=0.0
                     ),
                 # LC3:
                 LC(name='s_asym', category='imposed-load', file_name='LC3.csv',
                    exclusive_to=['s_sym', 'WF'],
                     psi_0=0.5, psi_1=0.2, psi_2=0.0
                     ),
                 # LC4:
                 LC(name='w_neg', category='imposed-load', file_name='LC4.csv',
                     exclusive_to=['w_pos', 'w_asym', 'w_int', 'WF'],
                     psi_0=0.6, psi_1=0.2, psi_2=0.0
                     ),
                 # LC5:
                 LC(name='w_pos', category='imposed-load', file_name='LC5.csv',
                     exclusive_to=['w_neg', 'w_asym', 'w_int', 'WF'],
                     psi_0=0.6, psi_1=0.2, psi_2=0.0
                     ),
                 # LC6:
                 # w_asym:
                 LC(name='w_asym', category='imposed-load', file_name='LC6.csv',
                     exclusive_to=['w_pos', 'w_neg', 'w_int', 'WF'],
                     psi_0=0.6, psi_1=0.2, psi_2=0.0
                     ),
                 # LC7:
                 LC(name='w_int', category='imposed-load', file_name='LC7.csv',
                     exclusive_to=['w_pos', 'w_neg', 'w_asym', 'WF', 'T_shrinkage'],
                     psi_0=0.6, psi_1=0.2, psi_2=0.0,
                     comment='Bauzustand'
                     ),
                 # LC8:
                 LC(name='Q_corner', category='imposed-load', file_name='LC8.csv',
                    exclusive_to=['Q_edge', 'WF'], psi_0=0.0, psi_1=0.0, psi_2=0.0,
                    comment='1 kN man load (corner)'
                    ),
                 # LC9:
                 LC(name='Q_edge', category='imposed-load', file_name='LC9.csv',
                     exclusive_to=['Q_corner', 'WF'], psi_0=0.0, psi_1=0.0, psi_2=0.0,
                     comment='1 kN man load (edge, center)'
                     ),
                 # LC10:
                 LC(name='T_pos', category='imposed-load', file_name='LC10.csv',
                     exclusive_to=['T_neg', 'WF'], psi_0=0.6, psi_1=0.5, psi_2=0.0,
                     comment='temperature (sommer)'
                     ),
                 # LC11:
                 LC(name='T_neg', category='imposed-load', file_name='LC11.csv',
                     exclusive_to=['T_pos', 'WF'], psi_0=0.6, psi_1=0.5, psi_2=0.0,
                     comment='temperature (winter)'
                     ),
                 # LC12:
                 LC(name='T_shrinkage', category='imposed-load', file_name='LC12.csv',
                     exclusive_to=['WF'], psi_0=0.8, psi_1=0.7, psi_2=0.5,
                     comment='shrinkage, combination coefficients taken from case "general imposed load"'
                     ),
               ]

#    #---------------------------------------------------------
#    # "aussergewoehnliche Bemessungssitauation":
#    #---------------------------------------------------------
#    #
#    lc_list = [
#                 # LC1:
#                 LC(name = 'g', category = 'dead-load',
#                     file_name = 'LC1.csv' ,
#                     gamma_unf = 1.00,
#                     comment = 'dead weight',
#                     ),
#                 # LC4:
#                 LC(name = 'w_neg', category = 'imposed-load',
#                     file_name = 'LC4.csv' ,
#                     exclusive_to = ['w_pos', 'w_asym', ],
#                     gamma_unf = 0.2, # = psi_1
#                     psi_0 = 0.0, # = psi_2
#                     ),
#                 # LC5:
#                 LC(name = 'w_pos', category = 'imposed-load',
#                     file_name = 'LC5.csv',
#                     exclusive_to = ['w_neg', 'w_asym'],
#                     gamma_unf = 0.2, # = psi_1
#                     psi_0 = 0.0      # = psi_2
#                     ),
#                 # LC6:
#                 LC(name = 'w_asym', category = 'imposed-load',
#                     file_name = 'LC6.csv',
#                     exclusive_to = ['w_pos', 'w_neg'],
#                     gamma_unf = 0.2, # = psi_1
#                     psi_0 = 0.0      # = psi_2
#                     ),
#                 # LC10:
#                 LC(name = 'T_pos', category = 'imposed-load',
#                     file_name = 'LC10.csv',
#                     exclusive_to = ['T_neg'],
#                     gamma_unf = 0.6, # = psi_0
#                     psi_0 = 0.5 / 0.6 # = psi_1
#                     ),
#                 # LC11:
#                 LC(name = 'T_neg', category = 'imposed-load',
#                     file_name = 'LC11.csv',
#                     exclusive_to = ['T_pos'],
#                     gamma_unf = 0.7, # = psi_1
#                     psi_0 = 0.5 / 0.7, # = psi_2
#                     comment = 'temperature (winter)'
#                     ),
#                 # LC12:
#                 LC(name = 'T_shrinkage', category = 'imposed-load',
#                     file_name = 'LC12.csv',
#                     gamma_unf = 0.8, # = psi_0
#                     psi_0 = 0.7 / 0.8, # = psi_1
#                     comment = 'shrinkage: combination coefficients taken from case "general imposed load"'
#                     ),
#                 # LC13:
#                 LC(name = 'WF', category = 'imposed-load',
#                     file_name = 'LC13.csv',
#                     gamma_unf = 1.00,
#                     psi_0 = 1.0,
#                     comment = 'water filling'
#                     ),
#                 # LC14_1-2:
#                 LC(name = 'EQ_1-2', category = 'imposed-load',
#                    file_name = 'LC14_1-2.csv',
#                    gamma_unf = 1.00,
#                    exclusive_to = ['EQ_3-4'],
#                    psi_0 = 1.0,
#                    comment = 'earth quake'
#                    ),
#                 # LC14_3-4:
#                 LC(name = 'EQ_3-4', category = 'imposed-load',
#                    file_name = 'LC14_3-4.csv',
#                    gamma_unf = 1.00,
#                    exclusive_to = ['EQ_1-2'],
#                    psi_0 = 1.0,
#                    comment = 'earthquake'
#                    ),
#               ]

#--------------------------------------------------------
# evaluation for eta_n-eta_m-interaction with varying angle of deflection
#--------------------------------------------------------

#    do = 'ZiE'
    do = 'QS'

    #--------------------------------------------------------
    # strength characteristics of the material:
    #--------------------------------------------------------
    # NOTE: the same strength in 90-direction has been assumed (safe side);

    if do == 'ZiE':
        #--------------------------------------------------------
        # (a) design values for SFB-demonstrator (ZiE)
        #--------------------------------------------------------
        # tensile strength [kN/m]
        #
        n_0_Rdt = n_90_Rdt = 412.  # = 57.7 (=F_tRd)/ 0.14 m ###  F_tm = 103.4 kN (mean value)

        # compressive strength [kN/m]
        #
        n_Rdc = 2200  # = ( 55 (=f_ck for C55/67) / 1.5 ) * 0.06 m * 1000 ### f_cm = 74.5 MPa (mean value)

        # bending strength [kNm/m]
        #
        m_0_Rd = m_90_Rd = 9.6  # = 1.93 (=M_Rd) / 0.20 m ### 3.5 kNm (mean value)

    if do == 'QS':
        #--------------------------------------------------------
        # (b) design values for SFB-demonstrator (QS)
        #--------------------------------------------------------
        # design values for quality tests of SFB-demonstrator (TT-SH2 and BT-SH4) on specimens with thickness 6 cm

        # tensile strength [kN/m]
        #
        # use strength reduction factor 'k_b' for discontinuities in the fabric layers due to butt joints:
        # i.e. k_b = 10/12 = 0.833333
        n_0_Rdt = n_90_Rdt = 0.833333 * 538.6  # = 75.4 (=F_tRd)/ 0.14 ### F_Rtm = 139.5 kN (mean value)

        # compressive strength [kN/m]
        #
        n_Rdc = 2200  # C55/67

        # bending strength [kNm/m]
        #
        m_0_Rd = m_90_Rd = 0.833333 * 8.3  # = 1.66 (=M_Rd) / 0.20 ### M_Rm = 3.1 kNm (mean value)

    #--------------------------------------------------------
    # ULS evaluation
    #--------------------------------------------------------
    lct = LCCTableULS(data_dir=data_dir,
                      data_filter=remove_midpoints,
                      lc_list=lc_list,
                      show_lc_characteristic=True,
                      strength_characteristics={'n_0_Rdt' : n_0_Rdt, 'm_0_Rd':m_0_Rd, 'n_Rdc' : n_Rdc,
                                                'n_90_Rdt' : n_90_Rdt, 'm_90_Rd':m_90_Rd},
                      k_alpha_min=False,  # NO simplification used for 'k_alpha' on the resistance side
                      )

    #--------------------------------------------------------------
    # 'combi_arr': array with indices of all loading case combinations
    #--------------------------------------------------------------
    #
    print('lct.combi_arr', lct.combi_arr.shape, '\n')
#     np.savetxt('combi_arr_LC1-12', lct.combi_arr, delimiter=';')

    #--------------------------------------------------------------
    # nm-interaction plot (normal force - bending moment)
    #--------------------------------------------------------------
    #
#    lct.plot_nm_interaction(save_fig_to_file='nm_interaction_LC1-12_' + do)

    #--------------------------------------------------------------
    # interaction plot of material usage 'eta_nm' (utilization ratio)
    #--------------------------------------------------------------
    #
#    lct.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_LC1-12_' + do)

    #--------------------------------------------------------------
    # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
    # (surrounding values of all loading cases)
    #--------------------------------------------------------------

#    lct.plot_assess_value('eta_nm_tot', scale_mode='scalar', scale_factor=0.8)

    #--------------------------------------------------------------
    # brows the loading case combinations within an interactive table view
    #--------------------------------------------------------------
    lct.configure_traits()

#--------------------------------------------------------
# SLS evaluation
#--------------------------------------------------------
#    lct = LCCTableSLS( data_dir = data_dir,
#                      data_filter = remove_midpoints,
#                       lc_list = lc_list,
#                       combination_SLS = 'rare',
# #                       combination_SLS = 'freq',
# #                       combination_SLS = 'perm',
# #                       show_lc_characteristic = True
#                        )
