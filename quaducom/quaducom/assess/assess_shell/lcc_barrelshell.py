if __name__ == '__main__':

    from matresdev.db.simdb import \
        SimDB

    import numpy as np

    import os

    from lcc_table import LCCTableULS, LC, LCCTableSLS

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    #------------------------
    # define filter:
    #------------------------

    def remove_support_elems(lcc_table, arr):
        '''filter for barrel shell used to remove the elements
        in direct contact to the support nodes
        (only valid for the specific discretization fines given,
        i.e. fine regular mesh (5cmx5cm-elements)
        as stored in '2cm-feines-Netz')
        '''
        elem_no = lcc_table.geo_data_orig['elem_no']

        # first row elements (direct contact to support node)
        #
        support_elems_list = [692, 693, 2324, 2325, 744, 745, 2376, 2377, \
                              # second row elements (= adjacent elements to first row elements)
                              #
                              691, 694, 2323, 2326, 743, 746, 2375, 2378]
        cond_arr = np.array(
            [elem_no != support_elem_number for support_elem_number in support_elems_list])
        cond_elem_active = np.product(cond_arr, axis=0)

        elem_active_idx = np.where(cond_elem_active)[0]
        if np.all(arr == lcc_table.geo_data_orig['t_elem_node_map'])\
                or np.all(arr == lcc_table.geo_data_orig['q_elem_node_map']) \
                or np.all(arr == lcc_table.geo_data_orig['t_idx'])\
                or np.all(arr == lcc_table.geo_data_orig['q_idx'])\
                or np.all(arr == lcc_table.lc_list[0].state_data_orig['ux'])\
                or np.all(arr == lcc_table.lc_list[0].state_data_orig['uy'])\
                or np.all(arr == lcc_table.lc_list[0].state_data_orig['uz'])\
                or np.all(arr == lcc_table.lc_list[0].state_data_orig['node_U']):
            return arr
        else:
            return arr[elem_active_idx]

    #------------------------
    # define directory:
    #------------------------

#    do = 'shell-test'
#    do = 'predimensioning'
    do = 'dimensioning'

    # specify weather to use strength characteristics of ZiE-test series or QS-test series
    #
    use_QS_values = True

    # experimental evaluation of shell strength
    # (surface load = 1kN/m^2,  along 10cm width, 30 cm free space from edges)
    #
    if do == 'shell-test':
        data_dir = os.path.join(simdb.simdb_dir,
                                'simdata',
                                'input_data_barrelshell',
                                '2cm-feines-Netz',
                                )

    # predimensioning of barrel shells for all main loading cases
    # (without superposition of temperature loading cases)
    #
    if do == 'predimensioning':
        data_dir = os.path.join(simdb.simdb_dir,
                                'simdata',
                                'input_data_barrelshell',
                                '2cm-feines-Netz',
                                )

    # dimensioning of barrel shells for all loading cases
    # (incl. superposition of temperature loading cases)
    #
    if do == 'dimensioning':
        data_dir = os.path.join(simdb.simdb_dir,
                                'simdata',
                                'input_data_barrelshell',
                                '2cm-feines-Netz-EC1',
                                )

    #------------------------
    # define loading cases:
    #------------------------

    #----------------------------------------------------------------------
    # (a) cable load (experimental test setup)
    #----------------------------------------------------------------------
    if do == 'shell-test':

        lc_list = [
            LC(name='cable-FL-m3', category='imposed-load', file_name='LC22.txt',
               # loading case only used for evaluation of experimental test setup
               # (estimation of load-bearing capacity)
               # scale load factor in order to get evaluation within
               # relevant values (>>0)
               gamma_unf=1000.0, gamma_fav=1.0,
               psi_0=1.0,
               ),
        ]

    #---------------------------------------------------------
    # (b) predimensioning
    #---------------------------------------------------------
    if do == 'predimensioning':

        lc_list = [

            #------------------------------------------------------------------
            # dead load
            #---------------------------------------------------------
            # LC1:
            LC(name='g', category='dead-load', file_name='LC1.txt'
               ),

            #---------------------------------------------------------
            # snow
            #---------------------------------------------------------

            # LC2:
            LC(name='s_hinten', category='imposed-load', file_name='LC2.txt',
               exclusive_to=[
                   's_feld', 's_vorne', 's_links', 's_rechts', 's_komplett'],
               psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),
            # LC3:
            LC(name='s_feld', category='imposed-load', file_name='LC3.txt',
               exclusive_to=[
                   's_hinten', 's_vorne', 's_links', 's_rechts', 's_komplett'],
               psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),
            # LC4:
            #                     LC(name='s_vorne', category='imposed-load', file_name='LC4.txt',
            #                        exclusive_to=['s_hinten', 's_feld', 's_links', 's_rechts', 's_komplett'],
            #                        psi_0=0.5, psi_1=0.2, psi_2=0.0
            #                        ),
            # LC5:
            LC(name='s_links', category='imposed-load', file_name='LC5.txt',
               exclusive_to=[
                   's_hinten', 's_feld', 's_vorne', 's_rechts', 's_komplett'],
               psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),
            # LC6:
            #                     LC(name='s_rechts', category='imposed-load', file_name='LC6.txt',
            #                        exclusive_to=['s_hinten', 's_feld', 's_vorne', 's_links', 's_komplett'],
            #                        psi_0=0.5, psi_1=0.2, psi_2=0.0
            #                        ),
            # LC7:
            LC(name='s_komplett', category='imposed-load', file_name='LC7.txt',
               exclusive_to=[
                   's_hinten', 's_feld', 's_vorne', 's_links', 's_rechts'],
               psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),

            #---------------------------------------------------------
            # man load (1 KN)
            #---------------------------------------------------------
            #                     # LC8:
            #                     LC(name = 'Q_hinten_mitte', category = 'imposed-load', file_name = 'LC8.txt',
            #                        exclusive_to = ['Q_feld_mitte', 'Q_feld_li','s_vorne', 's_hinten', 's_feld', 's_komplett', 's_links', 's_rechts','T_schwinden','w_vonlinks_komplett','w_vonrechts_komplett','w_druck_komplett','w_sog_komplett'],
            #                        psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
            #                        ),
            #                     # LC9:
            #                     LC(name = 'Q_feld_li', category = 'imposed-load', file_name = 'LC9.txt',
            #                        exclusive_to = ['Q_feld_mitte', 'Q_hinten_mitte','s_vorne', 's_hinten', 's_feld', 's_komplett', 's_links', 's_rechts','T_schwinden','w_vonlinks_komplett','w_vonrechts_komplett','w_druck_komplett','w_sog_komplett'],
            #                        psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
            #                        ),
            #                     # LC10:
            #                     LC(name = 'Q_feld_mitte', category = 'imposed-load', file_name = 'LC10.txt',
            #                        exclusive_to = ['Q_hinten_mitte', 'Q_feld_li','s_vorne', 's_hinten', 's_feld', 's_komplett', 's_links', 's_rechts','T_schwinden','w_vonlinks_komplett','w_vonrechts_komplett','w_druck_komplett','w_sog_komplett'],
            #                        psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
            #                        ),
            #                     # LC11:
            #                     LC(name = 'Q_feld_re', category = 'imposed-load', file_name = 'LC11.txt',
            #                        exclusive_to = ['Q_hinten_mitte', 'Q_feld_li', 'Q_feld_mitte', 'Q_vorne_mitte'],
            #                        psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
            #                        ),
            #                     # LC12:
            #                     LC(name = 'Q_vorne_mitte', category = 'imposed-load', file_name = 'LC12.txt',
            #                        exclusive_to = ['Q_hinten_mitte', 'Q_feld_li', 'Q_feld_mitte', 'Q_feld_re'],
            #                        psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
            #                        ),

            #---------------------------------------------------------
            # temperature
            #---------------------------------------------------------

            #                     # LC13:
            #                     LC(name = 'T_N_neg', category = 'imposed-load', file_name = 'LC13.txt',
            #                        exclusive_to = ['T_N_pos', 'T_uo_neg', 'T_uo_pos'],
            #                        psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0
            #                        ),
            #                     # LC14:
            #                     LC(name = 'T_N_pos', category = 'imposed-load', file_name = 'LC14.txt',
            #                        exclusive_to = ['T_N_neg', 'T_uo_neg', 'T_uo_pos'],
            #                        psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0
            #                        ),
            #                     # LC15:
            #                     LC(name = 'T_uo_neg', category = 'imposed-load', file_name = 'LC15.txt',
            #                        exclusive_to = ['T_N_neg', 'T_N_pos', 'T_uo_pos'],
            #                        psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0
            #                        ),
            #                     # LC16:
            #                     LC(name = 'T_uo_pos', category = 'imposed-load', file_name = 'LC16.txt',
            #                        exclusive_to = ['T_N_neg', 'T_N_pos', 'T_uo_neg'],
            #                        psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0
            #                        ),

            #---------------------------------------------------------
            # shrinkage
            #---------------------------------------------------------
            # LC17:
            #                     LC(name='T_schwinden', category='imposed-load', file_name='LC17.txt',
            #                        psi_0=0.8, psi_1=0.7, psi_2=0.5,
            #                        ),

            #---------------------------------------------------------
            # wind load
            #---------------------------------------------------------
            # LC42:
            LC(name='w_vonlinks_komplett', category='imposed-load', file_name='LC42.txt',
               exclusive_to=[
                   'w_vonrechts_komplett', 'w_sog_komplett', 'w_druck_komplett'],
               psi_0=0.6, psi_1=0.2, psi_2=0.0
               ),
            # LC43:
            #                     LC(name='w_vonrechts_komplett', category='imposed-load', file_name='LC43.txt',
            #                        exclusive_to=['w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett'],
            #                        psi_0=0.6, psi_1=0.2, psi_2=0.0
            #                        ),
            # LC44:
            LC(name='w_sog_komplett', category='imposed-load', file_name='LC44.txt',
               exclusive_to=[
                   'w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_druck_komplett'],
               psi_0=0.6, psi_1=0.2, psi_2=0.0
               ),
            # LC45:
            LC(name='w_druck_komplett', category='imposed-load', file_name='LC45.txt',
               exclusive_to=[
                   'w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_sog_komplett'],
               psi_0=0.6, psi_1=0.2, psi_2=0.0
               ),
        ]

    #----------------------------------------------------------------------
    # (c) dimensioning with loads according to DIN EN 1991-1-3 (snow) and DIN EN 1991-1-4 (wind)
    #----------------------------------------------------------------------

    if do == 'dimensioning':

        # NOTE: 'lc_list_Q' contains all imposed loading cases (except for temperature)
        # the loading cases for temperature are evaluated separately (see 'lc_list_T' below)
        # and superposed later
        #
        lc_list_Q = [

            #------------------------------------------------------------------
            # dead load
            #-------------------------------------------------------------
            # LC1:
            LC(name='g', category='dead-load', file_name='LC1.txt'
               ),

            #-------------------------------------------------------------
            # snow
            #-------------------------------------------------------------

            # LC2:
            LC(name='s_komplett', category='imposed-load', file_name='LC2.txt',
                    exclusive_to=[
                        's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld'],
                    psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),
            #                 # LC3:
            #                 LC(name='s_verweht_re', category='imposed-load', file_name='LC3.txt',
            #                    exclusive_to=['s_komplett', 's_scheddach_re', 's_hinten', 's_feld'],
            #                    psi_0=0.5, psi_1=0.2, psi_2=0.0
            #                    ),
            # LC4:
            LC(name='s_scheddach_re', category='imposed-load', file_name='LC4.txt',
                    exclusive_to=[
                        's_komplett', 's_verweht_re', 's_hinten', 's_feld'],
                    psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),
            #                 # LC5:
            #                 LC(name='s_hinten', category='imposed-load', file_name='LC5.txt',
            #                    exclusive_to=['s_komplett', 's_verweht_re', 's_scheddach_re'],
            #                    psi_0=0.5, psi_1=0.2, psi_2=0.0
            #                    ),
            #                 # LC6:
            #                 LC(name='s_feld', category='imposed-load', file_name='LC6.txt',
            #                    exclusive_to=['s_komplett', 's_verweht_re', 's_scheddach_re'],
            #                    psi_0=0.5, psi_1=0.2, psi_2=0.0
            #                    ),

            #-------------------------------------------------------------
            # wind load
            #-------------------------------------------------------------
            #                 # LC7:
            #                 LC(name='w_sog_komplett', category='imposed-load', file_name='LC7.txt',
            #                    exclusive_to=['w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_druck_komplett'],
            #                    psi_0=0.6, psi_1=0.2, psi_2=0.0
            #                    ),
            # LC8:
            LC(name='w_druck_komplett', category='imposed-load', file_name='LC8.txt',
                    exclusive_to=[
                        'w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_sog_komplett'],
                    psi_0=0.6, psi_1=0.2, psi_2=0.0
               ),
            #                 # LC9:
            #                 LC(name='w_vonlinks_komplett', category='imposed-load', file_name='LC9.txt',
            #                    exclusive_to=['w_vonrechts_komplett', 'w_sog_komplett', 'w_druck_komplett'],
            #                    psi_0=0.6, psi_1=0.2, psi_2=0.0
            #                    ),
            # LC10:
            LC(name='w_vonrechts_komplett', category='imposed-load', file_name='LC10.txt',
                    exclusive_to=[
                        'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett'],
                    psi_0=0.6, psi_1=0.2, psi_2=0.0
               ),

            #-------------------------------------------------------------
            # man load (1 KN)
            #-------------------------------------------------------------
            #                 # LC11:
            #                 LC(name='Q_hinten', category='imposed-load', file_name='LC11.txt',
            #                    exclusive_to=[ 'Q_links', 'Q_mitte', 's_komplett', 's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld', 'w_vonrechts_komplett', 'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett', 'T_schwinden'],
            #                    psi_0=0.0, psi_1=0.2, psi_2=0.0
            #                    ),
            #                 # LC12:
            #                 LC(name='Q_mitte', category='imposed-load', file_name='LC12.txt',
            #                    exclusive_to=[ 'Q_links', 'Q_hinten', 's_komplett', 's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld', 'w_vonrechts_komplett', 'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett', 'T_schwinden'],
            #                    psi_0=0.0, psi_1=0.2, psi_2=0.0
            #                    ),
            #                 # LC13:
            #                 LC(name='Q_links', category='imposed-load', file_name='LC13.txt',
            #                    exclusive_to=[ 'Q_hinten', 'Q_mitte', 's_komplett', 's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld', 'w_vonrechts_komplett', 'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett', 'T_schwinden'],
            #                    psi_0=0.0, psi_1=0.2, psi_2=0.0
            #                    ),
            #-------------------------------------------------------------
            # shrinkage
            #-------------------------------------------------------------
            # LC14:
            LC(name='T_schwinden', category='imposed-load', file_name='LC14.txt',
                    psi_0=0.8, psi_1=0.7, psi_2=0.5,
               ),

            # LC15:
            LC(name='T_N_pos', category='imposed-load', file_name='LC15.txt',
                    exclusive_to=['T_N_neg', 'T_uo_pos', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC16:
            LC(name='T_N_neg', category='imposed-load', file_name='LC16.txt',
                    exclusive_to=['T_N_pos', 'T_uo_pos', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC17:
            LC(name='T_uo_pos', category='imposed-load', file_name='LC17.txt',
                    exclusive_to=['T_N_pos', 'T_N_neg', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC18:
            LC(name='T_uo_neg', category='imposed-load', file_name='LC18.txt',
                    exclusive_to=['T_N_pos', 'T_N_neg', 'T_uo_pos'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),

        ]

        lc_list_T = [

            #------------------------------------------------------------------
            # temperature (superpose separately = corresponds to conservatively setting psi_0 = 1.0)
            #-------------------------------------------------------------

            # LC15:
            LC(name='T_N_pos', category='imposed-load', file_name='LC15.txt',
                    exclusive_to=['T_N_neg', 'T_uo_pos', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC16:
            LC(name='T_N_neg', category='imposed-load', file_name='LC16.txt',
                    exclusive_to=['T_N_pos', 'T_uo_pos', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC17:
            LC(name='T_uo_pos', category='imposed-load', file_name='LC17.txt',
                    exclusive_to=['T_N_pos', 'T_N_neg', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC18:
            LC(name='T_uo_neg', category='imposed-load', file_name='LC18.txt',
                    exclusive_to=['T_N_pos', 'T_N_neg', 'T_uo_pos'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
        ]

    #--------------------------------------------------------------
    # ULS evaluation
    #--------------------------------------------------------------

    if do == 'shell-test':
        #--------------------------------------------------------
        # strength characteristics for experimental (mean) values for barrelshell
        # (specimens thickness = 2 cm; specimn width = 10 cm; 6 layers carbon)
        #--------------------------------------------------------
        if 'use_QS_values':
            # design values for quality tests of barrelshell (TTb-bs4 and
            # BT_bs4) on specimens with thickness 2 cm

            # tensile strength [kN/m]
            #
            # = 25.7 (=F_tRd)/ 0.10 ### F_Rtm = 47.6 kN (mean value)
            n_0_Rdt = n_90_Rdt = 476.

            # compressive strength [kN/m]
            #
            n_Rdc = 1360  # f_cm=68 MPa

            # bending strength [kNm/m]
            #
            # = 0.18 (=M_Rd) / 0.10 ### M_Rm = 0.33 kNm (mean value)
            m_0_Rd = m_90_Rd = 3.3

            print 'mean values (QS-bs4) used for strength characteristics  (no k_b; no reduction for scatter; no gamma_M)'

        else:
            # tensile strength [kN/m]
            #
            n_0_Rdt = n_90_Rdt = 41.1 / 0.10  # [kN/m]

            # bending strength [kNm/m]
            #
            m_0_Rd = m_90_Rd = (3.5 * 0.46 / 4.) / 0.10  # [kNm/m]

            # compressive strength [kN/m]
            # (design value; f_cd = 37,5 MPa)
            #
            n_Rdc = 750.  # = 37,5 MPa * (100 cm * 2 cm) * 0.1

        lct = LCCTableULS(data_dir=data_dir,
                          reader_type='InfoCAD',
                          #                          data_filter=remove_support_elems,
                          lc_list=lc_list,
                          strength_characteristics={'n_0_Rdt': n_0_Rdt,
                                                    'm_0_Rd': m_0_Rd,
                                                    'n_Rdc': n_Rdc,
                                                    'n_90_Rdt': n_90_Rdt,
                                                    'm_90_Rd': m_90_Rd},
                          # NO simplification used for 'k_alpha' on the
                          # resistance side
                          k_alpha_min=False,
                          show_lc_characteristic=True
                          )

        #--------------------------------------------------------------
        # 'combi_arr': array with indices of all loading case combinations
        #--------------------------------------------------------------
        #
#        print 'lct.combi_arr', lct.combi_arr.shape
#        np.savetxt('combi_arr_wo_temp_LCs', lct.combi_arr, delimiter=';')

        #--------------------------------------------------------------
        # nm-interaction plot (normal force - bending moment)
        #--------------------------------------------------------------
        #
        lct.plot_nm_interaction(save_fig_to_file='nm_interaction_shell-test')

        #--------------------------------------------------------------
        # interaction plot of material usage 'eta_nm' (utilization ratio)
        #--------------------------------------------------------------
        #
#        lct.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_shell-test')

        #--------------------------------------------------------------
        # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
        # (surrounding values of all loading cases)
        #--------------------------------------------------------------
        #
        lct.plot_assess_value('eta_nm_tot',
                              save_fig_to_file='eta_nm_tot_shell-test')

        #--------------------------------------------------------------
        # brows the loading case combinations within an interactive table view
        #--------------------------------------------------------------
        lct.configure_traits()

    if do == 'predimensioning':
        #--------------------------------------------------------
        # strength characteristics (design) values for barrelshell
        # (specimens thickness = 2 cm; specimn width = 10 cm; 6 layers carbon)
        #--------------------------------------------------------

        # tensile strength [kN/m]
        #
        n_0_Rdt = n_90_Rdt = 223.  # [kN/m] # ZiE value

        # bending strength [kNm/m]
        #
        m_0_Rd = m_90_Rd = 1.7  # [kNm/m] # ZiE value

        # compressive strength [kN/m]
        # (design value; f_cd = 37,5 MPa)
        #
        n_Rdc = 750.  # = 37,5 MPa * (100 cm * 2 cm) * 0.1

        print 'design values used for strength characteristics'

        lct = LCCTableULS(data_dir=data_dir,
                          reader_type='InfoCAD',
                          data_filter=remove_support_elems,
                          lc_list=lc_list,
                          strength_characteristics={'n_0_Rdt': n_0_Rdt, 'm_0_Rd': m_0_Rd, 'n_Rdc': n_Rdc,
                                                    'n_90_Rdt': n_90_Rdt, 'm_90_Rd': m_90_Rd},
                          # NO simplification used for 'k_alpha' on the
                          # resistance side
                          k_alpha_min=True,
                          show_lc_characteristic=False
                          )

        #--------------------------------------------------------------
        # 'combi_arr': array with indices of all loading case combinations
        #--------------------------------------------------------------
        #
        print 'lct.combi_arr', lct.combi_arr.shape
#        np.savetxt('combi_arr_wo_temp_LCs', lct_Q.combi_arr, delimiter = ';')

        #--------------------------------------------------------------
        # nm-interaction plot (normal force - bending moment)
        #--------------------------------------------------------------
        #
        lct.plot_nm_interaction(save_fig_to_file='nm_interaction_LC22')

        #--------------------------------------------------------------
        # interaction plot of material usage 'eta_nm' (utilization ratio)
        #--------------------------------------------------------------
        #
        lct.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_LC22')

        #--------------------------------------------------------------
        # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
        # (surrounding values of all loading cases)
        #--------------------------------------------------------------
        #
        lct.plot_assess_value('eta_nm_tot')

        #--------------------------------------------------------------
        # brows the loading case combinations within an interactive table view
        #--------------------------------------------------------------
        lct.configure_traits()

    if do == 'dimensioning':
        #--------------------------------------------------------
        # strength characteristics (design) values for barrelshell
        # (specimens thickness = 2 cm; specimn width = 10 cm; 6 layers carbon)
        #--------------------------------------------------------
        if 'use_QS_values':
            # design values for quality tests of barrelshell (TTb-bs4 and
            # BT_bs4) on specimens with thickness 2 cm

            # tensile strength [kN/m]
            #
            # use strength reduction factor 'k_b' for discontinuities in the fabric layers due to butt joints
            # NOTE: (omit as no butt joints are located in midfield):
            # i.e. k_b = 5/6 = 0.833333
            # = 25.7 (=F_tRd)/ 0.10 ### F_Rtm = 47.6 kN (mean value)
            n_0_Rdt = n_90_Rdt = 0.83333 * 257

            # compressive strength [kN/m]
            #
            n_Rdc = 750  # f_cd=37.5 MPa

            # bending strength [kNm/m]
            #
            # = 0.18 (=M_Rd) / 0.10 ### M_Rm = 0.33 kNm (mean value)
            m_0_Rd = m_90_Rd = 0.83333 * 1.8

            print 'design values (QS-bs4) used for strength characteristics'

        else:
            # tensile strength [kN/m]
            #
            n_0_Rdt = n_90_Rdt = 223.  # [kN/m] # ZiE value

            # bending strength [kNm/m]
            #
            m_0_Rd = m_90_Rd = 1.7  # [kNm/m] # ZiE value

            # compressive strength [kN/m]
            # (design value; f_cd = 37,5 MPa)
            #
            n_Rdc = 750.  # = 37,5 MPa * (100 cm * 2 cm) * 0.1

            print 'design values (ZiE) used for strength characteristics'

        # LCCTable for imposed loads (without temperature)
        #
        lct_Q = LCCTableULS(data_dir=data_dir,
                            reader_type='InfoCAD',
                            data_filter=remove_support_elems,
                            lc_list=lc_list_Q,
                            strength_characteristics={'n_0_Rdt': n_0_Rdt, 'm_0_Rd': m_0_Rd, 'n_Rdc': n_Rdc,
                                                      'n_90_Rdt': n_90_Rdt, 'm_90_Rd': m_90_Rd},
                            # simplification: use the minimum value for k_alpha
                            # on the resistance side
                            k_alpha_min=False,
                            #                          show_lc_characteristic = True
                            )

        # LCCTable for temperature loading cases only
        #
        lct_T = LCCTableULS(data_dir=data_dir,
                            reader_type='InfoCAD',
                            data_filter=remove_support_elems,
                            lc_list=lc_list_T,
                            strength_characteristics={'n_0_Rdt': n_0_Rdt, 'm_0_Rd': m_0_Rd, 'n_Rdc': n_Rdc,
                                                      'n_90_Rdt': n_90_Rdt, 'm_90_Rd': m_90_Rd},
                            # k_alpha_min=True,  # simplification: use the
                            # minimum value for k_alpha on the resistance side
                            # NOTE: no simplification: use calculated values
                            # for k_alpha on the resistance side
                            k_alpha_min=False,
                            #       and evaluation in 'alpha_1' AND 'alpha_2' direction
                            #                          show_lc_characteristic = True
                            )

        #--------------------------------------------------------------
        # 'combi_arr': array with indices of all loading case combinations
        #--------------------------------------------------------------
        #
        print 'lct_Q.combi_arr', lct_Q.combi_arr.shape, '\n'
#        np.savetxt('combi_arr_wo_temp_LCs', lct_Q.combi_arr, delimiter=';')

        #--------------------------------------------------------------
        # nm-interaction plot (normal force - bending moment)
        #--------------------------------------------------------------
        #
#        lct_T.plot_nm_interaction(save_max_min_nm_to_file='max_min_nm_arr_LC15-18', save_fig_to_file='nm_interaction_LC15-18')
#        lct_Q.plot_nm_interaction(save_fig_to_file='nm_interaction_LC1-14')
#        lct_Q.plot_nm_interaction(add_max_min_nm_from_file='max_min_nm_arr_LC15-18', save_fig_to_file='nm_interaction_LC1-18')

        #--------------------------------------------------------------
        # interaction plot of material usage 'eta_nm' (utilization ratio)
        #--------------------------------------------------------------
        #
#        lct_T.plot_eta_nm_interaction(save_max_min_eta_nm_to_file='max_min_eta_nm_arr_LC15-18', save_fig_to_file='eta_nm_interaction_LC15-18')
#        lct_Q.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_LC1-14')
#        lct_Q.plot_eta_nm_interaction(add_max_min_eta_nm_from_file='max_min_eta_nm_arr_LC15-18', save_fig_to_file='eta_nm_interaction_LC1-18')

        #--------------------------------------------------------------
        # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
        # (surrounding values of all loading cases)
        #--------------------------------------------------------------
        #
#        lct_T.plot_assess_value('eta_nm_tot', save_assess_values_to_file='eta_nm_tot_LC15-18')
#        lct_Q.plot_assess_value('eta_nm_tot', add_assess_values_from_file='eta_nm_tot_LC15-18')

        #--------------------------------------------------------------
        # brows the loading case combinations within an interactive table view
        #--------------------------------------------------------------
        lct_Q.configure_traits()
#        lct_T.configure_traits()
=======
if __name__ == '__main__':

    from matresdev.db.simdb import \
        SimDB

    import numpy as np

    import os

    from lcc_table import LCCTableULS, LC, LCCTableSLS

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    #------------------------
    # define filter:
    #------------------------

    def remove_support_elems(lcc_table, arr):
        '''filter for barrel shell used to remove the elements
        in direct contact to the support nodes
        (only valid for the specific discretization fines given,
        i.e. fine regular mesh (5cmx5cm-elements)
        as stored in '2cm-feines-Netz')
        '''
        elem_no = lcc_table.geo_data_orig['elem_no']

        # first row elements (direct contact to support node)
        #
        support_elems_list = [692, 693, 2324, 2325, 744, 745, 2376, 2377, \
                              # second row elements (= adjacent elements to first row elements)
                              #
                              691, 694, 2323, 2326, 743, 746, 2375, 2378]
        cond_arr = np.array(
            [elem_no != support_elem_number for support_elem_number in support_elems_list])
        cond_elem_active = np.product(cond_arr, axis=0)

        elem_active_idx = np.where(cond_elem_active)[0]
        if np.all(arr == lcc_table.geo_data_orig['t_elem_node_map'])\
                or np.all(arr == lcc_table.geo_data_orig['q_elem_node_map']) \
                or np.all(arr == lcc_table.geo_data_orig['t_idx'])\
                or np.all(arr == lcc_table.geo_data_orig['q_idx'])\
                or np.all(arr == lcc_table.lc_list[0].state_data_orig['ux'])\
                or np.all(arr == lcc_table.lc_list[0].state_data_orig['uy'])\
                or np.all(arr == lcc_table.lc_list[0].state_data_orig['uz'])\
                or np.all(arr == lcc_table.lc_list[0].state_data_orig['node_U']):
            return arr
        else:
            return arr[elem_active_idx]

    #------------------------
    # define directory:
    #------------------------

    do = 'shell-test'
#    do = 'predimensioning'
#    do = 'dimensioning'

    # specify weather to use strength characteristics of ZiE-test series or QS-test series
    #
    use_QS_values = True

    # experimental evaluation of shell strength
    # (surface load = 1kN/m^2,  along 10cm width, 30 cm free space from edges)
    #
    if do == 'shell-test':
        data_dir = os.path.join(simdb.simdb_dir,
                                'simdata',
                                'input_data_barrelshell',
                                '2cm-feines-Netz',
                                )

    # predimensioning of barrel shells for all main loading cases
    # (without superposition of temperature loading cases)
    #
    if do == 'predimensioning':
        data_dir = os.path.join(simdb.simdb_dir,
                                'simdata',
                                'input_data_barrelshell',
                                '2cm-feines-Netz',
                                )

    # dimensioning of barrel shells for all loading cases
    # (incl. superposition of temperature loading cases)
    #
    if do == 'dimensioning':
        data_dir = os.path.join(simdb.simdb_dir,
                                'simdata',
                                'input_data_barrelshell',
                                '2cm-feines-Netz-EC1',
                                )

    #------------------------
    # define loading cases:
    #------------------------

    #----------------------------------------------------------------------
    # (a) cable load (experimental test setup)
    #----------------------------------------------------------------------
    if do == 'shell-test':

        lc_list = [
            LC(name='cable-FL-m3', category='imposed-load', file_name='LC22.txt',
               # loading case only used for evaluation of experimental test setup
               # (estimation of load-bearing capacity)
               # scale load factor in order to get evaluation within
               # relevant values (>>0)
               gamma_unf=100.0, gamma_fav=1.0,
               psi_0=1.0,
               ),
        ]

    #---------------------------------------------------------
    # (b) predimensioning
    #---------------------------------------------------------
    if do == 'predimensioning':

        lc_list = [

            #------------------------------------------------------------------
            # dead load
            #---------------------------------------------------------
            # LC1:
            LC(name='g', category='dead-load', file_name='LC1.txt'
               ),

            #---------------------------------------------------------
            # snow
            #---------------------------------------------------------

            # LC2:
            LC(name='s_hinten', category='imposed-load', file_name='LC2.txt',
               exclusive_to=[
                   's_feld', 's_vorne', 's_links', 's_rechts', 's_komplett'],
               psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),
            # LC3:
            LC(name='s_feld', category='imposed-load', file_name='LC3.txt',
               exclusive_to=[
                   's_hinten', 's_vorne', 's_links', 's_rechts', 's_komplett'],
               psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),
            # LC4:
            #                     LC(name='s_vorne', category='imposed-load', file_name='LC4.txt',
            #                        exclusive_to=['s_hinten', 's_feld', 's_links', 's_rechts', 's_komplett'],
            #                        psi_0=0.5, psi_1=0.2, psi_2=0.0
            #                        ),
            # LC5:
            LC(name='s_links', category='imposed-load', file_name='LC5.txt',
               exclusive_to=[
                   's_hinten', 's_feld', 's_vorne', 's_rechts', 's_komplett'],
               psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),
            # LC6:
            #                     LC(name='s_rechts', category='imposed-load', file_name='LC6.txt',
            #                        exclusive_to=['s_hinten', 's_feld', 's_vorne', 's_links', 's_komplett'],
            #                        psi_0=0.5, psi_1=0.2, psi_2=0.0
            #                        ),
            # LC7:
            LC(name='s_komplett', category='imposed-load', file_name='LC7.txt',
               exclusive_to=[
                   's_hinten', 's_feld', 's_vorne', 's_links', 's_rechts'],
               psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),

            #---------------------------------------------------------
            # man load (1 KN)
            #---------------------------------------------------------
            #                     # LC8:
            #                     LC(name = 'Q_hinten_mitte', category = 'imposed-load', file_name = 'LC8.txt',
            #                        exclusive_to = ['Q_feld_mitte', 'Q_feld_li','s_vorne', 's_hinten', 's_feld', 's_komplett', 's_links', 's_rechts','T_schwinden','w_vonlinks_komplett','w_vonrechts_komplett','w_druck_komplett','w_sog_komplett'],
            #                        psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
            #                        ),
            #                     # LC9:
            #                     LC(name = 'Q_feld_li', category = 'imposed-load', file_name = 'LC9.txt',
            #                        exclusive_to = ['Q_feld_mitte', 'Q_hinten_mitte','s_vorne', 's_hinten', 's_feld', 's_komplett', 's_links', 's_rechts','T_schwinden','w_vonlinks_komplett','w_vonrechts_komplett','w_druck_komplett','w_sog_komplett'],
            #                        psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
            #                        ),
            #                     # LC10:
            #                     LC(name = 'Q_feld_mitte', category = 'imposed-load', file_name = 'LC10.txt',
            #                        exclusive_to = ['Q_hinten_mitte', 'Q_feld_li','s_vorne', 's_hinten', 's_feld', 's_komplett', 's_links', 's_rechts','T_schwinden','w_vonlinks_komplett','w_vonrechts_komplett','w_druck_komplett','w_sog_komplett'],
            #                        psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
            #                        ),
            #                     # LC11:
            #                     LC(name = 'Q_feld_re', category = 'imposed-load', file_name = 'LC11.txt',
            #                        exclusive_to = ['Q_hinten_mitte', 'Q_feld_li', 'Q_feld_mitte', 'Q_vorne_mitte'],
            #                        psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
            #                        ),
            #                     # LC12:
            #                     LC(name = 'Q_vorne_mitte', category = 'imposed-load', file_name = 'LC12.txt',
            #                        exclusive_to = ['Q_hinten_mitte', 'Q_feld_li', 'Q_feld_mitte', 'Q_feld_re'],
            #                        psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
            #                        ),

            #---------------------------------------------------------
            # temperature
            #---------------------------------------------------------

            #                     # LC13:
            #                     LC(name = 'T_N_neg', category = 'imposed-load', file_name = 'LC13.txt',
            #                        exclusive_to = ['T_N_pos', 'T_uo_neg', 'T_uo_pos'],
            #                        psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0
            #                        ),
            #                     # LC14:
            #                     LC(name = 'T_N_pos', category = 'imposed-load', file_name = 'LC14.txt',
            #                        exclusive_to = ['T_N_neg', 'T_uo_neg', 'T_uo_pos'],
            #                        psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0
            #                        ),
            #                     # LC15:
            #                     LC(name = 'T_uo_neg', category = 'imposed-load', file_name = 'LC15.txt',
            #                        exclusive_to = ['T_N_neg', 'T_N_pos', 'T_uo_pos'],
            #                        psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0
            #                        ),
            #                     # LC16:
            #                     LC(name = 'T_uo_pos', category = 'imposed-load', file_name = 'LC16.txt',
            #                        exclusive_to = ['T_N_neg', 'T_N_pos', 'T_uo_neg'],
            #                        psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0
            #                        ),

            #---------------------------------------------------------
            # shrinkage
            #---------------------------------------------------------
            # LC17:
            #                     LC(name='T_schwinden', category='imposed-load', file_name='LC17.txt',
            #                        psi_0=0.8, psi_1=0.7, psi_2=0.5,
            #                        ),

            #---------------------------------------------------------
            # wind load
            #---------------------------------------------------------
            # LC42:
            LC(name='w_vonlinks_komplett', category='imposed-load', file_name='LC42.txt',
               exclusive_to=[
                   'w_vonrechts_komplett', 'w_sog_komplett', 'w_druck_komplett'],
               psi_0=0.6, psi_1=0.2, psi_2=0.0
               ),
            # LC43:
            #                     LC(name='w_vonrechts_komplett', category='imposed-load', file_name='LC43.txt',
            #                        exclusive_to=['w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett'],
            #                        psi_0=0.6, psi_1=0.2, psi_2=0.0
            #                        ),
            # LC44:
            LC(name='w_sog_komplett', category='imposed-load', file_name='LC44.txt',
               exclusive_to=[
                   'w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_druck_komplett'],
               psi_0=0.6, psi_1=0.2, psi_2=0.0
               ),
            # LC45:
            LC(name='w_druck_komplett', category='imposed-load', file_name='LC45.txt',
               exclusive_to=[
                   'w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_sog_komplett'],
               psi_0=0.6, psi_1=0.2, psi_2=0.0
               ),
        ]

    #----------------------------------------------------------------------
    # (c) dimensioning with loads according to DIN EN 1991-1-3 (snow) and DIN EN 1991-1-4 (wind)
    #----------------------------------------------------------------------

    if do == 'dimensioning':

        # NOTE: 'lc_list_Q' contains all imposed loading cases (except for temperature)
        # the loading cases for temperature are evaluated separately (see 'lc_list_T' below)
        # and superposed later
        #
        lc_list_Q = [

            #------------------------------------------------------------------
            # dead load
            #-------------------------------------------------------------
            # LC1:
            LC(name='g', category='dead-load', file_name='LC1.txt'
               ),

            #-------------------------------------------------------------
            # snow
            #-------------------------------------------------------------

            # LC2:
            LC(name='s_komplett', category='imposed-load', file_name='LC2.txt',
                    exclusive_to=[
                        's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld'],
                    psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),
            #                 # LC3:
            #                 LC(name='s_verweht_re', category='imposed-load', file_name='LC3.txt',
            #                    exclusive_to=['s_komplett', 's_scheddach_re', 's_hinten', 's_feld'],
            #                    psi_0=0.5, psi_1=0.2, psi_2=0.0
            #                    ),
            # LC4:
            LC(name='s_scheddach_re', category='imposed-load', file_name='LC4.txt',
                    exclusive_to=[
                        's_komplett', 's_verweht_re', 's_hinten', 's_feld'],
                    psi_0=0.5, psi_1=0.2, psi_2=0.0
               ),
            #                 # LC5:
            #                 LC(name='s_hinten', category='imposed-load', file_name='LC5.txt',
            #                    exclusive_to=['s_komplett', 's_verweht_re', 's_scheddach_re'],
            #                    psi_0=0.5, psi_1=0.2, psi_2=0.0
            #                    ),
            #                 # LC6:
            #                 LC(name='s_feld', category='imposed-load', file_name='LC6.txt',
            #                    exclusive_to=['s_komplett', 's_verweht_re', 's_scheddach_re'],
            #                    psi_0=0.5, psi_1=0.2, psi_2=0.0
            #                    ),

            #-------------------------------------------------------------
            # wind load
            #-------------------------------------------------------------
            #                 # LC7:
            #                 LC(name='w_sog_komplett', category='imposed-load', file_name='LC7.txt',
            #                    exclusive_to=['w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_druck_komplett'],
            #                    psi_0=0.6, psi_1=0.2, psi_2=0.0
            #                    ),
            # LC8:
            LC(name='w_druck_komplett', category='imposed-load', file_name='LC8.txt',
                    exclusive_to=[
                        'w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_sog_komplett'],
                    psi_0=0.6, psi_1=0.2, psi_2=0.0
               ),
            #                 # LC9:
            #                 LC(name='w_vonlinks_komplett', category='imposed-load', file_name='LC9.txt',
            #                    exclusive_to=['w_vonrechts_komplett', 'w_sog_komplett', 'w_druck_komplett'],
            #                    psi_0=0.6, psi_1=0.2, psi_2=0.0
            #                    ),
            # LC10:
            LC(name='w_vonrechts_komplett', category='imposed-load', file_name='LC10.txt',
                    exclusive_to=[
                        'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett'],
                    psi_0=0.6, psi_1=0.2, psi_2=0.0
               ),

            #-------------------------------------------------------------
            # man load (1 KN)
            #-------------------------------------------------------------
            #                 # LC11:
            #                 LC(name='Q_hinten', category='imposed-load', file_name='LC11.txt',
            #                    exclusive_to=[ 'Q_links', 'Q_mitte', 's_komplett', 's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld', 'w_vonrechts_komplett', 'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett', 'T_schwinden'],
            #                    psi_0=0.0, psi_1=0.2, psi_2=0.0
            #                    ),
            #                 # LC12:
            #                 LC(name='Q_mitte', category='imposed-load', file_name='LC12.txt',
            #                    exclusive_to=[ 'Q_links', 'Q_hinten', 's_komplett', 's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld', 'w_vonrechts_komplett', 'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett', 'T_schwinden'],
            #                    psi_0=0.0, psi_1=0.2, psi_2=0.0
            #                    ),
            #                 # LC13:
            #                 LC(name='Q_links', category='imposed-load', file_name='LC13.txt',
            #                    exclusive_to=[ 'Q_hinten', 'Q_mitte', 's_komplett', 's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld', 'w_vonrechts_komplett', 'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett', 'T_schwinden'],
            #                    psi_0=0.0, psi_1=0.2, psi_2=0.0
            #                    ),
            #-------------------------------------------------------------
            # shrinkage
            #-------------------------------------------------------------
            # LC14:
            LC(name='T_schwinden', category='imposed-load', file_name='LC14.txt',
                    psi_0=0.8, psi_1=0.7, psi_2=0.5,
               ),

            # LC15:
            LC(name='T_N_pos', category='imposed-load', file_name='LC15.txt',
                    exclusive_to=['T_N_neg', 'T_uo_pos', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC16:
            LC(name='T_N_neg', category='imposed-load', file_name='LC16.txt',
                    exclusive_to=['T_N_pos', 'T_uo_pos', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC17:
            LC(name='T_uo_pos', category='imposed-load', file_name='LC17.txt',
                    exclusive_to=['T_N_pos', 'T_N_neg', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC18:
            LC(name='T_uo_neg', category='imposed-load', file_name='LC18.txt',
                    exclusive_to=['T_N_pos', 'T_N_neg', 'T_uo_pos'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),

        ]

        lc_list_T = [

            #------------------------------------------------------------------
            # temperature (superpose separately = corresponds to conservatively setting psi_0 = 1.0)
            #-------------------------------------------------------------

            # LC15:
            LC(name='T_N_pos', category='imposed-load', file_name='LC15.txt',
                    exclusive_to=['T_N_neg', 'T_uo_pos', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC16:
            LC(name='T_N_neg', category='imposed-load', file_name='LC16.txt',
                    exclusive_to=['T_N_pos', 'T_uo_pos', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC17:
            LC(name='T_uo_pos', category='imposed-load', file_name='LC17.txt',
                    exclusive_to=['T_N_pos', 'T_N_neg', 'T_uo_neg'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
            # LC18:
            LC(name='T_uo_neg', category='imposed-load', file_name='LC18.txt',
                    exclusive_to=['T_N_pos', 'T_N_neg', 'T_uo_pos'],
                    psi_0=1.0, psi_1=0.5, psi_2=0.0
               ),
        ]

    #--------------------------------------------------------------
    # ULS evaluation
    #--------------------------------------------------------------

    if do == 'shell-test':
        #--------------------------------------------------------
        # strength characteristics for experimental (mean) values for barrelshell
        # (specimens thickness = 2 cm; specimn width = 10 cm; 6 layers carbon)
        #--------------------------------------------------------
        if 'use_QS_values':
            # design values for quality tests of barrelshell (TTb-bs4 and
            # BT_bs4) on specimens with thickness 2 cm

            # tensile strength [kN/m]
            #
            # = 25.7 (=F_tRd)/ 0.10 ### F_Rtm = 47.6 kN (mean value)
            n_0_Rdt = n_90_Rdt = 476.

            # compressive strength [kN/m]
            #
            n_Rdc = 1360  # f_cm=68 MPa

            # bending strength [kNm/m]
            #
            # = 0.18 (=M_Rd) / 0.10 ### M_Rm = 0.33 kNm (mean value)
            m_0_Rd = m_90_Rd = 3.3

            print 'mean values (QS-bs4) used for strength characteristics  (no k_b; no reduction for scatter; no gamma_M)'

        else:
            # tensile strength [kN/m]
            #
            n_0_Rdt = n_90_Rdt = 41.1 / 0.10  # [kN/m]

            # bending strength [kNm/m]
            #
            m_0_Rd = m_90_Rd = (3.5 * 0.46 / 4.) / 0.10  # [kNm/m]

            # compressive strength [kN/m]
            # (design value; f_cd = 37,5 MPa)
            #
            n_Rdc = 750.  # = 37,5 MPa * (100 cm * 2 cm) * 0.1

        lct = LCCTableULS(data_dir=data_dir,
                          reader_type='InfoCAD',
                          #                          data_filter=remove_support_elems,
                          lc_list=lc_list,
                          strength_characteristics={'n_0_Rdt': n_0_Rdt,
                                                    'm_0_Rd': m_0_Rd,
                                                    'n_Rdc': n_Rdc,
                                                    'n_90_Rdt': n_90_Rdt,
                                                    'm_90_Rd': m_90_Rd},
                          # NO simplification used for 'k_alpha' on the
                          # resistance side
                          k_alpha_min=False,
                          show_lc_characteristic=True
                          )

        #--------------------------------------------------------------
        # 'combi_arr': array with indices of all loading case combinations
        #--------------------------------------------------------------
        #
#        print 'lct.combi_arr', lct.combi_arr.shape
#        np.savetxt('combi_arr_wo_temp_LCs', lct.combi_arr, delimiter=';')

        #--------------------------------------------------------------
        # nm-interaction plot (normal force - bending moment)
        #--------------------------------------------------------------
        #
        lct.plot_nm_interaction(save_fig_to_file='nm_interaction_shell-test')

        #--------------------------------------------------------------
        # interaction plot of material usage 'eta_nm' (utilization ratio)
        #--------------------------------------------------------------
        #
        lct.plot_eta_nm_interaction(
            save_fig_to_file='eta_nm_interaction_shell-test')

        #--------------------------------------------------------------
        # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
        # (surrounding values of all loading cases)
        #--------------------------------------------------------------
        #
        lct.plot_assess_value('eta_nm_tot',
                              save_fig_to_file='eta_nm_tot_shell-test')

        #--------------------------------------------------------------
        # brows the loading case combinations within an interactive table view
        #--------------------------------------------------------------
        lct.configure_traits()

    if do == 'predimensioning':
        #--------------------------------------------------------
        # strength characteristics (design) values for barrelshell
        # (specimens thickness = 2 cm; specimn width = 10 cm; 6 layers carbon)
        #--------------------------------------------------------

        # tensile strength [kN/m]
        #
        n_0_Rdt = n_90_Rdt = 223.  # [kN/m] # ZiE value

        # bending strength [kNm/m]
        #
        m_0_Rd = m_90_Rd = 1.7  # [kNm/m] # ZiE value

        # compressive strength [kN/m]
        # (design value; f_cd = 37,5 MPa)
        #
        n_Rdc = 750.  # = 37,5 MPa * (100 cm * 2 cm) * 0.1

        print 'design values used for strength characteristics'

        lct = LCCTableULS(data_dir=data_dir,
                          reader_type='InfoCAD',
                          data_filter=remove_support_elems,
                          lc_list=lc_list,
                          strength_characteristics={'n_0_Rdt': n_0_Rdt, 'm_0_Rd': m_0_Rd, 'n_Rdc': n_Rdc,
                                                    'n_90_Rdt': n_90_Rdt, 'm_90_Rd': m_90_Rd},
                          # NO simplification used for 'k_alpha' on the
                          # resistance side
                          k_alpha_min=True,
                          show_lc_characteristic=False
                          )

        #--------------------------------------------------------------
        # 'combi_arr': array with indices of all loading case combinations
        #--------------------------------------------------------------
        #
        print 'lct.combi_arr', lct.combi_arr.shape
#        np.savetxt('combi_arr_wo_temp_LCs', lct_Q.combi_arr, delimiter = ';')

        #--------------------------------------------------------------
        # nm-interaction plot (normal force - bending moment)
        #--------------------------------------------------------------
        #
        lct.plot_nm_interaction(save_fig_to_file='nm_interaction_LC22')

        #--------------------------------------------------------------
        # interaction plot of material usage 'eta_nm' (utilization ratio)
        #--------------------------------------------------------------
        #
        lct.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_LC22')

        #--------------------------------------------------------------
        # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
        # (surrounding values of all loading cases)
        #--------------------------------------------------------------
        #
        lct.plot_assess_value('eta_nm_tot')

        #--------------------------------------------------------------
        # brows the loading case combinations within an interactive table view
        #--------------------------------------------------------------
        lct.configure_traits()

    if do == 'dimensioning':
        #--------------------------------------------------------
        # strength characteristics (design) values for barrelshell
        # (specimens thickness = 2 cm; specimn width = 10 cm; 6 layers carbon)
        #--------------------------------------------------------
        if 'use_QS_values':
            # design values for quality tests of barrelshell (TTb-bs4 and
            # BT_bs4) on specimens with thickness 2 cm

            # tensile strength [kN/m]
            #
            # use strength reduction factor 'k_b' for discontinuities in the fabric layers due to butt joints
            # NOTE: (omit as no butt joints are located in midfield):
            # i.e. k_b = 5/6 = 0.833333
            # = 25.7 (=F_tRd)/ 0.10 ### F_Rtm = 47.6 kN (mean value)
            n_0_Rdt = n_90_Rdt = 0.83333 * 257

            # compressive strength [kN/m]
            #
            n_Rdc = 750  # f_cd=37.5 MPa

            # bending strength [kNm/m]
            #
            # = 0.18 (=M_Rd) / 0.10 ### M_Rm = 0.33 kNm (mean value)
            m_0_Rd = m_90_Rd = 0.83333 * 1.8

            print 'design values (QS-bs4) used for strength characteristics'

        else:
            # tensile strength [kN/m]
            #
            n_0_Rdt = n_90_Rdt = 223.  # [kN/m] # ZiE value

            # bending strength [kNm/m]
            #
            m_0_Rd = m_90_Rd = 1.7  # [kNm/m] # ZiE value

            # compressive strength [kN/m]
            # (design value; f_cd = 37,5 MPa)
            #
            n_Rdc = 750.  # = 37,5 MPa * (100 cm * 2 cm) * 0.1

            print 'design values (ZiE) used for strength characteristics'

        # LCCTable for imposed loads (without temperature)
        #
        lct_Q = LCCTableULS(data_dir=data_dir,
                            reader_type='InfoCAD',
                            data_filter=remove_support_elems,
                            lc_list=lc_list_Q,
                            strength_characteristics={'n_0_Rdt': n_0_Rdt, 'm_0_Rd': m_0_Rd, 'n_Rdc': n_Rdc,
                                                      'n_90_Rdt': n_90_Rdt, 'm_90_Rd': m_90_Rd},
                            # simplification: use the minimum value for k_alpha
                            # on the resistance side
                            k_alpha_min=False,
                            #                          show_lc_characteristic = True
                            )

        # LCCTable for temperature loading cases only
        #
        lct_T = LCCTableULS(data_dir=data_dir,
                            reader_type='InfoCAD',
                            data_filter=remove_support_elems,
                            lc_list=lc_list_T,
                            strength_characteristics={'n_0_Rdt': n_0_Rdt, 'm_0_Rd': m_0_Rd, 'n_Rdc': n_Rdc,
                                                      'n_90_Rdt': n_90_Rdt, 'm_90_Rd': m_90_Rd},
                            # k_alpha_min=True,  # simplification: use the
                            # minimum value for k_alpha on the resistance side
                            # NOTE: no simplification: use calculated values
                            # for k_alpha on the resistance side
                            k_alpha_min=False,
                            #       and evaluation in 'alpha_1' AND 'alpha_2' direction
                            #                          show_lc_characteristic = True
                            )

        #--------------------------------------------------------------
        # 'combi_arr': array with indices of all loading case combinations
        #--------------------------------------------------------------
        #
        print 'lct_Q.combi_arr', lct_Q.combi_arr.shape, '\n'
#        np.savetxt('combi_arr_wo_temp_LCs', lct_Q.combi_arr, delimiter=';')

        #--------------------------------------------------------------
        # nm-interaction plot (normal force - bending moment)
        #--------------------------------------------------------------
        #
#        lct_T.plot_nm_interaction(save_max_min_nm_to_file='max_min_nm_arr_LC15-18', save_fig_to_file='nm_interaction_LC15-18')
#        lct_Q.plot_nm_interaction(save_fig_to_file='nm_interaction_LC1-14')
#        lct_Q.plot_nm_interaction(add_max_min_nm_from_file='max_min_nm_arr_LC15-18', save_fig_to_file='nm_interaction_LC1-18')

        #--------------------------------------------------------------
        # interaction plot of material usage 'eta_nm' (utilization ratio)
        #--------------------------------------------------------------
        #
#        lct_T.plot_eta_nm_interaction(save_max_min_eta_nm_to_file='max_min_eta_nm_arr_LC15-18', save_fig_to_file='eta_nm_interaction_LC15-18')
#        lct_Q.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_LC1-14')
#        lct_Q.plot_eta_nm_interaction(add_max_min_eta_nm_from_file='max_min_eta_nm_arr_LC15-18', save_fig_to_file='eta_nm_interaction_LC1-18')

        #--------------------------------------------------------------
        # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
        # (surrounding values of all loading cases)
        #--------------------------------------------------------------
        #
#        lct_T.plot_assess_value('eta_nm_tot', save_assess_values_to_file='eta_nm_tot_LC15-18')
#        lct_Q.plot_assess_value('eta_nm_tot', add_assess_values_from_file='eta_nm_tot_LC15-18')

        #--------------------------------------------------------------
        # brows the loading case combinations within an interactive table view
        #--------------------------------------------------------------
        lct_Q.configure_traits()
#        lct_T.configure_traits()
