if __name__ == '__main__':

    from etsproxy.mayavi import \
        mlab

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
        '''
        elem_no = lcc_table.geo_data_orig['elem_no']
        # support elements for fine regular mesh (5cmx5cm-elements)
        # as stored in '2cm-feines-Netz'

        # first row elements (direct contact to support node)
        #
        support_elems_list = [692, 693, 2324, 2325, 744, 745, 2376, 2377, \
        # second row elements (= adjacent elements to first row elements)
        #
                              691, 694, 2323, 2326, 743, 746, 2375, 2378 ]
        cond_arr = np.array([elem_no != support_elem_number for support_elem_number in support_elems_list])
        cond_elem_active = np.product(cond_arr, axis=0)

#        X = lcc_table.geo_data_orig['X']
#        Y = lcc_table.geo_data_orig['Y']
#        elem_size = Y[1] - Y[2]
#        print 'elem_size', elem_size
#        cond_X1 = X < -1.03 + 2. * elem_size
#        cond_X2 = X > 1.03 - 2. * elem_size
#        cond_Y1a = Y > -0.40 - 2. * elem_size
#        cond_Y1b = Y < -0.40 + 2. * elem_size
#        cond_Y2a = Y > -3.00 - 2. * elem_size
#        cond_Y2b = Y < -3.00 + 2. * elem_size
#        cond_X = cond_X1 + cond_X2
#        cond_Y = cond_Y1a * cond_Y1b + cond_Y2a * cond_Y2b
#        cond_elem_active = np.where( 1 - cond_X * cond_Y )[0]

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

#    do = 'predimensioning'
    do = 'dimensioning'

    # Vorbemessung
    #
    if do == 'predimensioning':
        data_dir = os.path.join(simdb.simdb_dir,
                                'simdata',
                                'input_data_barrelshell',
                                '2cm-feines-Netz',
                                )

    # Bemessung
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
    #
    #---------------------------------------------------------
    # "staendige und voruebergehende Bemessungssitauation":
    #---------------------------------------------------------
    if do == 'predimensioning':

        lc_list = [

                     #----------------------------------------------------------------------
                     # dead load
                     #----------------------------------------------------------------------
                     # LC1:
                     LC(name='g', category='dead-load', file_name='LC1.txt'
                        ),

                     #----------------------------------------------------------------------
                     # snow
                     #----------------------------------------------------------------------

                     # LC2:
                     LC(name='s_hinten', category='imposed-load', file_name='LC2.txt',
                        exclusive_to=['s_feld', 's_vorne', 's_links', 's_rechts', 's_komplett'],
                        psi_0=0.5, psi_1=0.2, psi_2=0.0
                        ),
                     # LC3:
                     LC(name='s_feld', category='imposed-load', file_name='LC3.txt',
                        exclusive_to=['s_hinten', 's_vorne', 's_links', 's_rechts', 's_komplett'],
                        psi_0=0.5, psi_1=0.2, psi_2=0.0
                        ),
                     # LC4:
#                     LC(name = 's_vorne', category = 'imposed-load', file_name = 'LC4.txt',
#                        exclusive_to = ['s_hinten', 's_feld', 's_links', 's_rechts', 's_komplett'],
#                        psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0
#                        ),
#                     # LC5:
                     LC(name='s_links', category='imposed-load', file_name='LC5.txt',
                        exclusive_to=['s_hinten', 's_feld', 's_vorne', 's_rechts', 's_komplett'],
                        psi_0=0.5, psi_1=0.2, psi_2=0.0
                        ),
                     # LC6:
                     LC(name='s_rechts', category='imposed-load', file_name='LC6.txt',
                        exclusive_to=['s_hinten', 's_feld', 's_vorne', 's_links', 's_komplett'],
                        psi_0=0.5, psi_1=0.2, psi_2=0.0
                        ),
                     # LC7:
                     LC(name='s_komplett', category='imposed-load', file_name='LC7.txt',
                        exclusive_to=['s_hinten', 's_feld', 's_vorne', 's_links', 's_rechts'],
                        psi_0=0.5, psi_1=0.2, psi_2=0.0
                        ),

                     #----------------------------------------------------------------------
                     # man load (1 KN)
                     #----------------------------------------------------------------------
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
    #                 # LC11:
    #                 LC(name = 'Q_feld_re', category = 'imposed-load', file_name = 'LC11.txt',
    #                    exclusive_to = ['Q_hinten_mitte', 'Q_feld_li', 'Q_feld_mitte', 'Q_vorne_mitte'],
    #                    psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
    #                    ),
    #                 # LC12:
    #                 LC(name = 'Q_vorne_mitte', category = 'imposed-load', file_name = 'LC12.txt',
    #                    exclusive_to = ['Q_hinten_mitte', 'Q_feld_li', 'Q_feld_mitte', 'Q_feld_re'],
    #                    psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
    #                    ),

                     #----------------------------------------------------------------------
                     # temperature
                     #----------------------------------------------------------------------

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

                     #----------------------------------------------------------------------
                     # shrinkage
                     #----------------------------------------------------------------------
    #                 # LC17:
                     LC(name='T_schwinden', category='imposed-load', file_name='LC17.txt',
                        psi_0=0.8, psi_1=0.7, psi_2=0.5,
                        ),

                     #----------------------------------------------------------------------
                     # wind load
                     #----------------------------------------------------------------------
                     # LC42:
                     LC(name='w_vonlinks_komplett', category='imposed-load', file_name='LC42.txt',
                        exclusive_to=['w_vonrechts_komplett', 'w_sog_komplett', 'w_druck_komplett'],
                        psi_0=0.6, psi_1=0.2, psi_2=0.0
                        ),
#                     # LC43:
                     LC(name='w_vonrechts_komplett', category='imposed-load', file_name='LC43.txt',
                        exclusive_to=['w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett'],
                        psi_0=0.6, psi_1=0.2, psi_2=0.0
                        ),
                     # LC44:
                     LC(name='w_sog_komplett', category='imposed-load', file_name='LC44.txt',
                        exclusive_to=['w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_druck_komplett'],
                        psi_0=0.6, psi_1=0.2, psi_2=0.0
                        ),
                     # LC45:
                     LC(name='w_druck_komplett', category='imposed-load', file_name='LC45.txt',
                        exclusive_to=['w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_sog_komplett'],
                        psi_0=0.6, psi_1=0.2, psi_2=0.0
                        ),
                   ]


    #----------------------------------------------------------------------
    # loads according to DIN EN 1991-1-3 (snow) and DIN EN 1991-1-4 (wind)
    #----------------------------------------------------------------------

    if do == 'dimensioning':

        # NOTE: 'lc_list_Q' contains all imposed loading cases (except for temperature)
        # the loading cases for temperature are evaluated separately (see 'lc_list_T' below)
        # and superposed later
        #
        lc_list_Q = [

                 #----------------------------------------------------------------------
                 # dead load
                 #----------------------------------------------------------------------
                 # LC1:
                 LC(name='g', category='dead-load', file_name='LC1.txt'
                    ),

                 #----------------------------------------------------------------------
                 # snow
                 #----------------------------------------------------------------------

                 # LC2:
                 LC(name='s_komplett', category='imposed-load', file_name='LC2.txt',
                    exclusive_to=['s_verweht_re', 's_scheddach_re', 's_hinten', 's_feld'],
                    psi_0=0.5, psi_1=0.2, psi_2=0.0
                    ),
                 # LC3:
                 LC(name='s_verweht_re', category='imposed-load', file_name='LC3.txt',
                    exclusive_to=['s_komplett', 's_scheddach_re', 's_hinten', 's_feld'],
                    psi_0=0.5, psi_1=0.2, psi_2=0.0
                    ),
                 # LC4:
                 LC(name='s_scheddach_re', category='imposed-load', file_name='LC4.txt',
                    exclusive_to=['s_komplett', 's_verweht_re', 's_hinten', 's_feld'],
                    psi_0=0.5, psi_1=0.2, psi_2=0.0
                    ),
                 # LC5:
                 LC(name='s_hinten', category='imposed-load', file_name='LC5.txt',
                    exclusive_to=['s_komplett', 's_verweht_re', 's_scheddach_re'],
                    psi_0=0.5, psi_1=0.2, psi_2=0.0
                    ),
                 # LC6:
                 LC(name='s_feld', category='imposed-load', file_name='LC6.txt',
                    exclusive_to=['s_komplett', 's_verweht_re', 's_scheddach_re'],
                    psi_0=0.5, psi_1=0.2, psi_2=0.0
                    ),

                 #----------------------------------------------------------------------
                 # wind load
                 #----------------------------------------------------------------------
                 # LC7:
                 LC(name='w_sog_komplett', category='imposed-load', file_name='LC7.txt',
                    exclusive_to=['w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_druck_komplett'],
                    psi_0=0.6, psi_1=0.2, psi_2=0.0
                    ),
                 # LC8:
                 LC(name='w_druck_komplett', category='imposed-load', file_name='LC8.txt',
                    exclusive_to=['w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_sog_komplett'],
                    psi_0=0.6, psi_1=0.2, psi_2=0.0
                    ),
                 # LC9:
                 LC(name='w_vonlinks_komplett', category='imposed-load', file_name='LC9.txt',
                    exclusive_to=['w_vonrechts_komplett', 'w_sog_komplett', 'w_druck_komplett'],
                    psi_0=0.6, psi_1=0.2, psi_2=0.0
                    ),
                 # LC10:
                 LC(name='w_vonrechts_komplett', category='imposed-load', file_name='LC10.txt',
                    exclusive_to=['w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett'],
                    psi_0=0.6, psi_1=0.2, psi_2=0.0
                    ),

                 #----------------------------------------------------------------------
                 # man load (1 KN)
                 #----------------------------------------------------------------------
                 # LC11:
                 LC(name='Q_hinten', category='imposed-load', file_name='LC11.txt',
                    exclusive_to=[ 'Q_links', 'Q_mitte', 's_komplett', 's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld', 'w_vonrechts_komplett', 'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett', 'T_schwinden'],
                    psi_0=0.0, psi_1=0.2, psi_2=0.0
                    ),
                 # LC12:
                 LC(name='Q_mitte', category='imposed-load', file_name='LC12.txt',
                    exclusive_to=[ 'Q_links', 'Q_hinten', 's_komplett', 's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld', 'w_vonrechts_komplett', 'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett', 'T_schwinden'],
                    psi_0=0.0, psi_1=0.2, psi_2=0.0
                    ),
                 # LC13:
                 LC(name='Q_links', category='imposed-load', file_name='LC13.txt',
                    exclusive_to=[ 'Q_hinten', 'Q_mitte', 's_komplett', 's_verweht_re', 's_scheddach_re', 's_hinten', 's_feld', 'w_vonrechts_komplett', 'w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett', 'T_schwinden'],
                    psi_0=0.0, psi_1=0.2, psi_2=0.0
                    ),

                 #----------------------------------------------------------------------
                 # shrinkage
                 #----------------------------------------------------------------------
                 # LC14:
                 LC(name='T_schwinden', category='imposed-load', file_name='LC14.txt',
                    psi_0=0.8, psi_1=0.7, psi_2=0.5,
                    ),

               ]


        lc_list_T = [

                 #----------------------------------------------------------------------
                 # temperature (superpose separately = corresponds to conservatively setting psi_0 = 1.0)
                 #----------------------------------------------------------------------

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


    if do == 'dimensioning':

        # LCCTable for imposed loads (without temperature)
        #
        lct_Q = LCCTableULS(data_dir=data_dir,
                          reader_type='InfoCAD',
                          data_filter=remove_support_elems,
                          lc_list=lc_list_Q,
                          k_alpha_min=True,  # simplification: use the minimum value for k_alpha on the resistance side
#                          show_lc_characteristic = True
                          )

        # LCCTable for temperature loading cases only
        #
        lct_T = LCCTableULS(data_dir=data_dir,
                          reader_type='InfoCAD',
                          data_filter=remove_support_elems,
                          lc_list=lc_list_T,
                          k_alpha_min=True,  # simplification: use the minimum value for k_alpha on the resistance side
#                          show_lc_characteristic = True
                          )

        #--------------------------------------------------------------
        # 'combi_arr': array with indices of all loading case combinations
        #--------------------------------------------------------------
        #
#        print 'lct_Q.combi_arr', lct_Q.combi_arr.shape
#        np.savetxt('combi_arr_wo_temp_LCs', lct_Q.combi_arr, delimiter = ';')

        #--------------------------------------------------------------
        # nm-interaction plot (normal force - bending moment)
        #--------------------------------------------------------------
        #
#        lct_T.plot_nm_interaction( save_max_min_nm_to_file = 'max_min_nm_arr_LC15-18', save_fig_to_file = 'nm_interaction_LC15-18')
#        lct_Q.plot_nm_interaction( save_fig_to_file = 'nm_interaction_LC1-14' )
#        lct_Q.plot_nm_interaction( add_max_min_nm_from_file = 'max_min_nm_arr_LC15-18', save_fig_to_file = 'nm_interaction_LC1-18' )

        #--------------------------------------------------------------
        # interaction plot of material usage 'eta_nm' (Ausnutzungsgrad)
        #--------------------------------------------------------------
        #
#        lct_T.plot_eta_nm_interaction( save_max_min_eta_nm_to_file = 'max_min_eta_nm_arr_LC15-18', save_fig_to_file = 'eta_nm_interaction_LC15-18' )
#        lct_Q.plot_eta_nm_interaction( save_fig_to_file = 'eta_nm_interaction_LC1-14')
#        lct_Q.plot_eta_nm_interaction( add_max_min_eta_nm_from_file = 'max_min_eta_nm_arr_LC15-18', save_fig_to_file = 'eta_nm_interaction_LC1-18')

        #--------------------------------------------------------------
        # plot of structure with color indication of material usage 'eta_nm' (Ausnutzungsgrad)
        # (surrounding values of all loading cases)
        #--------------------------------------------------------------
        #
#        lct_T.plot_assess_value( save_assess_values_to_file = 'eta_nm_tot_LC15-18' )
#        lct_Q.plot_assess_value( add_assess_values_from_file = 'eta_nm_tot_LC15-18' )

        #--------------------------------------------------------------
        # brows the loading case combinations within an interactive table view
        #--------------------------------------------------------------
#        lct_Q.configure_traits()
        lct_T.configure_traits()


