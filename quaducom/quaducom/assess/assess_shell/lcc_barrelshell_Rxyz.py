if __name__ == '__main__':

    from mayavi import \
        mlab

    from matresdev.db.simdb import \
        SimDB

    import numpy as np

    import os
    import string

    # changes concerning 'Rxyz': overload the default definitions of
    # 'LC' without changing the default class names
    # as those are used in the definition of LCCTable
    from .lcc_table_Rxyz import LCCTableULS, LC

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    #------------------------
    # define directory:
    #------------------------

#    do = 'predimensioning'
    do = 'dimensioning'

    if do == 'predimensioning':
        data_dir = os.path.join(simdb.simdb_dir,
                                'simdata',
                                'input_data_barrelshell',
                                '2cm-feines-Netz',
                                )

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
                     LC(name='s_vorne', category='imposed-load', file_name='LC4.txt',
                        exclusive_to=['s_hinten', 's_feld', 's_links', 's_rechts', 's_komplett'],
                        psi_0=0.5, psi_1=0.2, psi_2=0.0
                        ),
                     # LC5:
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
                     # LC17:
                     LC(name='T_schwinden', category='imposed-load', file_name='LC17.txt',
                        psi_0=0.8, psi_1=0.7, psi_2=0.5,
                        ),

                     #----------------------------------------------------------------------
                     # cable load (experimental test setup)
                     #----------------------------------------------------------------------
    #                 LC(name = 'cable-FL-m3', category = 'imposed-load', file_name = 'LC22.txt',
    #                    # loading case only used for evaluation of experimental test setup
    #                    # (estimation of load-bearing capacity)
    #                    gamma_unf = 1000.0, gamma_fav = 1.0,
    #                    psi_0 = 1.0,
    #                    ),


    #                 # LC18:
    #                 LC(name = 'cable-L-m0', category = 'imposed-load', file_name = 'LC18.txt',
    #                    # loading case only used for evaluation of experimental test setup
    #                    # (estimation of load-bearing capacity)
    #                    gamma_unf = 1.0, gamma_fav = 1.0,
    #                    psi_0 = 1.0,
    #                    ),

    #                 LC(name = 'cable-FL-m4', category = 'imposed-load', file_name = 'LC19.txt',
    #                    # loading case only used for evaluation of experimental test setup
    #                    # (estimation of load-bearing capacity)
    #                    gamma_unf = 1.0, gamma_fav = 1.0,
    #                    psi_0 = 1.0,
    #                    ),
    #                 LC(name = 'cable-FL-m2', category = 'imposed-load', file_name = 'LC21.txt',
    #                    # loading case only used for evaluation of experimental test setup
    #                    # (estimation of load-bearing capacity)
    #                    gamma_unf = 1.0, gamma_fav = 1.0,
    #                    psi_0 = 1.0,
    #                    ),
    #                 LC(name = 'cable-FL-m0', category = 'imposed-load', file_name = 'LC20.txt',
    #                    # loading case only used for evaluation of experimental test setup
    #                    # (estimation of load-bearing capacity)
    #                    gamma_unf = 1.0, gamma_fav = 1.0,
    #                    psi_0 = 1.0,
    #                    ),

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

    do = 'ULS'

    if do == 'ULS':

#        assess_name = 'max_Rx' # @todo: compare with shear resistance of the screw
#        assess_name = 'min_Rx'
#        assess_name = 'max_Ry'
#        assess_name = 'min_Ry'
#        assess_name = 'max_Rz'
#        assess_name = 'min_Rz' # @todo: compare with pull-out resistance of the screw
#        assess_name = 'max_Rres'
        assess_name = 'max_eta_R_tot'

        # LCCTable for imposed loads (without temperature)
        #
        lct_Q = LCCTableULS(data_dir=data_dir,
                            reader_type='InfoCADRxyz',
                            lc_list=lc_list_Q,
                            assess_name=assess_name,
#                            show_lc_characteristic = True
                           )

        # LCCTable for temperature loading cases only
        #
        lct_T = LCCTableULS(data_dir=data_dir,
                            reader_type='InfoCADRxyz',
                            lc_list=lc_list_T,
                            assess_name=assess_name,
#                            show_lc_characteristic = True
                           )

        #--------------------------------------------------------------
        # 'combi_arr': array with indices of all loading case combinations
        #--------------------------------------------------------------
        #
#        print 'lct_Q.combi_arr', lct_Q.combi_arr.shape
#        np.savetxt('combi_arr_wo_temp_LCs', lct_Q.combi_arr, delimiter = ';')

        #--------------------------------------------------------------
        # brows the loading case combinations within an interactive table view
        #--------------------------------------------------------------
#        lct_Q.configure_traits()
#        lct_T.configure_traits()

        #--------------------------------------------------------------
        # RxRz-interaction plot (tangential reaction force - radial reaction force)
        #--------------------------------------------------------------
        #
#        lct_T.plot_RxRz_interaction( save_max_min_RxRz_to_file = 'max_min_RxRz_arr_LC15-18', save_fig_to_file = 'RxRz_interaction_LC15-18')
#        lct_Q.plot_RxRz_interaction( save_fig_to_file = 'RxRz_interaction_LC1-14' )
        lct_Q.plot_RxRz_interaction(add_max_min_RxRz_from_file='max_min_RxRz_arr_LC15-18', save_fig_to_file='RxRz_interaction_LC1-18')
#        lct_Q.plot_RxRz_interaction( show_tension_only = True )

        #--------------------------------------------------------------
        # interaction plot of build-in screw usage 'eta_RxRz' (Ausnutzungsgrad)
        #--------------------------------------------------------------
        #
#        lct_T.plot_eta_RxRz_interaction(show_tension_only=True, save_fig_to_file='eta_RxRz_interaction_LC15-18')
#        lct_Q.plot_eta_RxRz_interaction( show_tension_only = True, save_fig_to_file = 'eta_RxRz_interaction_LC1-14')
        lct_Q.plot_eta_RxRz_interaction(show_tension_only=True, add_max_min_RxRz_from_file='max_min_RxRz_arr_LC15-18', save_fig_to_file='eta_RxRz_interaction_LC1-18')
