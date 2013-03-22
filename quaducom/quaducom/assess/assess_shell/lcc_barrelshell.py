


if __name__ == '__main__':

    from etsproxy.mayavi import \
        mlab

    from matresdev.db.simdb import \
        SimDB

    import os

    from lcc_table import LCCTableULS, LC, LCCTableSLS

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    #---------------------------------------------
    # 2 shells: 
    # new geometry with new loading cases plus waterfilling
    #---------------------------------------------

    data_dir = os.path.join(simdb.simdb_dir,
                            'simdata', 
                            'input_data_barrelshell',
#                            '2cm-feines-Netz', 
                            '2cm-Elastomer', 
#                            '2cm', 
#                            '2cm-Zugausfall-komplett', 
#                            '2cm-Zugausfall-vorne', 
#                            '2cm-Zugausfall-hinten', 
#                            '2cm-Zugausfall-links', 
#                            '2cm-Zugausfall-rechts'
#                            '2cm-Zugausfall-hinten-links'
#                            '2cm-Zugausfall-vorne-links'
#                            '2cm-Zugausfall-vorne-rechts'
#                            '2cm-Zugausfall-hinten-rechts'
#                            '2cm-Zugausfall-nie-hinten-links'
#                            '2cm-Zugausfall-nie-vorne-links'
#                            '2cm-Zugausfall-nie-vorne-rechts'
#                            '2cm-Zugausfall-nie-hinten-rechts'
#                            '2cm-Zugausfall-hinten-links-vorne-rechts'
#                            '2cm-Zugausfall-vorne-links-hinten-rechts'
#                            '2cm-fest', 
#                            '3cm', 
#                            '3-4cm', 
                            )

    #------------------------
    # define loading cases:
    #------------------------
    # NOTE:
    #
    #---------------------------------------------------------
    # "staendige und voruebergehende Bemessungssitauation":
    #---------------------------------------------------------
    lc_list = [

                 #----------------------------------------------------------------------
                 # dead load
                 #----------------------------------------------------------------------
#                 # LC1:
                 LC(name = 'g', category = 'dead-load', file_name = 'LC1.txt'
                    ),

                 #----------------------------------------------------------------------
                 # snow
                 #----------------------------------------------------------------------

                 # LC2:
                 LC(name = 's_hinten', category = 'imposed-load', file_name = 'LC2.txt',
                    exclusive_to = ['s_feld', 's_vorne', 's_links', 's_rechts', 's_komplett'],
                    psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0
                    ),
                 # LC3:
                 LC(name = 's_feld', category = 'imposed-load', file_name = 'LC3.txt',
                    exclusive_to = ['s_hinten', 's_vorne', 's_links', 's_rechts', 's_komplett'],
                    psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0
                    ),
                 # LC4:
                 LC(name = 's_vorne', category = 'imposed-load', file_name = 'LC4.txt',
                    exclusive_to = ['s_hinten', 's_feld', 's_links', 's_rechts', 's_komplett'],
                    psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0
                    ),
                 # LC5:
                 LC(name = 's_links', category = 'imposed-load', file_name = 'LC5.txt',
                    exclusive_to = ['s_hinten', 's_feld', 's_vorne', 's_rechts', 's_komplett'],
                    psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0
                    ),
                 # LC6:
                 LC(name = 's_rechts', category = 'imposed-load', file_name = 'LC6.txt',
                    exclusive_to = ['s_hinten', 's_feld', 's_vorne', 's_links', 's_komplett'],
                    psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0
                    ),
                 # LC7:
                 LC(name = 's_komplett', category = 'imposed-load', file_name = 'LC7.txt',
                    exclusive_to = ['s_hinten', 's_feld', 's_vorne', 's_links', 's_rechts'],
                    psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0
                    ),

                 #----------------------------------------------------------------------
                 # man load (1 KN)
                 #----------------------------------------------------------------------
#                 # LC8:
#                 LC(name = 'Q_hinten_mitte', category = 'imposed-load', file_name = 'LC8.txt',
#                    exclusive_to = ['Q_feld_li', 'Q_feld_mitte', 'Q_feld_re', 'Q_vorne_mitte'],
#                    psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
#                    ),
#                 # LC9:
#                 LC(name = 'Q_feld_li', category = 'imposed-load', file_name = 'LC9.txt',
#                    exclusive_to = ['Q_hinten_mitte', 'Q_feld_li', 'Q_feld_mitte', 'Q_vorne_mitte'],
#                    psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
#                    ),
#                 # LC10:
#                 LC(name = 'Q_feld_mitte', category = 'imposed-load', file_name = 'LC10.txt',
#                    exclusive_to = ['Q_hinten_mitte', 'Q_feld_li', 'Q_feld_mitte', 'Q_vorne_mitte'],
#                    psi_0 = 0.0, psi_1 = 0.2, psi_2 = 0.0
#                    ),
                 # LC11:
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

                 # LC13:
#                 LC(name = 'T_N_neg', category = 'imposed-load', file_name = 'LC13.txt',
#                    exclusive_to = ['T_N_pos', 'T_uo_neg', 'T_uo_pos'],
#                    psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0                 
#                    ),
                 # LC14:
#                 LC(name = 'T_N_pos', category = 'imposed-load', file_name = 'LC14.txt',
#                    exclusive_to = ['T_N_neg', 'T_uo_neg', 'T_uo_pos'],
#                    psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0              
#                    ),
                 # LC15:
                 LC(name = 'T_uo_neg', category = 'imposed-load', file_name = 'LC15.txt',
                    exclusive_to = ['T_N_neg', 'T_N_pos', 'T_uo_pos'],
                    psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0                   
                    ),
                 # LC16:
                 LC(name = 'T_uo_pos', category = 'imposed-load', file_name = 'LC16.txt',
                    exclusive_to = ['T_N_neg', 'T_N_pos', 'T_uo_neg'],
                    psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0                    
                    ),
                    
                 #----------------------------------------------------------------------
                 # shrinkage 
                 #----------------------------------------------------------------------
                 # LC17:
                 LC(name = 'T_schwinden', category = 'imposed-load', file_name = 'LC17.txt',
                    psi_0 = 0.8, psi_1 = 0.7, psi_2 = 0.5,
                    ),

                 #----------------------------------------------------------------------
                 # cable load (experimental test setup) 
                 #----------------------------------------------------------------------
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
#                 LC(name = 'cable-FL-m3', category = 'imposed-load', file_name = 'LC22.txt',
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
                 LC(name = 'w_vonlinks_komplett', category = 'imposed-load', file_name = 'LC42.txt',
                    exclusive_to = ['w_vonrechts_komplett', 'w_sog_komplett', 'w_druck_komplett'],
                    psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0
                    ),
                 # LC43:
                 LC(name = 'w_vonrechts_komplett', category = 'imposed-load', file_name = 'LC43.txt',
                    exclusive_to = ['w_vonlinks_komplett', 'w_sog_komplett', 'w_druck_komplett'],
                    psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0
                    ),
                 # LC44:
                 LC(name = 'w_sog_komplett', category = 'imposed-load', file_name = 'LC44.txt',
                    exclusive_to = ['w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_druck_komplett'],
                    psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0
                    ),
                 # LC45:
                 LC(name = 'w_druck_komplett', category = 'imposed-load', file_name = 'LC45.txt',
                    exclusive_to = ['w_vonlinks_komplett', 'w_vonrechts_komplett', 'w_sog_komplett'],
                    psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0
                    ),                                        
               ]

#    #---------------------------------------------------------
#    # "aussergewoehnliche Bemessungssituation":
#    #---------------------------------------------------------

#--------------------------------------------------------------

    do = 'ULS'

    if do == 'ULS':
        lct = LCCTableULS(data_dir = data_dir,
                          reader_type = 'InfoCAD',
                          lc_list = lc_list,
                          show_lc_characteristic = True
                          )

        lct.configure_traits()
#        lct.plot_n_tex()


    if do == 'SLS':
        lct = LCCTableSLS( data_dir = data_dir,
                           reader_type = 'InfoCAD',
                           lc_list = lc_list,
    #                       show_lc_characteristic = True,
                           show_lc_characteristic = False,

                           # NOTE: chose here the load-case combination for servicability limit state
                           # admissible stresses are defined for bending ('f_m') and tension ('f_ck') 
                           # in file 'ls_table' (set to approx. 4-5 MPa as observed in the tensile test 
                           # (composite first cracking strength = linear branch of the stress strain diagram) 
                           #
    #                       combination_SLS = 'rare',
                           combination_SLS = 'freq',
    #                       combination_SLS = 'perm',
                            )
        lct.configure_traits()




#--------------------------------------------------------------------------
#    mlab.figure(figure = "barrelshell",
#                 bgcolor = (1.0, 1.0, 1.0),
#                 fgcolor = (0.0, 0.0, 0.0))
#    lct.plot_geo(mlab)
#    lct.plot_sr(mlab)
#    mlab.show()

#    print 'lc_arr', lct.lc_arr
#    print 'lc_list[0].sr_arr.shape[0]', lct.lc_list[0].sr_arr.shape[0]
#    print 'lc_arr.shape', lct.lc_arr.shape
#    print 'combi_arr', lct.combi_arr
#    print 'combi_arr.shape', lct.combi_arr.shape
#    print 'lcc_arr', lct.lcc_arr
