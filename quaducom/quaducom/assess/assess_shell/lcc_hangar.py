if __name__ == '__main__':

    from matresdev.db.simdb import \
        SimDB

    import numpy as np

    import os

    from lcc_table import LCCTableULS, LC, LCCTableSLS

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    # dimensioning of barrel shells for all loading cases
    # (incl. superposition of temperature loading cases)
    #
    data_dir = os.path.join(simdb.simdb_dir,
                            'simdata',
                            'input_data_barrelshell',
                            '2cm-feines-Netz-EC1',
                            )

    #------------------------
    # define loading cases:
    #------------------------

    #----------------------------------------------------------------------
    # (c) dimensioning with loads according to DIN EN 1991-1-3 (snow) and DIN EN 1991-1-4 (wind)
    #----------------------------------------------------------------------

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
    ]

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

    print 'design values (ZiE) used for strength characteristics'

    # LCCTable for imposed loads (without temperature)
    #
    lct_Q = LCCTableULS(data_dir=data_dir,
                        reader_type='InfoCAD',
                        lc_list=lc_list_Q,
                        strength_characteristics={'n_0_Rdt': n_0_Rdt, 'm_0_Rd': m_0_Rd, 'n_Rdc': n_Rdc,
                                                  'n_90_Rdt': n_90_Rdt, 'm_90_Rd': m_90_Rd},
                        # simplification: use the minimum value for k_alpha
                        # on the resistance side
                        k_alpha_min=False,
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
#        lct_Q.plot_nm_interaction(save_fig_to_file='nm_interaction_LC1-14')
#        lct_Q.plot_nm_interaction(add_max_min_nm_from_file='max_min_nm_arr_LC15-18', save_fig_to_file='nm_interaction_LC1-18')

    #--------------------------------------------------------------
    # interaction plot of material usage 'eta_nm' (utilization ratio)
    #--------------------------------------------------------------
    #
#        lct_Q.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_LC1-14')
#        lct_Q.plot_eta_nm_interaction(add_max_min_eta_nm_from_file='max_min_eta_nm_arr_LC15-18', save_fig_to_file='eta_nm_interaction_LC1-18')

    #--------------------------------------------------------------
    # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
    # (surrounding values of all loading cases)
    #--------------------------------------------------------------
    #
#        lct_Q.plot_assess_value('eta_nm_tot', add_assess_values_from_file='eta_nm_tot_LC15-18')

    #--------------------------------------------------------------
    # brows the loading case combinations within an interactive table view
    #--------------------------------------------------------------
    lct_Q.configure_traits()
