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
                            'input_data_hangar',
                            'full',
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
        LC(name='g_excl.P', category='dead-load', file_name='LC1.txt'
           ),
        #LC10:
        LC(name='P', category='dead-load', file_name='LC10.txt', gamma_unf = 1.00,
           ),


        #-------------------------------------------------------------
        # snow
        # load cases: 's_sym.', 's_unsym.', 's_unsym_dreieck',
        #-------------------------------------------------------------
        # LC2:
       # LC(name='s_sym.', category='imposed-load', file_name='LC2.txt',
       #        exclusive_to=[
        #            's_unsym_dreieck', 's_unsym.'],
        #        psi_0=0.5, psi_1=0.2, psi_2=0.0
       #    ),
                 # LC3:
       # LC(name='s_unsym.', category='imposed-load', file_name='LC3.txt',
       #        exclusive_to=[
       #             's_sym.', 's_unsym_dreieck',],
        #        psi_0=0.5, psi_1=0.2, psi_2=0.0
        #   ),
                # LC4:
      #  LC(name='s_unsym_dreieck.', category='imposed-load', file_name='LC4.txt',
       #        exclusive_to=[
       #             's_unsym.', 's_sym.'],
        #        psi_0=0.5, psi_1=0.2, psi_2=0.0
         #  ), 
            #---------------------------------------------------------
            # wind load
            #load cases: 'w_X', 'w_Y+', w_Y-',
            #---------------------------------------------------------
            # LC5:
   #         LC(name='w_X', category='imposed-load', file_name='LC5.txt',
    #           exclusive_to=[
     #              'w_Y+', 'w_Y-'],
      #         psi_0=0.6, psi_1=0.2, psi_2=0.0
       #        ), 
            # LC6:
        #    LC(name='w_Y+', category='imposed-load', file_name='LC6.txt',
         #      exclusive_to=[
          #         'w_X', 'w_Y-'],
           #    psi_0=0.6, psi_1=0.2, psi_2=0.0
            #   ), 
            # LC7:
#            LC(name='w_Y-', category='imposed-load', file_name='LC7.txt',
 #              exclusive_to=[
  #                 'w_X', 'w_Y+'],
   #            psi_0=0.6, psi_1=0.2, psi_2=0.0
    #           ),  
            #---------------------------------------------------------
            # temperature
            #---------------------------------------------------------
            # LC11:
           # LC(name = 'T_N_neg', category = 'imposed-load', file_name = 'LC11.txt',
            #    psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0
            # ),
            #---------------------------------------------------------
            # shrinkage
            #---------------------------------------------------------
            #LC12:
           # LC(name='T_schwinden', category='imposed-load', file_name='LC12.txt',
            #    psi_0=0.8, psi_1=0.7, psi_2=0.5,
           #     ),                 
    ]

    #--------------------------------------------------------
    # strength characteristics (design) assumed
    # (specimens thickness = 40 cm; internal lever 35 cm,  specimen width = 100 cm; 8 layers Q142/142-CCE-38)
    #--------------------------------------------------------
    # tensile strength [kN/m]
    #
    n_0_Rdt = n_90_Rdt = 36960.  # [kN/m] # 2*3300*140/100 kN

    # bending strength [kNm/m]
    #
    m_0_Rd = m_90_Rd = 161.7 # [kNm/m] # 0.35*3300*140*/100 kNm/m

    # compressive strength [kN/m]
    # (C3-B2-HF-1-140-5 with characteristic value; f_ck = 100 MPa)
    # (design value; f_cd = 56,7 MPa)
    #
    n_Rdc = 22680.  # = 56,7 MPa * (100 cm * 40 cm) * 0.1

    print 'design values calculated by hand for strength characteristics'

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
#    lct_Q.plot_nm_interaction(save_fig_to_file='nm_interaction_LC1-14')
#    lct_Q.plot_nm_interaction(add_max_min_nm_from_file='max_min_nm_arr_LC15-18', save_fig_to_file='nm_interaction_LC1-18')

    #--------------------------------------------------------------
    # interaction plot of material usage 'eta_nm' (utilization ratio)
    #--------------------------------------------------------------
    #
#   lct_Q.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_LC1-14')
#   lct_Q.plot_eta_nm_interaction(add_max_min_eta_nm_from_file='max_min_eta_nm_arr_LC15-18', save_fig_to_file='eta_nm_interaction_LC1-18')

    #--------------------------------------------------------------
    # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
    # (surrounding values of all loading cases)
    #--------------------------------------------------------------
    #
#   lct_Q.plot_assess_value('eta_nm_tot', add_assess_values_from_file='eta_nm_tot_LC15-18')

    #--------------------------------------------------------------
    # brows the loading case combinations within an interactive table view
    #--------------------------------------------------------------
    lct_Q.configure_traits()
