if __name__ == '__main__':

    from matresdev.db.simdb import \
        SimDB

    import numpy as np

    import os

    from .lcc_table import LCCTableULS, LC, LCCTableSLS

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    # dimensioning of barrel shells for all loading cases
    # (incl. superposition of temperature loading cases)
    #
    data_dir = os.path.join(simdb.simdb_dir,
                            'simdata',
                            'input_data_hangar',
                            'bearing_sector',
                            )

    #------------------------
    # define loading cases:
    #------------------------

    #----------------------------------------------------------------------
    # (c) dimensioning with loads according to DIN EN 1991-1-3 (snow) and DIN EN 1991-1-4 (wind)
    #----------------------------------------------------------------------
    #---------------------------------------------------------
    # "staendige und voruebergehende Bemessungssitauation":
    #---------------------------------------------------------
    #
    lc_list_Q = [

        #------------------------------------------------------------------
        # dead load
        #-------------------------------------------------------------
        # LC1:
        LC(name='g_exclP', category='dead-load', file_name='LC1.txt',
           ),
        
        #LC10:
        LC(name='P', category='dead-load', file_name='LC10.txt', gamma_unf = 1.00,
           ),


        #-------------------------------------------------------------
        # snow
        # load cases: 's_sym.', 's_unsym.', 's_unsym_dreieck',
        #-------------------------------------------------------------
        # LC2:
        LC(name='s_sym', category='imposed-load', file_name='LC2.txt',
               exclusive_to=[
                    's_asym_dreieck', 's_asym'],
                psi_0=0.5, psi_1=0.2, psi_2=0.0
           ),
         # LC3:
        LC(name='s_unsym', category='imposed-load', file_name='LC3.txt',
               exclusive_to=[
                    's_sym', 's_asym_dreieck',],
                psi_0=0.5, psi_1=0.2, psi_2=0.0
           ),
          # LC4:
        LC(name='s_asym_dreieck', category='imposed-load', file_name='LC4.txt',
               exclusive_to=[
                    's_asym', 's_sym'],
                psi_0=0.5, psi_1=0.2, psi_2=0.0
           ), 
            #---------------------------------------------------------
            # wind load
            #load cases: 'w_X', 'w_Y+', w_Y-',
            #---------------------------------------------------------
            # LC5:
            LC(name='w_X', category='imposed-load', file_name='LC5.txt',
               exclusive_to=[
                   'w_Y+', 'w_Y-'],
               psi_0=0.6, psi_1=0.2, psi_2=0.0
               ), 
            # LC6:
            LC(name='w_Y+', category='imposed-load', file_name='LC6.txt',
               exclusive_to=[
                   'w_X', 'w_Y-'],
               psi_0=0.6, psi_1=0.2, psi_2=0.0
               ), 
            # LC7:
            LC(name='w_Y-', category='imposed-load', file_name='LC7.txt',
               exclusive_to=[
                   'w_X', 'w_Y+'],
               psi_0=0.6, psi_1=0.2, psi_2=0.0
               ), 
          #---------------------------------------------------------
            # temperature
         #---------------------------------------------------------
 
            # LC11:
            LC(name = 'T_sommer', category = 'imposed-load', file_name = 'LC11.txt',
                psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0
             ),
            #---------------------------------------------------------
            # shrinkage
            #---------------------------------------------------------
            #LC12:
            LC(name='T_schwinden', category='imposed-load', file_name='LC12.txt',
                psi_0=0.8, psi_1=0.7, psi_2=0.5,
                ), 
    ]

    #--------------------------------------------------------
    # strength characteristics (design) assumed
    #(specimens thickness = 0.65 m; internal lever 0.52m,  specimen width = 1m
    #--------------------------------------------------------
   ### tensile strength [kN/m]###
    # elected c-bar: carbopree gamma_cfk = 1.3, f_cfk,d = 2040 MN/m^2

    # A_cfk,0,vorh = 184,8 cm^2/m 
    n_0_Rdt = 37699. 
 
    # A_cfk,90,vorh = 55,44 cm^2/m
    n_90_Rdt =11310. #[kN/m]
     
    
    ### bending strength [kNm/m] ###
    ### aus iterativer QS-Bemessung mit epsilon,c=2,8 promille und epsilon,t= 9,74 promille,
    #   C90/105, f_cd= 51 MN/m^2, F_cd = 4,84 MN/m
    m_0_Rd = 2522.
    m_90_Rd = 2522.  
    

    ### compressive strength [kN/m] ###
    # C90/105, f_cd= 51 MN/m^2, high utilised = 0.45 m
    n_Rdc = 22950.  
    
    
    print('design values calculated by hand for strength characteristics')

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
    print('lct_Q.combi_arr', lct_Q.combi_arr.shape, '\n')
#        np.savetxt('combi_arr_wo_temp_LCs', lct_Q.combi_arr, delimiter=';')

    #--------------------------------------------------------------
    # nm-interaction plot (normal force - bending moment)
    #--------------------------------------------------------------
    #
    lct_Q.plot_nm_interaction(save_fig_to_file='nm_interaction_bearing')
   # lct_Q.plot_nm_interaction(add_max_min_nm_from_file='max_min_nm_arr_LC15-18', save_fig_to_file='nm_interaction_LC1-18')

    #--------------------------------------------------------------
    # interaction plot of material usage 'eta_nm' (utilization ratio)
    #--------------------------------------------------------------
    #
    #lct_Q.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_bearing_sv'),
    #lct_Q.plot_eta_nm_interaction(add_max_min_eta_nm_from_file='max_min_eta_nm', save_fig_to_file='eta_nm_interaction_bearing')

    #--------------------------------------------------------------
    # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
    # (surrounding values of all loading cases)
    #--------------------------------------------------------------
    #
   # lct_Q.plot_assess_value('eta_nm_tot', add_assess_values_from_file='eta_nm_tot'),

    #--------------------------------------------------------------
    # brows the loading case combinations within an interactive table view
    #--------------------------------------------------------------
    lct_Q.configure_traits()