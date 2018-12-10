if __name__ == '__main__':

    from mayavi import \
        mlab

    from matresdev.db.simdb import \
        SimDB

    import numpy as np

    import os

    from .lcc_table import LCCTableULS, LC, LCCTableSLS

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    #------------------------
    # define directory:
    #------------------------
    data_dir = os.path.join(simdb.simdb_dir,
                            'simdata',
                            'input_data_slabtest',
                            'slabtest_125x125x6cm_quarter_lineload_R=0.25MN'
                            )

    #------------------------
    # define loading case:
    #------------------------
    lc_list = [LC(name='1MN-center-load', category='dead-load', file_name='LC1.txt',
                  gamma_fav=1.0, gamma_unf=1.0)]

    #--------------------------------------------------------
    # mean values for SFB-demonstrator (QS)
    #--------------------------------------------------------
    # design values for quality tests of SFB-demonstrator
    # (TT-SH2 and BT-SH4) on specimens with thickness 6 cm

    # tensile strength [kN/m]
    #
    # = 139.5 kN / 0.14 ### F_Rtm = 139.5 kN (mean value)
    n_0_Rdt = n_90_Rdt = 996.4

    # compressive strength [kN/m]
    #
    n_Rdc = 2200  # C55/67

    # bending strength [kNm/m]
    #
    m_0_Rd = m_90_Rd = 15.5  # = 3.1 / 0.20 ### M_Rm = 3.1 kNm (mean value)

    # LCCTable for imposed loads
    #
    lct = LCCTableULS(data_dir=data_dir,
                      reader_type='InfoCAD',
                      lc_list=lc_list,
                      show_lc_characteristic=False,
                      strength_characteristics={'n_0_Rdt': n_0_Rdt,
                                                'm_0_Rd': m_0_Rd,
                                                'n_Rdc': n_Rdc,
                                                'n_90_Rdt': n_90_Rdt,
                                                'm_90_Rd': m_90_Rd},
                      # NO simplification used for 'k_alpha' on the resistance
                      # side
                      k_alpha_min=False,
                      )

    #--------------------------------------------------------------
    # 'combi_arr': array with indices of all loading case combinations
    #--------------------------------------------------------------
    #
#    print 'lct.combi_arr', lct.combi_arr.shape
#    np.savetxt('combi_arr_wo_temp_LCs', lct_Q.combi_arr, delimiter = ';')

    #--------------------------------------------------------------
    # nm-interaction plot (normal force - bending moment)
    #--------------------------------------------------------------
    #
#    lct.plot_nm_interaction(save_fig_to_file='nm_interaction_slabtest')

    #--------------------------------------------------------------
    # interaction plot of material usage 'eta_nm' (utilization ratio)
    #--------------------------------------------------------------
    #
#    lct.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_slabtest')

    #--------------------------------------------------------------
    # plot of structure with color indication of material usage 'eta_nm' (utilization ratio)
    # (surrounding values of all loading cases)
    #--------------------------------------------------------------
    #
#    lct.plot_assess_value('eta_nm_tot', scale_factor=0.04)

    #--------------------------------------------------------------
    # brows the loading case combinations within an interactive table view
    #--------------------------------------------------------------
    lct.configure_traits()
