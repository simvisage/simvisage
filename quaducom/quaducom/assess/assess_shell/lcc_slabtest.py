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
    # define directory:
    #------------------------
    data_dir = os.path.join(simdb.simdb_dir,
                            'simdata',
                            'input_data_slabtest',
                            'slabtest_125x125x6cm_quarter_lineload'
                            )

    #------------------------
    # define loading cases:
    #------------------------

    lc_list = [
                     # LC1:
                     LC(name='Q', category='dead-load', file_name='LC1.txt'
                        ),
              ]

    # LCCTable for imposed loads
    #
    lct = LCCTableULS(data_dir=data_dir,
                      reader_type='InfoCAD',
                      lc_list=lc_list,
                      k_alpha_min=False,  # simplification: use the minimum value for k_alpha on the resistance side
                      show_lc_characteristic=True
                      )

    #--------------------------------------------------------------
    # 'combi_arr': array with indices of all loading case combinations
    #--------------------------------------------------------------
    #
    print 'lct.combi_arr', lct.combi_arr.shape
#    np.savetxt('combi_arr_wo_temp_LCs', lct_Q.combi_arr, delimiter = ';')

    #--------------------------------------------------------------
    # nm-interaction plot (normal force - bending moment)
    #--------------------------------------------------------------
    #
#    lct.plot_nm_interaction(save_fig_to_file='nm_interaction_LC22')

    #--------------------------------------------------------------
    # interaction plot of material usage 'eta_nm' (Ausnutzungsgrad)
    #--------------------------------------------------------------
    #
#    lct.plot_eta_nm_interaction(save_fig_to_file='eta_nm_interaction_LC22')

    #--------------------------------------------------------------
    # plot of structure with color indication of material usage 'eta_nm' (Ausnutzungsgrad)
    # (surrounding values of all loading cases)
    #--------------------------------------------------------------
    #
#    lct.plot_assess_value()

    #--------------------------------------------------------------
    # brows the loading case combinations within an interactive table view
    #--------------------------------------------------------------
    lct.configure_traits()



