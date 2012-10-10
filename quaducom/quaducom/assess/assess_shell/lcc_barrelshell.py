


if __name__ == '__main__':

    from etsproxy.mayavi import \
        mlab

    from matresdev.db.simdb import \
        SimDB

    import os

    from lcc_table import LCCTableULS, LC

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    #---------------------------------------------
    # 2 shells: 
    # new geometry with new loading cases plus waterfilling
    #---------------------------------------------

    data_dir = os.path.join(simdb.simdb_dir,
                            'simdata', 'input_data_barrelshell')

    #------------------------
    # define loading cases:
    #------------------------
    # NOTE:
    #
    #---------------------------------------------------------
    # "staendige und voruebergehende Bemessungssitauation":
    #---------------------------------------------------------
    lc_list = [
                 # LC1:
                 LC(name = 'g', category = 'dead-load', file_name = 'LC1.txt'
                    ),
                 # LC2:
                 LC(name = 'snow', category = 'imposed-load', file_name = 'LC2.txt'
                    ),
               ]

#    #---------------------------------------------------------
#    # "aussergewoehnliche Bemessungssitauation":
#    #---------------------------------------------------------
#    #

#--------------------------------------------------------


    lct = LCCTableULS(data_dir = data_dir,
                      reader_type = 'InfoCAD',
                       lc_list = lc_list,

                       # remove only the lowest point = connection shell/column
                       # as this is a singularity of the FE-shell-model
                       #
#                       cut_z_fraction = 0.01,
#                       cut_z_fraction = 0.05, # corresponds to 50cm x 50cm
#                       cut_z_fraction = 0.10, # corresponds to 75cm x 75cm
                       cut_z_fraction = 0.15, # corresponds to 100cm x 100cm
                       show_lc_characteristic = False
                        )


    mlab.figure(figure = "SFB532Demo",
                 bgcolor = (1.0, 1.0, 1.0),
                 fgcolor = (0.0, 0.0, 0.0))

    lct.plot_geo(mlab)

    lct.plot_sr(mlab)


    mlab.show()
#    lct = LCCTableSLS( data_dir = data_dir,
#                       lc_list = lc_list,
#                       cut_z_fraction = 0.2,
#                       combination_SLS = 'rare',
##                       combination_SLS = 'freq',
##                       combination_SLS = 'perm',
##                       show_lc_characteristic = True
#                        )


#    lct.configure_traits()


#    print 'lc_arr', lct.lc_arr
#    print 'lc_list[0].sr_arr.shape[0]', lct.lc_list[0].sr_arr.shape[0]
#    print 'lc_arr.shape', lct.lc_arr.shape
#    print 'combi_arr', lct.combi_arr
#    print 'combi_arr.shape', lct.combi_arr.shape
#    print 'lcc_arr', lct.lcc_arr
