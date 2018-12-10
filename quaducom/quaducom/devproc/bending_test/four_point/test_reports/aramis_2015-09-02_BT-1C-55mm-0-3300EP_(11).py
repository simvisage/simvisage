'''
Created on Jul 1, 2014

Script evaluating the mxn interaction between
tensile and bending load
'''

if __name__ == '__main__':

    from .exp_bt_4pt_aramis2d import ExpBT4PTAramis2d
    from matresdev.db.simdb import SimDB
    from aramis_cdt import AramisUI
    simdb = SimDB()
    import os
    import numpy as np

    #-------------------------------------------
    # list of available aramis files and resolution keys
    #-------------------------------------------
    #                BT-1C-55mm-0-3300EP-V2_S3P2(11)-Aramis2d-Xf15s3-Yf10s1,
    #                BT-1C-55mm-0-3300EP-V2_S3P2(11)-Aramis2d-Xf15s13-Yf15s13,
    #                BT-1C-55mm-0-3300EP-V2_S4P2(13)-cyc-Aramis2d-Xf15s13-Yf15s13,
    #                BT-1C-55mm-0-3300EP-V2_S4P2(13)-cyc-Aramis2d-Xf15s3-Yf10s1,
    #                BT-1C-55mm-0-3300SBR-V2_S4P2(14)-Aramis2d-Xf15s13-Yf15s13,
    #                BT-1C-55mm-0-3300SBR-V2_S4P2(14)-Aramis2d-Xf15s3-Yf10s1,
    #                BT-1C-55mm-0-3300SBR-V3_S1P1(4)-Aramis2d-Xf15s13-Yf15s13,
    #                BT-1C-55mm-0-3300SBR-V3_S1P1(4)-Aramis2d-Xf15s3-Yf10s1,
    #                BT-1C-55mm-0-3300SBR-V3_S2P1(12)-cyc-Aramis2d-Xf15s13-Yf15s13,
    #                BT-1C-55mm-0-3300SBR-V3_S2P1(12)-cyc-Aramis2d-Xf15s3-Yf10s1

    #-------------------------------------------
    # select test file and resolution key
    #-------------------------------------------
    test_file = 'BT-1C-55mm-0-3300EP-V2_S3P2(11)-Aramis2d.DAT'

    test_file_path = os.path.join(simdb.exdata_dir,
                                  'bending_tests', 'four_point',
                                  '2015-09-02_BT-1C-55mm-0-3300SBR_cyc-Aramis2d',
                                  )

    res_key = 'Xf15s3-Yf10s1'
#     res_key = 'Xf15s13-Yf15s13'

    #-------------------------------------------
    # specify parameters for aramis_cdt
    #-------------------------------------------
    start_time_aramis = 0.
    crack_detaction_step = 162  # last step before rupture
    integ_radius = 3  # radius used for crack detection
    integ_radius_crack = 5  # radius used for calculation of displacement jumps (=crack width)
    dd_ux_avg_threshold = -1e-07
    dd_ux_threshold = -1e-07
    ddd_ux_avg_threshold = 0.00001
    ddd_ux_threshold = 0.00001

    #-------------------------------------------
    # specify parameters for 2d-scaling of coordinates (calculation see below)
    #-------------------------------------------
    n_px_between_load_introduction_points = 725  # [px] number of px between load introduction points for given png-file
    l_real = 350  # [mm] real distance between the load introduction points (=350 mm)

    #-------------------------------------------
    # setup experiment run with aramis data
    #-------------------------------------------
    e_run = ExpBT4PTAramis2d(data_file=os.path.join(test_file_path, test_file),
                       aramis_resolution_key=res_key,
                       start_time_aramis=start_time_aramis,
                       integ_radius=integ_radius,
                       integ_radius_crack=integ_radius_crack,
                       ddd_ux_avg_threshold=ddd_ux_avg_threshold,
                       ddd_ux_threshold=ddd_ux_threshold,
                       dd_ux_avg_threshold=dd_ux_avg_threshold,
                       dd_ux_threshold=dd_ux_threshold,
                       delta_t_aramis=0,
                       crack_detaction_current_step=False,
                       crack_detaction_step=crack_detaction_step
                       )

    e_run.process_source_data()
    print('crack filter', e_run.crack_filter_avg)

    #-------------------------------------------
    # non-deformed coordinates of the measuring field
    #-------------------------------------------
    x = e_run.aramis_field_data.x_arr_0
    y = e_run.aramis_field_data.y_arr_0
    z = e_run.aramis_field_data.z_arr_0

    #-------------------------------------------
    # dimensions of the measuring field (as read in from file)
    #-------------------------------------------
    lx_bottom = x[0, -1] - x[0, 0]
    lx_top = x[-1, -1] - x[-1, 0]
    ly_left = y[0, 0] - y[-1, 0]
    ly_right = y[0, -1] - y[-1, -1]
    lx_avg = (lx_bottom + lx_top) / 2.
    ly_avg = (ly_right + ly_left) / 2.
    print('Measuring field length: lx_bottom ', lx_bottom)
    print('Measuring field length: lx_top ', lx_top)
    print('Measuring field length: ly_left ', ly_left)
    print('Measuring field length: ly_right ', ly_right)
    print('Measuring field length: lx_avg ', lx_avg)
    print('Measuring field length: ly_avg ', ly_avg)

    #-------------------------------------------
    # resolution used for aramis evaluation
    #-------------------------------------------
    res_key_list = res_key.split('-')
    X_list = res_key_list[0].split('Xf')[1].split('s')
    Y_list = res_key_list[1].split('Yf')[1].split('s')
    Xf = int(X_list[0])
    Xs = int(X_list[1])
    Yf = int(Y_list[0])
    Ys = int(Y_list[1])
    print('res_key_list', res_key_list)
    print('X_list', X_list)
    print('Y_list', Y_list)

    #-------------------------------------------
    # average pixel size
    #-------------------------------------------
    dist_x = x[:, 1:] - x[:, :-1]
    lenght_px_x = np.mean(dist_x / Xs)
    dist_y = y[:, 1:] - y[:, :-1]
    lenght_px_y = np.mean(dist_y / Ys)
    print('Mean, std, min, max of facets center distance [mm]: x-direction', np.mean(dist_x), np.std(dist_x), np.min(dist_x), np.max(dist_x))
    print('Mean, std, min, max of facets center distance [mm]: y-direction', np.mean(dist_y), np.std(dist_y), np.min(dist_y), np.max(dist_y))
    print('Mean and std of pixel size [mm]: x-direction', np.mean(dist_x / Xs), np.std(dist_x / Xs))
    print('Mean and std of pixel size [mm]: y-direction', np.mean(dist_y / Ys), np.std(dist_y / Ys))

    #-------------------------------------------
    # calculate 2d-caling factor as ratio of read-in aramis coordinates and real length scale of tested specimens
    #-------------------------------------------
    l_default = lenght_px_x * n_px_between_load_introduction_points  # [mm]
    scale_data_factor = l_real / l_default
    lx_avg_real = lx_avg * scale_data_factor
    print('l_default', l_default)
    print('l_real ', l_real)
    print('scale_data_factor ', scale_data_factor)
    print('lx_avg_real', lx_avg_real)

    #-------------------------------------------
    # apply scale factor 2d to coordinates, displacements and crack field arrays
    #-------------------------------------------
    e_run.aramis_field_data.scale_data_factor = scale_data_factor

    #-------------------------------------------
    # cut off the boundary cells that can not be evaluated as they exceed the integ_radius-zone
    #-------------------------------------------
    left_i = e_run.aramis_field_data.left_i
    e_run.aramis_field_data.left_i = integ_radius_crack
    right_i = e_run.aramis_field_data.right_i
    e_run.aramis_field_data.right_i -= integ_radius_crack
    print('left_i', left_i)
    print('right_i', right_i)
    e_run.aramis_field_data.top_j += 9  # NOTE: cutoff 9 cells from top surface as purious displacement jumps are detected

    #-------------------------------------------
    # scale gylphs for mlab.points3d
    #-------------------------------------------
    overlap_factor = 1.25  # cubes overlap for smoother display without gaps
    range_l = max(x.shape[0], x.shape[1])
    range_t = min(x.shape[0], x.shape[1])
    dist_l = lx_avg * scale_data_factor
    dist_t = ly_avg * scale_data_factor
    glyph_l_length = dist_l / range_l
    glyph_l_length *= overlap_factor
    glyph_t_length = dist_t / range_t
    glyph_t_length *= overlap_factor
    print('range_l', range_l)
    print('range_t', range_t)
    print('dist_l', dist_l)
    print('dist_t', dist_t)
    print('glyph_l_length', glyph_l_length)
    print('glyph_t_length', glyph_t_length)

    #-------------------------------------------
    # start aramis user interface
    #-------------------------------------------
    aui = AramisUI(aramis_info=e_run.aramis_info,
                   aramis_data=e_run.aramis_field_data,
                   aramis_cdt=e_run.aramis_cdt,
                   glyph_x_length=glyph_l_length,
                   glyph_y_length=glyph_t_length,
                   glyph_z_length=0.0,
                   glyph_x_length_cr=glyph_l_length * 6.,  # extrude the plot_var in the out-of-plane direction
                   glyph_y_length_cr=glyph_l_length,
                   glyph_z_length_cr=glyph_t_length,
                   warp_factor=1.0)

    aui.configure_traits()

