'''
Created on 31.07.2012

@author: alexander
'''

from etsproxy.traits.api import \
    Float, Instance, Array, Int, Property, cached_property, on_trait_change, Bool, \
    HasTraits, File, Event, Trait, Str

from etsproxy.traits.ui.api import \
    View, Item, FileEditor, HSplit, Group, VSplit, \
    Handler

from etsproxy.traits.ui.menu import \
    Action, CloseAction, HelpAction, Menu, \
    MenuBar, NoButtons, Separator, ToolBar

import numpy as np

from mathkit.mfn import MFnLineArray

from sig_fl_calib import SigFlCalib

import pylab as p



if __name__ == '__main__':

    import pylab as p

    #------------------------------------------------
    # 1) CALIBRATION:
    # get 'eps_t' and the parameter of the effective 
    # crack bridge function 'var_a' for a given 'eps_c'
    #------------------------------------------------
    #
    print '\n'
    print 'setup SigFlCalib'
    print '\n'
    sig_fl_calib = SigFlCalib( # concrete strength after 9 days
                               #
                               f_ck = 49.9,

                               # measured strain at bending test rupture (0-dir)
                               #
                               eps_c_fail = 3.3 / 1000.,

                               # measured value in bending test [kNm]
                               # value per m: M = 5*3.49
                               #
                               M_fail = 3.49,
                           
                               # values for experiment beam with width = 0.20 m
                               #
                               width = 0.20,
                               n_roving = 23.,

#                               # values per m
#                               #
#                               width = 1.0,
#                               n_roving = 120,

                               # define shape of the crack-bridge law ('linear', 'bilinear' or 'quadratic')
                               #
                               calc_mode = 'calib',
#                               calib_config = 'linear',
#                               calib_config = 'quadratic',
                               calib_config = 'cubic',

                               # define shape of the concrete stress-strain-law ('block', 'bilinear' or 'quadratic')
                               #
                               sig_c_config = 'quadratic'
#                               sig_c_config = 'bilinear'
                              
                              )

    print '1 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

    sig_fl_calib.calib_sig_t_mfn()
    print 'eps_c_fail', sig_fl_calib.eps_c_fail
    print 'eps_t_fail', sig_fl_calib.eps_t_fail
    u_sol = sig_fl_calib.u_sol
#    sig_fl_calib.plot_sig_t_mfn( u_sol )
    
    
    #------------------------------------------------
    # 2) EVALUATION / VALIDATION:
    # get 'eps_lo', 'esp_up' for given/calibrated cb-law 
    #------------------------------------------------

    print '2 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    
    eps_lo = 0.005
    eps_up = -0.002
    
    N_internal, M_internal = sig_fl_calib.eval_N_M( eps_lo, eps_up )

    # select the right stress_mode for calculation
    #
    u_sol = [ abs(sig_fl_calib.eps_lo), abs(sig_fl_calib.eps_up) ]
    max_sig = sig_fl_calib.get_sig_max( u_sol )     
    print 'eps_up', abs( sig_fl_calib.eps_up )
    print 'eps_lo', abs( sig_fl_calib.eps_lo )
    print 'eps_c_fail', sig_fl_calib.eps_c_fail
    print 'eps_t_fail', sig_fl_calib.eps_t_fail
    print 'max_sig', sig_fl_calib.get_sig_max( u_sol )                      
    
    # plot response:
    # graph shows sig_comp at the height of the textil layer in [MPa] 
    #
#    sig_fl_calib.plot_sig_t_mfn( u_sol )


    #------------------------------------------------
    # 3 construct failure envelope for a given value of eps_n 
    #------------------------------------------------
    print '3 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    
    
    n_A = 20
    n_B = 20
    n_C = 20

    eps_c_fail = sig_fl_calib.eps_c_fail
    eps_t_fail = sig_fl_calib.eps_t_fail

    eps_lo_arr_A = np.linspace( eps_t_fail,   eps_t_fail, n_A )
    eps_up_arr_A = np.linspace( eps_t_fail, - eps_c_fail, n_A )
    
    eps_lo_arr_B = np.linspace( eps_t_fail,           0., n_B )
    eps_up_arr_B = np.linspace(-eps_c_fail, - eps_c_fail, n_B )
    
    eps_lo_arr_C = np.linspace(         0.,       -0.002, n_C )
    eps_up_arr_C = np.linspace(-eps_c_fail,       -0.002, n_C )
        
    eps_lo_arr = np.hstack([ eps_lo_arr_A, eps_lo_arr_B, eps_lo_arr_C ])
    eps_up_arr = np.hstack([ eps_up_arr_A, eps_up_arr_B, eps_up_arr_C ])

#    just one diagram with a specific value for n_layers and thickness is plotted
# ------------------------------------------------------------------------------------


#    eval_N_M_vectorized = np.frompyfunc( sig_fl_calib.eval_N_M, 2, 2 )
#    N_arr, M_arr = eval_N_M_vectorized( eps_lo_arr, eps_up_arr )
#    
#    from mpl_toolkits.axes_grid.axislines import SubplotZero
#    
#    fig = p.figure(1)
#
#    ax = SubplotZero(fig, 111)
#    fig.add_subplot(ax)
#
##    for direction in ["xzero", "yzero"]:
##        ax.axis[direction].set_axisline_style("-|>")
##        ax.axis[direction].set_visible(True)
##
##    for direction in ["left", "right", "bottom", "top"]:
##        ax.axis[direction].set_visible(False)
#
##   plot_list = [for (M, N) in zip(M_list, N_list) ]
#    ax.plot( M_arr, N_arr )
#
##    ax = p.gca()
#
#    ax.set_ylim(ax.get_ylim()[::-1])
#
#
#
#    p.show()

#    various diagrams with a specific value for n_layers and thickness are plotted
# ------------------------------------------------------------------------------------

    #n_layers_arr = np.array([8,9,10,11,12,13,14,15,16,17,18])
    n_layers_arr = np.array([8,12,18])
    thickness_arr = np.array ([0.04,0.06])
    #thickness_arr = np.array ([0.06])
    
    N_list = []
    M_list = []
    for n in n_layers_arr:   
        for k in thickness_arr:
            
            sig_fl_calib.n_layers = n        
            sig_fl_calib.thickness = k
            eval_N_M_vectorized = np.frompyfunc( sig_fl_calib.eval_N_M, 2, 2 )
            N_arr, M_arr = eval_N_M_vectorized( eps_lo_arr, eps_up_arr )
            N_list.append( N_arr )
            M_list.append( M_arr )
        


    from mpl_toolkits.axes_grid.axislines import SubplotZero
    
    fig = p.figure(1)

    ax = SubplotZero(fig, 111)
    fig.add_subplot(ax)


#    plot_list = [for (M, N) in zip(M_list, N_list) ]

#    plt_list = []
#    for n in range( len( M_list )):
#        plt_list.append(M_list[n])
#        plt_list.append(N_list[n])
#    ax.plot( plt_list[0] )
    #ax.plot (M_list[0], N_list[0],M_list[1], N_list[1],M_list[2], N_list[2],M_list[3], N_list[3],M_list[4], N_list[4],M_list[5], N_list[5],M_list[6], N_list[6],M_list[7], N_list[7],M_list[8], N_list[8],M_list[9], N_list[9],M_list[10], N_list[10],M_list[11], N_list[11],M_list[12], N_list[12],M_list[13], N_list[13],M_list[14], N_list[14],M_list[15], N_list[15],M_list[16], N_list[16],M_list[17], N_list[17],M_list[18], N_list[18],M_list[19], N_list[19],M_list[20], N_list[20],M_list[21], N_list[21])
    ax.plot (M_list[0], N_list[0],M_list[1], N_list[1],M_list[2], N_list[2],M_list[3], N_list[3],M_list[4], N_list[4],M_list[5], N_list[5])
    #for negative values at the top
    #ax.set_ylim(ax.get_ylim()[::-1])
    
    from pylab import *
    grid(True)
    from mpl_toolkits.axes_grid.axislines import SubplotZero
    import matplotlib.pyplot as plt
    

    for direction in ["xzero", "yzero"]:
        ax.axis[direction].set_axisline_style("-|>")
        ax.axis[direction].set_visible(True)

    for direction in ["left", "right", "bottom", "top"]:
        ax.axis[direction].set_visible(False)
    plt.axis([0, 6, 250 ,-700])
    # change thickenss of a line in the diagram
    plt.plot(M_list[0], N_list[0],M_list[1], N_list[1],M_list[2], N_list[2],M_list[3], N_list[3],M_list[4], N_list[4],M_list[5], N_list[5],color='black')
    plt.plot(M_list[3], N_list[3],color='blue', linewidth=2.0)
    #plt.tick_params(axis='x',direction='out', length=6, width=2, colors='r')
    plt.title('Interaction Diagram M/N')
    #plt.xlabel('Normal Force')
    #plt.ylabel('Moment')
    
    p.show()
    
#
#   normed Ny-My interaction diagram
# ------------------------------------------------------------------------------------

#    
#    n_layers_arr = np.array([8,12,18])
#    thickness_arr = np.array ([0.04,0.06])
#    #thickness_arr = np.array ([0.06])
#    
#
#    Ny_list = []
#    My_list = []
#    for n in n_layers_arr:   
#        for k in thickness_arr:
#            
#            sig_fl_calib.n_layers = n        
#            sig_fl_calib.thickness = k
#            eval_N_M_vectorized = np.frompyfunc( sig_fl_calib.eval_N_M, 2, 2 )
#            N_arr, M_arr = eval_N_M_vectorized( eps_lo_arr, eps_up_arr )
#            Ny_arr = N_arr / (1000. * sig_fl_calib.thickness * sig_fl_calib.width * sig_fl_calib.f_ck)
#            Ny_list.append( Ny_arr )
#            My_arr = M_arr / (1000. * sig_fl_calib.thickness **2. * sig_fl_calib.width  * sig_fl_calib.f_ck)
#            My_list.append( My_arr )
#            
#        
#
#
#    from mpl_toolkits.axes_grid.axislines import SubplotZero
#    
#    fig = p.figure(1)
#
#    ax = SubplotZero(fig, 111)
#    fig.add_subplot(ax)
#
#    ax.plot (My_list[0], Ny_list[0],My_list[1], Ny_list[1],My_list[2], Ny_list[2])
#    ax.set_ylim(ax.get_ylim()[::-1])
#    
#    
#  
#
#
#    p.show()

#
#
#
#
#
#
#
#
#    
#    
##    
##    kappa_arr = np.linspace( 0., 0.25, 20 )
##    print 'kappa_arr', kappa_arr
##
##    N_list = []
##    M_list = []
##    
##    for kappa in kappa_arr:
##
##        sig_fl_calib.set( eps_n = eps_n,
##                          kappa = kappa,
##                          sig_t_mfn = sig_t_mfn, 
##                          calc_mode = 'eval' )
##
##        u_sol = [ abs(sig_fl_calib.eps_lo), abs(sig_fl_calib.eps_up) ]
##        print 'u_sol', u_sol
##        
##        N, M = sig_fl_calib.get_lack_of_fit( u_sol, 0, 0 )
##    
###        max_sig = sig_fl_calib.get_sig_max( u_sol )                                   
###        print 'eps_up', abs(sig_fl_calib.eps_up)
###        print 'eps_lo', abs(sig_fl_calib.eps_lo)
###        print 'max_sig', sig_fl_calib.get_sig_max( u_sol )    
##
##        N_list.append( N )
##        M_list.append( M )
#
#    N_idx_max = np.argmax( N_list ) 
#    M_idx_max = np.argmax( M_list ) 
#    print 'N_idx_max', N_idx_max
#    print 'M_idx_max', M_idx_max
#
#    kappa_N_max = kappa_arr[ N_idx_max ]
#    kappa_M_max = kappa_arr[ M_idx_max ]
#    print 'kappa_N_max', kappa_N_max
#    print 'kappa_M_max', kappa_M_max
#
#    N_max = N_list[ N_idx_max ]
#    M_max = M_list[ M_idx_max ]
#    print 'N_max', N_max
#    print 'M_max', M_max
#
#    print 'kappa_arr.shape', kappa_arr.shape
#    print 'M_list.shape', len(M_list)
#    print 'M_arr', np.array(M_list)
