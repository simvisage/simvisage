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


#    various diagrams with a specific value for n_layers and thickness are plotted
# ------------------------------------------------------------------------------------

    import pylab as p
    from mpl_toolkits.axes_grid.axislines import SubplotZero

    fig = p.figure(1)
    fig.set_facecolor('w')

    ax = SubplotZero(fig, 111, axisbg='w')
    fig.add_subplot( ax )

    for direction in ["xzero", "yzero"]:
#        ax.axis[direction].set_axisline_style("-|>")
        ax.axis[direction].set_visible(True)

    for direction in ["left", "right", "bottom", "top"]:
        ax.axis[direction].set_visible(False)

    
    nu_list = []
    mu_list = []
    N_list = []
    M_list = []
    rho_list = []

    f_ck = sig_fl_calib.f_ck
    b = sig_fl_calib.width
    A_roving = sig_fl_calib.A_roving
    n_rovings = sig_fl_calib.n_rovings

    # constant 'n_layers':
    #
#    n_layers_list  = [  12,   12,    12,   12,   12]
#    thickness_list = [0.02, 0.04,  0.06, 0.08, 0.10]
    
    # constant 'thickness':
    #
    n_layers_list  = [   4,    8,    12,   16,   20]
    thickness_list = [0.06, 0.06,  0.06, 0.06, 0.06]

    # constant 'rho_tex':
    #
#    n_layers_list  = [   6,    8,    12,   18,   24]
#    thickness_list = [0.03, 0.04,  0.06, 0.09, 0.12]

    zip_list = zip( n_layers_list, thickness_list  )

    # varying layers for two different heights':
    #
#    n_layers_list = np.array([8,12,18])
#    thickness_list = np.array ([0.04,0.06])
#    zip_list =[]
#    for n in n_layers_list:
#        for t in thickness_list:
#            zip_list.append((n,t)) 


    for n,t in zip_list:
        #
        sig_fl_calib.n_layers = n        
        sig_fl_calib.thickness = t
        #
        eval_N_M_vectorized = np.frompyfunc( sig_fl_calib.eval_N_M, 2, 2 )
        N_arr, M_arr = eval_N_M_vectorized( eps_lo_arr, eps_up_arr )
        #
        c1 = 1000. * t* b * f_ck
        c2 = 1000. * t**2. * b * f_ck 
        #
        # @todo: why does an error occure when multiplying an array by a float value (caused by vectorized)
        nu_arr = N_arr / np.array([ c1 ]) 
        mu_arr = M_arr / np.array([ c2 ])
        rho_tex = n * A_roving * n_rovings / c1
        #
        N_list.append( N_arr )
        M_list.append( M_arr )
        nu_list.append( nu_arr )
        mu_list.append( mu_arr )
        rho_list.append( rho_tex )
        
        # plot normalized values nu and mu
        #
#        if t == 0.06 and n == 12:
#            p.plot( mu_arr, nu_arr,
#                    color='blue', 
#                    linewidth=2.0)
#        else:
#            p.plot( mu_arr, nu_arr )

        # plot absolute values for N and M
        #
        if t == 0.06 and n == 12:
            p.plot( M_arr, N_arr,
                    color='blue', 
                    linewidth=2.0)
        else:
            p.plot( M_arr, N_arr )

 
    # @todo: format plot and axes (legend, grid, font etc.)
#    p.legend( rho_list )
#    p.grid()
#    legend( (l2, l4), ('oscillatory', 'damped'), 'upper right', shadow=True)
#    p.axis([0, 6, 250 ,-700])
#    p.xticks(0.5*np.arange(12))
#    ax.set_xticks([0., 1., 2., 3., 4., 5., 6.])
    ax.set_ylim(ax.get_ylim()[::-1])

    from matplotlib.ticker import MaxNLocator
#    ax.xaxis.set_major_locator(MaxNLocator(10))
    ax.yaxis.set_major_locator(MaxNLocator(10))
        
#    p.tick_params(axis='x',direction='out', length=6, width=2, colors='r')
    
    p.title('Interaction Diagram M/N')
#    p.xlabel('moment [kNm]')
#    p.ylabel('normal force [kN]')

    p.xlabel('mu')
    p.ylabel('nu')

#    from matplotlib import rcParams
#    rcParams['text.usetex']=True
#    rcParams['text.latex.unicode']=True
#    p.xlabel(r"\textbf{$\mu$}")
#    p.ylabel(r"\textbf{$\nu$}")
    
    p.show()
    
