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
                               f_ck = 55.7,

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

                               # define shape of the crack-bridge law ('linear', 'bilinear' or 'quadratic')
                               #
#                               calib_config = 'root2',
#                               calib_config = 'quadratic',
#                               calib_config = 'linear',  
#                               calib_config = 'quadratic_eps_max', 
#                               calib_config = 'quadratic_monoton', 
#                               calib_config = 'quadratic_TT',  
#                               calib_config = 'plastic', 
                               calib_config = 'cubic',

                               # define shape of the concrete stress-strain-law ('block', 'bilinear' or 'quadratic')
                               #
                              sig_c_config = 'quadratic'
#                               sig_c_config = 'bilinear'
#                               sig_c_config = 'block'
                              )

    
    print '1 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

    sig_fl_calib.calib_sig_t_mfn()
    u_sol = sig_fl_calib.u_sol
    max_sig = sig_fl_calib.get_sig_max( u_sol )     
    print 'eps_c_fail', sig_fl_calib.eps_c_fail
    print 'eps_t_fail', sig_fl_calib.eps_t_fail
    print 'max_sig', sig_fl_calib.get_sig_max( u_sol )                      
#    sig_fl_calib.plot_sig_t_mfn( u_sol )
    
    
    #------------------------------------------------
    # 2) EVALUATION / VALIDATION:
    # get 'eps_lo', 'esp_up' for given/calibrated cb-law 
    #------------------------------------------------

    print '2 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    
    # reproduce the forces for the calibration test:
    #
    eps_lo =   sig_fl_calib.eps_t_fail
    eps_up = - sig_fl_calib.eps_c_fail
    N_internal, M_internal = sig_fl_calib.eval_N_M( eps_lo, eps_up )

    
    #------------------------------------------------
    # 3 construct failure envelope for a given value of eps_n 
    #------------------------------------------------
    print '3 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    
    n_A = 20
    n_B = 20
    n_C = 20

    # mainly compression strain:
    if  sig_fl_calib.sig_c_config== 'quadratic_2':    
        eps_compression = -(2. + 0.085 * ((sig_fl_calib.f_ck + 8.) - 50)**0.53) / 1000.
        
    elif sig_fl_calib.sig_c_config== 'bilinear':   
    
        eps_compression = -(1.75 + 0.55 * (((sig_fl_calib.f_ck + 8) - 50.) / 40.)) / 1000.
        
    elif sig_fl_calib.sig_c_config== 'quadratic':   
        
        eps_compression = -min(0.7 *( sig_fl_calib.f_ck + 8. ) ** 0.31, 2.8) / 1000. #eps_c1
        
    elif sig_fl_calib.sig_c_config== 'block':   
        
        eps_compression = -0.002
                
    eps_c_fail = 3.3/1000.
    eps_t_fail = sig_fl_calib.eps_t_fail
    print'eps_t_fail ',eps_t_fail
    # mainly tension
    
    eps_lo_arr_A1 = np.linspace( eps_t_fail,   eps_t_fail, n_A )
    eps_up_arr_A1 = np.linspace( eps_t_fail,            0, n_A )
    
    #  flexion2
    
    eps_lo_arr_A2 = np.linspace(eps_t_fail,   eps_t_fail, n_A )
    eps_up_arr_A2 = np.linspace(         0, - eps_c_fail, n_A )
    
    #  flexion3
    
    eps_lo_arr_B = np.linspace( eps_t_fail,           0., n_B )
    eps_up_arr_B = np.linspace(-eps_c_fail, - eps_c_fail, n_B )
    
    # pure compression
    
    eps_lo_arr_C = np.linspace(         0.,       eps_compression, n_C )
    eps_up_arr_C = np.linspace(-eps_c_fail,       eps_compression, n_C )

#    all stress cases

    eps_lo_arr = np.hstack([ eps_lo_arr_A1,eps_lo_arr_A2, eps_lo_arr_B, eps_lo_arr_C ])
    eps_up_arr = np.hstack([ eps_up_arr_A1,eps_up_arr_A2, eps_up_arr_B, eps_up_arr_C ])

# plot stress cases
    
    eps_lo_arr_psc = np.array  ([ 0.,          0., 0.,    eps_t_fail, 0. ,eps_t_fail, 0.])
    eps_up_arr_psc = np.array  ([ 0., -eps_c_fail, 0.,   -eps_c_fail, 0. ,        0., 0.])

#    all stress cases with classification of stress cases

#    eps_lo_arr = np.hstack([ eps_lo_arr_A1,eps_lo_arr_A2, eps_lo_arr_B, eps_lo_arr_C,eps_lo_arr_psc ])
#    eps_up_arr = np.hstack([ eps_up_arr_A1,eps_up_arr_A2, eps_up_arr_B, eps_up_arr_C,eps_up_arr_psc ])
#    all stress cases without classification of stress cases    
    eps_lo_arr = np.hstack([ eps_lo_arr_A1,eps_lo_arr_A2, eps_lo_arr_B, eps_lo_arr_C ])
    eps_up_arr = np.hstack([ eps_up_arr_A1,eps_up_arr_A2, eps_up_arr_B, eps_up_arr_C ])
    
    eps_lo_arr_psc = ([ eps_lo_arr_psc])
    eps_up_arr_psc = ([ eps_up_arr_psc])
#    mainly tension

#    eps_lo_arr = np.hstack([ eps_lo_arr_A1 ])
#    eps_up_arr = np.hstack([ eps_up_arr_A1 ])
    
    #    mainly compression

#    eps_lo_arr = np.hstack([ eps_lo_arr_C ])
#    eps_up_arr = np.hstack([ eps_up_arr_C ])
    
    #    flexion2

#    eps_lo_arr = np.hstack([ eps_lo_arr_A2 ])
#    eps_up_arr = np.hstack([ eps_up_arr_A2 ])
    
     #    flexion3

#    eps_lo_arr = np.hstack([ eps_lo_arr_B ])
#    eps_up_arr = np.hstack([ eps_up_arr_B ])




# ------------------------------------------------------------------------------------

    import pylab as p
    from mpl_toolkits.axes_grid.axislines import SubplotZero
    from matplotlib.ticker import MaxNLocator
    from matplotlib.ticker import AutoMinorLocator

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
#    zip_list = zip( n_layers_list, thickness_list  )
#    legend_list = thickness_list 


    # constant 'thickness':
    #
#    n_layers_list  = [   12]
#    thickness_list = [ 0.06]
#    n_layers_list  = [   4,    8,    12,   16,   20]
#    thickness_list = [0.06, 0.06,  0.06, 0.06, 0.06]
#    zip_list = zip( n_layers_list, thickness_list  )
#    legend_list = n_layers_list 
#    
    # constant 'rho_tex':
    
    n_layers_list  = [  6,    8,    12,   18,   24]
    thickness_list = [0.03, 0.04,  0.06, 0.09, 0.12]
    zip_list = zip( n_layers_list, thickness_list  )
    legend_list = n_layers_list 

    # varying layers for two different heights':
    # use plot with M,N instead of nu,mu !
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
        # @todo: check which value to be used for normalization if sig_c = parabolic is used
        # NOTE that for this case f_cm  is used not f_ck; NOTE that for pure compression eps = 2E-3 is used
        # which does not yield f_cm in the parabolic law; for the default concrete laws for pure
        # compression the value should evaluated to 1.0!
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
        
        
#    # coordinates for visualization stress cases
#
#        # mainly compression / pure tension (2)
#        eps_lo_1 = -eps_c_fail
#        eps_up_1 = 0.
#
#        N, M = eval_N_M_vectorized( eps_lo_1, eps_up_1 )    
#        nu_arr = N /  c1 
#        mu_arr = M /  c2
#        stress_case_1_arr = [(mu_arr, nu_arr),(0.,0.)] 
######
#
        if t == 0.06 and n == 12:
                p.plot( mu_arr, nu_arr,
                        color='blue', 
                        linewidth=2.0)
        else:
                p.plot( mu_arr, nu_arr )
        # plot absolute values for N and M
#        
#        if t == 0.06 and n == 12:
#            p.plot( M_arr, N_arr,
#                    color='blue', 
#                    linewidth=2.0)
#        else:
#            p.plot( M_arr, N_arr )
    
    
    
    ######### @todo: format plot and axes (legend, grid, font etc.) #########

    ### GRID
    #
    p.grid(b=None, which='both')


    ### LEGEND
    #
    p.legend( legend_list )
#    legend( (l2, l4), ('oscillatory', 'damped'), 'upper right', shadow=True)


    ### TICKS
    #
#    p.axis([0, 6, 250 ,-700])
#    p.axis([0, 0.2, 0.5 ,-1.2])

    #p.xticks(0.5*np.arange(12))
#    ax.set_xticks([0., 1., 2., 3., 4., 5., 6.])
    ax.set_ylim(ax.get_ylim()[::-1])
    
    ax.xaxis.set_major_locator(MaxNLocator(10))
    ax.yaxis.set_major_locator(MaxNLocator(10))
    minorLocator   = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minorLocator)
#    p.tick_params(axis='both',which='major',direction='out', length=15, width=2, colors='black')
#    p.tick_params(axis='both',which='minor',direction='in', length=100, width=1, colors='black')


    ### TITEL
    #
#    p.title(r'$\nu$ -$\mu$ Interaction Diagram ',fontsize=22,
#                                verticalalignment='baseline',
#                                horizontalalignment = 'center')
    #p.title('Interaction Diagram M/N')
#    p.xlabel('moment [kNm]')
#    p.ylabel('normal force [kN]')

    
    ### LABELS
    #
    #p.xlabel(r'$\mu$',fontsize='20',
#             verticalalignment = 'top',
#             horizontalalignment = 'left')
    #p.ylabel(r'$\nu$ ', fontsize= 20)
#    from matplotlib import rcParams
#    rcParams['text.usetex']=True
#    rcParams['text.latex.unicode']=True
#    p.xlabel(r"\textbf{$\mu$}")
 #  p.ylabel(r\textbf{$\nu$}")
    
    p.show()
    