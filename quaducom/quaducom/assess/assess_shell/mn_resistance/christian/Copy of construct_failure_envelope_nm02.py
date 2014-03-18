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

import Image



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
    sig_fl_calib = SigFlCalib(# concrete strength after 9 days
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

    #------------------------------------------------
    # 2) Calibration/ EVALUATION / VALIDATION:
    # get 'eps_lo', 'esp_up' for given/calibrated cb-law 
    #------------------------------------------------

    print '2 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    
# ------------------------------------------------------------------------------------

    import pylab as p
    from mpl_toolkits.axes_grid.axislines import SubplotZero
    from matplotlib.ticker import MaxNLocator
    from matplotlib.ticker import AutoMinorLocator

# create plot

    fig = p.figure()   
    fig.set_facecolor('w')
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('')
    ax.spines['left'].set_position(('data', 0))
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['top'].set_color('none')
#    ax.spines['left'].set_smart_bounds(True)
#    ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
#    ##### for image
#    
#    im = Image.open('C:\\Users\Christian Schmitz\Desktop\Dehnungsverteilung_Interaktionsdiagramm.png')
#    height = im.size[1]
#
#    # We need a float array between 0-1, rather than
#    # a uint8 array between 0-255
#    im = np.array(im).astype(np.float) / 255
#
#    fig = p.figure()
#
#    # With newer (1.0) versions of matplotlib, you can 
#    # use the "zorder" kwarg to make the image overlay
#    # the plot, rather than hide behind it... (e.g. zorder=10)
#    fig.figimage(im, 0, fig.bbox.ymax - height)
#
#    # (Saving with the same dpi as the screen default to
#    #  avoid displacing the logo image)    
#    fig.savefig('C:\\Users\Christian Schmitz\Desktop\Dehnungsverteilung_Interaktionsdiagramm.png', dpi=80)

    
#    for direction in ["xzero", "yzero"]:
#        ax.axes[direction].set_axisline_style("-|>")
#        ax.axes[direction].set_visible(True)
#
#    for direction in ["left", "right", "bottom", "top"]:
#        ax.spines[direction].set_visible(False)

# determine frame work requirements

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


    # constant 'thickness':
    #
    n_layers_list = [   12]
    thickness_list = [ 0.06]
#    n_layers_list  = [   4,    8,    12,   16,   20]
#    thickness_list = [0.06, 0.06,  0.06, 0.06, 0.06]
    zip_list = zip(n_layers_list, thickness_list)
    
#    
    # constant 'rho_tex':
    #
#    n_layers_list  = [   6,    8,    12,   18,   24]
#    thickness_list = [0.03, 0.04,  0.06, 0.09, 0.12]
#    zip_list = zip( n_layers_list, thickness_list  )
#    legend_list = n_layers_list 

    # varying layers for two different heights':
    # use plot with M,N instead of nu,mu !
    #
#    n_layers_list = np.array([8,12,18])
#    thickness_list = np.array ([0.04,0.06])
#    zip_list =[]
#    for n in n_layers_list:
#        for t in thickness_list:
#            zip_list.append((n,t)) 


    for n, t in zip_list:
###################################################################################  
        #decide which cbls shall be plottet
       #  possible options: 'linear','cubic','fbm','bilinear'
        
      
        for  calib_config in [ 'fbm']:   
           
            sig_fl_calib.calib_config = calib_config
            sig_fl_calib.calib_sig_t_mfn()
            u_sol = sig_fl_calib.u_sol
            max_sig = sig_fl_calib.get_sig_max( u_sol )     
##            
####################################################################################            
#        #
#        # decide which interaction diagram depending on cclaw shall be plottet
#        # possible options: 'bilinear','block','quadratic'
#        
#                
#        for  sig_c_config in ['bilinear', 'block', 'quadratic']:   
#            sig_fl_calib.sig_c_config = sig_c_config
#            sig_fl_calib.calib_config = 'fbm'
#            sig_fl_calib.calib_sig_t_mfn()
#            u_sol = sig_fl_calib.u_sol
          
#        
###################################################################################   
#       # compare effects in the interaction diagramm for bilinear cclaw
#
#        for   xd  in [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.]:     
#            
#            sig_fl_calib.xd = xd
#            sig_fl_calib.sig_c_config ='quadratic'
#            sig_fl_calib.calib_config = 'bilinear'
#            sig_fl_calib.calib_sig_t_mfn()
#            u_sol = sig_fl_calib.u_sol
#            sig_fl_calib
###           
#            
            
                         
        #    sig_fl_calib.plot_sig_t_mfn( u_sol )
            
            
            #------------------------------------------------
            # 2) EVALUATION / VALIDATION:
            # get 'eps_lo', 'esp_up' for given/calibrated cb-law 
            #------------------------------------------------
        
            print '2 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            
            # reproduce the forces for the calibration test:
            #
            eps_lo = sig_fl_calib.eps_t_fail
            
            
        
            
            #------------------------------------------------
            # 3 construct failure envelope for a given value of eps_n 
            #------------------------------------------------
            print '3 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            
            n_A = 100
            n_B = 100
            n_C = 100
            
            ##########################################################################
            #mainly compression strain:
            
            if  sig_fl_calib.sig_c_config == 'quadratic_2':    
                eps_compression = -(2. + 0.085 * ((f_ck + 8.) - 50) ** 0.53) / 1000.  #eps_c2
                
            elif sig_fl_calib.sig_c_config == 'bilinear':   
            
                eps_compression = -(1.75 + 0.55 * (((f_ck + 8) - 50.) / 40.)) / 1000. #eps_c3
                
            elif sig_fl_calib.sig_c_config == 'quadratic':   
                
                eps_compression = -min(0.7 * (f_ck + 8.) ** 0.31, 2.8) / 1000. #eps_c1
                
            elif sig_fl_calib.sig_c_config == 'block':   
                
                eps_compression = -0.002
            
            ##########################################################################
       
          
            
            eps_c_fail = 3.3 / 1000.
            eps_t_fail = sig_fl_calib.eps_t_fail
            print'eps_t_fail ', eps_t_fail
            # mainly tension
            
            eps_lo_arr_A1 = np.linspace(eps_t_fail, eps_t_fail, n_A)
            eps_up_arr_A1 = np.linspace(eps_t_fail, 0, n_A)
            
            #  flexion2
            
            eps_lo_arr_A2 = np.linspace(eps_t_fail, eps_t_fail, n_A)
            eps_up_arr_A2 = np.linspace(0, -eps_c_fail, n_A)
            
            #  flexion3
            
            eps_lo_arr_B = np.linspace(eps_t_fail, 0., n_B)
            eps_up_arr_B = np.linspace(-eps_c_fail, -eps_c_fail, n_B)
            
            # pure compression
            
            eps_lo_arr_C = np.linspace(0., eps_compression, n_C)
            eps_up_arr_C = np.linspace(-eps_c_fail, eps_compression, n_C)
            
            # plot stress cases
            
            eps_lo_arr_psc = np.array  ([ 0., 0., 0., eps_t_fail, 0. , eps_t_fail, 0.])
            eps_up_arr_psc = np.array  ([ 0., -eps_c_fail, 0., -eps_c_fail, 0. , 0., 0.])
        
        #    all stress cases 

            eps_lo_arr = np.hstack([ eps_lo_arr_A1, eps_lo_arr_A2, eps_lo_arr_B, eps_lo_arr_C ])
            eps_up_arr = np.hstack([ eps_up_arr_A1, eps_up_arr_A2, eps_up_arr_B, eps_up_arr_C ])
##            
            psc = 'True'
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
        
        #         
            sig_fl_calib.n_layers = n        
            sig_fl_calib.thickness = t
            #
            eval_N_M_vectorized = np.frompyfunc(sig_fl_calib.eval_N_M, 2, 2)
            N_arr, M_arr = eval_N_M_vectorized(eps_lo_arr, eps_up_arr)
            
            N_arr_psc, M_arr_psc = eval_N_M_vectorized(eps_lo_arr_psc, eps_up_arr_psc)
            
            #
            # @todo: check which value to be used for normalization if sig_c = parabolic is used
            # NOTE that for this case f_cm  is used not f_ck; NOTE that for pure compression eps = 2E-3 is used
            # which does not yield f_cm in the parabolic law; for the default concrete laws for pure
            # compression the value should evaluated to 1.0!
            #
            c1 = 1000. * t * b * f_ck
            c2 = 1000. * t ** 2. * b * f_ck 
            #
            # @todo: why does an error occure when multiplying an array by a float value (caused by vectorized)
            nu_arr = N_arr / np.array([ c1 ]) 
            mu_arr = M_arr / np.array([ c2 ])
            
            nu_arr_psc = N_arr_psc / np.array([ c1 ]) 
            mu_arr_psc = M_arr_psc / np.array([ c2 ])
         
            rho_tex = n * A_roving * n_rovings / c1
            #
            N_list.append(N_arr)
            M_list.append(M_arr)
            nu_list.append(nu_arr)
            mu_list.append(mu_arr)
            rho_list.append(rho_tex)
            
            #############################################
            #plot nu-mu-interaction
            
            if psc == 'True' and sig_fl_calib.calib_config == 'fbm':
                    p.plot (mu_arr_psc[0, :], nu_arr_psc[0, :], 'k--',
                                                    color= 'turquoise',
                                                  linewidth = 2.0)
            if psc == 'True' and sig_fl_calib.calib_config == 'cubic':
                    p.plot (mu_arr_psc[0, :], nu_arr_psc[0, :], 'k--',
                                                    color= 'red',
                                                  linewidth = 2.0)       
            if psc == 'True' and sig_fl_calib.calib_config == 'linear':
                    p.plot (mu_arr_psc[0, :], nu_arr_psc[0, :], 'k--',
                                                    color= 'green',
                                                  linewidth = 2.0)
            if psc == 'True' and sig_fl_calib.calib_config == 'bilinear':
                    p.plot (mu_arr_psc[0, :], nu_arr_psc[0, :], 'k--',
                                                    color= 'purple',
                                                  linewidth = 2.0) 

########################### distinction of stress cases all black

#            if psc == 'True':
#                    p.plot (mu_arr_psc[0, :], nu_arr_psc[0, :], 'k--',
#                                                    color= 'black',
#                                                  linewidth = 1.0)            
############################
    
            if t == 0.06 and n == 12:
                
                if sig_fl_calib.calib_config == 'linear':
                    p.plot(mu_arr, nu_arr,label='linear',
                            color='green', 
                            linewidth = 2.0)
                elif sig_fl_calib.calib_config == 'cubic':  
                    p.plot(mu_arr, nu_arr,label='cubic',
                            color='red', 
                            linewidth = 2.0)
                    
                elif sig_fl_calib.calib_config == 'bilinear':  
                    p.plot(mu_arr, nu_arr,label='bilinear',
                            color = 'purple',
                            linewidth = 2.0)
                    
                elif sig_fl_calib.calib_config == 'fbm':  
                    p.plot(mu_arr, nu_arr,label='fbm',
                            color = 'turquoise',
                            linewidth = 2.0)
                         
            else:
                
                if sig_fl_calib.calib_config == 'linear':
                    p.plot(mu_arr, nu_arr,label='linear',
                            #color='green', 
                            linewidth = 1.0)
                elif sig_fl_calib.calib_config == 'cubic':  
                    p.plot(mu_arr, nu_arr,label='cubic',
                            color = 'red',
                            linewidth = 1.0)
                    
                elif sig_fl_calib.calib_config == 'bilinear':  
                    p.plot(mu_arr, nu_arr,label='bilinear',
                            color = 'purple',
                            linewidth = 1.0)
                    
                elif sig_fl_calib.calib_config == 'fbm':  
                    p.plot(mu_arr, nu_arr,label='fbm',
                            color = 'turquoise',
                            linewidth = 1.0)
#                    
                        
   
                    
#    #    
#            # plot absolute values for N and M
#        #        
#            if psc == 'True':
#                        p.plot (M_arr_psc[0, :], N_arr_psc[0, :], 'k--',
#                                                          color = 'black',
#                                                          linewidth = 1.0)
#
#            if t == 0.06 and n == 12 :
##                    
#                    if sig_fl_calib.calib_config =='linear':
#                        p.plot( M_arr, N_arr,
#                               # label = 'linear' ,
#                                color='green', 
#                                linewidth=2.0)
#                    elif sig_fl_calib.calib_config == 'cubic':  
#                        p.plot( M_arr, N_arr,
#                                label= 'cubic',
#                                color='red',  
#                                linewidth=2.0)
#                        
#                    elif sig_fl_calib.calib_config == 'bilinear' and sig_fl_calib.xd==0.:  
#                        p.plot( M_arr, N_arr, 
#                                color='purple',
#                               # label= 'bilinear', 
#                                linewidth=2.0)    
#                    elif sig_fl_calib.calib_config == 'bilinear' and sig_fl_calib.xd==1.0:  
#                        p.plot( M_arr, N_arr, 
#                                color='green',
#                               # label= 'bilinear', 
#                                linewidth=2.0)        
#                    elif sig_fl_calib.calib_config == 'bilinear' and sig_fl_calib.xd==0.1:  
#                        p.plot( M_arr, N_arr, 
#                                color='yellow',
#                               # label= 'bilinear', 
#                                linewidth=2.0)       
#                        
#                    elif sig_fl_calib.calib_config == 'bilinear':  
#                        p.plot( M_arr, N_arr, 
#                               # color='purple',
#                               # label= 'bilinear', 
#                                linewidth=2.0)
#                        
#                    elif sig_fl_calib.calib_config == 'fbm':  
#                        p.plot( M_arr, N_arr,
#                                label= 'fbm',
#                                color='turquoise', 
#                                linewidth=2.0)
#             
#                           
#                           
#            else:
#                    if sig_fl_calib.calib_config =='linear':
#                        p.plot( M_arr, N_arr,
#                                color='green', 
#                                linewidth=1.0)
#                    elif sig_fl_calib.calib_config == 'cubic':  
#                        p.plot( M_arr, N_arr, 
#                                color='red', 
#                                linewidth=1.0)
#                        
#                    elif sig_fl_calib.calib_config == 'bilinear':  
#                        p.plot( M_arr, N_arr, 
#                                #color='purple', 
#                                linewidth=1.0)
#                        
#                    elif sig_fl_calib.calib_config == 'fbm':  
#                        p.plot( M_arr, N_arr, 
#                                color='turquoise', 
#                                linewidth=1.0)
#        
#                  
        
        
        ######### @todo: format plot and axes (legend, grid, font etc.) #########
    
        ### GRID
        #
        p.grid(b = None, which = 'major')
    
    
        
###################################################################################
     ### LEGEND
     
#        legend_list = [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.] 
#        legend_list = thickness_list 
#        legend_list = ['linear','cubic','fbm','bilinear']
#        legend_list = ['bilinear', 'block', 'quadratic']
#        legend_list = n_layers_list 
###################################################################################
#        p.legend( legend_list, 'upper right', shadow = True)
        p.legend()
        
    
    
        ### TICKS
        # nu-mu-interaction
        p.axis([0, 0.15, 0.5 , -1.25])
        ax.set_xticks([0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09,0.105,0.12,0.135,0.15])
        ax.set_yticks([0.4,0.2,0.,-0.2,-0.4,-0.6,-0.8,-1,-1.2])
        
        # n-m-interaction
        
#        p.axis([0, 6.4, 250 ,-900])
#        p.xticks([0.,0.8,1.6,2.4,3.2,4.,4.8,5.6,6.4])
#        p.yticks([300.,150.,0.,-150.,-300.,-450.,-600.,-750.])
 

    
# #      for automatic ticks
#        ax.xaxis.set_major_locator(MaxNLocator(10))
#        ax.yaxis.set_major_locator(MaxNLocator(10))
#        minorLocator = AutoMinorLocator()
#        ax.xaxis.set_minor_locator(minorLocator)

#        ax.tick_params(axis = 'both', which = 'major', direction = 'out', length = 6, width = 2, colors = 'black')
#        ax.tick_params(axis = 'both', which = 'minor', direction = 'out', length = 0, width = 0, colors = 'black')
    
    
        ### TITEL
        #
    #    p.title(r'$\nu$ -$\mu$ Interaction Diagram ',fontsize=22,
    #                                verticalalignment='baseline',
    #                                horizontalalignment = 'center')

        
        ### LABELS
        # M-N
        
#        p.xlabel('M',fontstyle='italic', fontsize = '15',
#                 verticalalignment = 'top',
#                 horizontalalignment = 'left')
#        p.ylabel('N ',fontstyle='italic', fontsize = '15')
        
        #mu-nu
        
        p.xlabel(r'$\mu$', fontsize = '20',
                 verticalalignment = 'top',
                 horizontalalignment = 'left')
        p.ylabel(r'$\nu$ ', fontsize = 20)
        p.show()
