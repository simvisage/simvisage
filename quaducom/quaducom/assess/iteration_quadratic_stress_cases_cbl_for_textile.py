'''
Created on Jun 23, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    Float, Instance, Array, Int, Property, cached_property, on_trait_change, Bool, \
    HasTraits, File, Event, Trait

from etsproxy.traits.ui.api import \
    View, Item, FileEditor, HSplit, Group, VSplit, \
    Handler

from etsproxy.traits.ui.menu import \
    Action, CloseAction, HelpAction, Menu, \
    MenuBar, NoButtons, Separator, ToolBar

import numpy as np

from mathkit.mfn import MFnLineArray

from scipy.optimize import fsolve

data_file_editor = FileEditor(filter = ['*.DAT'])

from matresdev.db.simdb import SimDB
simdb = SimDB()

class SigFlCalib(HasTraits):

    #---------------------------------------------------------------
    # material properties textile reinforcement
    #-------------------------------------------------------------------

    # security factor 'gamma_tex' for the textile reinforcement
    #
#    gamma_tex = Float(1.5, input = True)

    # reduction factor to drive the characteristic value from mean value
    # of the experiment (EN DIN 1990)
    #
#    beta = Float(0.81, input = True)

    #---------------------------------------------------------------
    # material properties concrete
    #-------------------------------------------------------------------

    # characteristic compressive force [MPa]
    #
    f_ck = Float(60.0, input = True)

#    # security factor 'gamma_c' for hight strength concrete (>C55/67)
#    #
#    gamma_c = Property
#    def _get_gamma_c( self ):
#        return 1.5 * ( 1 / ( 1.1 - self.f_ck / 500. ) )
#
#    # design value of the compressive force [MPa]
#    #
#    f_cd = Property
#    def _get_f_cd( self ):
#        return 0.85 * self.f_ck / self.gamma_c

    # compressive strain at the top at rupture [-]
    # (positive value is used)
    # value measured at the experiment (with DMS)
    # for 0: about 3.3
    # for 90: about 3.0
    # ultimute strain theoreticaly (Brockman): about 4.5
    # NOTE: strain was meassured at a distance of 5 cm
    #
    eps_c = Float(0.003, input = True) # float value corresponds to 3 promile

    #---------------------------------------------------------------
    # properties of the composite cross section
    #-------------------------------------------------------------------

    # thickness of reinforced cross section
    #
    thickness = Float(0.06, geo_input = True)

    # total number of reinforcement layers [-]
    # 
    n_layers = Float(12, geo_input = True)

    # spacing between the layers [m]
    #
    s_tex_z = Property(depends_on = '+geo_input')
    @cached_property
    def _get_s_tex_z(self):
        return self.thickness / (self.n_layers + 1)

    # distance from the top of each reinforcement layer [m]:
    #
    z_t_i_arr = Property(depends_on = '+geo_input')
    @cached_property
    def _get_z_t_i_arr(self):
        return np.array([ self.thickness - (i + 1) * self.s_tex_z for i in range(self.n_layers) ],
                      dtype = float)

    #---------------------------------------------------------------
    # properties of the bending specimens (tree point bending test):
    #---------------------------------------------------------------

    # width of the cross section [m]
    #
    width = Float(0.20, input = True)
    
    # number of rovings in 0-direction of one composite 
    # layer of the bending test [-]:
    #
    n_rovings = Int(23, input = True)

    # cross section of one roving [mm**2]:
    #
    A_roving = Float(0.461, input = True)


    # ultimate textile stress meassured in the tensile test [MPa]
    #
    sig_tex_fail = Float( 1216. )
    
    def calib_layer_response( self, u ):
        
        '''CALIBRATION: derive the unknown constitutive law of the layer
        (effective crack bridge law)
        '''
        eps_t, var_a = u

        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_t'
        # ------------------------------------------------------------------------                
        thickness = self.thickness
        z_t_i_arr = self.z_t_i_arr
        
        # ------------------------------------------------------------------------                
        # derived params depending on value for 'eps_t'
        # ------------------------------------------------------------------------                

        # heights of the compressive zone:
        #
        x = abs(self.eps_c) / (abs(self.eps_c) + abs(eps_t)) * thickness

        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = eps_t / (thickness - x) * (z_t_i_arr - x)

        # use a ramp function to consider only negative strains
        eps_c_i_arr = (-np.fabs(eps_i_arr) + eps_i_arr) / 2.0

        # use a ramp function to consider only positive strains
        eps_t_i_arr = (np.fabs(eps_i_arr) + eps_i_arr) / 2.0

        # ------------------------------------------------------------------------                
        # function used to approximate the effective crack bridge law
        # ------------------------------------------------------------------------                
        
        # linear cb-law using only the effective yarn modulus
        if self.calib_config == 'calibration-linear':
            E_yarn = var_a
            xdata = np.array([ 0., 1. ])
            ydata = np.array([ 0., E_yarn ])

        # linear cb-law up to reaching the rupture strain, behaves plastic afterwards
        elif self.calib_config == 'calibration-bilinear':
            E_yarn = var_a
            sig_fail = self.sig_tex_fail
            eps_fail = sig_fail / E_yarn
            xdata = np.array([0., eps_fail, 2.*eps_fail ])
            ydata = np.array([0., sig_fail, sig_fail ])

        # full plastic cb-law 
        elif self.calib_config == 'calibration-plastic':
            sig_fail = self.sig_tex_fail
            eps_fail = eps_t_i_arr[0]
            xdata = np.array([0., 0.01*eps_fail, eps_fail ])
            ydata = np.array([0., 0.99*sig_fail, sig_fail ])

        # quadratic cb-law up to reaching the rupture strain, behaves plastic afterwards
        elif self.calib_config == 'calibration-quadratic':
            # use strain at the lowest textile layer as rupture strain 
            eps_fail = eps_t_i_arr[0]
            sig_fail = self.sig_tex_fail
            eps_arr = np.arange(0, eps_fail, eps_fail / 50.) 
            var_b = (sig_fail - var_a * eps_fail **2) / eps_fail 
            sig_tex_arr = var_a * eps_arr ** 2 + var_b * eps_arr 
            xdata = np.hstack([eps_arr, 2.*eps_arr[-1]])  
            ydata = np.hstack([sig_tex_arr, sig_tex_arr[-1]]) 
            print 'var_b', var_b 

        sig_t_mfn = MFnLineArray( xdata = xdata, ydata = ydata )
        
        return x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn


    def eval_layer_response_f( self, u ):
        '''EVALUATION: using the calibrated constitutive law of the layer
        '''
        eps_t, eps_c = u
#        eps_t = eps_lo
#        eps_c = eps_up
        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_t'
        # ------------------------------------------------------------------------                

        thickness = self.thickness
        z_t_i_arr = self.z_t_i_arr
       
        # ------------------------------------------------------------------------                
        # derived params depending on value for 'eps_t'
        # ------------------------------------------------------------------------                

        # heights of the compressive zone:
        #
        x = abs(self.eps_c) / (abs(self.eps_c) + abs(eps_t)) * thickness

        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = eps_t / (thickness - x) * (z_t_i_arr - x)

        # use a ramp function to consider only negative strains
        eps_c_i_arr = (-np.fabs(eps_i_arr) + eps_i_arr) / 2.0

        # use a ramp function to consider only positive strains
        eps_t_i_arr = (np.fabs(eps_i_arr) + eps_i_arr) / 2.0
        
        sig_t_mfn = self.sig_t_mfn

        return x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn



    def eval_layer_response_t( self, u ):
        '''EVALUATION: using the calibrated constitutive law of the layer
        ''' 
        eps_lo, eps_up = u

        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_tl'
        # ------------------------------------------------------------------------                

        thickness = self.thickness
        z_t_i_arr = self.z_t_i_arr
   
        # ------------------------------------------------------------------------                
        # derived params depending on value for 'eps_tl'
        # ------------------------------------------------------------------------                

        # heights of the compressive zone: the whole cross-section is stressed with tension            #
        x = 0

        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = eps_lo - (eps_lo - eps_up ) / thickness * (thickness - z_t_i_arr)
        eps_t_i_arr = eps_i_arr
        eps_c_i_arr = (np.fabs(eps_i_arr) + eps_i_arr) / 2.0
        sig_t_mfn = self.sig_t_mfn

        return x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn
     
    def eval_layer_response_c( self, u ):
        '''EVALUATION: using the calibrated constitutive law of the layer
        ''' 
        eps_lo, eps_up = u

        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_cl'
        # ------------------------------------------------------------------------                

        thickness = self.thickness
        z_t_i_arr = self.z_t_i_arr
   
        # ------------------------------------------------------------------------                
        # derived params depending on value for 'eps_cl'
        # ------------------------------------------------------------------------                

        # heights of the compressive zone: the whole cross-section is stressed by compression            #
        x = thickness

        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = eps_lo - (eps_lo - eps_up) / thickness * z_t_i_arr
        eps_t_i_arr = (-np.fabs(eps_i_arr) + eps_i_arr) / 2.0
        eps_c_i_arr = eps_i_arr
        sig_t_mfn = self.sig_t_mfn

        return x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn



    def get_f_i_arr( self, u ):
        '''tensile force at the height of each reinforcement layer [kN]:
        '''
        x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn = self.layer_response( u )

        get_sig_t_i_arr = np.frompyfunc( sig_t_mfn.get_value, 1, 1 )
        sig_t_i_arr = get_sig_t_i_arr( eps_t_i_arr )

        # print 'sig_i_arr', sig_i_arr

        # tensile force of one reinforced composite layer [kN]:
        #
        n_rovings = self.n_rovings
        A_roving = self.A_roving
        f_t_i_arr = sig_t_i_arr * n_rovings * A_roving / 1000.
        # print 'f_t_i_arr', f_t_i_arr

        # compressive strain in the compression zone 'x' in each sub-area ('A_c' divided in 'n_c' parts)
        #
        sig_c_i_arr = np.array([self.get_sig_c(eps_c_i) for eps_c_i in eps_c_i_arr])
        #visualize signed sig_c_i_arr 
        sig_c_i_arr = -sig_c_i_arr
        f_c_i_arr = sig_c_i_arr * self.width * self.s_tex_z * 1000.

        return x, f_t_i_arr, f_c_i_arr


    def get_sig_comp_i_arr( self, u ):
        '''tensile stress at the height of each reinforcement layer [MPa]:
        '''
        x, f_t_i_arr, f_c_i_arr = self.get_f_i_arr( u )
        return (f_t_i_arr - f_c_i_arr) / self.width / self.s_tex_z / 1000.0


    #-----------------------------
    # for simplified constant stress-strain-diagram of the concrete (EC2)
    #-----------------------------
    sig_c_mfn_block = Property
    def _get_sig_c_mfn_block(self):
        #(for standard concrete)
        if self.f_ck <= 50:
            lamda = 0.8
            eta = 1.0  
            eps_cu3 = 0.0035
        # (for high strength concrete)
        #
        else:
            eta = 1.0 - ( self.f_ck / 50.) / 200.
        # factor [-] to calculate the height of the compressive zone  
            lamda = 0.8 - (self.f_ck - 50.) / 400.
            eps_cu3 = 2.6 + 35. * (90. - self.f_ck) **4 / 100000000 
#        xdata = np.array([0., (1. - lamda) * eps_cu3 - 0.0000001, (1 - lamda) * 0.0035, 0.0035,0.0035+0.0000001,0.0035*2]) 
#        ydata = np.array([0., 0., eta*self.f_ck, eta * self.f_ck,0.,0.])
        #concrete law without limit for eps_c
        xdata = np.array([0., (1. - lamda) * eps_cu3 - 0.0000001, (1 - lamda) * 0.0035, 0.0035]) 
        ydata = np.array([0., 0., eta*self.f_ck, eta * self.f_ck])
        return MFnLineArray(xdata = xdata, ydata = ydata)   


    #-----------------------------
    # for bilinear stress-strain-diagram of the concrete (EC2)
    #-----------------------------
    sig_c_mfn_bilinear = Property
    def _get_sig_c_mfn_bilinear(self):
        #(for standard concrete)
        if self.f_ck <= 50.:
            epsilon_c3 = 0.00175
            epsilon_cu3 = 0.0035
        #(for high strength concrete)
        else :
            epsilon_c3 = 1,75 + 0,55 * (self.f_ck - 50.) / 40.
            epsilon_cu3 = 2,6 + 35 * (90 - self.f_ck) ** 4. / 100000000.      
  
#        xdata = np.array( [0., epsilon_c3, epsilon_cu3,epsilon_cu3+0.00000001,epsilon_cu3*2] )
#        ydata = np.array( [0., self.f_ck, self.f_ck,0.,0.] )
        
        #concrete law without limit for eps_c
        xdata = np.array( [0., epsilon_c3, epsilon_cu3] )
        ydata = np.array( [0., self.f_ck, self.f_ck] )
        return MFnLineArray(xdata = xdata, ydata = ydata)


    #-----------------------------
    # for quadratic stress-strain-diagram of the concrete
    #-----------------------------
    sig_c_mfn_quadratic = Property
    def _get_sig_c_mfn_quadratic( self ):
        # (for all concretes up to f_cm=88 N/mm2) #max epislon_c1u
        f_cm = self.f_ck + 8
        E_tan = 9500 * (f_cm) ** (1/3) #SBT 17
        E_sec = f_cm / 0.0022
        epsilon_c1 = 0.0022 # min(0.7*f_cm**0.31, 2.8)/1000 #EC
        epsilon = np.arange(0, 0.0022, 0.000022)
        sigma_c = (E_tan / E_sec * epsilon / epsilon_c1 - (epsilon / epsilon_c1) **2) / (1 + (E_tan/E_sec - 2)*epsilon / epsilon_c1) * f_cm
#        print'epsilon', epsilon
#        print'sigma_c=', sigma_c
#        xdata = np.hstack([ epsilon, epsilon[-1] + 0.00000001, 2.*epsilon[-1]])  
#        ydata = np.hstack([sigma_c, 0, 0 ])
        #concrete law without limit for eps_c
        xdata = epsilon  
        ydata = sigma_c
      
#        xdata = epsilon
#        ydata = sigma_c
        return MFnLineArray( xdata = xdata, ydata = ydata )

    sig_c_config = Trait('quadratic',{
                          'quadratic' : 'sig_c_mfn_quadratic',
                          'bilinear'  : 'sig_c_mfn_bilinear',
                          'block'     : 'sig_c_mfn_block'
                          },
                         modified = True)

    sig_c_mfn = Property(depends_on = 'sig_c_config')
    def _get_sig_c_mfn(self):
        return getattr( self, self.sig_c_config_)

    def get_sig_c(self, eps_c):
        sig_c = self.sig_c_mfn.get_value(eps_c)
        return sig_c

    sig_c_mfn_vect = Property(depends_on = '+input')
    @cached_property
    def _get_sig_c_mfn_vect(self):
        return np.frompyfunc( self.sig_c_mfn.get_value, 1, 1 )

    get_sig_c_vectorized = np.frompyfunc( get_sig_c, 1, 1 )

    # number of subdivisions of the compressive zone
    #
    n_c = float(5)

    x_frac_i = Property(Array)
    @cached_property
    def _get_x_frac_i(self):
        '''subdivide the compression zone 'x' in 'n_c' sub-areas;
        'x_frac_i' giving the fraction of each distance of the sub-area from the upper surface 
        with respect to the compressive zone 'x'
        '''
#        print 'x_frac_i', np.arange( self.n_c ) / self.n_c + 1. / ( 2. * self.n_c )
        return np.arange(self.n_c) / self.n_c + 1. / (2. * self.n_c)

    # iteration counter
    #
    n = 0

    def get_lack_of_fit(self, u, M, N):
        '''Return the difference between 
        N_c (=compressive force of the compressive zone of the concrete) 
        N_t (=total tensile force of the reinforcement layers)

        NOTE: eps_t (=tensile strain at the bottom [MPa]) is the search parameter
        to be found iteratively!
        '''
        print '------------- iteration: %g ----------------------------' % ( self.n )

       
        
        if self.calc_mode == 'eval':
            stress_case = self.determine_stress_case( M, N )
#            print 'stress_case',stress_case
            self.eval_config = stress_case
      


        x, f_t_i_arr, f_c_i_arr = self.get_f_i_arr( u )

        # height of the cross section
        #
        thickness = self.thickness
#        print 'thickness', thickness

        z_t_i_arr = self.z_t_i_arr
#        print 'z_t_i_arr', z_t_i_arr

        # total tensile force of all layers of the composite tensile zone [kN]:
        # (characteristic value)
        #
        N_tk = sum(f_t_i_arr)
        # print 'N_tk', N_tk

        # total tensile force of all layers of the composite tensile zone [kN]:
        # (design value)
        #
#        N_td = N_tk / self.gamma_tex * self.beta
#        print 'N_td', N_td

#        k_exact = ( 1.74 * self.eps_c / 4.56 - ( self.eps_c / 4.56 ) ** 2 / ( 1 - 0.12 * self.eps_c / 4.56 ) )

#        # distance [m] of the resulting compressive 
#        # force from the top, i.e. a = k * x / 2
#        #
#        a = self.k * x / 2.
#        # print 'a', a
#
#        # total compressive force of the composite compressive zone [kN]:
#        # (characteristic value)
#        #
#        N_ck = 2.0 * a * self.width * self.chi * self.f_ck * 1000.

        # compressive strain in the compression zone 'x' in each sub-area ('A_c' divided in 'n_c' parts)
        #
        if self.calib_config == 'evaluation-flexion' or \
                                    'evaluation-tension':                            
             
                eps_c_i_arr = (1-self.x_frac_i) * self.eps_c
                
        elif  self.calib_config == 'evaluation-compression':
            eps_lo = self.eps_lo
            eps_up = self.eps_up
            
            if eps_lo < eps_up:
                 
                
                eps_c_i_arr = eps_lo + (1-self.x_frac_i) * (eps_up-self.eps_lo)
                
            else:
                
                eps_c_i_arr = eps_up + (self.x_frac_i) * (eps_lo-self.eps_up)
                
                
#        print 'eps_c_i_arr', eps_c_i_arr

        # @todo: get vectorized version running: 
#        sig_c_i_arr = self.sig_c_mfn_vect( eps_c_i_arr )
#        get_sig_c_i_arr = np.frompyfunc( self.sig_c_mfn.get_value, 1, 1 )
#        sig_c_i_arr = get_sig_c_i_arr( eps_c_i_arr )
        sig_c_i_arr = np.array([self.get_sig_c(eps_c_i) for eps_c_i in eps_c_i_arr])
        
#        print 'sig_c_i_arr', sig_c_i_arr
        f_c_i_arr = sig_c_i_arr * self.width * self.x_frac_i[0] * 2. * x * 1000.

        z_c_i_arr = x * self.x_frac_i
        

        N_ck = sum(f_c_i_arr)

#        # total compressive force of the composite compressive zone [kN]:
#        # (design value)
#        #
#        N_cd = x * self.width * self.chi * self.f_cd * 1000.
#        print 'N_cd', N_cd

        # absolute error (equilibrium of sum N = 0):
        #
        N_external = N
        N_internal = -N_ck + N_tk
        d_N = N_internal - N_external

#        print 'd_N', d_N

        # resistance moment of one reinforced composite layer [kNm]:
        # (characteristic value)
        #
        M_tk = np.dot(f_t_i_arr, z_t_i_arr)
        # print 'M_tk', M_tk

#        M_ck = -a * N_ck
        M_ck = np.dot(f_c_i_arr, z_c_i_arr)

        M_internal = M_tk + M_ck

        M_external = M + N_external * thickness / 2.0

        d_M = M_internal - M_external
#        print 'd_M', d_M

        self.n += 1
#        print 'n', self.n

        return np.array([ d_N, d_M ], dtype = float)

    
    # iteration counter
    #
    m = 0

    def fit_response( self, M, N ):
#                      elem_no = 0, mx = 0.0, my = 0.0, mxy = 0.0, nx = 0.0, ny = 0.0, nxy = 0.0, \
#                      sig1_up = 0, sig1_lo_sig_up = 0, sig1_lo = 0, sig1_up_sig_lo = 0, ):
        '''iterate 'eps_t' such that the lack of fit between the calculated
        normal forces in the tensile reinforcement and the compressive zone (concrete)
        is smaller then 'xtol' defined in function 'brentq'.
        NOTE: the method 'get_lack_of_fit' returns the relative error.
        '''
        self.m += 1
#        print '--- fit_response called --- %g' % ( self.m )

#      

        # use scipy-functionality to get the iterated value of 'eps_t'
            # NOTE: get_lack_of_fit must have a sign change as a requirement
            # for the function call 'brentq' to work property. 

            # The method brentq has optional arguments such as
            #   'xtol'    - absolut error (default value = 1.0e-12)
            #   'rtol'    - relative error (not supported at the time)
            #   'maxiter' - maximum numbers of iterations used
            #
        xtol = 1.0e-5

            #----------------
            # @todo: how can the rupture strain in the bending test be estimated realistically?
            # in which boundaries shall brentq search for a solution? (5*eps_u)?
            #----------------
#            E_comp = 30000. # MPa
#            eps_t_estimate = abs(sig_plus) / E_comp
#            eps_c_estimate = abs(sig_minus) / E_comp
#            u0 = [eps_t_estimate, eps_c_estimate]
        
        u0 = self.u0
        u_sol = fsolve(self.get_lack_of_fit, u0, args = (M, N), xtol = xtol)
#            print 'u_sol', u_sol

        
            # @todo: check if 'brenth' gives better fitting results; faster? 
    #            phi_new = brenth( self.get_lack_of_fit, 0., eps_t )
            
#        print 'u_sol', u_sol
#        print 'u_sol.shape', u_sol.shape
#        print 'type(u_sol)', type( u_sol )
#        return u_sol[0], u_sol[1]
        return u_sol
  
    def determine_stress_case(self, M, N):
#         print self.thickness
#        print self.width

        thickness = self.thickness

        W = thickness ** 2. * self.width / 6.0
        A = thickness * self.width

        sig_bending = M / W / 1000.0
        sig_normal = N / A / 1000.0

        sig_plus = sig_normal + sig_bending
        sig_minus = sig_normal - sig_bending
#        print 'M', M
#        print 'N', N
#        print 'W', W
#        print 'A', A
#        print 'M/W', sig_bending
#        print 'N/A', sig_normal
 
        print 'sig_plus', sig_plus
        print 'sig_minus', sig_minus
        
        if sig_plus * sig_minus < 0.0:
            stress_case = 'evaluation-flexion'

        elif sig_plus < 0.0 and sig_minus < 0.0:
            stress_case = 'evaluation-compression'
          
        elif sig_plus > 0.0 and sig_minus > 0.0:
            stress_case = 'evaluation-tension'
        
        print 'evaluation for stress case: ', stress_case
        return stress_case
    # set configuration for calibration or evaluation    
    
   
     
    eval_config = Trait('evaluation-flexion',
                          {
                           'evaluation-flexion' : ('eval_layer_response_f',
                                              np.array([ 0.01, 0.0033 ])),
                           'evaluation-tension' : ('eval_layer_response_t',
                                              np.array([ 0.02, 0.02 ])),
                           'evaluation-compression' : ('eval_layer_response_c',
                                             np.array([ 0.0033, 0.0033 ])) },
                         modified = True)
       
       
    calib_config = Trait('calibration-quadratic',
                          {'calibration-quadratic' : ('calib_layer_response',
                                              np.array([ 0.01, -50000 ])),
                           'calibration-linear' : ('calib_layer_response',
                                              np.array([ 0.01, 50000.0 ])),
                           'calibration-plastic' : ('calib_layer_response',
                                              np.array([ 0.01, 50000.0 ])),
                           'calibration-bilinear' : ('calib_layer_response',
                                              np.array([ 0.01, 50000.0 ])) },
                         modified = True)

    calc_mode = Trait('calib',
                    {'calib' : ('calib_config'),
                     'eval'  : ('eval_config')},
                                modified = True)
                              
    layer_response = Property(depends_on = 'calib_config')
    @cached_property
    def _get_layer_response(self):
        return getattr(self, self.calib_config_[0])

    u0 = Property(Array(float), depends_on = 'calib_config')
    @cached_property
    def _get_u0(self):
        return self.calib_config_[1]

    def get_sig_max(self, u):
        sig_max = np.max(self.get_sig_comp_i_arr( u ))
        print 'sig_max', sig_max
        return sig_max

    sig_t_mfn = MFnLineArray()



if __name__ == '__main__':

    import pylab as p

    #------------------------------
    # define input params
    #------------------------------
    #
    N = 0

    # measured value in bending test
    #
    M = 3.49

    # value per m
#    M = 5*3.49

#    M = 2.40048145726
#    N = -46.4211016241

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
                               eps_c = 3.3 / 1000.,

                               # values for experiment beam with width = 0.20 m
                               #
                               width = 0.20,
                               n_roving = 23,

#                               # values per m
#                               #
#                               width = 1.0,
#                               n_roving = 120,

                               # define shape of the crack-bridge law ('linear', 'bilinear' or 'quadratic')
                               #
                               calc_mode = 'calib',
                               calib_config = 'calibration-quadratic',
                               eval_config = 'evaluation-flexion' ,

                               # define shape of the concrete stress-strain-law ('block', 'bilinear' or 'quadratic')
                               #
                               sig_c_config = 'quadratic'
                              
                              )

    eps_t, var_a = sig_fl_calib.fit_response( M, N )
    u_sol = [ eps_t, var_a ]
    eps_c = sig_fl_calib.eps_c
    print 'eps_c', eps_c
    print 'eps_t', eps_t
    print 'u_sol', u_sol
    print 'var_a', var_a
    x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn = sig_fl_calib.layer_response( u_sol )

    
    #------------------------------------------------
    # 2) EVALUATION / VALIDATION:
    # get 'eps_lo', 'esp_up' for given/calibrated cb-law 
    #------------------------------------------------
    
    sig_fl_calib.set( sig_t_mfn = sig_t_mfn, calc_mode = 'eval',
                      eval_config = 'evaluation-flexion' )
    
    # select the right stress_mode for calculation
    #
    stress_case = sig_fl_calib.determine_stress_case( M, N)
    eval_config = stress_case
    print 'stress_case', stress_case


    eps_lo, eps_up = sig_fl_calib.fit_response( M, N )
    print 'eps_lo', eps_lo
    print 'eps_up', eps_up
    print 'max_sig', sig_fl_calib.get_sig_max( [ eps_lo, eps_up ] )                      
                              

    #------------------------------------------------
    # plot response
    #------------------------------------------------
    #
    layer_arr = np.arange(sig_fl_calib.n_layers)
    sig_comp_i_arr = sig_fl_calib.get_sig_comp_i_arr( [eps_lo, eps_up] )
    p.bar(layer_arr, sig_comp_i_arr, 0.2)
    p.show()

