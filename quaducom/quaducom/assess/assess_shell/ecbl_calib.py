'''
Created on Jun 23, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    Float, Instance, Array, Int, Property, cached_property, on_trait_change, Bool, \
    HasTraits, File, Event, Trait, Str, List

from etsproxy.traits.ui.api import \
    View, Item, FileEditor, HSplit, Group, VSplit, \
    Handler

from etsproxy.traits.ui.menu import \
    Action, CloseAction, HelpAction, Menu, \
    MenuBar, NoButtons, Separator, ToolBar

import math

import numpy as np

import pylab as p

from mathkit.mfn import MFnLineArray

from scipy.optimize import fsolve

data_file_editor = FileEditor(filter = ['*.DAT'])

from matresdev.db.simdb import SimDB
simdb = SimDB()

from math import exp, log

class ECBFBase(HasTraits):
    '''Base class for Effective Crack Bridge Laws.'''
    
class ECBLFBM(ECBFBase):
    '''Base class for Effective Crack Bridge Laws.'''
    
    sig_u_tt = Float
        
    def __call__(self, eps_u_tt, m, eps):
        return (self.sig_u_tt / eps_u_tt / exp(-pow(exp(-log(m) / m), 1.0 * m)) * 
                eps * np.exp(-np.power(eps / eps_u_tt * exp(-log(m) / m), 1.0 * m)))
        

class ECBLCalib(HasTraits):

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
    # material properties 
    #-------------------------------------------------------------------

    # characteristic compressive stress [MPa]
    #
    f_ck = Float(60.0, sig_c_modified = True)

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

    # ultimate textile stress measured in the tensile test [MPa]
    #
    sig_tex_fail = Float(1216., sig_t_modified = True)
    
    
    # compressive strain at the top at rupture [-]
    # (positive value is used)
    # value measured at the experiment (with DMS)
    # for 0: about 3.3
    # for 90: about 3.0
    # ultimate strain theoretically (Brockman): about 4.5
    # NOTE: strain was meassured at a distance of 5 cm
    #
    eps_c_fail = Float(0.0033, sig_t_modified = True) # float value corresponds to 3 promile

    # ultimate limit strain at the tensile surface of the cross section [-]:
    # value is derived based on the calibrated crack bridge law
    #
    eps_t_fail = Float
    sig_t_mfn = Instance(MFnLineArray)

    # rupture moment and normal force meassured in the calibration experiment
    # (three point bending test)
    #
    M_fail = Float(3.5) # [kNm]
    N_fail = Float(0.0) # [kN]

    #---------------------------------------------------------------
    # properties of the composite cross section
    #-------------------------------------------------------------------

    # thickness of reinforced cross section
    #
    thickness = Float(0.06, geo_input = True)

    # total number of reinforcement layers [-]
    # 
    n_layers = Int(12, geo_input = True)

    # spacing between the layers [m]
    #
    s_tex_z = Property(depends_on = '+geo_input')
    @cached_property
    def _get_s_tex_z(self):
        return self.thickness / (self.n_layers + 1)
    print 's_tex_z', s_tex_z
    
    # distance from the top of each reinforcement layer [m]:
    #
    z_t_i_arr = Property(depends_on = '+geo_input')
    @cached_property
    def _get_z_t_i_arr(self):
        return np.array([ self.thickness - (i + 1) * self.s_tex_z for i in range(self.n_layers) ],
                      dtype = float)

    zz_t_i_arr = Property
    def _get_zz_t_i_arr(self):
        return self.thickness - self.z_t_i_arr
    #---------------------------------------------------------------
    # properties of the bending specimens (tree point bending test):
    #---------------------------------------------------------------

    # width of the cross section [m]
    #
    width = Float(0.20, geo_input = True)
    
    # number of rovings in 0-direction of one composite 
    # layer of the bending test [-]:
    #
    n_rovings = Int(23, geo_input = True)
    
    # cross section of one roving [mm**2]:
    #
    A_roving = Float(0.461, geo_input = True)

    #===========================================================================
    # Derived input properties
    #===========================================================================
    #-----------------------------
    # Concrete constitutive law
    # for simplified constant stress-strain-diagram of the concrete (EC2)
    #-----------------------------
    sig_c_mfn_block = Property(depends_on = '+sig_c_modified')
    @cached_property
    def _get_sig_c_mfn_block(self):
        '''simplified constant stress-strain-diagram of the concrete (EC2)
        '''
        print 'sig_c: block'
        #(for standard concrete)
        f_ck = self.f_ck + 8. 
        if f_ck <= 50:
            lamda = 0.8
            eta = 1.0  
            eps_cu3 = 0.0035
        # (for high strength concrete)
        #
        else:
            eta = 1.0 - (f_ck / 50.) / 200.
        # factor [-] to calculate the height of the compressive zone  
            lamda = 0.8 - (f_ck - 50.) / 400.
            eps_cu3 = (2.6 + 35. * ((90. - f_ck) / 100) ** 4.) / 1000. 
            
        xdata = np.hstack([0., (1. - lamda) * eps_cu3 - 0.00001, (1 - lamda) * eps_cu3, eps_cu3 ]) 
        ydata = np.hstack([0., 0., eta * (f_ck), eta * (f_ck), ])
       
        return MFnLineArray(xdata = xdata, ydata = ydata)   

    #-----------------------------
    # for bilinear stress-strain-diagram of the concrete (EC2)
    #-----------------------------
    sig_c_mfn_bilinear = Property(depends_on = '+sig_c_modified')
    @cached_property
    def _get_sig_c_mfn_bilinear(self):
        '''bilinear stress-strain-diagram of the concrete
        '''
        print 'sig_c: bilinear'
        #(for standard concrete)
        f_ck = self.f_ck + 8.  
        if f_ck <= 50.:
            eps_c3 = 0.00175
            eps_cu3 = 0.0035
        #(for high strength concrete)
        else :
            eps_c3 = (1.75 + 0.55 * (f_ck - 50.) / 40.) / 1000.
            eps_cu3 = (2.6 + 35 * (90 - f_ck) ** 4.) / 1000.      
        # concrete law with limit for eps_c
       
        xdata = np.hstack([0., eps_c3, eps_cu3])
        ydata = np.hstack([0., (f_ck), (f_ck)])
        
        return MFnLineArray(xdata = xdata, ydata = ydata)


    #-----------------------------
    # for quadratic stress-strain-diagram of the concrete
    #-----------------------------
    sig_c_mfn_quadratic = Property(depends_on = '+sig_c_modified')
    @cached_property
    def _get_sig_c_mfn_quadratic(self):
        '''quadratic stress-strain-diagram of the concrete
        '''
        print 'sig_c: quadratic'
        # (for all concretes up to f_cm=88 N/mm2) #max epislon_c1u
        f_cm = self.f_ck + 8 
        E_tan = 22. * (f_cm / 10) ** 0.3 * 1000.
#        print 'E_tan', E_tan
        eps_c1 = min(0.7 * f_cm ** 0.31, 2.8) / 1000. #EC2
        # @todo: with constant value this yields negative values for strains close to 'eps_c1u'
#        eps_c1 = 0.0022 #Brockmann
#        print 'eps_c1',eps_c1  #EC 0.0022
        E_sec = f_cm / eps_c1   
       
        if self.f_ck <= 50.:
            eps_c1u = 0.0035
            eps_arr = np.linspace(0., eps_c1u, num = 100.)
        
        elif self.f_ck > 50.:   
            eps_c1u = (2.8 + 27. * (((98. - f_cm) / 100.) ** 4.)) / 1000.
            print 'eps_c1u', eps_c1u
            eps_arr = np.linspace(0., eps_c1u, num = 100.) 

        k = E_tan / E_sec
        sig_c_arr = ((k * eps_arr / eps_c1 - (eps_arr / eps_c1) ** 2.) / (1. + (k - 2.) * eps_arr / eps_c1)) * f_cm
        
        #print 'sig_c_arr[0]',sig_c_arr[0]
        xdata = eps_arr  
        ydata = sig_c_arr
        
        return MFnLineArray(xdata = xdata, ydata = ydata)


    sig_c_config = Trait('quadratic', {'quadratic' : 'sig_c_mfn_quadratic',
                                      'bilinear'  : 'sig_c_mfn_bilinear',
                                      'block'     : 'sig_c_mfn_block'
                                      }, sig_c_modified = True)

    sig_c_mfn = Property(depends_on = '+sig_c_modified')
    def _get_sig_c_mfn(self):
        return getattr(self, self.sig_c_config_)

    sig_c_mfn_vect = Property(depends_on = '+sig_c_modified')
    @cached_property
    def _get_sig_c_mfn_vect(self):
        return np.frompyfunc(self.sig_c_mfn.get_value, 1, 1)

    # number of subdivisions of the compressive zone
    #
    n_c = Float(20., sig_c_modified = True)

    xi_arr = Property(Array, depends_on = '+sig_c_modified')
    @cached_property
    def _get_xi_arr(self):
        '''subdivide the compression zone 'x' in 'n_c' sub-areas;
        'xi_arr' giving the fraction of each distance of the sub-area from the upper surface 
        with respect to the compressive zone 'x'
        '''
        return np.arange(self.n_c) / self.n_c + 1. / (2. * self.n_c)

    #===========================================================================
    # Calculation of stress resultants
    #===========================================================================
    def calib_layer_response(self, u):
        '''CALIBRATION: derive the unknown constitutive law of the layer
        (effective crack bridge law)
        '''
        print 'calib_layer_response called'
        
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
        x = abs(self.eps_c_fail) / (abs(self.eps_c_fail) + abs(eps_t)) * thickness

        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = math.fabs(eps_t) / (thickness - x) * (z_t_i_arr - x)

        # use a ramp function to consider only negative strains
        #eps_c_i_arr = -eps_i_arr[ np.where(eps_i_arr <= 0.0)]
        eps_c_i_arr = np.fabs(-np.fabs(eps_i_arr) + eps_i_arr) / 2.0

        # use a ramp function to consider only positive strains
        #eps_t_i_arr = eps_i_arr[ np.where(eps_i_arr <= 0.0)]
        eps_t_i_arr = (np.fabs(eps_i_arr) + eps_i_arr) / 2.0        
       
        # ------------------------------------------------------------------------                
        # function used to approximate the effective crack bridge law
        # ------------------------------------------------------------------------                
        # linear cb-law using only the effective yarn modulus
        if self.ecb_law == 'linear':
            E_yarn = abs(var_a)
            # with limit for eps_tex
            #
            eps_fail = eps_t_i_arr[0]
            sig_fail = E_yarn * eps_fail
        
            eps_arr = np.array([ 0., eps_fail])
            sig_tex_arr = E_yarn * eps_arr
            #E_yarn 131176.670505
            print 'E_yarn', E_yarn
            print 'sig_fail', sig_fail
         
        if self.ecb_law == 'fbm':
            ecbl_fbm = ECBLFBM(sig_u_tt = self.sig_tex_fail)
            eps_fail = eps_t_i_arr[0]
            eps_arr = np.linspace(0, eps_fail, num = 100.)            
            sig_tex_arr = ecbl_fbm(eps_fail, var_a, eps_arr)
         
        elif self.ecb_law == 'root1':
            # use strain at the lowest textile layer as rupture strain 
            eps_fail = eps_t_i_arr[0]
            sig_fail = np.sqrt(eps_t_i_arr[0] * var_a)
            
            eps_arr = np.linspace(0, eps_fail, num = 100.)
            sig_tex_arr = np.sqrt((var_a) * eps_arr) 

#        elif self.ecb_law == 'root2':
#            # use strain at the lowest textile layer as rupture strain 
#            eps_fail = eps_t_i_arr[0]
#            sig_fail = self.sig_tex_fail
#            eps_max = var_a
#            eps_arr = np.linspace(0, eps_fail, num = 100.)
#            var_c = (2. * sig_fail**2.) / (eps_max**2.)
#            var_b = (sig_fail **2 - var_c * eps_max **2) / eps_max
#            sig_tex_arr = np.sqrt((var_c) * eps_arr ** 2. + (var_b) * eps_arr) 
            
        elif self.ecb_law == 'root2':
            # use strain at the lowest textile layer as rupture strain 
            eps_fail = eps_t_i_arr[0]
            sig_fail = self.sig_tex_fail
            eps_arr = np.linspace(0, eps_fail, num = 100.)
            var_b = (sig_fail ** 2. - (var_a) * eps_fail) / eps_fail ** 2. 
            sig_tex_arr = np.sqrt((var_b) * eps_arr ** 2. + (var_a) * eps_arr) 
            
        elif self.ecb_law == 'root3':
            # use strain at the lowest textile layer as rupture strain 
            eps_fail = eps_t_i_arr[0]
            sig_fail = self.sig_tex_fail
            eps_arr = np.linspace(0, eps_fail, num = 100.)
            var_b = ((sig_fail + var_a) ** 2. - var_a ** 2.) / eps_fail
            sig_tex_arr = -var_a + np.sqrt(var_a ** 2. + var_b * eps_arr) 
                
        # full plastic cb-law 
        elif self.ecb_law == 'plastic':
            E_yarn = abs(var_a)
            eps_fail = eps_t_i_arr[0]
            sig_fail = E_yarn * eps_t_i_arr[0]
            eps_arr = np.hstack([0., 0.01 * eps_fail, eps_fail ])
            sig_tex_arr = np.hstack([0., 0.99 * sig_fail, sig_fail])
            print 'E_yarn', E_yarn
            
        # linear cb-law up to reaching the rupture strain, behaves plastic afterwards
        elif self.ecb_law == 'bilinear':
            E_yarn = abs(var_a)
            sig_fail = self.sig_tex_fail
            eps_fail = sig_fail / E_yarn
            eps_arr = np.array([0., eps_fail, 2.*eps_fail ])
            sig_tex_arr = np.array([0., sig_fail, sig_fail ])
            print 'E_yarn', E_yarn
            
        # quadratic cb-law up to reaching the rupture strain, behaves plastic afterwards
        elif self.ecb_law == 'quadratic':
            # use strain at the lowest textile layer as rupture strain 
            eps_fail = eps_t_i_arr[0]
            sig_fail = self.sig_tex_fail
            eps_arr = np.linspace(0, eps_fail, num = 100.)
            var_b = -eps_fail * 2 * var_a 
            sig_tex_arr = var_a * eps_arr ** 2 + var_b * eps_arr 
        
#        elif self.ecb_law == 'quadratic_eps_max':
#            # use strain at the lowest textile layer as rupture strain 
#            eps_fail = eps_t_i_arr[0]
#            eps_max = abs(var_a)
#            sig_fail = self.sig_tex_fail
#            eps_arr = np.linspace(0, eps_fail, num = 100.)
#            var_c = - sig_fail / (eps_max **2.)
#            var_b = -2. * eps_max * var_c 
#            sig_tex_arr = var_c * eps_arr ** 2 + var_b * eps_arr 
#            xdata = np.hstack([ eps_arr,eps_arr[-1] ])  
#            ydata = np.hstack([ sig_tex_arr, sig_tex_arr[-1] ]) 
#            print 'xdata',xdata
#            print 'ydata',ydata
#            print 'var_b', var_b 
#            print 'var_a', var_a
#            print 'var_c',var_c
        
        elif self.ecb_law == 'quadratic_eps_max':
            # use strain at the lowest textile layer as rupture strain 
            eps_fail = eps_t_i_arr[0]
            d_eps = var_a
            sig_fail = self.sig_tex_fail
            eps_arr = np.linspace(0, eps_fail, num = 100.)
            var_c = -sig_fail / ((eps_fail + var_a) ** 2.)
            var_b = -2. * (eps_fail + d_eps) * var_c 
            sig_tex_arr = var_c * eps_arr ** 2 + var_b * eps_arr 
                
        elif self.ecb_law == 'quadratic_monoton':
            # use strain at the lowest textile layer as rupture strain 
            eps_fail = eps_t_i_arr[0]
            sig_fail = self.sig_tex_fail
            eps_arr = np.linspace(0, eps_fail, num = 100.)
           
            var_b = sig_fail / eps_fail / (1 - eps_fail / (2. *(eps_fail + abs(var_a))))
            var_c = -var_b / (2 * (eps_fail + abs(var_a)))
           
            sig_tex_arr = var_c * eps_arr ** 2 + abs(var_b) * eps_arr 

        elif self.ecb_law == 'quadratic_TT':
            # use strain at the lowest textile layer as rupture strain 
            eps_fail = eps_t_i_arr[0]
            sig_fail = self.sig_tex_fail
            eps_arr = np.linspace(0, eps_fail, num = 100.)

            var_b = (sig_fail - var_a * eps_fail ** 2) / eps_fail 
            sig_tex_arr = var_a * eps_arr ** 2 + var_b * eps_arr

        if self.ecb_law == 'cubic':
            # use strain at the lowest textile layer as rupture strain 
            eps_fail = eps_t_i_arr[0]
            sig_fail = self.sig_tex_fail
          
            eps_arr = np.linspace(0, eps_fail, num = 100.)
            
            #for horizontal tangent at eps_fail
            var_b = -(sig_fail + 2. * var_a * eps_fail ** 3.) / eps_fail ** 2. 
            
            var_c = -3. * var_a * eps_fail ** 2. - 2. * var_b * eps_fail 
            #for horizontal tangent at 1.5 eps_fail
            #var_b = -( sig_fail + 5.75 * var_a * eps_fail**3.) / (2. * eps_fail**2.) 
            
            #var_c = -3. * var_a * (1.5 * eps_fail) ** 2. - 2. * var_b * (1.5 *eps_fail)
            #for horizontal tangent at 3 eps_fail
            #var_b = -( sig_fail + 26 * var_a * eps_fail**3.) / (5. * eps_fail**2.) 
            
            #var_c = -3. * var_a * (3 * eps_fail) ** 2. - 2. * var_b * (3 *eps_fail)
            sig_tex_arr = var_a * eps_arr ** 3. + var_b * eps_arr ** 2. + var_c * eps_arr 
        
            
        print'eps_t1', eps_t_i_arr[0]
        
        # function for cb-law 
        xdata = eps_arr
        ydata = sig_tex_arr
        sig_t_mfn = MFnLineArray(xdata = xdata, ydata = ydata)

#        print 'xdata',xdata
#        print 'ydata',ydata
#        print 'var_b', var_b 
#        print 'var_a', var_a

#        # get the slope of the cb-function to check for monotony
#        #
#        diff_arr = sig_t_mfn.get_diffs( sig_t_mfn.xdata )
#        #print 'diff_arr ',diff_arr 
#        p_min = min( np.min(diff_arr), 0. )
#        if p_min != 0.:
#            # introduce penalty factor if monotony is not obayed
#            #
#            k_penalty = 1.0
#            sig_t_mfn.ydata = k_penalty * sig_t_mfn.ydata
        
        return x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn, eps_t, var_a

    def get_f_i_arr(self, u):
        '''force at the height of each reinforcement layer [kN]:
        '''
        x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn, u1, u2 = self.calib_layer_response(u)
                
        get_sig_t_i_arr = np.frompyfunc(sig_t_mfn.get_value, 1, 1)
        sig_t_i_arr = get_sig_t_i_arr(eps_t_i_arr)

        # tensile force of one reinforced composite layer [kN]:
        #
        n_rovings = self.n_rovings
        A_roving = self.A_roving
        f_t_i_arr = sig_t_i_arr * n_rovings * A_roving / 1000. 
    
        # compressive stress in the compression zone 'x' in each sub-area ('A_c' divided in 'n_c' parts)
        #
        sig_c_i_arr = np.array([self.sig_c_mfn.get_value(eps_c_i) for eps_c_i in eps_c_i_arr])

        # visualize signed sig_c_i_arr in[kN]
        #
        f_c_i_arr = sig_c_i_arr * self.width * self.s_tex_z * 1000.

        return f_t_i_arr, f_c_i_arr

    def get_eps_c_arr(self, u):
        '''get compression strain at each integration layer of the compressive zone [-]:
        for 'stress_case' flexion
        '''
        x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn, eps_t, eps_c = self.calib_layer_response(u)

        # for calibration us measured compressive strain
        # @todo: use mapped traits instead
        #
        eps_c = self.eps_c_fail
        
        eps_c_arr = (1. - self.xi_arr) * abs(eps_c)
        return x, eps_c_arr

    def get_sig_c_arr(self, u):
        # @todo: get vectorized version running: 
#        sig_c_i_arr = self.sig_c_mfn_vect( eps_c_i_arr )
#        get_sig_c_i_arr = np.frompyfunc( self.sig_c_mfn.get_value, 1, 1 )
#        sig_c_i_arr = get_sig_c_i_arr( eps_c_i_arr )
        x, eps_c_arr = self.get_eps_c_arr(u)
        sig_c_arr = np.array([self.sig_c_mfn.get_value(eps_c_i) for eps_c_i in eps_c_arr])
        return self.thickness - x * self.xi_arr, sig_c_arr
                
    # iteration counter
    #
    n = 0
    def get_lack_of_fit(self, u):
        '''Return the difference between 'N_external' and 'N_internal' as well as 'M_external' and 'M_internal'
        N_c (=compressive force of the compressive zone of the concrete) 
        N_t (=total tensile force of the reinforcement layers) 

        NOTE: eps_t (=tensile strain at the bottom [MPa]) is the search parameter
        to be found iteratively!
        '''
        # all stress cases refer to a stress configuration in which M is positive 
        #(maximum eps_c at the top, maximum eps_t at the bottom of the cross-section 
        #
        print '------------- iteration: %g ----------------------------' % (self.n)

        M_external = math.fabs(self.M_fail)
        N_external = self.N_fail
        x, eps_c_arr = self.get_eps_c_arr(u)
        f_t_i_arr, f_c_i_arr = self.get_f_i_arr(u)
        z_t_i_arr = self.z_t_i_arr

        # total tensile force of all layers of the composite tensile zone [kN]:
        # (characteristic value)
        #
        N_tk = sum(f_t_i_arr)
        print 'N_tk', N_tk 

        # @todo: get vectorized version running: 
#        sig_c_i_arr = self.sig_c_mfn_vect( eps_c_i_arr )
#        get_sig_c_i_arr = np.frompyfunc( self.sig_c_mfn.get_value, 1, 1 )
#        sig_c_i_arr = get_sig_c_i_arr( eps_c_i_arr )
        sig_c_arr = np.array([self.sig_c_mfn.get_value(eps_c_i) for eps_c_i in eps_c_arr])
        f_c_arr = sig_c_arr * self.width * self.xi_arr[0] * 2. * x * 1000.
        
        # 'z_c_arr' meassured from the top surface to the resulting compressive forces in 'f_c_arr'
        #
        z_c_arr = x * self.xi_arr

        N_ck = sum(f_c_arr)
        print 'N_ck', N_ck 

#        # total compressive force of the composite compressive zone [kN]:
#        # (design value)
#        #
#        N_cd = x * self.width * self.chi * self.f_cd * 1000.
#        print 'N_cd', N_cd

        # absolute error (equilibrium of sum N = 0):
        #
        
        # resistance moment of one reinforced composite layer [kNm]:
        # (characteristic value)
        #
        M_tk = np.dot(f_t_i_arr, z_t_i_arr)
        print 'M_tk', M_tk

        M_ck = np.dot(f_c_arr, z_c_arr)
        print 'M_ck', M_ck
                  
        N_internal = -N_ck + N_tk
        
        # moment evaluated with respect to the upper layer
        # denoted by underline character
        #
        M_internal_ = M_tk - M_ck

        # moment evaluated with respect to the center line
        #
        M_internal = M_internal_ - N_internal * self.thickness / 2.
        print 'M_internal', M_internal

        z0_t_i_arr = z_t_i_arr - self.thickness / 2 
        z0_c_arr = z_c_arr - self.thickness / 2 
        M0_tk = np.dot(f_t_i_arr, z0_t_i_arr)
        print 'M0_tk', M0_tk
        M0_ck = np.dot(f_c_arr, z0_c_arr)
        print 'M0_ck', M0_ck
        M0_internal = M0_tk - M0_ck
        print 'M0', M0_internal
        
        # set iteration counter
        #
        self.n += 1
        
        d_N = N_internal - N_external
        print 'd_N', d_N
         
        d_M = M_internal - M_external
        print 'd_M', d_M

        return np.array([ d_N, d_M ], dtype = float)

    def get_lack_of_fit_dN(self, eps_fail, var_a):
        fn = np.vectorize(lambda eps_fail, var_a : self.get_lack_of_fit([eps_fail, var_a])[0])
        return fn(eps_fail, var_a) 

    def get_lack_of_fit_dM(self, eps_fail, var_a):
        fn = np.vectorize(lambda eps_fail, var_a : self.get_lack_of_fit([eps_fail, var_a])[1])
        return fn(eps_fail, var_a) 

    # solution vector returned by 'fit_response'
    #
    u_sol = Property(Array(Float))
    @cached_property
    def _get_u_sol(self):
        '''iterate 'eps_t' such that the lack of fit between the calculated
        normal forces in the tensile reinforcement and the compressive zone (concrete)
        is smaller then 'xtol' defined in function 'brentq'.
        NOTE: the method 'get_lack_of_fit' returns the relative error.
        '''

        # use scipy-functionality to get the iterated value of 'eps_t'
            # NOTE: get_lack_of_fit must have a sign change as a requirement
            # for the function call 'brentq' to work property. 

            # The method brentq has optional arguments such as
            #   'xtol'    - absolut error (default value = 1.0e-12)
            #   'rtol'    - relative error (not supported at the time)
            #   'maxiter' - maximum numbers of iterations used
            #
        xtol = 1.0e-5
        u0 = self.u0
        return fsolve(self.get_lack_of_fit, u0, xtol = xtol)

    ecb_law = Trait('quadratic',
                          {'linear'   : [ 0.01, 80000. ],
                           'fbm' : [0.014, 0.5 ],
                           'root1': [ 0.02, 90000. ],
                           'root2': [ 0.02, 500000 ],
                           'root3': [ 0.018, -80. ], # solved for: eps_tex= 0.13,var_a = 11000  plastic: var_a -520.290186234,var_b 26537508.6905,eps_t1 0.00803806678189
                           'plastic'  : [ 0.014, 50000. ],
                           'bilinear' : [ 0.010, 50000. ],
                           'quadratic': [ 0.010, 50000.],
                           'quadratic_monoton': [ 0.016, 0.016],
                           'quadratic_TT': [ 0.010, -900000.],
                           'quadratic_eps_max': [ 0.019, 0.001],
                           'cubic'    : [ 0.016, -5000000. ]},
                         config_modified = True)

    u0 = Property(Array(float), depends_on = '+config_modified')
    @cached_property
    def _get_u0(self):
        return np.array(self.ecb_law_, dtype = 'f')      

    #===========================================================================
    # Postprocessing for plotting
    #===========================================================================
    def get_sig_comp_i_arr(self):
        '''tensile stress at the height of each reinforcement layer [MPa]:
        '''
        f_t_i_arr, f_c_i_arr = self.get_f_i_arr(self.u_sol)
        return (f_t_i_arr - f_c_i_arr) / self.width / self.s_tex_z / 1000.

    def get_sig_max(self):
        return np.max(self.get_sig_comp_i_arr())

    def plot_sig_t_mfn(self):
        '''plot calibrated effective crack bridge law
        '''
        # graph shows sig_comp at the height of the textile layer in [MPa] 
        layer_arr = np.arange(self.n_layers)
        zz_t_arr = self.zz_t_i_arr
        sig_comp_i_arr = self.get_sig_comp_i_arr()
        p.bar(zz_t_arr, sig_comp_i_arr, 0.02 * self.thickness, align = 'center')
        zz_c_arr, sig_c_arr = self.get_sig_c_arr(self.u_sol)
        p.plot(zz_c_arr, -sig_c_arr, color = 'red')

    def plot_sig_c_mfn(self):
        '''plot concrete law
        '''
        # graph shows sig_c in [MPa] 
        eps_arr = self.sig_c_mfn.xdata
        sig_c_arr = self.sig_c_mfn.ydata
        p.plot(eps_arr, sig_c_arr)

    def plot_ecb_law(self, u = None):
        '''plot concrete law
        '''
        # graph shows sig_c in [MPa]
        if u == None:
            u = self.u_sol
        x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn, eps_t, eps_c = self.calib_layer_response(u)
        eps_arr = sig_t_mfn.xdata
        sig_t_arr = sig_t_mfn.ydata
        p.plot(eps_arr, sig_t_arr)

if __name__ == '__main__':

    #------------------------------------------------
    # 1) CALIBRATION:
    # get 'eps_t' and the parameter of the effective 
    # crack bridge function 'var_a' for a given 'eps_c_fail'
    #------------------------------------------------
    #
    print '\n'
    print 'setup ECBLCalib'
    print '\n'
    sig_fl_calib = ECBLCalib(# mean concrete strength after 9 days
                               # 7d: f_ck,cube = 62 MPa; f_ck,cyl = 62/1.2=52
                               # 9d: f_ck,cube = 66.8 MPa; f_ck,cyl = 55,7
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

                               # define shape of the crack-bridge law ('linear', 'bilinear' , 'quadratic','root3''root1','root2','cubic')
                               #
#                                ecb_law = 'root3',
#                               ecb_law = 'quadratic',
#                               ecb_law = 'linear',  
#                              ecb_law = 'quadratic_eps_max', # with k_penalty= 1.1; 1216 at eps_fail;with k_penalty= 1.05; 1216 at eps_fail
#                               ecb_law = 'quadratic_monoton', 
#                                ecb_law = 'quadratic_TT',  
#                               ecb_law = 'plastic', 
#                               ecb_law = 'cubic',
                               ecb_law = 'fbm',
                               # define shape of the concrete stress-strain-law ('block', 'bilinear' or 'quadratic')
                               #
#                              sig_c_config = 'block'       #cubic:   #eps_t_fail 0.0137706348812      
#                              sig_c_config = 'bilinear'              #eps_t_fail 0.0142981075345
                              sig_c_config = 'quadratic'             #eps_t_fail 0.0137279096658                              
                              )


    ## call calibration and get 'u_sol':
    
    print 'XXX CALIB XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    u_sol = sig_fl_calib.u_sol
    print 'u_sol', u_sol
    max_sig = sig_fl_calib.get_sig_max()     
    print 'u_sol', u_sol
    print 'eps_c_fail', sig_fl_calib.eps_c_fail
    print 'eps_t_fail', sig_fl_calib.eps_t_fail
    print 'max_sig', max_sig

    ### plot crack bridge law [MPa]:
    #
    p.subplot(2, 2, 1)
    sig_fl_calib.plot_ecb_law(u_sol)

    ### plot concrete law
#    sig_fl_calib.plot_sig_c_mfn()
    p.subplot(2, 2, 2)
    sig_fl_calib.plot_sig_t_mfn()

    u_grid = np.ogrid[0.01:0.02:30j, 0.1:2.0:30j]
    N_grid = sig_fl_calib.get_lack_of_fit_dN(u_grid[0], u_grid[1])
    M_grid = sig_fl_calib.get_lack_of_fit_dM(u_grid[0], u_grid[1])
    ones = np.ones_like(N_grid)
    p.subplot(2, 2, 3)
    p.pcolor(u_grid[0] * ones, u_grid[1] * ones, N_grid)
    p.colorbar()
    p.subplot(2, 2, 4)
    p.pcolor(u_grid[0] * ones, u_grid[1] * ones, M_grid)
    p.colorbar()
    p.show()
