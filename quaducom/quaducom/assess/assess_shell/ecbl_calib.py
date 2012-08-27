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
import etsproxy.mayavi.mlab as m

from mathkit.mfn import MFnLineArray

from scipy.optimize import fsolve

data_file_editor = FileEditor(filter = ['*.DAT'])

from matresdev.db.simdb import SimDB
simdb = SimDB()

from math import exp, log

class ECBLBase(HasTraits):
    '''Base class for Effective Crack Bridge Laws.'''
    
    u0 = List([0.0, 0.0])

class ECBLLinear(ECBLBase):
    '''Effective crack bridge Law with linear elastic response.'''
    
    sig_tex_u = Float

    u0 = List([ 0.01, 80000. ])
                                   
    def __call__(self, eps_tex_u, var_a):
        E_yarn = abs(var_a)

        # with limit for eps_tex
        #
        eps_tex_arr = np.array([ 0., eps_tex_u])
        sig_tex_arr = E_yarn * eps_tex_arr
        return eps_tex_arr, sig_tex_arr 
    
class ECBLFBM(ECBLBase):
    '''Effective crack bridge Law based on fiber-bundle-model.'''
    
    sig_tex_u = Float
        
    u0 = List([0.014, 0.5 ])
        
    def __call__(self, eps_tex_u, m):
        eps_tex_arr = np.linspace(0, eps_tex_u, num = 100.)
        sig_tex_arr = (self.sig_tex_u / eps_tex_u / exp(-pow(exp(-log(m) / m), 1.0 * m)) * 
                     eps_tex_arr * np.exp(-np.power(eps_tex_arr / eps_tex_u * exp(-log(m) / m), 1.0 * m)))            
        return eps_tex_arr, sig_tex_arr 
        
class ECBLCubic(ECBLBase):
    '''Effective crack bridge Law using a cubic polynomial.'''
    
    sig_tex_u = Float

    u0 = List([ 0.016, -5000000. ])
                                   
    def __call__(self, eps_tex_u, var_a):
        sig_tex_u = self.sig_tex_u
        eps_tex_arr = np.linspace(0, eps_tex_u, num = 100.)
        # for horizontal tangent at eps_tex_u
        var_b = -(sig_tex_u + 2. * var_a * eps_tex_u ** 3.) / eps_tex_u ** 2. 
        var_c = -3. * var_a * eps_tex_u ** 2. - 2. * var_b * eps_tex_u 
        sig_tex_arr = var_a * eps_tex_arr ** 3. + var_b * eps_tex_arr ** 2. + var_c * eps_tex_arr
        return eps_tex_arr, sig_tex_arr         

class ECBLPlastic(ECBLBase):
    '''Effective crack bridge Law using a cubic polynomial.'''
    
    sig_tex_u = Float

    u0 = List([ 0.014, 50000. ])
                                   
    def __call__(self, eps_tex_u, var_a):
        E_yarn = abs(var_a)
        sig_tex_u = self.sig_tex_u
        eps_tex_arr = np.hstack([0., 0.01 * eps_tex_u, eps_tex_u ])
        sig_tex_arr = np.hstack([0., 0.99 * sig_tex_u, sig_tex_u])
        print 'E_yarn', E_yarn
        return eps_tex_arr, sig_tex_arr         

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

    # ultimate textile stress measured in the tensile test [MPa]
    #
    sig_tex_u = Float(1216., sig_t_modified = True)
    
    # compressive strain at the top at rupture [-]
    # (positive value is used)
    # value measured at the experiment (with DMS)
    # for 0: about 3.3
    # for 90: about 3.0
    # ultimate strain theoretically (Brockman): about 4.5
    # NOTE: strain was meassured at a distance of 5 cm
    #
    eps_cu = Float(0.0033, sig_t_modified = True) # float value corresponds to 3 promile

    # ultimate limit strain at the tensile surface of the cross section [-]:
    # value is derived based on the calibrated crack bridge law
    #
    eps_tu = Float
    sig_tex_mfn = Instance(MFnLineArray)

    # rupture moment and normal force meassured in the calibration experiment
    # (three point bending test)
    #
    Mu = Float(3.5) # [kNm]
    Nu = Float(0.0) # [kN]

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
    
    # distance from the top of each reinforcement layer [m]:
    #
    z_ti_arr = Property(depends_on = '+geo_input')
    @cached_property
    def _get_z_ti_arr(self):
        return np.array([ self.thickness - (i + 1) * self.s_tex_z for i in range(self.n_layers) ],
                      dtype = float)

    zz_ti_arr = Property
    def _get_zz_ti_arr(self):
        return self.thickness - self.z_ti_arr
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
        # (for all concretes up to f_cm=88 N/mm2) #max epislon_c1u
        f_cm = self.f_ck + 8 
        E_tan = 22. * (f_cm / 10) ** 0.3 * 1000.
        eps_c1 = min(0.7 * f_cm ** 0.31, 2.8) / 1000. #EC2
        # @todo: with constant value this yields negative values for strains close to 'eps_c1u'
#        eps_c1 = 0.0022 #Brockmann
        E_sec = f_cm / eps_c1   
       
        if self.f_ck <= 50.:
            eps_c1u = 0.0035
            eps_arr = np.linspace(0., eps_c1u, num = 100.)
        
        elif self.f_ck > 50.:   
            eps_c1u = (2.8 + 27. * (((98. - f_cm) / 100.) ** 4.)) / 1000.
            eps_arr = np.linspace(0., eps_c1u, num = 100.) 

        k = E_tan / E_sec
        sig_c_arr = ((k * eps_arr / eps_c1 - (eps_arr / eps_c1) ** 2.) / (1. + (k - 2.) * eps_arr / eps_c1)) * f_cm
        
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
        return np.vectorize(self.sig_c_mfn.get_value)

    # number of subdivisions of the compressive zone
    #
    n_cj = Float(20., sig_c_modified = True)

    zeta_cj_arr = Property(Array, depends_on = '+sig_c_modified')
    @cached_property
    def _get_zeta_cj_arr(self):
        '''subdivide the compression zone 'x' in 'n_cj' sub-areas;
        'zeta_cj_arr' giving the fraction of each distance of the sub-area from the upper surface 
        with respect to the compressive zone 'x'
        '''
        return np.arange(self.n_cj) / self.n_cj + 1. / (2. * self.n_cj)

    #===========================================================================
    # Effective crack bridge law
    #===========================================================================
    ecbl_type = Trait('fbm', dict(fbm = ECBLFBM,
                                   cubic = ECBLCubic,
                                   linear = ECBLLinear,
                                   plastic = ECBLPlastic),
                      config_modified = True)
    
    ecbl = Property(depends_on = 'ecbl_type')
    @cached_property
    def _get_ecbl(self):
        return self.ecbl_type_(sig_tex_u = self.sig_tex_u)
    
    def get_ecbl_tex_arr(self, eps_tex_u, var_a):
        eps_tex_arr, sig_tex_arr = self.ecbl(eps_tex_u, var_a)    
        return eps_tex_arr, sig_tex_arr
    
    def get_ecbl_tex_mfn(self, eps_tex_u, var_a):
        eps_tex_arr, sig_tex_arr = self.get_ecbl_tex_arr(eps_tex_u, var_a)
        return MFnLineArray(xdata = eps_tex_arr, ydata = sig_tex_arr)

    #===========================================================================
    # Discretization conform to the tex layers
    #===========================================================================
    def get_eps_i_arr(self, eps_tu):
        '''CALIBRATION: derive the unknown constitutive law of the layer
        (effective crack bridge law)
        '''
        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_t'
        # ------------------------------------------------------------------------                
        thickness = self.thickness
        z_ti_arr = self.z_ti_arr
        
        # ------------------------------------------------------------------------                
        # derived params depending on value for 'eps_t'
        # ------------------------------------------------------------------------                

        # heights of the compressive zone:
        #
        x = abs(self.eps_cu) / (abs(self.eps_cu) + abs(eps_tu)) * thickness

        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = math.fabs(eps_tu) / (thickness - x) * (z_ti_arr - x)

        # use a ramp function to consider only negative strains
        #eps_ci_arr = -eps_i_arr[ np.where(eps_i_arr <= 0.0)]
        eps_ci_arr = np.fabs(-np.fabs(eps_i_arr) + eps_i_arr) / 2.0

        # use a ramp function to consider only positive strains
        #eps_ti_arr = eps_i_arr[ np.where(eps_i_arr <= 0.0)]
        eps_ti_arr = (np.fabs(eps_i_arr) + eps_i_arr) / 2.0     
        
        return x, eps_ti_arr, eps_ci_arr
        
    def get_f_ti_arr(self, eps_tu, var_a):
        '''force at the height of each reinforcement layer [kN]:
        '''
        x, eps_ti_arr, eps_ci_arr = self.get_eps_i_arr(eps_tu)
        eps_tex_u = eps_ti_arr[0]
                        
        sig_tex_mfn = self.get_ecbl_tex_mfn(eps_tex_u, var_a)
                
        get_sig_ti_arr = np.vectorize(sig_tex_mfn.get_value)
        sig_ti_arr = get_sig_ti_arr(eps_ti_arr)

        # tensile force of one reinforced composite layer [kN]:
        #
        n_rovings = self.n_rovings
        A_roving = self.A_roving
        return sig_ti_arr * n_rovings * A_roving / 1000. 

    #===========================================================================
    # Discretization of the compressive zone - tex layer conform
    #===========================================================================
    def get_f_ci_arr(self, eps_tu, var_a):
        '''Compressive stress in the compresive zone 'x' for each layer i.
        '''
        x, eps_ti_arr, eps_ci_arr = self.get_eps_i_arr(eps_tu)
        sig_ci_arr = np.array([self.sig_c_mfn.get_value(eps_ci) for eps_ci in eps_ci_arr])
        return sig_ci_arr * self.width * self.s_tex_z * 1000.

    #===========================================================================
    # Fine discretization of the compressive zone
    #===========================================================================
    def get_eps_cj_arr(self, eps_tu):
        '''get compressive strain at each integration layer of the compressive zone [-]:
        for 'stress_case' flexion
        '''
        x, eps_ti_arr, eps_ci_arr = self.get_eps_i_arr(eps_tu)

        # for calibration us measured compressive strain
        # @todo: use mapped traits instead
        #
        eps_c = self.eps_cu
        
        eps_cj_arr = (1. - self.zeta_cj_arr) * abs(eps_c)
        return x, eps_cj_arr

    def get_sig_cj_arr(self, eps_tu):
        x, eps_cj_arr = self.get_eps_cj_arr(eps_tu)
        return x, self.sig_c_mfn_vect(eps_cj_arr)

    def get_f_cj_arr(self, eps_tu):
        x, sig_cj_arr = self.get_sig_cj_arr(eps_tu)
        return sig_cj_arr * self.width * x / self.n_cj * 1000. 
                
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

        eps_tu, var_a = u
        M_external = math.fabs(self.Mu)
        N_external = self.Nu
        x, eps_cj_arr = self.get_eps_cj_arr(eps_tu)
        f_ti_arr = self.get_f_ti_arr(eps_tu, var_a)
        z_ti_arr = self.z_ti_arr

        # total tensile force of all layers of the composite tensile zone [kN]:
        # (characteristic value)
        #
        N_tk = sum(f_ti_arr)

        f_cj_arr = self.get_f_cj_arr(eps_tu)
        
        # 'z_c_arr' measured from the top surface to the resulting compressive forces in 'f_c_arr'
        #
        z_cj_arr = x * self.zeta_cj_arr

        N_ck = sum(f_cj_arr)

        # resistance moment of one reinforced composite layer [kNm]:
        # (characteristic value)
        #
        M_tk = np.dot(f_ti_arr, z_ti_arr)
        M_ck = np.dot(f_cj_arr, z_cj_arr)
                  
        N_internal = -N_ck + N_tk
        
        # moment evaluated with respect to the upper layer
        # denoted by underline character
        #
        M_internal_ = M_tk - M_ck

        # moment evaluated with respect to the center line
        #
        M_internal = M_internal_ - N_internal * self.thickness / 2.

        # set iteration counter
        #
        self.n += 1
        d_N = N_internal - N_external
        d_M = M_internal - M_external

        return np.array([ d_N, d_M ], dtype = float)

    def get_lack_of_fit_dN(self, eps_tu, var_a):
        fn = np.vectorize(lambda eps_tu, var_a : self.get_lack_of_fit([eps_tu, var_a])[0])
        return fn(eps_tu, var_a) 

    def get_lack_of_fit_dM(self, eps_tu, var_a):
        fn = np.vectorize(lambda eps_tu, var_a : self.get_lack_of_fit([eps_tu, var_a])[1])
        return fn(eps_tu, var_a) 

    # solution vector returned by 'fit_response'
    #
    u_sol = Property(Array(Float), depends_on = '+config_modified')
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
        return fsolve(self.get_lack_of_fit, self.u0, xtol = 1.0e-5)

    u0 = Property(Array(float), depends_on = '+config_modified')
    @cached_property
    def _get_u0(self):
        return self.ecbl.u0

    #===========================================================================
    # Postprocessing for plotting
    #===========================================================================
    def get_sig_comp_i_arr(self):
        '''tensile stress at the height of each reinforcement layer [MPa]:
        '''
        f_t_i_arr, f_c_i_arr = self.get_f_i_arr(*self.u_sol)
        return (f_t_i_arr - f_c_i_arr) / self.width / self.s_tex_z / 1000.

    def get_sig_max(self):
        return np.max(self.get_sig_comp_i_arr())

    def plot_sig_comp_i(self):
        '''plot calibrated effective crack bridge law
        '''
        # graph shows sig_comp at the height of the textile layer in [MPa] 
        zz_t_arr = self.zz_ti_arr
        sig_comp_i_arr = self.get_sig_comp_i_arr()
        p.bar(zz_t_arr, sig_comp_i_arr, 0.02 * self.thickness, align = 'center')
        zz_c_arr, sig_c_arr = self.get_sig_c_arr(self.u_sol[0])
        p.plot(zz_c_arr, -sig_c_arr, color = 'red')

    def plot_sig_c_mfn(self):
        '''plot concrete law
        '''
        # graph shows sig_c in [MPa] 
        eps_c_arr = self.sig_c_mfn.xdata
        sig_c_arr = self.sig_c_mfn.ydata
        p.plot(eps_c_arr, sig_c_arr)

    def plot_ecb_law(self, u = None):
        '''plot concrete law
        '''
        # graph shows sig_c in [MPa]
        if u == None:
            u = self.u_sol
        eps_tex_arr, sig_tex_arr = self.get_ecbl_tex_arr(*u)
        p.plot(eps_tex_arr, sig_tex_arr)

if __name__ == '__main__':

    #------------------------------------------------
    # 1) CALIBRATION:
    # get 'eps_t' and the parameter of the effective 
    # crack bridge function 'var_a' for a given 'eps_cu'
    #------------------------------------------------
    #
    print '\n'
    print 'setup ECBLCalib'
    print '\n'
    ecbl_calib = ECBLCalib(# mean concrete strength after 9 days
                             # 7d: f_ck,cube = 62 MPa; f_ck,cyl = 62/1.2=52
                             # 9d: f_ck,cube = 66.8 MPa; f_ck,cyl = 55,7
                             f_ck = 55.7,

                               # measured strain at bending test rupture (0-dir)
                               #
                               eps_cu = 3.3 / 1000.,

                               # measured value in bending test [kNm]
                               # value per m: M = 5*3.49
                               #
                               Mu = 3.49,
                           
                               # values for experiment beam with width = 0.20 m
                               #
                               width = 0.20,
                               n_roving = 23.,

                               ecbl_type = 'linear',
                               # define shape of the concrete stress-strain-law ('block', 'bilinear' or 'quadratic')
                               #
#                              sig_c_config = 'block'       #cubic:   #eps_tu 0.0137706348812      
#                              sig_c_config = 'bilinear'              #eps_tu 0.0142981075345
                              sig_c_config = 'quadratic'             #eps_tu 0.0137279096658                              
                              )


    ## call calibration and get 'u_sol':
    
    for ecbl_type in ['linear', 'cubic', 'fbm']:
        print 'XXX CALIB ', ecbl_type
        ecbl_calib.n = 0
        ecbl_calib.ecbl_type = ecbl_type
        ecbl_calib.plot_ecb_law()

    p.show()
#    ### plot crack bridge law [MPa]:
#    #
#    p.subplot(1, 2, 1)
#    sig_fl_calib.plot_ecb_law(u_sol)
#
#    ### plot concrete law
##    sig_fl_calib.plot_sig_c_mfn()
#    p.subplot(1, 2, 2)
#    sig_fl_calib.plot_sig_comp_i()
#


#    u_grid = np.ogrid[0.01:0.02:30j, 0.1:2.0:30j]
#    N_grid = sig_fl_calib.get_lack_of_fit_dN(u_grid[0], u_grid[1])
#    M_grid = sig_fl_calib.get_lack_of_fit_dM(u_grid[0], u_grid[1])
#    ones = np.ones_like(N_grid)
#    p.subplot(2, 2, 3)
#    p.pcolor(u_grid[0] * ones, u_grid[1] * ones, N_grid)
#    p.colorbar()
#    p.subplot(2, 2, 4)
#    p.pcolor(u_grid[0] * ones, u_grid[1] * ones, M_grid)
#    p.colorbar()
    p.show()
#
#    print 'N_grid', N_grid
#    print 'u_grid[0] ', u_grid[0]
#    print 'u_grid[1] ', u_grid[1]
#
#
