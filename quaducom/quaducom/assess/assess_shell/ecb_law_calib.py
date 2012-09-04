'''
Created on Jun 23, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    Float, Instance, Array, Property, cached_property, \
    HasTraits, Trait

import math
import numpy as np
import pylab as p

from mathkit.mfn import MFnLineArray

from scipy.optimize import fsolve

from ecb_cross_section import \
    ECBCrossSection
    
from ecb_law import \
    ECBLLinear, ECBLFBM, ECBLCubic, ECBLPlastic

from cc_law import \
    CCLawBlock, CCLawLinear, CCLawQuadratic

from matresdev.db.simdb import SimDB
simdb = SimDB()

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

    # rupture moment and normal force meassured in the calibration experiment
    # (three point bending test)
    #
    Mu = Float(3.5) # [kNm]
    Nu = Float(0.0) # [kN]

    #===========================================================================
    # Cross Section Specification (Geometry and Layout)
    #===========================================================================

    cs = Instance(ECBCrossSection)

    #===========================================================================
    # Compressive concrete constitutive law
    #===========================================================================
    cc_law_type = Trait('constant', dict(constant = CCLawBlock,
                                         linear = CCLawLinear,
                                         quadratic = CCLawQuadratic),
                      config_modified = True)
    
    cc_law = Property(depends_on = 'cc_law_type')
    @cached_property
    def _get_cc_law(self):
        '''Construct the compressive concrete law'''
        return self.cc_law_type_(f_ck = self.f_ck)

    sig_c_mfn = Property(depends_on = '+sig_c_modified')
    def _get_sig_c_mfn(self):
        return self.cc_law.mfn

    sig_c_mfn_vct = Property(depends_on = '+sig_c_modified')
    @cached_property
    def _get_sig_c_mfn_vct(self):
        return np.vectorize(self.sig_c_mfn.get_value)

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
        '''Get the arrays sigma - epsilon defining the crack bridge law.
        '''
        eps_tex_arr, sig_tex_arr = self.ecbl(eps_tex_u, var_a)    
        return eps_tex_arr, sig_tex_arr
    
    def get_ecbl_tex_mfn(self, eps_tex_u, var_a):
        '''Get the callable function for effective crack brige law.
        '''
        eps_tex_arr, sig_tex_arr = self.get_ecbl_tex_arr(eps_tex_u, var_a)
        return MFnLineArray(xdata = eps_tex_arr, ydata = sig_tex_arr)

    ecbl_mfn = Property(depends_on = '+config_modified')
    @cached_property
    def _get_ecbl_mfn(self):
        return self.get_ecbl_tex_mfn(*self.u_sol)
    
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
        thickness = self.cs.thickness
        z_ti_arr = self.cs.z_ti_arr
        
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

    #===========================================================================
    # Layer conform discretization of the tensile zone
    #===========================================================================
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
        n_rovings = self.cs.n_rovings
        A_roving = self.cs.A_roving
        return sig_ti_arr * n_rovings * A_roving / 1000. 

    #===========================================================================
    # Discretization of the compressive zone - tex layer conform
    #===========================================================================
    def get_f_ci_arr(self, eps_tu, var_a):
        '''Compressive stress in the compresive zone 'x' for each layer i.
        '''
        x, eps_ti_arr, eps_ci_arr = self.get_eps_i_arr(eps_tu)
        sig_ci_arr = np.array([self.sig_c_mfn.get_value(eps_ci) for eps_ci in eps_ci_arr])
        return sig_ci_arr * self.cs.width * self.cs.s_tex_z * 1000.

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
        
        eps_cj_arr = (1. - self.cs.zeta_cj_arr) * abs(eps_c)
        return x, eps_cj_arr

    def get_sig_cj_arr(self, eps_tu):
        x, eps_cj_arr = self.get_eps_cj_arr(eps_tu)
        return x, self.sig_c_mfn_vct(eps_cj_arr)

    def get_f_cj_arr(self, eps_tu):
        x, sig_cj_arr = self.get_sig_cj_arr(eps_tu)
        return sig_cj_arr * self.cs.width * x / self.cs.n_cj * 1000. 
                
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
        #print '------------- iteration: %g ----------------------------' % (self.n)

        eps_tu, var_a = u
        M_external = math.fabs(self.Mu)
        N_external = self.Nu
        x, eps_cj_arr = self.get_eps_cj_arr(eps_tu)
        f_ti_arr = self.get_f_ti_arr(eps_tu, var_a)
        z_ti_arr = self.cs.z_ti_arr

        # total tensile force of all layers of the composite tensile zone [kN]:
        # (characteristic value)
        #
        N_tk = sum(f_ti_arr)

        f_cj_arr = self.get_f_cj_arr(eps_tu)
        
        # 'z_c_arr' measured from the top surface to the resulting compressive forces in 'f_c_arr'
        #
        z_cj_arr = x * self.cs.zeta_cj_arr

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
        M_internal = M_internal_ - N_internal * self.cs.thickness / 2.

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
        f_ti_arr = self.get_f_ti_arr(*self.u_sol)
        f_ci_arr = self.get_f_ci_arr(*self.u_sol)
        return (f_ti_arr - f_ci_arr) / self.cs.width / self.cs.s_tex_z / 1000.

    def get_sig_max(self):
        return np.max(self.get_sig_comp_i_arr())

    def plot_sig_comp_i(self):
        '''plot calibrated effective crack bridge law
        '''
        # graph shows sig_comp at the height of the textile layer in [MPa] 
        zz_t_arr = self.cs.zz_ti_arr
        sig_comp_i_arr = self.get_sig_comp_i_arr()
        p.bar(zz_t_arr, sig_comp_i_arr, 0.02 * self.cs.thickness, align = 'center')
        zz_c_arr, sig_c_arr = self.get_sig_c_arr(self.u_sol[0])
        p.plot(zz_c_arr, -sig_c_arr, color = 'red')

    def plot_sig_c_mfn(self):
        '''plot concrete law
        '''
        # graph shows sig_c in [MPa] 
        eps_c_arr = self.sig_c_mfn.xdata
        sig_c_arr = self.sig_c_mfn.ydata
        p.plot(eps_c_arr, sig_c_arr)

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
    p.plot([0, 0], [0, 2.4e3])
    
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
                           cs = ECBCrossSection(
                                                width = 0.20,
                                                n_roving = 23.,
                                                ),
                           ecbl_type = 'linear',

                           # define shape of the concrete stress-strain-law ('block', 'bilinear' or 'quadratic')
                           #
#                              sig_c_config = 'block'       #cubic:   #eps_tu 0.0137706348812      
#                              sig_c_config = 'bilinear'              #eps_tu 0.0142981075345
                          sig_c_config = 'quadratic'             #eps_tu 0.0137279096658                              
                          )

    do = 'plot_cs_state'

    if do == 'plot_ecbl':
        for ecbl_type in ['linear', 'cubic', 'fbm', 'plastic']:
            print 'CALIB TYPE', ecbl_type
            ecbl_calib.n = 0
            ecbl_calib.ecbl_type = ecbl_type
            ecbl_calib.ecbl_mfn.plot(p)
            print 'E_yarn', ecbl_calib.ecbl_mfn.get_diff(0.00001)
            print 'INTEG', ecbl_calib.ecbl_mfn.integ_value

    elif do == 'plot_cs_state':    
    #    ### plot crack bridge law [MPa]:
    #    #
    #    p.subplot(1, 2, 1)
    #    sig_fl_calib.plot_ecb_law(u_sol)
    #
        ### plot concrete law
        ecbl_calib.plot_sig_c_mfn()
    #    p.subplot(1, 2, 2)
        ecbl_calib.plot_sig_comp_i()
    #
    elif do == 'plot_MN_grid':
    
        u_grid = np.ogrid[0.01:0.02:30j, 0.1:2.0:30j]
        N_grid = ecbl_calib.get_lack_of_fit_dN(u_grid[0], u_grid[1])
        M_grid = ecbl_calib.get_lack_of_fit_dM(u_grid[0], u_grid[1])
        ones = np.ones_like(N_grid)
        p.subplot(2, 2, 3)
        p.pcolor(u_grid[0] * ones, u_grid[1] * ones, N_grid)
        p.colorbar()
        p.subplot(2, 2, 4)
        p.pcolor(u_grid[0] * ones, u_grid[1] * ones, M_grid)
        p.colorbar()
    
        print 'N_grid', N_grid
        print 'u_grid[0] ', u_grid[0]
        print 'u_grid[1] ', u_grid[1]
    
    
    p.show()
        
