'''
Created on Jun 23, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    Float, Instance, Array, Property, cached_property, \
    HasTraits, DelegatesTo

import math
import numpy as np
import pylab as p

from scipy.optimize import fsolve

from ecb_cross_section import \
    ECBCrossSection
    
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

    # rupture moment and normal force measured in the calibration experiment
    # (three point bending test)
    #
    Mu = Float(3.5) # [kNm]
    Nu = Float(0.0) # [kN]

    #===========================================================================
    # Cross Section Specification (Geometry and Layout)
    #===========================================================================

    cs = Instance(ECBCrossSection)
    
    ecbl = DelegatesTo('cs')

    u0 = Property(Array(float), depends_on = '+tt_input, ecbl_modified')
    @cached_property
    def _get_u0(self):
        return self.ecbl.u0

    def get_NM_internal(self, eps_lo, eps_up):
        '''
        NOTE: eps_t (=tensile strain at the bottom [MPa]) is the search parameter
        to be found iteratively!
        '''
        # all stress cases refer to a stress configuration in which M is positive 
        #(maximum eps_c at the top, maximum eps_t at the bottom of the cross-section 
        #
        #print '------------- iteration: %g ----------------------------' % (self.n)
        
        self.cs.set(eps_lo = eps_lo, eps_up = eps_up)
        self.ecbl.eps_tex_u = self.cs.eps_ti_arr[0]
        
        N_tk = sum(self.cs.f_ti_arr)
        M_tk = np.dot(self.cs.f_ti_arr, self.cs.z_ti_arr)

        N_ck = sum(self.cs.f_cj_arr)
        M_ck = np.dot(self.cs.f_cj_arr, self.cs.z_cj_arr)
                  
        N_internal = -N_ck + N_tk
        M_internal_ = M_tk - M_ck

        # moment evaluated with respect to the center line
        #
        M_internal = M_internal_ - N_internal * self.cs.thickness / 2.
        
        # return the internal stress resultants
        return N_internal, M_internal

    # iteration counter
    #
    n = 0
    def get_lack_of_fit(self, u):
        '''Return the difference between 'N_external' and 'N_internal' as well as 'M_external' and 'M_internal'
        N_c (=compressive force of the compressive zone of the concrete) 
        N_t (=total tensile force of the reinforcement layers) 
        '''
        self.n += 1
        # set iteration counter
        #
        self.ecbl.set_cparams(*u)
        N_internal, M_internal = self.get_NM_internal(eps_lo = u[0],
                                                      eps_up = self.eps_cu)

        M_external = math.fabs(self.Mu)
        N_external = self.Nu
        
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
    u_sol = Property(Array(Float), depends_on = '+config_modified, ecbl_modified')
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

    #===========================================================================
    # Calibrated ecbl_mfn
    #===========================================================================
    ecbl_mfn = Property(depends_on = '+config_modified, ecbl_modified')
    @cached_property
    def _get_ecbl_mfn(self):
        self.ecbl.set_cparams(*self.u_sol)
        return self.ecbl.mfn

    #===========================================================================
    # Postprocessing for plotting
    #===========================================================================
    def get_sig_comp_i_arr(self):
        '''tensile stress at the height of each reinforcement layer [MPa]:
        '''
        
        f_ti_arr = self.get_f_ti_arr(eps_ti_arr)
        f_ci_arr = self.get_f_ci_arr(eps_ci_arr)
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

    do = 'plot_ecbl'

    if do == 'plot_ecbl':
        for sig_tex_u, color in zip([1200, 1300, 1400], ['red', 'green', 'blue', 'black', 'orange', 'brown']):

            for ecbl_type in ['linear', 'cubic', 'fbm']:
            #for ecbl_type in ['fbm']:
                print 'CALIB TYPE', ecbl_type
                ecbl_calib.n = 0
                ecbl_calib.ecbl_type = ecbl_type
                ecbl_calib.ecbl.sig_tex_u = sig_tex_u
                ecbl_calib.ecbl_mfn.plot(p, color = color, linewidth = 8)
                print 'E_yarn', ecbl_calib.ecbl_mfn.get_diff(0.00001)
                print 'INTEG', ecbl_calib.ecbl_mfn.integ_value

            ecbl_calib.ecbl_type = 'bilinear'
            ecbl_calib.ecbl.sig_tex_u = sig_tex_u
            for eps_el_fraction in [0.999]: # np.linspace(0.25, 0.99999, 4):
                ecbl_calib.n = 0
                ecbl_calib.ecbl.eps_el_fraction = eps_el_fraction
                ecbl_calib.ecbl_mfn.plot(p, color = color)
                print 'E_yarn', ecbl_calib.ecbl_mfn.get_diff(0.00001)
                print 'INTEG', ecbl_calib.ecbl_mfn.integ_value
        p.plot([0.0, 0.01], [0.0, 2400], color = 'black')
            
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
        
