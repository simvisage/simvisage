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

    # rupture moment and normal force measured in the calibration experiment
    # (three point bending test)
    #
    Mu = Float(3.5) # [kNm]
    Nu = Float(0.0) # [kN]

    #===========================================================================
    # Cross Section Specification (Geometry and Layout)
    #===========================================================================

    cs = Instance(ECBCrossSection)
    def _cs_default(self):
        return ECBCrossSection()

    ecbl_type = DelegatesTo('cs')
    ecbl = DelegatesTo('cs')

    u0 = Property(Array(float), depends_on = 'cs.modified')
    @cached_property
    def _get_u0(self):
        u0 = self.ecbl.u0
        eps_lo = self.cs.convert_eps_tex_u_2_lo(u0[0])
        return np.array([eps_lo, u0[1] ], dtype = 'float')

    # iteration counter
    #
    n = 0
    def get_lack_of_fit(self, u):
        '''Return the difference between 'N_external' and 'N_internal' as well as 'M_external' and 'M_internal'
        N_c (=compressive force of the compressive zone of the concrete) 
        N_t (=total tensile force of the reinforcement layers) 
        '''

        print '--------------------iteration', self.n, '------------------------'
        self.n += 1
        # set iteration counter
        #
        self.cs.set(eps_lo = u[0], eps_up = -self.eps_cu)   
        eps_tex_u = self.cs.convert_eps_lo_2_tex_u(u[0])
        self.cs.ecbl.set_cparams(eps_tex_u, u[1])
        
        N_internal = self.cs.N
        M_internal = self.cs.M
        
        d_N = N_internal - self.Nu
        d_M = M_internal - self.Mu

        return np.array([ d_N, d_M ], dtype = float)

    # solution vector returned by 'fit_response'
    #
    u_sol = Property(Array(Float), depends_on = 'cs.modified')
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
    ecbl_mfn = Property(depends_on = 'cs.modified')
    @cached_property
    def _get_ecbl_mfn(self):
        self.ecbl.set_cparams(*self.u_sol)
        return self.ecbl.mfn

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
    
    ec = ECBLCalib(# mean concrete strength after 9 days
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
                       )

    do = 'plot_ecbl'

    if do == 'debug':
        
        print ec.u0[0]
        ec.cs.set(eps_lo = ec.u0[0], eps_up = -ec.eps_cu)   
        eps_tex_u = ec.cs.convert_eps_lo_2_tex_u(ec.u0[0])
        print 'eps_tex_u', eps_tex_u
        ec.ecbl.set_cparams(eps_tex_u, ec.u0[1])
        
        print 'N'
        print ec.cs.N

        print 'M'
        print ec.cs.M
        
        ec.ecbl_mfn
                        
    elif do == 'plot_ecbl':
        for sig_tex_u, color in zip([1200, 1300, 1400], ['red', 'green', 'blue', 'black', 'orange', 'brown']):
        #for sig_tex_u, color in zip([1216], ['red']):

            #for ecbl_type in ['linear', 'cubic', 'fbm']:
            for ecbl_type in ['cubic']:
                print 'CALIB TYPE', ecbl_type
                ec.n = 0
                ec.cs.ecbl_type = ecbl_type
                ec.ecbl.sig_tex_u = sig_tex_u
                ec.get_lack_of_fit(ec.u0)
                ec.ecbl_mfn.plot(p, color = color, linewidth = 8)
                print 'E_yarn', ec.ecbl_mfn.get_diff(0.00001)
                print 'INTEG', ec.ecbl_mfn.integ_value

#            ec.ecbl_type = 'bilinear'
#            ec.ecbl.sig_tex_u = sig_tex_u
#            for eps_el_fraction in np.linspace(0.25, 0.99999, 4):
#                ec.n = 0
#                ec.ecbl.eps_el_fraction = eps_el_fraction
#                ec.ecbl_mfn.plot(p, color = color)
#                print 'E_yarn', ec.ecbl_mfn.get_diff(0.00001)
#                print 'INTEG', ec.ecbl_mfn.integ_value
        p.plot([0.0, 0.01], [0.0, 2400], color = 'black')
            
    elif do == 'plot_cs_state':    
    #    ### plot crack bridge law [MPa]:
    #    #
    #    p.subplot(1, 2, 1)
    #    sig_fl_calib.plot_ecb_law(u_sol)
    #
        ### plot concrete law
        ec.plot_sig_c_mfn()
    #    p.subplot(1, 2, 2)
        ec.plot_sig_comp_i()
    #
    elif do == 'plot_MN_grid':
    
        u_grid = np.ogrid[0.01:0.02:30j, 0.1:2.0:30j]
        N_grid = ec.get_lack_of_fit_dN(u_grid[0], u_grid[1])
        M_grid = ec.get_lack_of_fit_dM(u_grid[0], u_grid[1])
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
        
