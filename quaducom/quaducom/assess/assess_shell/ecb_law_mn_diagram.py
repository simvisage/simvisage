'''
Created on Sep 4, 2012

@author: rch
'''
from etsproxy.traits.api import \
    HasTraits, Float, Instance, Property, cached_property, DelegatesTo

from ecb_law_calib import \
    ECBLCalib
    
import numpy as np
    
class ECBLMNDiagram(HasTraits):

    # calibrator supplying the effective material law
    calib = Instance(ECBLCalib)
    
    # cross section
    cs = DelegatesTo('calib')

    def _get_eps_f(self, eps_lo, eps_up):
        '''EVALUATION for bending case: using the calibrated constitutive law of the layer
        '''
        # @todo: u = 'u_sol' is not needed for the evaluation of given stress at the top and bottom but must be
        # passed in general so that if calibration method is called 'layer_response' works correctly
        # use option u_sol = () instead of argument  

        print 'eval_layer_response_f called'

        eps_t = eps_lo  
        eps_c = abs(eps_up) 

        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_t'
        # ------------------------------------------------------------------------                

        thickness = self.thickness
        z_t_i_arr = self.cs.z_t_i_arr
       
        # ------------------------------------------------------------------------                
        # derived params depending on value for 'eps_t'
        # ------------------------------------------------------------------------                

        # heights of the compressive zone:
        #
        x = abs(eps_c) / (abs(eps_c) + abs(eps_t)) * thickness

        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = abs(eps_t) / (thickness - x) * (z_t_i_arr - x)

        # use a ramp function to consider only negative strains
        # NOTE: used only for plot at the height of the layers
        #
        eps_c_i_arr = abs(-np.fabs(eps_i_arr) + eps_i_arr) / 2.0

        # use a ramp function to consider only positive strains
        #
        eps_t_i_arr = (np.fabs(eps_i_arr) + eps_i_arr) / 2.0

        sig_t_mfn = self.sig_t_mfn

        return x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn, eps_t, eps_c

    def _get_eps_t(self, eps_lo, eps_up):
        '''EVALUATION for tension case: using the calibrated constitutive law of the layer
        ''' 
        print 'eval_layer_response_t called'

        eps_t_lo = abs(eps_lo)  
        eps_t_up = abs(eps_up) 
        
        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_tl'
        # ------------------------------------------------------------------------                

        thickness = self.thickness
        z_t_i_arr = self.z_t_i_arr
   
        # ------------------------------------------------------------------------                
        # derived params depending on value for 'eps_tl'
        # ------------------------------------------------------------------------                

        # heights of the compressive zone: the whole cross-section is stressed with tension            
        #
        x = 0.

        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = eps_t_lo - (eps_t_lo - eps_t_up) / thickness * (thickness - z_t_i_arr)
        eps_t_i_arr = eps_i_arr
        eps_c_i_arr = abs(-np.fabs(eps_i_arr) + eps_i_arr) / 2.0
        
        sig_t_mfn = self.sig_t_mfn

        return x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn, eps_t_lo, eps_t_up

    def _get_eps_c(self, eps_lo, eps_up):
        '''EVALUATION for compression case: using the calibrated constitutive law of the layer
        ''' 

        print 'eval_layer_response_c called'
        eps_c_lo = abs(eps_lo)  
        eps_c_up = abs(eps_up) 
        
        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_cl'
        # ------------------------------------------------------------------------                

        thickness = self.thickness
        z_t_i_arr = self.z_t_i_arr
       
        # ------------------------------------------------------------------------                
        # derived params depending on value for 'eps_cl'
        # ------------------------------------------------------------------------                

        # heights of the compressive zone: the whole cross-section is stressed by compression 
        #
        x = thickness
        
        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = -(abs(eps_c_lo) - (abs(eps_c_lo) - abs(eps_c_up)) / thickness * (x - z_t_i_arr))
        eps_t_i_arr = (np.fabs(eps_i_arr) + eps_i_arr) / 2.0
        eps_c_i_arr = abs(eps_i_arr)
        sig_t_mfn = self.sig_t_mfn

        return x, eps_t_i_arr, eps_c_i_arr, sig_t_mfn, eps_c_lo, eps_c_up

    def _get_eps_ij(self, eps_lo, eps_up):
        '''determine_stress_case
        '''
        if eps_up < 0 and eps_lo > 0:
            return self._get_eps_f(eps_lo, eps_up)

        elif self.eps_up <= 0 and self.eps_lo <= 0:
            return self._get_eps_c(eps_lo, eps_up)
          
        elif self.eps_up >= 0 and self.eps_lo >= 0:
            return self._get_eps_t(eps_lo, eps_up)
            
    def get_M_N_value(self, eps_lo, eps_up):
        '''evaluate the normal force and bending moment for given strains at the top and bottom
        using a calibrated crack bridge law;
        '''
        print "'eval_N_M' called"
        # reset iteration counter
#        self.n = 0
        self.eps_lo = eps_lo
        self.eps_up = eps_up
        self.calc_mode = 'eval'
        self.stress_case = self.determine_stress_case()
        u_sol = [ abs(self.eps_lo), abs(self.eps_up) ]

        # use method 'lack_of_fit' setting external forces to 0.
        # yielding 'N_internal = -N_ck + N_tk' and 'M_internal'
        # NOTE: the returned internal forces would correspond to an 
        # external loading yielding the given strain distribution
        #
        N_internal, M_internal = self.get_lack_of_fit(u_sol, 0., 0.)
        return N_internal, M_internal

    N_M_diagram = Property(depends_on = 'calib.+config_changed')
    @cached_property
    def _get_N_M_diagram(self, eps_lo, eps_up):
        pass

if __name__ == '__main__':
    c = ECBLCalib(f_ck = 55.7,
                  eps_cu = 3.3 / 1000.,
                  Mu = 3.49,
                  width = 0.20,
                  n_roving = 23.,
                  ecbl_type = 'linear',
                  sig_c_config = 'quadratic'             #eps_tu 0.0137279096658                              
                  )
                     
    ECBLMNDiagram(calib = c)
