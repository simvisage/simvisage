'''
Created on Aug 23, 2012

@author: rch
'''

from etsproxy.traits.api import \
    HasTraits, List, Float
    
import numpy as np

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
        sig_tex_u = self.sig_tex_u
        eps_tex_arr = np.hstack([0., 0.999999 * eps_tex_u, eps_tex_u ])
        sig_tex_arr = np.hstack([0., var_a * sig_tex_u, sig_tex_u])
        return eps_tex_arr, sig_tex_arr         

if __name__ == '__main__':
    ecbl = ECBLPlastic(sig_tex_u = 2e+3)
    print ecbl(0.04, 0.02)
    
