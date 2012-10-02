'''
Created on Apr 20, 2012

@author: rostar
'''
from material import Material
from etsproxy.traits.api import Float, Property, cached_property, Array, Int
from scipy.stats import weibull_min
import numpy as np
from math import e
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from scipy.integrate import odeint
import pickle

def H(self, x):
    return x > 0

class ACK(Material):
    '''blah eps_c, eps_f, eps_m, sigma_c, sigma_m, sigma_f'''
    
    nd = Int(100)
    max_sigma = Float
    
    sigma_arr = Property(Array, depends_on='sigma,nd')
    def _get_sigma_arr(self):
        return np.linspace(1e-10, sigma, self.nd)
    
    

    
    #Ack functions
  
        
    def w_c(self, Em, Ef, Vf, sigma, Lc):
        Kc = ((1 - Vf) * Em + Vf * Ef) 
        return sigma / Kc * Lc
    
    def w_f(self, Em, Ef, Vf, sigma, Lc):
        sigma_f = sigma / Vf
        Kf = Ef 
        return sigma_f / Kf * Lc
    
    def w(self):
        r_c = self.w_c(self.Em, self.Ef, self.Vf, self.sigma_arr, self.Lc)
        r_f = self.w_f(self.Em, self.Ef, self.Vf, self.sigma_arr, self.Lc)     
        return r_c * self.H(sigma_0 - self.sigma_arr) + r_f * self.H(self.sigma_arr - sigma_0)
    
    def d(self):
        r_c = self.w_c(self.Em, self.Ef, self.Vf, self.sigma_arr, self.Lc) / self.Lc * 1e3
        r_f = self.w_f(self.Em, self.Ef, self.Vf, self.sigma_arr, self.Lc) / self.Lc * 1e3
        return r_c * self.H(sigma_0 - self.sigma_arr) + r_f * self.H(self.sigma_arr - sigma_0)
      
if __name__ == '__main__':
    from matplotlib import pyplot as plt
    sigma = 10.
    sigma_0 = 5.
    a = ACK(max_sigma=sigma)
    #Em=
    #Ef=
    #Vf=
    #sigma=
    #Lc=
    #sigma_0=
    def wplot():
        plt.plot(a.w(), a.sigma_arr)
        plt.ylabel('sigma in [MPa]', fontsize=16)
        plt.xlabel('displacement in [mm]', fontsize=16)
        plt.title('ACK-Model ')
        plt.show()
        
    def dplot():
        plt.plot(a.w(), a.sigma_arr)
        plt.ylabel('sigma in [MPa]', fontsize=16)
        plt.xlabel(' $\epsilon$  in [10e-3]', fontsize=16)
        plt.title('ACK-Model')
        plt.show()
        
    #wplot()
    dplot()
    
