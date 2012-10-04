from material import Material
from etsproxy.traits.api import Float, Property, cached_property, Array, Int
import numpy as np
import array
def H( x ):
    return x >= 0

class ACK( Material ):
    '''
    ACK model with one matrix breaking strain all over the specimen
     
    material parameters:
    tension strength matrix:          sigma_mu[MPa]
    
    E-Modulus matrix:                 Em [MPa]                                
    E-Modulus fibers:                 Ef [MPa]
    reinforcement ratio:              Vf [-]
    
    program parameters :
    plot range:                       sigma_max [N/mm^2]
    '''
    
    sigma_max = Float
    sigma_mu = Float
    
    def eps_1( self, sigma ):
        #sigma<sigma_mu
        Kc = ( ( 1 - self.Vf ) * self.Em + self.Vf * self.Ef )
        eps_m_u = self.sigma_mu / Kc
        return eps_m_u

    def eps_2( self, sigma ):
        #sigma=sigma_mu
        alpha = self.Em * self.Vm / ( self.Ef * self.Vf )
        eps_m_c = ( 1 + 0.666 * alpha ) * self.sigma_mu / self.Em
        return eps_m_c
    
    def eps_3( self, sigma ):
        #sigma>sigma_mu
        eps_c_u = sigma / self.Ef / self.Vf - self.sigma_mu / self.Ef / self.Vf + self.eps_2( self.sigma_mu ) 
        return eps_c_u
        
  
    
    def plot_diagram( self ):
        eps_list = [0, self.eps_1( self.sigma_mu ), self.eps_2( self.sigma_mu ), self.eps_3( self.sigma_max )]
        sigma_list = [0, self.sigma_mu, self.sigma_mu, self.sigma_max]
        plt.plot( eps_list, sigma_list )
        plt.ylabel( '$\sigma$ in [MPa]', fontsize = 16 )
        plt.xlabel( '$\epsilon$ in [-]', fontsize = 16 )
        plt.title( 'ACK-Model ' )
        plt.show()
        
if __name__ == '__main__':
    from matplotlib import pyplot as plt
    a = ACK( sigma_max = 10.,
            V_f = 0.02,
            E_f = 200e3,
            E_m = 30e3,
            sigma_mu = 5.0 )
    
    a.plot_diagram()

 
    
