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
    breaking strain of the matrix:    eps_m_u [-]
    
    E-Modulus matrix:                 Em [N/mm^2]                                
    E-Modulus fibers:                 Ef [N/mm^2]
    reinforcement ratio:              Vf [-]
    
    program parameters :
    datapoints:                       nd
    plot range:                       sigma_plot [N/mm^2]
    '''
    
    nd = Int( 10 )
    sigma_plot = Float
    sigma_m_u = Float
    
    eps_m_u = Property( Float )
    def _get_eps_m_u( self ):
        return self.sigma_m_u / self.Em
    
    sigma_arr = Property( Array, depends_on = 'sigma,nd' )
    def _get_sigma_arr( self ):
        return np.linspace( 1e-10, sigma_plot, self.nd )

        
    def eps_c( self, Em, Ef, Vf, sigma ):
        Kc = ( ( 1 - Vf ) * Em + Vf * Ef )
        return sigma / Kc 
    
    def eps_f( self, Ef, Vf, sigma ):
        sigma_f = sigma / Vf
        Kf = Ef 
        return sigma_f / Kf 
    
    def d( self ):
        '''evaluates the strain and gives sigma arr'''
        
        #tension stiffening component
        alpha = self.Em * self.Vm / ( self.Ef * self.Vf )
        eps_p_b = self.eps_m_u * ( 1 + alpha / 2 )
        eps_p_c = self.eps_f( self.Ef, self.Vf, self.eps_m_u * self.Ec )
        eps_diff = eps_p_c - eps_p_b 

        #strain evaluation
        eps_c = self.eps_c( self.Em, self.Ef, self.Vf, self.sigma_arr )  
        eps_f = self.eps_f( self.Ef, self.Vf, self.sigma_arr ) - eps_diff
        eps_c_arr = eps_c * H( self.eps_m_u * self.Ec - self.sigma_arr ) + eps_f * H( self.sigma_arr - self.eps_m_u * self.Ec )
        
        #smoothing by inserting new datapoint
        index_l = list( H( self.eps_m_u * self.Ec - self.sigma_arr ) )
        ix = index_l.index( 0 )
        sigma_list = list( self.sigma_arr )
        sigma_list.insert( ix, sigma_list[ix - 1] )
        eps_c_list = list( eps_c_arr )
        eps_c_list.insert( ix , eps_f[ix - 1] )
        return  eps_c_list , sigma_list 
        
if __name__ == '__main__':
    from matplotlib import pyplot as plt
    
    '''must be set'''
    sigma_plot = 4.
    sigma_m_u = 3.
   
    '''parameters that can be set'''
    Vf = 0.02
    #Em=
    #Ef=
    #nd
    a = ACK( sigma_plot = sigma_plot, sigma_m_u = sigma_m_u, Vf = Vf )
    
    def dplot():
        eps, sigma = a.d()
        plt.plot( eps, sigma )
        plt.ylabel( '$\sigma$ in [MPa]', fontsize = 16 )
        plt.xlabel( '$\epsilon$ in [-]', fontsize = 16 )
        plt.title( 'ACK-Model ' )
        plt.show()

    dplot()
 
    
