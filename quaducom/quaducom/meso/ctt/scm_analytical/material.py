'''
Created on Apr 18, 2012

@author: rostar
'''
from etsproxy.traits.api import HasTraits, Float, Property, \
                                cached_property, Array
from etsproxy.traits.ui.api import View, Item, Tabbed, VGroup, \
                                VSplit, Group
from etsproxy.traits.ui.menu import OKButton, CancelButton

import numpy as np
from scipy.stats import weibull_min

from math import cos, pi as Pi, e

class Material( HasTraits ):
    '''
    Material parameters
    '''
    Ef = Float( 72e+3, auto_set = False, enter_set = True, # [N/mm^2]
                 desc = 'Modulus of elasticity of the fiber [MPa]',
                 modified = True )

    Em = Float( 30e+3, auto_set = False, enter_set = True, # [N/mm^2]
                 desc = 'Modulus of elasticity the matrix [MPa]',
                 modified = True )

    tau = Float( .1, auto_set = False, enter_set = True, # [N/mm^2]
                 desc = 'frictional stress between fiber and matrix [MPa]',
                 modified = True )

    r = Float( 5e-4, auto_set = False, enter_set = True, # [mm]
                 desc = 'radius of the fiber',
                 modified = True )

    m = Float( 5.3, auto_set = False, enter_set = True, # [-]
                 desc = 'Weibull shape parameter for the matrix tensile strength [-]',
                 modified = True )

    sigma_0 = Float( 6.0, auto_set = False, enter_set = True, # [N/mm^2]
                       desc = 'Weibull scale parameter for the matrix tensile strength [MPa]',
                       modified = True )

    Vf = Float(0.0175, auto_set = False, enter_set = True, # [-]
                       desc = 'reinforcement ratio [-]',
                       modified = True )

    l0 = Float(10., auto_set = False, enter_set = True, # [-]
                   desc = 'reference length for the matrix strength distribution [mm]',
                   modified = True )
    
    Lc = Float(1000., auto_set = False, enter_set = True, # [-]
                   desc = 'specimen length [mm]',
                   modified = True )

    sigma = Float(auto_set = False, enter_set = True, # [-]
                   desc = 'load [MPa]',
                   modified = True )

    Vm = Property( Float, depends_on = 'Vf' )
    @cached_property
    def _get_Vm( self ):
        return 1 - self.Vf

    Ec = Property( Float, depends_on = 'Em, Ef, Vf' )
    @cached_property
    def _get_Ec( self ):
        return self.Ef * self.Vf + self.Em * (1.-self.Vf)

    def sigma_mu_distr(self, L):
        return weibull_min(self.m, scale = self.sigma_0*(L/self.l0)**(-1./self.m))
    
if __name__ == '__main__':
    from matplotlib import pyplot as plt
    m = Material()
    x = np.linspace(0,20,200)
    y = m.sigma_mu_distr(20).cdf(x)
    plt.plot(x,y)
    plt.show()