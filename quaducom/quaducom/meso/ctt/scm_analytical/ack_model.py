'''
Created on Apr 20, 2012

@author: rostar
'''

from material import Material
from etsproxy.traits.api import Array, Property


class ACK(Material):
    
    sigma = Array
    eps = Property()
    def _get_eps(self):
        return 1.0