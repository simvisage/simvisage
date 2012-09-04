'''
Created on Sep 4, 2012

@author: rch
'''
from traits.api import \
    HasTraits, Float, Property, cached_property, Int
    
import numpy as np

class ECBCrossSection(HasTraits):
    
    # thickness of reinforced cross section
    #
    thickness = Float(0.06, geo_input = True)

    # total number of reinforcement layers [-]
    # 
    n_layers = Int(12, geo_input = True)

    # spacing between the layers [m]
    #
    s_tex_z = Property(depends_on = '+geo_input')
    @cached_property
    def _get_s_tex_z(self):
        return self.thickness / (self.n_layers + 1)

    # distance from the top of each reinforcement layer [m]:
    #
    z_ti_arr = Property(depends_on = '+geo_input')
    @cached_property
    def _get_z_ti_arr(self):
        return np.array([ self.thickness - (i + 1) * self.s_tex_z for i in range(self.n_layers) ],
                      dtype = float)

    zz_ti_arr = Property
    def _get_zz_ti_arr(self):
        return self.thickness - self.z_ti_arr
    #---------------------------------------------------------------
    # properties of the bending specimens (tree point bending test):
    #---------------------------------------------------------------

    # width of the cross section [m]
    #
    width = Float(0.20, geo_input = True)
    
    # number of rovings in 0-direction of one composite 
    # layer of the bending test [-]:
    #
    n_rovings = Int(23, geo_input = True)
    
    # cross section of one roving [mm**2]:
    #
    A_roving = Float(0.461, geo_input = True)

if __name__ == '__main__':
    ecs = ECBCrossSection()
    print ecs.z_t_i_arr
