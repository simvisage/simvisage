'''
Created on Sep 4, 2012

@author: rch
'''
from traits.api import \
    HasTraits, Float, Property, cached_property, Int, Array
    
import numpy as np

class ECBCrossSection(HasTraits):

    #---------------------------------------------------------------
    # Cross section characteristics needed for tensile specimens 
    #---------------------------------------------------------------
    
    # thickness of reinforced cross section
    #
    thickness = Float(0.06, geo_input = True)

    # total number of reinforcement layers [-]
    # 
    n_layers = Int(12, geo_input = True)

    #---------------------------------------------------------------
    # Cross section characteristics needed for bending specimens 
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

    # distance of reinforcement layers from the bottom 
    #
    zz_ti_arr = Property
    def _get_zz_ti_arr(self):
        return self.thickness - self.z_ti_arr
    
    # number of subdivisions of the compressive zone
    #
    n_cj = Float(20., sig_c_modified = True)

    zeta_cj_arr = Property(Array)
    @cached_property
    def _get_zeta_cj_arr(self):
        '''subdivide the compression zone 'x' in 'n_cj' sub-areas;
        'zeta_cj_arr' giving the fraction of each distance of the sub-area from the upper surface 
        with respect to the compressive zone 'x'
        '''
        return np.arange(self.n_cj) / self.n_cj + 1. / (2. * self.n_cj)
    
if __name__ == '__main__':
    ecs = ECBCrossSection()
    ecs.configure_traits()
    print ecs.z_ti_arr
