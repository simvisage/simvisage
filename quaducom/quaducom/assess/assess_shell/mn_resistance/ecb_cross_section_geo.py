'''
Created on Sep 4, 2012

@todo: introduce the dock feature for the views
@todo: classify the state changes and provide examples.

@author: rch
'''
from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, \
    Event, on_trait_change, Callable

from etsproxy.traits.ui.api import \
    View, Item, Group, HGroup

import numpy as np

class ECBCrossSectionGeo(HasStrictTraits):
    '''Cross section characteristics needed for tensile specimens.
    '''

    thickness = Float(0.06, auto_set=False, enter_set=True, geo_input=True)
    '''thickness of reinforced cross section
    '''

    n_layers = Int(12, auto_set=False, enter_set=True, geo_input=True)
    '''total number of reinforcement layers [-]
    '''

    #---------------------------------------------------------------
    # Cross section characteristics needed for bending specimens 
    #---------------------------------------------------------------

    width = Float(0.20, auto_set=False, enter_set=True, geo_input=True)
    '''width of the cross section [m]
    '''

    n_rovings = Int(23, auto_set=False, enter_set=True, geo_input=True)
    '''number of rovings in 0-direction of one composite layer of the
    bending test [-]:
    '''

    A_roving = Float(0.461, auto_set=False, enter_set=True, geo_input=True)
    '''cross section of one roving [mm**2]'''

    #===========================================================================
    # Derived geometric data
    #===========================================================================

    s_tex_z = Property(depends_on='+geo_input')
    '''property: spacing between the layers [m]
    '''
    @cached_property
    def _get_s_tex_z(self):
        return self.thickness / (self.n_layers + 1)

    z_ti_arr = Property(depends_on='+geo_input')
    '''property: distance from the top of each reinforcement layer [m]:
    '''
    @cached_property
    def _get_z_ti_arr(self):
        return np.array([ self.thickness - (i + 1) * self.s_tex_z
                         for i in range(self.n_layers) ],
                      dtype=float)

    zz_ti_arr = Property
    '''property: distance of reinforcement layers from the bottom
    '''
    def _get_zz_ti_arr(self):
        return self.thickness - self.z_ti_arr

    modified = Event
    @on_trait_change('+geo_input')
    def set_modified(self):
        self.modified = True

    view = View(HGroup(
                Group(Item('thickness', springy=True),
                      Item('width'),
                      Item('n_layers'),
                      Item('n_rovings'),
                      Item('A_roving'),
                      label='Geometry',
                      springy=True
                      ),
                springy=True,
                ),
                resizable=True,
                buttons=['OK', 'Cancel'])

if __name__ == '__main__':
    ecs = ECBCrossSectionGeo(n_layers=3, thickness=1)
    print 'zz_ti_arr', ecs.zz_ti_arr
    ecs.configure_traits()
