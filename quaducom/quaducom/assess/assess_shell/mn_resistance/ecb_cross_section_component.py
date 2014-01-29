'''
Created on Sep 4, 2012

@todo: introduce the dock feature for the views
@todo: classify the state changes and provide examples.

@author: rch
'''
from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, \
    Event, on_trait_change, Callable, Instance, WeakRef, Constant

from etsproxy.traits.ui.api import \
    View, Item, Group, HGroup

from ecb_cross_section_state import \
    ECBCrossSectionState

import numpy as np


ECB_COMPONENT_DEPEND = '+eps_input,+geo_input,+law_input'

class ECBCrossSectionComponent(HasStrictTraits):
    '''Cross section component supplying the normal force and moment..
    '''
    state = WeakRef(ECBCrossSectionState)
    '''Strain state of a cross section
    '''

    unit_conversion_factor = Constant(1000.0)

    #===========================================================================
    # State management
    #===========================================================================
    modified = Event
    @on_trait_change(ECB_COMPONENT_DEPEND)
    def set_modified(self):
        self.modified = True


    #===========================================================================
    # Cross-sectional stress resultants
    #===========================================================================

    N = Property(depends_on=ECB_COMPONENT_DEPEND)
    '''Get the resulting normal force.
    '''
    @cached_property
    def _get_N(self):
        return sum(self.f_ti_arr)

    M = Property(depends_on=ECB_COMPONENT_DEPEND)
    '''Get the resulting moment evaluated with respect to the center line
    '''
    @cached_property
    def _get_M(self):
        return np.dot(self.f_ti_arr, self.z_ti_arr)

if __name__ == '__main__':
    #ecs.configure_traits()
    pass
