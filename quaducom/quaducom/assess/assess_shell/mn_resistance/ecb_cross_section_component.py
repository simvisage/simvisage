'''
Created on Sep 4, 2012

@todo: introduce the dock feature for the views
@todo: classify the state changes and provide examples.

@author: rch
'''
from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, \
    Event, on_trait_change, DelegatesTo, Instance, WeakRef, Constant

from etsproxy.traits.ui.api import \
    View, Item, Group, HGroup

from ecb_cross_section_state import \
    ECBCrossSectionState

import numpy as np


ECB_COMPONENT_CHANGE = '+geo_input,+law_input'
ECB_COMPONENT_AND_EPS_CHANGE = 'eps_changed,+geo_input,+law_input'

class ECBCrossSectionComponent(HasStrictTraits):
    '''Cross section component supplying the normal force and moment..
    '''
    state = WeakRef(ECBCrossSectionState)
    '''Strain state of a cross section
    '''

    unit_conversion_factor = Constant(1000.0)

    eps_changed = Event
    '''State notifier that is set to true if the cross section state has changed
    upon modifications of eps_lo and eps_up
    '''

    @on_trait_change(ECB_COMPONENT_CHANGE)
    def notify_change(self):
        '''Propagate the change of the component geometry to the cross section state.
        '''
        if self.state:
            self.state.changed = True


    #===========================================================================
    # Cross-sectional stress resultants
    #===========================================================================

    N = Property(depends_on=ECB_COMPONENT_AND_EPS_CHANGE)
    '''Get the resulting normal force.
    '''
    @cached_property
    def _get_N(self):
        return sum(self.f_ti_arr)

    M = Property(depends_on=ECB_COMPONENT_AND_EPS_CHANGE)
    '''Get the resulting moment evaluated with respect to the center line
    '''
    @cached_property
    def _get_M(self):
        return np.dot(self.f_ti_arr, self.z_ti_arr)

if __name__ == '__main__':
    #ecs.configure_traits()
    pass
