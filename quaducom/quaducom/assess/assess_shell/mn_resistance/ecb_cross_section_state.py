'''
Created on Sep 4, 2012

@todo: introduce the dock feature for the views
@todo: classify the state changes and provide examples.


@author: rch
'''
from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property

class ECBCrossSectionState(HasStrictTraits):
    '''
    Cross section state is defined by the linear profile of strains
    with eps_up and eps_lo at the top and at the bottom of the cross section,
    respectively.
    '''

    eps_up = Float(-0.0033, auto_set=False, enter_set=True, eps_input=True)
    eps_lo = Float(0.0140, auto_set=False, enter_set=True, eps_input=True)

