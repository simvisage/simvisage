'''
Created on Oct 29, 2012

@author: rch
'''

from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, \
    Trait, Event, on_trait_change, Instance, Button, Callable

from ecb_cross_section import ECBCrossSection

class ECBCrossSectionHistory(ECBCrossSection):
    '''
    '''

