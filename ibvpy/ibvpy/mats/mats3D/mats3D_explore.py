
from etsproxy.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasStrictTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate, List, WeakRef

from util.traits.either_type import \
    EitherType

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import MATS3DMicroplaneDamage

from ibvpy.mats.matsXD.matsXD_explore import MATSXDExplore

class MATS3DExplore(MATSXDExplore):
    '''
    Simulate the loading histories of a material point in 2D space.
    '''

    mats_eval = EitherType(klasses=[MATS3DElastic,
                                       MATS3DMicroplaneDamage,
                                        ])

    def _mats_eval_default(self):
        return MATS3DElastic()
