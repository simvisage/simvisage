#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Oct 11, 2011 by: rch

from enthought.traits.api import HasTraits, Float, Property, cached_property, \
                                Instance, List, on_trait_change, Int, Tuple, Bool, \
                                Event, Button, Str
import os
import pylab as p

from numpy import \
    linspace, frompyfunc

from stats.spirrid.rf import RF
from brittle_filament import Filament
from cb_clamped_fiber import CBClampedFiber
from cb_infinite_fiber import CBInfiniteFiber
from cb_short_fiber import CBShortFiber
from po_clamped_fiber import POClampedFiber
from po_infinite_fiber import POInfiniteFiber
from po_short_fiber import POShortFiber

class WikiGen(HasTraits):

    rf_list = List([Filament,
                     CBClampedFiber,
                     CBInfiniteFiber,
                     CBShortFiber,
                     POClampedFiber,
                     POInfiniteFiber,
                     POShortFiber])

    max_x = Float(0.01, enter_set = True, auto_set = False, config_change = True)
    n_points = Int(200, enter_set = True, auto_set = False, config_change = True)

    def export_wiki(self):

        fname = os.path.join('wiki', 'rf.wiki')
        f = open(fname, 'w')

        for rf_class in self.rf_list:
            q = rf_class()

            f.write(str(q) + '\n')
            qname = q.__class__.__name__

            p.figure()
            q.plot(p, linewidth = 2, color = 'navy')

            fig_fname = os.path.join('wiki', qname + '.png')
            p.savefig(fig_fname)

if __name__ == '__main__':
    rf = WikiGen()
    rf.export_wiki()
