#!/usr/bin/env python
""" The entry point for an Envisage application. """

# Standard library imports.
import sys
import os.path
import logging

# Enthought library imports.
from etsproxy.mayavi.plugins.app import get_plugins, setup_logger
from etsproxy.traits.api import List
from etsproxy.envisage.api import Plugin, ServiceOffer
from etsproxy.envisage.ui.workbench.api import WorkbenchApplication
from etsproxy.pyface.workbench.api import Perspective, PerspectiveItem

###############################################################################
# `IBVPYPlugin` class.
###############################################################################
class TStepperPlugin(Plugin):

    # Extension points we contribute to.
    SERVICE_OFFERS = 'enthought.envisage.ui.workbench.service_offers'

    # The plugin's unique identifier.
    id = 'TStepper.TStepper'

    # The plugin's name (suitable for displaying to the user).
    name = 'IBVPY'

    # Services we contribute.
    service_offers = List(contributes_to = SERVICE_OFFERS)
    
    ######################################################################
    # Private methods.
    def _service_offers_default(self):
        """ Trait initializer. """
        tstepper_service_offer = ServiceOffer(
            protocol = 'ibvpy.plugins.tstepper_service.TStepperService',
            factory = 'ibvpy.plugins.tstepper_service.TStepperService'
        )

        return [tstepper_service_offer]
