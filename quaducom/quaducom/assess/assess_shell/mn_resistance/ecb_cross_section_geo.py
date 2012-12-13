'''
Created on Sep 4, 2012

@todo: introduce the dock feature for the views
@todo: classify the state changes and provide examples.


@author: rch
'''
from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, \
    Trait, Event, on_trait_change, Instance, Button, Callable, \
    WeakRef, Dict

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, VGroup, HGroup

from ecb_law import \
    ECBLBase, ECBLLinear, ECBLFBM, ECBLCubic, ECBLBilinear, ECBLPiecewiseLinear

from constitutive_law import \
    ConstitutiveLawModelView

from cc_law import \
    CCLawBase, CCLawBlock, CCLawLinear, CCLawQuadratic, CCLawBilinear, CCLawQuad

import numpy as np

class ECBCrossSectionGeo(HasStrictTraits):

    #---------------------------------------------------------------
    # Cross section characteristics needed for tensile specimens 
    #---------------------------------------------------------------

    # thickness of reinforced cross section
    #
    thickness = Float(0.06, auto_set = False, enter_set = True, geo_input = True)

    # total number of reinforcement layers [-]
    # 
    n_layers = Int(12, auto_set = False, enter_set = True, geo_input = True)

    #---------------------------------------------------------------
    # Cross section characteristics needed for bending specimens 
    #---------------------------------------------------------------

    # width of the cross section [m]
    #
    width = Float(0.20, auto_set = False, enter_set = True, geo_input = True)

    # number of rovings in 0-direction of one composite 
    # layer of the bending test [-]:
    #
    n_rovings = Int(23, auto_set = False, enter_set = True, geo_input = True)

    # cross section of one roving [mm**2]:
    #
    A_roving = Float(0.461, auto_set = False, enter_set = True, geo_input = True)

    notify_change = Callable
    #===========================================================================
    # State management
    #===========================================================================
    modified = Event
    @on_trait_change('+geo_input,+cc_input,+tt_input')
    def set_modified(self):
        self.modified = True
        if self.notify_change:
            self.notify_change()

    #===========================================================================
    # Distribution of reinforcement
    #===========================================================================

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
    n_cj = Int(20, auto_set = False, enter_set = True,
                 cc_input = True)

    #===========================================================================
    # Compressive concrete constitutive law
    #===========================================================================
    cc_law_type = Trait('constant', dict(constant = CCLawBlock,
                                         linear = CCLawLinear,
                                         bilinear = CCLawBilinear,
                                         quadratic = CCLawQuadratic,
                                         quad = CCLawQuad),
                        cc_input = True)

    cc_law_params = Dict

    cc_law = Property(Instance(CCLawBase), depends_on = '+cc_input')
    @cached_property
    def _get_cc_law(self):
        '''Construct the compressive concrete law'''
        kw = {'cs' : self }
        kw.update(self.cc_law_params[self.cc_law_type])
        return self.cc_law_type_(**kw)

    show_cc_law = Button
    def _show_cc_law_fired(self):
        cc_law_mw = ConstitutiveLawModelView(model = self.cc_law)
        cc_law_mw.edit_traits(kind = 'live')
        return

    cc_modified = Event

    #===========================================================================
    # Effective crack bridge law
    #===========================================================================
    ecb_law_type = Trait('fbm', dict(fbm = ECBLFBM,
                                  cubic = ECBLCubic,
                                  linear = ECBLLinear,
                                  bilinear = ECBLBilinear,
                                  piecewise_linear = ECBLPiecewiseLinear),
                      tt_input = True)

    ecb_law = Property(Instance(ECBLBase), depends_on = '+tt_input')
    @cached_property
    def _get_ecb_law(self):
        return self.ecb_law_type_(cs = self)

    show_ecb_law = Button
    def _show_ecb_law_fired(self):
        ecb_law_mw = ConstitutiveLawModelView(model = self.ecb_law)
        ecb_law_mw.edit_traits(kind = 'live')
        return

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor = 'white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    data_changed = Event

    replot = Button
    def _replot_fired(self):

        self.figure.clear()
        ax = self.figure.add_subplot(1, 2, 1)
        self.cc_law.plot(ax)

        ax = self.figure.add_subplot(1, 2, 2)
        self.ecb_law.plot(ax)

        self.data_changed = True

    view = View(HSplit(Group(
                HGroup(
                Group(Item('thickness', springy = True),
                      Item('width'),
                      Item('n_layers'),
                      Item('n_rovings'),
                      Item('A_roving'),
                      label = 'Geometry',
                      springy = True
                      ),
                springy = True,
                ),
                HGroup(
                Group(VGroup(
                      Item('cc_law_type', show_label = False, springy = True),
                      Item('cc_law', label = 'Edit', show_label = False, springy = True),
                      Item('show_cc_law', label = 'Show', show_label = False, springy = True),
                      springy = True
                      ),
                      Item('f_ck', label = 'Compressive strength'),
                      Item('n_cj', label = 'Discretization'),
                      label = 'Concrete',
                      springy = True
                      ),
                Group(VGroup(
                      Item('ecb_law_type', show_label = False, springy = True),
                      Item('ecb_law', label = 'Edit', show_label = False, springy = True),
                      Item('show_ecb_law', label = 'Show', show_label = False, springy = True),
                      springy = True,
                      ),
                      label = 'Reinforcement',
                      springy = True
                      ),
                springy = True,
                ),
                Group(Item('s_tex_z', label = 'vertical spacing', style = 'readonly'),
                      label = 'Layout',
                      ),
                scrollable = True,
                ),
                Group(Item('replot', show_label = False),
                      Item('figure', editor = MPLFigureEditor(),
                           resizable = True, show_label = False),
                      id = 'simexdb.plot_sheet',
                      label = 'plot sheet',
                      dock = 'tab',
                      ),
                       ),
                width = 0.8,
                height = 0.7,
                resizable = True,
                buttons = ['OK', 'Cancel'])

if __name__ == '__main__':
    ecs = ECBCrossSectionGeo(# 7d: f_ck,cube = 62 MPa; f_ck,cyl = 62/1.2=52
                           # 9d: f_ck,cube = 66.8 MPa; f_ck,cyl = 55,7
                           f_ck = 55.7,

                           ecb_law_type = 'fbm',
                           cc_law_type = 'quadratic'
                           )

    ecs.configure_traits()
