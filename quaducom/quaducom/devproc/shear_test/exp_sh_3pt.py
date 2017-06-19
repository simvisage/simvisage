#-------------------------------------------------------------------------
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
# Created on Feb 15, 2010 by: rch

from numpy import  fabs, where, copy, \
    argmax, unique, around
from traits.api import \
    Int, Float, \
    on_trait_change, Instance, \
    Array, Property, cached_property, \
    Event, implements, DelegatesTo
from traitsui.api import \
    View, Item, VGroup, \
    Group

from matresdev.db.exdb.ex_type import ExType
from matresdev.db.exdb.i_ex_type import IExType
from matresdev.db.matdb.trc.composite_cross_section import \
    CompositeCrossSection, plain_concrete
from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp


class ExpSH3PT(ExType):
    '''Experiment: Shear Test Three Point
    '''
#    label = Str('four point bending test')

    implements(IExType)

    file_ext = 'DAT'

    #--------------------------------------------------------------------
    # register a change of the traits with metadata 'input'
    #--------------------------------------------------------------------

    input_change = Event

    @on_trait_change('+input, ccs.input_change, +ironing_param')
    def _set_input_change(self):
        self.input_change = True

    #-------------------------------------------------------------------------
    # specify inputs:
    #-------------------------------------------------------------------------

    # effective length of the shear test specimen
    # (does not include the part at each side of the specimens that leaps over the support lines)
    #
    length = Float(0.12, unit='m', input=True, table_field=True,
                   auto_set=False, enter_set=True)
    length_left = Float(0.06, unit='m', input=True, table_field=True,
                        auto_set=False, enter_set=True)
    length_right = Float(0.06, unit='m', input=True, table_field=True,
                         auto_set=False, enter_set=True)
    width = Float(0.19, unit='m', input=True, table_field=True,
                  auto_set=False, enter_set=True)
    thickness = Float(0.03, unit='m', input=True, table_field=True,
                      auto_set=False, enter_set=True)

    # age of the concrete at the time of testing
    age = Int(242, unit='d', input=True, table_field=True,
              auto_set=False, enter_set=True)
    loading_rate = Float(1.0, unit='mm/min', input=True, table_field=True,
                         auto_set=False, enter_set=True)

    # additional own weight of load introduction
    weight_load_introduction = Float(0.0, unit='kN', input=True, table_field=True,
                                     auto_set=False, enter_set=True)

    #--------------------------------------------------------------------------
    # composite cross section
    #--------------------------------------------------------------------------

    ccs = Instance(CompositeCrossSection)

    def _ccs_default(self):
        '''default settings'
        '''
        # SFB 532 - demonstrator textil and concrete:
        fabric_layout_key = 'Q121/121-AAE-38'
        concrete_mixture_key = 'C3-HF2-155-5'
        orientation_fn_key = 'all0'
        n_layers = 1
        s_tex_z = 0.03 / (n_layers + 1)
        ccs = CompositeCrossSection(
            fabric_layup_list=[
                plain_concrete(s_tex_z * 0.5),
                FabricLayUp(
                    n_layers=n_layers,
                    orientation_fn_key=orientation_fn_key,
                    s_tex_z=s_tex_z,
                    fabric_layout_key=fabric_layout_key
                ),
                plain_concrete(s_tex_z * 0.5)
            ],
            concrete_mixture_key=concrete_mixture_key
        )
        return ccs

    #--------------------------------------------------------------------------
    # Get properties of the composite
    #--------------------------------------------------------------------------

    # E-modulus of the composite at the time of testing
    E_c = Property(
        Float, unit='MPa', depends_on='input_change', table_field=True)

    def _get_E_c(self):
        return self.ccs.get_E_c_time(self.age)

    # E-modulus of the composite after 28 days
    E_c28 = DelegatesTo('ccs', listenable=False)

    # reinforcement ration of the composite
    rho_c = DelegatesTo('ccs', listenable=False)

    #-------------------------------------------------------------------------
    # define processing
    #-------------------------------------------------------------------------

    # put this into the ironing procedure processor
    #
    jump_rtol = Float(0.9,
                      auto_set=False, enter_set=True,
                      ironing_param=True)

    def process_source_data(self):
        '''read in the measured data from file and assign
        attributes after array processing.
        If necessary modify the assigned data, i.e. change
        the sign or specify an offset for the specific test setup.
        '''
        print '*** process source data ***'

        super(ExpSH3PT, self).process_source_data()

        self._read_data_array()

        # set attributes:
        #
        self._set_array_attribs()

        # PEEKEL-measuring software:

        # convert units and change signs
        self.Kraft -= self.Kraft[0]
        self.Kraft *= -1

        # add weight of load introduction to force
        print 'add weight of steel traverse to force'
        self.Kraft += self.weight_load_introduction
        print 'force at initial state ', self.weight_load_introduction
        # @todo: interpolate an initial deformation based on the initial force and the initial stiffness
        #       measured in order to start the F-w-curve at the origin!

        # (reset displacement gauges by their initial values and change sign
        # in order to return a positive value for a displacement)

        # vertical displacements [mm]:
        self.Ver_re -= self.Ver_re[0]
        self.Ver_re *= -1
        
        self.Ver_li -= self.Ver_li[0]
        self.Ver_li *= -1        

    #-------------------------------------------------------------------------
    # plot templates
    #-------------------------------------------------------------------------

    plot_templates = {'force / gauge displacement': '_plot_force_displacement',
                      'force / deflection': '_plot_force_deflection',
                      'shearforce / deflection': '_plot_shearforce_deflection',
                      'force / deflection ascending': '_plot_force_deflection_asc',
                      'shearforce / deflection ascending':'_plot_shearforce_deflection_asc',
                      }

    default_plot_template = 'force / deflection'

    # get only the ascending branch of the response curve
    #
    max_force_idx = Property(Int)

    def _get_max_force_idx(self):
        '''get the index of the maximum force'''
        # NOTE: processed data returns positive values for force and
        # displacement
        return argmax(self.Kraft)

    def _plot_force_displacement(self, axes, offset_w=0., linestyle='-', linewidth=1., label=None):
        f = self.Kraft
        # add curves
        axes.plot(self.Ver_re, self.Kraft, color='b')
        axes.plot(self.Ver_li, self.Kraft, color='r')
        # add axes labels
        #
        xkey = 'displacement [mm]'
        ykey = 'force [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))

    def _plot_force_deflection(self, axes, offset_w=0., color='black', linestyle='-', linewidth=1., label=None):
        f = self.Kraft
        w = 0.5 * (self.Ver_re + self.Ver_li)

        # add curves
        #
        axes.plot(w, f, linewidth=linewidth, linestyle=linestyle, label=label, color=color)

        # add axes labels
        #
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))

    def _plot_shearforce_deflection(self, axes, offset_w=0., color='black', linestyle='-', linewidth=1., label=None):
        f = 0.5 * self.Kraft
        w = 0.5 * (self.Ver_re + self.Ver_li)

        # add curves
        #
        axes.plot(w, f, linewidth=linewidth, linestyle=linestyle, label=label, color=color)

        # add axes labels
        #
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))
        
    def _plot_force_deflection_asc(self, axes, offset_w=0., color='black', linestyle='-', linewidth=1., label=None):
        f = self.Kraft
        w = 0.5 * (self.Ver_re + self.Ver_li)

        # add curves
        #
        axes.plot(w[:self.max_force_idx + 1], f[:self.max_force_idx + 1], linewidth=linewidth, linestyle=linestyle, label=label, color=color)

        # add axes labels
        #
        xkey = 'deflection [mm]'
        ykey = 'force [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))

    def _plot_shearforce_deflection_asc(self, axes, offset_w=0., color='black', linestyle='-', linewidth=1., label=None):
        f = 0.5 * self.Kraft
        w = 0.5 * (self.Ver_re + self.Ver_li)

        # add curves
        #
        axes.plot(w[:self.max_force_idx + 1], f[:self.max_force_idx + 1], linewidth=linewidth, linestyle=linestyle, label=label, color=color)

        # add axes labels
        #
        xkey = 'deflection [mm]'
        ykey = 'shearforce [kN]'
        axes.set_xlabel('%s' % (xkey,))
        axes.set_ylabel('%s' % (ykey,))

    #-------------------------------------------------------------------------
    # view
    #-------------------------------------------------------------------------

    traits_view = View(VGroup(
        Group(
            Item('length', format_str="%.3f"),
            Item('length_left', format_str="%.3f"),
            Item('length_right', format_str="%.3f"),
            Item('width', format_str="%.3f"),
            Item('thickness', format_str="%.3f"),
            label='geometry'
        ),
        Group(
            Item('weight_load_introduction'),
            Item('loading_rate'),
            Item('age'),
            label='loading rate and age'
        ),
        Group(
            Item('E_c', show_label=True, style='readonly', format_str="%.0f"),
            Item('ccs@', show_label=False),
            label='composite cross section'
        )
    ),
        scrollable=True,
        resizable=True,
        height=0.8,
        width=0.6
    )

if __name__ == '__main__':

    from matresdev.db.exdb.ex_run_table import ExRunClassExt
    ex = ExRunClassExt(klass=ExpSH3PT)
    ex.configure_traits()
