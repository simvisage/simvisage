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
# Created on Jun 27, 2017 by: rch, li, jb


import os
from traits.api import \
    Int, Float, \
    on_trait_change, Instance, \
    Array, Property, cached_property, \
    Bool, Event, implements, \
    DelegatesTo, Callable
from traitsui.api \
    import View, Item, HSplit, Group, VSplit

from matresdev.db.exdb import ExRun
from matresdev.db.exdb.ex_run_table import ExRunClassExt
from matresdev.db.exdb.ex_type import ExType
from matresdev.db.exdb.i_ex_type import IExType
from matresdev.db.matdb.trc.composite_cross_section \
    import CompositeCrossSection, plain_concrete
from matresdev.db.matdb.trc.concrete_mixture \
    import ConcreteMixture
from matresdev.db.matdb.trc.fabric_layout \
    import FabricLayOut
from matresdev.db.matdb.trc.fabric_layup \
    import FabricLayUp
from matresdev.db.simdb.simdb import simdb
import numpy as np


class ExpSPO(ExType):

    '''Experiment: Single Sided Pullout-test with clamped specimen
    '''

    implements(IExType)

    file_ext = 'DAT'

    # --------------------------------------------------------------------
    # register a change of the traits with metadata 'input'
    # --------------------------------------------------------------------

    input_change = Event

    @on_trait_change('+input, ccs.input_change')
    def _set_input_change(self):
        print '*** raising input change in CTT'
        self.input_change = True

    avg_disp = Callable
    # -------------------------------------------------------------------------
    # specify inputs:
    # -------------------------------------------------------------------------

    width = Float(0.12, unit='m', input=True, table_field=True,
                  auto_set=False, enter_set=True)
    '''Width of the cross section.
    '''

    n_rovings = Int(1, unit='-', input=True, table_field=True,
                    auto_set=False, enter_set=True)
    '''Number of rovings within a cross section.
    '''

    A_roving = Float(0.000001, unit='mm2', input=True, table_field=True,
                     auto_set=False, enter_set=True)
    '''Cross sectional area of a roving.
    '''

    gauge_length = Float(0.25, unit='m', input=True, table_field=True,
                         auto_set=False, enter_set=True)
    '''Gauge length is just for information on the used test setup.
    For this test setup, only one single crack should form
    '''

    age = Int(28, unit='d', input=True, table_field=True,
              auto_set=False, enter_set=True)
    '''Age of the concrete at the time of testing.
    '''

    loading_rate = Float(1.0, unit='mm/min', input=True, table_field=True,
                         auto_set=False, enter_set=True)
    '''Applied loading rate.
    '''

    # --------------------------------------------------------------------------
    # composite cross section
    # --------------------------------------------------------------------------

    ccs = Instance(CompositeCrossSection)
    '''Link to the composite cross section object containing
    defining the layout and lay-up of the reinforcement fabrics.
    '''

    def _ccs_default(self):
        '''default settings correspond to
        setup '9u_MAG-07-03_PZ-0708-1'
        '''
        print 'ccs default used'
        fabric_layout_key = 'Q121/121-AAE-38'
        concrete_mixture_key = 'C3-HF2-155-5'
        orientation_fn_key = 'all0'
        n_layers = 1
        thickness = 0.021

        s_tex_z = thickness / (n_layers + 1)
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

  # --------------------------------------------------------------------------
    # Get properties of the composite
    # --------------------------------------------------------------------------

    E_c = Property(
        Float, unit='MPa', depends_on='input_change', table_field=True)
    '''E-modulus of the composite at the time of testing.
    '''

    @cached_property
    def _get_E_c(self):
        return self.ccs.get_E_c_time(self.age)

    E_m = Property(
        Float, unit='MPa', depends_on='input_change', table_field=True)
    '''E-modulus of the concrete matrix at the time of testing.
    '''

    @cached_property
    def _get_E_m(self):
        return self.ccs.get_E_m_time(self.age)

    A_c = Property(Float, unit='m^2', depends_on='input_change')
    '''Cross-sectional-area of the composite.
    '''

    @cached_property
    def _get_A_c(self):
        return self.width * self.ccs.thickness

    A_tex = Property(Float, unit='mm^2', depends_on='input_change')
    '''Total cross-sectional-area of the textile reinforcement.
    '''

    @cached_property
    def _get_A_tex(self):
        return self.ccs.a_tex * self.width

    E_c28 = DelegatesTo('ccs', listenable=False)
    '''E-modulus of the composite after 28 days.
    '''

    rho_c = DelegatesTo('ccs', listenable=False)
    '''Reinforcement ratio of the composite.
    '''

    # -------------------------------------------------------------------------
    # define processing
    # -------------------------------------------------------------------------

    def process_source_data(self):
        '''Extend the default data processing with
        the handling of gauges placed at the front and back side of
        the specimen. If necessary modify the assigned data, i. e. change
        the sign or specify an offset for the specific test setup.
        @todo: make the distinction in the subclasses.
        '''
        super(ExpSPO, self).process_source_data()

        if self.processing_done:
            return

        # NOTE: the small  test specimen (INSTRON) with width = 0.12 m have
        # only 2 or 3 displacement gauges
        #
        if hasattr(self, "W10_re") and hasattr(self, "W10_li") \
                and hasattr(self, "W10_vo"):
            self.W10_re -= self.W10_re[0]
            self.W10_re *= -1
            self.W10_li -= self.W10_li[0]
            self.W10_li *= -1
            self.W10_vo -= self.W10_vo[0]
            self.W10_vo *= -1

        if hasattr(self, "W10_re") and hasattr(self, "W10_li") \
                and not hasattr(self, "W10_vo"):
            self.W10_re -= self.W10_re[0]
            self.W10_re *= -1
            self.W10_li -= self.W10_li[0]
            self.W10_li *= -1

    '''Stress in the textile fabrics cross section.
    '''
    @cached_property
    def _get_sig_tex(self):
        # measured force:
        force = self.Kraft  # [kN]
        # cross sectional area of the reinforcement:
        A_tex = self.A_tex
        # calculated stresses:
        sig_tex = (force * 1000.) / A_tex  # [MPa]
        return sig_tex

    # -------------------------------------------------------------------------
    # plot templates
    # -------------------------------------------------------------------------

    plot_templates = {'force / machine displacement': '_plot_force_displacement_machine',
                      'force / gauge displacement': '_plot_force_displacement',
                      'per yarn force / gauge displacement': '_plot_yforce_displacement',
                      'textile stress / gauge displacement': '_plot_sigtex_displacement',
                      }

    default_plot_template = 'per yarn force / gauge displacement'

    def _plot_force_displacement_machine(self, axes, color='black',
                                         linewidth=1., linestyle='-',
                                         label=None):
        '''plot force-displacement diagram
        '''
        if hasattr(self, "Weg") and hasattr(self, "Kraft"):
            axes.plot(-self.Weg, self.Kraft, color=color,
                      linewidth=linewidth, linestyle=linestyle, label=label)
            axes.set_xlabel('displacement [mm]')
            axes.set_ylabel('force [kN]')

    def _plot_force_displacement(self, axes):
        '''plot force-displacement diagram
        '''
        if hasattr(self, "W10_re") and hasattr(self, "W10_li"):
            axes.plot(self.W10_re, self.Kraft)
            axes.plot(self.W10_li, self.Kraft)
            if hasattr(self, "W10_vo"):
                axes.plot(self.W10_vo, self.Kraft)
#            axes.set_xlabel('%s' % ('displacement [mm]',))
#            axes.set_ylabel('%s' % ('force [kN]',))

    def _plot_yforce_displacement(self, axes, color='black',
                                  linewidth=1., linestyle='-',
                                  label=None):
        '''plot per-yarn force-displacement diagram
        '''
        F_max_idx = np.argmax(self.Kraft)

        print 'PLOTTING'
        axes.plot(self.w, self.Kraft / self.n_rovings, color=color,
                  linewidth=linewidth, linestyle=linestyle, label=label)
        axes.set_xlabel('displacement [mm]')
        axes.set_ylabel('force [kN]')

    def _plot_sigtex_displacement(self, axes, color='black',
                                  linewidth=1., linestyle='-',
                                  label=None):
        pass

    # ---------------------------------
    # view
    # ---------------------------------

    traits_view = View(VSplit(
        HSplit(Group(
            Item('width', format_str="%.3f"),
            Item('gauge_length', format_str="%.3f"),
            springy=True,
            label='geometry',
            id='matresdev.db.exdb.ex_composite_tensile_test.geometry',
            dock='tab',
        ),
            Group(
            Item('loading_rate'),
            Item('age'),
            springy=True,
            label='loading rate and age',
            id='matresdev.db.exdb.ex_composite_tensile_test.loading',
            dock='tab',),
            id='matresdev.db.exdb.ex_composite_tensile_test.xxx',
            dock='tab',
        ),
        Group(
            Item('n_rovings'),
            Item('A_roving'),
            Item('ccs@', resizable=True, show_label=False),
            label='composite cross section'
        ),
        #                               label = 'input variables',
        #                               id = 'matresdev.db.exdb.ex_composite_tensile_test.vgroup.inputs',
        #                               dock = 'tab',
        #                               scrollable = True,
        #                               ),
        Group(
            Item('E_c', visible_when='derived_data_available',
                 style='readonly', show_label=True, format_str="%.0f"),
            #             Item('PYF_max', visible_when='derived_data_available',
            #                  style='readonly', emphasized=True, format_str="%.2f"),
            #             Item('sig_tex_max', visible_when='derived_data_available',
            # style='readonly', emphasized=True, format_str="%.2f"),
            label='output characteristics',
            id='matresdev.db.exdb.ex_composite_tensile_test.vgroup.outputs',
            dock='tab',
            scrollable=True,
        ),
        scrollable=True,
        id='matresdev.db.exdb.ex_composite_tensile_test.vgroup',
        dock='tab',
    ),
        id='matresdev.db.exdb.ex_composite_tensile_test',
        dock='tab',
        scrollable=True,
        resizable=True,
        height=0.8,
        width=0.5,
    )

ExpSPO.db = ExRunClassExt(klass=ExpSPO)

# if __name__ == '__main__':
# 
#     import pylab as p
#     fig = p.figure(facecolor='white', figsize=(12, 9))
# 
#     test_file_path = os.path.join(simdb.exdata_dir,
#                                   'double_pullout',
#                                   '2016-03-16_DPO-15mm-0-3300SBR_R4',
#                                   'raw_data',
#                                   'DPO-70cm-0-3300SBR-V2_R4.DAT'
#                                   )
#     exrun = ExRun(data_file=test_file_path)
# 
#     axes = p.subplot(111)
#     exrun.ex_type._plot_yforce_displacement(axes)
#     p.show()
