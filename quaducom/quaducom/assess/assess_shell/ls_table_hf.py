'''
Created on Jun 23, 2010

@author: alexander
'''

from traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, Trait

from traitsui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from mayavi import \
    mlab

# from etsproxy.mayavi.mlab import \
#    colorbar, show, points3d
#
# from etsproxy.mayavi.api import \
#    Engine


from traitsui.table_column import \
    ObjectColumn

from traitsui.menu import \
    OKButton, CancelButton

from traitsui.tabular_adapter \
    import TabularAdapter

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, sqrt, frompyfunc, \
                max as ndmax, ones

import numpy as np

from math import pi
from string import split
import os

from .lcc_reader import LCCReader, LCCReaderRFEM, LCCReaderInfoCAD

class LSArrayAdapter (TabularAdapter):

    columns = Property
    def _get_columns(self):
#        print 'GETTING COLUMNS', self.object.columns, self.object, self.object.__class__
        columns = self.object.columns
        return [ (name, idx) for idx, name in enumerate(columns) ]

    font = 'Courier 10'
    alignment = 'right'
    format = '%5.2f'  # '%g'
    even_bg_color = Color(0xE0E0FF)
    width = Float(80)

    # @todo: format columns using 'column_id'
#    adapter_column_map = Property(depends_on = 'adapters,columns')


class LS(HasTraits):
    '''Limit state class
    '''

    # backward link to the info shell to access the
    # input data when calculating
    # the limit-state-specific values
    #
    ls_table = WeakRef

    #-------------------------------
    # ls columns
    #-------------------------------
    # defined in the subclasses
    #
    ls_columns = List
    show_ls_columns = Bool(True)

    #-------------------------------
    # geo columns from info shell
    #-------------------------------

    geo_columns = List([ 'node_no' ])
    show_geo_columns = Bool(True)

    node_no = Property(Array)
    def _get_node_no(self):
        return self.ls_table.node_no

    #-------------------------------
    # state columns from info shell
    #-------------------------------

    state_columns = List(['Rx', 'Ry', 'Rz', 'Mx', 'My', 'Mz' ])

    show_state_columns = Bool(True)

    # forces at the supports [kN]
    #
    Rx = Property(Array)
    def _get_Rx(self):
        return self.ls_table.Rx

    Ry = Property(Array)
    def _get_Ry(self):
        return self.ls_table.Ry

    Rz = Property(Array)
    def _get_Rz(self):
        return self.ls_table.Rz

    #-------------------------------
    # ls table
    #-------------------------------

    # all columns associated with the limit state including the corresponding
    # stress resultants
    #
    columns = Property(List, depends_on='show_geo_columns, show_state_columns, show_ls_columns')
    @cached_property
    def _get_columns(self):
        columns = []

        if self.show_geo_columns:
            columns += self.geo_columns

        if self.show_state_columns:
            columns += self.state_columns

        if self.show_ls_columns:
            columns += self.ls_columns

        return columns

    # select column used for sorting the data in selected sorting order
    #
    sort_column = Enum(values='columns')
    def _sort_column_default(self):
        return self.columns[-1]

    sort_order = Enum('descending', 'ascending', 'unsorted')

    #-------------------------------------------------------
    # get the maximum value of the selected variable
    # 'max_in_column' of the current sheet (only one sheet)
    #-------------------------------------------------------

    # get the maximum value of the chosen column
    #
    max_in_column = Enum(values='columns')
    def _max_in_column_default(self):
        return self.columns[-1]

    max_value = Property(depends_on='max_in_column')
    def _get_max_value(self):
        col = getattr(self, self.max_in_column)[:, 0]
        return max(col)

    #-------------------------------------------------------
    # get the maximum value and the corresponding case of
    # the selected variable 'max_in_column' in all (!) sheets
    #-------------------------------------------------------

    max_value_all = Property(depends_on='max_in_column')
    def _get_max_value_all(self):
        return self.ls_table.max_value_and_case[ self.max_in_column ]['max_value']

    max_case = Property(depends_on='max_in_column')
    def _get_max_case(self):
        return self.ls_table.max_value_and_case[ self.max_in_column ]['max_case']

    #-------------------------------------------------------
    # get ls_table for View
    #-------------------------------------------------------

    # stack columns together for table used by TabularEditor
    #
    ls_array = Property(Array, depends_on='sort_column, sort_order, \
                                              show_geo_columns, \
                                              show_state_columns, \
                                              show_ls_columns')

    @cached_property
    def _get_ls_array(self):

        arr_list = [ getattr(self, col) for col in self.columns ]

        # get the array currently selected by the sort_column enumeration
        #
        sort_arr = getattr(self, self.sort_column)[:, 0]
        sort_idx = argsort(sort_arr)
        ls_array = hstack(arr_list)

        if self.sort_order == 'descending':
            return ls_array[ sort_idx[::-1] ]
        if self.sort_order == 'ascending':
            return ls_array[ sort_idx ]
        if self.sort_order == 'unsorted':
            return ls_array

    #---------------------------------
    # plot outputs in mlab-window
    #---------------------------------
    warp_factor = Float(1000., input=True)

    plot_column = Enum(values='columns')
    plot = Button
    def _plot_fired(self):

        plot_col = getattr(self, self.plot_column)[:, 0]

        mlab.figure(figure="SFB532Demo",
                     bgcolor=(1.0, 1.0, 1.0),
                     fgcolor=(0.0, 0.0, 0.0))

        gd = self.ls_table.geo_data
        sd = self.ls_table.state_data

        if self.plot_column == 'n_tex':
            plot_col = where(plot_col < 0, 0, plot_col)

        r = self.ls_table.reader
        r.plot_col(mlab, plot_col, gd, state_data=sd, warp_factor=self.warp_factor)

        mlab.scalarbar(title=self.plot_column, orientation='vertical')
        mlab.show

    # name of the trait that is used to assess the evaluated design
    #
    assess_name = Property(Str)
    def _get_assess_name(self):
        return self.ls_table.assess_name

    #-------------------------------
    # ls group
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed
    # does not work in connection with the LSArrayAdapter
    ls_group = VGroup(
                        HGroup(# Item( 'assess_name' ),
                                Item('max_in_column'),
                                Item('max_value', style='readonly', format_str='%6.2f'),
                              ),
                        HGroup(Item('sort_column'),
                                Item('sort_order'),
                                Item('show_geo_columns', label='show geo'),
                                Item('show_state_columns', label='show state'),
                                Item('show_ls_columns', label='show ls'),
                                Item('plot_column'),
                                Item('plot'),
                                Item('warp_factor')
                              ),
                     )


class SLSRxyz(LS):
    '''Serviceability limit state
    '''

    # ------------------------------------------------------------
    # SLS: outputs
    # ------------------------------------------------------------

    ls_columns = List([])

    ls_values = Property(depends_on='+input')
    @cached_property
    def _get_ls_values(self):
        '''get the outputs for SLS
        '''
        return {}

    assess_name = ''

    #-------------------------------
    # ls view
    #-------------------------------
    # @todo: the dynamic selection of the columns to be displayed
    # does not work in connection with the LSArrayAdapter
    traits_view = View(VGroup(
                            HGroup(Item(name='', label=''),
                                   Item(name='', label='')
                                   ),
                            VGroup(
                                Include('ls_group'),
                                Item('ls_array', show_label=False,
                                      editor=TabularEditor(adapter=LSArrayAdapter()))
                                  ),
                              ),
                      resizable=True,
                      scrollable=True,
                      height=1000,
                      width=1100
                      )

class ULSRxyz(LS):
    '''Ultimate limit state
    '''

    #--------------------------------------------------------
    # ULS: material parameters (Inputs)
    #--------------------------------------------------------

    # shear Resistance
    #
    Rx_Rd = Float(4.8, input=True)

    # pull-out Resistance
    #
    Rz_Rd = Float(4.7, input=True)

    # (unused as Ry = 0. for all cases)
    #
    Ry_Rd = Float(1., input=True)

    Mx_Rd = Float(1., input=True)
    My_Rd = Float(1., input=True)
    Mz_Rd = Float(1., input=True)

    # ------------------------------------------------------------
    # ULS: outputs
    # ------------------------------------------------------------

    ls_values = Property(depends_on='+input')
    @cached_property
    def _get_ls_values(self):
        '''get the outputs for ULS
        '''
        #---------------------------------------------------------
        # conditions for case distinction
        # (-- tension / compression reactions --)
        #---------------------------------------------------------
        # @todo: use this in order to sort out the cases with compression
        # up to now for cases with compression only the positive part is taken into account
        # leading to eta=0 for this direction and a possible non-zero value in the other direction

#        # reaction force tangential to the shell
#        #
#        cond_Rx_ge_0 = self.Rx >= 0. # positive value corresponds to shear in the screw
#        cond_Rx_le_0 = self.Rx <= 0. # negative value corresponds to compression
#
#        # reaction force radial to the shell
#        #
#        cond_Rz_ge_0 = self.Rz >= 0. # positive value corresponds to compression
#        cond_Rz_le_0 = self.Rz <= 0. # negative value corresponds to pull-out force in the screw

        #---------------------------------------------------------
        # resulting reaction forces and 'eta'
        #---------------------------------------------------------

        # evaluate resulting forces and moments
        #
        Rres = sqrt(self.Rx * self.Rx + self.Ry * self.Ry + self.Rz * self.Rz)
        Mres = sqrt(self.Mx * self.Mx + self.My * self.My + self.Mz * self.Mz)
        # note: positive values of 'Rx' correspond to shear forces for the support screw
        #       negative values are taken by the compression cushion at the support directly
        #       Therefore take only the positive part of support force 'Rx' into account
        #       for the evaluation of 'eta_Rx'
        Rx_pos = (abs(self.Rx) + self.Rx) / 2.

        # eta shear forces
        #
        eta_Rx = Rx_pos / self. Rx_Rd


        # note: negative values of 'Rz' correspond to pull-out forces for the support screw
        #       positive values are taken by the compression cushion at the support directly
        #       Therefore take only the negative values of the support force 'Rz' into account
        #       for the evaluation of 'eta_Rz'
        Rz_neg = (abs(self.Rz) - self.Rz) / 2.

        # eta pull-out
        #
        eta_Rz = Rz_neg / self. Rz_Rd

        # eta shear forces (unused as shear force in y-direction is always 0.)
        #
        eta_Ry = abs(self.Ry) / self. Ry_Rd

        # total eta for linear interaction:
        #
        eta_R_tot = eta_Rx + eta_Rz

        eta_Mx = self.Mx / self. Mx_Rd
        eta_My = self.My / self. My_Rd
        eta_Mz = self.Mz / self. Mz_Rd

        #------------------------------------------------------------
        # construct a dictionary containing the return values
        #------------------------------------------------------------

        return {
                 'Rres' : Rres,
                 'Mres' : Mres,
                 'eta_Rx' : eta_Rx,
                 'eta_Ry' : eta_Ry,
                 'eta_Rz' : eta_Rz,
                 'eta_R_tot' : eta_R_tot,
                 'eta_Mx' : eta_Mx,
                 'eta_My' : eta_My,
                 'eta_Mz' : eta_Mz,
               }

    #-----------------------------------------------
    # LS_COLUMNS: specify the properties that are displayed in the view
    #-----------------------------------------------

    # NOTE: the definition of ls_table.assess_name is given in constructor of 'LCCTable'
    #
#    assess_name = 'max_Rx' # @todo: compare with shear resistance of the screw
#    assess_name = 'min_Rx'
#    assess_name = 'max_Ry'
#    assess_name = 'min_Ry'
#    assess_name = 'max_Rz'
#    assess_name = 'min_Rz' # @todo: compare with pull-out resistance of the screw
#    assess_name = 'max_Rres'
    assess_name = 'max_eta_R_tot'

    ls_columns = List(['Rx', 'Ry', 'Rz', 'Rres',
                       'Mx', 'My', 'Mz', 'Mres',
                       'eta_Rx', 'eta_Ry', 'eta_Rz', 'eta_R_tot',
                       'eta_Mx', 'eta_My', 'eta_Mz'])

    Rres = Property(Array)
    def _get_Rres(self):
        return self.ls_values['Rres']

    Mres = Property(Array)
    def _get_Mres(self):
        return self.ls_values['Mres']

    eta_Rx = Property(Array)
    def _get_eta_Rx(self):
        return self.ls_values['eta_Rx']

    eta_Ry = Property(Array)
    def _get_eta_Ry(self):
        return self.ls_values['eta_Ry']

    eta_Rz = Property(Array)
    def _get_eta_Rz(self):
        return self.ls_values['eta_Rz']

    eta_R_tot = Property(Array)
    def _get_eta_R_tot(self):
        return self.ls_values['eta_R_tot']

    eta_Mx = Property(Array)
    def _get_eta_Mx(self):
        return self.ls_values['eta_Mx']

    eta_My = Property(Array)
    def _get_eta_My(self):
        return self.ls_values['eta_My']

    eta_Mz = Property(Array)
    def _get_eta_Mz(self):
        return self.ls_values['eta_Mz']

    #-------------------------------------------------
    # choose the assess parameter used for sorting
    # defined by the property name 'assess_name'
    #-------------------------------------------------

    max_Rx = Property(depends_on='+input')
    @cached_property
    def _get_max_Rx(self):
        return np.max(self.Rx)

    min_Rx = Property(depends_on='+input')
    @cached_property
    def _get_min_Rx(self):
        return np.min(self.Rx)

    max_Ry = Property(depends_on='+input')
    @cached_property
    def _get_max_Ry(self):
        return np.max(self.Ry)

    min_Ry = Property(depends_on='+input')
    @cached_property
    def _get_min_Ry(self):
        return np.min(self.Ry)

    max_Rz = Property(depends_on='+input')
    @cached_property
    def _get_max_Rz(self):
        return np.max(self.Rz)

    min_Rz = Property(depends_on='+input')
    @cached_property
    def _get_min_Rz(self):
        return np.min(self.Rz)

    max_Rres = Property(depends_on='+input')
    @cached_property
    def _get_max_Rres(self):
        return ndmax(self.Rres)

    max_eta_R_tot = Property(depends_on='+input')
    @cached_property
    def _get_max_eta_R_tot(self):
        return ndmax(self.eta_R_tot)

    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed
    # does not work in connection with the LSArrayAdapter
    traits_view = View(
                       VGroup(
                        HGroup(
                            VGroup(
                                Item(name='Rx_Rd', label='resistance R_xd [kN]', style='readonly', format_str="%.1f"),
                                Item(name='Ry_Rd', label='resistance R_yd [kN]', style='readonly', format_str="%.1f"),
                                Item(name='Rz_Rd', label='resistance R_zd [kN]', style='readonly', format_str="%.1f"),
                                label='material properties (longitudinal)'
                                  ),
                            VGroup(
                                Item(name='assess_name', label='assess_name', style='readonly', format_str="%s"),
                                label='sort rows according to'
                                  )
                             ),

                        VGroup(
                            Include('ls_group'),
                            Item('ls_array', show_label=False,
                                  editor=TabularEditor(adapter=LSArrayAdapter()))
                              ),
                            ),
                      resizable=True,
                      scrollable=True,
                      height=1000,
                      width=1100
                      )

LSLIST = [ SLSRxyz, ULSRxyz ]

class LSTableRxyz(HasTraits):
    '''Assessment tool
    '''

    is_id = Int(0)
    # geo data: coordinates and element thickness
    #
    geo_data = Dict

    node_no = Property(Array)
    def _get_node_no(self):
#        return self.state_data['node_no']
        return self.geo_data['node_no']

    # state data: stress resultants
    #
    state_data = Dict

    Rx = Property(Array)
    def _get_Rx(self):
        return self.state_data['Rx']

    Ry = Property(Array)
    def _get_Ry(self):
        return self.state_data['Ry']

    Rz = Property(Array)
    def _get_Rz(self):
        return self.state_data['Rz']

    Mx = Property(Array)
    def _get_Mx(self):
        return self.state_data['Mx']

    My = Property(Array)
    def _get_My(self):
        return self.state_data['My']

    Mz = Property(Array)
    def _get_Mz(self):
        return self.state_data['Mz']

    #------------------------------------------
    # combinations of limit states, stress resultants and directions
    #------------------------------------------

    ls = Trait('ULS',
                {'ULS' : ULSRxyz,
                 'SLS' : SLSRxyz })

    ls_class = Instance(LS)
    def _ls_class_default(self):
        '''ls instances, e.g. ULS()
        '''
        ls_class = self.ls_
        return ls_class(ls_table=self)

    assess_name = Str

    assess_value = Property
    def _get_assess_value(self):
        ls = self.ls_class
        return getattr(ls, self.assess_name)

    traits_view = View(Tabbed(
                            Item('ls_class@' , label="ls", show_label=False),
                            scrollable=False,
                         ),
                      resizable=True,
                      scrollable=True,
                      height=1000,
                      width=1100
                      )


