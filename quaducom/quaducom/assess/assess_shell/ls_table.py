'''
Created on Jun 23, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, Trait

from etsproxy.util.home_directory import \
    get_home_directory

from etsproxy.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from etsproxy.mayavi import \
    mlab

import numpy as np

from etsproxy.traits.ui.table_column import \
    ObjectColumn

from etsproxy.traits.ui.menu import \
    OKButton, CancelButton

from etsproxy.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, sqrt, frompyfunc, \
                max as ndmax, ones

from math import pi
from string import split
import os

from lcc_reader import LCCReader, LCCReaderRFEM, LCCReaderInfoCAD, LCCReaderInfoCADRxyz

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
    # geo data shown in 'geo_columns'
    #-------------------------------

    geo_columns = List([ 'elem_no', 'X', 'Y', 'Z', 'thickness' ])
    show_geo_columns = Bool(True)

    elem_no = Property(Array)
    def _get_elem_no(self):
        return self.ls_table.elem_no

    X = Property(Array)
    def _get_X(self):
        return self.ls_table.X

    Y = Property(Array)
    def _get_Y(self):
        return self.ls_table.Y

    Z = Property(Array)
    def _get_Z(self):
        return self.ls_table.Z

    thickness = Property(Array)
    def _get_thickness(self):
        return self.ls_table.thickness

    #-------------------------------
    # state data shown in 'state_columns'
    #-------------------------------

    state_columns = List([
                           'mx', 'my', 'mxy', 'nx', 'ny', 'nxy',
                           'sigx_lo', 'sigy_lo', 'sigxy_lo',
#                           'sig1_lo', 'sig1_up_sig_lo',
                           'm_sig1_lo', 'n_sig1_lo',
                           'm_sig2_lo', 'n_sig2_lo',
                           'sigx_up', 'sigy_up', 'sigxy_up',
#                           'sig1_up', 'sig1_lo_sig_up', 'alpha_sig1_up',
                           'm_sig1_up', 'n_sig1_up',
                           'm_sig2_up', 'n_sig2_up',
                            ])

    show_state_columns = Bool(True)

    mx = Property(Array)
    def _get_mx(self):
        return self.ls_table.mx

    my = Property(Array)
    def _get_my(self):
        return self.ls_table.my

    mxy = Property(Array)
    def _get_mxy(self):
        return self.ls_table.mxy

    nx = Property(Array)
    def _get_nx(self):
        return self.ls_table.nx

    ny = Property(Array)
    def _get_ny(self):
        return self.ls_table.ny

    nxy = Property(Array)
    def _get_nxy(self):
        return self.ls_table.nxy

    #-------------------------------
    # derived state date
    # (upper face):
    #-------------------------------

    # calculate global stresses from state data (upper face)
    #
    sigx_up = Property(Array)
    def _get_sigx_up(self):
        return self.ls_table.sigx_up

    sigy_up = Property(Array)
    def _get_sigy_up(self):
        return self.ls_table.sigy_up

    sigxy_up = Property(Array)
    def _get_sigxy_up(self):
        return self.ls_table.sigxy_up

    # evaluate principal stress-direction (upper face)
    #
    alpha_sig1_up = Property(Array)
    def _get_alpha_sig1_up(self):
        return self.ls_table.alpha_sig1_up

    alpha_sig2_up = Property(Array)
    def _get_alpha_sig2_up(self):
        return self.ls_table.alpha_sig2_up

    alpha_sig1_up_deg = Property(Array)
    def _get_alpha_sig1_up_deg(self):
        return self.ls_table.alpha_sig1_up_deg

    alpha_sig2_up_deg = Property(Array)
    def _get_alpha_sig2_up_deg(self):
        return self.ls_table.alpha_sig2_up_deg

    # evaluate principal stresses (upper face)
    #
    sig1_up = Property(Array)
    def _get_sig1_up(self):
        return self.ls_table.sig1_up

    sig2_up = Property(Array)
    def _get_sig2_up(self):
        return self.ls_table.sig2_up

    # evaluate n(alpha),m(alpha) in principal stress direction (upper face)
    #
    n_sig1_up = Property(Array)
    def _get_n_sig1_up(self):
        return self.ls_table.n_sig1_up

    m_sig1_up = Property(Array)
    def _get_m_sig1_up(self):
        return self.ls_table.m_sig1_up

    n_sig2_up = Property(Array)
    def _get_n_sig2_up(self):
        return self.ls_table.n_sig2_up

    m_sig2_up = Property(Array)
    def _get_m_sig2_up(self):
        return self.ls_table.m_sig2_up

    #-------------------------------
    # derived state date
    # (lower face):
    #-------------------------------

    # calculate global stresses from state data (lower face)
    #
    sigx_lo = Property(Float)
    def _get_sigx_lo(self):
        return self.ls_table.sigx_lo

    sigy_lo = Property(Float)
    def _get_sigy_lo(self):
        return self.ls_table.sigy_lo

    sigxy_lo = Property(Float)
    def _get_sigxy_lo(self):
        return self.ls_table.sigxy_lo

    # evaluate principal stress-direction (lower face)
    #
    alpha_sig1_lo = Property(Float)
    def _get_alpha_sig1_lo(self):
        return self.ls_table.alpha_sig1_lo

    alpha_sig2_lo = Property(Array)
    def _get_alpha_sig2_lo(self):
        return self.ls_table.alpha_sig2_lo

    alpha_sig1_lo_deg = Property(Float)
    def _get_alpha_sig1_lo_deg(self):
        return self.ls_table.alpha_sig1_lo_deg

    alpha_sig2_lo_deg = Property(Float)
    def _get_alpha_sig2_lo_deg(self):
        return self.ls_table.alpha_sig2_lo_deg

    # evaluate principal stresses (lower face)
    #
    sig1_lo = Property(Float)
    def _get_sig1_lo(self):
        return self.ls_table.sig1_lo

    sig2_lo = Property(Float)
    def _get_sig2_lo(self):
        return self.ls_table.sig2_lo

    # evaluate n(alpha),m(alpha) in principal stress direction (lower face)
    #
    n_sig1_lo = Property(Array)
    def _get_n_sig1_lo(self):
        return self.ls_table.n_sig1_lo

    m_sig1_lo = Property(Array)
    def _get_m_sig1_lo(self):
        return self.ls_table.m_sig1_lo

    n_sig2_lo = Property(Array)
    def _get_n_sig2_lo(self):
        return self.ls_table.n_sig2_lo

    m_sig2_lo = Property(Array)
    def _get_m_sig2_lo(self):
        return self.ls_table.m_sig2_lo

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
    warp_factor = Float(100., input=True)

    plot_column = Enum(values='columns')
    plot = Button
    def _plot_fired(self):

        plot_col = getattr(self, self.plot_column).flatten()
        if self.plot_column == 'n_tex':
            plot_col = where(plot_col < 0, 0, plot_col)

        mlab.figure(figure="SFB532Demo",
                     bgcolor=(1.0, 1.0, 1.0),
                     fgcolor=(0.0, 0.0, 0.0))

        gd = self.ls_table.geo_data
        sd = self.ls_table.state_data

        r = self.ls_table.reader
        # use plotting function defined by the specific LCCTableReader
        # extract global coordinates ('X','Y','Z') from 'geo_data' and
        # global displacements ('ux_elem','uy_elem','uz_elem') from 'state_data'
        # if this information is available (distinguished by the specific Reader)
        r.plot_col(mlab, plot_col, gd, state_data=sd, warp_factor=self.warp_factor)

        mlab.scalarbar(title=self.plot_column, orientation='vertical')
        mlab.show

    # name of the trait that is used to assess the evaluated design
    #
    assess_name = Str('')

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


class SLS(LS):
    '''Serviceability limit state
    '''
    # ------------------------------------------------------------
    # SLS: material parameters (Inputs)
    # ------------------------------------------------------------

    # tensile strength [MPa]
    f_ctk = Float(5.0, input=True)

    # flexural tensile strength [MPa]
    f_m = Float(8.0, input=True)

    # ------------------------------------------------------------
    # SLS: outputs
    # ------------------------------------------------------------

    ls_columns = List(['m_sig1_lo', 'n_sig1_lo',
                       'm_sig1_up', 'n_sig1_up',
                       'eta_n_sig1_lo', 'eta_m_sig1_lo', 'eta_tot_sig_lo',
                       'eta_n_sig1_up', 'eta_m_sig1_up', 'eta_tot_sig_up',
                       'eta_tot'])

    ls_values = Property(depends_on='+input')
    @cached_property
    def _get_ls_values(self):
        '''get the outputs for SLS
        '''
        f_ctk = self.f_ctk
        f_m = self.f_m
        A = self.ls_table.A  # [m**2/m]
        W = self.ls_table.W  # [m**3/m]

        n_sig1_lo = self.ls_table.n_sig1_lo  # [kN/m]
        m_sig1_lo = self.ls_table.m_sig1_lo  # [kNm/m]
        sig_n_sig1_lo = n_sig1_lo / A / 1000.  # [MPa]
        sig_m_sig1_lo = m_sig1_lo / W / 1000.  # [MPa]

        n_sig1_up = self.ls_table.n_sig1_up
        m_sig1_up = self.ls_table.m_sig1_up
        sig_n_sig1_up = n_sig1_up / A / 1000.
        sig_m_sig1_up = m_sig1_up / W / 1000.

        eta_n_sig1_lo = sig_n_sig1_lo / f_ctk
        eta_m_sig1_lo = sig_m_sig1_lo / f_m
        eta_tot_sig_lo = eta_n_sig1_lo + eta_m_sig1_lo

        eta_n_sig1_up = sig_n_sig1_up / f_ctk
        eta_m_sig1_up = sig_m_sig1_up / f_m
        eta_tot_sig_up = eta_n_sig1_up + eta_m_sig1_up

        eta_tot = ndmax(hstack([ eta_tot_sig_up, eta_tot_sig_lo]), axis=1)[:, None]

        return { 'sig_n_sig1_lo':sig_n_sig1_lo, 'sig_m_sig1_lo':sig_m_sig1_lo,
                 'sig_n_sig1_up':sig_n_sig1_up, 'sig_m_sig1_up':sig_m_sig1_up,
                 'eta_n_sig1_lo':eta_n_sig1_lo, 'eta_m_sig1_lo':eta_m_sig1_lo, 'eta_tot_sig_lo':eta_tot_sig_lo,
                 'eta_n_sig1_up':eta_n_sig1_up, 'eta_m_sig1_up':eta_m_sig1_up, 'eta_tot_sig_up':eta_tot_sig_up,
                 'eta_tot':eta_tot, }

    # evaluate stresses at lower side:
    #
    sig_n_sig1_lo = Property
    def _get_sig_n_sig1_lo(self):
        return self.ls_values['sig_n_sig1_lo']

    sig_m_sig1_lo = Property
    def _get_sig_m_sig1_lo(self):
        return self.ls_values['sig_m_sig1_lo']

    eta_n_sig1_lo = Property
    def _get_eta_n_sig1_lo(self):
        return self.ls_values['eta_n_sig1_lo']

    eta_m_sig1_lo = Property
    def _get_eta_m_sig1_lo(self):
        return self.ls_values['eta_m_sig1_lo']

    eta_tot_sig_lo = Property
    def _get_eta_tot_sig_lo(self):
        return self.ls_values['eta_tot_sig_lo']

    # evaluate stresses at upper side:
    #
    sig_n_sig1_up = Property
    def _get_sig_n_sig1_up(self):
        return self.ls_values['sig_n_sig1_up']

    sig_m_sig1_up = Property
    def _get_sig_m_sig1_up(self):
        return self.ls_values['sig_m_sig1_up']

    eta_n_sig1_up = Property
    def _get_eta_n_sig1_up(self):
        return self.ls_values['eta_n_sig1_up']

    eta_m_sig1_up = Property
    def _get_eta_m_sig1_up(self):
        return self.ls_values['eta_m_sig1_up']

    eta_tot_sig_up = Property
    def _get_eta_tot_sig_up(self):
        return self.ls_values['eta_tot_sig_up']

    # total eta (upper AND lower side)
    #
    eta_tot = Property
    def _get_eta_tot(self):
        return self.ls_values['eta_tot']

    assess_name = 'max_eta_tot'

    max_eta_tot = Property(depends_on='+input')
    @cached_property
    def _get_max_eta_tot(self):
        return ndmax(self.eta_tot)

    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed
    # does not work in connection with the LSArrayAdapter
    traits_view = View(VGroup(
                            HGroup(Item(name='f_ctk', label='Tensile strength concrete [MPa]: f_ctk '),
                                    Item(name='f_m', label='Flexural tensile trength concrete [MPa]: f_m ')
                                   ),
                            VGroup(
                                Include('ls_group'),

                                # @todo: currently LSArrayAdapter must be called both
                                #        in SLS and ULS separately to configure columns
                                #        arrangement individually
                                #
                                Item('ls_array', show_label=False,
                                      editor=TabularEditor(adapter=LSArrayAdapter()))
                                  ),
                              ),
                      resizable=True,
                      scrollable=True,
                      height=1000,
                      width=1100
                      )

class ULS(LS):
    '''Ultimate limit state
    '''
    #-----------------------------------
    # strength characteristics (ULS)
    #-----------------------------------
    ### 'eval_mode = eta_nm' ###

    # tensile strength (0-direction) [kN/m]
    #
    n_0_Rdt = Property(Float)
    def _get_n_0_Rdt(self):
        return self.ls_table.strength_characteristics['n_0_Rdt']

    # bending (0-direction) strength [kNm/m]
    #
    m_0_Rd = Property(Float)
    def _get_m_0_Rd(self):
        return self.ls_table.strength_characteristics['m_0_Rd']

    # tensile strength (90-direction) [kN/m]
    #
    n_90_Rdt = Property(Float)
    def _get_n_90_Rdt(self):
        return self.ls_table.strength_characteristics['n_90_Rdt']

    # bending (90-direction) strength [kNm/m]
    #
    m_90_Rd = Property(Float)
    def _get_m_90_Rd(self):
        return self.ls_table.strength_characteristics['m_90_Rd']

    # compressive strength [kN/m]
    #
    n_Rdc = Property(Float)
    def _get_n_Rdc(self):
        return self.ls_table.strength_characteristics['n_Rdc']

    # ------------------------------------------------------------
    # ULS: outputs
    # ------------------------------------------------------------

    ls_values = Property(depends_on='+input')
    @cached_property
    def _get_ls_values(self):
        '''get the outputs for ULS
        '''
        #------------------------------------------------------------
        # get angle of deflection of the textile reinforcement
        #------------------------------------------------------------
        # angle 'beta' denotes the absolute deflection between the textile reinforcement and the
        # 1st or 2nd principle stress direction; no distiction of the rotation direction is necessary
        # for orthogonal textiles (standard case);
        # By default the 0-direction coincides with the x-axis

        # @todo: generalization of orientation of the 0-direction != x-axis
        # introduce the angle 'alpha_tex' to define the angle between the global x-axis and the
        # the 0-direction; then beta = abs(alpha_sig1 - alpha_tex); make sure that the deflection
        # angle still lies within -/+ pi/2 otherwise rotate by substracting n*pi with n = int(beta/pi)

        #-------------------------------------------
        # formulas are valid for the general case R_0 != R_90
        #-------------------------------------------
        #----------------------------------------------------------
        # upper side
        #----------------------------------------------------------
        beta_up_deg = abs(self.alpha_sig1_up_deg)  # [degree]
        beta_up = abs(self.alpha_sig1_up)  # [rad]
        n_Rdt_up_1 = self.n_0_Rdt * cos(beta_up) * (1 - beta_up_deg / 90.) + \
                   self.n_90_Rdt * sin(beta_up) * (beta_up_deg / 90.)
        m_Rd_up_1 = self.m_0_Rd * cos(beta_up) * (1 - beta_up_deg / 90.) + \
                  self.m_90_Rd * sin(beta_up) * (beta_up_deg / 90.)
        k_alpha_up = cos(beta_up) * (1 - beta_up_deg / 90.) + \
                     sin(beta_up) * (beta_up_deg / 90.)
        n_Rdt_up_2 = self.n_90_Rdt * cos(beta_up) * (1 - beta_up_deg / 90.) + \
                   self.n_0_Rdt * sin(beta_up) * (beta_up_deg / 90.)
        m_Rd_up_2 = self.m_90_Rd * cos(beta_up) * (1 - beta_up_deg / 90.) + \
                  self.m_0_Rd * sin(beta_up) * (beta_up_deg / 90.)

        #----------------------------------------------------------
        # lower side
        #----------------------------------------------------------
        beta_lo_deg = abs(self.alpha_sig1_lo_deg)  # [degree]
        beta_lo = abs(self.alpha_sig1_lo)  # [rad]
        n_Rdt_lo_1 = self.n_0_Rdt * cos(beta_lo) * (1 - beta_lo_deg / 90.) + \
                   self.n_90_Rdt * sin(beta_lo) * (beta_lo_deg / 90.)
        m_Rd_lo_1 = self.m_0_Rd * cos(beta_lo) * (1 - beta_lo_deg / 90.) + \
                  self.m_90_Rd * sin(beta_lo) * (beta_lo_deg / 90.)
        k_alpha_lo = cos(beta_lo) * (1 - beta_lo_deg / 90.) + \
                     sin(beta_lo) * (beta_lo_deg / 90.)
        n_Rdt_lo_2 = self.n_90_Rdt * cos(beta_lo) * (1 - beta_lo_deg / 90.) + \
                   self.n_0_Rdt * sin(beta_lo) * (beta_lo_deg / 90.)
        m_Rd_lo_2 = self.m_90_Rd * cos(beta_lo) * (1 - beta_lo_deg / 90.) + \
                  self.m_0_Rd * sin(beta_lo) * (beta_lo_deg / 90.)
        #-------------------------------------------------
        # NOTE: the principle tensile stresses are used to evaluate the principle direction\
        # 'eta_nm_tot' is evaluated based on linear nm-interaction (derived from test results)
        # 1st-principle stress direction = maximum (tensile) stresses
        # 2nd-principle stress direction = minimum (compressive) stresses
        #-------------------------------------------------
        #
        print "NOTE: the principle tensile stresses are used to evaluate the deflection angle 'beta', i.e. k_alpha(beta)"
        print "      'eta_nm_tot' is evaluated based on linear nm-interaction (derived from test results)"
        print "      'eta_nm_tot' is evaluated both for 1st and 2nd principle direction"

        if self.ls_table.k_alpha_min == True:
            print "minimum value 'k_alpha_min'=0.707 has been used to evaluate resistance values"
            # NOTE: conservative simplification: k_alpha_min = 0.707 used
            #
            n_Rdt_lo_1 = n_Rdt_up_1 = n_Rdt_lo_2 = n_Rdt_up_2 = min(self.n_0_Rdt, self.n_90_Rdt) * 0.707 * np.ones_like(self.elem_no)
            m_Rd_lo_1 = m_Rd_up_1 = m_Rd_lo_2 = m_Rd_up_2 = min(self.m_0_Rd, self.m_90_Rd) * 0.707 * np.ones_like(self.elem_no)

        #-------------------------------------------------
        # compressive strength (independent from principle and evaluation direction)
        #-------------------------------------------------
        n_Rdc = self.n_Rdc * np.ones_like(self.elem_no)

        #----------------------------------------
        # caluclate eta_nm
        #----------------------------------------
        # NOTE: destinction of the sign of the normal force necessary

        #---------------
        # 1-direction:
        #---------------

        # initialize arrays to be filled based on case distinction
        #
        eta_n_up = np.zeros_like(self.n_sig1_up)
        eta_n_lo = np.zeros_like(self.n_sig1_lo)

        # cases with a tensile normal force
        #
        cond_nsu_ge_0 = self.n_sig1_up >= 0.  # tensile force in direction of principle stress at upper side
        cond_nsl_ge_0 = self.n_sig1_lo >= 0.  # tensile force in direction of principle stress at lower side

        # compare imposed tensile normal force with 'n_Rd,t' as obtained from tensile test
        #
        bool_arr = cond_nsu_ge_0
        eta_n_up[bool_arr] = self.n_sig1_up[bool_arr] / n_Rdt_up_1[bool_arr]

        bool_arr = cond_nsl_ge_0
        eta_n_lo[bool_arr] = self.n_sig1_lo[bool_arr] / n_Rdt_lo_1[bool_arr]

        # cases with a compressive normal force
        #
        cond_nsu_lt_0 = self.n_sig1_up < 0.  # compressive force in direction of principle stress at upper side
        cond_nsl_lt_0 = self.n_sig1_lo < 0.  # compressive force in direction of principle stress at lower side

        # compare imposed compressive normal force with 'n_Rdc' as obtained from compression test
        #
        bool_arr = cond_nsu_lt_0
        eta_n_up[bool_arr] = self.n_sig1_up[bool_arr] / n_Rdc[bool_arr]

        bool_arr = cond_nsl_lt_0
        eta_n_lo[bool_arr] = self.n_sig1_lo[bool_arr] / n_Rdc[bool_arr]

        # get 'eta_m' based on imposed moment compared with moment resistence
        # NOTE: use a linear increase factor for resistance moment based on reference thickness (= minimum thickness)
        #
        min_thickness = np.min(self.thickness)
        eta_m_lo = self.m_sig1_lo / (m_Rd_lo_1 * self.thickness / min_thickness)
        eta_m_up = self.m_sig1_up / (m_Rd_up_1 * self.thickness / min_thickness)

        # get total 'eta_mn' based on imposed normal force and moment
        # NOTE: if eta_n is negative (caused by a compressive normal force) take the absolute value
        # NOTE: if eta_m is negative (caused by a negative moment) take the absolute value
        #
        eta_nm_lo = np.abs(eta_n_lo) + np.abs(eta_m_lo)
        eta_nm_up = np.abs(eta_n_up) + np.abs(eta_m_up)

        # get maximum 'eta_mn' of both principle directions of upper and lower side
        #
        eta_nm1_tot = ndmax(hstack([ eta_nm_up, eta_nm_lo]), axis=1)[:, None]

        #---------------
        # 2-direction:
        #---------------

        # initialize arrays to be filled based on case distinction
        #
        eta_n2_up = np.zeros_like(self.n_sig2_up)
        eta_n2_lo = np.zeros_like(self.n_sig2_lo)

        # cases with a tensile normal force
        #
        cond_ns2u_ge_0 = self.n_sig2_up >= 0.  # tensile force in direction of principle stress at upper side
        cond_ns2l_ge_0 = self.n_sig2_lo >= 0.  # tensile force in direction of principle stress at lower side

        # compare imposed tensile normal force with 'n_Rd,t' as obtained from tensile test
        #
        bool_arr = cond_ns2u_ge_0
        eta_n2_up[bool_arr] = self.n_sig2_up[bool_arr] / n_Rdt_up_2[bool_arr]

        bool_arr = cond_ns2l_ge_0
        eta_n2_lo[bool_arr] = self.n_sig2_lo[bool_arr] / n_Rdt_lo_2[bool_arr]

        # cases with a compressive normal force
        #
        cond_ns2u_lt_0 = self.n_sig2_up < 0.  # compressive force in direction of principle stress at upper side
        cond_ns2l_lt_0 = self.n_sig2_lo < 0.  # compressive force in direction of principle stress at lower side

        # compare imposed compressive normal force with 'n_Rdc' as obtained from compression test
        #
        bool_arr = cond_ns2u_lt_0
        eta_n2_up[bool_arr] = self.n_sig2_up[bool_arr] / n_Rdc[bool_arr]

        bool_arr = cond_ns2l_lt_0
        eta_n2_lo[bool_arr] = self.n_sig2_lo[bool_arr] / n_Rdc[bool_arr]

        # get 'eta_m' based on imposed moment compared with moment resistence
        # NOTE: use a linear increase factor for resistance moment based on reference thickness (= minimum thickness)
        #
        eta_m2_lo = self.m_sig2_lo / (m_Rd_lo_2 * self.thickness / min_thickness)
        eta_m2_up = self.m_sig2_up / (m_Rd_up_2 * self.thickness / min_thickness)

        # get total 'eta_mn' based on imposed normal force and moment
        # NOTE: if eta_n is negative (caused by a compressive normal force) take the absolute value
        # NOTE: if eta_m is negative (caused by a negative moment) take the absolute value
        #
        eta_nm2_lo = np.abs(eta_n2_lo) + np.abs(eta_m2_lo)
        eta_nm2_up = np.abs(eta_n2_up) + np.abs(eta_m2_up)

        # get maximum 'eta_mn' of both principle directions of upper and lower side
        #
        eta_nm2_tot = ndmax(hstack([ eta_nm2_up, eta_nm2_lo]), axis=1)[:, None]

        # overall maximum eta_nm for 1st and 2nd principle direction
        #
        eta_nm_tot = ndmax(hstack([ eta_nm1_tot, eta_nm2_tot]), axis=1)[:, None]

        # overall maximum eta_n and eta_m distinguishing normal forces and bending moment influence:
        #
        eta_n_tot = ndmax(hstack([ eta_n_lo, eta_n2_lo, eta_n_up, eta_n2_up]), axis=1)[:, None]
        eta_m_tot = ndmax(hstack([ eta_m_lo, eta_m2_lo, eta_m_up, eta_m2_up]), axis=1)[:, None]

        #------------------------------------------------------------
        # construct a dictionary containing the return values
        #------------------------------------------------------------
        return {
                 'beta_up_deg':beta_up_deg,
                 'beta_lo_deg':beta_lo_deg,

                 'n_Rdt_up_1':n_Rdt_up_1,
                 'n_Rdt_lo_1':n_Rdt_lo_1,
                 'm_Rd_up_1':m_Rd_up_1,
                 'm_Rd_lo_1':m_Rd_lo_1,

                 'n_Rdt_up_2':n_Rdt_up_2,
                 'n_Rdt_lo_2':n_Rdt_lo_2,
                 'm_Rd_up_2':m_Rd_up_2,
                 'm_Rd_lo_2':m_Rd_lo_2,

                 'eta_n_up':eta_n_up,
                 'eta_m_up':eta_m_up,
                 'eta_nm_up':eta_nm_up,
                 'eta_n_lo':eta_n_lo,
                 'eta_m_lo':eta_m_lo,
                 'eta_nm_lo':eta_nm_lo,
                 'eta_nm_tot':eta_nm_tot,

                 'eta_n_tot':eta_n_tot,
                 'eta_m_tot':eta_m_tot,

                 'eta_n2_up':eta_n2_up,
                 'eta_m2_up':eta_m2_up,
                 'eta_nm2_up':eta_nm2_up,
                 'eta_n2_lo':eta_n2_lo,
                 'eta_m2_lo':eta_m2_lo,
                 'eta_nm2_lo':eta_nm2_lo,

                 'k_alpha_lo' : k_alpha_lo,
                 'k_alpha_up' : k_alpha_up
                 }

    #-----------------------------------------------
    # LS_COLUMNS: specify the properties that are displayed in the view
    #-----------------------------------------------

    # choose the assess parameter used for sorting
    # defined by the property name 'assess_name'
    #
    assess_name = 'max_eta_nm_tot'

    max_eta_nm_tot = Property(depends_on='+input')
    @cached_property
    def _get_max_eta_nm_tot(self):
        return ndmax(self.eta_nm_tot)

    ls_columns = List([
                       'alpha_sig1_up_deg',
                       'alpha_sig2_up_deg',
                       'alpha_sig1_lo_deg',
                       'alpha_sig2_lo_deg',

                       'beta_up_deg', 'k_alpha_up',
                       'beta_lo_deg', 'k_alpha_lo',

                       'n_Rdt_up_1',
                       'n_Rdt_lo_1',
                       'm_Rd_up_1',
                       'm_Rd_lo_1',

#                       'n_Rdt_up_2',
#                       'n_Rdt_lo_2',
#                       'm_Rd_up_2',
#                       'm_Rd_lo_2',

                       'eta_n_up', 'eta_m_up', 'eta_nm_up',
                       'eta_n_lo', 'eta_m_lo', 'eta_nm_lo',
                       'eta_n2_up', 'eta_m2_up', 'eta_nm2_up',
                       'eta_n2_lo', 'eta_m2_lo', 'eta_nm2_lo',
                       'eta_nm_tot'])

    # specify the material properties for the view:
    #
    plot_item_mpl = Item(name='n_0_Rdt', label='normal tensile strength [kN/m]:  n_0_Rdt ', style='readonly', format_str="%.1f"), \
                    Item(name='n_Rdc', label='normal compressive strength [kN/m]:  n_0_Rdc ', style='readonly', format_str="%.1f"), \
                    Item(name='m_0_Rd', label='bending strength [kNm/m]:  m_0_Rd ', style='readonly', format_str="%.1f")

    plot_item_mpt = Item(name='n_90_Rdt', label='normal tensile strength [kN/m]:  n_90_Rd ', style='readonly', format_str="%.1f"), \
                    Item(name='n_Rdc', label='normal compressive strength [kN/m]:  n_0_Rdc ', style='readonly', format_str="%.1f"), \
                    Item(name='m_90_Rd', label='bending strength [kNm/m]:  m_90_Rd ', style='readonly', format_str="%.1f")

    beta_up_deg = Property(Array)
    def _get_beta_up_deg(self):
        return self.ls_values['beta_up_deg']

    beta_lo_deg = Property(Array)
    def _get_beta_lo_deg(self):
        return self.ls_values['beta_lo_deg']

    #------------------------------------
    # eval == 'eta_nm'
    #------------------------------------

    n_Rdt_up_1 = Property(Array)
    def _get_n_Rdt_up_1(self):
        return self.ls_values['n_Rdt_up_1']

    n_Rdt_lo_1 = Property(Array)
    def _get_n_Rdt_lo_1(self):
        return self.ls_values['n_Rdt_lo_1']

    m_Rd_up_1 = Property(Array)
    def _get_m_Rd_up_1(self):
        return self.ls_values['m_Rd_up_1']

    m_Rd_lo_1 = Property(Array)
    def _get_m_Rd_lo_1(self):
        return self.ls_values['m_Rd_lo_1']

    n_Rdt_up_2 = Property(Array)
    def _get_n_Rdt_up_2(self):
        return self.ls_values['n_Rdt_up_2']

    n_Rdt_lo_2 = Property(Array)
    def _get_n_Rdt_lo_2(self):
        return self.ls_values['n_Rdt_lo_2']

    m_Rd_up_2 = Property(Array)
    def _get_m_Rd_up_2(self):
        return self.ls_values['m_Rd_up_2']

    m_Rd_lo_2 = Property(Array)
    def _get_m_Rd_lo_2(self):
        return self.ls_values['m_Rd_lo_2']

    k_alpha_up = Property(Array)
    def _get_k_alpha_up(self):
        return self.ls_values['k_alpha_up']

    k_alpha_lo = Property(Array)
    def _get_k_alpha_lo(self):
        return self.ls_values['k_alpha_lo']

    # 1st-principle direction
    #
    eta_n_up = Property(Array)
    def _get_eta_n_up(self):
        return self.ls_values['eta_n_up']

    eta_m_up = Property(Array)
    def _get_eta_m_up(self):
        return self.ls_values['eta_m_up']

    eta_nm_up = Property(Array)
    def _get_eta_nm_up(self):
        return self.ls_values['eta_nm_up']

    eta_n_lo = Property(Array)
    def _get_eta_n_lo(self):
        return self.ls_values['eta_n_lo']

    eta_m_lo = Property(Array)
    def _get_eta_m_lo(self):
        return self.ls_values['eta_m_lo']

    eta_nm_lo = Property(Array)
    def _get_eta_nm_lo(self):
        return self.ls_values['eta_nm_lo']

    # max from 1st and 2nd-principle direction
    #
    eta_nm_tot = Property(Array)
    def _get_eta_nm_tot(self):
        return self.ls_values['eta_nm_tot']

    # max from 1st and 2nd-principle direction only from normal forces
    #
    eta_n_tot = Property(Array)
    def _get_eta_n_tot(self):
        return self.ls_values['eta_n_tot']

    # max from 1st and 2nd-principle direction only from bending moments
    #
    eta_m_tot = Property(Array)
    def _get_eta_m_tot(self):
        return self.ls_values['eta_m_tot']

    # 2nd-principle direction
    #
    eta_n2_up = Property(Array)
    def _get_eta_n2_up(self):
        return self.ls_values['eta_n2_up']

    eta_m2_up = Property(Array)
    def _get_eta_m2_up(self):
        return self.ls_values['eta_m2_up']

    eta_nm2_up = Property(Array)
    def _get_eta_nm2_up(self):
        return self.ls_values['eta_nm2_up']

    eta_n2_lo = Property(Array)
    def _get_eta_n2_lo(self):
        return self.ls_values['eta_n2_lo']

    eta_m2_lo = Property(Array)
    def _get_eta_m2_lo(self):
        return self.ls_values['eta_m2_lo']

    eta_nm2_lo = Property(Array)
    def _get_eta_nm2_lo(self):
        return self.ls_values['eta_nm2_lo']

    #-------------------------------
    # ls view
    #-------------------------------
    # @todo: the dynamic selection of the columns to be displayed
    # does not work in connection with the LSArrayAdapter
    traits_view = View(
                       VGroup(
                        HGroup(
                            VGroup(
                                plot_item_mpl,
                                label='material Properties (longitudinal)'
                                  ),
                            VGroup(
                                plot_item_mpt,
                                label='material Properties (transversal)'
                                  ),
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

LSLIST = [ SLS, ULS ]

class LSTable(HasTraits):
    '''Assessment tool
    '''

    is_id = Int(0)

    # geo data: coordinates and element thickness
    #
    geo_data = Dict

    #------------------------------------------------------------
    # evaluation with conservative simplification for 'k_alpha'
    #------------------------------------------------------------
    # if flag is set to 'True' the resistance values 'n_Rdt' and 'm_Rd' below are
    # multiplied with the highest reduction factor 'k_alpha = 0.707', independently
    # of the true deflection angle 'beta_q' and 'beta_l'
    #
    k_alpha_min = Bool

    # specify the strength characteristics of the material;
    # necessary parameter in order to pass values from 'LCCTableULS' to 'ls_table' and as WeekRef to 'ULS';
    # this gives the possibility to specify strength values within the call of 'LCCTableULS';
    #
    strength_characteristics_dict = Dict

    #-----------------------------------
    # geo data
    #-----------------------------------
    elem_no = Property(Array)
    def _get_elem_no(self):
        return self.geo_data['elem_no']

    X = Property(Array)
    def _get_X(self):
        return self.geo_data['X']

    Y = Property(Array)
    def _get_Y(self):
        return self.geo_data['Y']

    Z = Property(Array)
    def _get_Z(self):
        return self.geo_data['Z']

    thickness = Property(Array)
    def _get_thickness(self):
        '''element thickness [m])'''
        return self.geo_data['thickness']

    #-----------------------------------
    # derived geometric values
    #-----------------------------------
    # area
    #
    A = Property(Array)
    def _get_A(self):
        return self.thickness * 1.

    # moment of inertia
    #
    W = Property(Array)
    def _get_W(self):
        return 1. * self.thickness ** 2 / 6.

    #-----------------------------------
    # state data: stress resultants
    #-----------------------------------
    #
    state_data = Dict

    mx = Property(Array)
    def _get_mx(self):
        return self.state_data['mx']

    my = Property(Array)
    def _get_my(self):
        return self.state_data['my']

    mxy = Property(Array)
    def _get_mxy(self):
        return self.state_data['mxy']

    nx = Property(Array)
    def _get_nx(self):
        return self.state_data['nx']

    ny = Property(Array)
    def _get_ny(self):
        return self.state_data['ny']

    nxy = Property(Array)
    def _get_nxy(self):
        return self.state_data['nxy']

    check_consistency_alpha = Bool(True)

    # ------------------------------------------------------------
    # Index sig: calculate principle direction of the stresses at
    # the lower and upper side and get the corresponding values of
    # the stresses at the opposite side. Also get the corresponding
    # values of the normal force and the moment in this direction
    # ------------------------------------------------------------

    princ_values_sig = Property(Dict, depends_on='data_file_stress_resultants')
    @cached_property
    def _get_princ_values_sig(self):
        '''principle value of the stresses for the lower ('lo') and upper ('up') face:
        '''
        # stress_resultants in global coordinates
        # --> moments in kNm
        # --> normal forces in kN
        # --> area in m**2
        # --> resisting moment in m**3
        # --> stresses in MPa
        #
        mx = self.mx
        my = self.my
        mxy = self.mxy
        nx = self.nx
        ny = self.ny
        nxy = self.nxy

        # geometrical properties:
        #
        A = self.A
        W = self.W

        # compare the formulae with the RFEM-manual p.290

        # stresses [MPa] upper face in global direction:
        #
        sigx_up = (nx / A - mx / W) / 1000.
        sigy_up = (ny / A - my / W) / 1000.
        sigxy_up = (nxy / A - mxy / W) / 1000.

        # stresses [MPa] lower face in global direction:
        #
        sigx_lo = (nx / A + mx / W) / 1000.
        sigy_lo = (ny / A + my / W) / 1000.
        sigxy_lo = (nxy / A + mxy / W) / 1000.

        #--------------
        # upper face:
        #--------------

        # principal stresses lower face:
        #
        sig1_up = 0.5 * (sigx_up + sigy_up) + 0.5 * sqrt((sigx_up - sigy_up) ** 2 + 4 * sigxy_up ** 2)
        sig2_up = 0.5 * (sigx_up + sigy_up) - 0.5 * sqrt((sigx_up - sigy_up) ** 2 + 4 * sigxy_up ** 2)

        # formula corresponds to mechanic formula with positive z-rotation for positive 'sigxy'
        # NOTE: underline score in the name indicates that the angle is still unsorted
        # formula returns an angle that lies within -/+90deg
        #
        alpha_sig1_up = pi / 2. * ones_like(sig1_up)
        bool_arr = sigx_up != sig2_up
        alpha_sig1_up[ bool_arr ] = arctan(sigxy_up[ bool_arr ] / (sigx_up[ bool_arr ] - sig2_up[ bool_arr ]))

        # angle of principle stresses (2-direction = minimum stresses (compression))
        # NOTE: formula guarantees that the 2nd principle angle returned for the first or second quadrant,
        #       i.e. abs(alpha_sig2_up_) <= pi/2
        #       i.e. if alpha_sig1_up_ lies in the first quadrant give second angle for second quadrant
        #            if alpha_sig1_up lies in the second quadrant give angle for first quadrant
        #            if alpha_sig1_up == 0. 'alpha_sig2_up' is set to pi/2
        #
        alpha_sig2_up = pi / 2. * ones_like(sig1_up)
        bool_arr = np.sign(alpha_sig1_up) != 0.
        alpha_sig2_up[bool_arr] = -pi / 2 * np.sign(alpha_sig1_up[bool_arr]) + alpha_sig1_up[bool_arr]

        # convert units from radiant to degree
        #
        alpha_sig1_up_deg = alpha_sig1_up * 180. / pi
        alpha_sig2_up_deg = alpha_sig2_up * 180. / pi

        # transform moments and normal forces in the direction of the principal stresses (1-direction)
        #
        m_sig1_up = 0.5 * (mx + my) + 0.5 * (mx - my) * cos(2 * alpha_sig1_up) + mxy * sin(2 * alpha_sig1_up)
        n_sig1_up = 0.5 * (nx + ny) + 0.5 * (nx - ny) * cos(2 * alpha_sig1_up) + nxy * sin(2 * alpha_sig1_up)

        # transform moments and normal forces in the direction of the principal stresses (2-direction)
        #
        m_sig2_up = 0.5 * (mx + my) - 0.5 * (mx - my) * cos(2 * alpha_sig1_up) - mxy * sin(2 * alpha_sig1_up)
        n_sig2_up = 0.5 * (nx + ny) - 0.5 * (nx - ny) * cos(2 * alpha_sig1_up) - nxy * sin(2 * alpha_sig1_up)

        #--------------
        # lower face:
        #--------------

        # principal stresses lower face:
        #
        sig1_lo = 0.5 * (sigx_lo + sigy_lo) + 0.5 * sqrt((sigx_lo - sigy_lo) ** 2 + 4 * sigxy_lo ** 2)
        sig2_lo = 0.5 * (sigx_lo + sigy_lo) - 0.5 * sqrt((sigx_lo - sigy_lo) ** 2 + 4 * sigxy_lo ** 2)

        # formula corresponds to mechanic formula with positive z-rotation for positive 'sigxy'
        # NOTE: underline score in the name indicates that the angle is still unsorted
        # formula returns an angle that lies within -/+90deg
        #
        alpha_sig1_lo = pi / 2. * ones_like(sig1_lo)
        bool_arr = sigx_lo != sig2_lo
        alpha_sig1_lo[ bool_arr ] = arctan(sigxy_lo[ bool_arr ] / (sigx_lo[ bool_arr ] - sig2_lo[ bool_arr ]))

        # angle of principle stresses (2-direction = minimum stresses (compression))
        # NOTE: formula guarantees that the 2nd principle angle returned for the first or second quadrant,
        #       i.e. abs(alpha_sig2_lo_) <= pi/2
        #       i.e. if alpha_sig1_lo_ lies in the first quadrant give second angle for second quadrant
        #            if alpha_sig1_lo lies in the second quadrant give angle for first quadrant
        #            if alpha_sig1_lo == 0. 'alpha_sig2_lo' is set to pi/2
        #
        alpha_sig2_lo = pi / 2. * ones_like(sig1_lo)
        bool_arr = np.sign(alpha_sig1_lo) != 0.
        alpha_sig2_lo[bool_arr] = -pi / 2 * np.sign(alpha_sig1_lo[bool_arr]) + alpha_sig1_lo[bool_arr]

        # convert units from radiant to degree
        #
        alpha_sig1_lo_deg = alpha_sig1_lo * 180. / pi
        alpha_sig2_lo_deg = alpha_sig2_lo * 180. / pi

        # transform moments and normal forces in the direction of the principal stresses (1-direction)
        #
        m_sig1_lo = 0.5 * (mx + my) + 0.5 * (mx - my) * cos(2 * alpha_sig1_lo) + mxy * sin(2 * alpha_sig1_lo)
        n_sig1_lo = 0.5 * (nx + ny) + 0.5 * (nx - ny) * cos(2 * alpha_sig1_lo) + nxy * sin(2 * alpha_sig1_lo)

        # transform moments and normal forces in the direction of the principal stresses (2-direction)
        #
        m_sig2_lo = 0.5 * (mx + my) - 0.5 * (mx - my) * cos(2 * alpha_sig1_lo) - mxy * sin(2 * alpha_sig1_lo)
        n_sig2_lo = 0.5 * (nx + ny) - 0.5 * (nx - ny) * cos(2 * alpha_sig1_lo) - nxy * sin(2 * alpha_sig1_lo)

        return {
                 'sigx_up' : sigx_up, 'sigy_up' : sigy_up, 'sigxy_up' : sigxy_up,
                 'sig1_up' : sig1_up, 'sig2_up' : sig2_up,
                 'alpha_sig1_up' : alpha_sig1_up, 'alpha_sig2_up' : alpha_sig2_up,
                 'alpha_sig1_up_deg' : alpha_sig1_up_deg, 'alpha_sig2_up_deg' : alpha_sig2_up_deg,
                 'm_sig1_up' : m_sig1_up, 'n_sig1_up' : n_sig1_up,
                 'm_sig2_up' : m_sig2_up, 'n_sig2_up' : n_sig2_up,

                 'sigx_lo' : sigx_lo, 'sigy_lo' : sigy_lo, 'sigxy_lo' : sigxy_lo,
                 'sig1_lo' : sig1_lo, 'sig2_lo' : sig2_lo,
                 'alpha_sig1_lo' : alpha_sig1_lo, 'alpha_sig2_lo' : alpha_sig2_lo,
                 'alpha_sig1_lo_deg' : alpha_sig1_lo_deg, 'alpha_sig2_lo_deg' : alpha_sig2_lo_deg,
                 'm_sig1_lo' : m_sig1_lo, 'n_sig1_lo' : n_sig1_lo,
                 'm_sig2_lo' : m_sig2_lo, 'n_sig2_lo' : n_sig2_lo,
                }

    # stresses upper face:
    #
    sigx_up = Property(Float)
    def _get_sigx_up(self):
        return self.princ_values_sig['sigx_up']

    sigy_up = Property(Float)
    def _get_sigy_up(self):
        return self.princ_values_sig['sigy_up']

    sigxy_up = Property(Float)
    def _get_sigxy_up(self):
        return self.princ_values_sig['sigxy_up']

    sig1_up = Property(Float)
    def _get_sig1_up(self):
        return self.princ_values_sig['sig1_up']

    sig2_up = Property(Float)
    def _get_sig2_up(self):
        return self.princ_values_sig['sig2_up']

    alpha_sig1_up = Property(Float)
    def _get_alpha_sig1_up(self):
        return self.princ_values_sig['alpha_sig1_up']

    alpha_sig2_up = Property(Float)
    def _get_alpha_sig2_up(self):
        return self.princ_values_sig['alpha_sig2_up']

    alpha_sig1_up_deg = Property(Float)
    def _get_alpha_sig1_up_deg(self):
        return self.princ_values_sig['alpha_sig1_up_deg']

    alpha_sig2_up_deg = Property(Float)
    def _get_alpha_sig2_up_deg(self):
        return self.princ_values_sig['alpha_sig2_up_deg']

    m_sig1_up = Property(Float)
    def _get_m_sig1_up(self):
        return self.princ_values_sig['m_sig1_up']

    n_sig1_up = Property(Float)
    def _get_n_sig1_up(self):
        return self.princ_values_sig['n_sig1_up']

    m_sig2_up = Property(Float)
    def _get_m_sig2_up(self):
        return self.princ_values_sig['m_sig2_up']

    n_sig2_up = Property(Float)
    def _get_n_sig2_up(self):
        return self.princ_values_sig['n_sig2_up']

    # stresses lower face:
    #
    sigx_lo = Property(Float)
    def _get_sigx_lo(self):
        return self.princ_values_sig['sigx_lo']

    sigy_lo = Property(Float)
    def _get_sigy_lo(self):
        return self.princ_values_sig['sigy_lo']

    sigxy_lo = Property(Float)
    def _get_sigxy_lo(self):
        return self.princ_values_sig['sigxy_lo']

    sig1_lo = Property(Float)
    def _get_sig1_lo(self):
        return self.princ_values_sig['sig1_lo']

    sig2_lo = Property(Float)
    def _get_sig2_lo(self):
        return self.princ_values_sig['sig2_lo']

    alpha_sig1_lo = Property(Float)
    def _get_alpha_sig1_lo(self):
        return self.princ_values_sig['alpha_sig1_lo']

    alpha_sig2_lo = Property(Float)
    def _get_alpha_sig2_lo(self):
        return self.princ_values_sig['alpha_sig2_lo']

    alpha_sig1_lo_deg = Property(Float)
    def _get_alpha_sig1_lo_deg(self):
        return self.princ_values_sig['alpha_sig1_lo_deg']

    alpha_sig2_lo_deg = Property(Float)
    def _get_alpha_sig2_lo_deg(self):
        return self.princ_values_sig['alpha_sig2_lo_deg']

    m_sig1_lo = Property(Float)
    def _get_m_sig1_lo(self):
        return self.princ_values_sig['m_sig1_lo']

    n_sig1_lo = Property(Float)
    def _get_n_sig1_lo(self):
        return self.princ_values_sig['n_sig1_lo']

    m_sig2_lo = Property(Float)
    def _get_m_sig2_lo(self):
        return self.princ_values_sig['m_sig2_lo']

    n_sig2_lo = Property(Float)
    def _get_n_sig2_lo(self):
        return self.princ_values_sig['n_sig2_lo']

    #------------------------------------------
    # combinations of limit states, stress resultants and directions
    #------------------------------------------

    ls = Trait('ULS',
                {'ULS' : ULS,
                 'SLS' : SLS })

    ls_class = Instance(LS)
    def _ls_class_default(self):
        '''ls instances, e.g. ULS()
        '''
        ls_class = self.ls_
        return ls_class(ls_table=self)

    assess_value = Property
    def _get_assess_value(self):
        ls = self.ls_class
        return getattr(ls, ls.assess_name)

    traits_view = View(Tabbed(
                            Item('ls_class@' , label="ls", show_label=False),
                            scrollable=False,
                         ),
                      resizable=True,
                      scrollable=True,
                      height=1000,
                      width=1100
                      )
