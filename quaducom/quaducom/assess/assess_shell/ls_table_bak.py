'''
Created on Jun 23, 2010

@author: alexander
'''

from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, Trait

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from enthought.mayavi import \
    mlab

#from enthought.mayavi.mlab import \
#    colorbar, show, points3d
#
#from enthought.mayavi.api import \
#    Engine


from enthought.traits.ui.table_column import \
    ObjectColumn

from enthought.traits.ui.menu import \
    OKButton, CancelButton

from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, sqrt, frompyfunc, \
                max as ndmax

from math import pi
from string import split
import os

from scipy.io import read_array


DIRLIST = ['x', 'y']
SRLIST = ['M', 'N']

class LSArrayAdapter ( TabularAdapter ):

    columns = Property
    def _get_columns( self ):
#        print 'GETTING COLUMNS', self.object.columns, self.object, self.object.__class__
        columns = self.object.columns
        return [ ( name, idx ) for idx, name in enumerate( columns ) ]

    font = 'Courier 10'
    alignment = 'right'
    format = '%5.2f'#'%g'
    even_bg_color = Color( 0xE0E0FF )
    width = Float( 80 )

    #@todo: format columns using 'column_id'
#    adapter_column_map = Property(depends_on = 'adapters,columns')


class LS( HasTraits ):
    '''Limit state class
    '''

    # backward link to the info shell to access the
    # input data when calculating 
    # the limit-state-specific values
    #
    ls_table = WeakRef

    # parameters of the limit state
    #
    dir = Enum( DIRLIST )
    stress_res = Enum( SRLIST )

    #-------------------------------
    # ls columns
    #-------------------------------
    # defined in the subclasses
    #
    ls_columns = List
    show_ls_columns = Bool( True )

    #-------------------------------
    # sr columns
    #-------------------------------

    # stress resultant columns - for ULS this is defined in the subclasses
    #
    sr_columns = List( ['m', 'n'] )
    show_sr_columns = Bool( True )

    # stress resultant columns - generated from the parameter combination
    # dir and stress_res - one of MX, NX, MY, NY
    #
    m_varname = Property( Str )
    def _get_m_varname( self ):
        # e.g. mx_N 
        appendix = self.dir + '_' + self.stress_res
        return 'm' + appendix

    n_varname = Property( Str )
    def _get_n_varname( self ):
        # e.g. nx_N 
        appendix = self.dir + '_' + self.stress_res
        return 'n' + appendix

    n = Property( Float )
    def _get_n( self ):
        return getattr( self.ls_table, self.n_varname )

    m = Property( Float )
    def _get_m( self ):
        return getattr( self.ls_table, self.m_varname )

    #-------------------------------
    # geo columns form info shell
    #-------------------------------

    geo_columns = List( [ 'elem_no', 'X', 'Y', 'Z', 'D_elem' ] )
    show_geo_columns = Bool( True )

    elem_no = Property( Float )
    def _get_elem_no( self ):
        return self.ls_table.elem_no

    X = Property( Float )
    def _get_X( self ):
        return self.ls_table.X

    Y = Property( Float )
    def _get_Y( self ):
        return self.ls_table.Y

    Z = Property( Float )
    def _get_Z( self ):
        return self.ls_table.Z

    D_elem = Property( Float )
    def _get_D_elem( self ):
        return self.ls_table.D_elem

    #-------------------------------
    # state columns form info shell
    #-------------------------------

#    state_columns = List( ['mx', 'my', 'mxy', 'nx', 'ny', 'nxy' ] )
    state_columns = List( ['mx', 'my', 'mxy', 'nx', 'ny', 'nxy',
                           'sigx_lo', 'sigy_lo', 'sigxy_lo', 'sig1_lo', 'sig2_lo', 'alpha_sig_lo',
                           'sigx_up', 'sigy_up', 'sigxy_up', 'sig1_up', 'sig2_up', 'alpha_sig_up',
                            ] )

    show_state_columns = Bool( True )

    mx = Property( Float )
    def _get_mx( self ):
        return self.ls_table.mx

    my = Property( Float )
    def _get_my( self ):
        return self.ls_table.my

    mxy = Property( Float )
    def _get_mxy( self ):
        return self.ls_table.mxy

    nx = Property( Float )
    def _get_nx( self ):
        return self.ls_table.nx

    ny = Property( Float )
    def _get_ny( self ):
        return self.ls_table.ny

    nxy = Property( Float )
    def _get_nxy( self ):
        return self.ls_table.nxy


    # evaluate principal stresses
    # upper face:
    #
    sigx_up = Property( Float )
    def _get_sigx_up( self ):
        return self.ls_table.sigx_up

    sigy_up = Property( Float )
    def _get_sigy_up( self ):
        return self.ls_table.sigy_up

    sigxy_up = Property( Float )
    def _get_sigxy_up( self ):
        return self.ls_table.sigxy_up

    sig1_up = Property( Float )
    def _get_sig1_up( self ):
        return self.ls_table.sig1_up

    sig2_up = Property( Float )
    def _get_sig2_up( self ):
        return self.ls_table.sig2_up

    alpha_sig_up = Property( Float )
    def _get_alpha_sig_up( self ):
        return self.ls_table.alpha_sig_up

    # lower face:
    #
    sigx_lo = Property( Float )
    def _get_sigx_lo( self ):
        return self.ls_table.sigx_lo

    sigy_lo = Property( Float )
    def _get_sigy_lo( self ):
        return self.ls_table.sigy_lo

    sigxy_lo = Property( Float )
    def _get_sigxy_lo( self ):
        return self.ls_table.sigxy_lo

    sig1_lo = Property( Float )
    def _get_sig1_lo( self ):
        return self.ls_table.sig1_lo

    sig2_lo = Property( Float )
    def _get_sig2_lo( self ):
        return self.ls_table.sig2_lo

    alpha_sig_lo = Property( Float )
    def _get_alpha_sig_lo( self ):
        return self.ls_table.alpha_sig_lo

    #-------------------------------
    # ls table
    #-------------------------------

    # all columns associated with the limit state including the corresponding
    # stress resultants
    #
    columns = Property( List, depends_on = 'show_geo_columns, show_state_columns,\
                                             show_sr_columns, show_ls_columns' )
    @cached_property
    def _get_columns( self ):
        columns = []

        if self.show_geo_columns:
            columns += self.geo_columns

        if self.show_state_columns:
            columns += self.state_columns

        if self.show_sr_columns:
            columns += self.sr_columns

        if self.show_ls_columns:
            columns += self.ls_columns

        return columns

    # select column used for sorting the data in selected sorting order 
    #
    sort_column = Enum( values = 'columns' )
    def _sort_column_default( self ):
        return self.columns[-1]

    sort_order = Enum( 'descending', 'ascending', 'unsorted' )

    #-------------------------------------------------------
    # get the maximum value of the selected variable 
    # 'max_in_column' of the current sheet (only one sheet)
    #-------------------------------------------------------

    # get the maximum value of the chosen column
    #
    max_in_column = Enum( values = 'columns' )
    def _max_in_column_default( self ):
        return self.columns[-1]

    max_value = Property( depends_on = 'max_in_column' )
    def _get_max_value( self ):
        col = getattr( self, self.max_in_column )[:, 0]
        return max( col )

    #-------------------------------------------------------
    # get the maximum value and the corresponding case of 
    # the selected variable 'max_in_column' in all (!) sheets
    #-------------------------------------------------------

    max_value_all = Property( depends_on = 'max_in_column' )
    def _get_max_value_all( self ):
        return self.ls_table.max_value_and_case[ self.max_in_column ]['max_value']

    max_case = Property( depends_on = 'max_in_column' )
    def _get_max_case( self ):
        return self.ls_table.max_value_and_case[ self.max_in_column ]['max_case']

    #-------------------------------------------------------
    # get ls_table for View
    #-------------------------------------------------------

    # stack columns together for table used by TabularEditor
    #
    ls_array = Property( Array, depends_on = 'sort_column, sort_order, show_geo_columns, \
                                              show_state_columns, show_sr_columns, show_ls_columns' )
    @cached_property
    def _get_ls_array( self ):

        arr_list = [ getattr( self, col ) for col in self.columns ]

        # get the array currently selected by the sort_column enumeration
        #
        sort_arr = getattr( self, self.sort_column )[:, 0]
        sort_idx = argsort( sort_arr )
        ls_array = hstack( arr_list )

        if self.sort_order == 'descending':
            return ls_array[ sort_idx[::-1] ]
        if self.sort_order == 'ascending':
            return ls_array[ sort_idx ]
        if self.sort_order == 'unsorted':
            return ls_array

    #---------------------------------
    # plot outputs in mlab-window 
    #---------------------------------

    plot_column = Enum( values = 'columns' )
    plot = Button
    def _plot_fired( self ):
        X = self.ls_table.X[:, 0]
        Y = self.ls_table.Y[:, 0]
        Z = self.ls_table.Z[:, 0]
        plot_col = getattr( self, self.plot_column )[:, 0]

        if self.plot_column == 'n_tex':
            plot_col = where( plot_col < 0, 0, plot_col )

        mlab.figure( figure = "SFB532Demo",
                     bgcolor = ( 1.0, 1.0, 1.0 ),
                     fgcolor = ( 0.0, 0.0, 0.0 ) )

        mlab.points3d( X, Y, ( -1.0 ) * Z, plot_col,
#                       colormap = "gist_rainbow",
#                       colormap = "Reds",
                       colormap = "YlOrBr",
                       mode = "cube",
                       scale_factor = 0.15 )

        mlab.scalarbar( title = self.plot_column, orientation = 'vertical' )

        mlab.show



    # name of the trait that is used to assess the evaluated design
    #
    assess_name = Str( '' )

    #-------------------------------
    # ls group
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    ls_group = VGroup( 
                        HGroup( #Item( 'assess_name' ),
                                Item( 'max_in_column' ),
                                Item( 'max_value', style = 'readonly', format_str = '%6.2f' ),
                                Item( 'max_value_all', style = 'readonly', format_str = '%6.2f' ),
                                Item( 'max_case', style = 'readonly', label = 'found in case: ' ),
                              ),
                        HGroup( Item( 'sort_column' ),
                                Item( 'sort_order' ),
                                Item( 'show_geo_columns', label = 'show geo' ),
                                Item( 'show_state_columns', label = 'show state' ),
                                Item( 'show_sr_columns', label = 'show sr' ),
                                Item( 'plot_column' ),
                                Item( 'plot' ),
                              ),
                     )


class SLS( LS ):
    '''Serviceability limit state
    '''

    # ------------------------------------------------------------
    # SLS: material parameters (Inputs)
    # ------------------------------------------------------------

    # tensile strength [MPa]
    f_ctk = Float( 4.0, input = True )

    # flexural tensile strength [MPa]
    f_m = Float( 5.0, input = True )

    # ------------------------------------------------------------
    # SLS - derived params:
    # ------------------------------------------------------------

    # area
    #
    A = Property( Float )
    def _get_A( self ):
        return self.ls_table.D_elem * 1.

    # moment of inertia
    #
    W = Property( Float )
    def _get_W( self ):
        return 1. * self.ls_table.D_elem ** 2 / 6.

    # ------------------------------------------------------------
    # SLS: outputs
    # ------------------------------------------------------------

    ls_columns = List( ['sig_n', 'sig_m', 'eta_n', 'eta_m', 'eta_tot', ] )

    ls_values = Property( depends_on = '+input' )
    @cached_property
    def _get_ls_values( self ):
        '''get the outputs for SLS
        '''
        n = self.n
        m = self.m

        A = self.A
        W = self.W
        f_ctk = self.f_ctk
        f_m = self.f_m

        sig_n = n / A / 1000.
        sig_m = abs( m / W ) / 1000.
        eta_n = sig_n / f_ctk
        eta_m = sig_m / f_m
        eta_tot = eta_n + eta_m

        return { 'sig_n':sig_n, 'sig_m':sig_m,
                 'eta_n':eta_n, 'eta_m':eta_m,
                 'eta_tot':eta_tot }

    sig_n = Property
    def _get_sig_n( self ):
        return self.ls_values['sig_n']

    sig_m = Property
    def _get_sig_m( self ):
        return self.ls_values['sig_m']

    eta_n = Property
    def _get_eta_n( self ):
        return self.ls_values['eta_n']

    eta_m = Property
    def _get_eta_m( self ):
        return self.ls_values['eta_m']

    eta_tot = Property
    def _get_eta_tot( self ):
        return self.ls_values['eta_tot']

    assess_name = 'max_eta_tot'

#    assess_name = 'max_sig1_up'
#
#    max_sig1_up = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_sig1_up( self ):
#        return ndmax( self.sig1_up )

    # @todo: make it possible to select the assess value:
    #
#    assess_name = Enum( values = 'columns' )
#    def _assess_name_default( self ):
#        return self.columns[-1]

    max_eta_tot = Property( depends_on = '+input' )
    @cached_property
    def _get_max_eta_tot( self ):
        return ndmax( self.eta_tot )


    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    traits_view = View( VGroup( 
                            HGroup( Item( name = 'f_ctk', label = 'Tensile strength concrete [MPa]: f_ctk ' ),
                                    Item( name = 'f_m', label = 'Flexural tensile trength concrete [MPa]: f_m ' )
                                   ),
                            VGroup( 
                                Include( 'ls_group' ),

                                # @todo: currently LSArrayAdapter must be called both 
                                #        in SLS and ULS separately to configure columns 
                                #        arrangement individually
                                #
                                Item( 'ls_array', show_label = False,
                                      editor = TabularEditor( adapter = LSArrayAdapter() ) )
                                  ),
                              ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )

class ULS( LS ):
    '''Ultimate limit state
    '''

    #--------------------------------------------------------
    # ULS: material parameters (Inputs)
    #--------------------------------------------------------

    # gamma-factor 
    gamma = Float( 1.5, input = True )

    # long term reduction factor
    beta = Float( 1.0, input = True )

    # INDEX l: longitudinal direction of the textile (MAG-02-02-06a)
    # characteristic tensile strength of the tensile specimen [N/mm2]
    f_tk_l = Float( 537, input = True )

    # design value of the tensile strength of the tensile specimen [N/mm2]
    # containing a gamma-factor of 1.5 and d long term reduction factor of 0.7
    # f_td_l = 251

    f_td_l = Property( Float, depends_on = '+input' )
    def _get_f_td_l( self ):
        return self.beta * self.f_tk_l / self.gamma

    # cross sectional area of the reinforcement [mm2/m]
    a_t_l = Float( 71.65, input = True )

    # INDEX q: orthogonal direction of the textile (MAG-02-02-06a)
    # characteristic tensile strength of the tensile specimen [N/mm2]
    f_tk_q = Float( 511, input = True )

    # design value of the tensile strength of the tensile specimen [kN/m]
    # f_td_q = 238
    f_td_q = Property( Float, depends_on = '+input' )
    def _get_f_td_q( self ):
        return self.beta * self.f_tk_q / self.gamma

    # cross sectional area of the reinforcement [mm2/m]
    a_t_q = Float( 53.31, input = True )

    # tensile strength of the textile reinforcement [kN/m]
    f_Rtex_l = Property( Float, depends_on = '+input' )
    def _get_f_Rtex_l( self ):
        return self.a_t_l * self.f_td_l / 1000.

    # tensile strength of the textile reinforcement [kN/m]
    f_Rtex_q = Property( Float )
    def _get_f_Rtex_q( self, depends_on = '+input' ):
        return self.a_t_q * self.f_td_q / 1000.

    # tensile strength of the textile reinforcement [kN/m]
    # as simplification the lower value of 0- and 90-direction is taken
    #
    # sig_composite,exp = 103.4 kN/14cm/6cm = 12.3 MPa
    # f_tex,exp = 103.4 kN/14cm/12layers = 61.5 kN/m
    # f_tex,k = f_tex,exp * 0,81 (EN-DIN 1990) = 0.82*61.5 kN/m = 50.4 kN/m
    # f_tex,d = f_tex,k / 1,5 = 33.6 kN/m
    #
    f_Rtex_l = f_Rtex_q = 34.3
    print 'NOTE: f_Rtex_l = f_Rtex_q = set to %g kN/m !' % ( f_Rtex_l )

    k_fl = 7.95 / 6.87
    print 'NOTE: k_fl = set to %g [-] !' % ( k_fl )


    # ------------------------------------------------------------
    # ULS - derived params:
    # ------------------------------------------------------------

    # Parameters for the cracked state (GdT):
    # assumptions!

    # (resultierende statische Nutzhoehe) 
    #
    d = Property( Float )
    def _get_d( self ):
        return 0.75 * self.ls_table.D_elem

    # (Abstand Schwereachse zur resultierende Bewehrungslage) 
    # chose the same amount of reinforcement at the top as at the bottom 
    # i.e. zs = zs1 = zs2
    #
    zs = Property( Float )
    def _get_zs( self ):
        return self.d - self.ls_table.D_elem / 2.

    # (Innerer Hebelarm) 
    #
    z = Property( Float )
    def _get_z( self ):
        return 0.9 * self.d

    # ------------------------------------------------------------
    # ULS: outputs
    # ------------------------------------------------------------

    ls_columns = List( [ 'e', 'm_Eds', 'f_t', 'f_t_sig', 'beta_l', 'beta_q', 'f_Rtex', 'k_fl_NM', 'n_tex' ] )

    sr_columns = ['m', 'n', 'alpha', 'd', 'zs', 'z']

    alpha_varname = Property()
    def _get_alpha_varname( self ):
        return 'alpha_' + self.stress_res

    alpha = Property
    def _get_alpha( self ):
        return getattr( self.ls_table, self.alpha_varname )

    ls_values = Property( depends_on = '+input' )
    @cached_property
    def _get_ls_values( self ):
        '''get the outputs for ULS
        '''
        #-------------------------------------------------
        # VAR 1:use simplified reinforced concrete approach
        #-------------------------------------------------
        n = self.n
        m = self.m
        alpha = self.alpha

        zs = self.zs
        z = self.z
        f_Rtex_l = self.f_Rtex_l
        f_Rtex_q = self.f_Rtex_q

        # (Exzentrizitaet)
        e = abs( m / n )
        e[ n == 0 ] = 1E9 # if normal force is zero set e to very large value

        # moment at the height of the resulting reinforcement layer:
        m_Eds = abs( m ) - zs * n

        # tensile force in the reinforcement for bending and compression
        f_t = m_Eds / z + n

        # check if the two conditions are true:
        cond1 = n > 0
        cond2 = e < zs
        bool_arr = cond1 * cond2
        # in case of pure tension in the cross section:
        f_t[ bool_arr ] = n[ bool_arr ] * ( zs[ bool_arr ] + e[ bool_arr ] ) / ( zs[ bool_arr ] + zs[ bool_arr ] )

        #-------------------------------------------------
        # VAR 2:use principal stresses to calculate the resulting tensile force
        #-------------------------------------------------
        #
        princ_stress_eval = True
#        princ_stress_eval = False

        f_t_sig = self.sig1_lo * self.D_elem * 1000.
        if princ_stress_eval == True:
            print "NOTE: the principle tensile stresses are used to evaluate 'n_tex'"
            # resulting tensile force of the composite cross section[kN]
            # the entire (!) cross section is used!
            # as the maximum value of the tensile stresses at the top or the bottom 
            # i.e. sig1_max = min( 0, max( self.sig1_up, self.sig1_lo ) )


            #---------------------------------------------------------
            # initialize arrays t be filled by case distinction:
            #---------------------------------------------------------
            #
            f_t_sig = zeros_like ( self.sig1_up )
            alpha = zeros_like ( self.sig1_up )

            k_fl_NM = ones_like ( self.sig1_up )
            # absolute value of the bending stress created by mx or my
            sig_b = ( abs( self.sig1_lo ) + abs( self.sig1_up ) ) / 2

            #---------------------------------------------------------
            # conditions for case distinction  
            #---------------------------------------------------------

            cond_t_up = self.sig1_up > 0. # compression stress upper side
            cond_t_lo = self.sig1_lo > 0. # compression stress lower side

            cond_u_gt_l = abs( self.sig1_up ) > abs( self.sig1_lo ) # absolute value of upper stress greater then lower stress
            cond_l_gt_u = abs( self.sig1_lo ) > abs( self.sig1_up ) # absolute value of lower stress greater then lower stress
            cond_u_eq_l = abs( self.sig1_up ) == abs( self.sig1_lo ) # absolute values of upper stress equals lower stress

            cond_b = self.sig1_up * self.sig1_lo < 0 # bending case
            cond_u_gt_0 = self.sig1_up > 0 # value of upper stress greater then 0.
            cond_l_gt_0 = self.sig1_lo > 0 # value of lower stress greater then 0.

            cond_c_up = self.sig1_up < 0. # compression stress upper side
            cond_c_lo = self.sig1_lo < 0. # compression stress lower side


            #---------------------------------------------------------
            # tension case:
            #---------------------------------------------------------

            # sig_up > sig_lo
            bool_arr = cond_t_up * cond_t_lo * cond_u_gt_l
            alpha[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            k_fl_NM[ bool_arr ] = 1.0
            f_t_sig[ bool_arr ] = self.sig1_up[ bool_arr ] * self.D_elem[ bool_arr ] * 1000.

            # sig_lo > sig_up
            bool_arr = cond_t_up * cond_t_lo * cond_l_gt_u
            alpha[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            k_fl_NM[ bool_arr ] = 1.0
            f_t_sig[ bool_arr ] = self.sig1_lo[ bool_arr ] * self.D_elem[ bool_arr ] * 1000.

            # sig_lo = sig_up
            bool_arr = cond_t_up * cond_t_lo * cond_u_eq_l
            alpha[ bool_arr ] = ( self.alpha_sig_lo[ bool_arr ] + self.alpha_sig_up[ bool_arr ] ) / 2.
            k_fl_NM[ bool_arr ] = 1.0
            f_t_sig[ bool_arr ] = self.sig1_lo[ bool_arr ] * self.D_elem[ bool_arr ] * 1000.


            #---------------------------------------------------------
            # bending case (= different signs)
            #---------------------------------------------------------

            # bending with tension at the lower side 
            # AND sig_N > 0 (--> bending and tension; sig1_lo results from bending and tension)
            bool_arr = cond_b * cond_l_gt_0 * cond_l_gt_u
            alpha[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            k_fl_NM[ bool_arr ] = 1.0 + ( self.k_fl - 1.0 ) * \
                                 ( 1.0 - ( sig_b[ bool_arr ] - abs( self.sig1_up[ bool_arr] ) ) / self.sig1_lo[ bool_arr ] )
            f_t_sig[ bool_arr ] = self.sig1_lo[ bool_arr ] * self.D_elem[ bool_arr ] * 1000. / k_fl_NM[ bool_arr ]

            # bending with tension at the lower side 
            # AND sig_N < 0 (--> bending and compression; sig1_lo results only from bending)
            bool_arr = cond_b * cond_l_gt_0 * cond_u_gt_l
            alpha[ bool_arr ] = self.alpha_sig_lo[ bool_arr ]
            k_fl_NM[ bool_arr ] = self.k_fl
            f_t_sig[ bool_arr ] = self.sig1_lo[ bool_arr ] * self.D_elem[ bool_arr ] * 1000. / k_fl_NM[ bool_arr ]

            # bending with tension at the upper side 
            # AND sig_N > 0 (--> bending and tension; sig1_up results from bending and tension)
            bool_arr = cond_b * cond_u_gt_0 * cond_u_gt_l
            alpha[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            k_fl_NM[ bool_arr ] = 1.0 + ( self.k_fl - 1.0 ) * \
                                 ( 1.0 - ( sig_b[ bool_arr ] - abs( self.sig1_lo[ bool_arr] ) ) / self.sig1_up[ bool_arr ] )
            f_t_sig[ bool_arr ] = self.sig1_up[ bool_arr ] * self.D_elem[ bool_arr ] * 1000. / k_fl_NM[ bool_arr ]

            # bending with tension at the upper side 
            # AND sig_N < 0 (--> bending and compression; sig1_up results only from bending)
            bool_arr = cond_b * cond_u_gt_0 * cond_l_gt_u
            alpha[ bool_arr ] = self.alpha_sig_up[ bool_arr ]
            k_fl_NM[ bool_arr ] = self.k_fl
            f_t_sig[ bool_arr ] = self.sig1_up[ bool_arr ] * self.D_elem[ bool_arr ] * 1000. / k_fl_NM[ bool_arr ]

            #---------------------------------------------------------
            # compression case:
            #---------------------------------------------------------

            bool_arr = cond_c_up * cond_c_lo
            # deflection angle is of minor interest use this simplification
            alpha[ bool_arr ] = ( self.alpha_sig_up[ bool_arr ] + self.alpha_sig_up[ bool_arr ] ) / 2.
            f_t_sig[ bool_arr ] = 0.
            k_fl_NM[ bool_arr ] = 1.0


        #------------------------------------------------------------
        # get angel of deflection of the textile reinforcement
        #------------------------------------------------------------

        # angel of deflection of the textile reinforcement for dimensioning in x-direction
        # distinguished between longitudinal (l) and transversal (q) direction

#        # ASSUMPTION: worst case angle used
#        # as first step use the worst case reduction due to deflection possible (at 55 degrees)
#        beta_l = 55. * pi / 180. * ones_like( alpha )
#        beta_q = ( 90. - 55. ) * pi / 180. * ones_like( alpha )

        # @todo: as second step use the value for an alternating layup (i.e. deflection angle)
        # @todo: get the correct formula for the demonstrator arrangement 
        # i.e. the RFEM coordinate system orientation     
#        print "NOTE: deflection angle is used to evaluate 'n_tex'"
        beta_l = 90 - abs( alpha ) # [degree]   
        beta_q = abs( alpha )      # [degree]   

        # resulting strength of the bi-directional textile considering the 
        # deflection of the reinforcement in the loading direction:
        f_Rtex = f_Rtex_l * cos( beta_l * pi / 180. ) * ( 1 - beta_l / 90. ) + \
                 f_Rtex_q * cos( beta_q * pi / 180. ) * ( 1 - beta_q / 90. )

        # f_Rtex= 11.65 kN/m corresponds to sig_Rtex = 585 MPa
        #
#        f_Rtex = 11.65 * ones_like( alpha )

        # f_Rtex= 10.9 kN/m corresponds to sig_Rtex = 678 MPa
        # with 1/3*( 679 + 690 + 667 ) = 678 MPa
        #
#        f_Rtex = 10.9 * ones_like( alpha )

        # f_Rtex = 10.00 kN/m corresponds to sig_Rtex = 610 MPa
        #
#        f_Rtex = 10.00 * ones_like( alpha )

#        print 'NOTE: f_Rtex set to %g kN/m !' % ( f_Rtex[0] )

        if princ_stress_eval == False:
            # necessary number of reinforcement layers
            n_tex = f_t / f_Rtex

            return { 'e':e, 'm_Eds':m_Eds, 'f_t':f_t, 'f_t_sig' : f_t_sig,
                     'beta_l':beta_l, 'beta_q':beta_q, 'f_Rtex':f_Rtex,
                     'n_tex':n_tex}

        elif princ_stress_eval == True:
            # NOTE: as the entire (!) cross section is used to calculate 'f_tex_sig'
            # the number of layers is evaluated also directly for the entire (!)
            # cross section! 
            #
            n_tex = f_t_sig / f_Rtex

            return { 'e':e, 'm_Eds':m_Eds, 'f_t':f_t, 'f_t_sig' : f_t_sig,
                     'beta_l':beta_l, 'beta_q':beta_q, 'f_Rtex':f_Rtex,
                     'n_tex':n_tex, 'k_fl_NM' : k_fl_NM }

    e = Property
    def _get_e( self ):
        return self.ls_values['e']

    m_Eds = Property
    def _get_m_Eds( self ):
        return self.ls_values['m_Eds']

    f_t = Property
    def _get_f_t( self ):
        return self.ls_values['f_t']

    f_t_sig = Property
    def _get_f_t_sig( self ):
        return self.ls_values['f_t_sig']

    beta_l = Property
    def _get_beta_l( self ):
        return self.ls_values['beta_l']

    beta_q = Property
    def _get_beta_q( self ):
        return self.ls_values['beta_q']

    f_Rtex = Property
    def _get_f_Rtex( self ):
        return self.ls_values['f_Rtex']

    k_fl_NM = Property
    def _get_k_fl_NM( self ):
        return self.ls_values['k_fl_NM']

    n_tex = Property
    def _get_n_tex( self ):
        return self.ls_values['n_tex']

    assess_name = 'max_n_tex'
    max_n_tex = Property( depends_on = '+input' )
    @cached_property
    def _get_max_n_tex( self ):
        return ndmax( self.n_tex )

#    assess_name = 'max_sig1_up'
#    max_sig1_up = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_sig1_up( self ):
#        return ndmax( self.sig1_up )

#    assess_name = 'max_sig1_lo'
#    max_sig1_lo = Property( depends_on = '+input' )
#    @cached_property
#    def _get_max_sig1_lo( self ):
#        return ndmax( self.sig1_lo )

    #-------------------------------
    # ls view
    #-------------------------------

    # @todo: the dynamic selection of the columns to be displayed 
    # does not work in connection with the LSArrayAdapter 
    traits_view = View( 
                       VGroup( 
                        HGroup( 
                            VGroup( 
                                Item( name = 'gamma', label = 'safety factor material [-]:  gamma ' ),
                                Item( name = 'beta', label = 'reduction long term durability [-]:  beta ' ),
                                label = 'safety factors'
                                  ),
                            VGroup( 
                                Item( name = 'f_tk_l', label = 'characteristic strength textil [MPa]:  f_tk_l ', format_str = "%.1f" ),
                                Item( name = 'f_td_l', label = 'design strength textil [MPa]:  f_td_l ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'a_t_l', label = 'cross sectional area textil [mm^2]:  a_t_l ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'f_Rtex_l', label = 'design strength textil [kN/m]:  f_Rtex_l ', style = 'readonly', format_str = "%.0f" ),
                                label = 'material properties (longitudinal)'
                                  ),
                            VGroup( 
                                Item( name = 'f_tk_q', label = 'characteristic strength textil [MPa]:  f_tk_q ', format_str = "%.1f" ),
                                Item( name = 'f_td_q', label = 'design strength textil [MPa]:  f_td_q ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'a_t_q', label = 'cross sectional area textil [mm^2]:  a_t_q ', style = 'readonly', format_str = "%.1f" ),
                                Item( name = 'f_Rtex_q', label = 'design strength textil [kN/m]: f_Rtex_q ', style = 'readonly', format_str = "%.0f" ),
                                label = 'material Properties (transversal)'
                                  ),
                             ),

                        VGroup( 
                            Include( 'ls_group' ),
                            Item( 'ls_array', show_label = False,
                                  editor = TabularEditor( adapter = LSArrayAdapter() ) )
                              ),
                            ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )

LSLIST = [ SLS, ULS ]

class LSTable( HasTraits ):
    '''Assessment tool
    '''

    is_id = Int( 0 )
    # geo data: coordinates and element thickness
    # 
    geo_data = Dict

    elem_no = Property( Array )
    def _get_elem_no( self ):
        return self.geo_data['elem_no']

    X = Property( Array )
    def _get_X( self ):
        return self.geo_data['X']

    Y = Property( Array )
    def _get_Y( self ):
        return self.geo_data['Y']

    Z = Property( Array )
    def _get_Z( self ):
        return self.geo_data['Z']

    D_elem = Property( Array )
    def _get_D_elem( self ):
        '''element thickness (units changed form [mm] to [m])'''
        return self.geo_data['thickness'] / 1000.

    # state data: stress resultants 
    # 
    state_data = Dict

    mx = Property( Array )
    def _get_mx( self ):
        return self.state_data['mx']

    my = Property( Array )
    def _get_my( self ):
        return self.state_data['my']

    mxy = Property( Array )
    def _get_mxy( self ):
        return self.state_data['mxy']

    nx = Property( Array )
    def _get_nx( self ):
        return self.state_data['nx']

    ny = Property( Array )
    def _get_ny( self ):
        return self.state_data['ny']

    nxy = Property( Array )
    def _get_nxy( self ):
        return self.state_data['nxy']


    # ------------------------------------------------------------
    # Index M: calculate principle moments with corresponding normal forces
    # ------------------------------------------------------------

    princ_values_M = Property( Dict, depends_on = 'data_file_stress_resultants' )
    @cached_property
    def _get_princ_values_M( self ):
        '''principle value of the moments forces:
        and principle angle of the moments forces:
        mx_M, my_M, nx_M, ny_M: transform the values in the principle direction
        '''
        # stress_resultants in global coordinates
        #
        mx = self.mx
        my = self.my
        mxy = self.mxy
        nx = self.nx
        ny = self.ny
        nxy = self.nxy

        # principal values
        #
        m1 = 0.5 * ( mx + my ) + 0.5 * sqrt( ( mx - my ) ** 2 + 4 * mxy ** 2 )
        m2 = 0.5 * ( mx + my ) - 0.5 * sqrt( ( mx - my ) ** 2 + 4 * mxy ** 2 )
        alpha_M = pi / 2. * ones_like( m1 )
        bool = m2 != mx
        alpha_M[ bool ] = arctan( mxy[ bool ] / ( m2[ bool ] - mx[ bool ] ) )
        alpha_M_deg = alpha_M * 180. / pi

        # transform to principal directions
        #
        mx_M = 0.5 * ( my + mx ) - 0.5 * ( my - mx ) * cos( 2 * alpha_M ) - mxy * sin( 2 * alpha_M )
        my_M = 0.5 * ( my + mx ) + 0.5 * ( my - mx ) * cos( 2 * alpha_M ) + mxy * sin( 2 * alpha_M )
        nx_M = 0.5 * ( ny + nx ) - 0.5 * ( ny - nx ) * cos( 2 * alpha_M ) - nxy * sin( 2 * alpha_M )
        ny_M = 0.5 * ( ny + nx ) + 0.5 * ( ny - nx ) * cos( 2 * alpha_M ) + nxy * sin( 2 * alpha_M )
        return { 'm1':m1, 'm2':m2, 'alpha_M':alpha_M_deg,
                 'mx_M':mx_M, 'my_M':my_M,
                 'nx_M':nx_M, 'ny_M':ny_M }

    m1 = Property( Float )
    def _get_m1( self ):
        return self.princ_values_M['m1']

    m2 = Property( Float )
    def _get_m2( self ):
        return self.princ_values_M['m2']

    alpha_M = Property( Float )
    def _get_alpha_M( self ):
        return self.princ_values_M['alpha_M']

    mx_M = Property( Float )
    def _get_mx_M( self ):
        return self.princ_values_M['mx_M']

    my_M = Property( Float )
    def _get_my_M( self ):
        return self.princ_values_M['my_M']

    nx_M = Property( Float )
    def _get_nx_M( self ):
        return self.princ_values_M['nx_M']

    ny_M = Property( Float )
    def _get_ny_M( self ):
        return self.princ_values_M['ny_M']

    # ------------------------------------------------------------
    # Index N: principle normal forces with corresponding moments
    # ------------------------------------------------------------

    princ_values_N = Property( Dict, depends_on = 'data_file_stress_resultants' )
    @cached_property
    def _get_princ_values_N( self ):
        '''principle value of the normal forces:
        and principle angle of the normal forces:
        mx_N, my_N, nx_N, ny_N: transform the values in the principle normal direction
        '''
        # stress_resultants in global coordinates
        #
        mx = self.mx
        my = self.my
        mxy = self.mxy
        nx = self.nx
        ny = self.ny
        nxy = self.nxy

        # principal values
        #
        n1 = 0.5 * ( nx + ny ) + 0.5 * sqrt( ( nx - ny ) ** 2 + 4 * nxy ** 2 )
        n2 = 0.5 * ( nx + ny ) - 0.5 * sqrt( ( nx - ny ) ** 2 + 4 * nxy ** 2 )
        alpha_N = pi / 2. * ones_like( n1 )
        bool = n2 != nx
        alpha_N[ bool ] = arctan( nxy[ bool ] / ( n2[ bool ] - nx[ bool ] ) )
        alpha_N_deg = alpha_N * 180. / pi

        # transform to principal directions
        mx_N = 0.5 * ( my + mx ) - 0.5 * ( my - mx ) * cos( 2 * alpha_N ) - mxy * sin( 2 * alpha_N )
        my_N = 0.5 * ( my + mx ) + 0.5 * ( my - mx ) * cos( 2 * alpha_N ) + mxy * sin( 2 * alpha_N )
        nx_N = 0.5 * ( ny + nx ) - 0.5 * ( ny - nx ) * cos( 2 * alpha_N ) - nxy * sin( 2 * alpha_N )
        ny_N = 0.5 * ( ny + nx ) + 0.5 * ( ny - nx ) * cos( 2 * alpha_N ) + nxy * sin( 2 * alpha_N )

        return{'n1' : n1, 'n2' : n2, 'alpha_N' : alpha_N_deg,
               'mx_N' : mx_N, 'my_N' : my_N,
               'nx_N' : nx_N, 'ny_N' : ny_N }

    n1 = Property( Float )
    def _get_n1( self ):
        return self.princ_values_N['n1']

    n2 = Property( Float )
    def _get_n2( self ):
        return self.princ_values_N['n2']

    alpha_N = Property( Float )
    def _get_alpha_N( self ):
        return self.princ_values_N['alpha_N']

    mx_N = Property( Float )
    def _get_mx_N( self ):
        return self.princ_values_N['mx_N']

    my_N = Property( Float )
    def _get_my_N( self ):
        return self.princ_values_N['my_N']

    nx_N = Property( Float )
    def _get_nx_N( self ):
        return self.princ_values_N['nx_N']

    ny_N = Property( Float )
    def _get_ny_N( self ):
        return self.princ_values_N['ny_N']

    # ------------------------------------------------------------
    # Index sig: calculate principle stresses 
    # ------------------------------------------------------------

    princ_values_sig = Property( Dict, depends_on = 'data_file_stress_resultants' )
    @cached_property
    def _get_princ_values_sig( self ):
        '''principle value of the stresses for the lower ('lo') and upper ('up') face:
        '''
        # stress_resultants in global coordinates
        #
        mx = self.mx
        my = self.my
        mxy = self.mxy
        nx = self.nx
        ny = self.ny
        nxy = self.nxy

        # geometrical properties:
        #
        A = self.D_elem * 1000.
        W = self.D_elem ** 2 / 6 * 1000.

        # compare the formulae with the RFEM-manual p.290

        # normal stresses upper face:
        #
        sigx_up = nx / A - mx / W
        sigy_up = ny / A - my / W
        sigxy_up = nxy / A - mxy / W

        # principal stresses upper face:
        #
        sig1_up = 0.5 * ( sigx_up + sigy_up ) + 0.5 * sqrt( ( sigx_up - sigy_up ) ** 2 + 4 * sigxy_up ** 2 )
        sig2_up = 0.5 * ( sigx_up + sigy_up ) - 0.5 * sqrt( ( sigx_up - sigy_up ) ** 2 + 4 * sigxy_up ** 2 )
        alpha_sig_up = pi / 2. * ones_like( sig1_up )
        bool = sig2_up != sigx_up
        alpha_sig_up[ bool ] = arctan( sigxy_up[ bool ] / ( sigy_up[ bool ] - sigx_up[ bool ] ) )
        alpha_sig_up_deg = alpha_sig_up * 180. / pi

        # normal stresses lower face:
        #
        sigx_lo = nx / A + mx / W
        sigy_lo = ny / A + my / W
        sigxy_lo = nxy / A + mxy / W

        # principal stresses lower face:
        #
        sig1_lo = 0.5 * ( sigx_lo + sigy_lo ) + 0.5 * sqrt( ( sigx_lo - sigy_lo ) ** 2 + 4 * sigxy_lo ** 2 )
        sig2_lo = 0.5 * ( sigx_lo + sigy_lo ) - 0.5 * sqrt( ( sigx_lo - sigy_lo ) ** 2 + 4 * sigxy_lo ** 2 )
        alpha_sig_lo = pi / 2. * ones_like( sig1_lo )
        bool = sig2_lo != sigx_lo
        alpha_sig_lo[ bool ] = arctan( sigxy_lo[ bool ] / ( sig2_lo[ bool ] - sigx_lo[ bool ] ) )
        alpha_sig_lo_deg = alpha_sig_lo * 180. / pi

        return {
                 'sigx_up' : sigx_up, 'sigy_up' : sigy_up, 'sigxy_up' : sigxy_up,
                 'sig1_up' : sig1_up, 'sig2_up' : sig2_up, 'alpha_sig_up' : alpha_sig_up_deg,

                 'sigx_lo' : sigx_lo, 'sigy_lo' : sigy_lo, 'sigxy_lo' : sigxy_lo,
                 'sig1_lo' : sig1_lo, 'sig2_lo' : sig2_lo, 'alpha_sig_lo' : alpha_sig_lo_deg,
               }

    # stresses upper face:
    #
    sigx_up = Property( Float )
    def _get_sigx_up( self ):
        return self.princ_values_sig['sigx_up']

    sigy_up = Property( Float )
    def _get_sigy_up( self ):
        return self.princ_values_sig['sigy_up']

    sigxy_up = Property( Float )
    def _get_sigxy_up( self ):
        return self.princ_values_sig['sigxy_up']

    sig1_up = Property( Float )
    def _get_sig1_up( self ):
        return self.princ_values_sig['sig1_up']

    sig2_up = Property( Float )
    def _get_sig2_up( self ):
        return self.princ_values_sig['sig2_up']

    alpha_sig_up = Property( Float )
    def _get_alpha_sig_up( self ):
        return self.princ_values_sig['alpha_sig_up']

    # stresses lower face:
    #
    sigx_lo = Property( Float )
    def _get_sigx_lo( self ):
        return self.princ_values_sig['sigx_lo']

    sigy_lo = Property( Float )
    def _get_sigy_lo( self ):
        return self.princ_values_sig['sigy_lo']

    sigxy_lo = Property( Float )
    def _get_sigxy_lo( self ):
        return self.princ_values_sig['sigxy_lo']

    sig1_lo = Property( Float )
    def _get_sig1_lo( self ):
        return self.princ_values_sig['sig1_lo']

    sig2_lo = Property( Float )
    def _get_sig2_lo( self ):
        return self.princ_values_sig['sig2_lo']

    alpha_sig_lo = Property( Float )
    def _get_alpha_sig_lo( self ):
        return self.princ_values_sig['alpha_sig_lo']

    #------------------------------------------
    # combinations of limit states, stress resultants and directions
    #------------------------------------------

    ls = Trait( 'ULS',
                {'ULS' : ULS,
                 'SLS' : SLS } )

    sr_tree = Dict
    def _sr_tree_default( self ):
        '''collect the ls instances in a dictionary tree: 
        sr_tree[ sr ][ dir ] = LSTable-object
        e.g. sr_tree[ 'M' ][ 'x' ] = ULS( dir = x, sr = M )
        '''
        dir_list = DIRLIST
        sr_list = SRLIST
        ls_class = self.ls_
        sr_dict = {}
        for sr in sr_list:
            dir_dict = {}
            for dir in dir_list:
                dir_dict[ dir ] = ls_class( ls_table = self, dir = dir, sr = sr )
            sr_dict[ sr ] = dir_dict
        return sr_dict

    ls_list = Property
    def _get_ls_list( self ):
        '''collect the ls instances in a list: 
        e.g. ls_list = [ sr_tree['M']['x'], sr_tree['M']['y'],
                         sr_tree['N']['x'], sr_tree['N']['y'] ]
        (NOTE: list order is undefined because it is constructed from a dictionary)
        '''
        ls_list = []
        for sr in self.sr_tree.values():
            for dir in sr.values():
                ls_list.append( dir )
        return ls_list

    #------------------------------------------
    # get max value and case for all cases:
    #------------------------------------------

    max_value_and_case = Dict
    def _max_value_and_case_default( self ):

        dir_list = DIRLIST
        sr_list = SRLIST
        sr_tree = self.sr_tree

        # Derive the column list from the first case.
        # NOTE: (all cases are structured identically, otherwise no 
        # overall comparison would be possible.
        #
        columns = sr_tree[ 'M' ][ 'x' ].columns

        # run through allvariables defined in the 'column'
        # attribute of the limit state class ('LS')
        #
        var_dir = {}
        for var in columns:
            # reset auxilary vairables
            #
            max_value_all = None
            max_case_all = None
            # run through all cases (e.g. 'Nx', 'Ny', 'Mx', 'My' )
            #
            for sr in sr_list:
                for dir in dir_list:
                    # set the 'max_in_column' attribute of the LS-class:
                    # and get the max value of the currently selected limit state:
                    # 
                    col = getattr( self.sr_tree[ sr ][ dir ], var ) [:, 0]
                    max_value_current = max( col )

                    # compare with the maximum reached so far in the investigation 
                    # of all cases:
                    #
                    if max_value_current >= max_value_all:
                        max_value_all = max_value_current
                        max_case_all = sr + dir
                        var_dir[ var ] = {'max_value' : max_value_all,
                                          'max_case' : max_case_all }
        return var_dir


    #------------------------------------------
    # get arrays for the TabularEditor:
    #------------------------------------------

    Mx = Property( Instance( LS ) )
    def _get_Mx( self ):
        return self.sr_tree['M']['x']

    My = Property( Instance( LS ) )
    def _get_My( self ):
        return self.sr_tree['M']['y']

    Nx = Property( Instance( LS ) )
    def _get_Nx( self ):
        return self.sr_tree['N']['x']

    Ny = Property( Instance( LS ) )
    def _get_Ny( self ):
        return self.sr_tree['N']['y']

    assess_value = Property
    def _get_assess_value( self ):
        return max( [ getattr( ls, ls.assess_name ) for ls in self.ls_list ] )

    # ------------------------------------------------------------
    # View 
    # ------------------------------------------------------------

    traits_view = View( Tabbed( 
                            Item( 'Nx@' , label = "NX", show_label = False ),
                            Item( 'Ny@' , label = "NY", show_label = False ),
                            Item( 'Mx@' , label = "MX", show_label = False ),
                            Item( 'My@' , label = "MY", show_label = False ),
                            scrollable = False,
                         ),
                      resizable = True,
                      scrollable = True,
                      height = 1000,
                      width = 1100
                      )


