'''
Created on Jun 29, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, DelegatesTo, Callable

from etsproxy.util.home_directory import \
    get_home_directory

from etsproxy.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from etsproxy.mayavi import \
    mlab

from etsproxy.traits.ui.table_column import \
    ObjectColumn

from etsproxy.traits.ui.menu import \
    OKButton, CancelButton

from etsproxy.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import \
    array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
    vstack, savetxt, hstack, argsort, fromstring, zeros_like, shape, \
    copy, c_, newaxis, argmax, where, sqrt, frompyfunc, sum, \
    ones, transpose, shape, append, argmin, fabs, identity, unique, vdot, \
    max as ndmax, min as ndmin

from ls_table_hf import \
    LSTable, ULS, SLS

from math import pi
from string import split
import os

from promod.simdb import \
    SimDB

import pickle
import string
import csv
from os.path import join


# Access to the top level directory of the database
#
simdb = SimDB()




class LC( HasTraits ):
    '''Loading case class
    '''
    # name of the file containing the hinge forces 
    # or the displacements
    #
    file_name = Str( input = True )

    # data filter (used to hide unwanted values, e.g. high sigularities etc.) 
    #
    data_filter = Callable( input = True )
    def _data_filter_default( self ):
        return lambda x: x    #  - do nothing by default

    # name of the loading case
    #
    name = Str( input = True )

    # category of the loading case
    #
    category = Enum( 'dead-load', 'additional dead-load', 'imposed-load', input = True )

    # list of keys specifying the names of the loading cases 
    # that can not exist at the same time, i.e. which are exclusive to each other
    # 
    exclusive_to = List( Str, input = True )
    def _exclusive_to_default( self ):
        return []

    # combination factors (need to be defined in case of imposed loads)
    # 
    psi_0 = Float( input = True )
    psi_1 = Float( input = True )
    psi_2 = Float( input = True )

    # security factors ULS
    #
    gamma_fav = Float( input = True )
    def _gamma_fav_default( self ):
        if self.category == 'dead-load':
            return 1.00
        if self.category == 'additional dead-load':
            return 0.00
        if self.category == 'imposed-load':
            return 0.00

    gamma_unf = Float( input = True )
    def _gamma_unf_default( self ):
        if self.category == 'dead-load':
            return 1.35
        if self.category == 'additional dead-load':
            return 1.35
        if self.category == 'imposed-load':
            return 1.50

    # security factors SLS:
    # (used to distinguish combinations where imposed-loads
    # or additional-dead-loads are favorable or unfavorable.) 
    #
    gamma_fav_SLS = Float( input = True )
    def _gamma_fav_SLS_default( self ):
        if self.category == 'dead-load':
            return 1.00
        elif self.category == 'additional dead-load' or \
             self.category == 'imposed-load':
            return 0.00

    gamma_unf_SLS = Float( input = True )
    def _gamma_unf_SLS_default( self ):
        return 1.00


    def _read_data( self, file_name ):
        '''read state data and geo data from csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
        print '*** read state data from file: %s ***' % ( file_name )

        # get the column headings defined in the first row 
        # of the csv hinge force input file
        # column_headings = array(["elem_no", "N_ip", "V_ip", "V_op"])
        #
        file = open( file_name, 'r' )
        lines = file.readlines()
        column_headings = lines[0].split( ';' )
        # remove '\n' from last string element in list
        #
        column_headings[-1] = column_headings[-1][:-1]
        column_headings_arr = array( column_headings )

        # data_hf:
        #
#        elem_no_idx = where( 'Nr.' == column_headings_arr )[0]


        # state_data:
        #
        N_ip_idx = where( 'N_ip[MN]' == column_headings_arr )[0]
        V_ip_idx = where( 'V_ip[MN]' == column_headings_arr )[0]
        V_op_idx = where( 'V_op[MN]' == column_headings_arr )[0]

        X_hf_idx = where( 'x[m]' == column_headings_arr )[0]
        Y_hf_idx = where( 'y[m]' == column_headings_arr )[0]
        Z_hf_idx = where( 'z[m]' == column_headings_arr )[0]

        U_x_idx = where( 'u-x[m]' == column_headings_arr )[0]
        U_y_idx = where( 'u-y[m]' == column_headings_arr )[0]
        U_z_idx = where( 'u-z[m]' == column_headings_arr )[0]
        dU_y_idx = where( 'du-y[m]' == column_headings_arr )[0]
        dU_z_idx = where( 'du-z[m]' == column_headings_arr )[0]

        X_u_idx = where( 'X[m]' == column_headings_arr )[0]
        Y_u_idx = where( 'Y[m]' == column_headings_arr )[0]
        Z_u_idx = where( 'Z[m]' == column_headings_arr )[0]

        file.close()

        input_arr = loadtxt( file_name , delimiter = ';', skiprows = 1 )

        # hinge forces [kN]
        #

        X_hf = input_arr[:, X_hf_idx]
        Y_hf = input_arr[:, Y_hf_idx]
        Z_hf = input_arr[:, Z_hf_idx]


        N_ip = 1000 * input_arr[:, N_ip_idx]
        V_ip = 1000 * input_arr[:, V_ip_idx]
        V_op = 1000 * input_arr[:, V_op_idx]


        # displacement data
        #       

        U_x = input_arr[:, U_x_idx]
        U_y = input_arr[:, U_y_idx]
        U_z = input_arr[:, U_z_idx]
        dU_y = input_arr[:, dU_y_idx]
        dU_z = input_arr[:, dU_z_idx]


        X_u = input_arr[:, X_u_idx]
        Y_u = input_arr[:, Y_u_idx]
        Z_u = input_arr[:, Z_u_idx]

        return {'X_hf' : X_hf, 'Y_hf' : Y_hf, 'Z_hf' : Z_hf,
                'N_ip' : N_ip, 'V_ip' : V_ip, 'V_op' : V_op,
                'X_u' : X_u, 'Y_u' : Y_u, 'Z_u' : Z_u,
                'U_x' : U_x, 'U_y' : U_y, 'U_z' : U_z, 'dU_y':dU_y, 'dU_z':dU_z
                }

    # original data (before filtering)
    #
    data_orig = Property( Dict, depends_on = 'file_name' )
    @cached_property
    def _get_data_orig( self ):
        return self._read_data( self.file_name )

    # data (after filtering)
    #
    data_dict = Property( Dict, depends_on = 'file_name, data_filter' )
    @cached_property
    def _get_data_dict( self ):
        d = {}
        for k, arr in self.data_orig.items():
            d[k] = self.data_filter( arr )
        return d


    #---------------------------------------------
    # NEEDS TO BE CHANGED FOR DISPLACEMENT OR FORCES DEPENDING ON OPTION
    #---------------------------------------------

    # use this line to evaluate the displacement files using LCCTable as combination tool
    # (for the files with spring elements the displacement in the .csv-files is only given for U_x, U_y and U_z)

    # if do == hf (hinge force evaluation)
    #
    sr_columns = List( ['N_ip', 'V_ip', 'V_op'] )
#    sr_columns = List( ['U_x', 'U_y', 'U_z', 'dU_y', 'dU_z'] )
    geo_columns = List( ['X_hf', 'Y_hf', 'Z_hf'] )

    # if do == u (edge displacement evaluation)
    #
#    sr_columns = List( ['U_x', 'U_y', 'U_z'] )
#    geo_columns = List( ['X_u', 'Y_u', 'Z_u'] )

    sr_arr = Property( Array )
    def _get_sr_arr( self ):
        '''return the stress resultants of the loading case 
        as stack of all sr-column arrays.
        '''
        data_dict = self.data_dict
        return hstack( [ data_dict[ sr_key ] for sr_key in self.sr_columns ] )

    geo_data_dict = Property( Array )
    def _get_geo_data_dict( self ):
        '''Dict of coords as sub-Dict of read in data dict
        '''
        data_dict = self.data_dict
        geo_data_dict = {}
        for geo_key in self.geo_columns:
            geo_data_dict[ geo_key ] = data_dict[ geo_key ]
        return geo_data_dict


    def _read_data_u( self, file_name_u ):
        '''read state data and geo date from csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter for displacement.
        '''
        print '*** read state data from file: %s ***' % ( file_name_u )

        # get the column headings defined in the first row 
        # of the csv hinge force input file
        # column_headings = array(["elem_no", "N_ip", "V_ip", "V_op"])
        #
        file = open( file_name_u, 'r' )
        lines = file.readlines()
        column_headings = lines[0].split( ';' )
        # remove '\n' from last string element in list
        #
        column_headings[-1] = column_headings[-1][:-1]
        column_headings_arr = array( column_headings )

        # geo_data:
        #
        X_idx = where( 'X[m]' == column_headings_arr )[0]
        Y_idx = where( 'Y[m]' == column_headings_arr )[0]
        Z_idx = where( 'Z[m]' == column_headings_arr )[0]

        # state_data:
        #
        U_x = where( 'u-x[m]' == column_headings_arr )[0]
        U_y = where( 'u-y[m]' == column_headings_arr )[0]
        U_z = where( 'u-z[m]' == column_headings_arr )[0]
        dU_y_idx = where( 'du-y[m]' == column_headings_arr )[0]
        dU_z_idx = where( 'du-z[m]' == column_headings_arr )[0]

        file.close()
        input_arr = loadtxt( file_name_u , delimiter = ';', skiprows = 1 )

        # global coords [m]
        #
        self.X_u = input_arr[:, X_idx]
        self.Y_u = input_arr[:, Y_idx]
        self.Z_u = input_arr[:, Z_idx]

        # hinge forces [kN]
        #
        self.U_x = input_arr[:, dU_y_idx]
        self.U_y = input_arr[:, dU_y_idx]
        self.U_z = input_arr[:, dU_y_idx]

        self.dU_y_idx = input_arr[:, dU_y_idx]
        self.dU_z_idx = input_arr[:, dU_z_idx]




class LCC( HasTraits ):

    lcc_id = Int

    #lcc_table = WeakRef()

    ls_table = Instance( LSTable )

    assess_value = Property()
    def _get_assess_value( self ):
        return self.ls_table.assess_value

    traits_view = View( Item( 'ls_table@', show_label = False ),
                        resizable = True,
                        scrollable = True
                        )

# The definition of the demo TableEditor:
lcc_list_editor = TableEditor( 
    columns_name = 'lcc_table_columns',
    editable = False,
    selection_mode = 'row',
    selected = 'object.lcc',
    show_toolbar = True,
    auto_add = False,
    configurable = True,
    sortable = True,
    reorderable = False,
    sort_model = False,
    auto_size = False,
    )

class LCCTable( HasTraits ):
    '''Loading Case Manager.
    Generates and sorts the loading case combinations
    of all specified loading cases.
    '''

    # define ls
    #
    ls = Trait( 'ULS',
                {'ULS' : ULS,
                 'SLS' : SLS } )

    # lcc-instance for the view 
    #
    lcc = Instance( LCC )

    #-------------------------------
    # Define loading cases:
    #-------------------------------

    # path to the directory containing the state data files
    #
    data_dir = Directory

    # list of load cases
    #
    lc_list_ = List( Instance( LC ) )

    lc_list = Property( List, depends_on = '+filter' )
    def _set_lc_list( self, value ):
        self.lc_list_ = value
    def _get_lc_list( self ):
#        for lc in self.lc_list_:
#            if lc.data_filter != self.data_filter:
#                lc.data_filter = self.data_filter
        return self.lc_list_

    lcc_table_columns = Property( depends_on = 'lc_list_, +filter' )
    def _get_lcc_table_columns( self ):
        return [ ObjectColumn( label = 'Id', name = 'lcc_id' ) ] + \
               [ ObjectColumn( label = lc.name, name = lc.name )
                for idx, lc in enumerate( self.lc_list ) ] + \
                [ ObjectColumn( label = 'assess_value', name = 'assess_value' ) ]

    geo_columns = Property( List( Str ), depends_on = 'lc_list_, +filter' )
    def _get_geo_columns( self ):
        '''derive the order of the geo columns
        from the first element in 'lc_list'. The internal
        consistency is checked separately in the
        'check_consistency' method.
        '''
        return self.lc_list[0].geo_columns

    sr_columns = Property( List( Str ), depends_on = 'lc_list_, +filter' )
    def _get_sr_columns( self ):
        '''derive the order of the stress resultants
        from the first element in 'lc_list'. The internal
        consistency is checked separately in the
        'check_consistency' method.
        '''
        return self.lc_list[0].sr_columns

    #-------------------------------
    # check consistency
    #-------------------------------

    def _check_for_consistency( self ):
        ''' check input files for consitency:
        '''
        return True

    #-------------------------------
    # lc_arr
    #-------------------------------

    lc_arr = Property( Array )
    def _get_lc_arr( self ):
        '''stack stress resultants arrays of all loading cases together.
        This yields an array of shape ( n_lc, n_elems, n_sr )
        '''
        sr_arr_list = [ lc.sr_arr for lc in self.lc_list ]
#        for x in sr_arr_list:
#            print x.shape

        return array( sr_arr_list )

    #-------------------------------
    # Array dimensions:
    #-------------------------------

    n_sr = Property( Int )
    def _get_n_sr( self ):
        return len( self.sr_columns )

    n_lc = Property( Int )
    def _get_n_lc( self ):
        return len( self.lc_list )

    n_lcc = Property( Int )
    def _get_n_lcc( self ):
        return self.combi_arr.shape[0]

    n_elems = Property( Int )
    def _get_n_elems( self ):
        return self.lc_list[0].sr_arr.shape[0]

    #-------------------------------
    # auxilary method for get_combi_arr 
    #-------------------------------

    def _product( self, args ):
        """
        Get all possible permutations of the security factors
        without changing the order of the loading cases.
        The method corresponds to the build-in function 'itertools.product'.
        Instead of returning a generator object a list of all 
        possible permutations is returned. As argument a list of list 
        needs to be defined. In the original version of 'itertools.product' 
        the function takes a tuple as argument ("*args"). 
        """
        pools = map( tuple, args ) #within original version args defined as *args
        result = [[]]
        for pool in pools:
            result = [x + [y] for x in result for y in pool]
        return result

    # ------------------------------------------------------------
    # 'combi_arr' - array containing indices of all loading case combinations:
    # ------------------------------------------------------------

    # list of indices of the position of the imposed loads in 'lc_list'
    #
#    imposed_idx_list = Property( List, depends_on = 'lc_list_, lc_list_.+input' )
    imposed_idx_list = Property( List, depends_on = 'lc_list_' )
    @cached_property
    def _get_imposed_idx_list( self ):
        '''list of indices for the imposed loads
        '''
        imposed_idx_list = []
        for i_lc, lc in enumerate( self.lc_list ):
            cat = lc.category
            if cat == 'imposed-load':
                imposed_idx_list.append( i_lc )
        return imposed_idx_list

    # array containing the psi with name 'psi_key' for the specified
    # loading cases defined in 'lc_list'. For dead-loads no value for
    # psi exists. In this case a value of 1.0 is defined. 
    # This yields an array of shape ( n_lc, )
    #
    def _get_psi_arr( self, psi_key ):
        '''psi_key must be defined as: 
        'psi_0', 'psi_1', or 'psi_2'
        Returns an 1d-array of shape ( n_lc, )
        '''
        # get list of ones (used for dead-loads):
        #
        psi_list = [1] * len( self.lc_list )

        # overwrite ones with psi-values in case of imposed-loads:
        #
        for imposed_idx in self.imposed_idx_list:
            psi_value = getattr( self.lc_list[ imposed_idx ], psi_key )
            psi_list[ imposed_idx ] = psi_value

        return array( psi_list, dtype = 'float_' )

    # list containing names of the loading cases
    #
    lc_name_list = Property( List, depends_on = 'lc_list_' )
    @cached_property
    def _get_lc_name_list( self ):
        '''list of names of all loading cases
        '''
        return [ lc.name for lc in self.lc_list ]

    show_lc_characteristic = Bool( True )

    # combination array:
    #
    combi_arr = Property( Array, depends_on = 'lc_list_, combination_SLS' )
    @cached_property
    def _get_combi_arr( self ):
        '''array containing the security and combination factors
        corresponding to the specified loading cases.
        This yields an array of shape ( n_lcc, n_lc )

        Properties defined in the subclasses 'LCCTableULS', 'LCCTableSLS':
        - 'gamma_list' = list of security factors (gamma) 
        - 'psi_lead' = combination factors (psi) of the leading imposed load 
        - 'psi_non_lead' = combination factors (psi) of the non-leading imposed loads 
        '''
        # printouts:
        #
        if self.ls == 'ULS':
            print '*** load case combinations for limit state ULS ***'
        else:
            print '*** load case combinations for limit state SLS ***'
            print '*** SLS combination used: % s ***' % ( self.combination_SLS )

        #---------------------------------------------------------------
        # get permutations of safety factors ('gamma')
        #---------------------------------------------------------------
        #
        permutation_list = self._product( self.gamma_list )

        combi_arr = array( permutation_list )

        # check if imposed loads are defined 
        # if not no further processing of 'combi_arr' is necessary:
        #
        if self.imposed_idx_list == []:

            # if option is set to 'True' the loading case combination table 
            # is enlarged with an identity matrix in order to see the 
            # characteristic values of each loading case.
            #
            if self.show_lc_characteristic:
                combi_arr = vstack( [ identity( self.n_lc ), combi_arr ] )

            return combi_arr

        #---------------------------------------------------------------
        # get leading and non leading combination factors ('psi')
        #---------------------------------------------------------------
        # go through all possible cases of leading imposed loads
        # For the currently investigated imposed loading case the
        # psi value is taken from 'psi_leading_arr' for all other 
        # imposed loads the psi value is taken from 'psi_non_lead_arr'

        # Properties are defined in the subclasses
        #
        psi_lead_arr = self.psi_lead_arr
        psi_non_lead_arr = self.psi_non_lead_arr

        # for SLS limit state case 'rare' all imposed loads are multiplied
        # with 'psi_2'. In this case no distinction between leading or 
        # non-leading imposed loads needs to be performed.  
        #
        if all( psi_lead_arr == psi_non_lead_arr ):
            combi_arr_psi = combi_arr * psi_lead_arr

        # generate a list or arrays obtained by multiplication
        # with the psi-factors.
        # This yields a list of length = number of imposed-loads.
        # 
        else:
            combi_arr_psi_list = []
            for imposed_idx in self.imposed_idx_list:
                # copy in order to preserve initial state of the array
                # and avoid in place modification
                psi_arr = copy( psi_non_lead_arr )
                psi_arr[imposed_idx] = psi_lead_arr[imposed_idx]
                combi_arr_lead_i = combi_arr[where( combi_arr[:, imposed_idx] != 0 )] * psi_arr
                combi_arr_psi_list.append( combi_arr_lead_i )

            combi_arr_psi_no_0 = vstack( combi_arr_psi_list )

            # missing cases without any dead load have to be added 
            # get combinations with all!! imposed = 0 
            #
            lcc_all_imposed_zero = where( ( combi_arr[:, self.imposed_idx_list] == 0 )
                                          .all( axis = 1 ) )

            # add to combinations
            #
            combi_arr_psi = vstack( ( combi_arr[lcc_all_imposed_zero], combi_arr_psi_no_0 ) )

        #---------------------------------------------------------------
        # get exclusive loading cases ('exclusive_to')
        #---------------------------------------------------------------

        # get a list of lists containing the indices of the loading cases
        # that are defined exclusive to each other. 
        # The list still contains duplicates, e.g. [1,2] and [2,1]
        #
        exclusive_list = []
        for i_lc, lc in enumerate( self.lc_list ):

            # get related load case number
            #
            for exclusive_name in lc.exclusive_to:
                if exclusive_name in self.lc_name_list:
                    exclusive_idx = self.lc_name_list.index( exclusive_name )
                    exclusive_list.append( [ i_lc, exclusive_idx ] )

        # eliminate the duplicates in 'exclusive_list'
        #
        exclusive_list_unique = []
        for exclusive_list_entry in exclusive_list:
            if sorted( exclusive_list_entry ) not in exclusive_list_unique:
                exclusive_list_unique.append( sorted( exclusive_list_entry ) )

        # delete the rows in combination array that contain
        # loading case combinations with imposed-loads that have been defined
        # as exclusive to each other. 
        # 
        combi_arr_psi_exclusive = combi_arr_psi
#        print 'combi_arr_psi_exclusive', combi_arr_psi_exclusive
        for exclusive_list_entry in exclusive_list_unique:
            # check where maximum one value of the exclusive load cases is unequal to one 
            #              LC1  LC2  LC3  (all LCs are defined as exclusive to each other)
            #
            # e.g.         1.5  0.9  0.8  (example of 'combi_arr_psi') 
            #              1.5  0.0  0.0  
            #              0.0  0.0  0.0  (combination with all imposed loads = 0 after multiplication wit psi and gamma)
            #              ...  ...  ...
            #
            # this would yield the following mask_arr (containing ones or zeros):
            # e.g.         1.0  1.0  1.0  --> sum = 3 --> true combi --> accepted combination
            #              1.0  0.0  0.0  --> sum = 1 --> false combi --> no accepted combination
            # e.g.         0.0  0.0  0.0  --> sum = 0 --> true combi --> accepted combination (only body-loads)
            #              ...  ...  ...
            #
            mask_arr = where( combi_arr_psi_exclusive[ :, exclusive_list_entry ] != 0, 1.0, 0.0 )
#            print 'mask_arr', mask_arr
            true_combi = where( sum( mask_arr, axis = 1 ) <= 1.0 )
#            print 'true_combi', true_combi
            combi_arr_psi_exclusive = combi_arr_psi_exclusive[ true_combi ]

        #---------------------------------------------------------------
        # create array with only unique load case combinations 
        #---------------------------------------------------------------
        # If the psi values of an imposed-load are defined as zero this
        # may led to zero entries in 'combi_arr'. This would yield rows
        # in 'combi_arr' which are duplicates. Those rows are removed.

        # Add first row in 'combi_arr_psi_exclusive' to '_unique' array 
        # This array must have shape (1, n_lc) in order to use 'axis'-option
        #
        combi_arr_psi_exclusive_unique = combi_arr_psi_exclusive[0][None, :]

        for row in combi_arr_psi_exclusive:
            # Check if all factors in one row are equal to the rows in 'unique' array.
            # If this is not the case for any row the combination is added to 'unique'.
            # Broadcasting is used for the bool evaluation:
            #
            if ( row == combi_arr_psi_exclusive_unique ).all( axis = 1.0 ).any() == False:
                combi_arr_psi_exclusive_unique = vstack( ( combi_arr_psi_exclusive_unique, row ) )

        # if option is set to 'True' the loading case combination table 
        # is enlarged with an identity matrix in order to see the 
        # characteristic values of each loading case.
        #
#        if self.show_lc_characteristic:
#            combi_arr_psi_exclusive_unique = vstack( [ identity( self.n_lc ), combi_arr_psi_exclusive_unique ] )

        return combi_arr_psi_exclusive_unique

    #-------------------------------
    # lcc_arr 
    #-------------------------------

    lcc_arr = Property( Array, depends_on = 'lc_list_' )
    @cached_property
    def _get_lcc_arr( self ):
        '''Array of all loading case combinations following the
        loading cases define in 'lc_list' and the combinations
        defined in 'combi_arr'.
        This yields an array of shape ( n_lcc, n_elems, n_sr )
        '''
        self._check_for_consistency()

        combi_arr = self.combi_arr

        # 'combi_arr' is of shape ( n_lcc, n_lc )
        # 'lc_arr' is of shape ( n_lc, n_elems, n_sr )
        #
        lc_arr = self.lc_arr

        # Broadcasting is used to generate array containing the multiplied lc's
        # yielding an array of shape ( n_lcc, n_lc, n_elems, n_sr )
        #
        lc_combi_arr = lc_arr[ None, :, :, :] * combi_arr[:, :, None, None ]

        # Then the sum over index 'n_lc' is evaluated yielding 
        # an array of all loading case combinations.
        # This yields an array of shape ( n_lcc, n_elem, n_sr ) 
        #
        lcc_arr = sum( lc_combi_arr, axis = 1 )

        return lcc_arr


    #-------------------------------
    # geo_arr 
    #-------------------------------

    geo_data_dict = Property( Dict, depends_on = 'lc_list_' )
    @cached_property
    def _get_geo_data_dict( self ):
        '''Array of global coords derived from the first loading case defined in lc_list.
        Coords are identical for all LC's.
        '''
        return self.lc_list[0].geo_data_dict


    #-------------------------------
    # min/max-values 
    #-------------------------------

    def get_min_max_state_data( self ):
        ''' get the surrounding curve of all 'lcc' values
        '''
        lcc_arr = self.lcc_arr
        min_arr = ndmin( lcc_arr, axis = 0 )
        max_arr = ndmax( lcc_arr, axis = 0 )
        return min_arr, max_arr

    max_sr_grouped_dict = Property( Dict )
    @cached_property
    def _get_max_sr_grouped_dict( self ):
        ''' get the surrounding curve for each stress resultant
            shape lcc_array ( n_lcc, n_elems, n_sr )
        '''
        sr_columns = self.sr_columns
        lcc_arr = self.lcc_arr
        dict = {}
        for i, sr in enumerate( self.sr_columns ):
            idx_1 = argmax( abs( lcc_arr[:, :, i] ), axis = 0 )
            idx_2 = arange( 0, idx_1.shape[0], 1 )
            dict[sr] = lcc_arr[idx_1, idx_2, :]
        return dict

    # choose linking type (in-plane shear dof blocked or not)
    #
    link_type = Enum( 'exc_V_ip', 'inc_V_ip' )

    # length of the shell (needed to plot the hinge forces plots correctly)
    #
    length_xy_quarter = 3.5 # m

    def export_hf_max_grouped( self, filename ):
        """exports with one leading value
        """
        from matplotlib import pyplot
        sr_columns = self.sr_columns
        dict = self.max_sr_grouped_dict
        length_xy_quarter = self.length_xy_quarter

        def save_bar_plot( x, y, filename = 'bla', title = 'Title',
                            xlabel = 'xlabel', ylabel = 'ylavel',
                            width = 0.1, xmin = 0, xmax = 1000 ,
                            ymin = -1000, ymax = 1000, figsize = [10, 5] ):
            fig = pyplot.figure( facecolor = "white", figsize = figsize )
            ax1 = fig.add_subplot( 1, 1, 1 )
            ax1.bar( x , y , width = width, align = 'center', color = 'green' )
            ax1.set_xlim( xmin, xmax )
            ax1.set_ylim( ymin, ymax )
            ax1.set_xlabel( xlabel, fontsize = 22 )
            ax1.set_ylabel( ylabel, fontsize = 22 )
            if title == 'N_ip max':
                title = 'max $N_{ip}$'
            if title == 'V_ip max':
                title = 'max $V_{ip}$'
            if title == 'V_op max':
                title = 'max $V_{op}$'
            ax1.set_title( title )
            fig.savefig( filename, orientation = 'portrait', bbox_inches = 'tight' )
            pyplot.clf()


        X = array( self.geo_data_dict['X_hf'] )
        Y = array( self.geo_data_dict['Y_hf'] )

        # symmetric axes
        #
        idx_sym = where( abs( Y[:, 0] - 2.0 * length_xy_quarter ) <= 0.0001 )
        X_sym = X[idx_sym].reshape( -1 )
        idx_r0_r1 = where( abs( X[:, 0] - 2.0 * length_xy_quarter ) <= 0.0001 )
        X_r0_r1 = Y[idx_r0_r1].reshape( -1 )

        for sr in sr_columns:
            F_int = dict[sr]    #first row N_ip, second V_ip third V_op
            F_sym = F_int[idx_sym, :].reshape( -1, len( sr_columns ) )
            F_r0_r1 = F_int[idx_r0_r1, :].reshape( -1, len( sr_columns ) )


            save_bar_plot( X_sym, F_sym[:, 0].reshape( -1 ),
                          xlabel = '$X$ [m]', ylabel = '$N_{ip}$ [kN]',
                          filename = filename + 'N_ip' + '_sym_' + sr + '_max',
                          title = sr + ' max',
                          xmin = 0.0, xmax = 4.0 * length_xy_quarter, figsize = [10, 5], ymin = -40, ymax = +40 )
            if self.link_type == 'inc_V_ip':
                save_bar_plot( X_sym, F_sym[:, 1].reshape( -1 ),
                              xlabel = '$X$ [m]', ylabel = '$V_{ip}$ [kN]',
                              filename = filename + 'V_ip' + '_sym_' + sr + '_max',
                              title = sr + ' max',
                              xmin = 0.0, xmax = 4.0 * length_xy_quarter, figsize = [10, 5], ymin = -40, ymax = +40 )


            save_bar_plot( X_sym, F_sym[:, 2].reshape( -1 ),
                          xlabel = '$X$ [m]', ylabel = '$V_{op}$ [kN]',
                          filename = filename + 'V_op' + '_sym_' + sr + '_max',
                          title = sr + ' max',
                          xmin = 0.0, xmax = 2.0 * length_xy_quarter, figsize = [10, 5], ymin = -10, ymax = +10 )

            # r0_r1
            #
            save_bar_plot( X_r0_r1, F_r0_r1[:, 0].reshape( -1 ),
                          xlabel = '$Y$ [m]', ylabel = '$N_{ip}$ [kN]',
                          filename = filename + 'N_ip' + '_r0_r1_' + sr + '_max',
                           title = sr + ' max',
                          xmin = 0.0, xmax = 2.0 * length_xy_quarter, figsize = [5, 5], ymin = -40, ymax = +40 )
            if self.link_type == 'inc_V_ip':
                save_bar_plot( X_r0_r1, F_r0_r1[:, 1].reshape( -1 ),
                              xlabel = '$Y$ [m]', ylabel = '$V_{ip}$ [kN]',
                              filename = filename + 'V_ip' + '_r0_r1_' + sr + '_max',
                              title = sr + ' max',
                              xmin = 0.0, xmax = 2.0 * length_xy_quarter, figsize = [5, 5], ymin = -40, ymax = +40 )
            save_bar_plot( X_r0_r1, F_r0_r1[:, 2].reshape( -1 ),
                          xlabel = '$Y$ [m]', ylabel = '$V_{op}$ [kN]',
                          filename = filename + 'V_op' + '_r0_r1_' + sr + '_max',
                          title = sr + ' max',
                          xmin = 0.0, xmax = 2.0 * length_xy_quarter, figsize = [5, 5], ymin = -10, ymax = +10 )

    def export_u( self, filename ):
        """exports u values, maximum and minimum of each lc and combination
        """
        length_xy_quarter = self.length_xy_quarter

        n_lc = self.n_lc
        n_sr = self.n_sr
        # max value charakteristic
        #
        max_export_arr = zeros( ( n_lc, n_sr ) )
        min_export_arr = zeros( ( n_lc, n_sr ) )


        for i in range( 0, n_lc ):
            for j in range( 0, n_sr ):
                max_export_arr[i, j] = ndmax( self.lc_arr[i, :, j] )
                min_export_arr[i, j] = ndmin( self.lc_arr[i, :, j] )

        # from combinated values
        #
        # from outer_edges (= r0_r1_bottom, r0_left, r1_right)
        #
        idx_bottom = where( abs( self.geo_data_dict['Y_u'].flatten() ) <= 0.001 )
        idx_left = where( abs( self.geo_data_dict['X_u'].flatten() ) <= 0.001 )
        idx_right = where( abs( self.geo_data_dict['X_u'].flatten() - 4.0 * length_xy_quarter ) <= 0.001 )
        idx_out = unique( hstack( [idx_bottom, idx_left, idx_right] ) )
#        print 'idx_out', idx_out

        # for internal edge
        #
        idx_sym = where( abs( self.geo_data_dict['Y_u'].reshape( -1 ) - 2.0 * length_xy_quarter ) <= 0.001 )
        idx_r0_r1 = where( abs( self.geo_data_dict['X_u'].reshape( -1 ) - 2.0 * length_xy_quarter ) <= 0.001 )
        idx_int = unique( hstack( [idx_sym, idx_left, idx_r0_r1] ) )

        # look in lcc_arr (= array with the combined values)
        # for the combination with the maximum or minimum combined value 
        # NOTE: lcc_arr.shape = ( n_lcc, n_elems, n_sr )
        #
        max_komb_out = zeros( ( n_sr ) )
        min_komb_out = zeros( ( n_sr ) )
        for j in range( 0, n_sr ):
            max_komb_out[j] = ndmax( lct.lcc_arr[:, idx_out, j].reshape( -1 ), axis = 0 )
            min_komb_out[j] = ndmin( lct.lcc_arr[:, idx_out, j].reshape( -1 ), axis = 0 )

            print"\n"
            print"-------------------"
            print self.sr_columns[j] + " - MAX"
            print"-------------------"
            idx_max = where ( lct.lcc_arr[:, :, j] == max_komb_out[j] )
            print 'idx_max', idx_max
            idx_comb = idx_max[0][0]
            print 'idx_comb', idx_comb
            idx_point = idx_max[1][0]
            print 'idx_point', idx_point
            idx_sr = j
            print 'max_komb_out', max_komb_out[j]
            print 'combi_arr', self.combi_arr[idx_comb, :]
            print 'lc_arr', self.lc_arr[:, idx_point, j]
            print self.lc_arr[:, idx_point, j] * self.combi_arr[idx_comb, :] / max_komb_out[j]
            print "design value ", vdot( self.lc_arr[:, idx_point, j], self.combi_arr[idx_comb, :] )

            print 'at position X,Y,Z', self.geo_data_dict['X_u'][idx_point], self.geo_data_dict['Y_u'][idx_point], self.geo_data_dict['Z_u'][idx_point]

            print"-------------------"
            print self.sr_columns[j] + " - MIN"
            print"-------------------"
            idx_min = where ( lct.lcc_arr[:, :, j] == min_komb_out[j] )
            print 'idx_min', idx_min
            idx_comb = idx_min[0][0]
            print 'idx_comb', idx_comb
            idx_point = idx_min[1][0]
            print 'idx_point', idx_point
            idx_sr = j

            print min_komb_out[j]
            print self.combi_arr[idx_comb, :]
            print self.lc_arr[:, idx_point, j]
            print self.lc_arr[:, idx_point, j] * self.combi_arr[idx_comb, :] / min_komb_out[j]
            print "design value ", vdot( self.combi_arr[idx_comb, :], self.lc_arr[:, idx_point, j] )
            print 'at position X,Y,Z', self.geo_data_dict['X_u'][idx_point], self.geo_data_dict['Y_u'][idx_point], self.geo_data_dict['Z_u'][idx_point]

        max_komb_int = zeros( ( n_sr ) )
        min_komb_int = zeros( ( n_sr ) )
        for j in range( 0, n_sr ):

            max_komb_int[j] = ndmax( lct.lcc_arr[:, idx_int, j].reshape( -1 ), axis = 0 )
            min_komb_int[j] = ndmin( lct.lcc_arr[:, idx_int, j].reshape( -1 ), axis = 0 )



        max_export_arr = vstack( ( max_export_arr, max_komb_out, max_komb_int ) )
        min_export_arr = vstack( ( min_export_arr, min_komb_out, min_komb_int ) )

        def csv_u( data, filename = 'U_data.csv' ):
            '''exports X_U_export and U_export data to csv - worksheet
            '''
            file = open( filename, 'w' )

            writer = csv.writer( file, delimiter = ";", lineterminator = "\n" )

            #first row
            #
            writer.writerow( ['-'] + self.sr_columns )

            # not combinated rows
            #
            for i in range( 0, self.n_lc ):
                a = [lc_list[i].name] + list( data[i] )
                writer.writerow( a )

            # combinated rows
            #
            writer.writerow( ['komb_out'] + list( data[-2] ) )
            writer.writerow( ['komb_int'] + list( data[-1] ) )


            file = file.close()

        csv_u( max_export_arr, filename = filename + 'max_U' + '.csv' )
        csv_u( min_export_arr, filename = filename + 'min_U' + '.csv' )


    def plot_interaction_tp( self ):
        """get interaction pairs max"""

        lcc_arr = self.lcc_arr

        N_tp_d = ( lcc_arr[:, :, 0] ** 2 ) ** 0.5
        V_tp_d = lcc_arr[:, :, 2]

        # s6cm
        #
        N_ck = 66.8
        V_ck = 8.5

        gamma_s = 3.0

        beta_N = N_tp_d / ( N_ck / gamma_s )
        beta_V = abs( V_tp_d / ( V_ck / gamma_s ) )
        beta_inter = ( beta_N ) ** 1.5 + ( beta_V ) ** 1.5

        idx_max_hinge = beta_inter.argmax( axis = 0 )
        idx_hinge = arange( 0, len( idx_max_hinge ), 1 )
        plot_beta_N = beta_N[idx_max_hinge, idx_hinge]
        plot_beta_V = beta_V[idx_max_hinge, idx_hinge]
        self.interaction_plot( plot_beta_N, plot_beta_V )



    def plot_interaction_s6cm( self ):
        """get interaction pairs max"""

        lcc_arr = self.lcc_arr

        N_s6cm_d = lcc_arr[:, :, 2] * 1.5
        V_s6cm_d = ( ( lcc_arr[:, :, 0] / 2 ) ** 2 + ( lcc_arr[:, :, 1] * 1.5 ) ** 2 ) ** 0.5

        # s6cm
        #
        N_ck = 29.2
        V_ck = 67.5

        gamma_s = 3.0

        beta_N = N_s6cm_d / ( N_ck / gamma_s )
        beta_V = abs( V_s6cm_d / ( V_ck / gamma_s ) )
        beta_inter = ( beta_N ) ** 1.5 + ( beta_V ) ** 1.5

        idx_max_hinge = beta_inter.argmax( axis = 0 )
        idx_hinge = arange( 0, len( idx_max_hinge ), 1 )
        plot_beta_N = beta_N[idx_max_hinge, idx_hinge]
        plot_beta_V = beta_V[idx_max_hinge, idx_hinge]
        self.interaction_plot( plot_beta_N, plot_beta_V )



    def interaction_plot( self, eta_N, eta_V ):
            from matplotlib import pyplot
            fig = pyplot.figure( facecolor = "white", figsize = [10, 10] )
            ax1 = fig.add_subplot( 1, 1, 1 )
            x = arange( 0, 1.01, 0.01 )
            y = ( 1 - x ** 1.5 ) ** ( 1 / 1.5 )
            limit = eta_N ** 1.5 + eta_V ** 1.5


            ax1.set_xlabel( '$V_{Ed}/V_{Rd}$' , fontsize = 24 )
            ax1.set_ylabel( '$N_{Ed}/N_{Rd}$', fontsize = 24 )
            ax1.plot( x , y, '--', color = 'black'
                      , linewidth = 2.0 )
            ax1.plot( eta_V[where( limit < 1 )] , eta_N[where( limit < 1 )], 'o', markersize = 8 )
            ax1.plot( eta_V[where( limit > 1 )] , eta_N[where( limit > 1 )], 'o', color = 'red', markersize = 8 )

            for xlabel_i in ax1.get_xticklabels():
                xlabel_i.set_fontsize( 18 )

            for ylabel_i in ax1.get_yticklabels():
                ylabel_i.set_fontsize( 18 )


    #        ax1.plot( x , 1 - x, '--', color = 'black', label = 'lineare Interaktion' )

            ax1.set_xlim( 0, 2.0 )
            ax1.set_ylim( 0, 2.0 )
            ax1.legend()
            pyplot.show()
            pyplot.clf()



    def export_hf_lc( self ):
        """exports with one leading value
        """

        from matplotlib import pyplot
        sr_columns = self.sr_columns
        dict = self.max_sr_grouped_dict
        length_xy_quarter = self.length_xy_quarter

        def save_bar_plot( x, y, filename = 'bla',
                              xlabel = 'xlabel', ylabel = 'ylavel', ymin = -10 , ymax = 10,
                              width = 0.1, xmin = 0, xmax = 1000 , figsize = [10, 5] ):
            fig = pyplot.figure( facecolor = "white", figsize = figsize )
            ax1 = fig.add_subplot( 1, 1, 1 )
            ax1.bar( x , y , width = width, align = 'center', color = 'blue' )
            ax1.set_xlim( xmin, xmax )
            ax1.set_ylim( ymin, ymax )
            ax1.set_xlabel( xlabel, fontsize = 22 )
            ax1.set_ylabel( ylabel, fontsize = 22 )
            fig.savefig( filename, orientation = 'portrait', bbox_inches = 'tight' )
            pyplot.clf()


        X = array( self.geo_data_dict['X_hf'] )
        Y = array( self.geo_data_dict['Y_hf'] )

        # symmetric axes
        #
        idx_sym = where( abs( Y[:, 0] - 2.0 * length_xy_quarter ) <= 0.0001 )
        X_sym = X[idx_sym].reshape( -1 )
        idx_r0_r1 = where( abs( X[:, 0] - 2.0 * length_xy_quarter ) <= 0.0001 )
        X_r0_r1 = Y[idx_r0_r1].reshape( -1 )
        F_int = self.lc_arr

        for i, lc_name in enumerate( self.lc_name_list ):
            filename = self.lc_list[i].plt_export

            max_N_ip = max( int( ndmax( F_int[i, :, 0], axis = 0 ) ) + 1, 1 )
            max_V_ip = max( int( ndmax( F_int[i, :, 1], axis = 0 ) ) + 1, 1 )
            max_V_op = max( int( ndmax( F_int[i, :, 2], axis = 0 ) ) + 1, 1 )


            F_int_lc = F_int[i, :, :]   #first row N_ip, second V_ip third V_op
            F_sym = F_int_lc[idx_sym, :].reshape( -1, len( sr_columns ) )
            F_r0_r1 = F_int_lc[idx_r0_r1, :].reshape( -1, len( sr_columns ) )


            save_bar_plot( X_sym, F_sym[:, 0].reshape( -1 ),
                          xlabel = '$X$ [m]', ylabel = '$N_{ip}$ [kN]',
                          filename = filename + 'N_ip' + '_sym',
                          xmin = 0.0, xmax = 4.0 * length_xy_quarter,
                          ymin = -max_N_ip, ymax = max_N_ip, figsize = [10, 5] )

            save_bar_plot( X_sym, F_sym[:, 1].reshape( -1 ),
                          xlabel = '$X$ [m]', ylabel = '$V_{ip}$ [kN]',
                          filename = filename + 'V_ip' + '_sym',
                          xmin = 0.0, xmax = 4.0 * length_xy_quarter,
                          ymin = -max_V_ip, ymax = max_V_ip, figsize = [10, 5] )


            save_bar_plot( X_sym, F_sym[:, 2].reshape( -1 ),
                          xlabel = '$X$ [m]', ylabel = '$V_{op}$ [kN]',
                          filename = filename + 'V_op' + '_sym',
                          xmin = 0.0, xmax = 4.0 * length_xy_quarter,
                          ymin = -max_V_op, ymax = max_V_op, figsize = [10, 5] )


            # r0_r1
            #
            save_bar_plot( X_r0_r1, F_r0_r1[:, 0].reshape( -1 ),
                          xlabel = '$Y$ [m]', ylabel = '$N_{ip}$ [kN]',
                          filename = filename + 'N_ip' + '_r0_r1',
                          xmin = 0.0, xmax = 2.0 * length_xy_quarter,
                          ymin = -max_N_ip, ymax = max_N_ip, figsize = [5, 5] )
            save_bar_plot( X_r0_r1, F_r0_r1[:, 1].reshape( -1 ),
                          xlabel = '$Y$ [m]', ylabel = '$V_{ip}$ [kN]',
                          filename = filename + 'V_ip' + '_r0_r1',
                          xmin = 0.0, xmax = 2.0 * length_xy_quarter,
                          ymin = -max_V_ip, ymax = max_V_ip, figsize = [5, 5] )
            save_bar_plot( X_r0_r1, F_r0_r1[:, 2].reshape( -1 ),
                          xlabel = '$Y$ [m]', ylabel = '$V_{op}$ [kN]',
                          filename = filename + 'V_op' + '_r0_r1',
                          xmin = 0.0, xmax = 2.0 * length_xy_quarter,
                          ymin = -max_V_op, ymax = max_V_op, figsize = [5, 5] )

    #-------------------------------
    # lcc_lists 
    #-------------------------------

    lcc_list = Property( List, depends_on = 'lc_list_' )
    @cached_property
    def _get_lcc_list( self ):
        '''list of loading case combinations (instances of LCC)
        '''
        combi_arr = self.combi_arr
        lcc_arr = self.lcc_arr
        sr_columns = self.sr_columns
        geo_columns = self.geo_columns

        n_lcc = self.n_lcc

#        print 'combi_arr', combi_arr
#        print 'lcc_arr', lcc_arr
#        print 'sr_columns', sr_columns
#        print 'geo_columns', geo_columns
#        print 'n_lcc', n_lcc

        # return a dictionary of the stress resultants
        # this is used by LSTable to determine the stress 
        # resultants of the current limit state 
        #
        lcc_list = []
        for i_lcc in range( n_lcc ):

            state_data_dict = {}
            for i_sr, name in enumerate( sr_columns ):
                state_data_dict[ name ] = lcc_arr[ i_lcc, :, i_sr ][:, None]

            geo_data_dict = self.geo_data_dict

            lcc = LCC( #lcc_table = self,
                       factors = combi_arr[ i_lcc, : ],
                       lcc_id = i_lcc,
                       ls_table = LSTable( geo_data = geo_data_dict,
                                           state_data = state_data_dict,
                                           ls = self.ls )
                       )

            for idx, lc in enumerate( self.lc_list ):
                lcc.add_trait( lc.name, Int( combi_arr[ i_lcc, idx ] ) )

            lcc_list.append( lcc )

        return lcc_list





    # ------------------------------------------------------------
    # View 
    # ------------------------------------------------------------

    traits_view = View( VGroup( 

                        VSplit( 
                                    Item( 'lcc_list', editor = lcc_list_editor,
                                          show_label = False ),
                                    Item( 'lcc@', show_label = False ),
                                    ),
                        ),
                      resizable = True,
                      scrollable = True,
                      height = 1.0,
                      width = 1.0
                      )


class LCCTableULS( LCCTable ):
    '''LCCTable for ultimate limit state
    '''

    # set limit state to 'ULS'
    # (attribute is used by 'LSTable')
    #
    ls = 'ULS'

    # 'gamma' - safety factors
    #
    gamma_list = Property( List, depends_on = 'lc_list_' )
    @cached_property
    def _get_gamma_list( self ):
        return [[ lc.gamma_fav, lc.gamma_unf ] for lc in self.lc_list ]

    # 'psi' - combination factors (psi) for leading
    # and non leading load cases
    #
    psi_non_lead_arr = Property( Array, depends_on = 'lc_list_' )
    @cached_property
    def _get_psi_non_lead_arr( self ):
        return self._get_psi_arr( 'psi_0' )

    psi_lead_arr = Property( Array, depends_on = 'lc_list_' )
    @cached_property
    def _get_psi_lead_arr( self ):
        return ones( len( self.lc_list ) )


class LCCTableSLS( LCCTable ):
    '''LCCTable for serviceability limit state
    '''

    # set limit state to 'SLS'
    # (attribute is used by 'LSTable')
    #
    ls = 'SLS'

    # possible definitions of the serviceability limit state    
    # are: 'rare', 'freq', 'perm'
    #
    combination_SLS = Enum( 'rare', 'freq', 'perm' )
    def _combination_SLS_default( self ):
        return 'rare'

    # 'gamma' - safety factors
    #
    gamma_list = Property( List, depends_on = 'lc_list_' )
    @cached_property
    def _get_gamma_list( self ):

        # generate [1.0]-entry in case of body-loads:
        #
        gamma_list = [[ 1.0 , 3.0]] * len( self.lc_list ) # for creeping
#        gamma_list = [[ 1.0]] * len( self.lc_list ) # for creeping

        # overwrite those in case of imposed-loads:
        #
        for imposed_idx in self.imposed_idx_list:
            gamma_fav_SLS = getattr( self.lc_list[ imposed_idx ], 'gamma_fav_SLS' )
            gamma_unf_SLS = getattr( self.lc_list[ imposed_idx ], 'gamma_unf_SLS' )
            gamma_list[ imposed_idx ] = [ gamma_unf_SLS, gamma_fav_SLS ]

        return gamma_list

    # 'psi' - combination factors
    #
    psi_lead_dict = Property( Dict )
    def _get_psi_lead_dict( self ):
        return {'rare' : ones_like( self._get_psi_arr( 'psi_0' ) ) ,
                'freq' : self._get_psi_arr( 'psi_1' ),
                'perm' : self._get_psi_arr( 'psi_2' )}

    psi_non_lead_dict = Property( Dict )
    def _get_psi_non_lead_dict( self ):
        return {'rare' : self._get_psi_arr( 'psi_0' ) ,
                'freq' : self._get_psi_arr( 'psi_2' ),
                'perm' : self._get_psi_arr( 'psi_2' )}

    # combination factors (psi) for leading
    # and non leading load cases
    #
    psi_lead_arr = Property( Array, depends_on = 'lc_list_, combination_SLS' )
    @cached_property
    def _get_psi_lead_arr( self ):
        return self.psi_lead_dict[ self.combination_SLS ]

    psi_non_lead_arr = Property( Array, depends_on = 'lc_list_, combination_SLS' )
    @cached_property
    def _get_psi_non_lead_arr( self ):
        return self.psi_non_lead_dict[ self.combination_SLS ]


if __name__ == '__main__':

    #---------------------------------------------
    # 2 shells: 
    # new geometry 7m x 7m
    #---------------------------------------------

    #------------------------
    # evaluate the combination for the displpacements (SLS) or the hinge forces (ULS)
    #------------------------

    # choose!
    #
#    do = 'hf'
    do = 'u'

    if do == 'u':
        # NOTE: switch those values from hf to u in the source code directly!
        # !!!!!!!    sr_columns and geo_columns need to be changed for u and hf options
        #            for hf : N_IP..... 
        #            for u: u_z .....
        #            as state in quellcode above

        # NOTE: switch those values from hf to u in the source code directly!
        sr_columns = ['U_x', 'U_y', 'U_z']
#        sr_columns = ['U_x', 'U_y', 'U_z', 'dU_y', 'dU_z']
        geo_columns = ['X_u', 'Y_u', 'Z_u']

    if do == 'hf':
        sr_columns = ['N_ip', 'V_ip', 'V_op']
        geo_columns = ['X_hf', 'Y_hf', 'Z_hf']

    #------------------------
    # specify linking option:
    #------------------------

#    link_case = 'equal_88cm'
    link_case = 'equal_100cm_7m'

    link_type = 'exc_V_ip'
#        link_type = 'inc_V_ip'

    spring = 'no_spring'
#        spring = 'spring'

    if spring == 'spring':
        spring_type = 'spring_1_'
    else:
        spring_type = ''

    #------------------------
    # specify directory containig csv-files of the hinge forces and the displacements:
    #------------------------

    data_dir = os.path.join( simdb.simdb_dir,
                             'simdata', 'output_data_mushroof',
                             'hinge_forces_and_displacement_2shells',
                             spring,
                             link_type,
                             link_case )

    #------------------------
    # specify loading cases (names, gamma, psi):
    #------------------------

    lc_list = [
                 # NOTE: in order to define the dead-loads and additional dead-loads exclusive to each other
                 # it is necessary to define these loading cases as "imposed load" and specify the gamma-values 
                 # explicitly for dead loads and additional dead-loads

                 # own weight, and additional loads:
                 LC( name = 'G', category = 'dead-load',
                     file_name = os.path.join( data_dir , spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_g.csv' ),
                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_g_',
                     #
#                     gamma_fav = 1.0, gamma_unf = 1.35,
                     gamma_fav_SLS = 1.0, gamma_unf_SLS = 3.0,
                     psi_0 = 1.0, psi_1 = 1.0, psi_2 = 1.0,
                     exclusive_to = ['g_own_weight'], #, 'g_edge_load' , 'g_tol_asym', 'g_tol_sym', 'g_surf_load'],
                     #  
                     sr_columns = sr_columns,
                     geo_columns = geo_columns
                     ),

#                 # own weight without additional loads:
#                 # NOTE: defined as "imposed-load" with predefined 'psi' and 'gamma' values
#                 # in order to define load as "exclusive_to"
#                 #
#                 LC( name = 'g_own_weight', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_g_own_weight.csv' ),
#                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_g_own_weight_',
#                     #
#                     gamma_fav = 1.0, gamma_unf = 1.35,
#                     gamma_fav_SLS = 1.0, gamma_unf_SLS = 3.0,
#                     psi_0 = 1.0, psi_1 = 1.0, psi_2 = 1.0,
#                     exclusive_to = ['G', 'g_edge_load', 'g_tol_asym', 'g_tol_sym', 'g_surf_load'],
#                     #
#                     sr_columns = sr_columns,
#                     geo_columns = geo_columns
#                     ),
##
#                 # additional thickness (symmetric):
#                 #
#                 LC( name = 'g_tol_asym', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_g_tol_asym.csv' ),
#                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_g_tol_asym_',
#                     #
#                     gamma_fav = 0.0, gamma_unf = 1.35,
#                     gamma_fav_SLS = 1.0, gamma_unf_SLS = 3.0,
#                     psi_0 = 1.0, psi_1 = 1.0, psi_2 = 1.0,
#                     exclusive_to = ['G', 'g_edge_load', 'g_own_weight', 'g_tol_sym', 'g_surf_load'],
#                     #
#                     sr_columns = sr_columns,
#                     geo_columns = geo_columns
#                     ),
#
#                 # additional thickness (asymmetric):
#                 #
#                 LC( name = 'g_tol_sym', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_g_tol_sym.csv' ),
#                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_g_tol_sym_',
#                     #
#                     gamma_fav = 0.0, gamma_unf = 1.35,
#                     gamma_fav_SLS = 1.0, gamma_unf_SLS = 3.0,
#                     psi_0 = 1.0, psi_1 = 1.0, psi_2 = 1.0,
#                     exclusive_to = ['G', 'g_edge_load', 'g_own_weight', 'g_tol_asym', 'g_surf_load'],
#                     #
#                     sr_columns = sr_columns,
#                     geo_columns = geo_columns
#                     ),
#
#                 # permanent fassade loads (edge):
#                 #
#                 LC( name = 'g_edge_load', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_g_edge_load.csv' ),
#                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_g_edge_load_',
#                     #
#                     gamma_fav = 1.0, gamma_unf = 1.35,
#                     gamma_fav_SLS = 1.0, gamma_unf_SLS = 3.0,
#                     psi_0 = 1.0, psi_1 = 1.0, psi_2 = 1.0,
#                     exclusive_to = ['G'],
#                     #
#                     sr_columns = sr_columns,
#                     geo_columns = geo_columns
#                     ),

#                 # surface loads (roof):
#                 #
#                 LC( name = 'g_surf_load', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_g_surf_load.csv' ),
#                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_g_surf_load_',
#                     #
#                     gamma_fav = 1.0, gamma_unf = 1.35,
#                     gamma_fav_SLS = 1.0, gamma_unf_SLS = 3.0,
#                     psi_0 = 1.0, psi_1 = 1.0, psi_2 = 1.0,
#                     exclusive_to = ['G', 'g_own_weight', 'g_tol_asym', 'g_tol_sym', 'g_edge_load', 'g_surf_load'],
#                     #
#                     sr_columns = sr_columns,
#                     geo_columns = geo_columns
#                     ),

                 # snow symmetric:
                 #
                 LC( name = 'S_sym', category = 'imposed-load',
                     file_name = os.path.join( data_dir , spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_s_sym.csv' ),
                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_s_sym_' ,
                     #
                     psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0 ,
                     exclusive_to = ['S_asym'],
                     #
                     sr_columns = sr_columns,
                     geo_columns = geo_columns
                     ),

                 # snow asymmetric:
                 #
                 LC( name = 'S_asym', category = 'imposed-load',
                     file_name = os.path.join( data_dir, spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_s_asym.csv' ),
                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_s_asym_' ,
                     #
                     exclusive_to = ['S_sym'],
                     psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0 ,
                     #
                     sr_columns = sr_columns,
                     geo_columns = geo_columns
                     ),

                 # wind asymmetric
                 #
                 LC( name = 'W_asym', category = 'imposed-load',
                     file_name = os.path.join( data_dir, spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_w_asym.csv' ),
                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_w_asym_' ,
                     #
                     exclusive_to = ['W_pos', 'W_neg'],
                     psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0 ,
                     #
                     sr_columns = sr_columns,
                     geo_columns = geo_columns
                     ),

                 # wind pos
                 #
                 LC( name = 'W_pos', category = 'imposed-load',
                     file_name = os.path.join( data_dir, spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_w_pos.csv' ),
                     plt_export = link_type + '_' + link_case + '_plt_' + 'lc_w_pos_' ,
                     #
                     exclusive_to = ['W_asym', 'W_neg'],
                     psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0 ,
                     #
                     sr_columns = sr_columns,
                     geo_columns = geo_columns
                     ),

                 # wind neg
                 #
                 LC( name = 'W_neg', category = 'imposed-load',
                     file_name = os.path.join( data_dir, spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_w_neg.csv' ),
                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_w_neg_' ,
                     #
                     exclusive_to = ['W_pos', 'W_asym'], psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0 ,
                     #
                     sr_columns = sr_columns,
                     geo_columns = geo_columns
                     ),

                # shrinkage
                #s
                LC( name = 'T_shrinkage', category = 'imposed-load',
                     file_name = os.path.join( data_dir, spring_type + link_type + '_' + link_case + '_' + do + '_' + 'lc_shrink.csv' ),
                     plt_export = spring_type + link_type + '_' + link_case + '_plt_' + 'lc_shrink_' ,
                     #
                     exclusive_to = [], psi_0 = 0.8, psi_1 = 0.7, psi_2 = 0.5 ,
                     #
                     sr_columns = sr_columns,
                     geo_columns = geo_columns
                     ),
               ]


    #-------------------------------------------------------
    # evaluate the displacement files csv:
    #-------------------------------------------------------
    #
    if do == 'u':

        lct = LCCTableSLS( data_dir = data_dir,
                            lc_list = lc_list,
        #                       cut_z_fraction = 0.2,
                            combination_SLS = 'rare',
        ##                       combination_SLS = 'freq',
        ##                       combination_SLS = 'perm',
        ##                       show_lc_characteristic = True
                            )

        # export of the max edge displacement 
        #
        lct.export_u( filename = spring + "_" + link_type + "_" + link_case )

        # LCC-TABLE
#        lct.configure_traits()


    #-------------------------------------------------------
    # evaluate the "hinge forces"-csv-files
    #-------------------------------------------------------
    #
    elif do == 'hf':

        # combination
        #
        lct = LCCTableULS( data_dir = data_dir,
                           lc_list = lc_list,
    #                       cut_z_fraction = 0.05,

                           # remove only the lowest point = connection shell/column
                           # as this is a singularity of the FE-shell-model
                           #
    #                       cut_z_fraction = 0.0000001,
                           show_lc_characteristic = True
                            )

        # INTERACTION plot

#        lct.plot_interaction_s6cm()
#        lct.plot_interaction_tp()

        # LCC-TABLE
#        lct.configure_traits()

        # Internal Force EXPORT
        #
        lct.link_type = link_type
        lct.export_hf_max_grouped( link_type + '_' + link_case + '_plt_' )
        lct.export_hf_lc()


#----------------------------------------
# run study for more link_cases
#----------------------------------------

#        spring_list = ['no_spring', 'spring']
#
#
#        link_case_list = [
#                          'equal_25cm',
#                          'equal_50cm',
#                          'equal_100cm',
#                          'equal_200cm',
#                          'corner_2to4_25cm',
#                          'gap_100cm',
#                          'middle_0to2_25cm',
#                          ]
#
#        link_type_list = [
#                          'exc_V_ip',
#                          'inc_V_ip',
#                          ]

#        spring_list = ['no_spring']
#        link_type_list = ['exc_V_ip']
##        link_case_list = ['equal_100cm_7m']
#        link_case_list = ['equal_88cm']
#
#        spring = spring_list[0]
##        link_type = link_type_list[0]
##        link_case = link_case_list[0]
##        spring_type = 'spring_1'
#
#        for link_case in link_case_list:
#            for link_type in link_type_list:
#
#                # file reading and export definition
#                #
#                data_dir = os.path.join( simdb.simdb_dir,
#                                         'simdata', 'output_data_mushroof',
#                                         'hinge_forces_and_displacement_2shells',
#                                         spring, link_type, link_case )
#
#                # combination
#                #
#                lct = LCCTableULS( data_dir = data_dir,
#                                   lc_list = lc_list,
#            #                       cut_z_fraction = 0.05,
#
#                                   # remove only the lowest point = connection shell/column
#                                   # as this is a singularity of the FE-shell-model
#                                   #
#            #                       cut_z_fraction = 0.0000001,
#                                   show_lc_characteristic = True
#                                    )
#
#
#                # INTERACTION plot
#
##                lct.plot_interaction_s6cm()
##                lct.plot_interaction_tp()
#
#                # LCC-TABLE
##                lct.configure_traits()
#
#
#                # Internal Force EXPORT
#                #
#                lct.link_type = link_type
#                lct.export_hf_max_grouped( link_type + '_' + link_case + '_plt_' )
#                lct.export_hf_lc()
#
#                # Print Maximum values
##                sr_list = ['N_ip', 'V_ip', 'V_op' ]
##
##                print link_case, link_type
##                for sr in sr_list:
##                    max_array = lct.max_sr_grouped_dict[sr]
##                    idx_sr = where( array( lct.sr_columns ) == sr )[0]
##                    if sr == 'N_ip':
##                        idx = argmax( max_array[:, idx_sr] )
##                    else:
##                       idx = argmax( abs( max_array[:, idx_sr] ) )
##                    print 'max', sr
##                    print max_array[idx]
