'''
Created on Jun 29, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, DelegatesTo, Callable, \
    Trait, Event

from etsproxy.util.home_directory import \
    get_home_directory

from etsproxy.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, \
    TableEditor, Group, ListEditor, VSplit, HSplit, VGroup, HGroup, Spring, \
    Include

from etsproxy.mayavi import \
    mlab
    
import pylab as p

from etsproxy.traits.ui.table_column import \
    ObjectColumn

from etsproxy.traits.ui.menu import \
    OKButton, CancelButton

from etsproxy.traits.ui.tabular_adapter \
    import TabularAdapter

from numpy import \
    array, loadtxt, ones_like, \
    vstack, hstack, \
    copy, where, sum, \
    ones, fabs, identity, \
    max as ndmax, min as ndmin

import numpy as np

import os.path

from ls_table import \
    ULS, SLS

from ls_table_Rxyz import \
    LSTableRxyz, ULSRxyz, SLSRxyz

from lcc_reader import LCCReader, LCCReaderRFEM, LCCReaderInfoCAD, LCCReaderInfoCADRxyz

class LC(HasTraits):
    '''Loading case class
    '''

    reader = WeakRef

    lcc_table = WeakRef

    # name of the file containing the stress resultants
    #
    file_name = Str(input = True)

    # data filter (used to hide unwanted values, e.g. high sigularities etc.) 
    #
    data_filter = Callable(input = True)

    # name of the loading case
    #
    name = Str(input = True)

    # category of the loading case
    #
    category = Enum('dead-load', 'additional dead-load', 'imposed-load', input = True)

    # list of keys specifying the names of the loading cases 
    # that can not exist at the same time, i.e. which are exclusive to each other
    # 
    exclusive_to = List(Str, input = True)
    def _exclusive_to_default(self):
        return []

    # combination factors (need to be defined in case of imposed loads)
    # 
    psi_0 = Float(input = True)
    psi_1 = Float(input = True)
    psi_2 = Float(input = True)

    # security factors ULS
    #
    gamma_fav = Float(input = True)
    def _gamma_fav_default(self):
        if self.category == 'dead-load':
            return 1.00
        if self.category == 'additional dead-load':
            return 0.00
        if self.category == 'imposed-load':
            return 0.00

    gamma_unf = Float(input = True)
    def _gamma_unf_default(self):
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
    gamma_fav_SLS = Float(input = True)
    def _gamma_fav_SLS_default(self):
        if self.category == 'dead-load':
            return 1.00
        elif self.category == 'additional dead-load' or \
             self.category == 'imposed-load':
            return 0.00

    gamma_unf_SLS = Float(input = True)
    def _gamma_unf_SLS_default(self):
        return 1.00

    # original state data (before filtering)
    #
    state_data_orig = Property(Dict, depends_on = 'file_name, lcc_table.data_dir')
    @cached_property
    def _get_state_data_orig(self):
        return self.reader.read_state_data(self.file_name)

    # state data (after filtering)
    #
    state_data_dict = Property(Dict, depends_on = 'file_name, data_filter, lcc_table.data_dir')
    @cached_property
    def _get_state_data_dict(self):
        d = {}
        for k, arr in self.state_data_orig.items():
            d[k] = self.data_filter(self.lcc_table, arr)
        return d

    # if reader == 'RFEM':
    #sr_columns = List(['mx', 'my', 'mxy', 'nx', 'ny', 'nxy'])
    # if reader == 'InfoCAD':
    #sr_columns = List(['mx', 'my', 'mxy', 'nx', 'ny', 'nxy', 'ux_elem', 'uy_elem', 'uz_elem'])
    # if reader == 'InfoCADRxyz':
    sr_columns = List(['Rx', 'Ry', 'Rz', 'Mx', 'My', 'Mz'])

    sr_arr = Property(Array)
    def _get_sr_arr(self):
        '''return the stress resultants of the loading case
        as stack of all sr-column arrays.
        '''
        sd_dict = self.state_data_dict
        return hstack([ sd_dict[ sr_key ] for sr_key in self.sr_columns ])


class LCC(HasTraits):

    lcc_id = Int

    #lcc_table = WeakRef()

#    ls_table = Instance(LSTable)
    # changes concerning 'Rxyz'
    ls_table = Instance(LSTableRxyz)
    
    assess_value = Property()
    def _get_assess_value(self):
        return self.ls_table.assess_value

    traits_view = View(Item('ls_table@', show_label = False),
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

class LCCTable(HasTraits):
    '''Loading Case Manager.
    Generates and sorts the loading case combinations
    of all specified loading cases.
    '''

    # define ls
    #
    ls = Trait('ULS',
                {'ULS' : ULS,
                 'SLS' : SLS })

    # manager of input files with state and geometry data 
    # 
    reader_type = Trait('RFEM', dict(RFEM = LCCReaderRFEM,
                                     InfoCAD = LCCReaderInfoCAD,
                                     # changes concerning 'Rxyz'
                                     InfoCADRxyz = LCCReaderInfoCADRxyz))

    reader = Property(Instance(LCCReader), depends_on = 'reader_type')
    @cached_property
    def _get_reader(self):
        return self.reader_type_(lcc_table = self)

    # lcc-instance for the view 
    #
    lcc = Instance(LCC)

    #-------------------------------
    # Define loading cases:
    #-------------------------------

    # path to the directory containing the state data files
    #
    data_dir = Directory

    # specify assess name of the parameter used to sort the rows in LCCTable in .configure_traits
    #
    assess_name = Str

    # list of load cases
    #
    lc_list_ = List(Instance(LC))
    lc_list = Property(List, depends_on = '+filter, data_dir')
    def _set_lc_list(self, value):
        self.lc_list_ = value
    def _get_lc_list(self):
        for lc in self.lc_list_:
            lc.reader = self.reader
            lc.lcc_table = self
            if lc.data_filter != self.data_filter:
                lc.data_filter = self.data_filter
        return self.lc_list_

    lcc_table_columns = Property(depends_on = 'lc_list_, +filter, data_dir')
    def _get_lcc_table_columns(self):
        return [ ObjectColumn(label = 'Id', name = 'lcc_id') ] + \
               [ ObjectColumn(label = lc.name, name = lc.name)
                for idx, lc in enumerate(self.lc_list) ] + \
                [ ObjectColumn(label = 'assess_value', name = 'assess_value') ]

    sr_columns = Property(List(Str), depends_on = 'lc_list_, +filter, data_dir')
    def _get_sr_columns(self):
        '''derive the order of the stress resultants
        from the first element in 'lc_list'. The internal
        consistency is checked separately in the
        'check_consistency' method.
        '''
        return self.lc_list[0].sr_columns

    #-------------------------------
    # check consistency
    #-------------------------------

    def _check_for_consistency(self):
        ''' check input files for consitency:
        '''
        pass
#        return self.reader.check_for_consistency(self.lc_list, self.geo_data_dict)

    #-------------------------------
    # lc_arr
    #-------------------------------

    lc_arr = Property(Array)
    def _get_lc_arr(self):
        '''stack stress resultants arrays of all loading cases together.
        This yields an array of shape ( n_lc, n_elems, n_sr )
        '''
        sr_arr_list = [ lc.sr_arr for lc in self.lc_list ]
        return array(sr_arr_list)

    #-------------------------------
    # Array dimensions:
    #-------------------------------

    n_sr = Property(Int)
    def _get_n_sr(self):
        return len(self.sr_columns)

    n_lc = Property(Int)
    def _get_n_lc(self):
        return len(self.lc_list)

    n_lcc = Property(Int)
    def _get_n_lcc(self):
        return self.combi_arr.shape[0]

    n_elems = Property(Int)
    def _get_n_elems(self):
        return self.lc_list[0].sr_arr.shape[0]

    #-------------------------------
    # auxilary method for get_combi_arr 
    #-------------------------------

    def _product(self, args):
        """
        Get all possible permutations of the security factors
        without changing the order of the loading cases.
        The method corresponds to the build-in function 'itertools.product'.
        Instead of returning a generator object a list of all
        possible permutations is returned. As argument a list of list
        needs to be defined. In the original version of 'itertools.product'
        the function takes a tuple as argument ("*args").
        """
        pools = map(tuple, args) #within original version args defined as *args
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
    imposed_idx_list = Property(List, depends_on = 'lc_list_')
    @cached_property
    def _get_imposed_idx_list(self):
        '''list of indices for the imposed loads
        '''
        imposed_idx_list = []
        for i_lc, lc in enumerate(self.lc_list):
            cat = lc.category
            if cat == 'imposed-load':
                imposed_idx_list.append(i_lc)
        return imposed_idx_list

    # array containing the psi with name 'psi_key' for the specified
    # loading cases defined in 'lc_list'. For dead-loads no value for
    # psi exists. In this case a value of 1.0 is defined. 
    # This yields an array of shape ( n_lc, )
    #
    def _get_psi_arr(self, psi_key):
        '''psi_key must be defined as:
        'psi_0', 'psi_1', or 'psi_2'
        Returns an 1d-array of shape ( n_lc, )
        '''
        # get list of ones (used for dead-loads):
        #
        psi_list = [1] * len(self.lc_list)

        # overwrite ones with psi-values in case of imposed-loads:
        #
        for imposed_idx in self.imposed_idx_list:
            psi_value = getattr(self.lc_list[ imposed_idx ], psi_key)
            psi_list[ imposed_idx ] = psi_value

        return array(psi_list, dtype = 'float_')

    # list containing names of the loading cases
    #
    lc_name_list = Property(List, depends_on = 'lc_list_')
    @cached_property
    def _get_lc_name_list(self):
        '''list of names of all loading cases
        '''
        return [ lc.name for lc in self.lc_list ]

    show_lc_characteristic = Bool(True)

    # combination array:
    #
    combi_arr = Property(Array, depends_on = 'lc_list_, combination_SLS')
    @cached_property
    def _get_combi_arr(self):
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
            print '*** SLS combination used: % s ***' % (self.combination_SLS)

        #---------------------------------------------------------------
        # get permutations of safety factors ('gamma')
        #---------------------------------------------------------------
        #
        permutation_list = self._product(self.gamma_list)

        combi_arr = array(permutation_list)

        # check if imposed loads are defined 
        # if not no further processing of 'combi_arr' is necessary:
        #
        if self.imposed_idx_list == []:

            # if option is set to 'True' the loading case combination table 
            # is enlarged with an identity matrix in order to see the 
            # characteristic values of each loading case.
            #
            if self.show_lc_characteristic:
                combi_arr = vstack([ identity(self.n_lc), combi_arr ])

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
        if all(psi_lead_arr == psi_non_lead_arr):
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
                psi_arr = copy(psi_non_lead_arr)
                psi_arr[imposed_idx] = psi_lead_arr[imposed_idx]
                combi_arr_lead_i = combi_arr[where(combi_arr[:, imposed_idx] != 0)] * psi_arr
                combi_arr_psi_list.append(combi_arr_lead_i)

            combi_arr_psi_no_0 = vstack(combi_arr_psi_list)

            # missing cases without any dead load have to be added 
            # get combinations with all!! imposed = 0 
            #
            lcc_all_imposed_zero = where((combi_arr[:, self.imposed_idx_list] == 0)
                                          .all(axis = 1))

            # add to combinations
            #
            combi_arr_psi = vstack((combi_arr[lcc_all_imposed_zero], combi_arr_psi_no_0))

        #---------------------------------------------------------------
        # get exclusive loading cases ('exclusive_to')
        #---------------------------------------------------------------

        # get a list of lists containing the indices of the loading cases
        # that are defined exclusive to each other. 
        # The list still contains duplicates, e.g. [1,2] and [2,1]
        #
        exclusive_list = []
        for i_lc, lc in enumerate(self.lc_list):

            # get related load case number
            #
            for exclusive_name in lc.exclusive_to:
                if exclusive_name in self.lc_name_list:
                    exclusive_idx = self.lc_name_list.index(exclusive_name)
                    exclusive_list.append([ i_lc, exclusive_idx ])

        # eliminate the duplicates in 'exclusive_list'
        #
        exclusive_list_unique = []
        for exclusive_list_entry in exclusive_list:
            if sorted(exclusive_list_entry) not in exclusive_list_unique:
                exclusive_list_unique.append(sorted(exclusive_list_entry))

        # delete the rows in combination array that contain
        # loading case combinations with imposed-loads that have been defined
        # as exclusive to each other. 
        # 
        combi_arr_psi_exclusive = combi_arr_psi
        for exclusive_list_entry in exclusive_list_unique:
            # check where maximum one value of the exclusive load cases is unequal to one 
            #              LC1  LC2  LC3  (all LCs are defined as exclusive to each other)
            #
            # e.g.         1.5  0.9  0.8  (example of 'combi_arr_psi') 
            #              1.5  0.0  0.0  
            #              0.0  0.0  0.0  (combination with all imposed loads = 0 after multiplication wit hpsi and gamma)
            #              ...  ...  ...
            #
            # this would yield the following mask_arr (containing ones or zeros):
            # e.g.         1.0  1.0  1.0  --> sum = 3 --> true combi --> accepted combination
            #              1.0  0.0  0.0  --> sum = 1 --> false combi --> no accepted combination
            # e.g.         0.0  0.0  0.0  --> sum = 0 --> true combi --> accepted combination (only body-loads)
            #              ...  ...  ...
            #
            mask_arr = where(combi_arr_psi_exclusive[ :, exclusive_list_entry ] != 0, 1.0, 0.0)
            true_combi = where(sum(mask_arr, axis = 1) <= 1.0)
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
            if (row == combi_arr_psi_exclusive_unique).all(axis = 1.0).any() == False:
                combi_arr_psi_exclusive_unique = vstack((combi_arr_psi_exclusive_unique, row))

        # if option is set to 'True' the loading case combination table 
        # is enlarged with an identity matrix in order to see the 
        # characteristic values of each loading case.
        #
        if self.show_lc_characteristic:
            combi_arr_psi_exclusive_unique = vstack([ identity(self.n_lc), combi_arr_psi_exclusive_unique ])

        return combi_arr_psi_exclusive_unique

    #-------------------------------
    # lcc_arr 
    #-------------------------------
    
    lcc_arr = Property(Array, depends_on = 'lc_list_, data_dir')
    @cached_property
    def _get_lcc_arr(self):
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
        lcc_arr = sum(lc_combi_arr, axis = 1)

        return lcc_arr

    lcc_arr_new = Array
    
    def _set_lcc_arr(self):
        print '_set_lcc_arr called!'
        self.lcc_arr = self.lcc_arr_new
#        return self.lcc_arr_new

#    lcc_changed = Event
#    @on_trait_change('lcc_arr_new')
#    def _set_lcc_changed(self):
#        setattr(self, 'lcc_arr', self.lcc_arr_new )
#        #self.lcc_arr = self.lcc_arr_new



    #-------------------------------
    # min/max-values - for verification only!
    #-------------------------------

    def get_min_max_state_data(self):
        ''' get the surrounding curve of all 'lcc' values
        '''
        lcc_arr = self.lcc_arr

        min_arr = ndmin(lcc_arr, axis = 0)
        max_arr = ndmax(lcc_arr, axis = 0)

        return min_arr, max_arr

    #------------------------------------------
    # read the geometry data from file 
    # (element number and thickness):
    #------------------------------------------

    # thickness_data input file:
    #
    geo_data_file = Property
    def _get_geo_data_file(self):
#        return os.path.join(self.data_dir, 'thickness.csv')
        # use the same file as used for state data as both
        # information is stored in the same file
        #
        pass 

    #------------------------------------------
    # get thickness data:
    #------------------------------------------

    # coordinates and element thickness read from file:
    # 
    geo_data_orig = Property(Dict, depends_on = 'geo_data_file')
    @cached_property
    def _get_geo_data_orig(self):
        # changes concerning 'Rxyz': get node numbers of the support from the same file as state data
        return self.lc_list[0].state_data_dict

    # parameter that defines for which z-coordinate values
    # the read in data is not evaluated (filtered out)
    # 
    cut_z_fraction = Float(0.0, filter = True)

    # construct a filter 
    #
    data_filter = Callable
    def _data_filter_default(self):
        return lambda lcc_table, x: x    #  - do nothing by default

    geo_data_dict = Property(Dict, depends_on = 'geo_data_file, +filter')
    @cached_property
    def _get_geo_data_dict(self):
        d = {}
        for k, arr in self.geo_data_orig.items():
            d[ k ] = self.data_filter(self, arr)
        return d

    #-------------------------------
    # lcc_lists 
    #-------------------------------

    lcc_list = Property(List, depends_on = 'lc_list_, data_dir')
    @cached_property
    def _get_lcc_list(self):
        '''list of loading case combinations (instances of LCC)
        '''
        combi_arr = self.combi_arr
        lcc_arr = self.lcc_arr
        sr_columns = self.sr_columns
        n_lcc = self.n_lcc

        # return a dictionary of the stress resultants
        # this is used by LSTable to determine the stress 
        # resultants of the current limit state 
        #
        lcc_list = []
        for i_lcc in range(n_lcc):

            state_data_dict = {}
            for i_sr, name in enumerate(sr_columns):
                state_data_dict[ name ] = lcc_arr[ i_lcc, :, i_sr ][:, None]

            lcc = LCC(#lcc_table = self,
                       factors = combi_arr[ i_lcc, : ],
                       lcc_id = i_lcc,
                       # changes concerning 'Rxyz': use LSTable for support force evaluation
                       ls_table = LSTableRxyz( reader = self.reader,
                                           reader_type = self.reader_type,
                                           assess_name = self.assess_name,
                                           geo_data = self.geo_data_dict,
                                           state_data = state_data_dict,
                                           ls = self.ls)
                       )

            for idx, lc in enumerate(self.lc_list):
            # @todo: use 'Float' instead of 'Int'    
                lcc.add_trait(lc.name, Int(combi_arr[ i_lcc, idx ]))

            lcc_list.append(lcc)

        return lcc_list

    def plot_geo(self, mlab):
        self.reader.plot_mesh(mlab, self.geo_data_dict)

    def plot_RxRz_interaction(self, save_fig_to_file = None, show_tension_only = False, add_max_min_RxRz_from_file = None, save_max_min_RxRz_to_file = None):
        '''plot the normal/shear-reaction forces as interaction diagram for all loading case combinations
        '''

        # get the list of all loading case combinations:
        #
        lcc_list = self.lcc_list

        #----------------------------------------------
        # run trough all loading case combinations:
        #----------------------------------------------

        Rx_Ed_list = [] # normal reaction force (tangential to the shell geometry)
        Rz_Ed_list = [] # shear reaction force (radial to the shell geometry)
        for lcc in lcc_list:

            # get the ls_table object and retrieve its 'ls_class'
            # (= LSTable_ULS-object)
            #
            ls_class = lcc.ls_table.ls_class
 
            # get Rx_Ed and Rz_Ed 
            #
            Rx_Ed = getattr(ls_class, 'Rx')
            Rz_Ed = getattr(ls_class, 'Rz')
            
            # add read in saved values to be superposed with currently read in values
            #
            if add_max_min_RxRz_from_file != None:
                max_min_RxRz_arr  = np.loadtxt( add_max_min_RxRz_from_file )
                max_Rx_Ed_arr = max_min_RxRz_arr[:,0][:,None]
                min_Rx_Ed_arr = max_min_RxRz_arr[:,1][:,None]
                max_Rz_Ed_arr = max_min_RxRz_arr[:,2][:,None]
                min_Rz_Ed_arr = max_min_RxRz_arr[:,3][:,None]
    
                # Rx_Ed
                #
                cond_Rx_Ed_ge_0 = Rx_Ed >= 0. # positive Rx value = shear force for the build-in srew 
                bool_arr = cond_Rx_Ed_ge_0
                Rx_Ed[bool_arr] += max_Rx_Ed_arr[bool_arr]

                cond_Rx_Ed_lt_0 = Rx_Ed < 0. # negative Rx value = compressive force in tangential direction at support (normal force) 
                bool_arr = cond_Rx_Ed_lt_0
                Rx_Ed[bool_arr] += min_Rx_Ed_arr[bool_arr]

                # Rz_Ed
                #
                cond_Rz_Ed_ge_0 = Rz_Ed >= 0. # positive Rz value = compressive force in radial direction at the support 
                bool_arr = cond_Rz_Ed_ge_0
                Rz_Ed[bool_arr] += max_Rz_Ed_arr[bool_arr]

                cond_Rz_Ed_lt_0 = Rz_Ed < 0. # negative Rz value = pull-out force for the build-in srew 
                bool_arr = cond_Rz_Ed_lt_0
                Rz_Ed[bool_arr] += min_Rz_Ed_arr[bool_arr]            
            
            Rx_Ed_list.append( Rx_Ed )
            Rz_Ed_list.append( Rz_Ed )

        # stack the list to an array in order to use plot-function
        #
        Rx_Ed_arr = hstack( Rx_Ed_list )
        Rz_Ed_arr = hstack( Rz_Ed_list )
        print 'Rz_Ed_arr.shape', Rz_Ed_arr.shape

        # get Rx_Rd, Rz_Rd 
        #
        Rx_Rd = getattr( ls_class, 'Rx_Rd' ) # shear resistance
        Rz_Rd = getattr( ls_class, 'Rz_Rd' ) # pull-out resistance

        # save min- and max-values to file in order to superpose them later
        #
        if save_max_min_RxRz_to_file != None:
                
            # get maximum or minimum values for superposition depending on
            # weather Rx or Rz are positive or negative values             
            
            # positive value for Rx
            #
            max_Rx_Ed_arr = np.max( Rx_Ed_arr, axis = 1) 

            # negative value for Rx
            #
            min_Rx_Ed_arr = np.min( Rx_Ed_arr, axis = 1) 

            # positive value for Rz
            #
            max_Rz_Ed_arr = np.max( Rz_Ed_arr, axis = 1) 

            # negative value for Rz
            #
            min_Rz_Ed_arr = np.min( Rz_Ed_arr, axis = 1) 

            # stack as four column array
            #
            max_min_RxRz_arr = np.hstack([max_Rx_Ed_arr[:,None], min_Rx_Ed_arr[:,None], max_Rz_Ed_arr[:,None], min_Rz_Ed_arr[:,None]])

            # save max and min values to file
            #
            np.savetxt( save_max_min_RxRz_to_file, max_min_RxRz_arr )
            print 'max_min_RxRz_arr saved to file %s' %(save_max_min_RxRz_to_file)

        #----------------------------------------------
        # plot
        #----------------------------------------------
        #
        p.figure(facecolor = 'white') # white background

        p.plot(Rx_Ed_arr, Rz_Ed_arr, 'wo', markersize=3) # blue dots
        print 'Rx_Ed_arr', Rx_Ed_arr
        print 'Rz_Ed_arr', Rz_Ed_arr
        x = np.array([0, Rx_Rd])
        y1 = np.array([ -Rz_Rd, 0. ])

#        p.title('$RxRz$-Interaktionsdiagramm')
    
        ax = p.gca()
        if show_tension_only == False:
#            ax.set_xticks([0., 0.2, 0.4, 0.6, 0.8, 1., 1.2])
#            ax.set_yticks([140., 120, 100, 80., 60., 40., 20., 0.])
#            p.axis([0., 1.05 * Rx_Rd, -1.2 * Rz_Rd, 0.]) # set plotting range for axis
            print 'show_tension_only == False'

        if show_tension_only == True:
#            ax.set_xticks([0., 0.2, 0.4, 0.6, 0.8, 1., 1.2])
#            ax.set_yticks([140., 120, 100, 80., 60., 40., 20., 0.])
            p.axis([0., 1.05 * Rx_Rd, -1.2 * Rz_Rd, 0.]) # set plotting range for axis
            print 'show_tension_only == True'

        ax.spines['left'].set_position(('data', 0))
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position(('data', 0))
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        p.plot(x,y1,'k--', linewidth=2.0) # black dashed line
        p.grid(True)

        p.xlabel('$R_{x,Ed}$ [kN]                                             ', fontsize=14, verticalalignment = 'top', horizontalalignment = 'right')
        p.ylabel('$R_{z,Ed}$ [kN]                                                              ', fontsize=14)

        # save figure as png-file
        #
        if save_fig_to_file != None:
            print 'figure saved to file %s' %(save_fig_to_file)
            p.savefig( save_fig_to_file, format='png' )
        
        p.show()
        

    #--------------------
    # eta_RxRz- interaction (new method with sign-consistend superposition of the values)
    #--------------------
    #
    # NOTE: old method superposed the eta_Rx and eta_Rz values too conservatively as inconsistend pairs are considered
    # that means that independently of the sign for Rx or Rz resulting from LC1-14 the most conservative eta of all temperature loading cases
    # is used. This corresponds to a translation of the obtained values for LC1-14 by a fixed horizontal and vertical shift.
    # a more realistic superposition needs to consider the signs of the (dominant) values obtained from LS1-14 and enlarge the eta-value only if the
    # resulting force of LC1-14 affects the in-build srew at all.   
    #
    def plot_eta_RxRz_interaction(self, show_tension_only = False, save_fig_to_file = None, add_max_min_RxRz_from_file = None):
        '''plot the eta_RxRz-interaction for all loading case combinations
        NOTE: the same superposition assumptions are made as used for method 'plot_RxRz_interaction'
        '''
        print 'plot_eta_RxRz_interaction; superposition based on max/min values of Rx and Rz'

        # get the list of all loading case combinations:
        #
        lcc_list = self.lcc_list

        #----------------------------------------------
        # run trough all loading case combinations:
        #----------------------------------------------
        
        eta_Rx_list = []
        eta_Rz_list = []

        for lcc in lcc_list:

            # get the ls_table object and retrieve its 'ls_class'
            # (= LSTable_ULS-object)
            #
            ls_class = lcc.ls_table.ls_class
 
            # get Rx_Ed and Rz_Ed 
            #
            Rx_Ed = getattr(ls_class, 'Rx')
            Rz_Ed = getattr(ls_class, 'Rz')

            # get Rx_Rd and Rz_Rd 
            #
            Rx_Rd = getattr(ls_class, 'Rx_Rd')
            Rz_Rd = getattr(ls_class, 'Rz_Rd')
            
            if add_max_min_RxRz_from_file == None:
                # get 'eta_Rx' and 'eta_Rz' 
                #
                eta_Rx = getattr(ls_class, 'eta_Rx')
                eta_Rz = getattr(ls_class, 'eta_Rz')
 
            # add read in saved values to be superposed with currently read in values
            #
            else:
                max_min_RxRz_arr  = np.loadtxt( add_max_min_RxRz_from_file )
                max_Rx_Ed_arr = max_min_RxRz_arr[:,0][:,None]
                min_Rx_Ed_arr = max_min_RxRz_arr[:,1][:,None]
                max_Rz_Ed_arr = max_min_RxRz_arr[:,2][:,None]
                min_Rz_Ed_arr = max_min_RxRz_arr[:,3][:,None]
    
                # Rx_Ed
                #
                cond_Rx_Ed_ge_0 = Rx_Ed >= 0. # positive Rx value = shear force for the build-in srew 
                bool_arr = cond_Rx_Ed_ge_0
                Rx_Ed[bool_arr] += max_Rx_Ed_arr[bool_arr]

                cond_Rx_Ed_lt_0 = Rx_Ed < 0. # negative Rx value = compressive force in tangential direction at support (normal force) 
                bool_arr = cond_Rx_Ed_lt_0
                Rx_Ed[bool_arr] += min_Rx_Ed_arr[bool_arr]

                # Rz_Ed
                #
                cond_Rz_Ed_ge_0 = Rz_Ed >= 0. # positive Rz value = compressive force in radial direction at the support 
                bool_arr = cond_Rz_Ed_ge_0
                Rz_Ed[bool_arr] += max_Rz_Ed_arr[bool_arr]

                cond_Rz_Ed_lt_0 = Rz_Ed < 0. # negative Rz value = pull-out force for the build-in srew 
                bool_arr = cond_Rz_Ed_lt_0
                Rz_Ed[bool_arr] += min_Rz_Ed_arr[bool_arr]         
                
                # note: positive values of 'Rx' correspond to shear forces for the support screw
                #       negative values are taken by the compression cushion at the support directly 
                #       Therefore take only the positive part of support force 'Rx' into account
                #       for the evaluation of 'eta_Rx'  
                Rx_pos = ( abs( Rx_Ed ) + Rx_Ed ) / 2.
        
                # eta shear forces 
                #
                eta_Rx = Rx_pos / Rx_Rd
        
                # note: negative values of 'Rz' correspond to pull-out forces for the support screw
                #       positive values are taken by the compression cushion at the support directly 
                #       Therefore take only the negative values of the support force 'Rz' into account
                #       for the evaluation of 'eta_Rz'  
                Rz_neg = ( abs( Rz_Ed ) - Rz_Ed ) / 2.
                
                # eta pull-out
                #
                eta_Rz = Rz_neg / Rz_Rd
            
            eta_Rx_list.append( eta_Rx )
            eta_Rz_list.append( eta_Rz )

        # stack the list to an array in order to use plot-function
        #
        eta_Rx_arr = hstack( eta_Rx_list )
        eta_Rz_arr = hstack( eta_Rz_list )
        print 'eta_Rx_arr.shape', eta_Rx_arr.shape

        #----------------------------------------------
        # plot
        #----------------------------------------------
        #
        p.figure(facecolor = 'white') # white background

        p.plot(eta_Rx_arr, eta_Rz_arr, 'wo', markersize=3) # blue dots
        x = np.array([0, 1. ])
        y1 = np.array([ -1., 0. ])
        y2 = np.array([  1., 0. ])
    
#        p.title('Ausnutzungsgrad $\eta$')

        ax = p.gca()
        if show_tension_only == False:
#            ax.set_xticks([0., 0.2, 0.4, 0.6, 0.8, 1., 1.2])
#            ax.set_yticks([140., 120, 100, 80., 60., 40., 20., 0.])
#            p.axis([0., 1.05 * Rx_Rd, -1.2 * Rz_Rd, 0.]) # set plotting range for axis
            print 'show_tension_only == False'

        if show_tension_only == True:
            ax.set_xticks([0., 0.2, 0.4, 0.6, 0.8, 1.])
            ax.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.])
            p.axis([0., 1., 1., 0.]) # set plotting range for axis

            print 'show_tension_only == True'

        p.plot(x,y1,'k--', linewidth=2.0) # black dashed line
        p.plot(x,y2,'k--', linewidth=2.0) # black dashed line
        p.grid(True)

        p.xlabel('$\eta_{Rx}$ [-] (Abscherkraft)', fontsize=14)
        p.ylabel('$\eta_{Rz}$ [-] (Auszugskraft)', fontsize=14)

        # save figure as png-file
        #
        if save_fig_to_file != None:
            print 'figure saved to file %s' %(save_fig_to_file)
            p.savefig( save_fig_to_file, format='png' )

        p.show()


    # ------------------------------------------------------------
    # View 
    # ------------------------------------------------------------

    traits_view = View(VGroup(

                        Item('geo_data_file',
                              label = 'Evaluated input file for geo_data ',
                               style = 'readonly', emphasized = True),
                        VSplit(
                                    Item('lcc_list', editor = lcc_list_editor,
                                          show_label = False),
                                    Item('lcc@', show_label = False),
                                    ),
                        ),
                      resizable = True,
                      scrollable = True,
                      height = 1.0,
                      width = 1.0
                      )

class LCCTableULS(LCCTable):
    '''LCCTable for ultimate limit state
    '''

    # set limit state to 'ULS'
    # (attribute is used by 'LSTable')
    #
    ls = 'ULS'

    # 'gamma' - safety factors
    #
    gamma_list = Property(List, depends_on = 'lc_list_')
    @cached_property
    def _get_gamma_list(self):
        return [[ lc.gamma_fav, lc.gamma_unf ] for lc in self.lc_list ]

    # 'psi' - combination factors (psi) for leading
    # and non leading load cases
    #
    psi_non_lead_arr = Property(Array, depends_on = 'lc_list_')
    @cached_property
    def _get_psi_non_lead_arr(self):
        return self._get_psi_arr('psi_0')

    psi_lead_arr = Property(Array, depends_on = 'lc_list_')
    @cached_property
    def _get_psi_lead_arr(self):
        return ones(len(self.lc_list))

class LCCTableSLS(LCCTable):
    '''LCCTable for serviceability limit state
    '''

    # set limit state to 'SLS'
    # (attribute is used by 'LSTable')
    #
    ls = 'SLS'

    # possible definitions of the serviceability limit state    
    # are: 'rare', 'freq', 'perm'
    #
    combination_SLS = Enum('rare', 'freq', 'perm')
    def _combination_SLS_default(self):
        return 'rare'

    # 'gamma' - safety factors
    #
    gamma_list = Property(List, depends_on = 'lc_list_')
    @cached_property
    def _get_gamma_list(self):

        # generate [1.0]-entry in case of body-loads:
        #
        gamma_list = [[ 1.0 ]] * len(self.lc_list)

        # overwrite those in case of imposed-loads:
        #
        for imposed_idx in self.imposed_idx_list:
            gamma_fav_SLS = getattr(self.lc_list[ imposed_idx ], 'gamma_fav_SLS')
            gamma_unf_SLS = getattr(self.lc_list[ imposed_idx ], 'gamma_unf_SLS')
            gamma_list[ imposed_idx ] = [ gamma_unf_SLS, gamma_fav_SLS ]

        return gamma_list

    # 'psi' - combination factors
    #
    psi_lead_dict = Property(Dict)
    def _get_psi_lead_dict(self):
        return {'rare' : ones_like(self._get_psi_arr('psi_0')) ,
                'freq' : self._get_psi_arr('psi_1'),
                'perm' : self._get_psi_arr('psi_2')}

    psi_non_lead_dict = Property(Dict)
    def _get_psi_non_lead_dict(self):
        return {'rare' : self._get_psi_arr('psi_0') ,
                'freq' : self._get_psi_arr('psi_2'),
                'perm' : self._get_psi_arr('psi_2')}

    # combination factors (psi) for leading
    # and non leading load cases
    #
    psi_lead_arr = Property(Array, depends_on = 'lc_list_, combination_SLS')
    @cached_property
    def _get_psi_lead_arr(self):
        return self.psi_lead_dict[ self.combination_SLS ]

    psi_non_lead_arr = Property(Array, depends_on = 'lc_list_, combination_SLS')
    @cached_property
    def _get_psi_non_lead_arr(self):
        return self.psi_non_lead_dict[ self.combination_SLS ]

