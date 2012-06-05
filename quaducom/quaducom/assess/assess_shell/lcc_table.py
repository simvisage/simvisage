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
    array, loadtxt, ones_like, \
    vstack, hstack, \
    copy, where, sum, \
    ones, fabs, identity, \
    max as ndmax, min as ndmin

import os.path

from ls_table import \
    LSTable, ULS, SLS

from matresdev.db.simdb import \
    SimDB

# Access to the top level directory of the database
#
simdb = SimDB()

class LC(HasTraits):
    '''Loading case class
    '''
    # name of the file containing the stress resultants
    #
    file_name = Str(input = True)

    # data filter (used to hide unwanted values, e.g. high sigularities etc.) 
    #
    data_filter = Callable(input = True)
    def _data_filter_default(self):
        return lambda x: x    #  - do nothing by default

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

    def _read_state_data(self, file_name):
        '''to read the stb-stress resultants save the xls-worksheet 
        to a csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''
        print '*** read state data from file: %s ***' % (file_name)

        # get the column headings defined in the second row 
        # of the csv soliciotations input file
        # column_headings = array(["Nr.","Punkt","X","Y","Z","mx","my","mxy","vx","vy","nx","ny","nxy"])
        #
        file = open(file_name, 'r')
        lines = file.readlines()
        column_headings = lines[1].split(';')
        # remove '\n' from last string element in list
        #
        column_headings[-1] = column_headings[-1][:-1]
        column_headings_arr = array(column_headings)

        elem_no_idx = where('Nr.' == column_headings_arr)[0]
        X_idx = where('X' == column_headings_arr)[0]
        Y_idx = where('Y' == column_headings_arr)[0]
        Z_idx = where('Z' == column_headings_arr)[0]
        mx_idx = where('mx' == column_headings_arr)[0]
        my_idx = where('my' == column_headings_arr)[0]
        mxy_idx = where('mxy' == column_headings_arr)[0]
        nx_idx = where('nx' == column_headings_arr)[0]
        ny_idx = where('ny' == column_headings_arr)[0]
        nxy_idx = where('nxy' == column_headings_arr)[0]

        file.close()

        # define arrays containing the information from the raw input file
        #
        input_arr = loadtxt(file_name , delimiter = ';', skiprows = 2)

        # element number:
        #
        elem_no = input_arr[:, elem_no_idx]

        # coordinates [m]:
        #
        X = input_arr[:, X_idx]
        Y = input_arr[:, Y_idx]
        Z = input_arr[:, Z_idx]

        # moments [kNm/m]
        #
        mx = input_arr[:, mx_idx]
        my = input_arr[:, my_idx]
        mxy = input_arr[:, mxy_idx]

        # normal forces [kN/m]:
        #
        nx = input_arr[:, nx_idx]
        ny = input_arr[:, ny_idx]
        nxy = input_arr[:, nxy_idx]

        return { 'elem_no' : elem_no, 'X' : X, 'Y' : Y, 'Z' : Z,
                 'mx' : mx, 'my' : my, 'mxy' : mxy,
                 'nx' : nx, 'ny' : ny, 'nxy' : nxy,
               }

    # original state data (before filtering)
    #
    state_data_orig = Property(Dict, depends_on = 'file_name')
    @cached_property
    def _get_state_data_orig(self):
        return self._read_state_data(self.file_name)

    # state data (after filtering)
    #
    state_data_dict = Property(Dict, depends_on = 'file_name, data_filter')
    @cached_property
    def _get_state_data_dict(self):
        d = {}
        for k, arr in self.state_data_orig.items():
            d[k] = self.data_filter(arr)
        return d

    sr_columns = List(['mx', 'my', 'mxy', 'nx', 'ny', 'nxy'])

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

    ls_table = Instance(LSTable)

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

    # lcc-instance for the view 
    #
    lcc = Instance(LCC)

    #-------------------------------
    # Define loading cases:
    #-------------------------------

    # path to the directory containing the state data files
    #
    data_dir = Directory

    # list of load cases
    #
    lc_list_ = List(Instance(LC))
    lc_list = Property(List, depends_on = '+filter')
    def _set_lc_list(self, value):
        self.lc_list_ = value
    def _get_lc_list(self):
        for lc in self.lc_list_:
            if lc.data_filter != self.data_filter:
                lc.data_filter = self.data_filter
        return self.lc_list_

    lcc_table_columns = Property(depends_on = 'lc_list_, +filter')
    def _get_lcc_table_columns(self):
        return [ ObjectColumn(label = 'Id', name = 'lcc_id') ] + \
               [ ObjectColumn(label = lc.name, name = lc.name)
                for idx, lc in enumerate(self.lc_list) ] + \
                [ ObjectColumn(label = 'assess_value', name = 'assess_value') ]

    sr_columns = Property(List(Str), depends_on = 'lc_list_, +filter')
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

        for lc in self.lc_list:
            # check internal LC-consitency: 
            # (compare coords-values of first LC with all other LC's in 'lc_list')
            #
            if not all(lc_list[0].state_data_dict['X'] == lc.state_data_dict['X']) and \
                not all(lc_list[0].state_data_dict['Y'] == lc.state_data_dict['Y']) and \
                not all(lc_list[0].state_data_dict['Z'] == lc.state_data_dict['Z']):
                raise ValueError, "coordinates in loading case '%s' and loading case '%s' are not identical. Check input files for consistency!" \
                        % (self.lc_list[0].name, lc.name)
                return False

            # check external consistency:
            # (compare 'elem_no' in 'thickness.csv' and in all loading-cases 
            # input files (e.g. 'LC1.csv') defined in 'lc_list')
            #
            if not all(self.geo_data_dict['elem_no'] == lc.state_data_dict['elem_no']):
                raise ValueError, "element numbers in loading case '%s' and loading case '%s' are not identical. Check input files for consistency!" \
                        % (self.lc_list[0].name, lc.name)
                return False

        print '*** input files checked for consistency (OK) ***'
        return True

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

    lcc_arr = Property(Array, depends_on = 'lc_list_')
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
        return os.path.join(self.data_dir, 'thickness.csv')

    def _read_geo_data(self, file_name):
        '''to read the stb - thickness save the xls - worksheet
        to a csv - file using ';' as filed delimiter and ' ' ( blank )
        as text delimiter.
        '''
        print '*** read geo data from file: %s ***' % (file_name)


        # coordinates [m]:
        # (NOTE: corrds are taken from the state data file of the first loading case) 

        # the column headings are defined in the first/second row 
        # of the csv thickness input file
        # Flaeche;;;Material;Dicke;;Exzentrizitaet;Integrierte Objekte;;;

        # (NOTE: the headings contain non-ascii characters. Therefore the
        #       column indices can not be derived using the 'where'-method.)

        # read the float data:
        #
        input_arr = loadtxt(file_name, usecols = (0, 5), delimiter = ';', skiprows = 2)
        elem_no_idx = 0
        thickness_idx = 1

        # element number:
        # (NOTE: column array must be of shape (n_elems, 1)
        #
        elem_no = input_arr[:, elem_no_idx][:, None]

        # element thickness [mm]:
        # (NOTE: column array must be of shape (n_elems, 1)
        #
        thickness = input_arr[:, thickness_idx][:, None]

        # coordinates [m]:
        # (NOTE: corrds are taken from the state data file of the first loading case) 
        #
        X = self.lc_list[0].state_data_orig['X']
        Y = self.lc_list[0].state_data_orig['Y']
        Z = self.lc_list[0].state_data_orig['Z']

        return  {'elem_no':elem_no,
                 'X':X, 'Y':Y, 'Z':Z,
                 'thickness':thickness }

    #------------------------------------------
    # get thickness data:
    #------------------------------------------

    # coordinates and element thickness read from file:
    # 
    geo_data_orig = Property(Dict, depends_on = 'geo_data_file')
    @cached_property
    def _get_geo_data_orig(self):
        return self._read_geo_data(self.geo_data_file)

    # parameter that defines for which z-coordinate values
    # the read in data is not evaluated (filtered out)
    # 
    cut_z_fraction = Float(0.0, filter = True)

    # construct a filter 
    #
    data_filter = Property(Callable, depends_on = 'geo_data_file, +filter')
    def _get_data_filter(self):

        def remove_midpoints(arr):
            # remove the center points of the shells 
            Zp = fabs(self.geo_data_orig['Z'])
            max_Z = max(Zp)
            min_Z = min(Zp)
            h = max_Z - min_Z
            z_active_idx = where(Zp >= (min_Z + self.cut_z_fraction * h))[0]
            return arr[ z_active_idx ]

        return remove_midpoints

    geo_data_dict = Property(Dict, depends_on = 'geo_data_file, +filter')
    @cached_property
    def _get_geo_data_dict(self):
        d = {}
        for k, arr in self.geo_data_orig.items():
            d[ k ] = self.data_filter(arr)
        return d

    #-------------------------------
    # lcc_lists 
    #-------------------------------

    lcc_list = Property(List, depends_on = 'lc_list_')
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
                       ls_table = LSTable(geo_data = self.geo_data_dict,
                                           state_data = state_data_dict,
                                           ls = self.ls)
                       )

            for idx, lc in enumerate(self.lc_list):
                lcc.add_trait(lc.name, Int(combi_arr[ i_lcc, idx ]))

            lcc_list.append(lcc)

        return lcc_list





    # ------------------------------------------------------------
    # View 
    # ------------------------------------------------------------

    traits_view = View(VGroup(

                        Item('geo_data_file',
                              label = 'Evaluated input file for thicknesses ',
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


if __name__ == '__main__':

    #---------------------------------------------
    # 2 shells: 
    # new geometry with new loading cases plus waterfilling
    #---------------------------------------------

    data_dir = os.path.join(simdb.simdb_dir,
                             'simdata', 'input_data_mushroof_stb',
                             'ZiE_state_data_2shells_delta_h_865mm_2011-08-16')


    #------------------------
    # define loading cases:
    #------------------------
    # NOTE:
    #
    #---------------------------------------------------------
    # "staendige und voruebergehende Bemessungssitauation":
    #---------------------------------------------------------
    lc_list = [
                 # LC1:
                 # own weight:
                 #
                 LC(name = 'g', category = 'dead-load',
                     file_name = os.path.join(data_dir, 'LC1.csv'),
                     ),

                 # LC2:
                 # s_sym:
                 #
                 LC(name = 's_sym', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC2.csv'),
                     exclusive_to = ['s_asym', 'WF'],
                     psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0
                     ),

                 # LC3:
                 # s_asym:
                 #
                 LC(name = 's_asym', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC3.csv'),
                     exclusive_to = ['s_sym', 'WF'],
                     psi_0 = 0.5, psi_1 = 0.2, psi_2 = 0.0
                     ),

                 # LC4:
                 # w_neg:
                 #
                 LC(name = 'w_neg', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC4.csv'),
                     exclusive_to = ['w_pos', 'w_asym', 'w_int', 'WF'],
                     psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0
                     ),

                 # LC5:
                 # w_pos:
                 #
                 LC(name = 'w_pos', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC5.csv'),
                     exclusive_to = ['w_neg', 'w_asym', 'w_int', 'WF'],
                     psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0
                     ),

                 # LC6:
                 # w_asym:
                 #
                 LC(name = 'w_asym', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC6.csv'),
                     exclusive_to = ['w_pos', 'w_neg', 'w_int', 'WF'],
                     psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0
                     ),

                 # LC7:
                 # w_int: Bauzustand
                 #
                 LC(name = 'w_int', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC7.csv'),
                     exclusive_to = ['w_pos', 'w_neg', 'w_asym', 'WF', 'T_shrinkage'],
                     psi_0 = 0.6, psi_1 = 0.2, psi_2 = 0.0
                     ),

                 # LC8: 
                 # 1 kN man load (corner):
                 #
                 LC(name = 'Q_corner', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC8.csv'),
                     exclusive_to = ['Q_edge', 'WF'], psi_0 = 0.0, psi_1 = 0.0, psi_2 = 0.0),

                 # LC9: 
                 # 1 kN man load (edge, center):
                 #
                 LC(name = 'Q_edge', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC9.csv'),
                     exclusive_to = ['Q_corner', 'WF'], psi_0 = 0.0, psi_1 = 0.0, psi_2 = 0.0),

                 # LC10:
                 # temperature (sommer):
                 #
                 LC(name = 'T_pos', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC10.csv'),
                     exclusive_to = ['T_neg', 'WF'], psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0),

                 # LC11:
                 # temperature (winter):
                 #
                 LC(name = 'T_neg', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC11.csv'),
                     exclusive_to = ['T_pos', 'WF'], psi_0 = 0.6, psi_1 = 0.5, psi_2 = 0.0),

                 # LC12:
                 # shrinkage:
                 # combination coefficients taken from case 'general imposed load'
                 # 
                 LC(name = 'T_shrinkage', category = 'imposed-load',
                     file_name = os.path.join(data_dir, 'LC12.csv'),
                     exclusive_to = ['WF'], psi_0 = 0.8, psi_1 = 0.7, psi_2 = 0.5),
               ]


    #---------------------------------------------------------
    # "aussergewoehnliche Bemessungssitauation":
    #---------------------------------------------------------
    #
#    lc_list = [
#                 # LC1:
#                 # own weight:
#                 #
#                 LC( name = 'g', category = 'dead-load',
#                     file_name = os.path.join( data_dir, 'LC1.csv' ),
#                     gamma_unf = 1.00
#                     ),
#
#                 # LC4:
#                 # w_neg:
#                 #
#                 LC( name = 'w_neg', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, 'LC4.csv' ),
#                     exclusive_to = ['w_pos', 'w_asym', ],
#                     gamma_unf = 0.2, # = psi_1
#                     psi_0 = 0.0      # = psi_2
#                     ),
#
#                 # LC5:
#                 # w_pos:
#                 #
#                 LC( name = 'w_pos', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, 'LC5.csv' ),
#                     exclusive_to = ['w_neg', 'w_asym'],
#                     gamma_unf = 0.2, # = psi_1
#                     psi_0 = 0.0      # = psi_2
#                     ),
#
#                 # LC6:
#                 # w_asym:
#                 #
#                 LC( name = 'w_asym', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, 'LC6.csv' ),
#                     exclusive_to = ['w_pos', 'w_neg'],
#                     gamma_unf = 0.2, # = psi_1
#                     psi_0 = 0.0      # = psi_2
#                     ),
#
#                 # LC10:
#                 # temperature (sommer):
#                 #
#                 LC( name = 'T_pos', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, 'LC10.csv' ),
#                     exclusive_to = ['T_neg'],
#                     gamma_unf = 0.6, # = psi_0
#                     psi_0 = 0.5 / 0.6 # = psi_1
#                     ),
#
#                 # LC11:
#                 # temperature (winter):
#                 #
#                 LC( name = 'T_neg', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, 'LC11.csv' ),
#                     exclusive_to = ['T_pos'],
#                     gamma_unf = 0.7, # = psi_1
#                     psi_0 = 0.5 / 0.7 # = psi_2
#                     ),
#
#                 # LC12:
#                 # shrinkage:
#                 # combination coefficients taken from case 'general imposed load'
#                 # 
#                 LC( name = 'T_shrinkage', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, 'LC12.csv' ),
#                     gamma_unf = 0.8, # = psi_0
#                     psi_0 = 0.7 / 0.8 # = psi_1
#                     ),
#
##                 # LC13:
##                 # water filling:
##                 #
##                 LC( name = 'WF', category = 'imposed-load',
##                     file_name = os.path.join( data_dir, 'LC13.csv' ),
##                     gamma_unf = 1.00,
##                     psi_0 = 1.0 ),
#
#                 # LC14_1-2:
#                 # earth quake:
#                 #
#                 LC( name = 'EQ_1-2', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, 'LC14_1-2.csv' ),
#                     gamma_unf = 1.00,
#                     exclusive_to = ['EQ_3-4'],
#                     psi_0 = 1.0 ),
#
#                 # LC14_3-4:
#                 # earth quake:
#                 #
#                 LC( name = 'EQ_3-4', category = 'imposed-load',
#                     file_name = os.path.join( data_dir, 'LC14_3-4.csv' ),
#                     gamma_unf = 1.00,
#                     exclusive_to = ['EQ_1-2'],
#                     psi_0 = 1.0 ),
#
#               ]

#--------------------------------------------------------


    lct = LCCTableULS(data_dir = data_dir,
                       lc_list = lc_list,

                       # remove only the lowest point = connection shell/column
                       # as this is a singularity of the FE-shell-model
                       #
#                       cut_z_fraction = 0.01,
#                       cut_z_fraction = 0.05, # corresponds to 50cm x 50cm
#                       cut_z_fraction = 0.10, # corresponds to 75cm x 75cm
#                       cut_z_fraction = 0.15, # corresponds to 100cm x 100cm
                       show_lc_characteristic = False
                        )

#    lct = LCCTableSLS( data_dir = data_dir,
#                       lc_list = lc_list,
#                       cut_z_fraction = 0.2,
#                       combination_SLS = 'rare',
##                       combination_SLS = 'freq',
##                       combination_SLS = 'perm',
##                       show_lc_characteristic = True
#                        )


#    lct.configure_traits()

    #----------------------------------------
    # script to get the maximum number of reinforcement ('n_tex')
    # at all given coordinate points for all possible loading
    # case combinations:
    #----------------------------------------
    # set option to "True" for surrounding 
    # evaluation of necessary layers "n_tex"
    # needed for all loading cases
    #
    plot_n_tex_all_lcc = True
#    plot_n_tex_all_lcc = False

    if plot_n_tex_all_lcc == True:

        print 'lct.combi_arr', lct.combi_arr

        # get the list of all loading case combinations:
        #
        lcc_list = lct.lcc_list

        #----------------------------------------------
        # run trough all loading case combinations:
        #----------------------------------------------

        n_tex_list = []
        for lcc in lcc_list:

            # get the ls_table object and retrieve its 'ls_class'
            # (= LSTable_ULS-object)
            #
            ls_class = lcc.ls_table.ls_class

            # get 'n_tex'-column array 
            #
            n_tex = ls_class.n_tex
#            n_tex = ls_class.n_tex_up
#            n_tex = ls_class.n_tex_lo
            n_tex_list.append(n_tex)


        # stack the list to an array in order to use ndmax-function
        #
        n_tex_arr = hstack(n_tex_list)

        #----------------------------------------------
        # get the overall maximum values:
        #----------------------------------------------

        n_tex_max = ndmax(n_tex_arr, axis = 1)[:, None]

        #----------------------------------------------
        # plot
        #----------------------------------------------
        #
        X = lcc_list[0].ls_table.X[:, 0]
        Y = lcc_list[0].ls_table.Y[:, 0]
        Z = lcc_list[0].ls_table.Z[:, 0]
        plot_col = n_tex_max[:, 0]

        # if n_tex is negative plot 0 instead:
        #
        plot_col = where(plot_col < 0, 0, plot_col)

        mlab.figure(figure = "SFB532Demo",
                     bgcolor = (1.0, 1.0, 1.0),
                     fgcolor = (0.0, 0.0, 0.0))

        mlab.points3d(X, Y, (-1.0) * Z, plot_col,
                       colormap = "YlOrBr",
                       mode = "cube",
                       scale_factor = 0.10)

        mlab.scalarbar(title = 'n_tex (all LCs)', orientation = 'vertical')

        mlab.show()



#    print 'lc_arr', lct.lc_arr
#    print 'lc_list[0].sr_arr.shape[0]', lct.lc_list[0].sr_arr.shape[0]
#    print 'lc_arr.shape', lct.lc_arr.shape
#    print 'combi_arr', lct.combi_arr
#    print 'combi_arr.shape', lct.combi_arr.shape
#    print 'lcc_arr', lct.lcc_arr
