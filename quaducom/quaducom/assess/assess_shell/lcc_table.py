'''
Created on Jun 29, 2010

@author: alexander
'''

from etsproxy.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, Enum, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Color, Bool, DelegatesTo, Callable, \
    Trait

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

import numpy as np

from numpy import \
    array, loadtxt, ones_like, \
    vstack, hstack, \
    copy, where, sum, \
    ones, fabs, identity, \
    max as ndmax, min as ndmin

import os.path

from matresdev.db.simdb import \
    SimDB

import os

from ls_table import \
    LSTable, ULS, SLS

from lcc_reader import LCCReader, LCCReaderRFEM, LCCReaderInfoCAD

simdb = SimDB()

class LC(HasTraits):
    '''Loading case class
    '''

    reader = WeakRef

    lcc_table = WeakRef

    # name of the file containing the stress resultants
    #
    file_name = Str(input=True)

    # data filter (used to hide unwanted values, e.g. high sigularities etc.)
    #
    data_filter = Callable(input=True)

    # name of the loading case
    #
    name = Str(input=True)

    # category of the loading case
    #
    category = Enum('dead-load', 'additional dead-load', 'imposed-load', input=True)

    # list of keys specifying the names of the loading cases
    # that can not exist at the same time, i.e. which are exclusive to each other
    #
    exclusive_to = List(Str, input=True)
    def _exclusive_to_default(self):
        return []

    # combination factors (need to be defined in case of imposed loads)
    #
    psi_0 = Float(input=True)
    psi_1 = Float(input=True)
    psi_2 = Float(input=True)

    # security factors ULS
    #
    gamma_fav = Float(input=True)
    def _gamma_fav_default(self):
        if self.category == 'dead-load':
            return 1.00
        if self.category == 'additional dead-load':
            return 0.00
        if self.category == 'imposed-load':
            return 0.00

    gamma_unf = Float(input=True)
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
    gamma_fav_SLS = Float(input=True)
    def _gamma_fav_SLS_default(self):
        if self.category == 'dead-load':
            return 1.00
        elif self.category == 'additional dead-load' or \
             self.category == 'imposed-load':
            return 0.00

    gamma_unf_SLS = Float(input=True)
    def _gamma_unf_SLS_default(self):
        return 1.00

    # original state data (before filtering)
    #
    state_data_orig = Property(Dict, depends_on='file_name')
    @cached_property
    def _get_state_data_orig(self):
        return self.reader.read_state_data(self.file_name)

    # state data (after filtering)
    #
    state_data_dict = Property(Dict, depends_on='file_name, +filter')
    @cached_property
    def _get_state_data_dict(self):
        d = {}
        for k, arr in self.state_data_orig.items():
            d[ k ] = self.data_filter(self.lcc_table, arr)
        return d

    sr_columns = Property(List)
    def _get_sr_columns(self):
        '''return the list of the stress resultants to be use within the combinations of LCC.
        '''
        # if LCCTableReaderInfoCAD is used:
        # use displacement stored in 'state_data' within 'plot_col' method of the reader
        # sr_columns = List(['mx', 'my', 'mxy', 'nx', 'ny', 'nxy', 'ux_elem', 'uy_elem', 'uz_elem'])
        # if LCCTableReaderRFEM is used:
        # no displacement is available yet
        # sr_columns = List(['mx', 'my', 'mxy', 'nx', 'ny', 'nxy'])
        return self.reader.sr_columns

    sr_arr = Property(Array)
    def _get_sr_arr(self):
        '''return the stress resultants of the loading case
        as stack of all sr-column arrays.
        '''
        sd_dict = self.state_data_dict
        return hstack([ sd_dict[ sr_key ] for sr_key in self.sr_columns ])

#    # deformation data
#    #
#    u_data_dict = Property(Dict)
#    @cached_property
#    def _get_u_data_dict(self):
#        return self.reader.read_u_data(self.file_name)
#
#    u_arr = Property(Array)
#    def _get_u_arr(self):
#        '''return the element deformation of the loading case
#        as stack of all u-column arrays.
#        '''
#        u_dict = self.u_data_dict
#        return hstack([ u_dict[ u_key ] for u_key in u_dict.keys() ])

class LCC(HasTraits):

    lcc_id = Int

    # lcc_table = WeakRef()

    ls_table = Instance(LSTable)

    assess_value = Property()
    def _get_assess_value(self):
        return self.ls_table.assess_value

    traits_view = View(Item('ls_table@', show_label=False),
                        resizable=True,
                        scrollable=True
                        )

# The definition of the demo TableEditor:
lcc_list_editor = TableEditor(
    columns_name='lcc_table_columns',
    editable=False,
    selection_mode='row',
    selected='object.lcc',
    show_toolbar=True,
    auto_add=False,
    configurable=True,
    sortable=True,
    reorderable=False,
    sort_model=False,
    auto_size=False,
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
    reader_type = Trait('RFEM', dict(RFEM=LCCReaderRFEM,
                                     InfoCAD=LCCReaderInfoCAD))

    reader = Property(Instance(LCCReader), depends_on='reader_type')
    @cached_property
    def _get_reader(self):
        return self.reader_type_(lcc_table=self)

    # lcc-instance for the view
    #
    lcc = Instance(LCC)

    #-------------------------------
    # Define loading cases:
    #-------------------------------

    # path to the directory containing the state data files
    #
    data_dir = Directory

    # specify weather exact evaluation 'k_alpha' is to be used or a lower bound 'k_min_alpha = 0.707' as simplification instead
    #
    k_alpha_min = Bool(False)

    # specify the strength characteristics of the material;
    # necessary parameter in order to pass values from 'LCCTableULS' to 'ls_table' and as WeekRef to 'ULS';
    # this gives the possibility to specify strength values within the call of 'LCCTableULS';
    #
    strength_characteristics = Dict

    # list of load cases
    #
    lc_list_ = List(Instance(LC))
    lc_list = Property(List, depends_on='+filter')
    def _set_lc_list(self, value):
        self.lc_list_ = value
    def _get_lc_list(self):
        for lc in self.lc_list_:
            lc.reader = self.reader
            lc.lcc_table = self
            if lc.data_filter != self.data_filter:
                lc.data_filter = self.data_filter
        return self.lc_list_

    lcc_table_columns = Property(depends_on='lc_list_, +filter')
    @cached_property
    def _get_lcc_table_columns(self):
        return [ ObjectColumn(label='Id', name='lcc_id') ] + \
               [ ObjectColumn(label=lc.name, name=lc.name)
                for idx, lc in enumerate(self.lc_list) ] + \
                [ ObjectColumn(label='assess_value', name='assess_value') ]

    sr_columns = Property(List(Str), depends_on='lc_list_, +filter')
    @cached_property
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

        return self.reader.check_for_consistency(self.lc_list, self.geo_data_dict)

    #-------------------------------
    # lc_arr
    #-------------------------------

    lc_arr = Property(Array)
    @cached_property
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
        pools = map(tuple, args)  # within original version args defined as *args
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
    imposed_idx_list = Property(List, depends_on='lc_list_')
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

        return array(psi_list, dtype='float_')

    # list containing names of the loading cases
    #
    lc_name_list = Property(List, depends_on='lc_list_')
    @cached_property
    def _get_lc_name_list(self):
        '''list of names of all loading cases
        '''
        return [ lc.name for lc in self.lc_list ]

    show_lc_characteristic = Bool(True)

    # combination array:
    #
    combi_arr = Property(Array, depends_on='lc_list_, combination_SLS')
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
                                          .all(axis=1))

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
            true_combi = where(sum(mask_arr, axis=1) <= 1.0)
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
            if (row == combi_arr_psi_exclusive_unique).all(axis=1.0).any() == False:
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

    lcc_arr = Property(Array, depends_on='lc_list_')
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
        lcc_arr = sum(lc_combi_arr, axis=1)

        return lcc_arr

    #-------------------------------
    # min/max-values - for verification only!
    #-------------------------------

    def get_min_max_state_data(self):
        ''' get the surrounding curve of all 'lcc' values
        '''
        lcc_arr = self.lcc_arr

        min_arr = ndmin(lcc_arr, axis=0)
        max_arr = ndmax(lcc_arr, axis=0)

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

    #------------------------------------------
    # get thickness data:
    #------------------------------------------

    # coordinates and element thickness read from file:
    #
    geo_data_orig = Property(Dict, depends_on='geo_data_file')
    @cached_property
    def _get_geo_data_orig(self):
        return self.reader.read_geo_data(self.geo_data_file)

    # parameter that defines for which z-coordinate values
    # the read in data is not evaluated (filtered out)
    #
    cut_z_fraction = Float(0.0, filter=True)

    # construct a filter
    #
    data_filter = Callable(filter=True)
    def _data_filter_default(self):
        return lambda lcc_table, x: x  #  - do nothing by default

    geo_data_dict = Property(Dict, depends_on='geo_data_file, +filter')
    @cached_property
    def _get_geo_data_dict(self):
        d = {}
        for k, arr in self.geo_data_orig.items():
            d[ k ] = self.data_filter(self, arr)
        return d

    #-------------------------------
    # lcc_lists
    #-------------------------------

    lcc_list = Property(List, depends_on='lc_list_')
    @cached_property
    def _get_lcc_list(self):
        '''list of loading case combinations (instances of LCC)
        '''
        print '################### constructing lcc_list ##################'

        combi_arr = self.combi_arr
        lcc_arr = self.lcc_arr
        sr_columns = self.sr_columns
        n_lcc = self.n_lcc

        print '################### arrays accessed ##################'

        # return a dictionary of the stress resultants
        # this is used by LSTable to determine the stress
        # resultants of the current limit state
        #
        lcc_list = []
        for i_lcc in range(n_lcc):

            state_data_dict = {}
            for i_sr, name in enumerate(sr_columns):
                state_data_dict[ name ] = lcc_arr[ i_lcc, :, i_sr ][:, None]

            print '################### lcc %s ##################' % str(i_lcc)

            lcc = LCC(# lcc_table = self,
                       factors=combi_arr[ i_lcc, : ],
                       lcc_id=i_lcc,
                       # changes necessary to distinguish plotting functionality defined in 'ls_table'
                       # e.i. plot deformation when 'LCCReaderInfoCAD' is used
                       ls_table=LSTable(reader=self.reader,
                                        geo_data=self.geo_data_dict,
                                        state_data=state_data_dict,
                                        strength_characteristics=self.strength_characteristics,
                                        k_alpha_min=self.k_alpha_min,
                                        ls=self.ls)
                       )

            for idx, lc in enumerate(self.lc_list):
            # @todo: use 'Float' instead of 'Int'
                lcc.add_trait(lc.name, Int(combi_arr[ i_lcc, idx ]))

            lcc_list.append(lcc)

        print '################### returning lcc_list ##################'

        return lcc_list

    def plot_geo(self, mlab):
        self.reader.plot_mesh(mlab, self.geo_data_dict)

    def plot_sr(self, mlab, lc=0, sr=0):
        lc = self.lc_arr[lc]
        gd = self.geo_data_dict
        mlab.points3d(gd['X'], gd['Y'], gd['Z'], lc[:, sr],
                       # colormap = "YlOrBr",
                       mode="cube",
                       scale_factor=0.1)

    def plot_assess_value(self, assess_name, scale_factor=0.1, scale_mode='none', azimuth=None, elevation=None, distance=None, focalpoint=None, title=None, save_fig_to_file=None,
                          add_assess_values_from_file=None, save_assess_values_to_file=None):
        '''plot-3d the assess value for all loading case combinations as structure plot with color legend
           options: - plot options for mlab
                    - save figure with specified name within "/simdb/simdata/lcc_table/output_images/save_fig_to_file.png"
                    - superpose data values, i.e. from imposed loads and temperature loads
        '''
        print '################### plotting assess_value ##################'

        #----------------------------------------
        # script to get the maximum values of 'assess_value'
        # at all given coordinate points for all possible loading
        # case combinations:
        #----------------------------------------

        # get the list of all loading case combinations:
        #
        lcc_list = self.lcc_list

        #----------------------------------------------
        # run trough all loading case combinations:
        #----------------------------------------------

        assess_value_list = []
        for lcc in lcc_list:

            # get the ls_table object and retrieve its 'ls_class'
            # (= LSTable_ULS-object)
            #
            ls_class = lcc.ls_table.ls_class

            # get 'assess_name'-column array
            assess_value = getattr(ls_class, assess_name)

            assess_value_list.append(assess_value)

        # stack the list to an array in order to use ndmax-function
        #
        assess_value_arr = hstack(assess_value_list)
        print 'assess_value_arr.shape', assess_value_arr.shape

        #----------------------------------------------
        # get the overall maximum values:
        #----------------------------------------------

        assess_value_max = ndmax(assess_value_arr, axis=1)[:, None]

        #----------------------------------------------
        # plot
        #----------------------------------------------
        #
        X = lcc_list[0].ls_table.X[:, 0]
        Y = lcc_list[0].ls_table.Y[:, 0]
        Z = lcc_list[0].ls_table.Z[:, 0]
        plot_col = assess_value_max[:, 0]

        # save assess values to file in order to superpose them later
        #
        if save_assess_values_to_file != None:
            simdata_dir = os.path.join(simdb.simdata_dir, 'lcc_table')
            filename = os.path.join(simdata_dir, save_assess_values_to_file)
            assess_value_arr = plot_col
            np.savetxt(filename, assess_value_arr)
            print 'assess_values saved to file %s' % (filename)

        # add read in saved assess values to be superposed with currently read in assess values
        #
        if add_assess_values_from_file != None:
            simdata_dir = os.path.join(simdb.simdata_dir, 'lcc_table')
            filename = os.path.join(simdata_dir, add_assess_values_from_file)
            assess_value_arr = np.loadtxt(filename)
            plot_col += assess_value_arr
            print 'superpose assess_value_arr with values read in from file %s' % (filename)

#        # if n_tex is negative plot 0 instead:
#        #
#        plot_col = where(plot_col < 0, 0, plot_col)

        mlab.figure(figure=title,
                     bgcolor=(1.0, 1.0, 1.0),
                     fgcolor=(0.0, 0.0, 0.0))

        mlab.points3d(X, Y, Z, plot_col,
                       colormap="YlOrBr",
                       mode="cube",
                       scale_mode=scale_mode,
                       scale_factor=scale_factor)

        mlab.scalarbar(title=assess_name + ' (all LCs)', orientation='vertical')

        mlab.view(azimuth=azimuth, elevation=elevation, distance=distance, focalpoint=focalpoint)

        if save_fig_to_file != None:
            simdata_dir = os.path.join(simdb.simdata_dir, 'lcc_table')
            img_dir = os.path.join(simdata_dir, 'output_images')
            filename = os.path.join(img_dir, save_fig_to_file + '.png')
            mlab.savefig(filename)  # , format='png')
            print 'figure saved to file %s' % (filename)

        mlab.show()

    def plot_nm_interaction(self, save_fig_to_file=None, show_tension_only=False, add_max_min_nm_from_file=None, save_max_min_nm_to_file=None):
        '''plot the nm-interaction for all loading case combinations
        '''

        # get the list of all loading case combinations:
        #
        lcc_list = self.lcc_list

        #----------------------------------------------
        # run trough all loading case combinations:
        #----------------------------------------------

        m_sig_lo_list = []
        n_sig_lo_list = []
        m_sig_up_list = []
        n_sig_up_list = []
        m_Ed_list = []
        n_Ed_list = []
        for lcc in lcc_list:

            # get the ls_table object and retrieve its 'ls_class'
            # (= LSTable_ULS-object)
            #
            ls_class = lcc.ls_table.ls_class

            # get n_Ed and m_Ed
            #
            m_sig_lo = np.copy(getattr(ls_class, 'm_sig_lo'))
            n_sig_lo = np.copy(getattr(ls_class, 'n_sig_lo'))
            m_sig_up = np.copy(getattr(ls_class, 'm_sig_up'))
            n_sig_up = np.copy(getattr(ls_class, 'n_sig_up'))

            # add read in saved values to be superposed with currently read in values
            #
            if add_max_min_nm_from_file != None:
                simdata_dir = os.path.join(simdb.simdata_dir, 'lcc_table')
                filename = os.path.join(simdata_dir, add_max_min_nm_from_file)
                max_min_nm_arr = np.loadtxt(filename)
                max_n_arr = max_min_nm_arr[:, 0][:, None]
                min_n_arr = max_min_nm_arr[:, 1][:, None]
                max_m_arr = max_min_nm_arr[:, 2][:, None]
                min_m_arr = max_min_nm_arr[:, 3][:, None]

                # n_sig_lo
                #
                cond_n_sig_lo_ge_0 = n_sig_lo >= 0.  # tensile normal force
                bool_arr = cond_n_sig_lo_ge_0
                n_sig_lo[bool_arr] += max_n_arr[bool_arr]

                cond_n_sig_lo_lt_0 = n_sig_lo < 0.  # compressive normal force
                bool_arr = cond_n_sig_lo_lt_0
                n_sig_lo[bool_arr] += min_n_arr[bool_arr]

                # n_sig_up
                #
                cond_n_sig_up_ge_0 = n_sig_up >= 0.  # tensile normal force
                bool_arr = cond_n_sig_up_ge_0
                n_sig_up[bool_arr] += max_n_arr[bool_arr]

                cond_n_sig_up_lt_0 = n_sig_up < 0.  # compressive normal force
                bool_arr = cond_n_sig_up_lt_0
                n_sig_up[bool_arr] += min_n_arr[bool_arr]

                # m_sig_lo
                #
                cond_m_sig_lo_ge_0 = m_sig_lo >= 0.  # positive bending moment
                bool_arr = cond_m_sig_lo_ge_0
                m_sig_lo[bool_arr] += max_m_arr[bool_arr]

                cond_m_sig_lo_lt_0 = m_sig_lo < 0.  # compressive normal force
                bool_arr = cond_m_sig_lo_lt_0
                m_sig_lo[bool_arr] += min_m_arr[bool_arr]

                # m_sig_up
                #
                cond_m_sig_up_ge_0 = m_sig_up >= 0.  # positive bending moment
                bool_arr = cond_m_sig_up_ge_0
                m_sig_up[bool_arr] += max_m_arr[bool_arr]

                cond_m_sig_up_lt_0 = m_sig_up < 0.  # compressive normal force
                bool_arr = cond_m_sig_up_lt_0
                m_sig_up[bool_arr] += min_m_arr[bool_arr]

            # NOTE: after the superposition only the absolute values of m_Ed are used for the plot
            #
            m_sig_lo_list.append(abs(m_sig_lo))
            m_sig_up_list.append(abs(m_sig_up))
            n_sig_lo_list.append(n_sig_lo)
            n_sig_up_list.append(n_sig_up)

            m_Ed_list = m_sig_lo_list + m_sig_up_list
            n_Ed_list = n_sig_lo_list + n_sig_up_list

        # stack the list to an array in order to use plot-function
        #
        m_Ed_arr = hstack(m_Ed_list)
        n_Ed_arr = hstack(n_Ed_list)
        print 'm_Ed_arr.shape', m_Ed_arr.shape

        # get n_tRd, n_cRd, m_Rd
        #
        m_0_Rd = getattr(ls_class, 'm_0_Rd')
        m_90_Rd = getattr(ls_class, 'm_90_Rd')
        n_0_Rdt = getattr(ls_class, 'n_0_Rdt')
        n_90_Rdt = getattr(ls_class, 'n_90_Rdt')
        n_Rdc = -getattr(ls_class, 'n_Rdc')

        # use simplification with minimum value for k_alpha = 0.707
        # and lower resistance of 0- or 90-degree direction
        #
        print 'simplification with k_alpha,min = 0.707 has been used for plot'
        m_Rd = min(m_0_Rd, m_90_Rd) * 0.707
        n_Rdt = min(n_0_Rdt, n_90_Rdt) * 0.707

        # save min- and max-values to file in order to superpose them later
        #
        if save_max_min_nm_to_file != None:

            # get maximum values for superposition if n_Ed is a positive number or minimum if it is a negative number

            # (positive) tensile force
            #
            max_n_arr = np.max(n_Ed_arr, axis=1)

            # (negative) compression force
            #
            min_n_arr = np.min(n_Ed_arr, axis=1)

            # positive bending moment
            #
            max_m_arr = np.max(m_Ed_arr, axis=1)

            # negative bending moment
            #
            min_m_arr = np.min(m_Ed_arr, axis=1)

            # stack as three column array
            #
            max_min_nm_arr = np.hstack([max_n_arr[:, None], min_n_arr[:, None], max_m_arr[:, None], min_m_arr[:, None]])

            # save max and min values to file
            #
            simdata_dir = os.path.join(simdb.simdata_dir, 'lcc_table')
            filename = os.path.join(simdata_dir, save_max_min_nm_to_file)
            np.savetxt(filename, max_min_nm_arr)
            print 'max_min_nm_arr saved to file %s' % (filename)

        #----------------------------------------------
        # plot
        #----------------------------------------------
        #
        p.figure(facecolor='white')  # white background

        p.plot(m_Ed_arr, n_Ed_arr, 'wo', markersize=3)  # blue dots
        x = np.array([0, m_Rd])
        y1 = np.array([ n_Rdc, 0. ])
        y2 = np.array([ n_Rdt, 0. ])

#        p.title('$nm$-Interaktionsdiagramm')

#        ax = p.gca()
#        if show_tension_only == False:
#            ax.set_xticks([0., 0.2, 0.4, 0.6, 0.8, 1., 1.2])
#            ax.set_yticks([200., 0., -200, -400, -600, -800])
#            p.axis([0., 1.05 * m_Rd, 1.2 * n_Rdt, 1.03 * n_Rdc])  # set plotting range for axis
#
#        if show_tension_only == True:
#            ax.set_xticks([0., 0.2, 0.4, 0.6, 0.8, 1., 1.2])
#            ax.set_yticks([140., 120, 100, 80., 60., 40., 20., 0.])
#            p.axis([0., 1.2, 140., 0.])  # set plotting range for axis

        p.plot(x, y1, 'k--', linewidth=2.0)  # black dashed line
        p.plot(x, y2, 'k--', linewidth=2.0)  # black dashed line
        p.grid(True)

#        ax.spines['left'].set_position(('data', 0))
#        ax.spines['right'].set_color('none')
#        ax.spines['bottom'].set_position(('data', 0))
#        ax.spines['top'].set_color('none')
#        ax.xaxis.set_ticks_position('bottom')
#        ax.yaxis.set_ticks_position('left')
        p.xlabel('$m_{Ed}$ [kNm/m]', fontsize=14, verticalalignment='top', horizontalalignment='right')
        p.ylabel('$n_{Ed}$ [kN/m]', fontsize=14)

        # save figure as png-file
        #
        if save_fig_to_file != None:
            simdata_dir = os.path.join(simdb.simdata_dir, 'lcc_table')
            img_dir = os.path.join(simdata_dir, 'output_images')
            filename = os.path.join(img_dir, save_fig_to_file)
            print 'figure saved to file %s' % (filename)
            p.savefig(filename, format='png')

        p.show()


    def plot_eta_nm_interaction(self, save_fig_to_file=None, show_tension_only=False, save_max_min_eta_nm_to_file=None, add_max_min_eta_nm_from_file=None):
        '''plot the eta_nm-interaction for all loading case combinations
        '''

        # get the list of all loading case combinations:
        #
        lcc_list = self.lcc_list

        #----------------------------------------------
        # run trough all loading case combinations:
        #----------------------------------------------

        # eta for 1st principal direction:
        #
        eta_n_lo_list = []
        eta_m_lo_list = []
        eta_n_up_list = []
        eta_m_up_list = []

        # eta for 2nd principal direction:
        #
        eta_n2_lo_list = []
        eta_m2_lo_list = []
        eta_n2_up_list = []
        eta_m2_up_list = []

        # get envelope over all loading case combinations (lcc)
        #
        for lcc in lcc_list:

            # get the ls_table object and retrieve its 'ls_class'
            # (= LSTable_ULS-object)
            #
            ls_class = lcc.ls_table.ls_class

            # get 'eta_n' and 'eta_m' (1-principle direction)
            #
            eta_m_lo = np.copy(getattr(ls_class, 'eta_m_lo'))
            eta_m_up = np.copy(getattr(ls_class, 'eta_m_up'))
            eta_n_lo = np.copy(getattr(ls_class, 'eta_n_lo'))
            eta_n_up = np.copy(getattr(ls_class, 'eta_n_up'))

            # get 'eta_n' and 'eta_m' (2-principle direction)
            #
            eta_m2_lo = np.copy(getattr(ls_class, 'eta_m2_lo'))
            eta_m2_up = np.copy(getattr(ls_class, 'eta_m2_up'))
            eta_n2_lo = np.copy(getattr(ls_class, 'eta_n2_lo'))
            eta_n2_up = np.copy(getattr(ls_class, 'eta_n2_up'))

            # add read in saved values to be superposed with currently read in values
            #
            if add_max_min_eta_nm_from_file != None:
                print "superpose max values for 'eta_n' and 'eta_m' with currently loaded values"
                simdata_dir = os.path.join(simdb.simdata_dir, 'lcc_table')
                filename = os.path.join(simdata_dir, add_max_min_eta_nm_from_file)
                max_min_eta_nm_arr = np.loadtxt(filename)
                max_eta_n_arr = max_min_eta_nm_arr[:, 0][:, None]
                min_eta_n_arr = max_min_eta_nm_arr[:, 1][:, None]
                max_eta_m_arr = max_min_eta_nm_arr[:, 2][:, None]
                min_eta_m_arr = max_min_eta_nm_arr[:, 3][:, None]

                # eta_n_lo
                #
                cond_eta_nlo_ge_0 = eta_n_lo >= 0.  # eta caused by tensile normal force
                bool_arr = cond_eta_nlo_ge_0
                eta_n_lo[bool_arr] += max_eta_n_arr[bool_arr]

                cond_eta_nlo_lt_0 = eta_n_lo < 0.  # eta caused by compressive normal force
                bool_arr = cond_eta_nlo_lt_0
                eta_n_lo[bool_arr] += min_eta_n_arr[bool_arr]

                # eta_n_up
                #
                cond_eta_nup_ge_0 = eta_n_up >= 0.  # eta caused by tensile normal force
                bool_arr = cond_eta_nup_ge_0
                eta_n_up[bool_arr] += max_eta_n_arr[bool_arr]

                cond_eta_nup_lt_0 = eta_n_up < 0.  # eta caused by compressive normal force
                bool_arr = cond_eta_nup_lt_0
                eta_n_up[bool_arr] += min_eta_n_arr[bool_arr]

                # eta_m_lo
                #
                cond_eta_mlo_ge_0 = eta_m_lo >= 0.  # eta caused by positive bending moment
                bool_arr = cond_eta_mlo_ge_0
                eta_m_lo[bool_arr] += max_eta_m_arr[bool_arr]

                cond_eta_mlo_lt_0 = eta_m_lo < 0.  # eta caused by negative bending moment
                bool_arr = cond_eta_mlo_lt_0
                eta_m_lo[bool_arr] += min_eta_m_arr[bool_arr]

                # eta_m_up
                #
                cond_eta_mup_ge_0 = eta_m_up >= 0.  # eta caused by positive bending moment
                bool_arr = cond_eta_mup_ge_0
                eta_m_up[bool_arr] += max_eta_m_arr[bool_arr]

                cond_eta_mup_lt_0 = eta_m_up < 0.  # eta caused by negative bending moment
                bool_arr = cond_eta_mup_lt_0
                eta_m_up[bool_arr] += min_eta_m_arr[bool_arr]

            eta_n_lo_list.append(eta_n_lo)
            eta_n_up_list.append(eta_n_up)
            eta_n2_lo_list.append(eta_n2_lo)
            eta_n2_up_list.append(eta_n2_up)

            # NOTE: after superposition take only the absolute values for the plot
            #
            eta_m_lo_list.append(abs(eta_m_lo))
            eta_m_up_list.append(abs(eta_m_up))
            eta_m2_lo_list.append(abs(eta_m2_lo))
            eta_m2_up_list.append(abs(eta_m2_up))

            eta_n_list = eta_n_lo_list + eta_n_up_list + eta_n2_lo_list + eta_n2_up_list
            eta_m_list = eta_m_lo_list + eta_m_up_list + eta_m2_lo_list + eta_m2_up_list

        # stack the list to an array in order to use plot-function
        #
        eta_n_arr = np.hstack(eta_n_list)
        eta_m_arr = np.hstack(eta_m_list)
        print 'eta_n_arr.shape', eta_n_arr.shape

        # save max-values to file in order to superpose them later
        #
        if save_max_min_eta_nm_to_file != None:

            # get maximum values if eta_n is a positive number or minimum if it is a negative number
            #
            max_eta_n_arr = np.max(eta_n_arr, axis=1)  # eta caused by (positive) tensile force

            # eta caused by (negative) compression force
            #
            min_eta_n_arr = np.min(eta_n_arr, axis=1)

            # eta_m cause by a positive bending moment
            #
            max_eta_m_arr = np.max(eta_m_arr, axis=1)

            # eta_m cause by a negative bending moment
            #
            min_eta_m_arr = np.min(eta_m_arr, axis=1)

            # stack as four column array
            #
            max_min_eta_nm_arr = np.hstack([max_eta_n_arr[:, None], min_eta_n_arr[:, None], max_eta_m_arr[:, None], min_eta_m_arr[:, None]])
            # save max values to file
            #
            simdata_dir = os.path.join(simdb.simdata_dir, 'lcc_table')
            filename = os.path.join(simdata_dir, save_max_min_eta_nm_to_file)
            np.savetxt(filename, max_min_eta_nm_arr)
            print 'max_min_eta_nm_arr saved to file %s' % (filename)

        #----------------------------------------------
        # plot
        #----------------------------------------------
        #
        p.figure(facecolor='white')  # white background

        p.plot(eta_m_arr, eta_n_arr, 'wo', markersize=3)  # blue dots
        x = np.array([0, 1. ])
        y1 = np.array([ -1., 0. ])
        y2 = np.array([  1., 0. ])

#        p.title('Ausnutzungsgrad $\eta_{nm}$')

        ax = p.gca()

        if show_tension_only == True:
            ax.set_xticks([0., 0.2, 0.4, 0.6, 0.8, 1.])
            ax.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.])
            p.axis([0., 1., 1., 0.])  # set plotting range for axis

        if show_tension_only == False:
            ax.set_xticks([0., 0.2, 0.4, 0.6, 0.8, 1.])
            ax.set_yticks([1., 0.8, 0.6, 0.4, 0.2, 0, -0.2, -0.4, -0.6, -0.8, -1.])
            p.axis([0., 1., 1., -1.])  # set plotting range for axis

        p.plot(x, y1, 'k--', linewidth=2.0)  # black dashed line
        p.plot(x, y2, 'k--', linewidth=2.0)  # black dashed line
        p.grid(True)

#        ax.spines['left'].set_position(('data', 0))
#        ax.spines['right'].set_color('none')
#        ax.spines['bottom'].set_position(('data', 0))
#        ax.spines['top'].set_color('none')
#        ax.xaxis.set_ticks_position('bottom')
#        ax.yaxis.set_ticks_position('left')
        p.xlabel('$\eta_m$ [-]', fontsize=14)
        p.ylabel('$\eta_n$ [-]', fontsize=14)

        # save figure as png-file
        #
        if save_fig_to_file != None:
            simdata_dir = os.path.join(simdb.simdata_dir, 'lcc_table')
            img_dir = os.path.join(simdata_dir, 'output_images')
            filename = os.path.join(img_dir, save_fig_to_file)
            print 'figure saved to file %s' % (filename)
            p.savefig(filename, format='png')

        p.show()

    # ------------------------------------------------------------
    # View
    # ------------------------------------------------------------

    traits_view = View(VGroup(

                        Item('geo_data_file',
                              label='Evaluated input file for thicknesses ',
                               style='readonly', emphasized=True),
                        VSplit(
                                    Item('lcc_list', editor=lcc_list_editor,
                                          show_label=False),
                                    Item('lcc@', show_label=False),
                                    ),
                        ),
                      resizable=True,
                      scrollable=True,
                      height=1.0,
                      width=1.0
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
    gamma_list = Property(List, depends_on='lc_list_')
    @cached_property
    def _get_gamma_list(self):
        return [[ lc.gamma_fav, lc.gamma_unf ] for lc in self.lc_list ]

    # 'psi' - combination factors (psi) for leading
    # and non leading load cases
    #
    psi_non_lead_arr = Property(Array, depends_on='lc_list_')
    @cached_property
    def _get_psi_non_lead_arr(self):
        return self._get_psi_arr('psi_0')

    psi_lead_arr = Property(Array, depends_on='lc_list_')
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
    gamma_list = Property(List, depends_on='lc_list_')
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
    psi_lead_arr = Property(Array, depends_on='lc_list_, combination_SLS')
    @cached_property
    def _get_psi_lead_arr(self):
        return self.psi_lead_dict[ self.combination_SLS ]

    psi_non_lead_arr = Property(Array, depends_on='lc_list_, combination_SLS')
    @cached_property
    def _get_psi_non_lead_arr(self):
        return self.psi_non_lead_dict[ self.combination_SLS ]

