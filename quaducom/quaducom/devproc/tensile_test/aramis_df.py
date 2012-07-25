'''
Created on Apr 25, 2012

@author: rch
'''

from etsproxy.traits.api import \
    HasTraits, Float, Property, cached_property, Int, Directory, Array, Bool, Str, Trait
    
import numpy as np
import numpy.ma as ma
import etsproxy.mayavi.mlab as m

from matresdev.db.simdb import \
    SimDB

import os

import os.path
import pylab as p
import string


# Access to the toplevel directory of the database
#
simdb = SimDB()

def get_d(u_arr, integ_radius):
    ir = integ_radius
    du_arr = np.zeros_like(u_arr)
    du_arr[ir:-ir, :] = (u_arr[2 * ir:, :] - 
                         u_arr[:-2 * ir, :])
    return du_arr

class CrackTracer(HasTraits):

    
    # time step used for the evaluation
    # (if defined as non-zero only the defined time step 
    # and the 0 time step are read into 'data_arr'
    #
    evaluated_time_step = Int(0)
    
    # size of a facette in x-direction defined in numbers of pixels 
    #
    n_px_facette_size_x = Int(19)

    # distance between the mid-points of two facettes in x-direction  
    #
    n_px_facette_distance_x = Int(15)
    
    # size of a facette in y-direction defined in numbers of pixels 
    #
    n_px_facette_size_y = Int(19)

    # distance between the mid-points of two facettes in y-direction  
    #
    n_px_facette_distance_y = Int(15)
    
    # size of the steps between the evaluated time steps
    # all time steps found in the read in directory are taken starting
    # with time_step = 0 untill the last time_step found in the directory in steps of
    # 'time_step_size'  
    #
    time_step_size = Int(430)
    
    # time step used to determine the crack pattern
    #
    w_detect_step = Int(-1)
    
    # integration radius for the non-local average of the measured strain
    # defined as integer value of the number of facettes (or elements) 
    # 
    # the value should correspond to ceil( float( n_px_facette_size / n_px_facette_distance ) )
    integ_radius = Int(7)

    #===============================================================================
    # Read data
    #===============================================================================
    #
    # 1) get all files from the directory
    # 2) loop over the files - read also the load level for each time step
    # 3) stack the arrays over time steps an array with the dimensions 
    # input_arr[time,node,values]
    
    
    data_dir = Directory
    def _data_dir_default(self):
        return os.path.join(simdb.exdata_dir, 'tensile_tests',
                            'dog_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4', 'ARAMIS',
                            'Probe-1-Ausschnitt-Xf15a1-Yf5a4')
    
    input_list = Property(Array(float), depends_on = 'data_dir')
    @cached_property
    def _get_input_list(self):
        
        fn_list = os.listdir(self.data_dir)

        # remove hidden directory files from the file name list 
        if '.directory' in fn_list :
            fn_list.remove('.directory')

        n_files = len(fn_list)
        print 'n_files', n_files
#        print 'fn_list', fn_list

        fn_decomposed = [ name + '-' for name in fn_list[0].split('-')[:-1]]
        basename = ''.join(fn_decomposed)

        print 'basename', basename
        self.basename = basename

        input_list = []
        min_shape = 1e10

        # if a particular time step is selected for evaluation 
        # construct a data_arr consiting of time_step = 0 and time_step = evaluated_time_step
        #
        if self.evaluated_time_step != 0:
            print 'NOTE: single time step has been selected for evaluation: t = %g' % (self.evaluated_time_step) 
            time_step_stop = self.evaluated_time_step + 1
            time_step_size = self.evaluated_time_step
        else:
            time_step_stop = n_files
            time_step_size = self.time_step_size

        for fn_idx in range(0, time_step_stop, time_step_size):
#        for fn_idx in range(0, n_files, self.time_step_size):
            fname = os.path.join(self.data_dir, '%s%d.txt' % (basename, fn_idx))
            print 'reading', fn_idx, '...'
            input_arr = np.loadtxt(fname,
                                   skiprows = 14,
                                   usecols = [0, 1, 2, 3, 4, 8, 9, 10])
            # handle the case with unrecognized facets
            # find the maximum row and node number.
            
            min_shape = min(input_arr.shape[0], min_shape)
#            print 'min_shape', min_shape
            input_list.append(input_arr)
        
        # @todo: find the real dimension of the input_arr
        return [input_arr[:min_shape, :] for input_arr in input_list ]

    data_arr = Property()
    @cached_property
    def _get_data_arr(self):
#        print 'first row', input_arr[0, 0, :]
        
        # Identify the points that do not belong to the specimen
        # 
        return np.array(self.input_list)


    idx_maps_t = Property()
    @cached_property
    def _get_idx_maps_t(self):

        data_arr = self.data_arr
        print 'data_arr.shape', data_arr.shape
        x_idx_t = np.array(data_arr[:, :, 0], dtype = int)
        y_idx_t = np.array(data_arr[:, :, 1], dtype = int)
#        print 'x_idx_t: ',x_idx_t
#        print 'y_idx_t: ',y_idx_t
        
        # find the shape of the point grid
        # by taking aramis indices stored in
        # first two columns
        # use minimum and maximum indices of all time steps to derive the grid size
        x_idx_min_t, x_idx_max_t = np.min(x_idx_t, axis = 1), np.max(x_idx_t, axis = 1)
        y_idx_min_t, y_idx_max_t = np.min(y_idx_t, axis = 1), np.max(y_idx_t, axis = 1)
#        print 'x_idx_min_t: ',x_idx_min_t
#        print 'y_idx_min_t: ',y_idx_min_t
#        print 'x_idx_max_t: ',x_idx_max_t
#        print 'y_idx_max_t: ',y_idx_max_t
        
        # make sure that indices start at 0
        x_idx_t -= x_idx_min_t[:, np.newaxis]
        y_idx_t -= y_idx_min_t[:, np.newaxis]
#        print 'x_idx_t_new: ',x_idx_t
#        print 'y_idx_t_new: ',y_idx_t
#        print 'x_idx_max_t_new: ',x_idx_max_t
#        print 'y_idx_max_t_new: ',y_idx_max_t
        
        # number of grid points 
        n_x = np.max(x_idx_max_t - x_idx_min_t) + 1
        n_y = np.max(y_idx_max_t - y_idx_min_t) + 1
        print 'nx', n_x
        print 'ny', n_y
        
        # number of time steps
        n_t = data_arr.shape[0]
        print 'nt', n_t
        t_idx = (np.arange(n_t)[:, np.newaxis] * 
                 np.ones((data_arr.shape[1],), dtype = int)[np.newaxis, :])
#        print 't_idx',t_idx

        # construct the mask for elements to be ignored
        grid_mask_t = np.zeros((n_t, n_x, n_y), dtype = bool)
        grid_mask_t[:, :, :] = True
        grid_mask_t[(t_idx.flatten(), x_idx_t.flatten(), y_idx_t.flatten())] = False
        
        return grid_mask_t, t_idx, x_idx_t, y_idx_t

    grid_mask_t = Property
    def _get_grid_mask_t(self):
        return self.idx_maps_t[0]
    
    grid_mask_w = Property
    def _get_grid_mask_w(self):
        grid_mask_w = self.grid_mask_t[self.w_detect_step, :, :]
        print 'number of missing facettes at crack_detect_step: ', np.sum(grid_mask_w)
        return 1.0 * grid_mask_w
    
    txy_idx = Property
    def _get_txy_idx(self):
        return self.idx_maps_t[1], self.idx_maps_t[2], self.idx_maps_t[3]

    idx_maps = Property()
    @cached_property
    def _get_idx_maps(self):

        data_arr = self.data_arr
#        x_idx = np.array(data_arr[-1, :, 0], dtype = int)
#        y_idx = np.array(data_arr[-1, :, 1], dtype = int)

        # derive the x_idx and y_idx from the undeformed state, 
        # e.g. the first time step (first slice index = 0)
        #
        x_idx = np.array(data_arr[0, :, 0], dtype = int)
        y_idx = np.array(data_arr[0, :, 1], dtype = int)
        
        # find the shape of the point grid
        # by taking aramis indices stored in
        # first two columns
        #
        x_idx_min, x_idx_max = np.min(x_idx), np.max(x_idx)
        y_idx_min, y_idx_max = np.min(y_idx), np.max(y_idx)

        # make sure that indices start at 0
        x_idx -= x_idx_min
        y_idx -= y_idx_min
        
        n_x = x_idx_max - x_idx_min + 1
        n_y = y_idx_max - y_idx_min + 1
        
#        n_t = data_arr.shape[0]
        # construct the mask for elements to be ignored
#        grid_mask = np.zeros((n_t, n_x, n_y), dtype = bool)
#        grid_mask[:, :, :] = True
#        grid_mask[(slice(None), x_idx, y_idx)] = False
        grid_mask = np.zeros((n_x, n_y), dtype = bool)
        grid_mask[:, :] = True
        grid_mask[(x_idx, y_idx)] = False

#        # length of the meassured field in x-direction, pixel size, facette size, facette distance
#        #
#        l_x = self.x_arr[0,0] - self.x_arr[-1,0]
#        px_size_x = l_x / n_x / self.n_px_facette_distance_x
#        facette_size_x = self.n_px_facette_size_x * px_size_x
#        facette_distance_x = self.n_px_facette_distance_x * px_size_x
#        print 'l_x',l_x 
#        print 'self.x_arr[0,0]',self.x_arr[0,0]
##        print 'self.x_arr[0,-1]',np.average(self.x_arr[-1,0]
#        print 'n_x',n_x 
#        print 'px_size_x',px_size_x 
#        print 'facette_size_x',facette_size_x 
#        print 'facette_distance_x',facette_distance_x 
#
#        l_y = self.y_arr[0,0] - self.y_arr[0,-1]
#        px_size_y = l_y / n_y / self.n_px_facette_distance_y
#        facette_size_y = self.n_px_facette_size_y * px_size_y
#        facette_distance_y = self.n_px_facette_distance_y * px_size_y
#        print 'l_y',l_y 
#        print 'n_y',n_y 
#        print 'px_size_y',px_size_y 
#        print 'facette_size_y',facette_size_y 
#        print 'facette_distance_y',facette_distance_y 

        return grid_mask, x_idx, y_idx, n_x, n_y


    grid_mask = Property
    def _get_grid_mask(self):
        return self.idx_maps[0]
    
    xy_idx = Property
    def _get_xy_idx(self):
        return self.idx_maps[1], self.idx_maps[2]
        
    data_t_orig = Property()
    @cached_property
    def _get_data_t_orig(self):
        '''original input data before coordinate transformation'''
        data_arr = self.data_arr
        
        # total number of data values associated with a single point
        n_values = data_arr.shape[2] - 2
    
        # define a three dimensional array
        # first two indices are grid indices, third
        # index selects one of the values associated with a point
        daf = np.zeros(self.grid_mask_t.shape + (n_values,), dtype = float)
        
        # define an array of indexes performing the mapping
        # from the data_arr into daf
        select_idx_arr = (self.txy_idx[0], self.txy_idx[1], self.txy_idx[2], slice(None))

        daf[select_idx_arr] = data_arr[:, :, 2:] 

        #data_t = ma.masked_array(daf, mask = self.grid_mask)

        print 'data_t.shape', daf.shape
        return daf

    transform_data = 'true'
    
    data_t = Property()
    @cached_property
    def _get_data_t(self):
        '''input data after coordinate transformation'''    
        
        if self.transform_data == 'true':    
            # Coordinate transformation:
            # derive transformation direction from the first time step (t = 0)
            #
            daf_0 = self.data_t_orig[0]
            daf = self.data_t_orig
            daf_new = np.copy(daf)
            
            # @todo: select a facette with none-zero coordinate values
            # e.g. select a facette from the last middle axis and in a list of the first 10 or last 10 facettes
            # 
    #        n_x = self.idx_maps[3]
    #        n_y = self.idx_maps[4]
    #        idx_x = self.idx_maps[1]
    #        idx_y = self.idx_maps[2]
    ##        n_x_middle = int( n_x/2 )
    ##        n_y_middle = int( n_y/2 )
    #        x_idx_middle = int( idx_x.shape[0] / 2 )
    #        y_idx_middle = int( idx_y.shape[0] / 2 )
    #        n_x_middle = idx_x[ x_idx_middle ]
    #        n_y_middle = idx_y[ y_idx_middle ]
    #        print 'n_x_middle', n_x_middle
    #        print 'n_y_middle', n_y_middle
            n_x_middle = 10
            n_y_middle = 10
    
            x_vec_ = daf_0[-1, n_x_middle, :3] - daf_0[0, n_x_middle, :3]
            y_vec_ = daf_0[n_y_middle, -1, :3] - daf_0[n_y_middle, 0, :3]
    
    #        x_vec_ = daf_0[-1, 10, :3] - daf_0[0, 10, :3]
    #        y_vec_ = daf_0[10, -1, :3] - daf_0[10, 0, :3]
    
            # base vectors (normed) of the local coordinate system x_, y_, z_
            #
            x_ = x_vec_ / np.math.sqrt(np.dot(x_vec_, x_vec_))                         
            y_ = y_vec_ / np.math.sqrt(np.dot(y_vec_, y_vec_))        
            z_ = np.cross(x_, y_)
            
            # base vectors (normed) of the global carthesian coordinate system
            #
            x = np.array([1, 0, 0])
            y = np.array([0, 1, 0])
            z = np.array([0, 0, 1])
            
            # get the direction cosines:
            #
            cos_xx_ = np.dot(x_, x)
            cos_yx_ = np.dot(x_, y)
            cos_zx_ = np.dot(x_, z)
            cos_xy_ = np.dot(y_, x)
            cos_yy_ = np.dot(y_, y)
            cos_zy_ = np.dot(y_, z)
            cos_xz_ = np.dot(z_, x)
            cos_yz_ = np.dot(z_, y)
            cos_zz_ = np.dot(z_, z)
             
            # rotation using transformation matrix T_mtx 
            # (cf. Zienkiewicz, 6th edition, p.192, (6.18):
            # T_mtx = np.array([[cos_xx_, cos_yx_, cos_zx_],
            #                   [cos_xy_, cos_yy_, cos_zy_],
            #                   [cos_xz_, cos_yz_, cos_zz_]])        
            # notation:
            # x,y,z = global cartesian coordinates
            # x_,y_,z_ = rotated coordinate system
            # cos_xz_ = direction cosine between axis 'x' and 'z_'
            
            # rotate the coordinate values (x,y,z)    
            #
            daf_new[:, :, :, 0] = daf_new[:, :, :, 0] * cos_xx_ + daf_new[:, :, :, 1] * cos_yx_ + daf_new[:, :, :, 2] * cos_zx_
            daf_new[:, :, :, 1] = daf_new[:, :, :, 0] * cos_xy_ + daf_new[:, :, :, 1] * cos_yy_ + daf_new[:, :, :, 2] * cos_zy_
            daf_new[:, :, :, 2] = daf_new[:, :, :, 0] * cos_xz_ + daf_new[:, :, :, 1] * cos_yz_ + daf_new[:, :, :, 2] * cos_zz_
    
            # rotate the displacement values (ux,uy,uz)    
            #
            daf_new[:, :, :, 3] = daf_new[:, :, :, 3] * cos_xx_ + daf_new[:, :, :, 4] * cos_yx_ + daf_new[:, :, :, 5] * cos_zx_
            daf_new[:, :, :, 4] = daf_new[:, :, :, 3] * cos_xy_ + daf_new[:, :, :, 4] * cos_yy_ + daf_new[:, :, :, 5] * cos_zy_
            daf_new[:, :, :, 5] = daf_new[:, :, :, 3] * cos_xz_ + daf_new[:, :, :, 4] * cos_yz_ + daf_new[:, :, :, 5] * cos_zz_
            
            # translation of the coordinates into the origin:
            # 'x_0_vec' derived from the first time step
            # = distance between origin and the position of the first facette
            # @todo: make sure that facette (0,0) is non-zero due to lost facette!
            #
            x_0_vec = daf[0, 0, 0, :3]
            print 'x_0_vec', x_0_vec
            
    #        daf_new[:,:,:,0] = daf_new[:,:,:,0] - x_0_vec[0]
    #        daf_new[:,:,:,1] = daf_new[:,:,:,1] - x_0_vec[1]
    #        daf_new[:,:,:,2] = daf_new[:,:,:,2] - x_0_vec[2]
    
            daf_new[:, :, :, :3] = daf_new[:, :, :, :3] - x_0_vec
    
            return daf_new
        else:
            return self.data_t_orig
    
    #===========================================================================
    # Geometry arrays
    #===========================================================================
    x_arr = Property
    @cached_property    
    def _get_x_arr(self):
        # take the x-coords from the w_detect_step time step
        return self.data_t[self.w_detect_step, :, :, 0] 

    x_arr_avg = Property
    @cached_property    
    def _get_x_arr_avg(self):
        # use average to eliminate errors in measuring of single points (yields a 1d-array)
        return np.average(self.data_t[self.w_detect_step, :, :, 0], axis = 1)

    x_idx_arr = Property
    @cached_property    
    def _get_x_idx_arr(self):        
        return np.arange(self.x_arr.shape[0], dtype = int)

    y_arr = Property
    @cached_property    
    def _get_y_arr(self):
        # take the y-coords derived from the w_detect_step time step
        return self.data_t[self.w_detect_step, :, :, 1]

    y_arr_avg = Property
    @cached_property    
    def _get_y_arr_avg(self):
        # use average to eliminate errors in measuring of single points (yields a 1d-array)
        return np.average(self.data_t[self.w_detect_step, :, :, 1], axis = 0)

    y_idx_arr = Property
    @cached_property    
    def _get_y_idx_arr(self):        
        return np.arange(self.y_arr.shape[0], dtype = int)

    l_x = Property
    @cached_property    
    def _get_l_x(self):
        return abs(self.x_arr_avg[-1] - self.x_arr_avg[0])

    z_arr = Property
    @cached_property    
    def _get_z_arr(self):
        return self.data_t[0, :, :, 2]

    t_arr = Property
    @cached_property
    def _get_t_arr(self):
        # get the load factor
        return np.arange(self.data_t.shape[0], dtype = float)

    t_idx_arr = Property
    @cached_property    
    def _get_t_idx_arr(self):        
        return np.arange(self.t_arr.shape[0], dtype = int)

    #===========================================================================
    # Fields in space and time
    #===========================================================================
    ux_t = Property
    @cached_property
    def _get_ux_t(self):        
        ux_arr = self.data_t[:, :, :, 3]
        # use an individual grid_mask for each time step
        # 
        ux_masked = ma.masked_array(ux_arr, mask = self.grid_mask_t)

        ux_avg = ma.average(ux_masked, axis = 2)
        t_idx_zeros, x_idx_zeros, y_idx_zeros = np.where(self.grid_mask_t)
        ux_arr[t_idx_zeros, x_idx_zeros, y_idx_zeros] = ux_avg[t_idx_zeros, x_idx_zeros]
        
        return ux_arr
    
    d_ux_t = Property
    @cached_property
    def _get_d_ux_t(self):
        # generate the elem_node_map
        
        ir = self.integ_radius
        ux_arr = self.ux_t
        d_ux_arr = np.zeros_like(ux_arr)
        d_ux_arr[:, ir:-ir, :] = (ux_arr[:, 2 * ir:, :] - 
                                   ux_arr[:, :-2 * ir, :])
        
        # cutoff the negative strains - noise
        d_ux_arr[ d_ux_arr < 0.0 ] = 0.0
        
        return d_ux_arr
    
    d_ux_t_avg = Property
    @cached_property
    def _get_d_ux_t_avg(self):
        return np.average(self.d_ux_t, axis = 2)

    #===========================================================================
    # Fields in the time step used for crack detection
    #===========================================================================
    ux_w = Property
    @cached_property
    def _get_ux_w(self):
        # generate the elem_node_map
        return self.ux_t[self.w_detect_step]
    
    ux_w_avg = Property
    @cached_property
    def _get_ux_w_avg(self):
        # generate the elem_node_map
        return np.average(self.ux_w, axis = 1)
    
    # tolerance to distinguish cracks from noise
    # i.e. use only displacement
    #
    d_ux_tol = Float(0.00)
    
    d_ux_w = Property
    @cached_property
    def _get_d_ux_w(self):
        d_ux_w = get_d(self.ux_w, self.integ_radius)
        # cutoff the negative strains - noise
        d_ux_w[ d_ux_w < 0.0 ] = 0.0
        # eliminate displacement jumps smaller then specified tolerance
        d_ux_w[ d_ux_w < self.d_ux_tol ] = 0.0
        return d_ux_w
        
    d_ux_avg_w = Property
    @cached_property
    def _get_d_ux_avg_w(self):
        return np.average(self.d_ux_w, axis = 1)
        
    dd_ux_w = Property
    @cached_property
    def _get_dd_ux_w(self):
        return get_d(self.d_ux_w, self.integ_radius)

    dd_ux_avg_w = Property
    @cached_property
    def _get_dd_ux_avg_w(self):
        return np.average(self.dd_ux_w, axis = 1)
        
    ddd_ux_w = Property
    @cached_property
    def _get_ddd_ux_w(self):
        return get_d(self.dd_ux_w, self.integ_radius)

    ddd_ux_avg_w = Property
    @cached_property
    def _get_ddd_ux_avg_w(self):
        return np.average(self.ddd_ux_w, axis = 1)
        
    #===========================================================================
    # Crack detection
    #===========================================================================
    crack_filter = Property
    @cached_property
    def _get_crack_filter(self):
        dd_ux_w = self.dd_ux_w 
        ddd_ux_w = self.ddd_ux_w
        return ((dd_ux_w[1:, :] * dd_ux_w[:-1, :] < 0.0) * 
                ((ddd_ux_w[1:, :] + ddd_ux_w[:-1, :]) / 2.0 < -0.005)) 
    
    crack_filter_avg = Property
    @cached_property
    def _get_crack_filter_avg(self):
        dd_ux_avg_w = self.dd_ux_avg_w
        ddd_ux_avg_w = self.ddd_ux_avg_w
        crack_filter_avg = ((dd_ux_avg_w[1:] * dd_ux_avg_w[:-1] < 0.0) * 
                           ((ddd_ux_avg_w[1:] + ddd_ux_avg_w[:-1]) / 2.0 < -0.01))
        print "number of cracks determined by 'crack_filter_avg': ", np.sum(crack_filter_avg)
        return crack_filter_avg

    crack_spacing_avg = Property
    @cached_property
    def _get_crack_spacing_avg(self):
        n_cr_avg = np.sum(self.crack_filter_avg)
        s_cr_avg = self.l_x / n_cr_avg
        print "average crack spacing [mm]: %.1f" % (s_cr_avg)
        return s_cr_avg 
        
    crack_arr_w = Property
    @cached_property
    def _get_crack_arr_w(self):
        return self.d_ux_w[np.where(self.crack_filter)]

    crack_field_w = Property
    @cached_property
    def _get_crack_field_w(self):    
        cf_w = np.zeros_like(self.d_ux_w)
        cf_w[np.where(self.crack_filter)] = self.crack_arr_w
        return cf_w

    crack_arr_t = Property
    @cached_property
    def _get_crack_arr_t(self):
        return self.d4_ux_t[np.where(self.crack_filter[None, :, :])]

    crack_field_t = Property
    @cached_property
    def _get_crack_field_t(self):    
        '''
        '''
        crack_idx_arr = np.where(self.crack_filter[None, :, :])
        crack_arr_t = self.d_ux_t[crack_idx_arr]
        cf_t = np.zeros_like(self.d_ux_t)
        cf_t[np.where(self.crack_filter[None, :, :])] = crack_arr_t
        return cf_t

    #===========================================================================
    # Plotting methods
    #===========================================================================

    def plot2d(self):    
        p.subplot(2, 2, 1)
        p.plot(self.x_idx_arr[:-1], self.crack_filter_avg * 0.1,
               color = 'magenta', linewidth = 2)
        p.plot(self.x_idx_arr, self.d_ux_avg_w, color = 'black')

        p.subplot(2, 2, 2)
        p.plot(self.x_idx_arr, self.ux_w, color = 'green')
        p.plot(self.x_idx_arr, self.ux_w_avg, color = 'red')
        
        p.subplot(2, 2, 3)
        p.plot(self.x_idx_arr[:-1], self.crack_filter_avg * 0.1,
               color = 'magenta', linewidth = 2)
        p.plot(self.x_idx_arr, self.dd_ux_avg_w, color = 'black')
        p.plot(self.x_idx_arr, self.ddd_ux_avg_w, color = 'blue')
        
        p.subplot(2, 2, 4)
        p.hist(self.crack_arr_w, bins = 40)
        
        # figure title 
        # uses directory path by default
        #
        if plot_title == 'true':
            m.title(os.path.join(self.data_dir, self.basename))
        
        p.show()
    
    def plot_w_hist(self):
        font = {'family' : 'times new roman',
                #'weight' : 'bold',
                'size'   : 20 }
        
        p.rc('font', **font)  # pass in the font dict as kwargs
        p.hist(self.crack_arr_w, normed = True, range = (0, 0.25),
               bins = 40, histtype = 'bar', color = 'gray')
        p.xlabel('w [mm]')
        p.ylabel('frequency [%]')
        p.twinx()
        p.hist(self.crack_arr_w, normed = True, range = (0, 0.25),
               histtype = 'step', color = 'black',
               cumulative = True, linewidth = 3, bins = 40)
        p.ylim(0, 1)
        p.ylabel('fraction [-]')
        # figure title 
        # uses directory path by default
        #
        if plot_title == 'true':
            m.title(os.path.join(self.data_dir, self.basename))
        
        p.show()
    
    
    # variable to be plotted in 3d-plot
    # set variable to 'd_ux_w', 'd_ux_w'
    #
    plot3d_var = Trait('d_ux_w', {'ux_w [mm]':('ux_w', None),
                                  'd_ux_w [mm]':('d_ux_w', 0.2),
                                  'w_cr [mm]': ('crack_field_w', 0.2),
                                  'grid_mask_w': ('grid_mask_w', 1.0)})

    # flag for figure title
    #
    plot_title = 'false'    

    # variables to configure glyphs for 'plot3d_points'
    #
    glyph_x_length = Float(0.015)
    glyph_y_length = Float(0.060)
    glyph_z_length = Float(0.000)

    def plot3d_points(self):

        # set background color to white and forground color to black
        #
        m.figure(fgcolor = (0, 0, 0), bgcolor = (1, 1, 1), size = (900, 600))                

        # NOTE: the missing facettes have all the coordinate position (0,0,0); 
        # the missing facettes therefore are displayed in points3d as missing elements (gaps)
        
        plot3d_var = getattr(self, self.plot3d_var_[0])
        vmax = self.plot3d_var_[1]
        
        m.points3d(self.x_arr, self.y_arr, self.z_arr, plot3d_var, mode = 'cube', colormap = "blue-red", scale_mode = 'none', vmax = vmax)

        # switch coordinate order to display glyphs at their head/tail insted of their center
#        m.points3d(self.z_arr, self.x_arr, self.y_arr, plot3d_var, mode = 'cube', colormap = "blue-red", scale_mode ='none', vmax = vmax)

        # plot scalarbar
        #
        m.scalarbar(orientation = 'horizontal', title = self.plot3d_var)
        # @todo: predefine scalarbar position
#        m.scalarbar().position = [0., 0.]

        # plot axes
        #
        m.axes() 
        
        # figure title 
        # uses directory path by default
        #
        if plot_title == 'true':
            m.title(os.path.join(self.data_dir, self.basename))

        # get mayavi engine 
        #
        engine = m.get_engine()
        
        # scale glyphs 
        #
        glyph = engine.scenes[0].children[0].children[0].children[0]
        glyph.glyph.glyph_source.glyph_source.x_length = self.glyph_x_length
        glyph.glyph.glyph_source.glyph_source.y_length = self.glyph_y_length
        glyph.glyph.glyph_source.glyph_source.z_length = self.glyph_z_length
        glyph.glyph.glyph_source.glyph_position = 'tail'

        # rotate scene 
        #
        scene = engine.scenes[0]
        scene.scene.parallel_projection = True
        scene.scene.z_plus_view()

        # preset position of scalarbar
        #
        module_manager = engine.scenes[0].children[0].children[0]
        module_manager.scalar_lut_manager.scalar_bar_representation.maximum_size = np.array([100000, 100000])
        module_manager.scalar_lut_manager.scalar_bar_representation.minimum_size = np.array([1, 1])
        module_manager.scalar_lut_manager.scalar_bar_representation.position2 = np.array([ 0.20603604, 0.16827586])
        module_manager.scalar_lut_manager.scalar_bar_representation.moving = 0
        module_manager.scalar_lut_manager.scalar_bar_representation.position = np.array([ 0.40858859, 0.09310345])
        module_manager.scalar_lut_manager.scalar_bar_representation.maximum_size = np.array([100000, 100000])

#        # @todo: use mlab-functionality instead
#        view = mlab.view()
#        roll = mlab.roll()
#        # Reposition the camera
#        mlab.view(*view)
#        mlab.roll(roll)
#        
#        f = mlab.gcf()
#        camera = f.scene.camera
#        cam.parallel_scale = 9
        
        m.show()
        
    # variables to configure glyphs for 'plot3d_cracks'
    #
    glyph_x_length_cr = Float(3.000)
    glyph_y_length_cr = Float(0.120)
    glyph_z_length_cr = Float(0.120)

    def plot3d_cracks(self):

        # set background color to white and forground color to black
        #
        m.figure(fgcolor = (0, 0, 0), bgcolor = (1, 1, 1), size = (900, 600))                

        #-----------------------------------
        # plot crack width ('crack_field_w')
        #-----------------------------------
        
        # flatten z_arr
        # @todo: use instead: extent = [xmin, xmax, ymin, ymax, zmin, zmax], 
        #
        z_arr = np.zeros_like(self.z_arr)
        # NOTE: coordinate order is switched in order to display glyphs at their tail instead of their center
        # (only possible in mayavis's x-direction)    
        #
        plot3d_var = getattr(self, self.plot3d_var_[0])
#        plot3d_var = getattr(self, 'crack_field_w')
        vmax = self.plot3d_var_[1]
        m.points3d(z_arr, self.x_arr, self.y_arr, plot3d_var, mode = 'cube', colormap = "blue-red", scale_mode = 'scalar', vmax = vmax)

        # plot scalarbar
        #
        m.scalarbar(orientation = 'horizontal', title = self.plot3d_var)

        # plot axes
        #
#        m.axes() 

        # figure title 
        # uses directory path by default
        #
        if plot_title == 'true':
            m.title(os.path.join(self.data_dir, self.basename))
        
        # get mayavi engine
        #
        engine = m.get_engine()

        # scale glyphs
        #
        glyph = engine.scenes[0].children[0].children[0].children[0]
        glyph.glyph.glyph_source.glyph_position = 'tail'
        glyph.glyph.glyph_source.glyph_source.x_length = self.glyph_x_length_cr
        glyph.glyph.glyph_source.glyph_source.y_length = self.glyph_y_length_cr
        glyph.glyph.glyph_source.glyph_source.z_length = self.glyph_z_length_cr
        
        #-----------------------------------
        # plot displacement jumps ('d_ux_w')
        #-----------------------------------

#        plot3d_var = getattr(self, 'd_ux_w')
        plot3d_var = getattr(self, 'crack_field_w')
        m.points3d(z_arr, self.x_arr, self.y_arr, plot3d_var, mode = 'cube', colormap = "blue-red", scale_mode = 'none', vmax = vmax)
        glyph1 = engine.scenes[0].children[1].children[0].children[0]
        # switch order of the scale_factor corresponding to the order of the 
        glyph1.glyph.glyph_source.glyph_source.x_length = self.glyph_z_length
        glyph1.glyph.glyph_source.glyph_source.y_length = self.glyph_x_length
        glyph1.glyph.glyph_source.glyph_source.z_length = self.glyph_y_length
        
        # rotate scene
        #
        scene = engine.scenes[0]
        scene.scene.parallel_projection = False
        scene.scene.camera.position = [616.49832063929375, 638.29074243238438, 514.06081220962164]
        scene.scene.camera.focal_point = [11.259753942489624, 11.990119934082031, 9.7502956390380859]
        scene.scene.camera.view_angle = 30.0
        scene.scene.camera.view_up = [0.79246046934594205, -0.54446721176015089, -0.27488517574095594]
        scene.scene.camera.clipping_range = [595.49259262137014, 1526.1362976999562]
        scene.scene.camera.compute_view_plane_normal()
        scene.scene.render()

        glyph.glyph.glyph_source.glyph_position = 'head'
        glyph.glyph.glyph_source.glyph_position = 'tail'
        
        module_manager = engine.scenes[0].children[0].children[0]
        module_manager.scalar_lut_manager.scalar_bar_representation.minimum_size = np.array([1, 1])
        module_manager.scalar_lut_manager.scalar_bar_representation.position2 = np.array([ 0.20603604, 0.16827586])
        module_manager.scalar_lut_manager.scalar_bar_representation.moving = 0
        module_manager.scalar_lut_manager.scalar_bar_representation.position = np.array([ 0.53971972, 0.19931035])
        module_manager.scalar_lut_manager.scalar_bar_representation.maximum_size = np.array([100000, 100000])
        
        m.show()
        
        
    # variable to reduce plot to non-missing facettes (mask)
    #
    use_mask = 'false'
    
    # variable scale function value in z-dirction of the plot
    #
    warp_scale = 1.
    
    def plot3d_surf(self):

        # set background color to white and forground color to black
        #
        m.figure(fgcolor = (0, 0, 0), bgcolor = (1, 1, 1), size = (1600, 900))                

        #----------------------------------------------------    
        # plot 'plot_var' as 3d-surface plot 
        #----------------------------------------------------    
        # @todo: is this comment still relavant? "why is moved in the y direction?"
        
        if self.use_mask == 'true':
            mask_arr = self.grid_mask_w
        else:
            mask_arr = None

        plot3d_var = getattr(self, self.plot3d_var_[0])
        vmax = self.plot3d_var_[1]

        m.surf(self.x_arr, self.y_arr, plot3d_var, vmax = vmax, mask = mask_arr, warp_scale = warp_scale)

        # use average coordinate values (1D-arrays)
        #    
#        m.surf(self.x_arr_avg, self.y_arr_avg, plot3d_var, vmax = vmax, mask = mask_arr, warp_scale = warp_scale)


        # plot scalarbar
        #
        m.scalarbar(orientation = 'horizontal', title = self.plot3d_var)

        # plot axes
        #
#        m.axes() 

        # figure title 
        # uses directory path by default
        #
        if plot_title == 'true':
            m.title(os.path.join(self.data_dir, self.basename))

        # get mayavi engine
        #
        engine = m.get_engine()

        # rotate scene
        #
        scene = engine.scenes[0]
        scene.scene.parallel_projection = True
        if self.warp_scale == 1:
            scene.scene.z_plus_view()
        else:
            scene.scene.isometric_view()
            scene.scene.camera.position = [-319.07443227647047, -469.29101319337502, 434.40818470843612]
            scene.scene.camera.focal_point = [-66.983721840860625, -79.447548413824947, 1.1328650920308507]
            scene.scene.camera.view_angle = 30.0
            scene.scene.camera.view_up = [0.49290441759866288, 0.48449293877207339, 0.72271144130401255]
            scene.scene.camera.clipping_range = [325.28180250321782, 995.97136098229157]
            scene.scene.camera.compute_view_plane_normal()
            scene.scene.render()
        
        # predefine the position of the scalarbar
        #
        module_manager = engine.scenes[0].children[0].children[0].children[0].children[0]
        module_manager.scalar_lut_manager.scalar_bar_representation.minimum_size = np.array([1, 1])
        module_manager.scalar_lut_manager.scalar_bar_representation.position2 = np.array([ 0.27780226, 0.11638842])
        module_manager.scalar_lut_manager.scalar_bar_representation.position = np.array([ 0.59707606, 0.16105125])
        module_manager.scalar_lut_manager.scalar_bar_representation.maximum_size = np.array([100000, 100000])        
        m.show()
        
    # @todo: 
    # sensitivity with respect to the integ_range
        
if __name__ == '__main__':

    aramis_dir = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4', 'ARAMIS')

#-----------------------------------------------------------------------------------
# select plot typ
#-----------------------------------------------------------------------------------

#    plot_type = '3d-surf'
#    plot_type = '3d-points'
#    plot_type = '3d-cracks'
#    plot_type = '2d'
    plot_type = 'w_hist'

#-----------------------------------------------------------------------------------
# select variable to be plotted
#-----------------------------------------------------------------------------------
    
#    plot3d_var = 'ux_w [mm]'
#    plot3d_var = 'd_ux_w [mm]'
    plot3d_var = 'w_cr [mm]'
#    plot3d_var = 'grid_mask_w'

#-----------------------------------------------------------------------------------
# set mask flag for 'plot3d-surf' and 'warp_scale'
#-----------------------------------------------------------------------------------
    
#    use_mask = 'true'
    use_mask = 'false'
    warp_scale = 100.
#    warp_scale = 1.
    
#-----------------------------------------------------------------------------------
# plot title
#-----------------------------------------------------------------------------------
    
    plot_title = 'false'
#    plot_title = 'true'

#-----------------------------------------------------------------------------------
# complete displacement field; varying resolutions, quadratic elements
#-----------------------------------------------------------------------------------

#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz','f5a3'),
#                     evaluated_time_step = 22,
#                     integ_radius = 1,
#                     w_detect_step = -1,
#                     transform_data = 'false',
#                     plot3d_var = plot3d_var,
#                     plot_title = plot_title,
#                     use_mask = use_mask,
#                     warp_scale = warp_scale,
#                     glyph_x_length = 0.05,
#                     glyph_y_length = 0.05
#                     )

#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz','f7a3'),
#                     evaluated_time_step = 3,
#                     integ_radius = 2,
#                     w_detect_step = -1,
#                     transform_data = 'false',
#                     plot3d_var = plot3d_var,
#                     plot_title = plot_title,
#                     use_mask = use_mask,
#                     warp_scale = warp_scale,
#                     glyph_x_length = 0.05,
#                     glyph_y_length = 0.05
#                    )

#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz','f8a4'),
#                     evaluated_time_step = 3,
#                     integ_radius = 2,
#                     w_detect_step = -1,
#                     transform_data = 'false',
#                     plot3d_var = plot3d_var,
#                     plot_title = plot_title,
#                     use_mask = use_mask,
#                     warp_scale = warp_scale,
#                     glyph_x_length = 0.05,
#                     glyph_y_length = 0.05
#                     )
  
    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz', 'f15a3'),
                     evaluated_time_step = 3,
                     integ_radius = 3,
                     w_detect_step = -1,
                     transform_data = 'true',
                     plot3d_var = plot3d_var,
                     plot_title = plot_title,
                     use_mask = use_mask,
                     warp_scale = warp_scale,
                     glyph_x_length = 0.05,
                     glyph_y_length = 0.05,
                     glyph_x_length_cr = 2.000,
                     glyph_y_length_cr = 0.120,
                     glyph_z_length_cr = 0.200
                     )

#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz','f15a13'),
#                     evaluated_time_step = 3,
#                     integ_radius = 2,
#                     w_detect_step = -1,
#                     transform_data = 'false',
#                     plot3d_var = plot3d_var,
#                     plot_title = plot_title,
#                     use_mask = use_mask,
#                     warp_scale = warp_scale,
#                     glyph_x_length = 0.12,
#                     glyph_y_length = 0.12
#                     )

#-------------------------------------------------------------------------------------------------
# complete displacement field; varying resolutions, rectangular elements with fine x-dircetion)
#-------------------------------------------------------------------------------------------------

#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz','Xf15a1-Yf5a3'),
#                     evaluated_time_step = 3,
#                     integ_radius = 8,
#                     w_detect_step = -1,
#                     transform_data = 'false',
#                     plot3d_var = plot3d_var,
#                     plot_title = plot_title,
#                     use_mask = use_mask,
#                     warp_scale = warp_scale,
#                     glyph_x_length = 0.013,
#                     glyph_y_length = 0.045
#                     )
#
#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz','Xf15a1-Yf5a4'),
#                     evaluated_time_step = 3,
#                     integ_radius = 8,
#                     w_detect_step = -1,
#                     transform_data = 'false',
#                     plot3d_var = plot3d_var,
#                     plot_title = plot_title,
#                     use_mask = use_mask,
#                     warp_scale = warp_scale,
#                     glyph_x_length = 0.013,
#                     glyph_y_length = 0.045
#                     )


#-----------------------------------------------------------------------------------
# 1/4 of the displacement field; fine resolution x-direction (Xf15a1-Yf5a4)
#-----------------------------------------------------------------------------------

#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'Probe-1-Ausschnitt-Xf15a1-Yf5a4'),
#                     evaluated_time_step = 432,
#                     integ_radius = 8,
#                     w_detect_step = -1,
#                     transform_data = 'false',
#                     plot3d_var = plot3d_var,
#                     plot_title = plot_title,
#                     use_mask = use_mask,
#                     warp_scale = warp_scale,
#                     glyph_x_length = 0.025,
#                     glyph_y_length = 0.120,
#                     glyph_z_length = 0.120,
#                     glyph_x_length_cr = 3.,
#                     glyph_y_length_cr = 0.120,
#                     glyph_z_length_cr = 0.120
#                     )

#-----------------------------------------------------------------------------------
# complete displacement field; rough resolution (f19a15)
#-----------------------------------------------------------------------------------
   
#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1'),
#                     evaluated_time_step = 430,
#                     integ_radius = 1,
#                     w_detect_step = -1,
#                     transform_data = 'true',
#                     plot3d_var = plot3d_var,
#                     plot_title = plot_title,
#                     use_mask = use_mask,
#                     warp_scale = warp_scale,
#                     glyph_x_length = 0.2,
#                     glyph_y_length = 0.2
#                     )
#    
#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V2'),
#                     evaluated_time_step = 430,
#                     integ_radius = 1,
#                     w_detect_step = -1,
#                     transform_data = 'false',
#                     plot3d_var = plot3d_var,
#                     plot_title = plot_title,
#                     use_mask = use_mask,
#                     warp_scale = warp_scale,
#                     glyph_x_length = 0.2,
#                     glyph_y_length = 0.2
#                     )
#
#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V3'),
#                     evaluated_time_step = 430,
#                     integ_radius = 1,
#                     w_detect_step = -1,
#                     transform_data = 'false',
#                     plot3d_var = plot3d_var,
#                     plot_title = plot_title,
#                     use_mask = use_mask,
#                     warp_scale = warp_scale,
#                     glyph_x_length = 0.2,
#                     glyph_y_length = 0.2
#                     )
#

#    print ct.data_t

    ct.grid_mask_w
    ct.crack_filter_avg
    ct.crack_spacing_avg
    print 'l_x', ct.l_x
    
    if plot_type == '3d-surf':
        ct.plot3d_surf()
    elif plot_type == '3d-points':
        ct.plot3d_points()
    if plot_type == '3d-cracks':
        ct.plot3d_cracks()
    elif plot_type == '2d':
        ct.plot2d()
    elif plot_type == 'w_hist':
        ct.plot_w_hist()

