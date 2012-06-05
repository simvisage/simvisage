'''
Created on Apr 25, 2012

@author: rch
'''

from etsproxy.traits.api import \
    HasTraits, Float, Property, cached_property, Int, Directory, Array
    
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
    time_step_eval = Int(0)
    
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
                            'dog_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4','ARAMIS',
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

        input_list = []
        min_shape = 1e10

        # if a particular time step is selected for evaluation 
        # construct a data_arr consiting of time_step = 0 and time_step = time_step_eval
        #
        if self.time_step_eval != 0:
            print 'NOTE: single time step has been selected for evaluation: t = %g' %(self.time_step_eval) 
            time_step_stop = self.time_step_eval + 1
            time_step_size = self.time_step_eval
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
        print 'data_arr.shape',data_arr.shape
        x_idx_t = np.array(data_arr[:, :, 0], dtype = int)
        y_idx_t = np.array(data_arr[:, :, 1], dtype = int)
#        print 'x_idx_t: ',x_idx_t
#        print 'y_idx_t: ',y_idx_t
        
        # find the shape of the point grid
        # by taking aramis indices stored in
        # first two columns
        x_idx_min_t, x_idx_max_t = np.min(x_idx_t, axis = 1), np.max(x_idx_t, axis = 1)
        y_idx_min_t, y_idx_max_t = np.min(y_idx_t, axis = 1), np.max(y_idx_t, axis = 1)
#        print 'x_idx_min_t: ',x_idx_min_t
#        print 'y_idx_min_t: ',y_idx_min_t
        
#        print 'x_idx_max_t: ',x_idx_max_t
#        print 'y_idx_max_t: ',y_idx_max_t
        
        # make sure that indices start at 0
        x_idx_t -= x_idx_min_t[:,np.newaxis]
        y_idx_t -= y_idx_min_t[:,np.newaxis]
#        print 'x_idx_t_new: ',x_idx_t
#        print 'y_idx_t_new: ',y_idx_t
#        print 'x_idx_max_t_new: ',x_idx_max_t
#        print 'y_idx_max_t_new: ',y_idx_max_t
        
        # number of grid points 
        n_x = np.max(x_idx_max_t - x_idx_min_t) + 1
        n_y = np.max(y_idx_max_t - y_idx_min_t) + 1
        print 'nx',n_x
        print 'ny',n_y
        
        # number of time steps
        n_t = data_arr.shape[0]
        print 'nt',n_t
        t_idx = (np.arange(n_t)[:, np.newaxis] * 
                 np.ones((data_arr.shape[1],), dtype = int)[np.newaxis, :])
#        print 't_idx',t_idx

        # construct the mask for elements to be ignored
        grid_mask_t = np.zeros((n_t, n_x, n_y), dtype = bool)
        grid_mask_t[:, :, :] = True
        grid_mask_t[(t_idx.flatten(), x_idx_t.flatten(), y_idx_t.flatten())] = False
        
        print 'grid_mask_t', grid_mask_t.shape

        return grid_mask_t, t_idx, x_idx_t, y_idx_t

    grid_mask_t = Property
    def _get_grid_mask_t(self):
        return self.idx_maps_t[0]
    

    
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
        grid_mask[( x_idx, y_idx)] = False

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
#        print 'number of missing facettes at time step = 0: ', np.sum(self.idx_maps[0])
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

        
    data_t = Property()
    @cached_property
    def _get_data_t(self):
        '''input data after coordinate transformation'''    
            
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
        x_ = x_vec_ / np.math.sqrt( np.dot( x_vec_, x_vec_ ))                         
        y_ = y_vec_ / np.math.sqrt( np.dot( y_vec_, y_vec_ ))        
        z_ = np.cross( x_, y_)
        
        # base vectors (normed) of the global carthesian coordinate system
        #
        x = np.array([1,0,0])
        y = np.array([0,1,0])
        z = np.array([0,0,1])
        
        # get the direction cosines:
        #
        cos_xx_ = np.dot(x_,x)
        cos_yx_ = np.dot(x_,y)
        cos_zx_ = np.dot(x_,z)
        cos_xy_ = np.dot(y_,x)
        cos_yy_ = np.dot(y_,y)
        cos_zy_ = np.dot(y_,z)
        cos_xz_ = np.dot(z_,x)
        cos_yz_ = np.dot(z_,y)
        cos_zz_ = np.dot(z_,z)
         
        # rotatation using transformation matrix T_mtx 
        # (cf. Zienkiewicz, 6th edition, p.192, (6.18):
        # T_mtx = np.array([[cos_xx_, cos_yx_, cos_zx_],
        #                   [cos_xy_, cos_yy_, cos_zy_],
        #                   [cos_xz_, cos_yz_, cos_zz_]])        
        #
        daf_new[:,:,:,0] = daf_new[:,:,:,0] * cos_xx_ + daf_new[:,:,:,1] * cos_yx_ + daf_new[:,:,:,2] * cos_zx_
        daf_new[:,:,:,1] = daf_new[:,:,:,0] * cos_xy_ + daf_new[:,:,:,1] * cos_yy_ + daf_new[:,:,:,2] * cos_zy_
        daf_new[:,:,:,2] = daf_new[:,:,:,0] * cos_xz_ + daf_new[:,:,:,1] * cos_yz_ + daf_new[:,:,:,2] * cos_zz_

        daf_new[:,:,:,3] = daf_new[:,:,:,3] * cos_xx_ + daf_new[:,:,:,4] * cos_yx_ + daf_new[:,:,:,5] * cos_zx_
        daf_new[:,:,:,4] = daf_new[:,:,:,3] * cos_xy_ + daf_new[:,:,:,4] * cos_yy_ + daf_new[:,:,:,5] * cos_zy_
        daf_new[:,:,:,5] = daf_new[:,:,:,3] * cos_xz_ + daf_new[:,:,:,4] * cos_yz_ + daf_new[:,:,:,5] * cos_zz_
        
        # translation of the coordinates into the origin:
        # 'x_0_vec' derived from the first time step
        # = distance between origin and the position of the first facette
        #
        x_0_vec = daf[0,0,0,:3]
        
#        daf_new[:,:,:,0] = daf_new[:,:,:,0] - x_0_vec[0]
#        daf_new[:,:,:,1] = daf_new[:,:,:,1] - x_0_vec[1]
#        daf_new[:,:,:,2] = daf_new[:,:,:,2] - x_0_vec[2]

        daf_new[:,:,:,:3] = daf_new[:,:,:,:3] - x_0_vec

        return daf_new
       
    
    #===========================================================================
    # Geometry arrays
    #===========================================================================
    x_arr = Property
    @cached_property    
    def _get_x_arr(self):
        # use average to eliminate errors in measuring of single points
#        print 'np.average(self.data_t[self.w_detect_step, :, :, 0], axis = 1)',np.average(self.data_t[self.w_detect_step, :, :, 0], axis = 1)
#        print 'np.average(self.data_t[self.w_detect_step, :, :, 0], axis = 1).shape',np.average(self.data_t[self.w_detect_step, :, :, 0], axis = 1).shape
#        return np.average(self.data_t[self.w_detect_step, :, :, 0], axis = 1)
        # take the x-coords derived from the first time step
        # transform the x-coords to start in the origin
#        print 'self.data_t[0, 0, 0, 0]', self.data_t[0, 0, 0, 0]
#        print 'x_min', np.min(self.data_t[0, :, :, 0])
        return self.data_t[0, :, :, 0] # - self.data_t[0, 0, 0, 0]

    x_idx_arr = Property
    @cached_property    
    def _get_x_idx_arr(self):        
        return np.arange(self.x_arr.shape[0], dtype = int)

    y_arr = Property
    @cached_property    
    def _get_y_arr(self):
        # use average to eliminate errors in measuring of single points
#        print 'np.average(self.data_t[self.w_detect_step, :, :, 1], axis = 0)',np.average(self.data_t[self.w_detect_step, :, :, 1], axis = 0)
#        return np.average(self.data_t[self.w_detect_step, :, :, 1], axis = 0)
        # take the y-coords derived from the first time step
        # transform the y-coords to start in the origin
#        print 'self.data_t[0, 0, 0, 1]', self.data_t[0, 0, 0, 1]
#        print 'self.data_t[0, :, :, 1]', self.data_t[0, :, :, 1]
#        print 'y_min', np.min(self.data_t[0, :, :, 1])
#        print 'self.data_t[0,:, :, 1]- self.data_t[0, 0, 0, 1]', self.data_t[0,:, :, 1] - self.data_t[0, 0, 0, 1]
#        print 'y_min', np.min(self.data_t[0, :, :, 1]) - self.data_t[0, 0, 0, 1]
        return self.data_t[0, :, :, 1] # - self.data_t[0, 0, 0, 1]

    z_arr = Property
    @cached_property    
    def _get_z_arr(self):
        return self.data_t[0, :, :, 2]

    
    y_idx_arr = Property
    @cached_property    
    def _get_y_idx_arr(self):        
        return np.arange(self.y_arr.shape[0], dtype = int)

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
    
    d_ux_w = Property
    @cached_property
    def _get_d_ux_w(self):
        d_ux_w = get_d(self.ux_w, self.integ_radius)
        # cutoff the negative strains - noise
        d_ux_w[ d_ux_w < 0.0 ] = 0.0
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
        return ((dd_ux_avg_w[1:] * dd_ux_avg_w[:-1] < 0.0) * 
                ((ddd_ux_avg_w[1:] + ddd_ux_avg_w[:-1]) / 2.0 < -0.01))

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
    def plot(self):    
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
        
        p.show()
    
    glyph_x_length = Float(0.015)
    glyph_y_length = Float(0.060)
    glyph_z_length = Float(0.000)
    
    def plot3d(self):
        #qargs = [ daf[:, :, i, np.newaxis] for i in range(0, 6)]
        
        #m.quiver3d(*qargs)

#        barr = np.zeros_like(self.w_t)
#        s = m.contour3d(self.t_idx_arr[:, np.newaxis, np.newaxis] + barr,
#                       self.x_idx_arr[np.newaxis, :, np.newaxis] + barr,
#                       self.y_idx_arr[np.newaxis, np.newaxis, :] + barr,
#                       self.d_ux_t, contours = 30)
        
        ## print a surface - why is moved in the y direction?
#        m.surf(self.x_idx_arr, self.y_idx_arr, self.d_ux_w)#, mask = mask_idx_array)
#        m.surf(self.x_arr, self.y_arr, self.d_ux_w)#, mask = mask_idx_array)
        print 'self.x_arr', self.x_arr
        print 'self.y_arr', self.y_arr
#        m.surf(self.x_arr, self.y_arr, self.crack_field_w )#, mask = self.grid_mask)

#        m.points3d(self.x_arr, self.y_arr, self.z_arr, 1.0*self.grid_mask, mode = 'cube', scale_factor = 15., colormap = "blue-red", scale_mode ='none')
        m.points3d(self.x_arr, self.y_arr, self.z_arr, self.d_ux_w, mode = 'cube', scale_factor = 15., colormap = "blue-red", scale_mode ='none')
#        m.points3d(self.x_arr, self.y_arr, self.z_arr, self.crack_field_w, mode = 'cube', scale_factor = 15., colormap = "blue-red", scale_mode ='none')

        engine = m.get_engine()
        glyph = engine.scenes[0].children[0].children[0].children[0]
        glyph.glyph.glyph_source.glyph_source.x_length = self.glyph_x_length
        glyph.glyph.glyph_source.glyph_source.y_length = self.glyph_y_length
        glyph.glyph.glyph_source.glyph_source.z_length = self.glyph_z_length
        scene = engine.scenes[0]
        scene.scene.parallel_projection = True
        scene.scene.z_plus_view()
                        
#        m.show_pipeline(self, engine=engine, rich_view=True)
         
        # plot the missing facettes in a 3d-plot 
        # (manual scaling of the glyphs in the corresponding directions is necessary)
        #
#        facette_gaps = 1.0 * self.grid_mask_t[-1,:,:]
#        print 'facette_gaps', facette_gaps
#        m.points3d(self.x_arr, self.y_arr, self.z_arr, facette_gaps, mode = 'cube', scale_factor = 0.2, colormap = "blue-red", scale_mode ='none')
       
        m.axes() 
        m.show()
        
        # @todo: 
        # sensitivity with respect to the integ_range
        # 
        
if __name__ == '__main__':

    aramis_dir = os.path.join(simdb.exdata_dir, 'tensile_tests', 'dog_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4','ARAMIS')

#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz','Xf15a1-Yf5a3'),
#                     time_step_eval = 3,
#                     integ_radius = 8,
#                     w_detect_step = -1)

#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz','Xf15a1-Yf5a4'),
#                     time_step_eval = 3,
#                     integ_radius = 8,
#                     w_detect_step = -1)

#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'Probe-1-Ausschnitt-Xf15a1-Yf5a4'),
#                     time_step_eval = 428,
#                     integ_radius = 8,
#                     w_detect_step = -1)




#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz','f15a13'),
#                     time_step_eval = 3,
#                     integ_radius = 8,
#                     w_detect_step = -1)
#    
#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1_kurz','f15a3'),
#                     time_step_eval = 3,
#                     integ_radius = 3,
#                     w_detect_step = -1,
#                     glyph_x_length = 0.045,
#                     glyph_y_length = 0.045,
#                     glyph_z_length = 0.000)


    
    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V1'),
                     time_step_eval = 3,
                     integ_radius = 8,
                     w_detect_step = -1)
#    
#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V2'),
#                     time_step_eval = 3,
#                     integ_radius = 8,
#                     w_detect_step = -1)
#
#    ct = CrackTracer(data_dir = os.path.join(aramis_dir, 'V3'),
#                     time_step_eval = 3,
#                     integ_radius = 8,
#                     w_detect_step = -1)
#

    print ct.data_t

#    ct.plot3d()
#    ct.plot()

