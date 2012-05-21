'''
Created on Apr 25, 2012

@author: rch
'''

from enthought.traits.api import \
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

def get_d4(u_arr, integ_radius):
    ir = integ_radius
    du_arr = np.zeros_like(u_arr)
    du_arr[ir:-ir, :] = (u_arr[2 * ir:, :] - 
                         u_arr[:-2 * ir, :])
    return du_arr

class CrackTracer(HasTraits):

    time_step_size = Int(430)
    
    w_detect_step = Int(-1)
    
    integ_radius = Int(7)

    #===============================================================================
    # Read data
    #===============================================================================
    #
    # 1) get all files from the directory
    # 2) loop over the files - read also the load level for each load step
    # 3) stack the arrays over load steps an array with the dimensions 
    # input_arr[time,node,values]
    # 
    
    
    data_dir = Directory
    def _data_dir_default(self):
        return os.path.join(simdb.exdata_dir, 'tensile_tests',
                            'dog_bone', '2012-04-12_TT-12c-6cm-0-TU_SH4',
                            'ARAMIS', 'V1_kurz')
    
    input_list = Property(Array(float), depends_on = 'data_dir')
    @cached_property
    def _get_input_list(self):
        
        fn_list = os.listdir(self.data_dir)
        n_files = len(fn_list)
        
        fn_decomposed = [ name + '-' for name in fn_list[0].split('-')[:-1]]
        basename = ''.join(fn_decomposed)
        input_list = []
        min_shape = 1e10
        for fn_idx in range(0, n_files, self.time_step_size):
            fname = os.path.join(self.data_dir, '%s%d.txt' % (basename, fn_idx))
            print 'reading', fn_idx, '...'
            input_arr = np.loadtxt(fname,
                                   skiprows = 14,
                                   usecols = [0, 1, 2, 3, 4, 8, 9, 10])
            # handle the case with unrecognized facets
            # find the maximum row and node number.
            
            min_shape = min(input_arr.shape[0], min_shape)
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
        x_idx = np.array(data_arr[:, :, 0], dtype = int)
        y_idx = np.array(data_arr[:, :, 1], dtype = int)
        
        # find the shape of the point grid
        # by taking aramis indices stored in
        # first two columns
        x_min, x_max = np.min(x_idx, axis = 1), np.max(x_idx, axis = 1)
        y_min, y_max = np.min(y_idx, axis = 1), np.max(y_idx, axis = 1)
        
        n_x = np.max(x_max - x_min)
        n_y = np.max(y_max - y_min)
        
        n_t = data_arr.shape[0]
        t_idx = (np.arange(n_t)[:, np.newaxis] * 
                 np.ones((data_arr.shape[1],), dtype = int)[np.newaxis, :])

        # construct the mask for elements to be ignored
        grid_mask = np.zeros((n_t, n_x + 1, n_y + 1), dtype = bool)
        grid_mask[:, :, :] = True
        grid_mask[(t_idx.flatten(), x_idx.flatten(), y_idx.flatten())] = False
        return grid_mask, t_idx, x_idx, y_idx

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
        x_idx = np.array(data_arr[-1, :, 0], dtype = int)
        y_idx = np.array(data_arr[-1, :, 1], dtype = int)
        
        # find the shape of the point grid
        # by taking aramis indices stored in
        # first two columns
        x_min, x_max = np.min(x_idx), np.max(x_idx)
        y_min, y_max = np.min(y_idx), np.max(y_idx)

        x_idx -= x_min
        y_idx -= y_min
        
        n_x = x_max - x_min
        n_y = y_max - y_min
        
        n_t = data_arr.shape[0]
        # construct the mask for elements to be ignored
        grid_mask = np.zeros((n_t, n_x + 1, n_y + 1), dtype = bool)
        grid_mask[:, :, :] = True
        grid_mask[(slice(None), x_idx, y_idx)] = False
        return grid_mask, x_idx, y_idx

    grid_mask = Property
    def _get_grid_mask(self):
        return self.idx_maps[0]
    
    xy_idx = Property
    def _get_xy_idx(self):
        return self.idx_maps[1], self.idx_maps[2]
        
    data_t = Property()
    @cached_property
    def _get_data_t(self):
    
        data_arr = self.data_arr
        
        # total number of data values associated with a single point
        n_values = data_arr.shape[2] - 2
    
        # define a three dimensional array
        # first two indices are grid indices, third
        # index selects one of the values associated with a point
        daf = np.zeros(self.grid_mask.shape + (n_values,), dtype = float)
        
        # define an array of indexes performing the mapping
        # from the data_arr into daf
        select_idx_arr = (self.txy_idx[0], self.txy_idx[1], self.txy_idx[2], slice(None))

        daf[select_idx_arr] = data_arr[:, :, 2:] 

        #data_t = ma.masked_array(daf, mask = self.grid_mask)

        return daf

    #===========================================================================
    # Geometry arrays
    #===========================================================================
    x_arr = Property
    @cached_property    
    def _get_x_arr(self):
        return np.average(self.data_t[self.w_detect_step, :, :, 0], axis = 1)

    x_idx_arr = Property
    @cached_property    
    def _get_x_idx_arr(self):        
        return np.arange(self.x_arr.shape[0], dtype = int)

    y_arr = Property
    @cached_property    
    def _get_y_arr(self):
        return np.average(self.data_t[self.w_detect_step, :, :, 1], axis = 0)
    
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
        ux_masked = ma.masked_array(ux_arr, mask = self.grid_mask)

        ux_avg = ma.average(ux_masked, axis = 2)
        t_idx_zeros, x_idx_zeros, y_idx_zeros = np.where(self.grid_mask)
        ux_arr[t_idx_zeros, x_idx_zeros, y_idx_zeros] = ux_avg[t_idx_zeros, x_idx_zeros]
        
        return ux_arr
    
    d4_ux_t = Property
    @cached_property
    def _get_d4_ux_t(self):
        # generate the elem_node_map
        
        ir = self.integ_radius
        ux_arr = self.ux_t
        d4_ux_arr = np.zeros_like(ux_arr)
        d4_ux_arr[:, ir:-ir, :] = (ux_arr[:, 2 * ir:, :] - 
                                   ux_arr[:, :-2 * ir, :])
        
        # cutoff the negative strains - noise
        d4_ux_arr[ d4_ux_arr < 0.0 ] = 0.0
        
        return d4_ux_arr
    
    d4_ux_t_avg = Property
    @cached_property
    def _get_d4_ux_t_avg(self):
        return np.average(self.d4_ux_t, axis = 2)

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
    
    d4_ux_w = Property
    @cached_property
    def _get_d4_ux_w(self):
        d4_ux_w = get_d4(self.ux_w, self.integ_radius)
        # cutoff the negative strains - noise
        d4_ux_w[ d4_ux_w < 0.0 ] = 0.0
        return d4_ux_w
        
    d4_ux_avg_w = Property
    @cached_property
    def _get_d4_ux_avg_w(self):
        return np.average(self.d4_ux_w, axis = 1)
        
    dd44_ux_w = Property
    @cached_property
    def _get_dd44_ux_w(self):
        return get_d4(self.d4_ux_w, self.integ_radius)

    dd44_ux_avg_w = Property
    @cached_property
    def _get_dd44_ux_avg_w(self):
        return np.average(self.dd44_ux_w, axis = 1)
        
    ddd444_ux_w = Property
    @cached_property
    def _get_ddd444_ux_w(self):
        return get_d4(self.dd44_ux_w, self.integ_radius)

    ddd444_ux_avg_w = Property
    @cached_property
    def _get_ddd444_ux_avg_w(self):
        return np.average(self.ddd444_ux_w, axis = 1)
        
    #===========================================================================
    # Crack detection
    #===========================================================================
    crack_filter = Property
    @cached_property
    def _get_crack_filter(self):
        dd44_ux_w = self.dd44_ux_w 
        ddd444_ux_w = self.ddd444_ux_w
        return ((dd44_ux_w[1:, :] * dd44_ux_w[:-1, :] < 0.0) * 
                ((ddd444_ux_w[1:, :] + ddd444_ux_w[:-1, :]) / 2.0 < -0.005)) 
    
    crack_filter_avg = Property
    @cached_property
    def _get_crack_filter_avg(self):
        dd44_ux_avg_w = self.dd44_ux_avg_w
        ddd444_ux_avg_w = self.ddd444_ux_avg_w
        return ((dd44_ux_avg_w[1:] * dd44_ux_avg_w[:-1] < 0.0) * 
                ((ddd444_ux_avg_w[1:] + ddd444_ux_avg_w[:-1]) / 2.0 < -0.01))

    crack_arr = Property
    @cached_property
    def _get_crack_arr(self):
        return self.d4_ux_w[np.where(self.crack_filter)]

    crack_field_w = Property
    @cached_property
    def _get_crack_field_w(self):    

        cf_w = np.zeros_like(self.d4_ux_w)
        cf_w[np.where(self.crack_filter)] = self.crack_arr
        return cf_w

    crack_field_t = Property
    @cached_property
    def _get_crack_field_t(self):    
        '''
        '''
        crack_idx_arr = np.where(self.crack_filter[None, :, :])
        crack_arr_t = self.d4_ux_t[crack_idx_arr]
        cf_t = np.zeros_like(self.d4_ux_t)
        cf_t[np.where(self.crack_filter[None, :, :])] = crack_arr_t
        return cf_t

    #===========================================================================
    # Plotting methods
    #===========================================================================
    def plot(self):    
        p.subplot(2, 2, 1)
        p.plot(self.x_idx_arr[:-1], self.crack_filter_avg * 0.1,
               color = 'magenta', linewidth = 2)
        p.plot(self.x_idx_arr, self.d4_ux_avg_w, color = 'black')

        p.subplot(2, 2, 3)
        p.plot(self.x_idx_arr[:-1], self.crack_filter_avg * 0.1,
               color = 'magenta', linewidth = 2)
        p.plot(self.x_idx_arr, self.dd44_ux_avg_w, color = 'black')
        p.plot(self.x_idx_arr, self.ddd444_ux_avg_w, color = 'blue')
        
        p.subplot(2, 2, 2)
        p.plot(self.x_idx_arr, self.ux_w, color = 'green')
        p.plot(self.x_idx_arr, self.ux_w_avg, color = 'red')
        
        p.subplot(2, 2, 4)
        p.hist(self.crack_arr, bins = 40)
        
        p.show()
    
    def plot3d(self):
        #qargs = [ daf[:, :, i, np.newaxis] for i in range(0, 6)]
        
        #m.quiver3d(*qargs)

#        barr = np.zeros_like(self.w_t)
#        s = m.contour3d(self.t_idx_arr[:, np.newaxis, np.newaxis] + barr,
#                       self.x_idx_arr[np.newaxis, :, np.newaxis] + barr,
#                       self.y_idx_arr[np.newaxis, np.newaxis, :] + barr,
#                       self.d4_ux_t, contours = 30)
        
        ## print a surface - why is moved in the y direction?
#        m.surf(self.x_idx_arr, self.y_idx_arr, self.d4_ux_w)#, mask = mask_idx_array)
#        m.surf(self.x_arr, self.y_arr, self.d4_ux_w)#, mask = mask_idx_array)
        m.surf(self.x_arr, self.y_arr, self.crack_field_w)#, mask = mask_idx_array)
#        
        m.show()
        
        # @todo: 
        # sensitivity with respect to the integ_range
        # 
        
if __name__ == '__main__':
    ct = CrackTracer(time_step_size = 40,
                     integ_radius = 10,
                     w_detect_step = -1)
    ct.plot3d()
