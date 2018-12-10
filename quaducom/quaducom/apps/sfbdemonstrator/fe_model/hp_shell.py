'''
Created on Jun 16, 2010

@author: alexander
'''


from etsproxy.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, Enum, \
    Dict, Bool, Int, List

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, mgrid, arange, \
    where, sum, sin, cos, vstack, hstack, argmax, newaxis, size, \
    shape, sqrt, frompyfunc, ones_like, zeros_like, ones, any, all, \
    sort, argsort, concatenate, add

from .rsurface_reader import \
    read_rsurface, normalize_rsurfaces

from os.path import join

from math import pi
from matplotlib.pyplot import *

# Interpolation
from scipy.interpolate import Rbf

def delete_second_rows(arr, nx = 20, ny = 20):
    '''Remove the second and second to last column and row from the regular grid of points.
    for lowerface_4x4m.robj and upperface_4x4m.robj data file
    '''
    va = arr
    va = va.reshape((nx, ny, 3))
    va_a = vstack([ va[ :1, :, : ], va[ 2:-2, :, : ], va[-1:, :, : ] ])
    va_b = hstack([ va_a[ :, :1, : ], va_a[ :, 2:-2, : ], va_a[ :, -1:, : ] ])
    n_size = va_b.shape[0] * va_b.shape[1]
    va = va_b.reshape((n_size, 3))
    return va

class HPShell(HasTraits):
    '''Geometry definition of a hyperbolic parabolid shell.
    '''

    #-----------------------------------------------------------------
    # geometric parameters of the shell
    #-----------------------------------------------------------------

    # dimensions of the shell for one quarter of mush_roof
    #
    length_xy_quarter = Float(3.5, input = True) # [m]
    length_z = Float(0.927, input = True) # [m]

    # corresponds to the delta in the geometry .obj-file with name '4x4m' as a cut off
    #
    delta_h = Float(0.865, input = True) # [m]

    # scale factors for geometric dimensions
    # NOTE: parameters can be scaled separately, i.e. scaling of 'delta_h' (inclination of the shell)
    # does not effect the scaling of the thickness
    #
    scalefactor_delta_h = Float(1.00, input = True) # [-] 
    scalefactor_length_xy = Float(1.00, input = True) # [-] 

    # thickness of the shell
    # NOTE: only used by the option 'const_reinf_layer'
    #
    t_shell = Float(0.06, input = True) # [m]

    width_top_col = Float(0.45, input = True) # [m]

    # factor shifting the z coordinates at the middle nodes of the edges upwards
    # in order to include the effect of imperfection in the formwork
    z_imperfection_factor = Float(0.0)

    #-----------------------------------------------------------------
    # specify the relation of the total structure (in the 'mushroof'-model)
    # with respect to a quarter of one shell defined in 'HPShell'
    #-----------------------------------------------------------------

    # @todo: "four" is not supported by "n_elems_xy_dict"in mushroff_model
    mushroof_part = Enum('one', 'quarter', 'four', input = True)

    # 'scale_size' parameter is used as scale factor for different mushroof parts
    # Defines the proportion between the length of the total model
    # with respect to the length of a quarter shell as the 
    # basic substructure of which the model consists of.
    # @todo: add "depends_on" or remove "cached_property"
    # @todo: move to class definition of "mushroof_model" not in "HPShell"
    #        (also see there "n_elems_dict" with implicit "scale_factor")
    #
    scale_size = Property(Float, depends_on = 'mushroof_part')
    @cached_property
    def _get_scale_size(self):
        scale_dict = {'quarter': 1.0 , 'one' : 2.0, 'four' : 4.0 }
        return scale_dict[ self.mushroof_part ]

    # origin of the shell
    #
    X0 = Array(float, input = True)
    def _X0_default(self):
        return array([ 0., 0., 0. ])


    #-----------------------------------------------------------------
    # discretisation
    #-----------------------------------------------------------------

    # number of element used for the discretisation ( dimensions of the entire model)
    #
    n_elems_xy_quarter = Int(5, input = True)
    n_elems_z = Int(3, input = True)

    n_elems_xy = Property(Int, depends_on = 'n_elems_xy_quarter, +input')
    @cached_property
    def _get_n_elems_xy(self):
        return  self.n_elems_xy_quarter * self.scale_size


    #-----------------------------------------------------------------
    # option: 'shift_elems'
    #-----------------------------------------------------------------

    # shift of column elements
    # if set to "True" (default) the information defined in 'shift_array' is used.
    #
    shift_elems = Bool(True, input = True)

    # 'shift_array' is used to place element corners at a defined global 
    # position in order to connect the shell with the corner nodes of the column. 
    # [x_shift, y_shift, number of element s between the coordinate position]
    # NOTE: 'shift_array' needs to have shape (:,3)!
    #
    shift_array = Array (float, input = True)
    def _shift_array_default(self):
        return array([[self.width_top_col / 2 ** 0.5, self.width_top_col / 2 ** 0.5, 1]])


    #-----------------------------------------------------------------
    # option: 'const_reinf_layer'
    #-----------------------------------------------------------------
    # 'const_reinf_layer' - parameter is used only for the non-linear analysis,
    # where an element layer with a constant thickness is needed to simulate the 
    # reinforced concrete at the top and bottom of the TRC-shell.
    #
    const_reinf_layer_elem = Bool(False, input = True)
    t_reinf_layer = Float (0.03, input = True) # [m]
    n_elems_reinf_layer = Int(1, input = True) #number of dofs used for edge refinement  


    #-----------------------------------------------------------------
    # read vertice points of the shell and derive a normalized 
    # RBF-function for the shell approximation
    #-----------------------------------------------------------------

    # "lowerface_cut_off" - option replaces constant height for the coordinates 
    # which connect to the column (this cuts of the shell geometry horizontally 
    # at the bottom of the lower face of the shell geometry.
    #
    cut_off_lowerface = Bool(True, input = True)

    # choose geometric file (obj-data file)
    #
    geo_input_name = Enum('350x350cm', '4x4m', '02', input = True)

    # filter for '4x4m' file needs to be done to have regular grid
    # in order to rbf-function leading to stable solution without oscilation
    #
    geo_filter = Dict({'4x4m' : delete_second_rows })
#                        ,'350x350cm' : delete_second_rows} )

    def _read_arr(self, side = 'lowerface_'):
        '''read the robj-file saved in the subdirectory 
        'geometry_files'
        '''
        file_name = side + self.geo_input_name + '.robj'
        file_path = join('geometry_files', file_name)
        # get an array with the vertice coordinates
        #
        v_arr = read_rsurface(file_path)
#        print 'v_arr before filtering \n', v_arr
#        print 'v_arr.shape before filtering \n', v_arr.shape
        filter = self.geo_filter.get(self.geo_input_name, None)
        if filter != None:
            v_arr = list(filter(v_arr))
#        print 'v_arr after filtering \n', v_arr
#        print 'v_arr.shape after filtering \n', v_arr.shape
        return v_arr

    # array of the vertex positions in global 
    # x,y,z-coordinates defining the lower surface of the shell 
    #
    vl_arr = Property(Array(float), depends_on = 'geo_input_name')
    @cached_property
    def _get_vl_arr(self):
        vl_arr = self._read_arr('lowerface_')
        if self.cut_off_lowerface == True:
            print('--- lower face z-coords cut off ---')

            # z-values of the coords from the lower face are cut off. 
            # From the highest z-coordinate of the lower face the vertical
            # distance is 'delta h (for 4x4m: delta_h = 1.0m).
            # At this limit the lower face is cut off. 
            # NOTE: the global z coordinate is assumed to point up 
            # and must be given in the same unite as 'delta_h', i.e. in [m].
            #
            vl_z_max = max(vl_arr[:, 2])
            if self.geo_input_name == '4x4m':
                # NOTE: the global z-coordinates are given in the geo data file in [m]
                # no conversion of unites necessary (self.delta_h is given in [m])
                delta_h = self.delta_h
            elif self.geo_input_name == '350x350cm':
                # NOTE: the global z-coordinates are given in the geo data file in [cm]
                # convert delta_h from [m] to [cm]
                #
                delta_h = self.delta_h * 100.
            vl_z_min = vl_z_max - delta_h
            vl_arr_z = where(vl_arr[:, 2] < vl_z_min, vl_z_min, vl_arr[:, 2])
            vl_arr = c_[vl_arr[:, 0:2], vl_arr_z]

        return vl_arr

    # array of the vertex positions in global 
    # x,y,z-coordinates defining the upper surface of the shell
    #
    vu_arr = Property(Array(float), depends_on = 'geo_input_name')
    @cached_property
    def _get_vu_arr(self):
        return self._read_arr('upperface_')

    # normalized coordinates of the vertices for the lowerface
    # NOTE: the underline character indicates a normalized value 
    # @todo: 'normalize_rsurfaces' is called twice for 'vl_arr_' and 'vu_arr_' 
    #
    vl_arr_ = Property(Array(float), depends_on = 'geo_input_name')
    @cached_property
    def _get_vl_arr_(self):
        vl_arr_, vu_arr_ = normalize_rsurfaces(self.vl_arr, self.vu_arr)
        return vl_arr_

    # normalized coordinates of the vertices for the lowerface
    # NOTE: the underline character indicates a normalized value 
    #
    vu_arr_ = Property(Array(float), depends_on = 'geo_input_name')
    @cached_property
    def _get_vu_arr_(self):
        vl_arr_, vu_arr_ = normalize_rsurfaces(self.vl_arr, self.vu_arr)
        return vu_arr_


    rbf_l_ = Property(Instance(Rbf), depends_on = 'geo_input_name')
    @cached_property
    def _get_rbf_l_(self):
        # use a radial basis function approximation (rbf) (i.e. interpolation of
        # scattered data) based on the normalized vertex points of the lower face
        #
        xl_ = self.vl_arr_[:, 0]
        yl_ = self.vl_arr_[:, 1]
        zl_ = self.vl_arr_[:, 2]

        # flip the orientation of the local coordinate axis 
        # depending on the geometry file used
        #
        if self.geo_input_name == '350x350cm':
            xl_ = 1 - self.vl_arr_[:, 0]
        if self.geo_input_name == '4x4m':
            yl_ = 1 - self.vl_arr_[:, 1]

        rbf_l_ = Rbf(xl_, yl_, zl_, function = 'cubic')
#        rbf_l_ = Rbf( xl_, yl_, zl_, function = 'linear' )
        return rbf_l_


    rbf_u_ = Property(Instance(Rbf), depends_on = 'geo_input_name')
    @cached_property
    def _get_rbf_u_(self):
        # use a radial basis function approximation (rbf) (i.e. interpolation of
        # scattered data) based on the normalized vertex points of the upper face
        #
        xu_ = self.vu_arr_[:, 0]
        yu_ = self.vu_arr_[:, 1]
        zu_ = self.vu_arr_[:, 2]

        # flip the orientation of the local coordinate axis 
        # depending on the geometry file used
        #
        if self.geo_input_name == '350x350cm':
            xu_ = 1 - self.vu_arr_[:, 0]
        if self.geo_input_name == '4x4m':
            yu_ = 1 - self.vu_arr_[:, 1]

        rbf_u_ = Rbf(xu_, yu_, zu_, function = 'cubic')
#        rbf_u_ = Rbf( xu_, yu_, zu_, function = 'linear' )
        return rbf_u_


    #------------------------------------------------------------------------------ 
    # hp_shell geometric transformation
    # NOTE: returns the global coordinates of the shell based on the supplied local
    #       grid points  
    #------------------------------------------------------------------------------ 

    def __call__(self, points):
        '''Return the global coordinates of the supplied local points.
        '''

        # number of local grid points for each coordinate direction
        # NOTE: values must range between 0 and 1
        #
        xi_, yi_, zi_ = points[:, 0], points[:, 1], points[:, 2]

        # insert imperfection (shift the middle node of the shell upwards)
        imp = self.z_imperfection_factor
        zi_ += imp * xi_ + imp * yi_ - 2 * imp * xi_ * yi_

        # size of total structure
        #
        # @todo: move to class definition of "mushroof_model" and send to "__call__"
        scale_size = self.scale_size
        # @todo: add "_quarter" (see above)
        length_xy_tot = self.length_xy_quarter * scale_size
        n_elems_xy_quarter = self.n_elems_xy_quarter
#        print 'HPShell n_elems_xy_quarter', n_elems_xy_quarter
        # distance from origin for each mushroof_part
        #
        def d_origin_fn(self, coords):
            if self.mushroof_part == 'quarter':
                return coords
            if self.mushroof_part == 'one':
                return abs(2.0 * coords - 1.0)
            # @todo: corresponding "scale_factor" needs to be added 
            #        in order for this to work
            if self.mushroof_part == 'four':
                return  where(coords < 0.5, abs(4 * coords - 1), abs(-4 * coords + 3))

        # element at column shift
        #
        if self.shift_elems == True:

            # define the origin for each model part
            #
            def origin_fn(self, coords):
                if self.mushroof_part == 'quarter':
                    return zeros_like(coords)
                if self.mushroof_part == 'one':
                    return ones_like(xi_) * 0.5
                if self.mushroof_part == 'four':
                    return where(coords < 0.5, 0.25, 0.75)

            def piecewise_linear_fn(x, x_fix_arr_, y_fix_arr_):
                '''creates a piecewise linear_fn going through the fix_points
                values need to be normed running between 0..1
                and values have to be unique'''
                x_fix_arr_ = hstack((0, x_fix_arr_, 1))
                y_fix_arr_ = hstack((0, y_fix_arr_, 1))
                rbf_fn_ = Rbf(x_fix_arr_, y_fix_arr_, function = 'linear') #rbf has to be linear
                return rbf_fn_(x)

            # define origin for quarter 
            #
            xi_origin_arr_ = origin_fn(self, xi_)
            yi_origin_arr_ = origin_fn(self, yi_)
#            print 'xi_origin_arr_', xi_origin_arr_

            # delta towards origin
            #
            xi_delta_arr_ = (xi_ - xi_origin_arr_) * scale_size
            yi_delta_arr_ = (yi_ - yi_origin_arr_) * scale_size
#            print 'xi_delta_arr_', xi_delta_arr

            # define sign  
            #
            xi_sign_arr = where(xi_delta_arr_ == 0., 0., xi_delta_arr_ / abs(xi_delta_arr_))
            yi_sign_arr = where(yi_delta_arr_ == 0., 0., yi_delta_arr_ / abs(yi_delta_arr_))
#            print 'xi_sign_arr', xi_sign_arr

            # fix points defined in shift array as normelized values
            #
            x_fix_ = self.shift_array[:, 0] / self.length_xy_quarter
#            print 'x_fix_', x_fix_

            y_fix_ = self.shift_array[:, 1] / self.length_xy_quarter
            n_fix_ = add.accumulate(self.shift_array[:, 2]) / n_elems_xy_quarter

#            print 'add.accumulate( self.shift_array[:, 2] )', add.accumulate( self.shift_array[:, 2] )
#            print 'n_fix_', n_fix_
#            print 'piecewise_linear_fn', piecewise_linear_fn( abs( xi_delta_arr_ ),
#                                                                   n_fix_,
#                                                                   x_fix_ ) / scale_size

            # new xi_
            #
            xi_ = xi_origin_arr_ + xi_sign_arr * piecewise_linear_fn(abs(xi_delta_arr_),
                                                                    n_fix_,
                                                                    x_fix_) / scale_size

#            print 'xi_new', xi_

            # new yi
            #
            yi_ = yi_origin_arr_ + yi_sign_arr * piecewise_linear_fn(abs(yi_delta_arr_),
                                                                    n_fix_,
                                                                    y_fix_) / scale_size

            #-------------------------------------------------------------------------------------


        # values are used to calculate the z-coordinate using RBF-function of the quarter
        # (= values of the distance to the origin as absolute value)
        #
        xi_rbf_ = d_origin_fn(self, xi_)
#        print 'xi_rbf_', xi_rbf_
        yi_rbf_ = d_origin_fn(self, yi_)

        # get the z-value at the supplied local grid points
        # of the lower face
        #
        zi_lower_ = self.rbf_l_(xi_rbf_, yi_rbf_)

        # get the z-value at the supplied local grid points
        # of the upper face
        # 
        zi_upper_ = self.rbf_u_(xi_rbf_, yi_rbf_)

        # constant edge element transformation
        #
        if self.const_reinf_layer_elem == True:
            # arrange and check data
            #
            if self.t_reinf_layer > self.t_shell / 2. or self.n_elems_z < 3:
                print('--- constant edge element transformation canceled ---')
                print('the following condition needs to be fullfilled: \n')
                print('self.t_reinf_layer <= self.t_shell/2 and self.n_elems_z >= 3')
            else:
                n_elems_z = float(self.n_elems_z)
                # normed thickness will evaluate as t_reinf_layer at each element
                t_reinf_layer_ = self.t_reinf_layer / self.length_z / (zi_upper_ - zi_lower_)

                # zi_old set off from top which needs to be shifted
                delta_ = self.n_elems_reinf_layer / n_elems_z

                # get upper, lower and internal coordinates, that need to be shifted
                zi_lower = where(zi_ <= delta_)
                zi_upper = where(abs(1 - zi_) <= delta_ + 1e-10)
                zi_inter = where(abs(zi_ - 0.5) < 0.5 - (delta_ + 1e-10))

                # narrowing of coordinates
                zi_[zi_lower] = zi_[zi_lower] * t_reinf_layer_[zi_lower] / delta_
                zi_[zi_upper] = 1 - (1 - zi_[zi_upper]) * t_reinf_layer_[zi_upper] / delta_
                zi_[zi_inter] = t_reinf_layer_[zi_inter] + \
                                (zi_[zi_inter] - delta_) / (1 - 2 * delta_)\
                                 * (1 - 2 * t_reinf_layer_[zi_inter])
                print('--- constant edge elements transformation done ---')

        # thickness is multiplied by the supplied zi coordinate
        #
        z_ = (zi_lower_ * self.scalefactor_delta_h + (zi_upper_ - zi_lower_) * zi_)

        # coordinates of origin
        #
        X_0, Y_0, Z_0 = self.X0

        print('--- geometric transformation done ---')

        # multiply the local grid points with the real dimensions in order to obtain the 
        # global coordinates of the mushroof_part:
        #
        return c_[ X_0 + (xi_ * length_xy_tot) * self.scalefactor_length_xy,
                   Y_0 + (yi_ * length_xy_tot) * self.scalefactor_length_xy,
                   Z_0 + z_ * self.length_z ]





if __name__ == '__main__':

    from numpy import mgrid, c_, hstack, vstack, shape
    from etsproxy.mayavi import mlab

    hp = HPShell()

    # discretization
    #
    hp.n_elems_xy_quarter = 10
    hp.n_elems_z = 1

    hp.shift_elems = True
    hp.cut_off_lowerface = True
    hp.geo_input_name = '350x350cm'
    hp.mushroof_part = 'one'

#    # scale ( 8m x 8m x 1m ) - geometry to ( 7m x 7m x 0.85m ) - geometry
#    # without changing thickness (remains 6 cm)
#    #
#    scalefactor_delta_h = Float( 0.85 ) # [-] 
#    scalefactor_length_xy = Float( 7. / 8 ) # [-] 

#    # used only for constant element thickness
#    # for the reinforcement layer in a non-linear calculation
#    #
#    hp.const_reinf_layer_elem = True
#    hp.n_elems_reinf_layer = 1
#    hp.t_reinf_layer = 0.01


    X, Y, Z = mgrid[0:1:complex(0, hp.n_elems_xy + 1),
                    0:1:complex(0, hp.n_elems_xy + 1),
                    0:1:complex(0, hp.n_elems_z + 1)]

    gpoints = c_[ X.flatten(), Y.flatten(), Z.flatten() ]

    fp1 = hp(gpoints)
    #print shape( fp1 )
    #print fp1

    mlab.points3d(fp1[:, 0], fp1[:, 1], fp1[:, 2],
                   scale_factor = 0.05 ,
#                   mode = "cube",
                   resolution = 8)

    mlab.show()
