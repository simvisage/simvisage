'''
MUSHROOF - For double symmetric loading cases, one Roof is sufficient 

TODO: @ Andreas
     - split dead load cases, include additional dead load case
'''

from traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool, HasTraits, Enum, Dict, Str

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic
from ibvpy.mats.mats3D.mats3D_tensor import \
    map3d_sig_eng_to_mtx

from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U
from ibvpy.mesh.fe_grid import \
    FEGrid

from ibvpy.mesh.fe_spring_array import \
     FESpringArray

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, shape, \
    cos, sin, arctan, where, abs, all, any, diag, \
    argsort, sum, transpose, vstack, append, arange, \
    hstack, reshape, shape, size, zeros_like, cumsum, \
    dsplit, copy, dtype, sort, ones_like, max, argmax

# Interpolation
from scipy.interpolate import Rbf

from math import \
    sqrt, asin, acos, pi as Pi

from matplotlib.pyplot import \
    bar, show, axhline, ion, ioff, xlabel, ylabel, title, figure, savefig

import csv


from matresdev.simiter.sim_pstudy import\
    ISimModel, SimOut, SimPStudy, SimArray, SimArrayView


from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

from geo_column import GEOColumn

#from hp_shell import HPShell
from hp_shell_conn_detail import HPShell

from geo_hp_shell import GeoHPShell

from mathkit.geo.transform.square2circle import GeoSquare2Circle

from mush_roof_model import MushRoofModel


class MRDetail(MushRoofModel):

    implements(ISimModel)

    # parameter used by hp_shell
    # if set to "one" hp shell uses symmetry of one roof to construct geometry
    # otherwise only one quarter is assumed.
    #
    mushroof_part = 'one'

    # shift of elements in hp_shell
    #
    shift_elems = True

    material_density_roof = Float(-22.4e-3)    # [MN/m^3]

    #----------------------------------------------------
    # elements
    #----------------------------------------------------

    n_elems_xy_quarter = Int(5, input = True)
    n_elems_z = Int(1, input = True)

    n_elems_col_z = Int(10 , input = True)
    n_elems_col_xy = Int(2 , input = True)

    vtk_r = Float(0.9, input = True)

    # fets used for roof
    #
    fe_roof = Instance((FETSEval), depends_on = '+ps_levels, +input')
    def _fe_roof_default(self):
        # fe_quad_serendipity_roof is defined in base class 
        # connected with material properties of the roof
        #
        #fets = self.fe_linear_roof
        fets = self.fe_quad_serendipity_roof
        fets.vtk_r *= self.vtk_r
        return fets

    # fets used for steel plate
    #
    fe_plate = Instance((FETSEval), depends_on = '+ps_levels, +input')
    def _fe_plate_default (self):
        #fets = self.fe_linear_plate
        fets = self.fe_quad_serendipity_plate
        # default integration scheme defined as 2x2x2
        # use higher order integration needed to get better 
        # results at connection points with the column
        #
        fets.ngp_r = 3
        fets.ngp_s = 3
        fets.ngp_t = 3
        fets.vtk_r *= self.vtk_r
        return fets

    #----------------------------------------------------
    # geometric transformation
    #----------------------------------------------------

    hp_shell = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_hp_shell(self):
        # list of local origins of each roof defined in global coordinates
        # (for MRDetail two roofs are considered due to symmetry)
        #
        X0 = [0, 0, 0]
        return HPShell(t_shell = self.t_shell,
                         X0 = X0,
                         delta_h_scalefactor = self.delta_h_scalefactor,
                         n_elems_xy = self.n_elems_xy,
                         n_elems_z = self.n_elems_z,
                         shift_elems = self.shift_elems,
                         const_edge_elem = self.const_edge_elem,
                         shift_array = self.shift_array,
                         mushroof_part = self.mushroof_part)

    shift_array = Array(value = [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1], ], input = True)

    plate = Property(Instance(GEOColumn) , depends_on = '+ps_levels, +input')
    @cached_property
    def _get_plate(self):
        # column transformation can be used to obtain the geometry of the plate
        #
        X0 = [ 0., 0., -self.t_plate ]
        return GEOColumn(width_top = self.width_top_col,
                          width_bottom = self.width_top_col,
                          X0 = X0,
                          h_col = self.t_plate)

    #----------------------------------------------------
    # fe-grids 
    #----------------------------------------------------

    fe_grid_roof = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_grid_roof(self):
        hp_shell = GeoHPShell(thickness = 0.6)

        square_to_circle = GeoSquare2Circle(post_transform = hp_shell,
                                             circle_radius = 0.15,
                                             #Scircle_center = [4.0, 4.0, 0.0],
                                             square_edge = 0.6)
        fe_grid = FEGrid(coord_min = (-1.0, -1.0, 0.0),
                          coord_max = (1.0, 1.0, 1.0),
                          geo_transform = square_to_circle,
                          shape = (self.n_elems_xy, self.n_elems_xy, self.n_elems_z),
                          fets_eval = self.fe_roof)

        mid_idx = self.mid_idx
        idx_min, idx_max = self.idx_min, self.idx_max
        print 'idx_min', idx_min
        print 'idx_max', idx_max
        interior_elems = fe_grid[ idx_min:idx_max, idx_min:idx_max, :, :, :, : ].elems
        fe_grid.inactive_elems = list(interior_elems)
        print 'elems', interior_elems
        return fe_grid

    n_hollow_elems = Int(1)

    mid_idx = Property()
    def _get_mid_idx(self):
        return self.n_elems_xy / 2

    idx_min = Property
    def _get_idx_min(self):
        return (self.mid_idx - self.n_hollow_elems)

    idx_max = Property
    def _get_idx_max(self):
        return (self.mid_idx + self.n_hollow_elems)

    bc_fix_dofs_in_the_hollow = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_fix_dofs_in_the_hollow(self):

        idx_min, idx_max = self.idx_min, self.idx_max
        roof = self.fe_grid_roof

        link_slice = roof[ idx_min, idx_min, 0, 0, 0, 0 ]

        s_list = [ roof[ idx_min:idx_max - 1, idx_min:idx_max - 1, :, 1:, 1:, : ],
                   roof[ idx_max - 1, idx_min:idx_max - 1, :, 1:-1, 1:, : ],
                   roof[ idx_min:idx_max - 1, idx_max - 1, :, 1:, 1:-1, : ],
                   roof[ idx_max - 1, idx_max - 1, :, 1:-1, 1:-1, : ] ]

        return [ BCSlice(var = 'u'  , dims = [0, 1, 2],
                           slice = s,
                           link_slice = link_slice,
                           link_coeffs = [0],
                           value = 1.0)
                           for s in s_list ]

    fe_grid_plate = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_grid_plate(self):
        return FEGrid(coord_min = (0.0, 0.0, 0.0),
                       coord_max = (1.0, 1.0, 1.0),
                       geo_transform = self.plate,
                       shape = (self.n_elems_col_xy, self.n_elems_col_xy, 2),
                       fets_eval = self.fe_plate)

    #----------------------------------------------------
    # ps_study
    #----------------------------------------------------
    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''

        # DISPLACEMENT
        #
        U = self.tloop.eval()

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [  SimOut(name = '$u_z$', unit = '[mm]'), ]


    #----------------------------------------------------
    # response tracer
    #----------------------------------------------------

    rtrace_list = List
    def _rtrace_list_default(self):
        return [  self.max_princ_stress, self.sig_app, self.u]

    #----------------------------------------------------
    # time loop
    #----------------------------------------------------

    tloop = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_tloop(self):

        roof = self.fe_grid_roof

                  # self.bc_roof_deadweight + \
        bc_list = self.bc_roof_list + self.bc_outer_boundary + \
                  self.bc_fix_dofs_in_the_hollow

        ts = TS(sdomain = [ roof ],
                 dof_resultants = True,
                 bcond_list = bc_list,
                 rtrace_list = self.rtrace_list

               )

        # Add the time-loop control
        #

        tloop = TLoop(tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline)

        return tloop

    #----------------------------------------------------
    # boundaries
    #----------------------------------------------------

    bc_roof_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_roof_list(self):
        '''
        links all plate corner nodes of each elements to the adjacent elements of the roof
        '''

        bc_roof_list = []
        n_from_center = 3

        roof, plate = self.fe_grid_roof, self.fe_grid_plate
        slice_1 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter - n_from_center ,
                                           self.n_elems_xy_quarter, 0,
                                           0, 0, 0 ],
#                              link_slice = plate[ 0 , 0 , -1, 0, 0, -1], link_coeffs = [1.0],
                              value = 0.)]
        slice_2 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter,
                                           self.n_elems_xy_quarter - n_from_center, 0,
                                           0, 0, 0 ],
#                              link_slice = plate[ -1, 0, -1, -1, 0, -1], link_coeffs = [1.0],
                              value = 0.)]

        slice_3 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter + n_from_center,
                                           self.n_elems_xy_quarter, 0,
                                           0, 0, 0 ],
#                              link_slice = plate[ -1 , -1 , -1, -1, -1, -1], link_coeffs = [1.0],
                              value = 0.)]

        slice_4 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter ,
                                           self.n_elems_xy_quarter + n_from_center, 0,
                                           0, 0, 0 ],
#                              link_slice = plate[ 0 , -1 , -1, 0, -1, -1], link_coeffs = [1.0],
                              value = 0.)]
        bc_roof_list = bc_roof_list + \
                          slice_1 + slice_2 + slice_3 + slice_4

        return bc_roof_list

    bc_plate_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_plate_list(self):
        '''
        links all plate corner nodes of each elements to the adjacent elements of the roof
        '''

        plate = self.fe_grid_plate
        slice_1 = BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = plate[ : , : , 0, :, :, 0],
                              value = 0.)

        return [ slice_1 ]

    bc_roof_deadweight = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_roof_deadweight(self):
        '''
        links all plate corner nodes of each elements to the adjacent elements of the roof
        '''

        roof = self.fe_grid_roof
        slice_1 = BCSlice(var = 'f', value = self.material_density_roof, dims = [2],
                          integ_domain = 'global',
                          slice = roof[:, :, :, :, :, :])

        return [ slice_1 ]

    bc_outer_boundary = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_outer_boundary(self):
        '''
        links all plate corner nodes of each elements to the adjacent elements of the roof
        '''

        roof = self.fe_grid_roof
        slice_1 = BCSlice(var = 'u', value = -0.01, dims = [2],
                          integ_domain = 'global',
                          slice = roof[0, ::2, 0, 0, -1, 0])
        slice_2 = BCSlice(var = 'u', value = -0.01, dims = [2],
                          integ_domain = 'global',
                          slice = roof[0, ::2, -1, 0, -1, -1])
        slice_3 = BCSlice(var = 'u', value = -0.01, dims = [2],
                          integ_domain = 'global',
                          slice = roof[0, 0, 0, 0, -1, 0])
        slice_4 = BCSlice(var = 'u', value = -0.01, dims = [2],
                          integ_domain = 'global',
                          slice = roof[0, 0, -1, 0, -1, -1])

        return [ slice_1, slice_2, slice_3, slice_4 ]

if __name__ == '__main__':

    sim_model = MRDetail(n_elems_z = 2,
                          n_elems_xy_quarter = 11
                           )

    do = 'ui'

    if do == 'ui':
        print '*** ui ***'

        sim_model.lc = 'lc_w_asym'

        # element of shif_array that should not be linked
        # (corresponds to 'equal_100cm')
        #
#        sim_model.not_linked_elem = 0

        # evaluation
        #
        sim_model.peval()

        # visualisation
        # 
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource = sim_model)
        app.main()


