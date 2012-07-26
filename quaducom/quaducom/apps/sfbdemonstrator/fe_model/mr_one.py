
#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Dec 9, 2009 by: rch

'''
MUSHROOF - Demostrator SFB532

@todos 
1) geo_transform for the bottom bar [Done]

'''

from etsproxy.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool, HasTraits, Enum

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel


from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic


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
from mathkit.mfn import MFnLineArray

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, shape, \
    cos, sin, arctan, where, abs, all, any, diag, vstack, transpose, \
    ones_like, size, arange, append, min, max, hstack

from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos, pi as Pi

from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

# Interpolation
from scipy.interpolate import Rbf

from geo_column import GEOColumn

from simiter.sim_pstudy import\
    ISimModel, SimOut, SimPStudy, SimArray, SimArrayView

from hp_shell import HPShell

from mush_roof_model import MushRoofModel

from time import time

from matplotlib.pyplot import bar, show, axhline

import csv

class MRone(MushRoofModel):

    implements(ISimModel)
    mushroof_part = 'one'
    #===============================================================================
    # fe_grid
    #===============================================================================

    n_elems_xy_quarter = Int(10, input = True)#, ps_levels = [4, 16, 5] )
    n_elems_z = Int(1, input = True)#, ps_levels = [1, 2, 1] )
    n_elems_col_z = Int(10  , input = True, ps_levels = [5, 20, 3 ])
    n_elems_col_xy = Int(2 , input = True, ps_levels = [2, 4, 1])

    shift_elems = True

    vtk_r = Float(0.90)

    #default roof
    fe_roof = Instance((FETSEval), depends_on = '+ps_levels, +input')
    def _fe_roof_default(self):
        fets = self.fe_quad_serendipity_roof
        fets.vtk_r *= 0.9
        return fets

    #default plate
    fe_plate = Instance((FETSEval), depends_on = '+ps_levels, +input')
    def _fe_plate_default (self):
        fets = self.fe_quad_serendipity_plate
        fets.ngp_r = 3
        fets.ngp_s = 3
        fets.ngp_t = 3
        fets.vtk_r *= 0.9
        return fets


    # shell
    #
    hp_shell = Property(Instance(HPShell) , depends_on = '+ps_levels, +input')
    @cached_property
    def _get_hp_shell(self):
        return HPShell(length_xy_quarter = self.length_xy_quarter,
                        length_z = self.length_z,
                        n_elems_xy_quarter = self.n_elems_xy_quarter,
                        n_elems_z = self.n_elems_z,
                        scalefactor_delta_h = self.scalefactor_delta_h,
                        const_reinf_layer_elem = self.const_reinf_layer_elem,
                        width_top_col = self.width_top_col,
                        mushroof_part = self.mushroof_part,
                        shift_array = self.shift_array,
                        X0 = self.X0)

    # plate
    #
    plate = Property(Instance(GEOColumn) , depends_on = '+ps_levels, +input')
    @cached_property
    def _get_plate(self):
        return GEOColumn(width_top = self.width_top_col,
                          width_bottom = self.width_top_col,
                          X0 = [ 3.5, 3.5, -self.t_plate ], # - 0.25],
                          h_col = self.t_plate)

    # column
    #
#    X0_column = Array( [ 4., 4., -3.] )
    column = Property(Instance(GEOColumn) , depends_on = '+ps_levels, +input')
    @cached_property
    def _get_column(self):
        return GEOColumn(width_top_col = self.width_top_col,
                          width_bottom_col = self.width_bottom_col,
                          h_col = self.h_col - self.t_plate,
                          #r_pipe = self.r_pipe,
                          X0 = [ 3.5, 3.5, -(self.h_col) ]) # - 0.5] )


    #default column
    fe_column = Instance((FETSEval), transient = True , depends_on = '+ps_levels, +input')
    def _fe_column_default(self):
        fets = self.fe_quad_serendipity_column
        fets.vtk_r *= 0.9
        return fets


    fe_grid_roof = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_grid_roof(self):
        return FEGrid(coord_min = (0.0, 0.0, 0.0),
                       coord_max = (1.0, 1.0, 1.0),
                       geo_transform = self.hp_shell,
                       shift_array = self.shift_array,
                       shape = (self.n_elems_xy, self.n_elems_xy, self.n_elems_z),
                       fets_eval = self.fe_roof)

    fe_grid_column = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_grid_column(self):
        return  FEGrid(coord_min = (0.0, 0.0, 0.0),
                        coord_max = (1.0, 1.0, 1.0),
                        geo_transform = self.column,
                        shape = (self.n_elems_col_xy, self.n_elems_col_xy, self.n_elems_col_z),
                        fets_eval = self.fe_column)

    fe_grid_plate = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_grid_plate(self):
        return  FEGrid(coord_min = (0.0, 0.0, 0.0),
                        coord_max = (1.0, 1.0, 1.0),
                        geo_transform = self.plate,
                        shape = (self.n_elems_col_xy, self.n_elems_col_xy, 2),
                        fets_eval = self.fe_plate)

    #===============================================================================
    # ps_study
    #===============================================================================
    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()

        U_edge = U[self.edge_corner_1_dof][0, 0, 2]

        F_int = self.tloop.tstepper.F_int

        F_int_slice_x = F_int[self.edge_roof_right]
        F_int_slice_y = F_int[self.edge_roof_top]

        #bring dofs into right order for plot
        #
        F_hinge_in_order_x = self.sort_by_dofs(self.edge_roof_top, F_int_slice_x)
        F_hinge_in_order_y = self.sort_by_dofs(self.edge_roof_top, F_int_slice_y)
        F_hinge_x = append(F_hinge_in_order_x[:, :-1, 0], F_hinge_in_order_x[-1, -1, 0])
        F_hinge_y = append(F_hinge_in_order_y[:, :-1, 1], F_hinge_in_order_y[-1, -1, 1])
        F_hinge_y_sum = sum(F_hinge_y.flatten())
        F_hinge_x_sum = sum(F_hinge_x.flatten())
#
#        self.visual_force_bar( F_hinge_x.flatten()
#                               , y_label = "internal force x [MN]"
#                               , Title = 'F_Hinge_x_shrinkage' )
#        self.visual_force_bar( F_hinge_y.flatten()
#                               , y_label = "internal force y [MN]"
#                               , Title = 'F_Hinge_y_shrinkage' )
        print "u_edge", U_edge
        print "n_elems_xy_col", self.n_elems_col_xy
        print "n_elems_z_col", self.n_elems_col_z
        print "n_elems_xy_quarter", self.n_elems_xy_quarter
        print "n_elems_z", self.n_elems_z



        return array([ U_edge,
#                        u_x_corner2,
#                       F_hinge_y_sum] )
#                        u_z_corner2,
#                        max_princ_stress ]
                       ], dtype = 'float_')

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name = 'U', unit = 'm'),
#                 SimOut( name = 'u_x_corner2', unit = 'm' ),
#                 SimOut( name = 'N Gelenk', unit = 'MN' ), ]
#                 SimOut( name = 'u_z_corner2', unit = 'm' ),
#                 SimOut( name = 'maximum principle stress', unit = 'MPa' )
                 ]

    #===============================================================================
    # response tracer
    #===============================================================================

    rtrace_list = List
    def _rtrace_list_default(self):
        return [  self.max_princ_stress, self.sig_app, self.u, self.f_dof]


    shift_array = Array(value = [[0.45 / 2 ** 0.5, 0.45 / 2 ** 0.5, 1], ], input = True)


    #===============================================================================
    # boundary conditions
    #===============================================================================

    bc_plate_roof_link_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_plate_roof_link_list(self):
        '''
        links all plate corner nodes of each elements to the adjacent elements of the roof
        '''
        roof = self.fe_grid_roof
        plate = self.fe_grid_plate
        bc_col_link_list = []

        slice_1 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter - 1 ,
                                           self.n_elems_xy_quarter, 0,
                                           0, 0, 0 ],
                              link_slice = plate[ 0 , 0 , -1, 0, 0, -1], link_coeffs = [1.0],
                              value = 0.)]
        slice_2 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter,
                                           self.n_elems_xy_quarter - 1, 0,
                                           0, 0, 0 ],
                              link_slice = plate[ -1, 0, -1, -1, 0, -1], link_coeffs = [1.0],
                              value = 0.)]

        slice_3 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter + 1,
                                           self.n_elems_xy_quarter, 0,
                                           0, 0, 0 ],
                              link_slice = plate[ -1 , -1 , -1, -1, -1, -1], link_coeffs = [1.0],
                              value = 0.)]

        slice_4 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                            slice = roof[self.n_elems_xy_quarter ,
                                         self.n_elems_xy_quarter + 1, 0,
                                         0, 0, 0 ],
                            link_slice = plate[ 0 , -1 , -1, 0, -1, -1], link_coeffs = [1.0],
                            value = 0.)]

        slice_5 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                            slice = roof[self.n_elems_xy_quarter ,
                                         self.n_elems_xy_quarter , 0,
                                         0, 0, 0 ],
                            link_slice = plate[ self.n_elems_col_xy / 2.0 , self.n_elems_col_xy / 2.0 , -1,
                                                0, 0, -1], link_coeffs = [1.0],
                            value = 0.)]


        bc_plate_roof_link_list = slice_1 + slice_2 + slice_3 + slice_4 + slice_5

        return bc_plate_roof_link_list


    bc_roof_top_roof_low_link_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_roof_top_roof_low_link_list(self):
        '''
        links all plate corner nodes of each elements to the adjacent elements of the roof
        '''
        roof = self.fe_grid_roof
        plate = self.fe_grid_plate
        bc_roof_top_roof_low_link_list = []

        slice_1 = [BCSlice(var = 'u'  , dims = [ 2],
                              link_slice = roof[self.n_elems_xy_quarter - 1 ,
                                           self.n_elems_xy_quarter, 0,
                                           0, 0, 0 ],
                              slice = roof[self.n_elems_xy_quarter - 1 ,
                                                self.n_elems_xy_quarter, -1,
                                                0, 0, -1 ], link_coeffs = [1.0],
                              value = 0.)]
        slice_2 = [BCSlice(var = 'u'  , dims = [ 2],
                              link_slice = roof[self.n_elems_xy_quarter,
                                           self.n_elems_xy_quarter - 1, 0,
                                           0, 0, 0 ],
                              slice = roof[self.n_elems_xy_quarter ,
                                                self.n_elems_xy_quarter - 1, -1,
                                                0, 0, -1 ], link_coeffs = [1.0],
                              value = 0.)]

        slice_3 = [BCSlice(var = 'u'  , dims = [ 2],
                              link_slice = roof[self.n_elems_xy_quarter + 1,
                                           self.n_elems_xy_quarter, 0,
                                           0, 0, 0 ],
                              slice = roof[self.n_elems_xy_quarter + 1 ,
                                                self.n_elems_xy_quarter, -1,
                                                0, 0, -1 ], link_coeffs = [1.0],
                              value = 0.)]

        slice_4 = [BCSlice(var = 'u'  , dims = [ 2],
                            link_slice = roof[self.n_elems_xy_quarter ,
                                         self.n_elems_xy_quarter + 1, 0,
                                         0, 0, 0 ],
                            slice = roof[self.n_elems_xy_quarter  ,
                                                self.n_elems_xy_quarter + 1, -1,
                                                0, 0, -1 ], link_coeffs = [1.0],
                            value = 0.)]

        slice_5 = [BCSlice(var = 'u'  , dims = [ 2],
                            link_slice = roof[self.n_elems_xy_quarter ,
                                         self.n_elems_xy_quarter , 0,
                                         0, 0, 0 ],
                            slice = roof[self.n_elems_xy_quarter  ,
                                                self.n_elems_xy_quarter, -1,
                                                0, 0, -1 ], link_coeffs = [1.0],
                            value = 0.)]


        bc_roof_top_roof_low_link_list = slice_1 + slice_2 + slice_3 + slice_4 + slice_5

        return bc_roof_top_roof_low_link_list



    bc_plate_column_link_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_plate_column_link_list(self):
        '''
        links all column nodes to plate nodes
        '''
        column = self.fe_grid_column
        plate = self.fe_grid_plate

        slice_1 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = plate[:, :, 0, -1, -1, 0 ],
                              link_slice = column[ :, :, -1 , -1, -1, -1], link_coeffs = [1.0],
                              value = 0.)]

        slice_2 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = plate[:, :, 0, 0, 0, 0 ],
                              link_slice = column[ :, :, -1 , 0, 0, -1], link_coeffs = [1.0],
                              value = 0.)]

        slice_3 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = plate[:, :, 0, 0, -1, 0 ],
                              link_slice = column[ :, :, -1 , 0, -1, -1], link_coeffs = [1.0],
                              value = 0.)]

        slice_4 = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                              slice = plate[:, :, 0, -1, 0, 0 ],
                              link_slice = column[ :, :, -1 , -1, 0, -1], link_coeffs = [1.0],
                              value = 0.)]

        return slice_1 + slice_2 + slice_3 + slice_4
#        return [BCSlice( var = 'u'  , dims = [0, 1, 2],
#                         slice = plate[:,:,0,:,:, 0 ],
#                         link_slice = column[ :,:,-1 ,:,:,-1], link_coeffs = [1.0], value = 0. )]

    link_edge_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_link_edge_list(self):
        '''
        links all edge nodes to one node, for this node boundary conditions are applied,
        the complete force within the edge hinge can therefore be evaluated at one node
        '''
        roof = self.fe_grid_roof
        dof_constraint_0 = [BCSlice(var = 'u', dims = [1],
                                     slice = roof[ : , -1, -1, :, -1, -1],
                                    value = 0.0)]
        dof_constraint_1 = [BCSlice(var = 'u', dims = [0],
                                     slice = roof[ -1 , :, -1, -1, :, -1],
                                    value = 0.0)]
        link_edge_list = dof_constraint_0 + dof_constraint_1
        return link_edge_list

    bc_col_clamped_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_col_clamped_list(self):
        column = self.fe_grid_column
        constraint = [ BCSlice(var = 'u', dims = [0, 1, 2],
                               slice = column[ :, :, 0, :, :, 0 ],
                               value = 0.0) ]
        return constraint

    bc_col_hinge_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_col_hinge_list(self):
        constraint = []
        column = self.fe_grid_column
        for i in range(0, self.n_elems_col):
            dof_const = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                         slice = column[i , 0 , 0, 0, 0, 0 ],
                         link_slice = column[ -1 - i , -1, 0, -1 , -1, 0], link_coeffs = [-1.0],
                         value = 0.0)]
            constraint = constraint + dof_const
        for i in range(0, self.n_elems_col):
            dof_const = [BCSlice(var = 'u'  , dims = [0, 1, 2],
                         slice = column[0 , -1 - i , 0, 0, -1, 0 ],
                         link_slice = column[ -1, i , 0, -1, 0 , 0], link_coeffs = [-1.0],
                         value = 0.0)]
            constraint = constraint + dof_const

        return constraint

    #===============================================================================
    # loading cases for mr_one only symmetric 
    #===============================================================================

    lc_g_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_lc_g_list(self):
        #slices   
        roof = self.fe_grid_roof
        column = self.fe_grid_column
        upper_surf = roof[:, :, -1, :, :, -1]
        bottom_edge_roof = roof[:, 0, -1, :, 0, -1]
        left_edge_roof = roof[0, :, -1, 0, :, -1]

        # loads in global z- direction
        material_density_roof = -22.4e-3 # [MN/m^3]
        material_density_column = -26e-3 # [MN/m^3]
        additional_surface_load = -0.20e-3 # [MN/m^2]
        additional_t_constr = -0.02 * 22.4e-3
        edge_load = -0.35e-3 # [MN/m]
        return [ BCSlice(var = 'f', value = material_density_roof, dims = [2],
                          integ_domain = 'global',
                          slice = roof[:, :, :, :, :, :]),
                 BCSlice(var = 'f', value = material_density_column, dims = [2],
                          integ_domain = 'global',
                          slice = column[:, :, :, :, :, :]), ]
#                 BCSlice( var = 'f', value = additional_surface_load + additional_t_constr,
#                          dims = [2], integ_domain = 'global',
#                          slice = upper_surf ),
#                 BCSlice( var = 'f', value = edge_load, dims = [2],
#                          integ_domain = 'global',
#                          slice = bottom_edge_roof ),
#                 BCSlice( var = 'f', value = edge_load, dims = [2],
#                          integ_domain = 'global',
#                          slice = left_edge_roof )]

    lc_s_list = Property(List, depends_on = '+ps_levels, +input')
    @cached_property
    def _get_lc_s_list(self):
        # slices
        roof = self.fe_grid_roof
        upper_surf = roof[:, :, -1, :, :, -1]
        # loads in global z- direction
        snow_load = -0.85e-3
        return [ BCSlice(var = 'f', value = snow_load, dims = [2],
                          integ_domain = 'global',
                          slice = upper_surf)]

    lc_shrink_list = Property(List, depends_on = '+ps_levels, +input')
    def _get_lc_shrink_list (self):
        self.initial_strain_roof = True
        self.initial_strain_col = True
        self.t_up = -100
        self.t_lo = -100

    #===============================================================================
    # time loop
    #===============================================================================

    tloop = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_tloop(self):
        roof = self.fe_grid_roof
        column = self.fe_grid_column
        plate = self.fe_grid_plate


        ts = TS(sdomain = [roof, plate, column],
                 dof_resultants = True,
                 bcond_list = 

                              # boundary conditions
                              #
                              self.bc_roof_top_roof_low_link_list + 
                              self.bc_plate_column_link_list + 
                              self.bc_plate_roof_link_list + 
                              self.link_edge_list + 
                              self.bc_col_clamped_list + 

                              # loading
                              #
                              self.lc_g_list
                              ,
                 rtrace_list = self.rtrace_list
               )

        # Add the time-loop control
        tloop = TLoop(tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline)
        self.edge_corner_1_dof = roof[0, 0, 0, 0, 0, 0].dofs
        self.edge_corner_2_dof = roof[-1, 0, -1, -1, 0, -1].dofs
        self.dof = roof[ -1, 0, -1, -1, 0, -1 ].dofs[0][0][0]
        self.edge_roof_top = roof[ : , -1, -1, :, -1, -1].dofs
        self.edge_roof_right = roof[-1, :, -1, -1, :, -1].dofs

        return tloop


if __name__ == '__main__':

    sim_model = MRone(#n_elems_xy_quarter = 8,
                       #n_elems_z = 1,
                       #n_elems_col = 2,
                       #vtk_r = 1.0,
                       #n_elems_per_column_height = 6,
                       )

#    interior_elems = sim_model.fe_grid_column[ 1:-1, 1:-1, :, :, :, : ].elems
#    sim_model.fe_grid_column.inactive_elems = list( interior_elems )

    do = 'ui'

    if do == 'eval':
        print 'eval', sim_model.peval()

    if do == 'ui':

#        sim_model.peval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource = sim_model)
        app.main()

    elif do == 'ps':

        sim_ps = SimPStudy(sim_model = sim_model)
        sim_ps.configure_traits()

    elif do == 'cs':

#        sim_ps = SimPStudy( sim_model = sim_model )
        pstudy = SimArray(sim_model = sim_model)
        pstudy_view = SimArrayView(model = pstudy)
#        pstudy_view._start_study()
        pstudy_view.configure_traits()


        factor_list = pstudy.factor_list


        n_levels_list = [ factor.get_n_levels() for factor in factor_list ]
#        print 'get array0', pstudy.output_array.shape
#        print 'get array1', pstudy.output_array
        from matplotlib.pyplot import *
        output_array = pstudy.output_array



        cs = 'shell'
        choosen_x , choosen_z = 10, 1


        if cs == 'shell':
            print output_array.shape
            y_0 = output_array[:, 0, :]
            y_1 = output_array[:, 1, :]
#            y_2 = output_array[:,2,:]


        if cs == 'column':

            y_0 = output_array
#            y_1 = output_array[1,:,:]
#            y_2 = output_array[2,:,:]


        # shape of input [kombinations , number of ps_level]
        #   
        input_array = pstudy.input_table

        #choosen output
        #

#        cs = Enum('column','shell') 


        idx_x_axes = Int
        idx_z_axes = Int

        for i, name in enumerate(pstudy.factor_names):
            if name == 'n_elems_col_xy':
               idx_z_axes = i
            if name == 'n_elems_col_z':
               idx_x_axes = i
            if name == 'n_elems_xy_quarter':
               idx_x_axes = i
            if name == 'n_elems_z':
               idx_z_axes = i


        def n_dof_serendipity_grid(x, y, z):
            return ((2 * x + 1) * (2 * y + 1) * (2 * z + 1) - 6 * x * y * z \
                   + ((x - 1) * y * z + (y - 1) * x * z + (z - 1) * x * y) - x * y * z) * 3.0

        if cs == 'shell':
            n_dofs = n_dof_serendipity_grid(input_array[:, 0] * 2,
                                            input_array[:, 0] * 2,
                                            input_array[:, 1])
            n_dofs_0 = n_dofs[where(input_array[:, idx_z_axes] == 1)]
            n_dofs_1 = n_dofs[where(input_array[:, idx_z_axes] == 2)]
            n_dofs_2 = n_dofs[where(input_array[:, idx_z_axes] == 3)]

        if cs == 'column':
            n_dofs = n_dof_serendipity_grid(ones_like(input_array) * 2,
                                            ones_like(input_array) * 2,
                                            input_array)
            n_dofs_0 = n_dofs
#            n_dofs_0 = n_dofs[where(input_array[:,idx_z_axes]==2)]
#            n_dofs_1 = n_dofs[where(input_array[:,idx_z_axes]==4)]
#            n_dofs_2 = n_dofs[where(input_array[:,idx_z_axes]==4)]





        idx_1 = where(input_array[:, idx_x_axes] == choosen_x)
        idx_2 = where(input_array[:, idx_z_axes] == choosen_z)

        for i in  idx_1[0]:
            for j in idx_2[0]:
                if i == j:
                    idx_choosen = i

#        for i in  idx_1[0]:
#            idx_choosen = i

        n_dofs_choosen = n_dofs[idx_choosen]

        y_choosen = c_[y_0,
                       y_1,
#                       y_2
                       ][where(c_[n_dofs_0,
                       n_dofs_1,
#                       n_dofs_2
                       ] == n_dofs_choosen)]

        fig = figure(facecolor = 'white')
        ax1 = fig.add_subplot(1, 1, 1)

#        from matplotlib import rc
#        rc( 'text', usetex = True )
#        rc( 'font', **{'family':'serif', 'serif':'Times'} )
#        rc( 'text.latex',  )
#        ax1.set_xticks((50000,100000,300000))
#        ax1.set_xticklabels(("50000","100000", '300000'))

        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20)

        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20)



        if cs == 'shell':

            ax1.plot(n_dofs_0 , y_0, color = 'b', label = '1', linewidth = 1.5)
            ax1.plot(n_dofs_1 , y_1, color = 'g', label = '2', linewidth = 1.5)
#            ax1.plot(n_dofs_2 ,y_2, color = 'r', label = '3', linewidth=1.5)
            ax1.plot(n_dofs_choosen.reshape(-1, 1), y_choosen.reshape(-1, 1), marker = 'o', color = 'b', markerfacecolor = 'b', markersize = 8)#
            xlim (0, 15000)

        if cs == 'column':
            xlim (0, 2000)
            ax1.plot(n_dofs_0 , y_0, color = 'b', label = '2', linewidth = 1.5)
#            ax1.plot(n_dofs_1 ,y_1, color = 'g', label = '4', linewidth=1.5)
#            ax1.plot(n_dofs_2 ,y_2, color = 'r', label = '4', linewidth=1.5)
            ax1.plot(n_dofs_choosen.reshape(-1, 1), y_choosen.reshape(-1, 1), marker = 'o', color = 'b', markerfacecolor = 'b', markersize = 8)#

        ylim (-0.0031, -0.0029)

#        ax1.plot( [0.25,0.75],[0.45/2**0.5,2], 'bo')
        ax1.set_xlabel('Freiheitsgrade ', fontsize = 24)
        ax1.set_ylabel('U [m]', fontsize = 24)

        legend()
        show()

        filename = 'cs_' + cs + '.csv'
#        
        print '*** writing study data to file,', filename, ' ***'

        if cs == 'column':
            X_data = vstack((hstack((n_dofs_0.reshape(-1, 1), ones_like(n_dofs_0.reshape(-1, 1)) * 2, y_0))))
        if cs == 'shell':
            X_data = vstack((hstack((n_dofs_0.reshape(-1, 1), ones_like(n_dofs_0.reshape(-1, 1)) * 1, y_0)),
                            hstack((n_dofs_1.reshape(-1, 1), ones_like(n_dofs_1.reshape(-1, 1)) * 2, y_1))))



        file = open(filename, 'w')
#        
        writer = csv.writer(file, delimiter = ";", lineterminator = "\n")
        writer.writerow(['nodes' + cs , 'n_elems_t', 'U_z'])
        writer.writerows(X_data)

        file = file.close()





#        pyplot.plot
#        sim_array.clear_cache()
#        print 'array_content', sim_array.output_array
#        sim_array.configure_traits()

#        sim_ps.configure_traits()

    elif do == 'pickle':

        import pickle
        filename = '/tmp/sim.pickle'
        file = open(filename, 'w')
        pickle.dump(sim_model, file)
        file.close()
        file = open(filename, 'r')
        sm = pickle.load(file)
        file.close()

