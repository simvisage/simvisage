
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

from enthought.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool, HasTraits

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
    cos, sin, arctan, where, abs, all, any, diag

from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos, pi as Pi

from rsurface_reader import \
    read_rsurface, normalize_rsurfaces

# Interpolation
from scipy.interpolate import Rbf

from geo_column import GEOColumn

from simiter.sim_pstudy import ISimModel, SimOut, SimPStudy

from hp_shell import HPShell

from mush_roof_model import MushRoofModel

from time import time

class MRone_free( MushRoofModel ):

    implements( ISimModel )
    mushroof_part = 'one'
    #===============================================================================
    # fe_grid
    #===============================================================================
    
    n_elems_xy_quarter = Int(4, ps_levels = [3,15,5])
    n_elems_z = Int(2, ps_levels = [1,4,2] )
    n_elems_per_column_height = Int( 3, ps_levels = [3,6,2]   )
    
    vtk_r = Float(0.95)
    #default column
    fe_column = Instance( FETSEval, transient = True )
    def _fe_column_default( self ):
        fets = self.fe_quad_serendipity
        fets.vtk_r *= self.vtk_r
        return fets
    #default roof
    fe_roof = Instance( FETSEval, ps_levels = ['fe_linear',
                                                'fe2d5_quad_serendipity',
                                                'fe_quad_serendipity',
                                                'fe_quad_lagrange' ], depends_on='+vtk_r, +initial_strain' )
    def _fe_roof_default( self ):
        fets = self.fe_quad_serendipity
        fets.vtk_r *= self.vtk_r
        return fets

    fe_grid_roof = Property( Instance( FEGrid ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_fe_grid_roof( self ):
        return FEGrid( coord_min = ( 0.0, 0.0, 0.0 ),
                       coord_max = ( 1.0, 1.0, 1.0 ),
                       geo_transform = self.hp_shell,
                       shape = ( self.n_elems_xy, self.n_elems_xy, self.n_elems_z ),
                       fets_eval = self.fe_roof )

    fe_grid_column = Property( Instance( FEGrid ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_fe_grid_column( self ):
        return  FEGrid( coord_min = ( 0.0, 0.0, 0.0 ),
                        coord_max = ( 1.0, 1.0, 1.0 ),
                        geo_transform = self.column,
                        shape = ( self.n_elems_col,self.n_elems_col, self.n_elems_per_column_height ),
                        fets_eval = self.fe_column )

    #===============================================================================
    # ps_study
    #===============================================================================
    def peval( self ):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()
        u_x_corner1 = U[ self.edge_corner_1_dof ][0, 0, 0]
        #u_x_corner2 = U[ self.edge_corner_2_dof ][0, 0, 0]
        u_z_corner1 = U[ self.edge_corner_1_dof ][0, 0, 2]
        #u_z_corner2 = U[ self.edge_corner_2_dof ][0, 0, 2]
        #max_princ_stress = max( self.max_princ_stress._get_field_data().flatten() )

        return array( [ u_x_corner1,
                      # u_x_corner2,
                       u_z_corner1,])
                       #u_z_corner2,
                       #max_princ_stress ], dtype = 'float_' )

    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut( name = 'u_x', unit = 'm' ),
                # SimOut( name = 'u_x_corner2', unit = 'm' ),
                 SimOut( name = 'u_z', unit = 'm' ),]
                 #SimOut( name = 'u_z_corner2', unit = 'm' ),
                 #SimOut( name = 'maximum principle stress', unit = 'MPa' ), ]

    link_column_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_link_column_list(self):

        roof = self.fe_grid_roof
        column = self.fe_grid_column
        link_column_list =[]
        for i in range (0,self.n_elems_col):
            slice = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter - self.n_elems_col + i,
                                           self.n_elems_xy_quarter +i, 0,
                                           0, 0, 0 ],
                              link_slice = column[ 0 , i, -1, 0, 0, -1], link_coeffs = [1.0],
                              value = 0. )]
            link_column_list = link_column_list  + slice
        for i in range (0,self.n_elems_col):
            slice = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter  + i +1 ,
                                           self.n_elems_xy_quarter -self.n_elems_col +i +1 , 0,
                                           0, 0, 0 ],
                              link_slice = column[ -1 , i, -1, -1, -1, -1],
                              link_coeffs = [1.0],
                              value = 0. )]
            link_column_list = link_column_list  + slice
        
        for i in range (0,self.n_elems_col):
            slice = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter  + i  ,
                                           self.n_elems_xy_quarter + self.n_elems_col -i , 0,
                                           0, 0, 0 ],
                              link_slice = column[ i, -1 , -1, 0, -1, -1],
                              link_coeffs = [1.0],
                              value = 0. )]
            link_column_list = link_column_list  + slice
        
        for i in range (0,self.n_elems_col):
            slice = [BCSlice( var = 'u'  , dims = [0, 1, 2],
                              slice = roof[self.n_elems_xy_quarter - self.n_elems_col + i +1 ,
                                           self.n_elems_xy_quarter - i -1 , 0,
                                           0, 0, 0 ],
                              link_slice = column[ i , 0 , -1, -1, 0, -1],
                              link_coeffs = [1.0],
                              value = 0. )]
            link_column_list = link_column_list  + slice
       
        return link_column_list
    
    #----------------------------------------------------
    # response tracer
    #----------------------------------------------------
        
    rtrace_list = List
    def _rtrace_list_default( self ):
        return [  self.max_princ_stress, self.sig_app, self.u]
    
              
    tloop = Property( depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_tloop( self ):
        roof = self.fe_grid_roof
        column = self.fe_grid_column

        self.edge_corner_1_dof = roof[0, 0, 0, 0, 0, 0].dofs
        self.edge_corner_2_dof = roof[-1, 0, -1, -1, 0, -1].dofs
        self.dof = roof[ -1 , 0, -1, -1, 0, -1].dofs[0][0][0]

        #----------------------------------------------------
        # loading cases (LC):
        #----------------------------------------------------

        #--- LC1: dead load
        # g = 22.4 kN/m^3 
        # orientation: global z-direction; 
        material_density_roof = -0.0224 # [MN/m^3]
        material_density_column = -0.026 # [MN/m^3]

        #--- LC2 additional dead load 
        # gA = 0,20 kN/m^2 
        # orientation: global z-direction (following the curved structure); 
        surface_load_gA = -0.20e-3 # [MN/m^2]

        #--- LC3 snow
        # s = 0,79 kN/m^2 
        # orientation: global z-direction (projection); 
        surface_load_full_s = -0.84e-3 # [MN/m^2]
        surface_load_half_s = -0.84e-3 * 0.5 # [MN/m^2]

        #--- LC4 wind (pressure) 
        # w = 0,13 kN/m^2 
        # orientation: local t-direction (surface normal); 
        surface_load_front_edge_w = 0.91e-3   # [MN/m^2]
#        surface_load_back_edge_w = 0.39e-3 # [MN/m^2]
        surface_load_x_pos_w = -0.13e-3
        surface_load_x_neg_w = 0.89e-3
        surface_load_side = 1.14e-3

#        surface_load_back_edge_w = # [MN/m^2]
#        surface_load_left_edge_w = # [MN/m^2]
#        surface_load_right_edge_w = # [MN/m^2]

        # NOTE: additional line-loads at the edge of the roof need to be considered!  

        upper_surface = roof[:, :, -1, :, :, -1]
        whole_roof = roof[:, :, :, :, :, :]

        #edges
        #
        front_edge_roof = roof[0, :, -1, 0, :, -1]
#        back_edge_roof = roof[-1, :, -1, -1, :, -1]
        side_edge_roof = roof[:, 0, -1, :, 0, -1]


        #sides
        #
        x_pos_surface = roof[self.n_elems_xy_quarter:, :, -1, :, :, -1]
        x_neg_surface = roof[:self.n_elems_xy_quarter, :, -1, :, :, -1]
        whole_column = column[:, :, :, :, :, :]


        force_bc = [
                     # LC1: dead load
#                     BCSlice( var = 'f', value = material_density_roof, dims = [2],
#                              integ_domain = 'global',
#                              slice = whole_roof ),
#                     BCSlice( var = 'f', value = material_density_column, dims = [2],
#                              integ_domain = 'global',
#                              slice = whole_column ),

##                     # LC2: additional dead load
                     BCSlice( var = 'f', value = -1.0e-3, dims = [2],
                              integ_domain = 'global',
                              slice = upper_surface ),
##
##                     # LC3: snow load         
#                     BCSlice( var = 'f', value = surface_load_full_s, dims = [2],
#                              integ_domain = 'global',
#                              slice = upper_surface ),
##                     BCSlice( var = 'f', value = surface_load_full_s, dims = [2],
##                              integ_domain = 'global',
##                              slice = x_pos_surface ),
##                      #LC4: wind load
##
##                     #front edge
###                     BCSlice( var = 'f', value = surface_load_front_edge_w, dims = [0],
###                              integ_domain = 'global',
###                              slice = front_edge_roof ),
#####                     #back edge
####                     BCSlice( var = 'f', value = surface_load_back_edge_w, dims = [0],
####                              integ_domain = 'global',
####                              slice = back_edge_roof ),
###                     #top
#                      BCSlice( var = 'f', value = surface_load_x_pos_w, dims = [2],
#                              integ_domain = 'global',
#                              slice = whole_roof ),
###                     #side
###                     BCSlice( var = 'f', value = surface_load_side, dims = [1],
###                              integ_domain = 'global',
###                              slice = side_edge_roof ),


                    ]




        bc_list = [# column
                        BCSlice( var = 'u', dims = [0, 1, 2],
                                     slice = column[ :, :, 0, :, :, 0 ],
                                     value = 0.0 ),]
                                

        #bc_corner_load   = BCSlice( var = 'f', value = -nodal_load, dims = [2], slice = roof[-1,-1,-1,-1,-1,-1] )
        #bc_topface_load  = BCSlice( var = 'f', value = -nodal_load, dims = [2], slice = roof[:,:,-1,:,:,-1] )

        w_z = roof[-1, -1, -1, -1, -1, -1].dofs[0, 0, 2]

#        self.f_w_diagram = RTraceGraph( name = 'load - corner deflection',
#                                           var_x = 'U_k', idx_x = w_z,
#                                           var_y = 'time', idx_y = 0,
#                                           record_on = 'update' )
#
#        rtrace_list = [ self.f_w_diagram ] + self.rtrace_list

        ts = TS( sdomain = [roof, column],
                 dof_resultants = True,
                 bcond_list = bc_list + force_bc +self.link_column_list,
#                 bcond_list = [bc_column_support,
#                               ] + force_bc,

                 rtrace_list = self.rtrace_list
               )

        # Add the time-loop control
        tloop = TLoop( tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline )

        return tloop

if __name__ == '__main__':

    sim_model = MRone_free( n_elems_xy_quarter = 5,
                            n_elems_z = 2,
                            shift_elems_column = True,
                            n_elems_col = 2,
                            vtk_r = 1.0,
                            n_elems_per_column_height = 5)

#    interior_elems = sim_model.fe_grid_column[ 1:-1, 1:-1, :, :, :, : ].elems
#    sim_model.fe_grid_column.inactive_elems = list( interior_elems )

    do = 'ui'

    if do == 'eval':
        print 'eval', sim_model.peval()

    if do == 'ui':


       sim_model.peval()
       from ibvpy.plugins.ibvpy_app import IBVPyApp
       app = IBVPyApp( ibv_resource = sim_model )
       app.main()

    elif do == 'ps':

        sim_ps = SimPStudy( sim_model = sim_model )
        sim_ps.configure_traits()

    elif do == 'pickle':

        import pickle
        filename = '/tmp/sim.pickle'
        file = open( filename, 'w' )
        pickle.dump( sim_model, file )
        file.close()
        file = open( filename, 'r' )
        sm = pickle.load( file )
        file.close()
