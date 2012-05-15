

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

from time import time

def temperature_strain( X_pnt, x_pnt ):

    alpha = 1.3e-5;
    t_up = -0.5
    t_lo = +0.5
    delta_t = t_lo + ( t_up - t_lo ) * x_pnt[2]
    epsilon_0 = alpha * delta_t

    return diag( [ epsilon_0 for i in range( 3 ) ] )

class SFBMushRoofModel( IBVModel ):
    '''SFB - Demonstrator model specification.
    '''
    implements( ISimModel )

    # dimensions of one quarter of the shell structure [m]
    #
    length_xy = Float( 4. )
    length_z = Float( 1.062 )
    t_shell = Float( 0.06 )

    # choose model and discretization:
    #
    mushroof_part = 'one'
    n_elems_xy_quarter = Int( 5 )
    n_elems_z = Int( 1 )

    # grid parameters:
    #
    shift_elems_column = Bool( False, input = True )
    const_edge_elem = Bool ( False, input = True )
    width_column = Float( 0.45, input = True ) # [m]
    n_elems_col = Int( 1, input = True )
    delta_h_scalefactor = Float( 1.00, ps_levels = ( 1.0, 1.6 , 4 ) )

    # Material properties: Youngs-modulus, poision ratio 
    #
    E_roof = Float( 28700 ) # [MPa]
    E_column = Float( 37000 ) # [MPa] Annahme C60
    nu = Float( 0.2 ) # [-]

    n_elems_xy = Property( Int , depends_on = '+ps_levels' )
    def _get_n_elems_xy( self ):
        return int ( self.n_elems_xy_quarter * 2 )


    # variable type of the finite element
    fets = Instance( FETSEval, ps_levels = ['fe_linear',
                                            'fe2d5_quad_serendipity',
                                            'fe_quad_serendipity',
                                            'fe_quad_lagrange' ] )
    def _fets_default( self ):
#        return self.fe_quad_serendipity
        return self.fe_linear

    mats_roof = Instance( MATS3DElastic )
    def _mats_roof_default( self ):
        return MATS3DElastic( E = self.E_roof, nu = self.nu , initial_strain = temperature_strain )

    mats_column = Instance( MATS3DElastic )
    def _mats_column_default( self ):
        return MATS3DElastic( E = self.E_column, nu = self.nu )

    fe_linear = Instance( FETSEval, transient = True )
    def _fe_linear_default( self ):
        return FETS3D8H( mats_eval = self.mats_roof )

    fe_quad_serendipity = Instance( FETSEval, transient = True )
    def _fe_quad_serendipity_default( self ):
        return FETS3D8H20U( mats_eval = self.mats_roof )

    fe2d5_quad_serendipity = Instance( FETSEval, transient = True )
    def _fe2d5_quad_serendipity_default( self ):
        return FETS2D58H20U( mats_eval = self.mats_roof )

    fe_quad_lagrange = Instance( FETSEval, transient = True )
    def _fe_quad_lagrange_default( self ):
        return FETS3D8H27U( mats_eval = self.mats_roof )

    fe_column = Instance( FETSEval, transient = True )
    def _fe_column_default( self ):
        fets = FETS3D8H20U( mats_eval = self.mats_column )
        fets.vtk_r *= 1.0
        return fets

    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut( name = 'u_x_corner1', unit = 'm' ),
                 SimOut( name = 'u_x_corner2', unit = 'm' ),
                 SimOut( name = 'u_z_corner1', unit = 'm' ),
                 SimOut( name = 'u_z_corner2', unit = 'm' ),
                 SimOut( name = 'maximum principle stress', unit = 'MPa' ), ]

    def peval( self ):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()

        u_x_corner1 = U[ self.edge_corner_1_dof ][0, 0, 0]
        u_x_corner2 = U[ self.edge_corner_1_dof ][0, 0, 0]
        u_z_corner1 = U[ self.edge_corner_1_dof ][0, 0, 2]
        u_z_corner2 = U[ self.edge_corner_2_dof ][0, 0, 2]



        max_princ_stress = max( self.max_princ_stress._get_field_data().flatten() )

        return array( [ u_x_corner1,
                       u_x_corner2,
                       u_z_corner1,
                       u_z_corner2,
                       max_princ_stress ], dtype = 'float_' )


    tline = Instance( TLine )
    def _tline_default( self ):
        return TLine( min = 0.0, step = 1.0, max = 1.0 )

    max_princ_stress = Instance( RTraceDomainListField )
    def _max_princ_stress_default( self ):
        return RTraceDomainListField( name = 'max principle stress' , idx = 0,
                                      var = 'max_principle_sig', warp = True,
                                      record_on = 'update', )

    rtrace_list = List
    def _rtrace_list_default( self ):
        return [  self.max_princ_stress, self.sig_app, self.u  ]

    sig_app = Property( Instance( RTraceDomainListField ), depends_on = '+ps_levels' )
    @cached_property
    def _get_sig_app( self ):
        return RTraceDomainListField( name = 'sig_app' ,
    #                                  position = 'int_pnts',
                                      var = 'sig_app',
                                      record_on = 'update', )

    u = Property( Instance( RTraceDomainListField ), depends_on = '+ps_levels' )
    @cached_property
    def _get_u( self ):
        return RTraceDomainListField( name = 'displacement' ,
                                      var = 'u', warp = True,
                                      record_on = 'update', )


    #----------------------------------------------------------------------------------
    # HP-Shell specification
    #----------------------------------------------------------------------------------

    hp_shell = Property( Instance( HPShell ) , depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_hp_shell( self ):
        return HPShell( n_elems_xy = self.n_elems_xy,
                        n_elems_z = self.n_elems_z,
                        shift_elems_column = self.shift_elems_column,
                        delta_h_scalefactor = self.delta_h_scalefactor,
                        const_edge_elem = self.const_edge_elem,
                        width_column = self.width_column,
                        n_elems_col = self.n_elems_col,
                        mushroof_part = self.mushroof_part )

    fe_grid_roof = Property( Instance( FEGrid ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_fe_grid_roof( self ):
        return FEGrid( coord_min = ( 0.0, 0.0, 0.0 ),
                       coord_max = ( 1.0, 1.0, 1.0 ),
                       geo_transform = self.hp_shell,
                       shape = ( self.n_elems_xy, self.n_elems_xy, self.n_elems_z ),
                       fets_eval = self.fets )


    #----------------------------------------------------------------------------------
    # Column specification
    #----------------------------------------------------------------------------------

    column_width_top = Float( 0.45, unit = 'm' )
    column_width_bottom = Float( 0.3, ps_levels = ( 0.30, 0.45, 4 ) )
    column_height = Float( 3.0, unit = 'm' )

    column = Property( Instance( HPShell ) , depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_column( self ):
        return GEOColumn( width_top = self.column_width_top,
                       width_bottom = self.column_width_bottom,
                       height = self.column_height / self.delta_h_scalefactor ,
                       X0 = [ 4., 4., -3. / self.delta_h_scalefactor ] )

    fe_grid_column = Property( Instance( FEGrid ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_fe_grid_column( self ):
        fe_grid_column = FEGrid( coord_min = ( 0.0, 0.0, 0.0 ),
                              coord_max = ( 1.0, 1.0, 1.0 ),
                              geo_transform = self.column,
                              shape = ( 3, 3, 4 ),
                              fets_eval = self.fe_column )
#        interior_elems = fe_grid_column[ 1:-1, 1:-1, :, :, :, : ].elems
#        fe_grid_column.inactive_elems = list( interior_elems )
        return fe_grid_column

    # time loop
    tloop = Property( depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_tloop( self ):
        self.fets.vtk_r *= .95

        roof = self.fe_grid_roof
        column = self.fe_grid_column

        self.edge_corner_1_dof = roof[-1, -1, -1, -1, -1, -1].dofs
        self.edge_corner_2_dof = roof[-1, 0, -1, -1, 0, -1].dofs

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

        bc_link_col01 = BCSlice( var = 'u'  , dims = [0, 1, 2],
                                    slice = roof[self.n_elems_xy_quarter - 1,
                                                   self.n_elems_xy_quarter - 1,
                                                   0,
                                                   0, -1, 0 ],
                                    link_slice = column[ 0, 0, -1, 0, 0, -1],
                                    link_coeffs = [1.0],
                                    value = 0. )
        bc_link_col02 = BCSlice( var = 'u'  , dims = [0, 1, 2],
                                    slice = roof[self.n_elems_xy_quarter - 1,
                                                   self.n_elems_xy_quarter - 1 ,
                                                   0,
                                                   - 1, 0, 0],
                                    link_slice = column[ -1, 0, -1, -1, 0, -1],
                                    link_coeffs = [1.0],
                                    value = 0. )
        bc_link_col03 = BCSlice( var = 'u'  , dims = [0, 1, 2],
                                    slice = roof[self.n_elems_xy_quarter,
                                                   self.n_elems_xy_quarter,
                                                   0,
                                                   - 1, 0, 0 ],
                                    link_slice = column[ -1, -1, -1, -1, -1, -1],
                                    link_coeffs = [1.0],
                                    value = 0. )
        bc_link_col04 = BCSlice( var = 'u'  , dims = [0, 1, 2],
                                    slice = roof[self.n_elems_xy_quarter,
                                                   self.n_elems_xy_quarter,
                                                   0,
                                                   0, -1, 0 ],
                                    link_slice = column[ 0, -1, -1, 0, -1, -1],
                                    link_coeffs = [1.0],
                                    value = 0. )

        bc_column_support = BCSlice( var = 'u', dims = [0, 1, 2],
                                     slice = column[ :, :, 0, :, :, 0 ],
                                     value = 0.0 )

        #bc_corner_load   = BCSlice( var = 'f', value = -nodal_load, dims = [2], slice = roof[-1,-1,-1,-1,-1,-1] )
        #bc_topface_load  = BCSlice( var = 'f', value = -nodal_load, dims = [2], slice = roof[:,:,-1,:,:,-1] )

        w_z = roof[-1, -1, -1, -1, -1, -1].dofs[0, 0, 2]

        self.f_w_diagram = RTraceGraph( name = 'load - corner deflection',
                                           var_x = 'U_k', idx_x = w_z,
                                           var_y = 'time', idx_y = 0,
                                           record_on = 'update' )

        rtrace_list = [ self.f_w_diagram ] + self.rtrace_list

        ts = TS( sdomain = [roof, column],
                 dof_resultants = True,
                 bcond_list = [bc_link_col01,
                               bc_link_col02,
                               bc_link_col03,
                               bc_link_col04,
                               bc_column_support,
                               ], # + force_bc,
                 rtrace_list = rtrace_list
               )

        # Add the time-loop control
        tloop = TLoop( tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline )

        return tloop

if __name__ == '__main__':

    sim_model = SFBMushRoofModel( n_elems_xy_quarter = 5,
                                  n_elems_z = 1,
                                  const_edge_elem = False,
                                  shift_elems_column = True,
                                  width_column = 0.45,
                                  n_elems_col = 1 )



    do = 'ui'

    if do == 'eval':
        print 'eval', sim_model.peval()

    if do == 'ui':

#        sim_model.peval()
        print sim_model.n_elems_col
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
