'''
Created on Jul 30, 2010

@author: abach
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


class MRfour( MushRoofModel ):

    #----------------------------------------------------
    # elements
    #----------------------------------------------------
    vtk_r = Float( 1.0 )
    #default column
    fe_column = Instance( FETSEval, transient = True )
    def _fe_column_default( self ):
        fets = self.fe_quad_serendipity_column
        fets.vtk_r *= self.vtk_r
        return fets
    #default roof
    fe_roof = Instance( FETSEval, ps_levels = ['fe_linear',
                                                'fe2d5_quad_serendipity',
                                                'fe_quad_serendipity_roof',
                                                'fe_quad_serendipity_column',
                                                'fe_quad_lagrange' ] )
    def _fe_roof_default( self ):
        fets = self.fe_quad_serendipity_roof
        fets.vtk_r *= self.vtk_r
        return fets


    #----------------------------------------------------
    # geometric transformation
    #----------------------------------------------------

    hp_shells = Property( List( HPShell ) , depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_hp_shells( self ):
        X_list = [ [0, 0, 0], [7, 0, 0], [0, 7, 0], [7, 7, 0] ]
        return [
                 HPShell( 
                    length_xy_quarter = self.length_xy_quarter,
                    length_z = self.length_z,
                    t_shell = self.t_shell,
                    X0 = X,
                    n_elems_xy_quarter = self.n_elems_xy_quarter,
                    n_elems_z = self.n_elems_z,
#                    shift_elems_column = self.shift_elems_column,
                    scalefactor_delta_h = self.scalefactor_delta_h,
                    const_reinf_layer_elem = self.const_reinf_layer_elem,
                    width_top_col = self.width_top_col,
#                    n_elems_col = self.n_elems_col,
                    mushroof_part = self.mushroof_part )
                 for X in X_list ]

    columns = Property( Instance( GEOColumn ) , depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_columns( self ):
        ch = self.h_col / self.scalefactor_delta_h
        X_list = [[ 3.5, 3.5, -ch ],
                  [ 10.5, 3.5, -ch ],
                  [ 3.5, 10.5, -ch ],
                  [ 10.5, 10.5, -ch ]]
        return [ GEOColumn( width_top = self.width_top_col,
                            X0 = X,
                            width_bottom = self.width_bottom_col,
                            h_col = ch )
                            for X in X_list ]

    #----------------------------------------------------
    # grid 
    #----------------------------------------------------
    fe_grid_roofs = Property( Instance( FEGrid ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_fe_grid_roofs( self ):
        return [
                FEGrid( coord_min = ( 0.0, 0.0, 0.0 ),
                       coord_max = ( 1.0, 1.0, 1.0 ),
                       geo_transform = hp_shell,
                       shape = ( self.n_elems_xy_quarter, self.n_elems_xy_quarter, self.n_elems_z ),
                       fets_eval = self.fe_roof )
                for hp_shell in self.hp_shells ]



    fe_grid_columns = Property( Instance( FEGrid ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_fe_grid_columns( self ):
        return [ FEGrid( coord_min = ( 0.0, 0.0, 0.0 ),
                         coord_max = ( 1.0, 1.0, 1.0 ),
                         geo_transform = column,
                         shape = ( 1, 1, 1 ),
                         fets_eval = self.fe_column )
                         for column in self.columns ]


    #----------------------------------------------------
    # ps_study
    #----------------------------------------------------

    def peval( self ):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()

        U_z = U[self.corner_dof_r00][0, 0, 2]

#        max_princ_stress = max( self.max_princ_stress._get_field_data().flatten() )

        return array( [ U_z ],
#                       max_princ_stress ],
                        dtype = 'float_' )
    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut( name = 'u_z_free_corner', unit = 'm' )]
#               SimOut( name = 'maximum principle stress', unit = 'MPa' ), ]


     #----------------------------------------------------
    # loading and boundaries
    #----------------------------------------------------
#    lc_symm_dict = Property( Dict )
#    def _get_lc_symm_dict ( self ):
#        return {'g_mat_dens_roof':Float( -0.0224, unit = 'MN/m^3' ),
#                'g_mat_dens_column':Float( -0.024, unit = 'MN/m^3' ),
#                'g_surf_load_gA':Float( -0.20e-3, unit = 'MN/m^2' ),
#                's_surf_load':Float( -0.84e-3, unit = 'MN/m^2' ),
#                'w_surf_load':Float( -0.13e-3, unit = 'MN/m^2' ),
#                }
    #--- LC0: verifying load for ps_study
    # g = 1.0 kN/m^2 
    # orientation: global z-direction; 
    f_z = Float( -1e-3, unit = 'MN/m^3' )

    #--- LC1: dead load
    # g = 22.4 kN/m^3 
    # orientation: global z-direction; 
    material_density_roof = Float( -0.0224, unit = 'MN/m^3' )
    material_density_column = Float( -0.0224, unit = 'MN/m^3' )

    #--- LC2 additional dead load 
    # gA = 0,20 kN/m^2 
    # orientation: global z-direction (following the curved structure); 
    surface_load_gA = Float( -0.20e-3, unit = 'MN/m^2' )

    #--- LC3 snow
    # s = 0,79 kN/m^2 
    # orientation: global z-direction (projection); 
    surface_load_s = Float( -0.84e-3, unit = 'MN/m^2' )

    #--- LC4 wind (pressure) 
    # w = 0,13 kN/m^2 
    # orientation: local t-direction (surface normal); 
    surface_load_w = Float( -0.13e-3, unit = 'MN/m^2' )

    force_bc_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_force_bc_list( self ):
        # NOTE: additional line-loads at the edge of the roof need to be considered!  

        force_bc_list = []
        for roof, column in zip( self.fe_grid_roofs, self.fe_grid_columns ):
            upper_surface = roof[:, :, -1, :, :, -1]
#            whole_roof = roof[:, :, :, :, :, :]
#            whole_column = column[:, :, :, :, :, :]

            force_bc = [
#                         # LC0: dead load
                         BCSlice( var = 'f', value = self.f_z, dims = [2],
                                  integ_domain = 'global',
                                  slice = upper_surface ),

#                         # LC1: dead load
#                         BCSlice( var = 'f', value = self.material_density_roof, dims = [2],
#                                  integ_domain = 'global',
#                                  slice = whole_roof ),
#                         BCSlice( var = 'f', value = self.material_density_column, dims = [2],
#                                  integ_domain = 'global',
#                                  slice = whole_column ),
                         # LC2: additional dead load
#                         BCSlice( var = 'f', value = self.surface_load_gA, dims = [2],
#                                  integ_domain = 'global',
#                                  slice = upper_surface ),
#                         # LC3: snow load         
#                         BCSlice( var = 'f', value = self.surface_load_s, dims = [2],
#                                  integ_domain = 'global',
#                                  slice = upper_surface )
                       ]
            force_bc_list += force_bc

        return force_bc_list

    displ_bc_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_displ_bc_list( self ):

        displ_bc_list = []
        for roof, column in zip( self.fe_grid_roofs, self.fe_grid_columns ):
            bc_list = [
                       BCSlice( var = 'u'  , dims = [0, 1, 2],
                                slice = roof[self.n_elems_xy_quarter - 1,
                                               self.n_elems_xy_quarter - 1,
                                               0,
                                               0, -1, 0 ],
                                link_slice = column[ 0, 0, -1, 0, 0, -1],
                                link_coeffs = [1.0],
                                value = 0. ),
                        BCSlice( var = 'u'  , dims = [0, 1, 2],
                                slice = roof[self.n_elems_xy_quarter - 1,
                                               self.n_elems_xy_quarter - 1 ,
                                               0,
                                               - 1, 0, 0],
                                link_slice = column[ -1, 0, -1, -1, 0, -1],
                                link_coeffs = [1.0],
                                value = 0. ),
                        BCSlice( var = 'u'  , dims = [0, 1, 2],
                                slice = roof[self.n_elems_xy_quarter,
                                               self.n_elems_xy_quarter,
                                               0,
                                               - 1, 0, 0 ],
                                link_slice = column[ -1, -1, -1, -1, -1, -1],
                                link_coeffs = [1.0],
                                value = 0. ),
                        BCSlice( var = 'u'  , dims = [0, 1, 2],
                                    slice = roof[self.n_elems_xy_quarter,
                                                   self.n_elems_xy_quarter,
                                                   0,
                                                   0, -1, 0 ],
                                    link_slice = column[ 0, -1, -1, 0, -1, -1],
                                    link_coeffs = [1.0],
                                    value = 0. ),
                        BCSlice( var = 'u', dims = [0, 1, 2],
                                 slice = column[ :, :, 0, :, :, 0 ],
#                                 slice = column[ 0, 0, 0, -1, -1, 0 ],
                                 value = 0.0 )
                        ]
            displ_bc_list += bc_list
        return displ_bc_list


    wind_bc_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_wind_bc_list( self ):

        w_left_edge = 0.91e-3 # [MN/m]
        w_right_edge = 0.39e-3 # [MN/m]
        w_bottom_edge_0 = 1.14e-3 # [Mn/m]
        w_bottom_edge_1 = 0.65e-3 # [Mn/m]
        w_top_edge_0 = -w_bottom_edge_0 # [Mn/m]
        w_top_edge_1 = -w_bottom_edge_1 # [Mn/m]
        w_face_0_left = 0.89e-3 # MN/m^2
        w_face_0_right = -0.13e-3 # MN/m^2
        w_face_1_left = 0.72e-3 # MN/m^2
        w_face_1_right = w_face_0_right

        r00, r10, r01, r11 = self.fe_grid_roofs

        left_edge_0 = r00[ 0, :, -1, 0, :, -1 ]
        left_edge_1 = r01[ 0, :, -1, 0, :, -1 ]
        right_edge_0 = r10[ -1, :, -1, 0, :, -1 ]
        right_edge_1 = r11[ -1, :, -1, 0, :, -1 ]
        bottom_edge_0 = r00[ :, 0, -1, :, 0, -1 ]
        bottom_edge_1 = r10[ :, 0, -1, :, 0, -1 ]
        top_edge_0 = r01[ :, -1, -1, :, -1, -1 ]
        top_edge_1 = r11[ :, -1, -1, :, -1, -1 ]

        n_e_q = self.n_elems_xy_quarter
        face_00_left = r00[ :n_e_q, :, -1, :, :, -1 ]
        face_00_right = r00[ n_e_q:, :, -1, :, :, -1 ]
        face_01_left = r01[ :n_e_q, :, -1, :, :, -1 ]
        face_01_right = r01[ n_e_q:, :, -1, :, :, -1 ]
        face_10_left = r10[ :n_e_q, :, -1, :, :, -1 ]
        face_10_right = r10[ n_e_q:, :, -1, :, :, -1 ]
        face_11_left = r11[ :n_e_q, :, -1, :, :, -1 ]
        face_11_right = r11[ n_e_q:, :, -1, :, :, -1 ]

        wind_list = [
                   # left edge - x direction
                   BCSlice( var = 'f'  , dims = [0],
                            slice = left_edge_0,
                            value = w_left_edge ),
                   BCSlice( var = 'f'  , dims = [0],
                            slice = left_edge_1,
                            value = w_left_edge ),
                   # right edge - x direction
                   BCSlice( var = 'f'  , dims = [0],
                            slice = right_edge_0,
                            value = w_right_edge ),
                   BCSlice( var = 'f'  , dims = [0],
                            slice = right_edge_1,
                            value = w_right_edge ),
                   # bottom edge - y direction
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_0,
                            value = w_bottom_edge_0 ),
                   BCSlice( var = 'f'  , dims = [1],
                            slice = bottom_edge_1,
                            value = w_bottom_edge_1 ),
                   # top edge - y direction
                   BCSlice( var = 'f'  , dims = [1],
                            slice = top_edge_0,
                            value = w_top_edge_0 ),
                   BCSlice( var = 'f'  , dims = [1],
                            slice = top_edge_1,
                            value = w_top_edge_1 ),
                   # upper face left - left 
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_00_left,
                            value = w_face_0_left ),
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_01_left,
                            value = w_face_0_left ),
                   # upper face left - right 
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_00_right,
                            value = w_face_0_right ),
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_01_right,
                            value = w_face_0_right ),
                   # upper face right - left 
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_10_left,
                            value = w_face_1_left ),
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_11_left,
                            value = w_face_1_left ),
                   # upper face right - right 
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_10_right,
                            value = w_face_1_right ),
                   BCSlice( var = 'f'  , dims = [2],
                            slice = face_11_right,
                            value = w_face_1_right ),
                        ]
        return wind_list

    link_bc_list = Property( List, depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_link_bc_list( self ):
        r00, r10, r01, r11 = self.fe_grid_roofs
        link_bc_list = [
                       BCSlice( var = 'u'  , dims = [0, 1, 2],
                                slice = r00[-1, :, -1, -1, :-1, -1 ],
                                link_slice = r10[ 0, :, -1, 0, :-1, -1],
                                link_coeffs = [1.0],
                                value = 0. ),
                        BCSlice( var = 'u'  , dims = [0, 1, 2],
                                slice = r01[:, 0, -1, :-1, 0, -1],
                                link_slice = r00[ :, -1, -1, :-1, -1, -1],
                                link_coeffs = [1.0],
                                value = 0. ),
                        BCSlice( var = 'u'  , dims = [0, 1, 2],
                                slice = r11[ 0, :, -1, 0, 1:, -1 ],
                                link_slice = r01[ -1, :, -1, -1, 1:, -1],
                                link_coeffs = [1.0],
                                value = 0. ),
                        BCSlice( var = 'u'  , dims = [0, 1, 2],
                                slice = r10[ :, -1, -1, 1:, -1, -1 ],
                                link_slice = r11[ :, 0, -1, 1:, 0, -1],
                                link_coeffs = [1.0],
                                value = 0. ),
                        ]
        return link_bc_list

    #----------------------------------------------------
    # response tracer
    #----------------------------------------------------

    rtrace_list = List
    def _rtrace_list_default( self ):
        return [  self.max_princ_stress, self.sig_app, self.u]



    # time loop
    tloop = Property( depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_tloop( self ):
        #self.fets.vtk_r *= 0.95

        roofs = self.fe_grid_roofs

        columns = self.fe_grid_columns

        r00, r10, r01, r11 = roofs

        self.corner_dof_r00 = r00[0, 0, 0, 0, 0, 0].dofs

        rtrace_list = self.rtrace_list

        ts = TS( sdomain = roofs + columns,
                 dof_resultants = True,
#                 bcond_list = self.wind_bc_list + self.displ_bc_list + self.link_bc_list,
                 bcond_list = self.force_bc_list + self.displ_bc_list + self.link_bc_list,
#                 bcond_list = self.displ_bc_list + self.link_bc_list,
                 rtrace_list = rtrace_list
               )

        # Add the time-loop control
        tloop = TLoop( tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline )

        return tloop





if __name__ == '__main__':

    sim_model = MRfour( n_elems_xy_quarter = 20,
                        mushroof_part = 'one',
                        n_elems_z = 20,
                        const_edge_elem = False,
                        shift_elems_column = True,
                        width_column = 0.45,
                        n_elems_col = 1 )

#    interior_elems = sim_model.fe_grid_column[ 1:-1, 1:-1, :, :, :, : ].elems
#    sim_model.fe_grid_column.inactive_elems = list( interior_elems )

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

