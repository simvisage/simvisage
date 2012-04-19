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
# Created on Jul 15, 2010 by: rch

# This script shows how to impose initial strain for 1-2-3 D domains
#

from etsproxy.traits.api import \
    HasTraits, Float, Int, Property

from numpy import \
    array, zeros, arange, array_equal, hstack, dot, sqrt, argsort, diag, shape

from ibvpy.api import \
    TStepper as TS, TLoop, \
    TLine, BCSlice, RTraceDomainListField

from ibvpy.plugins.ibvpy_app import \
    IBVPyApp

from ibvpy.mesh.fe_grid import \
    FEGrid

#--------------------------------------------------------------------------------
# Tracers for all three examples
#--------------------------------------------------------------------------------
sig_trace = RTraceDomainListField( name = 'Stress' ,
                                   var = 'sig_app', warp = False,
 #                                  position = 'int_pnts',
                                   record_on = 'update' )
eps_trace = RTraceDomainListField( name = 'Epsilon' ,
                                   var = 'eps_app', warp = True,
                                   record_on = 'update' )
eps0_trace = RTraceDomainListField( name = 'Epsilon 0' ,
                                   var = 'eps0_app', warp = True,
                                   record_on = 'update' )
eps1t_trace = RTraceDomainListField( name = 'Epsilon 1-t' ,
                                   var = 'eps1t_app', warp = True,
                                   record_on = 'update' )
max_p_sig_trace = RTraceDomainListField( name = 'max principle stress' , idx = 0,
                                      var = 'max_principle_sig', warp = True,
                                      record_on = 'update', )
u_trace = RTraceDomainListField( name = 'Displacement' ,
                                   var = 'u', warp = True,
                                   record_on = 'update' )

class TemperatureLinFn( HasTraits ):
    '''
    Callable class defining linear profile of initial strain. 
    '''
    alpha = Float( 1e-3 )
    T_right = Float( 50 )
    T_left = Float( 50 )
    T_ref = Float( 0 )

    strain_right = Property
    def _get_strain_right( self ):
        return self.alpha * ( self.T_right - self.T_ref )

    strain_left = Property
    def _get_strain_left( self ):
        return self.alpha * ( self.T_left - self.T_ref )

    length = Float( 10.0 )

    offset = Float( 0.0 )

    n_dims = Int( 1 )

    def __call__( self, X_pnt, x_pnt = None ):

        x = x_pnt[0]
        delta_strain = self.strain_right - self.strain_left
        if x < self.offset:
            eps = 0.0
        else:
            eps = delta_strain / self.length * ( x - self.offset ) + self.strain_left
        # return the initial volumetric strain tensor with n_dims  
        return diag( [ eps for i in range( self.n_dims ) ] )

#---------------------------------------------------------------------------------------
# 1D example
#---------------------------------------------------------------------------------------

from matplotlib import pylab as p

def plot_sig( sig_trace, text, legend ):
    '''Helper plotting function for 1D postprocessing in matplotlib.
    '''
    sig_field = sig_trace.subfields[0]
    xdata = sig_field.vtk_X[:, 0]
    ydata = sig_field.field_arr[:, :, 0][:, 0]
    idata = argsort( xdata )
    p.plot( xdata[idata], ydata[idata] )
    legend.append( text )

def demo1d():

    # Geometry
    #
    length = 1.0

    # Material and FE Formulation
    #
    from ibvpy.fets.fets1D import FETS1D2L, FETS1D2L3U
    from ibvpy.mats.mats1D import MATS1DElastic

    mats_eval = MATS1DElastic( E = 100.,
                               initial_strain = TemperatureLinFn( length = length, n_dims = 1, offset = 0.5 ) )
    fets_eval = FETS1D2L3U( mats_eval = mats_eval )
    fets_eval.vtk_r *= 0.99

    # Discretization
    #
    domain = FEGrid( coord_max = ( length, 0., 0. ),
                     n_elems = ( 10, ),
                     fets_eval = fets_eval )

    bcond_list = [BCSlice( var = 'u', dims = [0], slice = domain[0, 0], value = 0 ),
                  #BCSlice( var = 'u', dims = [0], slice = domain[-1, -1], value = 0 )
                  ]

    ts = TS( sdomain = domain,
             bcond_list = bcond_list,
             rtrace_list = [ sig_trace, eps_trace, eps0_trace, eps1t_trace ] )

    # Time integration
    #
    tloop = TLoop( tstepper = ts,
                        tline = TLine( min = 0.0, step = 1, max = 1.0 ) )

    tloop.eval()

    # Postprocessing
    #
    legend = []
    plot_sig( eps_trace, 'eps', legend )
    plot_sig( eps0_trace, 'eps0', legend )
    plot_sig( eps1t_trace, 'eps1t', legend )
    p.legend( legend )
    p.show()

#---------------------------------------------------------------------------------------
# 2D example
#---------------------------------------------------------------------------------------

def demo2d():

    # Geometry
    #
    length = 1.0

    from ibvpy.fets.fets2D import FETS2D4Q, FETS2D4Q8U, FETS2D4Q12U
    from ibvpy.mats.mats2D import MATS2DElastic

    # Material and FE Formulation
    #
    lin_x_temperature = TemperatureLinFn( length = length, n_dims = 2 )
    fets_eval = FETS2D4Q12U( mats_eval = MATS2DElastic( E = 30e5, nu = 0.2,
                                                        initial_strain = lin_x_temperature ) )
    fets_eval.vtk_r *= 0.99

    # Discretization
    #
    domain = FEGrid( coord_max = ( length, length, 0. ),
                                shape = ( 10, 10 ),
                                fets_eval = fets_eval )

    bcond_list = [BCSlice( var = 'u', dims = [0, 1], slice = domain[0, 0, 0, 0], value = 0 ),
                  BCSlice( var = 'u', dims = [1], slice = domain[0, -1, 0, -1], value = 0 ),
                  ]
    rtrace_list = [ sig_trace, eps_trace, eps0_trace, eps1t_trace, u_trace ]
    ts = TS( sdomain = domain,
             bcond_list = bcond_list,
             rtrace_list = rtrace_list,
             )

    # Time integration
    #
    tloop = TLoop( tstepper = ts,
                   tline = TLine( min = 0.0, step = 1, max = 1.0 ) )
    tloop.eval()


    # Postprocessing
    #
    app = IBVPyApp( ibv_resource = tloop )
    app.main()

#---------------------------------------------------------------------------------------
# 3D example
#---------------------------------------------------------------------------------------

def demo3d():

    # Geometry
    #
    length = 1.0

    from ibvpy.fets.fets3D import FETS3D8H, FETS3D8H20U, FETS3D8H27U, FETS3D8H20U
    from ibvpy.mats.mats3D import MATS3DElastic

    # Material and FE Formulation
    #
    lin_x_temperature = TemperatureLinFn( length = length, n_dims = 3, offset = 0.5 )
    fets_eval = FETS3D8H20U( mats_eval = MATS3DElastic( E = 30e3, nu = 0.2,
                                                       initial_strain = lin_x_temperature ) )
    fets_eval.vtk_r *= 0.99

    # Discretization
    #
    domain = FEGrid( coord_max = ( length, length, length ),
                     shape = ( 6, 3, 3 ),
                     fets_eval = fets_eval )

    bcond_list = [BCSlice( var = 'u', dims = [0, 1, 2], slice = domain[0, 0, 0, 0, 0, 0], value = 0 ),
                  BCSlice( var = 'u', dims = [0, 1], slice = domain[0, 0, -1, 0, 0, -1], value = 0 ),
                  BCSlice( var = 'u', dims = [0], slice = domain[0, -1, 0, 0, -1, 0], value = 0 ),
                  ]
    rtrace_list = [ sig_trace, eps_trace, eps0_trace, eps1t_trace, max_p_sig_trace, u_trace ]
    for rtrace in rtrace_list:
        rtrace.position = 'int_pnts'
        rtrace.warp = False

    corner_dof = domain[-1, -1, -1, -1, -1, -1].dofs[0, 0, 2]


    ts = TS( sdomain = domain,
             bcond_list = bcond_list,
             rtrace_list = rtrace_list
             )

    # Time integration
    #

    tloop = TLoop( tstepper = ts,
                   tline = TLine( min = 0.0, step = 3, max = 1.0 ) )
    tloop.eval()

    # Postprocessing
    #
    app = IBVPyApp( ibv_resource = tloop )
    app.main()


if __name__ == '__main__':

    #demo1d()
    demo2d()
    #demo3d()
