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
# Created on Mar 2, 2010 by: rch

from enthought.traits.api import \
    Float, Instance, Property, Int, Event, Callable, Button, on_trait_change, \
    cached_property, Array, HasTraits

from enthought.traits.ui.api import \
    View, Item, ToolBar, Action, \
    HSplit, VGroup, VSplit, OKButton, Group

from enthought.pyface.api import \
    ImageResource

from ibvpy.mats.mats1D.mats1D_damage.mats1D_damage import \
    MATS1DDamage

from ibvpy.mats.mats1D.mats1D_plastic.mats1D_plastic import \
    MATS1DEval, MATS1DPlastic

from ibvpy.mats.mats1D5.mats1D5_bond import \
    MATS1D5Bond

from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import \
    MATS1DElastic

from ibvpy.fets.fets1D5.fets1D52l4uLRH import \
    FETS1D52L4ULRH

from ibvpy.fets.fets1D5.fets1D52l6uLRH import \
    FETS1D52L6ULRH

from ibvpy.fets.fets1D5.fets1D52l8uLRH import \
    FETS1D52L8ULRH

from ibvpy.core.tstepper import \
    TStepper

from ibvpy.fets.fets1D.fets1D2l import \
    FETS1D2L

from ibvpy.fets.fets1D.fets1D2l3u import \
    FETS1D2L3U

from ibvpy.dots.dots_eval_avg import \
    DOTSEvalAvg

from numpy import \
    array, argsort, dot, unique, copy, ones_like, linspace, frompyfunc, \
    sign, fabs as nfabs, where, max

from math import \
    exp, fabs, pi, sqrt

from ibvpy.api import \
    RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval, IBVModel, \
    FETSEval, BCSlice, BCDofGroup

from ibvpy.mesh.fe_grid import \
    FEGrid

from ibvpy.mesh.fe_refinement_grid import \
    FERefinementGrid

from ibvpy.mesh.fe_domain import \
    FEDomain

from matplotlib.figure import \
    Figure

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from ibvpy.mats.mats1D.mats1D_damage.mats1D_damage_view import \
    MATS1DDamageView

from ibvpy.cntl.displ_avg.rt_nonlocal_averaging import \
    RTNonlocalAvg, QuarticAF, RTUAvg

from ibvpy.cntl.displ_avg.bar1d_damage import \
    MATS1DDamageWithFlaw

from matplotlib import \
    pylab as p


class FlawCenteredGeoTransform( HasTraits ):
    '''Slcale the length of the element to be maximum at center
    '''

    R = Float( 0.3 )

    C = Float( 0.5 )

    h = Callable

    def _h_default( self ):
        def h( x ):
            if fabs( x - self.C ) >= self.R:
                return 0.0
            else:
                R = self.R
                C = self.C
                return - ( R ** 2 - x ** 2 + 2 * x * C - C ** 2 ) ** 3 * ( x - C ) / ( ( R ** 2 - C ** 2 ) ** 2 * ( -7 * C ** 2 + R ** 2 ) )
                return - ( R ** 2 - x ** 2 + 2 * x * C - C ** 2 ) ** 3 * ( x - C ) / ( ( R ** 2 - C ** 2 ) ** 2 * ( -7 * C ** 2 + R ** 2 ) )
                return - ( R ** 2 - x ** 2 + 2 * x * C - C ** 2 ) ** 3 * ( x - C ) / ( ( R ** 2 - C ** 2 ) ** 2 * ( -7 * C ** 2 + R ** 2 ) )
                #return ( 1 - ( x - self.C ) ** 2 / self.R ** 2 ) ** 3 * x
        return h

    def __call__( self, points ):
        '''Return the global coordinates of the supplied local points.
        '''

        # number of local grid points for each coordinate direction
        # values must range between 0 and 1
        # xi, yi, zi = points[:,0], points[:,1], points[:,2]

        X = points[:, 0]

        dist = X - self.C

        R = X[ 1 ] - X[ 0 ]

        c = 0.1

        alpha = c + ( 1 - c ) / R * dist

        Xn = X - sign( dist ) * alpha * dist

        X_idx = where( nfabs( dist ) > self.R )[0]
        Xn[ X_idx ] = X[ X_idx ]

        p = copy( points )
        p[:, 0 ] = Xn
        return points


fcw = FlawCenteredGeoTransform()
fcw_v = frompyfunc( fcw.h, 1, 1 )

points = linspace( 0, 1, 11 )
print 'pppppppppp', fcw( points[:, None] )



class SimCrackLoc( IBVModel ):
    '''Model assembling the components for studying the restrained crack localization.
    '''

    geo_transform = Instance( FlawCenteredGeoTransform )
    def _geo_transform_default( self ):
        return FlawCenteredGeoTransform()

    shape = Int( 10, desc = 'Number of finite elements',
                   ps_levsls = ( 10, 40, 4 ), input = True )


    length = Float( 1000, desc = 'Length of the simulated region', unit = 'mm', input = True )

    flaw_position = Float( 500, input = True, unit = 'mm' )

    flaw_radius = Float( 100, input = True, unit = 'mm' )

    reduction_factor = Float( 0.9, input = True )

    elastic_fraction = Float( 0.9, input = True )

    avg_radius = Float( 400, input = True, unit = 'mm' )

    # tensile strength of concrete
    f_m_t = Float( 3.0, input = True, unit = 'MPa' )

    epsilon_0 = Property( unit = '-' )
    def _get_epsilon_0( self ):
        return self.f_m_t / self.E_m

    epsilon_f = Float( 10, input = True, unit = '-' )

    h_m = Float( 10, input = True, unit = 'mm' )

    b_m = Float( 8, input = True, unit = 'mm' )

    A_m = Property( unit = 'm^2' )
    def _get_A_m( self ):
        return self.b_m * self.h_m

    E_m = Float( 30.0e5, input = True, unit = 'MPa' )

    E_f = Float( 70.0e6, input = True, unit = 'MPa' )

    A_f = Float( 1.0, input = True, unit = 'mm^2' )

    s_crit = Float( 0.009, input = True, unit = 'mm' )

    P_f = Property( depends_on = '+input' )
    @cached_property
    def _get_P_f( self ):
        return sqrt( 4 * self.A_f * pi )

    K_b = Property( depends_on = '+input' )
    @cached_property
    def _get_K_b( self ):
        return self.T_max / self.s_crit

    tau_max = Float( 0.0, input = True, unit = 'MPa' )

    T_max = Property( depends_on = '+input', unit = 'N/mm' )
    @cached_property
    def _get_T_max( self ):
        return self.tau_max * self.P_f

    rho = Property( depends_on = '+input' )
    @cached_property
    def _get_rho( self ):
        return self.A_f / ( self.A_f + self.A_m )

    #-------------------------------------------------------
    # Material model for the matrix
    #-------------------------------------------------------
    mats_m = Property( Instance( MATS1DDamageWithFlaw ), depends_on = '+input' )
    @cached_property
    def _get_mats_m( self ):
        mats_m = MATS1DDamageWithFlaw( E = self.E_m * self.A_m,
                                       flaw_position = self.flaw_position,
                                       flaw_radius = self.flaw_radius,
                                       reduction_factor = self.reduction_factor,
                                       epsilon_0 = self.epsilon_0,
                                       epsilon_f = self.epsilon_f )
        return mats_m

    mats_f = Instance( MATS1DElastic )
    def _mats_f_default( self ):
        mats_f = MATS1DElastic( E = self.E_f * self.A_f )
        return mats_f

    mats_b = Instance( MATS1DEval )
    def _mats_b_default( self ):
        mats_b = MATS1DElastic( E = self.K_b )
        mats_b = MATS1DPlastic( E = self.K_b,
                                sigma_y = self.T_max,
                                K_bar = 0.,
                                H_bar = 0. ) # plastic function of slip
        return mats_b

    mats_fb = Property( Instance( MATS1D5Bond ), depends_on = '+input' )
    @cached_property
    def _get_mats_fb( self ):

        # Material model construction
        return MATS1D5Bond( mats_phase1 = MATS1DElastic( E = 0 ),
                            mats_phase2 = self.mats_f,
                            mats_ifslip = self.mats_b,
                            mats_ifopen = MATS1DElastic( E = 0 )   # elastic function of open - inactive
                            )

    #-------------------------------------------------------
    # Finite element type
    #-------------------------------------------------------
    fets_m = Property( depends_on = '+input' )
    @cached_property
    def _get_fets_m( self ):
        fets_eval = FETS1D2L( mats_eval = self.mats_m )
        #fets_eval = FETS1D2L3U( mats_eval = self.mats_m )
        return fets_eval

    fets_fb = Property( depends_on = '+input' )
    @cached_property
    def _get_fets_fb( self ):
        return FETS1D52L4ULRH( mats_eval = self.mats_fb )
        #return FETS1D52L6ULRH( mats_eval = self.mats_fb )
        #return FETS1D52L8ULRH( mats_eval = self.mats_fb )

    #--------------------------------------------------------------------------------------
    # Mesh integrator
    #--------------------------------------------------------------------------------------
    fe_domain_structure = Property( depends_on = '+input' )
    @cached_property
    def _get_fe_domain_structure( self ):
        '''Root of the domain hierarchy
        '''
        elem_length = self.length / float( self.shape )

        fe_domain = FEDomain()

        fe_m_level = FERefinementGrid( name = 'matrix domain',
                                       domain = fe_domain, fets_eval = self.fets_m )

        fe_grid_m = FEGrid( name = 'matrix grid',
                            coord_max = ( self.length, ),
                            shape = ( self.shape, ),
                            level = fe_m_level,
                            fets_eval = self.fets_m,
                            geo_transform = self.geo_transform )

        fe_fb_level = FERefinementGrid( name = 'fiber bond domain',
                                             domain = fe_domain, fets_eval = self.fets_fb )

        fe_grid_fb = FEGrid( coord_min = ( 0., length / 5. ),
                             coord_max = ( length, 0. ),
                             shape = ( self.shape, 1 ),
                             level = fe_fb_level,
                             fets_eval = self.fets_fb,
                             geo_transform = self.geo_transform )

        return fe_domain, fe_grid_m, fe_grid_fb, fe_m_level, fe_fb_level

    fe_domain = Property
    def _get_fe_domain( self ):
        return self.fe_domain_structure[0]

    fe_grid_m = Property
    def _get_fe_grid_m( self ):
        return self.fe_domain_structure[1]

    fe_grid_fb = Property
    def _get_fe_grid_fb( self ):
        return self.fe_domain_structure[2]

    fe_m_level = Property
    def _get_fe_m_level( self ):
        return self.fe_domain_structure[3]

    fe_fb_level = Property
    def _get_fe_fb_level( self ):
        return self.fe_domain_structure[4]
    #---------------------------------------------------------------------------
    # Load scaling adapted to the elastic and inelastic regime
    #---------------------------------------------------------------------------
    final_displ = Property( depends_on = '+input' )
    @cached_property
    def _get_final_displ( self ):
        damage_onset_displ = self.mats_m.epsilon_0 * self.length
        return damage_onset_displ / self.elastic_fraction

    step_size = Property( depends_on = '+input' )
    @cached_property
    def _get_step_size( self ):
        n_steps = self.n_steps
        return 1.0 / float( n_steps )

    time_function = Property( depends_on = '+input' )
    @cached_property
    def _get_time_function( self ):
        '''Get the time function so that the elastic regime 
        is skipped in a single step.
        '''
        step_size = self.step_size

        elastic_value = self.elastic_fraction * 0.98 * self.reduction_factor
        inelastic_value = 1.0 - elastic_value

        def ls( t ):
            if t <= step_size:
                return ( elastic_value / step_size ) * t
            else:
                return elastic_value + ( t - step_size ) * ( inelastic_value ) / ( 1 - step_size )

        return ls

    def plot_time_function( self, p ):
        '''Plot the time function.
        '''
        n_steps = self.n_steps
        mats = self.mats
        step_size = self.step_size

        ls_t = linspace( 0, step_size * n_steps, n_steps + 1 )
        ls_fn = frompyfunc( self.time_function, 1, 1 )
        ls_v = ls_fn( ls_t )

        p.subplot( 321 )
        p.plot( ls_t, ls_v, 'ro-' )

        final_epsilon = self.final_displ / self.length

        kappa = linspace( mats.epsilon_0, final_epsilon, 10 )
        omega_fn = frompyfunc( lambda kappa: mats._get_omega( None , kappa ), 1, 1 )
        omega = omega_fn( kappa )
        kappa_scaled = ( step_size + ( 1 - step_size ) * ( kappa - mats.epsilon_0 ) /
                   ( final_epsilon - mats.epsilon_0 ) )
        xdata = hstack( [array( [0.0], dtype = float ),
                         kappa_scaled ] )
        ydata = hstack( [array( [0.0], dtype = float ),
                         omega] )
        p.plot( xdata, ydata, 'g' )
        p.xlabel( 'regular time [-]' )
        p.ylabel( 'scaled time [-]' )

    run = Button

    @on_trait_change( 'run' )
    def peval( self ):
        '''Evaluation procedure.
        '''
        #mv = MATS1DDamageView( model = mats_eval )
        #mv.configure_traits()

        right_dof_m = self.fe_grid_m[-1, -1].dofs[0, 0, 0]

        right_dof_fb = self.fe_grid_fb[-1, -1, -1, -1].dofs[0, 0, 0]
        # Response tracers
        A = self.A_m + self.A_f
        self.sig_eps_m = RTraceGraph( name = 'F_u_m' ,
                                   var_y = 'F_int', idx_y = right_dof_m,
                                   var_x = 'U_k', idx_x = right_dof_m,
                                   transform_y = 'y / %g' % A )

        # Response tracers
        self.sig_eps_f = RTraceGraph( name = 'F_u_f' ,
                                   var_y = 'F_int', idx_y = right_dof_fb,
                                   var_x = 'U_k', idx_x = right_dof_fb,
                                   transform_y = 'y / %g' % A )

        self.eps_m_field = RTraceDomainListField( name = 'eps_m' ,
                                                  position = 'int_pnts',
                                                  var = 'eps_app',
                                                  warp = False )

        self.eps_f_field = RTraceDomainListField( name = 'eps_f' ,
                                                  position = 'int_pnts',
                                                  var = 'mats_phase2_eps_app',
                                                  warp = False )
        # Response tracers
        self.sig_m_field = RTraceDomainListField( name = 'sig_m' ,
                                                  position = 'int_pnts',
                                                  var = 'sig_app' )

        self.sig_f_field = RTraceDomainListField( name = 'sig_f' ,
                                                  position = 'int_pnts',
                                                  var = 'mats_phase2_sig_app' )

        self.omega_m_field = RTraceDomainListField( name = 'omega_m' ,
                                                    position = 'int_pnts',
                                                    var = 'omega',
                                                    warp = False )


        self.shear_flow_field = RTraceDomainListField( name = 'shear flow' ,
                                                    position = 'int_pnts',
                                                    var = 'shear_flow',
                                                    warp = False )

        self.slip_field = RTraceDomainListField( name = 'slip' ,
                                                position = 'int_pnts',
                                                var = 'slip',
                                                warp = False )

        avg_processor = None
        if self.avg_radius > 0.0:
            n_dofs = self.fe_domain.n_dofs
            avg_processor = RTUAvg( sd = self.fe_m_level, n_dofs = n_dofs,
                                    avg_fn = QuarticAF( radius = self.avg_radius ) )

        ts = TStepper( 
                      u_processor = avg_processor,
                      dof_resultants = True,
                      sdomain = self.fe_domain,
                      bcond_list = [# define the left clamping 
                                BCSlice( var = 'u', value = 0., dims = [0],
                                         slice = self.fe_grid_fb[ 0, 0, 0, :] ),
#                                BCSlice( var = 'u', value = 0., dims = [0], slice = self.fe_grid_m[ 0, 0 ] ),
                                # loading at the right edge
#                                 BCSlice( var = 'f', value = 1, dims = [0], slice = domain[-1, -1, -1, 0],
#                                         time_function = ls ),
                                 BCSlice( var = 'u', value = self.final_displ, dims = [0],
                                          slice = self.fe_grid_fb[-1, -1, -1, :],
                                          time_function = self.time_function ),
#                                 BCSlice( var = 'u', value = self.final_displ, dims = [0], slice = self.fe_grid_m[-1, -1],
#                                         time_function = self.time_function ),
                                # fix horizontal displacement in the top layer
#                                 BCSlice( var = 'u', value = 0., dims = [0], slice = domain[:, -1, :, -1] ),
                                # fix the vertical displacement all over the domain
                                 BCSlice( var = 'u', value = 0., dims = [1], slice = self.fe_grid_fb[ :, :, :, :] ),
#                            # Connect bond and matrix domains
                                 BCDofGroup( var = 'u', value = 0., dims = [0],
                                             get_link_dof_method = self.fe_grid_fb.get_bottom_dofs,
                                             get_dof_method = self.fe_grid_m.get_all_dofs,
                                             link_coeffs = [1.] )
                                ],
                 rtrace_list = [ self.sig_eps_m,
                                 self.sig_eps_f,
                                  self.eps_m_field,
                                  self.eps_f_field,
                                  self.sig_m_field,
                                  self.sig_f_field,
                                  self.omega_m_field,
                                  self.shear_flow_field,
                                  self.slip_field,
                                  ]
                )

        # Add the time-loop control
        tloop = TLoop( tstepper = ts, KMAX = 300, tolerance = 1e-5,
                       debug = False,
                       verbose_iteration = True,
                       verbose_time = False,
                       tline = TLine( min = 0.0, step = self.step_size, max = 1.0 ) )

        tloop.on_accept_time_step = self.plot

        U = tloop.eval()

        self.sig_eps_f.refresh()
        max_sig_m = max( self.sig_eps_m.trace.ydata )
        return array( [ U[right_dof_m], max_sig_m ], dtype = 'float_' )

    #--------------------------------------------------------------------------------------
    # Tracers
    #--------------------------------------------------------------------------------------

    def plot_sig_eps( self, p ):
        p.set_xlabel( 'control displacement [mm]' )
        p.set_ylabel( 'stress [MPa]' )

        self.sig_eps_m.refresh()
        self.sig_eps_m.trace.plot( p, 'o-' )

        self.sig_eps_f.refresh()
        self.sig_eps_f.trace.plot( p, 'o-' )

        p.plot( self.sig_eps_m.trace.xdata,
                self.sig_eps_m.trace.ydata + self.sig_eps_f.trace.ydata, 'o-' )

    def plot_eps( self, p ):
        eps_m = self.eps_m_field.subfields[0]
        xdata = eps_m.vtk_X[:, 0]
        ydata = eps_m.field_arr[:, 0, 0]
        idata = argsort( xdata )
        p.plot( xdata[idata], ydata[idata], 'o-' )

        eps_f = self.eps_f_field.subfields[1]
        xdata = eps_f.vtk_X[:, 0]
        ydata = eps_f.field_arr[:, 0, 0]
        idata = argsort( xdata )
        p.plot( xdata[idata], ydata[idata], 'o-' )

        p.set_ylim( ymin = 0 )
        p.set_xlabel( 'bar axis [mm]' )
        p.set_ylabel( 'strain [-]' )

    def plot_omega( self, p ):
        omega_m = self.omega_m_field.subfields[0]
        xdata = omega_m.vtk_X[:, 0]
        ydata = omega_m.field_arr[:]
        idata = argsort( xdata )
        p.fill( xdata[idata], ydata[idata], facecolor = 'gray', alpha = 0.2 )

        print 'max omega', max( ydata[idata] )

        p.set_ylim( ymin = 0, ymax = 1.0 )
        p.set_xlabel( 'bar axis [mm]' )
        p.set_ylabel( 'omega [-]' )

    def plot_sig( self, p ):
        sig_m = self.sig_m_field.subfields[0]
        xdata = sig_m.vtk_X[:, 0]
        ydata = sig_m.field_arr[:, 0, 0]
        idata = argsort( xdata )
        ymax = max( ydata )
        p.plot( xdata[idata], ydata[idata], 'o-' )

        sig_f = self.sig_f_field.subfields[1]
        xdata = sig_f.vtk_X[:, 0]
        ydata = sig_f.field_arr[:, 0, 0]
        idata = argsort( xdata )
        p.plot( xdata[idata], ydata[idata], 'o-' )

        xdata = sig_f.vtk_X[:, 0]
        ydata = sig_f.field_arr[:, 0, 0] + sig_m.field_arr[:, 0, 0]
        p.plot( xdata[idata], ydata[idata], 'ro-' )

        p.set_ylim( ymin = 0 ) # , ymax = 1.2 * ymax )
        p.set_xlabel( 'bar axis [mm]' )
        p.set_ylabel( 'stress [MPa]' )

    def plot_shear_flow( self, p ):
        shear_flow = self.shear_flow_field.subfields[1]
        xdata = shear_flow.vtk_X[:, 0]
        ydata = shear_flow.field_arr[:, 0] / self.P_f
        idata = argsort( xdata )
        ymax = max( ydata )
        p.plot( xdata[idata], ydata[idata], 'o-' )

        p.set_xlabel( 'bar axis [mm]' )
        p.set_ylabel( 'shear flow [N/m]' )

    def plot_slip( self, p ):
        slip = self.slip_field.subfields[1]
        xdata = slip.vtk_X[:, 0]
        ydata = slip.field_arr[:, 0]
        idata = argsort( xdata )
        ymax = max( ydata )
        p.plot( xdata[idata], ydata[idata], 'ro-' )

        p.set_xlabel( 'bar axis [mm]' )
        p.set_ylabel( 'slip [N/m]' )

    def plot_tracers( self, p = p ):

        p.subplot( 221 )
        self.plot_sig_eps( p )

        p.subplot( 223 )
        self.plot_eps( p )

        p.subplot( 224 )
        self.plot_sig( p )

    #---------------------------------------------------------------
    # PLOT OBJECT
    #-------------------------------------------------------------------
    figure_ld = Instance( Figure )
    def _figure_ld_default( self ):
        figure = Figure( facecolor = 'white' )
        figure.add_axes( [0.12, 0.13, 0.85, 0.74] )
        return figure

    #---------------------------------------------------------------
    # PLOT OBJECT
    #-------------------------------------------------------------------
    figure_eps = Instance( Figure )
    def _figure_eps_default( self ):
        figure = Figure( facecolor = 'white' )
        figure.add_axes( [0.12, 0.13, 0.85, 0.74] )
        return figure

    #---------------------------------------------------------------
    # PLOT OBJECT
    #-------------------------------------------------------------------
    figure_shear_flow = Instance( Figure )
    def _figure_shear_flow_default( self ):
        figure = Figure( facecolor = 'white' )
        figure.add_axes( [0.12, 0.13, 0.85, 0.74] )
        return figure

    #---------------------------------------------------------------
    # PLOT OBJECT
    #-------------------------------------------------------------------
    figure_sig = Instance( Figure )
    def _figure_sig_default( self ):
        figure = Figure( facecolor = 'white' )
        figure.add_axes( [0.12, 0.13, 0.85, 0.74] )
        return figure

    #---------------------------------------------------------------
    # PLOT OBJECT
    #-------------------------------------------------------------------
    figure_shear_flow = Instance( Figure )
    def _figure_shear_flow_default( self ):
        figure = Figure( facecolor = 'white' )
        figure.add_axes( [0.12, 0.13, 0.85, 0.74] )
        return figure

    def plot( self ):

        self.figure_ld.clear()
        ax = self.figure_ld.gca()
        self.plot_sig_eps( ax )

        self.figure_eps.clear()
        ax = self.figure_eps.gca()
        self.plot_eps( ax )
        ax2 = ax.twinx()
        self.plot_omega( ax2 )

        self.figure_sig.clear()
        ax = self.figure_sig.gca()
        self.plot_sig( ax )
        ax2 = ax.twinx()
        self.plot_omega( ax2 )

        self.figure_shear_flow.clear()
        ax = self.figure_shear_flow.gca()
        self.plot_shear_flow( ax )
        ax2 = ax.twinx()
        self.plot_slip( ax2 )

        self.data_changed = True

    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut( name = 'right end displacement', unit = 'm' ),
                 SimOut( name = 'peak load', uni = 'MPa' ) ]

    data_changed = Event

    toolbar = ToolBar( 
                      Action( name = "Run",
                             tooltip = 'Start computation',
                             image = ImageResource( 'kt-start' ),
                             action = "start_study" ),
                      Action( name = "Pause",
                             tooltip = 'Pause computation',
                             image = ImageResource( 'kt-pause' ),
                             action = "pause_study" ),
                      Action( name = "Stop",
                             tooltip = 'Stop computation',
                             image = ImageResource( 'kt-stop' ),
                             action = "stop_study" ),
                      image_size = ( 32, 32 ),
                      show_tool_names = False,
                      show_divider = True,
                      name = 'view_toolbar' ),

    traits_view = View( 
                    HSplit( 
                        VSplit( 
                                Item( 'run', show_label = False ),
                            VGroup( 
                                 Item( 'shape' ),
                                 Item( 'n_steps' ),
                                 Item( 'length' ),
                                 label = 'parameters',
                                 id = 'crackloc.viewmodel.factor.geometry',
                                 dock = 'tab',
                                 scrollable = True,
                                 ),
                            VGroup( 
                                Item( 'E_m' ),
                                Item( 'f_m_t' ),
                                Item( 'avg_radius' ),
                                Item( 'h_m' ),
                                Item( 'b_m' ),
                                Item( 'A_m', style = 'readonly' ),
                                Item( 'mats_m', show_label = False ),
                                 label = 'Matrix',
                                 dock = 'tab',
                                 id = 'crackloc.viewmodel.factor.matrix',
                                 scrollable = True
                                 ),
                            VGroup( 
                                Item( 'E_f' ),
                                Item( 'A_f' ),
                                Item( 'mats_f', show_label = False ),
                                 label = 'Fiber',
                                 dock = 'tab',
                                 id = 'crackloc.viewmodel.factor.fiber',
                                 scrollable = True
                                 ),
                            VGroup( 
                                Group( 
                                 Item( 'tau_max' ),
                                 Item( 's_crit' ),
                                 Item( 'P_f', style = 'readonly', show_label = True ),
                                 Item( 'K_b', style = 'readonly', show_label = True ),
                                 Item( 'T_max', style = 'readonly', show_label = True ),
                                 ),
                                 Item( 'mats_b', show_label = False ),
                                 label = 'Bond',
                                 dock = 'tab',
                                 id = 'crackloc.viewmodel.factor.bond',
                                 scrollable = True
                                ),
                            VGroup( 
                                 Item( 'rho', style = 'readonly', show_label = True ),
                                 label = 'Composite',
                                 dock = 'tab',
                                 id = 'crackloc.viewmodel.factor.jcomposite',
                                 scrollable = True
                                ),
                            id = 'crackloc.viewmodel.left',
                            label = 'studied factors',
                            layout = 'tabbed',
                            dock = 'tab',
                            ),
                               VSplit( 
                                    VGroup( 
                                        Item( 'figure_ld', editor = MPLFigureEditor(),
                                             resizable = True, show_label = False ),
                                             label = 'stress-strain',
                                            id = 'crackloc.viewmode.figure_ld_window',
                                            dock = 'tab',
                                    ),
                                    VGroup( 
                                        Item( 'figure_eps', editor = MPLFigureEditor(),
                                             resizable = True, show_label = False ),
                                             label = 'strains profile',
                                            id = 'crackloc.viewmode.figure_eps_window',
                                            dock = 'tab',
                                        ),
                                    VGroup( 
                                        Item( 'figure_sig', editor = MPLFigureEditor(),
                                             resizable = True, show_label = False ),
                                             label = 'stress profile',
                                            id = 'crackloc.viewmode.figure_sig_window',
                                            dock = 'tab',
                                        ),
                                    VGroup( 
                                        Item( 'figure_shear_flow', editor = MPLFigureEditor(),
                                             resizable = True, show_label = False ),
                                             label = 'bond shear and slip profiles',
                                            id = 'crackloc.viewmode.figure_shear_flow_window',
                                            dock = 'tab',
                                        ),
                                   id = 'crackloc.viewmodel.right',
                                 ),
                            id = 'crackloc.viewmodel.splitter',
                        ),
                        title = 'SimVisage Component: Crack localization',
                        id = 'crackloc.viewmodel',
                        dock = 'tab',
                        resizable = True,
                        height = 0.8, width = 0.8,
                        buttons = [OKButton] )

if __name__ == '__main__':
    # define a bar

    shape = 151
    n_steps = 30
    length = 1000
    flaw_radius = length / float( shape ) / 2.0

    sim_model = SimCrackLoc( shape = shape,
                             length = length,
                             elastic_fraction = 0.70,
                             n_steps = n_steps,
                             E_m = 28700, # [MPa]
                             h_m = 5, # [mm]
                             b_m = 5, # [mm]
                             f_t = 3.0, # [MPa]
                             epsilon_f = 0.0002,
                             E_f = 70000, # [MPa]
                             A_f = 1.0, # [mm^2]
                             s_crit = 0.009, # [ mm ]
                             tau_max = 8.96, # [ MPa*mm = N/mm ]
                             flaw_position = 500,
                             reduction_factor = 0.99,
                             flaw_radius = flaw_radius,
                             avg_radius = 10.0
                             )

    # if the composite strain is 3 MPa and the characteristic length about 0.1m
    # what bond level / reinforcement stiffness leads to the increase
    # the matrix breaking stress? 
    # stress transfer length ad 3 MPa must be 3 = tau * 0.1, i.e. tau = 3 * 10 = 30 N/m

    do = 'ui'

    if do == 'eval':
        u = sim_model.peval()
        print u
        sim_model.plot_tracers()
        p.show()
    if do == 'ui':
        sim_model.configure_traits()
