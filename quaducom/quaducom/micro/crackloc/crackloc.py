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

from etsproxy.traits.api import \
    Float, Instance, Property, Int, Event, Callable, Button, on_trait_change, \
    cached_property, Array

from etsproxy.traits.ui.api import \
    View, Item, ToolBar, Action, \
    HSplit, VGroup, VSplit, OKButton

from etsproxy.pyface.api import \
    ImageResource

from ibvpy.mats.mats1D.mats1D_damage.mats1D_damage import \
    MATS1DDamage

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
    array, argsort, dot

from math import \
    exp

from ibvpy.api import \
    RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval, IBVModel, \
    FETSEval, BCSlice

from ibvpy.mesh.fe_grid import \
    FEGrid

from matplotlib.figure import \
    Figure

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from ibvpy.mats.mats1D.mats1D_damage.mats1D_damage_view import \
    MATS1DDamageView

class MATS1DCrackLoc( MATS1DDamage ):
    '''Specialized damage model.
    
    The damage model is driven by a single damaage varialbe omega_0
    at the point x = 0. The material points are damage according
    to the nonlocal distribution function alpha implemnetd
    in the get_alpha procedure.
    
    The implementation reuses the standard MATS1DDamage but replaces
    the state variables at each material point by the single shared
    variable omega_0.
    '''

    R = Float( 0.2, desc = 'Interaction radius',
               enter_set = True, auto_set = False, modified = True )

    # state variable - damage in the discrete driving crack at x = 0
    # the value of the damage parameter in the vicinity of the 
    # driving crack is driven by the spatial damage distribution 
    # function alpha.
    # 
    _omega_0_trial = Float( 0.0 )
    _omega_0 = Float( 0.0 )
    _d_omega_0_trial = Float( 0.0 )

    # maximum strain at x = 0 achieved during the loading history.
    # used for the evaluation of the damage rule. 
    # 
    _kappa_0_trial = Float( 0.0 )
    _kappa_0 = Float( 0.0 )

    # default parameters of the damage rule 
    # onset of damage 
    epsilon_0 = 0.001
    # slope of the damage
    epsilon_f = 0.005

    #stiffness = 'algorithmic'

    traits_view = View( Item( 'R' ),
                        resizable = True
                         )

    def reset_state( self ):
        self._omega_0_trial = 0.0
        self._omega_0 = 0.0
        self._d_omega_0_trial = 0.0
        self._kappa_0 = 0.0
        self._kappa_0_trial = 0.0

    def get_alpha( self, X ):
        '''Quartic polynomial for spatial distribution of damage centered at x = 0'''
        if X <= self.R:
            return ( 1 - X ** 2 / self.R ** 2 ) ** 2
        else:
            return 0.0

    #-----------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------

    def get_corr_pred( self, sctx, eps_app_eng, d_eps, tn, tn1, eps_avg = None ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        if eps_avg != None:
            pass
        else:
            eps_avg = eps_app_eng
        #print "IN mtrl corr_pred"
        E = self.E
        D_el = array( [E] )

        if sctx.update_state_on:
            eps_n = eps_avg - d_eps

            e_max, omega = self._get_state_variables( sctx, eps_n, update_on = True )

            sctx.mats_state_array[:] = array( [e_max, omega] )

        e_max, omega = self._get_state_variables( sctx, eps_app_eng )

        if self.stiffness == "algorithmic":
            D_e_dam = array( [self._get_alg_stiffness( sctx, eps_app_eng,
                                                        e_max,
                                                        omega )] )
        else:
            D_e_dam = array( [( 1 - omega ) * D_el] )

        sigma = dot( array( [( 1 - omega ) * D_el] ), eps_app_eng )

        # print the stress you just computed and the value of the apparent E
        X_pnt = sctx.fets_eval.get_X_pnt( sctx )
        alpha = self.get_alpha( X_pnt[0] )
        print 'X', X_pnt, 'alpha', alpha, 'sigma', sigma
        return  sigma, D_e_dam

    def _get_state_variables( self, sctx, eps_app_eng, update_on = False ):
        '''Overloaded function of MATS1DDamage.
        
        Instead of using the state variables in the material point
        use the shared omega_0_trial. 
        '''
        e_max, omega = sctx.mats_state_array

        #print 'e_max', e_max, 'omega', omega

        X_pnt = sctx.fets_eval.get_X_pnt( sctx )
        alpha = self.get_alpha( X_pnt[0] )
        omega_0 = 0.5
        e_max = abs( eps_app_eng )
        return e_max, alpha * omega_0

        f_trial = abs( eps_app_eng ) - e_max
        if f_trial > 0 :
            e_max = abs( eps_app_eng )

            X_pnt = sctx.fets_eval.get_X_pnt( sctx )
            alpha = self.get_alpha( X_pnt[0] )

            if update_on:
                omega_0 = self._omega_0
            else:
                omega_0 = self._omega_0_trial
            omega = alpha * omega_0

        return e_max, omega

    def set_trial_state( self, eps_app_eng ):
        '''Set the state variables for the current trial step.
        
        Prior to calling the standard corrector predictor
        set the kappa_0_trial and omega_0_trial respecting
        the damage parameters and supplied strain state 
        in eps_app_eng.
        '''
        epsilon_0 = self.epsilon_0
        epsilon_f = self.epsilon_f

        # get the maximum strain achieved so far by comparing the
        # current and stored strains.
        eps_0 = eps_app_eng[0]
        kappa = max( self._kappa_0, eps_0 )
        self._kappa_0_trial = kappa

        if kappa >= epsilon_0 :
            self._omega_0_trial = 1. - epsilon_0 / kappa * exp( -1 * ( kappa - epsilon_0 ) / ( epsilon_f - epsilon_0 ) )
        else:
            self._omega_0_trial = 0.

    def update_state( self ):
        '''Update state.
        
        Make the trials to the accepted state variables. 
        '''
        self._omega_0 = self._omega_0_trial
        self._kappa_0 = self._kappa_0_trial


    def _get_alg_stiffness( self, sctx, eps_app_eng, e_max, omega ):
        '''Not working at the moment - probably caused by the term dodk
        
        that should be the derivative along the alpha * omega_0 - but actually
        omega_0 the dependence on epsilon_x is not explicitly available
        - thus, the determination is not possible 
        '''
        E_el = self.E
        epsilon_0 = self.epsilon_0
        epsilon_f = self.epsilon_f

        e_max_old, omega = sctx.mats_state_array
        e_max_old = max( e_max_old, epsilon_0 )
        if eps_app_eng < e_max_old: # epsilon_0:
            return E_el
        else:

            X_pnt = sctx.fets_eval.get_X_pnt( sctx )
            alpha = self.get_alpha( X_pnt[0] )

            dodk = epsilon_0 / ( e_max ** 2 ) * exp( -( e_max - epsilon_0 ) / ( epsilon_f - epsilon_0 ) ) + \
                   epsilon_0 / e_max / ( epsilon_f - epsilon_0 ) * exp( -( e_max - epsilon_0 ) / ( epsilon_f - epsilon_0 ) )
            E_inel = ( 1 - omega ) * E_el - E_el * eps_app_eng * alpha * dodk
            return E_inel


class TSCrackLoc( TStepper ):
    '''Modified time stepper with adjusted predictor and corrector.
    
    The management of the shared variable omega_0 is done before
    calling the standard corrector predictor. The strain state
    at the x = 0 is evaluated and provided to the MATS1DCrackLoc
    material model in order to evaluate the next trial damage state.
    
    During the corrpred call, the material model then uses these
    variables during the stress and stiffness integration.
    '''
    mats_eval = Property( Instance( MATS1DCrackLoc ) )
    def _get_mats_eval( self ):
        return self.subdomains[0].fets_eval.mats_eval.mats_phase1

    on_update = Callable

    fe_grid = Property( Instance( FEGrid ) )
    @cached_property
    def _get_fe_grid( self ):
        # get the single subgrid in the first subdomain.
        return self.subdomains[0].fe_subgrids[0]

    concrete_dofs = Property( Array( Float ) )
    @cached_property
    def _get_concrete_dofs( self ):
        # get the dofs associated with the bottom layer - concrete
        return self.fe_grid[:, 0, :, 0].dofs[:, :, 0]

    # specialize the 
    #
    def eval( self, step_flag, U_k, d_U, t_n, t_n1 ):

        print 'concrete dofs'
        print self.concrete_dofs

        if step_flag == 'corrector':

            # when in the corrector - switch the trial to accepted
            #
            self.mats_eval.update_state()
            self.on_update()

            #print 'E', self.subdomains[0].dots.state_array

        #----------------------------------
        # calculate the epxilon_0 at x = 0 
        #----------------------------------
        fe_grid = self.subdomains[0]
        # pick the first element
        elem_0 = fe_grid.elements[0]
        fets_eval = fe_grid.fets_eval
        # get the displacement vector in the first element
        u = U_k[ elem_0.get_dof_map() ]
        # and the coordinates
        # set the spatial context
        self.sctx.X = elem_0.get_X_mtx()
        self.sctx.x = elem_0.get_x_mtx()
        self.sctx.loc = array( [0], dtype = 'float_' )
        # get epsilon
        eps = fets_eval.get_eps_eng( self.sctx, u )
        # set the trial state for the obtained strain
        self.mats_eval.set_trial_state( eps )

        print 'BEGIN TS'
        # run the standard calculation
        return super( TSCrackLoc, self ).eval( step_flag, U_k, d_U, t_n, t_n1 )

class SimCrackLoc( IBVModel ):
    '''Model assembling the components for studying the restrained crack localization.
    '''
    shape = Int( 1, desc = 'Number of finite elements',
                   ps_levsls = ( 10, 40, 4 ) )
    length = Float( 1, desc = 'Length of the simulated region' )

    #-------------------------------------------------------
    # Material model for the matrix
    #-------------------------------------------------------
    mats_m = Instance( MATS1DCrackLoc )
    def _mats_m_default( self ):
        E_m = 1
        mats_m = MATS1DCrackLoc( E = E_m,
                                 R = 0.5,
                                 epsilon_0 = 0.001,
                                 epsilon_f = 0.01 )
        return mats_m

    mats_f = Instance( MATS1DElastic )
    def _mats_f_default( self ):
        E_f = 0
        mats_f = MATS1DElastic( E = E_f )
        return mats_f

    mats_b = Instance( MATS1DElastic )
    def _mats_b_default( self ):
        E_b = 0
        mats_b = MATS1DElastic( E = E_b )
        return mats_b

    mats = Instance( MATS1D5Bond )
    def _mats_default( self ):

        # Material model construction
        mats = MATS1D5Bond( mats_phase1 = self.mats_m,
                            mats_phase2 = self.mats_f,
                            mats_ifslip = self.mats_b,
                            mats_ifopen = MATS1DElastic( E = 0 )   # elastic function of open - inactive
                            )
        return mats
    #-------------------------------------------------------
    # Finite element type
    #-------------------------------------------------------
    fets = Instance( FETSEval )
    def _fets_default( self ):
        #return FETS1D52L8ULRH( mats_eval = self.mats )        
        return FETS1D52L4ULRH( mats_eval = self.mats )
        #return FETS1D2L3U( mats_eval = self.mats ) 

    run = Button

    @on_trait_change( 'run' )
    def peval( self ):
        '''Evaluation procedure.
        '''
        #mv = MATS1DDamageView( model = mats_eval )
        #mv.configure_traits()

        self.mats_m.reset_state()
        # Discretization
        #
        length = self.length
        domain = FEGrid( coord_min = ( 0., length / 5. ),
                          coord_max = ( length, 0. ),
                          shape = ( self.shape, 1 ),
                          fets_eval = self.fets )

        right_dofs = domain[-1, -1, -1, :].dofs[0, :, 0]
        print 'concrete_dofs', id( domain ), domain[:, 0, :, 0].dofs
        # Response tracers
        self.stress_strain = RTraceGraph( name = 'Fi,right over u_right (iteration)' ,
                                   var_y = 'F_int', idx_y = right_dofs[0],
                                   var_x = 'U_k', idx_x = right_dofs[0] )
        self.eps_m_field = RTraceDomainListField( name = 'eps_m' , position = 'int_pnts',
                                           var = 'eps1', idx = 0,
                                           warp = True )

        self.eps_f_field = RTraceDomainListField( name = 'eps_f' , position = 'int_pnts',
                                           var = 'eps2', idx = 0,
                                           warp = True )

        # Response tracers
        self.sig_m_field = RTraceDomainListField( name = 'sig_m' , position = 'int_pnts',
                                              var = 'mats_phase1_sig_app', idx = 0 )
        self.sig_f_field = RTraceDomainListField( name = 'sig_f' , position = 'int_pnts',
                                              var = 'mats_phase2_sig_app', idx = 0 )
        self.omega_m_field = RTraceDomainListField( name = 'omega_m' , position = 'int_pnts',
                                           var = 'mats_phase1_omega', idx = 0,
                                           warp = True )
        # 

        damage_onset_displ = self.mats_m.epsilon_0 * self.length
        go_behind = 1.5
        finish_displ = go_behind * damage_onset_displ
        n_steps = 20
        step_size = ( finish_displ - damage_onset_displ ) / n_steps
        tmax = 1 + n_steps

        def ls( t ):
            if t <= 1:
                return t
            else:
                return 1.0 + ( t - 1.0 ) / n_steps * ( go_behind - 1 )

        ts = TSCrackLoc( 
                 dof_resultants = True,
                 on_update = self.plot,
                 sdomain = domain,
                 bcond_list = [# define the left clamping 
                                BCSlice( var = 'u', value = 0., dims = [0], slice = domain[ 0, 0, 0, :] ),
                                # loading at the right edge
                                 BCSlice( var = 'f', value = 1, dims = [0], slice = domain[-1, -1, -1, 0],
                                         time_function = ls ),
#                                 BCSlice(var='u', value = finish_displ, dims = [0], slice = domain[-1,-1,-1, 0],
#                                         time_function = ls ),
                                # fix horizontal displacement in the top layer
                                 BCSlice( var = 'u', value = 0., dims = [0], slice = domain[:, -1, :, -1] ),
                                # fix the vertical displacement all over the domain
                                 BCSlice( var = 'u', value = 0., dims = [1], slice = domain[ :, :, :, :] )
                                ],
                 rtrace_list = [ self.stress_strain,
                                  self.eps_m_field,
                                  self.eps_f_field,
                                  self.sig_m_field,
                                  self.sig_f_field,
                                  self.omega_m_field ]
                )

        # Add the time-loop control
        tloop = TLoop( tstepper = ts, KMAX = 200, debug = True, tolerance = 1e-5,
                       tline = TLine( min = 0.0, step = 1.0, max = 1.0 ) )

        print ts.rte_dict.keys()
        U = tloop.eval()

        self.plot()

        return array( [ U[right_dofs[-1]] ], dtype = 'float_' )

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
    figure_omega = Instance( Figure )
    def _figure_omega_default( self ):
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

    def plot( self ):
        self.stress_strain.refresh()
        t = self.stress_strain.trace
        ax = self.figure_ld.axes[0]
        ax.clear()
        ax.plot( t.xdata, t.ydata, linewidth = 2, color = 'blue' )

        ax = self.figure_eps.axes[0]
        ax.clear()

        ef = self.eps_m_field.subfields[0]
        xdata = ef.vtk_X[:, 0]
        ydata = ef.field_arr[:, 0]
        idata = argsort( xdata )
        ax.plot( xdata[idata], ydata[idata], linewidth = 2, color = 'grey' )

        ef = self.eps_f_field.subfields[0]
        xdata = ef.vtk_X[:, 0]
        ydata = ef.field_arr[:, 0]
        idata = argsort( xdata )
        ax.plot( xdata[idata], ydata[idata], linewidth = 2, color = 'red' )
        ax.legend( ['eps_m', 'eps_f'] )

        ax = self.figure_sig.axes[0]
        ax.clear()

        ef = self.sig_m_field.subfields[0]
        xdata = ef.vtk_X[:, 0]
        ydata = ef.field_arr[:, 0]
        idata = argsort( xdata )
        ax.plot( xdata[idata], ydata[idata], linewidth = 2, color = 'grey' )

        ef = self.sig_f_field.subfields[0]
        xdata = ef.vtk_X[:, 0]
        ydata = ef.field_arr[:, 0]
        idata = argsort( xdata )
        ax.plot( xdata[idata], ydata[idata], linewidth = 2, color = 'red' )
        ax.legend( ['sig_m', 'sig_f'] )


        ax = self.figure_omega.axes[0]
        ax.clear()

        ef = self.omega_m_field.subfields[0]
        xdata = ef.vtk_X[:, 0]
        ydata = ef.field_arr[:, 0]
        idata = argsort( xdata )
        ax.plot( xdata[idata], ydata[idata], linewidth = 2, color = 'brown' )
        ax.legend( ['omega_m'] )


        self.data_changed = True

    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut( name = 'right end displacement', unit = 'm' ) ]

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
                                 Item( 'n_sjteps' ),
                                 Item( 'length' ),
                                 label = 'parameters',
                                 id = 'crackloc.viewmodel.factor.geometry',
                                 dock = 'tab',
                                 scrollable = True,
                                 ),
                            VGroup( 
                                 Item( 'mats_m@', show_label = False ),
                                 label = 'Matrix',
                                 dock = 'tab',
                                 id = 'crackloc.viewmodel.factor.matrix',
                                 scrollable = True
                                 ),
                            VGroup( 
                                 Item( 'mats_f@', show_label = False ),
                                 label = 'Fiber',
                                 dock = 'tab',
                                 id = 'crackloc.viewmodel.factor.fiber',
                                 scrollable = True
                                 ),
                            VGroup( 
                                 Item( 'mats_b@', show_label = False ),
                                 label = 'Bond',
                                 dock = 'tab',
                                 id = 'crackloc.viewmodel.factor.bond',
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
                                        Item( 'figure_omega', editor = MPLFigureEditor(),
                                             resizable = True, show_label = False ),
                                             label = 'omega profile',
                                            id = 'crackloc.viewmode.figure_omega_window',
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
    sim_model = SimCrackLoc()

    do = 'ui'

    if do == 'eval':
        print 'eval', sim_model.peval()
    if do == 'ui':
        sim_model.configure_traits()

