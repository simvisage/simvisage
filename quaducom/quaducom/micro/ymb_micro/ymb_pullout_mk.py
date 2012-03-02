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
# Created on Dec 20, 2010 by: rch

from enthought.traits.api import Float, Property, cached_property, Int, \
    on_trait_change, Interface, implements, HasTraits, Enum, Instance, Event
from enthought.traits.ui.api import Item, View, Group, Handler, HSplit, VSplit, \
    HGroup, HSplit, VGroup, Tabbed, Label, Spring
from enthought.traits.ui.menu import OKButton, CancelButton

from mathkit.mfn import MFnLineArray
from matplotlib.figure import Figure
from numpy import array, linspace, frompyfunc, mean, sort, trapz, sqrt
from stats.pdistrib.pdistrib import IPDistrib
from util.traits.editors.mpl_figure_editor import MPLFigureEditor
from ymb_data import YMBSource, YMBCutData, YMBSlider, IYMBData
from ymb_pullout import YarnPullOut
from stats.pdistrib.pdistrib import PDistrib
from quaducom.crackbridge.yarn_symmetrical import DoublePulloutSym

class MKYarnPDistrib( HasTraits ):

    implements( IPDistrib )

    n_int = Int( 50, auto_set = False, enter_set = True, input = True )
    p_s = Float( 0.0, auto_set = False, enter_set = True, input = True )
    p_c = Float( 0.0, auto_set = False, enter_set = True, input = True )
    k_s = Float( 0.0, auto_set = False, enter_set = True, input = True )
    k_c = Float( 0.0, auto_set = False, enter_set = True, input = True )

    x_array = Property( depends_on = '+input' )
    @cached_property
    def _get_x_array( self ):
        a = min( self.p_s, self.p_c )
        b = max( self.p_s, self.p_c )
        return linspace( a, b, self.n_int )

    polynom = Property( depends_on = '+input' )
    @cached_property
    def _get_polynom( self ):
        p_c = self.p_c
        p_s = self.p_s
        k_s = self.k_s
        k_c = self.k_c
        a = min( self.p_s, self.p_c )
        b = max( self.p_s, self.p_c )
        xi = linspace( -0.2, 1.2, 1000 )
        x = ( p_c - p_s ) * ( ( k_s + k_c - 2 ) * xi ** 3
                              + ( 3 - 2 * k_s - k_c ) * xi ** 2 + k_s * xi ) + p_s
        x = sort( x )
        xi = sort( xi )
        return MFnLineArray( xdata = x, ydata = xi )

    cdf_array = Property( depends_on = '+input' )
    @cached_property
    def _get_cdf_array( self ):
        line = self.polynom
        X = self.x_array
        pyCDF = frompyfunc( line.get_value, 1, 1 )
        return array( pyCDF( X ), dtype = 'float_' )

    pdf_array = Property( depends_on = '+input' )
    @cached_property
    def _get_pdf_array( self ):
        line = self.polynom
        X = self.x_array
        pyf = frompyfunc( line.get_diff, 1, 1 )
        PDF = array( pyf( X ), dtype = 'float_' )
        return PDF

    dx = Property()
    def _get_dx( self ):
        return abs( self.p_c - self.p_s ) / ( self.n_int - 1 )

    def get_pdf_array( self, x_array ):
        line = self.polynom
        X = x_array
        pyf = frompyfunc( line.get_diff, 1, 1 )
        PDF = array( pyf( X ), dtype = 'float_' )
        return PDF

    def integ( self ):
        return trapz( self.pdf_array, x = self.x_array )

    data_mean = Property( Float, depends_on = '+input' )
    @cached_property
    def _get_data_mean( self ):
        X = self.x_array
        x = linspace( X[0], X[-1], 300 )
        PDF = self.get_pdf_array( x )
        data_mean = trapz( x * PDF, x = x )
        return data_mean

    data_stdev = Property( Float, depends_on = '+input' )
    @cached_property
    def _get_data_stdev( self ):
        X = self.x_array
        x = linspace( X[0], X[-1], 300 )
        PDF = self.get_pdf_array( x )
        data_stdev = sqrt( trapz( x ** 2 * PDF, x = x ) - self.data_mean ** 2 )
        return data_stdev

    figure = Instance( Figure )
    def _figure_default( self ):
        figure = Figure( facecolor = 'white' )
        return figure

    data_changed = Event
    @on_trait_change( '+input' )
    def refresh( self ):
        figure = self.figure
        figure.clear()
        axes = figure.gca()
        # plot PDF and CDF
        X = self.x_array
        x = linspace( X[0], X[-1], 300 )
        PDF = self.get_pdf_array( x )
        line = self.polynom
        pyCDF = frompyfunc( line.get_value, 1, 1 )
        CDF = pyCDF( x )
        axes.plot( x, PDF, lw = 1.0, color = 'blue', \
                  label = 'PDF' )
        axes2 = axes.twinx()
        # plot CDF on a separate axis (tick labels left)
        axes2.plot( x, CDF, lw = 2, color = 'red', \
                  label = 'CDF' )
        # fill the unity area given by integrating PDF along the X-axis
        axes.fill_between( x, 0, PDF, color = 'lightblue',
                           alpha = 0.8, linewidth = 2 )

        # plot mean
        mean = self.data_mean
        axes.plot( [mean, mean], [0.0, self.get_pdf_array( mean )],
                   lw = 1.5, color = 'black', linestyle = '-' )

        # plot stdev
        stdev = self.data_stdev
        axes.plot( [mean - stdev, mean - stdev],
                   [0.0, self.get_pdf_array( mean - stdev )],
                   lw = 1.5, color = 'black', linestyle = '--' )
        axes.plot( [mean + stdev, mean + stdev],
                   [0.0, self.get_pdf_array( mean + stdev )],
                   lw = 1.5, color = 'black', linestyle = '--' )

        axes.legend( loc = 'upper center' )
        axes2.legend( loc = 'upper right' )
        axes.ticklabel_format( scilimits = ( -3., 4. ) )
        axes2.ticklabel_format( scilimits = ( -3., 4. ) )

        # plot limits on X and Y axes
        axes.set_ylim( 0.0, max( PDF ) * 1.15 )
        axes2.set_ylim( 0.0, 1.15 )
        axes.set_xlim( X[0],
                      X[-1] )
        axes2.set_xlim( X[0],
                      X[-1] )

        self.data_changed = True

    traits_view = View( HSplit( VGroup( Item( 'p_c' ),
                                     Item( 'p_s' ),
                                     Item( 'k_c' ),
                                     Item( 'k_s' ),
                                     Item( 'data_mean', label = 'mean' , style = 'readonly' ),
                                     Item( 'data_stdev', label = 'stdev', style = 'readonly' )
                                     ),
                              VGroup( Item( 'figure',
                                            editor = MPLFigureEditor(),
                                            show_label = False,
                                            resizable = True ),
                                id = 'mk.figure.view'
                                       ), ),
                                id = 'mk.view',
                                buttons = [OKButton, CancelButton],
                                resizable = True,
                                width = 600, height = 400
                        )

from promod.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

class MKPullOutParamDistribs( SimDBClass ):

    phi = Instance( IPDistrib )
    theta = Instance( IPDistrib )
    ell = Instance( IPDistrib )

    traits_view = View( Item( 'phi@' ),
                        Item( 'theta@' ),
                        Item( 'ell@' ) )

# Setup the database class extension 
#
MKPullOutParamDistribs.db = SimDBClassExt( 
            klass = MKPullOutParamDistribs,
            constants = {
                '1' : MKPullOutParamDistribs( 
                        phi = MKYarnPDistrib( n_int = 50, p_c = 0.,
                                            p_s = 1.,
                                            k_c = 2.109,
                                            k_s = 1.308 ),
                        theta = MKYarnPDistrib( n_int = 3, p_c = 0.04321,
                                            p_s = 0.0,
                                            k_c = 3.336,
                                            k_s = 0.1 ),
                        ell = MKYarnPDistrib( n_int = 3, p_c = 9.093,
                                            p_s = 0.0,
                                            k_c = 2.109,
                                            k_s = 1.308 )
                                           ),
                '2' : MKPullOutParamDistribs( 
                        phi = MKYarnPDistrib( n_int = 50, p_c = 0.,
                                            p_s = 1.,
                                            k_c = 2.233,
                                            k_s = 1.612 ),
                        theta = MKYarnPDistrib( n_int = 3, p_c = 0.04420,
                                            p_s = 0.0,
                                            k_c = 4.808,
                                            k_s = 0.1 ),
                        ell = MKYarnPDistrib( n_int = 3, p_c = 10.579,
                                            p_s = 0.0,
                                            k_c = 2.109,
                                            k_s = 1.612 )
                                           ),
                '3' : MKPullOutParamDistribs( 
                        phi = MKYarnPDistrib( n_int = 50, p_c = 0.,
                                            p_s = 1.,
                                            k_c = 2.668,
                                            k_s = 1.492 ),
                        theta = MKYarnPDistrib( n_int = 3, p_c = 0.00023,
                                            p_s = 0.0,
                                            k_c = 2.708,
                                            k_s = 0.1 ),
                        ell = MKYarnPDistrib( n_int = 3, p_c = 8.664,
                                            p_s = 0.0,
                                            k_c = 2.668,
                                            k_s = 1.492 )
                                           ),
                '4' : MKPullOutParamDistribs( 
                        phi = MKYarnPDistrib( n_int = 50, p_c = 0.,
                                            p_s = 1.,
                                            k_c = 1.929,
                                            k_s = 1.756 ),
                        theta = MKYarnPDistrib( n_int = 3, p_c = 0.00714,
                                            p_s = 0.0,
                                            k_c = 3.319,
                                            k_s = 0.1 ),
                        ell = MKYarnPDistrib( n_int = 3, p_c = 9.803,
                                            p_s = 0.0,
                                            k_c = 1.929,
                                            k_s = 1.756 )
                                           ),
                '5' : MKPullOutParamDistribs( 
                        phi = MKYarnPDistrib( n_int = 50, p_c = 0.,
                                            p_s = 1.,
                                            k_c = 1.834,
                                            k_s = 0.369 ),
                        theta = MKYarnPDistrib( n_int = 3, p_c = 0.00036,
                                            p_s = 0.0,
                                            k_c = 2.820,
                                            k_s = 0.1 ),
                        ell = MKYarnPDistrib( n_int = 3, p_c = 5.546,
                                            p_s = 0.0,
                                            k_c = 1.834,
                                            k_s = 0.369 )
                                           ),
                                           }
            )


class MKPullOut( HasTraits ):

    rf = Instance( DoublePulloutSym )
    def _rf_default( self ):
        return DoublePulloutSym( tau_fr = 3.14, l = 0.0, d = 25.5e-3, E_mod = 70.0e3,
                                 theta = 0.0, xi = 0.0179, phi = 1., L = 30.0 )

    #--------------------------------------------------------------------------------
    # select the concrete mixture from the concrete database:
    #--------------------------------------------------------------------------------

    param_distribs_key = Enum( MKPullOutParamDistribs.db.keys(),
                                 simdb = True, input = True,
                                 auto_set = False, enter_set = True )

    param_distribs_ref = Property( Instance( SimDBClass ),
                                     depends_on = 'param_distribs_key' )
    @cached_property
    def _get_param_distribs_ref( self ):
        return MKPullOutParamDistribs.db[ self.param_distribs_key ]

    po = Property( Instance( YarnPullOut ), depends_on = 'param_distribs_key' )
    @cached_property
    def _get_po( self ):
        pd = self.param_distribs_ref
        return YarnPullOut( rf = self.rf,
                            pdf_phi = pd.phi, pdf_phi_on = True,
                            pdf_theta = pd.theta, pdf_theta_on = True,
                            pdf_l = pd.ell, pdf_ell_on = True,
                            n_f = 1743,
                            w_max = 1.2 )

    #--------------------------------------------------------------------------------
    # view
    #--------------------------------------------------------------------------------

    traits_view = View( VSplit( 
                        VGroup( 
                            Spring(),
                            Item( 'param_distribs_key', label = 'parameters' ),
                            Spring(),
                            Item( 'param_distribs_ref@', show_label = False ),
                            Spring(),
                            id = 'mkpo.split.pd',
                            label = 'parameter distributions',
                            dock = 'tab',
                            ),
                        HSplit( 
                              Item( 'po@', show_label = False ),
                              label = 'Pull Out',
                              id = 'mkpo.split.pu',
                              scrollable = True,
                              dock = 'tab',
                              ),
                        dock = 'tab',
                        id = 'mkpo.split',
                        orientation = 'vertical',
                        ),
                        dock = 'tab',
                        id = 'mkpo',
                        scrollable = True,
                        resizable = True,
                        height = 0.4,
                        width = 0.5,
                        buttons = ['OK', 'Cancel'],
                        )


if __name__ == '__main__':

#    ell = PDistrib( n_segments = 100, distr_choice = 'uniform' )
#    ell.distr_type.set( loc = 0, scale = 0.25 )

#    phi = PDistrib( n_segments = 100, distr_choice = 'uniform' )
#    phi.distr_type.set( loc = 0, scale = 1.0 )

    #po.run = True
    po = MKPullOut()
    po.configure_traits()

