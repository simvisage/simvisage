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
# Created on Dec 13, 2010 by: kelidas

#import wxversion
#wxversion.select( '2.8' )

from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Tuple, Interface, implements, Any, Button, Str, List
from enthought.traits.trait_types import DelegatesTo
from enthought.traits.ui.api import Item, View, HGroup, RangeEditor, Group, HSplit, \
                                 HGroup, VGroup, VSplit, VGrid, Handler
from matplotlib.figure import Figure
from util.traits.editors.mpl_figure_editor import MPLFigureEditor
from enthought.traits.trait_types import DelegatesTo
import os
from stats.pdistrib.pdistrib import PDistrib, IPDistrib

from ymb_data import YMBSource, IYMBData, YMBCutData, YMBSegmentData, YMBData
from ymb_data import YMBSlider
from ymb_hist import YMBHist
from ymb_auto_correl import YMBAutoCorrelView, YMBAutoCorrel
from ymb_cross_correl import YMBCrossCorrel
from ymb_view3d import YMBView3D
from ymb_view2d import YMBView2D
from ymb_pullout import YMBPullOut
from ymb_pdistrib import YMBDistrib
from ymb_micro import YMBMicro

from promod.simdb import SimDB
from os.path import join

simdb = SimDB()




from matplotlib import rc, RcParams
rc( 'text', usetex = False )
rc( 'font', family = 'serif', style = 'normal', \
     variant = 'normal', stretch = 'normal', weight = '100', size = 12 ) # serif = 'Times New Roman'
rc( 'axes', labelsize = 14 )
#rc( 'legend', fontsize = 16 )

start_tex = '''
\documentclass[a4paper,10pt,oneside,final,titlepage,onecolumn]{report}

\usepackage{ucs}
\usepackage[utf8x]{inputenc}
\usepackage{amsmath}
\usepackage[english]{babel}
\usepackage[T1]{fontenc} 
\usepackage[pdftex]{graphicx} 

\usepackage[pdftex]{hyperref} 
\usepackage{graphicx}
 
\date{}

\\begin{document}
'''

fig_tex = '''
\\begin{figure}[!htb]
\includegraphics[width=\\textwidth]{figs/%s}
\caption{%s}
\end{figure}
'''

end_tex = '''\n\end{document} '''


def file_test( filepath ):
    return os.path.isfile( filepath )

class histPlot( HasTraits ):
    bins = Int( 20 )

class TitleHandler( Handler ):
    """ Change the title on the UI.
    """

    def object_yarn_changed( self, info ):
        """ Called whenever the "yarn" attribute changes on the handled
        object.
        """
        info.ui.title = info.object.yarn


class YMBReport( HasTraits ):

    body_tex = Str()

    data = Instance( YMBData, changed = True )

    yarn = Property( Str, depends_on = '+changed' )
    @cached_property
    def _get_yarn( self ):
        return self.data.source.yarn_type

    data_dir = Property( Str, depends_on = '+changed' )
    def _get_data_dir( self ):
        return join( simdb.exdata_dir, 'trc', 'yarn_structure', self.yarn, 'raw_data' )
    tex_dir = Property( Str, depends_on = '+changed' )
    def _get_tex_dir( self ):
        return join( simdb.exdata_dir, 'trc', 'yarn_structure', self.yarn, 'report' )
    fig_dir = Property( Str, depends_on = '+changed' )
    def _get_fig_dir( self ):
        return join( simdb.exdata_dir, 'trc', 'yarn_structure', self.yarn, 'report', 'figs' )

    # select data for report
    plot_3d_on = Bool( True )
    hist_rad_on = Bool( True )
    hist_cf_on = Bool( True )
    hist_bfl_on = Bool( True )
    hist_slack_on = Bool( True )
    corr_plot_rad_on = Bool( True )
    corr_plot_cf_on = Bool( True )
    corr_plot_bfl_on = Bool( True )
    corr_plot_slack_on = Bool( True )
    spirrid_on = Bool( False )

    hist_axes_adjust = List( [0.12, 0.17, 0.68, 0.68] )
    corr_axes_adjust = List( [0.15, 0.17, 0.75, 0.68] )
    n_bins = Int( 40 )


    #################################
    # BUILD REPORT
    #################################

    build = Button( label = 'build report' )
    def _build_fired ( self ):
        self._directory_test()

        # get the yarn type
        yt = self.data.source.yarn_type

        if self.plot_3d_on == True:
            self._save_plot3d_fired()
        if self.hist_rad_on == True:
            self._save_rad_fired()
        if self.hist_cf_on == True:
            self._save_cf_fired()
        if self.hist_bfl_on == True:
            self._save_bfl_fired()
        if self.hist_slack_on == True:

            self._save_slack_fired()
        if self.corr_plot_rad_on == True:
            self._save_corr_plot_rad_fired()
        if self.corr_plot_cf_on == True:
            self._save_corr_plot_cf_fired()
        if self.corr_plot_bfl_on == True:
            self._save_corr_plot_bfl_fired()
        if self.corr_plot_slack_on == True:
            self._save_corr_plot_slack_fired()

        print '================'
        print 'Figure(s) saved'
        print '================'


    #################################
    # BUILD REPORT
    #################################

        filename = 'ymb_report_' + yt
        texfile = join( self.tex_dir, filename + '.tex' )
        pdffile = join( self.tex_dir, filename + '.pdf' )

        bodyfile = join( self.tex_dir, texfile )
        body_out = open( bodyfile, 'w' )

        self.body_tex += 'Yarn with contact fraction limit = %s\n' % self.data.cf_limit

        if self.plot_3d_on == True:
            self.body_tex += self.plot3d_tex
        if self.hist_rad_on == True:
            self.body_tex += self.rad_tex
        if self.hist_cf_on == True:
            self.body_tex += self.cf_tex
        if self.hist_bfl_on == True:
            self.body_tex += self.bfl_tex
        if self.hist_slack_on == True:
            self.body_tex += self.slack_tex
        if self.corr_plot_rad_on == True:
            self.body_tex += self.corr_plot_rad_tex
        if self.corr_plot_cf_on == True:
            self.body_tex += self.corr_plot_cf_tex
        if self.corr_plot_bfl_on == True:
            self.body_tex += self.corr_plot_bfl_tex
        if self.corr_plot_slack_on == True:
            self.body_tex += self.corr_plot_slack_tex

        body_out.write( start_tex )
        body_out.write( '\section*{ Yarn %s }' % ( self.yarn ) )
        body_out.write( self.body_tex )
        body_out.write( end_tex )
        body_out.close()
        os.system( 'cd ' + self.tex_dir + ';pdflatex -shell-escape ' + texfile )
        print '=============================='
        print 'Report written to %s', texfile
        print '=============================='
        os.system( 'acroread ' + pdffile + ' &' )

    #################################
    # 3d PLOT
    #################################

    plot3d = Property( Instance( YMBView3D ) )
    @cached_property
    def _get_plot3d( self ):
        plot3d = YMBView3D( data = self.data,
                             color_map = 'binary' ) # black-white
        return plot3d

    save_plot3d = Button( label = 'save plot3d figure' )
    def _save_plot3d_fired( self ):
        self._directory_test()
        filename = join( self.fig_dir, 'plot3d.png' )
        self.plot3d.scene.save( filename, ( 1000, 800 ) )

    plot3d_tex = Property( Str )
    @cached_property
    def _get_plot3d_tex( self ):
        filename = 'plot3d.png'
        if file_test( join( self.fig_dir, filename ) ) == True:
            return fig_tex % ( filename, '3D yarn plot' )
        else:
            self._save_plot3d_fired()
            return fig_tex % ( filename, '3D yarn plot' )


    #################################
    # HISTOGRAM PLOT AND SAVE
    #################################

    hist_rad = Property( Instance( YMBHist ) )#Instance(Figure )
    @cached_property
    def _get_hist_rad( self ):
        histog = YMBHist()
        histog.set( edge_color = 'black', face_color = '0.75', axes_adjust = self.hist_axes_adjust )
        slider = YMBSlider( var_enum = 'radius', data = self.data )
        return histog.set( slider = slider, bins = self.n_bins, normed_on = True )

    save_rad = Button( label = 'save histogram of radius' )
    def _save_rad_fired( self ):
        self._directory_test()
        filename = join( self.fig_dir, 'radius.pdf' )
        self.hist_rad.figure.savefig( filename, format = 'pdf' )

    rad_tex = Property( Str )
    @cached_property
    def _get_rad_tex( self ):
        filename = 'radius.pdf'
        if file_test( join( self.fig_dir, filename ) ) == True:
            return fig_tex % ( filename, 'Histogram of filament radius' )
        else:
            self._save_rad_fired()
            return fig_tex % ( filename, 'Histogram of filament radius' )


    hist_cf = Property( Instance( YMBHist ) )
    @cached_property
    def _get_hist_cf( self ):
        histog = YMBHist()
        histog.set( edge_color = 'black', face_color = '0.75', axes_adjust = self.hist_axes_adjust )
        slider = YMBSlider( var_enum = 'contact fraction', data = self.data )
        return histog.set( slider = slider, bins = self.n_bins, normed_on = True )

    save_cf = Button( label = 'save histogram of contact fraction' )
    def _save_cf_fired( self ):
        self._directory_test()
        filename = join( self.fig_dir, 'contact_fraction.pdf' )
        self.hist_cf.figure.savefig( filename, format = 'pdf' )

    cf_tex = Property( Str )
    @cached_property
    def _get_cf_tex( self ):
        filename = 'contact_fraction.pdf'
        if file_test( join( self.fig_dir, filename ) ) == True:
            return fig_tex % ( filename, 'Histogram of contact fraction' )
        else:
            self._save_cf_fired()
            return fig_tex % ( filename, 'Histogram of contact fraction' )


    hist_bfl = Property( Instance( YMBHist ) )
    @cached_property
    def _get_hist_bfl( self ):
        histog = YMBHist()
        histog.set( edge_color = 'black', face_color = '0.75', axes_adjust = self.hist_axes_adjust )
        slider = YMBSlider( var_enum = 'bond free length', data = self.data )
        return histog.set( slider = slider, bins = self.n_bins, normed_on = True )

    save_bfl = Button( label = 'save histogram of bond free length' )
    def _save_bfl_fired( self ):
        self._directory_test()
        filename = 'bond_free_length.pdf'
        self.hist_bfl.figure.savefig( join( self.fig_dir, filename ), format = 'pdf' )

    bfl_tex = Property( Str )
    @cached_property
    def _get_bfl_tex( self ):
        filename = 'bond_free_length.pdf'
        if file_test( join( self.fig_dir, filename ) ) == True:
            return fig_tex % ( filename, 'Histogram of bond free length' )
        else:
            self._save_bfl_fired()
            return fig_tex % ( filename, 'Histogram of bond free length' )



    hist_slack = Property( Instance( YMBHist ) )
    @cached_property
    def _get_hist_slack( self ):
        histog = YMBHist()
        histog.set( edge_color = 'black', face_color = '0.75', axes_adjust = self.hist_axes_adjust )
        slider = YMBSlider( var_enum = 'slack', data = self.data )
        return histog.set( slider = slider, bins = self.n_bins, normed_on = True,
                           xlimit_on = True, xlimit = 0.03 )

    save_slack = Button( label = 'save histogram of slack' )
    def _save_slack_fired( self ):
        self._directory_test()
        filename = join( self.fig_dir, 'slack.pdf' )
        self.hist_slack.figure.savefig( filename, format = 'pdf' )

    slack_tex = Property( Str )
    @cached_property
    def _get_slack_tex( self ):
        filename = 'slack.pdf'
        if file_test( join( self.fig_dir, filename ) ) == True:
            return fig_tex % ( filename, 'Histogram of slack' )
        else:
            self._save_slack_fired()
            return fig_tex % ( filename, 'Histogram of slack' )


    corr_plot_rad = Property( Instance( YMBAutoCorrel ) )
    @cached_property
    def _get_corr_plot_rad( self ):
        plot = YMBAutoCorrelView()
        plot.set( color = 'black', axes_adjust = self.corr_axes_adjust )
        return plot.set( correl_data = YMBAutoCorrel( data = data, var_enum = 'radius' ) )

    save_corr_plot_rad = Button( label = 'save correlation plot of radius' )
    def _save_corr_plot_rad_fired( self ):
        self._directory_test()
        filename = join( self.fig_dir, 'corr_plot_rad.pdf' )
        self.corr_plot_rad.figure.savefig( filename, format = 'pdf' )

    corr_plot_rad_tex = Property( Str )
    @cached_property
    def _get_corr_plot_rad_tex( self ):
        filename = 'corr_plot_rad.pdf'
        if file_test( join( self.fig_dir, filename ) ) == True:
            return fig_tex % ( filename, 'Autocorrelation of radius' )
        else:
            self._save_corr_plot_rad_fired()
            return fig_tex % ( filename, 'Autocorrelation of radius' )


    corr_plot_cf = Property( Instance( YMBAutoCorrel ) )
    @cached_property
    def _get_corr_plot_cf( self ):
        plot = YMBAutoCorrelView()
        plot.set( color = 'black', axes_adjust = self.corr_axes_adjust )
        return plot.set( correl_data = YMBAutoCorrel( data = data, var_enum = 'contact fraction' ) )



    save_corr_plot_cf = Button( label = 'save correlation plot of contact fraction' )
    def _save_corr_plot_cf_fired( self ):
        self._directory_test()
        filename = join( self.fig_dir, 'corr_plot_cf.pdf' )
        self.corr_plot_cf.figure.savefig( filename, format = 'pdf' )

    corr_plot_cf_tex = Property( Str )
    @cached_property
    def _get_corr_plot_cf_tex( self ):
        filename = 'corr_plot_cf.pdf'
        if file_test( join( self.fig_dir, filename ) ) == True:
            return fig_tex % ( filename, 'Autocorrelation of contact fraction' )
        else:
            self._save_corr_plot_cf_fired()
            return fig_tex % ( filename, 'Autocorrelation of contact fraction' )



    corr_plot_bfl = Property( Instance( YMBAutoCorrel ) )
    @cached_property
    def _get_corr_plot_bfl( self ):
        plot = YMBAutoCorrelView()
        plot.set( color = 'black', axes_adjust = self.corr_axes_adjust )
        return plot.set( correl_data = YMBAutoCorrel( data = data, var_enum = 'bond free length' ) )



    save_corr_plot_bfl = Button( label = 'save corr plot of bond free length' )
    def _save_corr_plot_bfl_fired( self ):
        self._directory_test()
        filename = join( self.fig_dir, 'corr_plot_bfl.pdf' )
        self.corr_plot_bfl.figure.savefig( filename, format = 'pdf' )

    corr_plot_bfl_tex = Property( Str )
    @cached_property
    def _get_corr_plot_bfl_tex( self ):
        filename = 'corr_plot_bfl.pdf'
        if file_test( join( self.fig_dir, filename ) ) == True:
            return fig_tex % ( filename, 'Autocorrelation of bond free length' )
        else:
            self._save_corr_plot_bfl_fired()
            return fig_tex % ( filename, 'Autocorrelation of bond free length' )


    corr_plot_slack = Property( Instance( YMBAutoCorrel ) )
    @cached_property
    def _get_corr_plot_slack( self ):
        plot = YMBAutoCorrelView()
        plot.set( color = 'black', axes_adjust = self.corr_axes_adjust )
        return plot.set( correl_data = YMBAutoCorrel( data = data, var_enum = 'slack' ) )


    save_corr_plot_slack = Button( label = 'save corr plot of slack' )
    def _save_corr_plot_slack_fired( self ):
        self._directory_test()
        filename = join( self.fig_dir, 'corr_plot_slack.pdf' )
        self.corr_plot_slack.figure.savefig( filename, format = 'pdf' )

    corr_plot_slack_tex = Property( Str )
    @cached_property
    def _get_corr_plot_slack_tex( self ):
        filename = 'corr_plot_slack.pdf'
        if file_test( join( self.fig_dir, filename ) ) == True:
            return fig_tex % ( filename, 'Autocorrelation of slack' )
        else:
            self._save_corr_plot_slack_fired()
            return fig_tex % ( filename, 'Autocorrelation of slack' )


    #################################
    # CORRELATION PLOT AND TABLE
    #################################

    def corr_plot( self ):
        return 0

    #################################
    # SPIRRID PLOT
    #################################

    spirrid_plot = Property( Instance( YMBPullOut ) )
    @cached_property
    def _get_spirrid_plot( self ):
        return YMBPullOut( data = self.data )

    pdf_theta = Property( Instance( IPDistrib ) )
    def _get_pdf_theta( self ):
        theta = self.spirrid_plot.pdf_theta
        #theta.set( face_color = '0.75', edge_color = 'black', axes_adjust = [0.13, 0.18, 0.8, 0.7] )
        return theta

    pdf_l = Property( Instance( IPDistrib ) )
    def _get_pdf_l( self ):
        return self.spirrid_plot.pdf_l

    pdf_phi = Property( Instance( IPDistrib ) )
    def _get_pdf_phi( self ):
        return self.spirrid_plot.pdf_phi

    pdf_xi = Property( Instance( IPDistrib ) )
    def _get_pdf_xi( self ):
        return self.spirrid_plot.pdf_xi.figure

    def _directory_test( self ):
        if os.access( self.tex_dir, os.F_OK ) == False:
            os.mkdir( self.tex_dir )
        if os.access( self.fig_dir, os.F_OK ) == False:
            os.mkdir( self.fig_dir )


    traits_view = View( 
                HSplit( 
                Group( 
                Item( 'data@', show_label = False ),
                HGroup( 
                VGrid( 
                       Item( 'plot_3d_on', label = 'plot3d', ),
                       Item( 'save_plot3d', show_label = False,
                              visible_when = 'plot_3d_on == True' ),
                        Item( '_' ),
                       Item( 'hist_rad_on', label = 'radius_hist' ),
                       Item( 'save_rad', show_label = False,
                              visible_when = 'hist_rad_on == True' ),
                       Item( 'hist_cf_on', label = 'cf_hist' ),
                       Item( 'save_cf', show_label = False,
                              visible_when = 'hist_cf_on == True' ),
                       Item( 'hist_bfl_on', label = 'bfl_hist' ),
                       Item( 'save_bfl', show_label = False,
                              visible_when = 'hist_bfl_on == True' ),
                       Item( 'hist_slack_on', label = 'slack_hist' ),
                       Item( 'save_slack', show_label = False,
                              visible_when = 'hist_slack_on == True' ),
                        Item( '_' ),
                       Item( 'corr_plot_rad_on', label = 'corr_plot_rad' ),
                       Item( 'save_corr_plot_rad', show_label = False,
                              visible_when = 'hist_slack_on == True' ),
                       Item( 'corr_plot_cf_on', label = 'corr_plot_cf' ),
                       Item( 'save_corr_plot_cf', show_label = False,
                              visible_when = 'hist_slack_on == True' ),
                       Item( 'corr_plot_bfl_on', label = 'corr_plot_bfl' ),
                       Item( 'save_corr_plot_bfl', show_label = False,
                              visible_when = 'hist_slack_on == True' ),
                       Item( 'corr_plot_slack_on', label = 'corr_plot_slack' ),
                       Item( 'save_corr_plot_slack', show_label = False,
                              visible_when = 'hist_slack_on == True' ),
                       ),
                VGrid( 
                      Item( 'spirrid_on' ),
                      ),
                ),
                Item( 'build' ),
                label = 'report',
                id = 'report.bool',
                ),
                ),
                VGroup( 
                HGroup( 
                Item( 'hist_rad@', show_label = False ),
                Item( 'hist_cf@', show_label = False ),
                ),
                HGroup( 
                Item( 'hist_bfl@', show_label = False ),
                Item( 'hist_slack@', show_label = False ),
                ), label = 'histograms', id = 'report.hist', ),
                VGroup( 
                HGroup( 
                Item( 'corr_plot_rad@', show_label = False ),
                Item( 'corr_plot_cf@', show_label = False ),
                ),
                HGroup( 
                Item( 'corr_plot_bfl@', show_label = False ),
                Item( 'corr_plot_slack@', show_label = False ),
                ), label = 'correlation plot', id = 'report.corr_plot', ),


#                HGroup( 
#                Group( 
#                Item( 'pdf_theta@', show_label = False ),
#                Item( 'pdf_l@', show_label = False ),
#                ),
#                Group( 
#                Item( 'pdf_phi@', show_label = False ),
#                Item( 'pdf_xi@', show_label = False ),
#                ),
#                label = 'pdf',
#                id = 'report.pdf',
#                #scrollable = True,
#                ),
                Group( Item( 'plot3d@', show_label = False ),
                      label = 'plot3d',
                      id = 'report.plot3d' ),
                resizable = True,
                title = u"Yarn name",
                handler = TitleHandler(),
                id = 'report.main',
               )





if __name__ == '__main__':


    source = YMBSource( yarn_type = 'MAG' )

    data = YMBCutData( source = source, cf_limit = 0.5 )

    ymbreport = YMBReport( data = data )
    ymbreport.configure_traits()




