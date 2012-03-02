'''
Created on Nov 19, 2010

@author: kelidas
'''
#import wxversion

#wxversion.select( '2.8' )

from enthought.mayavi.core.api import PipelineBase
from enthought.mayavi.core.ui.api import MayaviScene, SceneEditor, \
    MlabSceneModel
from enthought.mayavi.modules.axes import Axes

from enthought.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str
from enthought.traits.ui.api import View, Item, Group
from numpy import ones_like, array
from ymb_data import IYMBData, YMBSegmentData, YMBCutData, YMBSource, var_dict
from enthought.mayavi import mlab
from enthought.mayavi.core.api import Engine

class YMBView3D( HasTraits ):
    data = Instance( IYMBData )

    var_enum = Trait( 'contact fraction', var_dict )
    scalar_arr = Property( depends_on = 'var_enum' )
    @cached_property
    def _get_scalar_arr( self ):
        return getattr( self.data, self.var_enum_ )

    zero = Constant( 0 )

    slider_max = Property
    def _get_slider_max( self ):
        return self.data.n_filaments - 1

    n_fiber = Property( depends_on = 'data' )
    start_fib = Range( 'zero', 'slider_max', 0, mode = 'slider' )
    end_fib = Range( 'zero', 'slider_max', 'slider_max', mode = 'slider' )

    color_map = Str( 'blue-red' )

    scene = Instance( MlabSceneModel, () )
    plot = Instance( PipelineBase )

    # When the scene is activated, or when the parameters are changed, we
    # update the plot.
    @on_trait_change( 'scalar_arr, start_fib, end_fib, scene.activated' )
    def update_plot( self ):

        x_arrr, y_arrr, z_arrr = self.data.cut_data[0:3]
        scalar_arrr = self.scalar_arr

        x_arr = x_arrr.filled()#[self.start_fib:self.end_fib + 1, :].filled()
        y_arr = y_arrr.filled()#[self.start_fib:self.end_fib + 1, :].filled()
        z_arr = z_arrr.filled()#[self.start_fib:self.end_fib + 1, :].filled()
        scalar_arr = scalar_arrr.filled()#[self.start_fib:self.end_fib + 1, :].filled()

        mask = y_arr > -1

        x = x_arr[mask]
        y = y_arr[mask]
        z = z_arr[mask]
        scalar = scalar_arr[mask]

        connections = -ones_like( x_arrr )
        mesk = x_arrr.filled() > -1
        connections[mesk] = range( 0, len( connections[mesk] ) )
        connections = connections[self.start_fib:self.end_fib + 1, :].filled()
        connection = connections.astype( int ).copy()
        connection = connection.tolist()

        # TODO: better
        for i in range( 0, self.data.n_cols + 1 ):
            for item in connection:
                try:
                    item.remove( -1 )
                except:
                    pass

        if self.plot is None:
            print 'plot 3d -- 1'
            #self.scene.parallel_projection = False
            pts = self.scene.mlab.pipeline.scalar_scatter( array( x ), array( y ),
                                                        array( z ), array( scalar ) )
            pts.mlab_source.dataset.lines = connection
            self.plot = self.scene.mlab.pipeline.surface( 
                    self.scene.mlab.pipeline.tube( 
#                        fig.scene.mlab.pipeline.stripper( 
                            pts, figure = self.scene.mayavi_scene ,
#                        ),
                        tube_sides = 10, tube_radius = 0.015,
                    ),
                )
            self.plot.actor.mapper.interpolate_scalars_before_mapping = True
            self.plot.module_manager.scalar_lut_manager.show_scalar_bar = True
            self.plot.module_manager.scalar_lut_manager.show_legend = True
            self.plot.module_manager.scalar_lut_manager.shadow = True
            self.plot.module_manager.scalar_lut_manager.label_text_property.italic = False


            self.plot.module_manager.scalar_lut_manager.scalar_bar.orientation = 'horizontal'
            self.plot.module_manager.scalar_lut_manager.scalar_bar_representation.position2 = array( [ 0.61775334, 0.17 ] )
            self.plot.module_manager.scalar_lut_manager.scalar_bar_representation.position = array( [ 0.18606834, 0.08273163] )
            self.plot.module_manager.scalar_lut_manager.scalar_bar.width = 0.17000000000000004

            self.plot.module_manager.scalar_lut_manager.lut_mode = self.color_map#'black-white'
            self.plot.module_manager.scalar_lut_manager.data_name = self.var_enum
            self.plot.module_manager.scalar_lut_manager.label_text_property.font_family = 'times'
            self.plot.module_manager.scalar_lut_manager.label_text_property.shadow = True
            self.plot.module_manager.scalar_lut_manager.title_text_property.color = ( 0.0, 0.0, 0.0 )
            self.plot.module_manager.scalar_lut_manager.label_text_property.color = ( 0.0, 0.0, 0.0 )
            self.plot.module_manager.scalar_lut_manager.title_text_property.font_family = 'times'
            self.plot.module_manager.scalar_lut_manager.title_text_property.shadow = True


            #fig.scene.parallel_projection = True
            self.scene.scene.background = ( 1.0, 1.0, 1.0 )
            self.scene.scene.camera.position = [16.319534155794827, 10.477447863842627, 6.1717943847883232]
            self.scene.scene.camera.focal_point = [3.8980860486356859, 2.4731178194274621, 0.14856957086692035]
            self.scene.scene.camera.view_angle = 30.0
            self.scene.scene.camera.view_up = [-0.27676100729835512, -0.26547169369097656, 0.92354107904740446]
            self.scene.scene.camera.clipping_range = [7.7372124315754673, 26.343575352248056]
            self.scene.scene.camera.compute_view_plane_normal()
            #fig.scene.reset_zoom()

            axes = Axes()
            self.scene.engine.add_filter( axes, self.plot )
            axes.label_text_property.font_family = 'times'
            axes.label_text_property.shadow = True
            axes.title_text_property.font_family = 'times'
            axes.title_text_property.shadow = True
            axes.property.color = ( 0.0, 0.0, 0.0 )
            axes.title_text_property.color = ( 0.0, 0.0, 0.0 )
            axes.label_text_property.color = ( 0.0, 0.0, 0.0 )
            axes.axes.corner_offset = .1
            axes.axes.x_label = 'x'
            axes.axes.y_label = 'y'
            axes.axes.z_label = 'z'
        else:
            print 'plot 3d -- 2'
            #self.plot.mlab_source.dataset.reset()
            #self.plot.mlab_source.set( x = x, y = y, z = z, scalars = scalar )
            #self.plot.mlab_source.dataset.points = array( [x, y, z] ).T
            self.plot.mlab_source.scalars = scalar
            self.plot.mlab_source.dataset.lines = connection
            self.plot.module_manager.scalar_lut_manager.data_name = self.var_enum

    # The layout of the dialog created
    view = View( Item( 'scene', editor = SceneEditor( scene_class = MayaviScene ),
                     height = 250, width = 300, show_label = False ),
                Group( 
                        '_', 'start_fib', 'end_fib', 'var_enum',
                     ),
                resizable = True,
                )



if __name__ == '__main__':

    yarn_type = 'MAG'
    my_model = YMBView3D( data = YMBCutData( source = YMBSource( yarn_type = yarn_type ) ) )
    my_model.configure_traits()
