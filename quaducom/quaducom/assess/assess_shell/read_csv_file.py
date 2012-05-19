from numpy import *
from promod.simdb import \
    SimDB

import os
import pickle
import string
from os.path import join

# Access to the top level directory of the database
#
simdb = SimDB()

data_dir = os.path.join( simdb.simdb_dir,
                         'simdata', 'input_data_E4_stb',
                         'state_data_lattice_shell' )

file_name = os.path.join( data_dir, 'LC1.csv' )

print '*** read state data from file: %s ***' % ( file_name )

# get the column headings defined in the second row 
# of the .csv-input file
# column_headings = array(['N' 'Vy' 'Vz' 'MT' 'My' 'Mz'])
#
file = open( file_name, 'r' )
lines = file.readlines()
column_headings = lines[1].split( ';' )
column_headings_arr = array( column_headings )
print column_headings_arr[4:10]

# NOTE: the RFEM lattice state data needs to be converted into csv-files using ";" as field delimiter and ""(blank) as text delimiter.
#       The float numbers need to be separated by "." instead of ",". (if necessary use "replace" in text editor)
#       Line numbers in first column must be integers (remove "," and replace by "" (none)
#       In Excel change format of the cells to float with 3 digits in order to save float values properly
#
#arr = loadtxt( file_name, skiprows = 2, delimiter = ";", usecols = [4, 5, 6, 7, 8, 9] )
arr = loadtxt( file_name, skiprows = 2, delimiter = ";", usecols = [0, 4, 5, 6, 7, 8, 9] )
state_data = arr[::14, :]
print state_data
print state_data.shape



#from numpy import array, random
#from enthought.tvtk.api import tvtk
#
## The numpy array data.
#points = array( [[0, -0.5, 0], [1.5, 0, 0], [0, 1, 0], [0, 0, 0.5],
#                [-1, -1.5, 0.1], [0, -1, 0.5], [-1, -0.5, 0],
#                [1, 0.8, 0]], 'f' )
#triangles = array( [[0, 1, 3], [1, 2, 3], [1, 0, 5],
#                   [2, 3, 4], [3, 0, 4], [0, 5, 4], [2, 4, 6],
#                    [2, 1, 7]] )
#scalars = random.random( points.shape )
#
## The TVTK dataset.
#mesh = tvtk.PolyData( points = points, polys = triangles )
#mesh.point_data.scalars = scalars
#mesh.point_data.scalars.name = 'scalars'
#


# Author: Prabhu Ramachandran <prabhu at aero dot iitb dot ac dot in>
# Copyright (c) 2007, Enthought, Inc.
# License: BSD style.

#from numpy import array
#from enthought.tvtk.api import tvtk
#from enthought.mayavi.scripts import mayavi2
#
## The numpy array data.
#points = array( [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], 'f' )
#triangles = array( [[0, 1, 3], [0, 3, 2], [1, 2, 3], [0, 2, 1]] )
#temperature = array( [10., 20., 30., 40.] )
#
## The TVTK dataset.
#mesh = tvtk.PolyData( points = points, polys = triangles )
#mesh.point_data.scalars = temperature
#mesh.point_data.scalars.name = 'Temperature'
#
## Uncomment the next two lines to save the dataset to a VTK XML file.
##w = tvtk.XMLPolyDataWriter(input=mesh, file_name='polydata.vtp')
##w.write()
#
## Now view the data.
#@mayavi2.standalone
#def view():
#    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
#    from enthought.mayavi.modules.surface import Surface
#
#    mayavi.new_scene()
#    src = VTKDataSource( data = mesh )
#    mayavi.add_source( src )
#    s = Surface()
#    mayavi.add_module( s )
#
#if __name__ == '__main__':
#    view()

#from enthought.mayavi import mlab
#import numpy as np
#
### Generate some random data along a straight line in the x-direction
##num = 100
##x = np.arange( num )
##y, z = np.ones( num ), np.ones( num )
##
##s = np.arange( num )
#
#
#arr = np.array( [
#       [0, 0, 0, 1],
#       [1, 2, 3, 22],
#       [5, 5, 5, 5]] )
#
#x = arr[0, :]
#y = arr[1, :]
#z = arr[2, :]
#s = arr[2, :]
#print x
#print y
#print z
#print s
#
#
## Plot using mayavi's mlab api
#fig = mlab.figure()
#
## First we need to make a line source from our data
##line = mlab.pipeline.line_source( x, y, z, s )
##mlab.pipeline.vectors( line )
##mlab.pipeline.glyph.glyph.color_mode = 'color_by_scalar'

#mlab.figure( figure = "SFB532Demo",
#             bgcolor = ( 1.0, 1.0, 1.0 ),
#             fgcolor = ( 0.0, 0.0, 0.0 ) )

#mlab.points3d( x, y, z, s,
#               colormap = "copper",
#               mode = "cube",
#                )
#mlab.outline()
#
##mlab.scalarbar( title = self.plot_column, orientation = 'vertical' )
#
#mlab.show





#x = arr[1:3, 0]
#y = arr[1:3, 1]
#z = arr[1:3, 2]
#s = arr[1:3, 3]
#print x
#
## First we need to make a line source from our data
#line = mlab.pipeline.line_source( x, y, z, s )
#mlab.pipeline.vectors( line )

#vectors = Vectors()
#vtk_data_source = engine.scenes[0].children[0]
#engine.add_filter( vectors, vtk_data_source )
#engine.add_filter( < enthought.mayavi.filters.tube.Tube object at 0x73de650 > , vtk_data_source )


# Then we apply the "tube" filter to it, and vary the radius by "s"
#tube = mlab.pipeline.tube( line, tube_sides = 20, tube_radius = 1.0 )

# Now we display the tube as a surface
#mlab.pipeline.surface( tube )

# And finally visualize the result
#mlab.show()
