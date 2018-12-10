'''
Created on Jul 7, 2010

@author: alexander
'''

from .hp_shell import HPShell
from ibvpy.mesh.fe_grid import FEGrid
from ibvpy.fets.fets3D.fets3D8h import FETS3D8H
from numpy import savetxt, array

if __name__ == '__main__':

    hp_shell = HPShell()

    fets = FETS3D8H()



    dens_xy = [10, 20, 50]
    dens_z = [1, 3]

    for xy in dens_xy:
        for z in dens_z:
            fe_grid = FEGrid( coord_min = ( 0.0, 0.0, 0.0 ),
                              coord_max = ( 1.0, 1.0, 1.0 ),
                              geo_transform = hp_shell,
                              shape = ( xy, xy, z ),
                              fets_eval = fets
                              )

            name = '%dx%dx%d.dat' % ( xy, xy, z )
            savetxt( 'nodes_' + name, fe_grid.geo_grid.cell_grid.points )
            savetxt( 'elems_' + name, fe_grid.geo_grid.cell_grid.cell_node_map, fmt = '%d' )

