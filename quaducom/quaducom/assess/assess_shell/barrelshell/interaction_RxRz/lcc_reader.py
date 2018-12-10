
from traits.api import \
    HasTraits, Directory, \
    Property, WeakRef

from numpy import array, random
from tvtk.api import tvtk

import os

import numpy as np

from io import StringIO

class LCCReader(HasTraits):
    '''Base class for LCC Readers.'''

    _dd = Directory

    data_dir = Property
    def _get_data_dir(self):
        if self._dd:
            return self._dd
        else:
            return self.lcc_table.data_dir
    def _set_data_dir(self, data_dir):
        self._dd = data_dir

    def check_for_consistency(self, lc_list, geo_data_dict):
        return True

    lcc_table = WeakRef

class LCCReaderRFEM(LCCReader):

    def read_state_data(self, f_name):

        '''to read the stb-stress resultants save the xls-worksheet
        to a csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''

        file_name = os.path.join(self.data_dir, f_name)

        print('*** read state data from file: %s ***' % (file_name))

        # get the column headings defined in the second row 
        # of the csv soliciotations input file
        # column_headings = np.array(["Nr.","Punkt","X","Y","Z","mx","my","mxy","vx","vy","nx","ny","nxy"])
        #
        file_ = open(file_name, 'r')
        lines = file_.readlines()
        column_headings = lines[1].split(';')
        # remove '\n' from last string element in list
        #
        column_headings[-1] = column_headings[-1][:-1]
        column_headings_arr = np.array(column_headings)

        elem_no_idx = np.where('Nr.' == column_headings_arr)[0]
        X_idx = np.where('X' == column_headings_arr)[0]
        Y_idx = np.where('Y' == column_headings_arr)[0]
        Z_idx = np.where('Z' == column_headings_arr)[0]
        mx_idx = np.where('mx' == column_headings_arr)[0]
        my_idx = np.where('my' == column_headings_arr)[0]
        mxy_idx = np.where('mxy' == column_headings_arr)[0]
        nx_idx = np.where('nx' == column_headings_arr)[0]
        ny_idx = np.where('ny' == column_headings_arr)[0]
        nxy_idx = np.where('nxy' == column_headings_arr)[0]

        file_.close()

        # define np.arrays containing the information from the raw input file
        #
        input_arr = np.loadtxt(file_name , delimiter = ';', skiprows = 2)

        # element number:
        #
        elem_no = input_arr[:, elem_no_idx]

        # coordinates [m]:
        #
        X = input_arr[:, X_idx]
        Y = input_arr[:, Y_idx]
        Z = input_arr[:, Z_idx]

        # moments [kNm/m]
        #
        mx = input_arr[:, mx_idx]
        my = input_arr[:, my_idx]
        mxy = input_arr[:, mxy_idx]

        # normal forces [kN/m]:
        #
        nx = input_arr[:, nx_idx]
        ny = input_arr[:, ny_idx]
        nxy = input_arr[:, nxy_idx]

        return { 'elem_no' : elem_no, 'X' : X, 'Y' : Y, 'Z' : Z,
                 'mx' : mx, 'my' : my, 'mxy' : mxy,
                 'nx' : nx, 'ny' : ny, 'nxy' : nxy,
               }

    # stress resultants to be multiplied within the LCC combinations
    #
    sr_columns = ['mx', 'my', 'mxy', 'nx', 'ny', 'nxy']

    def read_geo_data(self, file_name):
        '''to read the stb - thickness save the xls - worksheet
        to a csv - file using ';' as filed delimiter and ' ' ( blank )
        as text delimiter.
        '''
        print('*** read geo data from file: %s ***' % (file_name))


        # coordinates [m]:
        # (NOTE: corrds are taken from the state data file of the first loading case) 

        # the column headings are defined in the first/second row 
        # of the csv thickness input file
        # Flaeche;;;Material;Dicke;;Exzentrizitaet;Integrierte Objekte;;;

        # (NOTE: the headings contain non-ascii characters. Therefore the
        #       column indices can not be derived using the 'np.where'-method.)

        # read the float data:
        #
        input_arr = np.loadtxt(file_name, usecols = (0, 5), delimiter = ';', skiprows = 2)
        elem_no_idx = 0
        thickness_idx = 1

        # element number:
        # (NOTE: column np.array must be of shape (n_elems, 1)
        #
        elem_no = input_arr[:, elem_no_idx][:, None]

        # element thickness [mm]:
        # (NOTE: column np.array must be of shape (n_elems, 1)
        #
        thickness = input_arr[:, thickness_idx][:, None]
        
        # convert thickness to [m]
        thickness = thickness / 1000.

        # coordinates [m]:
        # (NOTE: corrds are taken from the state data file of the first loading case) 
        #
        X = self.lcc_table.lc_list[0].state_data_orig['X']
        Y = self.lcc_table.lc_list[0].state_data_orig['Y']
        Z = self.lcc_table.lc_list[0].state_data_orig['Z']

        return  {'elem_no':elem_no,
                 'X':X, 'Y':Y, 'Z':Z,
                 'thickness':thickness }

    def plot_col(self, mlab, plot_col, geo_data, state_data, warp_factor = 1.):
        '''
        plot the chosen plot_col array at the center of gravity of the elements;
        method is used by 'ls_table' to plot the selected plot variable
        NOTE: if RFEM-Reader is used no values for the deformed state are read in yet
        '''
        gd = geo_data

        # element coordinates of the undeformed shape 
        # (2d column arrays)
        #
        X = gd['X'].flatten()
        Y = gd['Y'].flatten()
        # switch orientation of the z-axis
        Z = (-1.0) * gd['Z'].flatten()

        # plot state data in the deformed geometry  
        #
        mlab.points3d(X, Y, Z, plot_col,
                           colormap = "YlOrBr",
                           mode = "cube",
                           scale_mode = 'none',
                           scale_factor = 0.15)


    def check_for_consistency(self, lc_list, geo_data_dict):

        for lc in lc_list:
            # check internal LC-consitency: 
            # (compare coords-values of first LC with all other LC's in 'lc_list')
            #
            if not all(lc_list[0].state_data_dict['X'] == lc.state_data_dict['X']) and \
                not all(lc_list[0].state_data_dict['Y'] == lc.state_data_dict['Y']) and \
                not all(lc_list[0].state_data_dict['Z'] == lc.state_data_dict['Z']):
                raise ValueError("coordinates in loading case '%s' and loading case '%s' are not identical. Check input files for consistency!" \
                        % (self.lc_list[0].name, lc.name))
                return False

            # check external consistency:
            # (compare 'elem_no' in 'thickness.csv' and in all loading-cases 
            # input files (e.g. 'LC1.csv') defined in 'lc_list')
            #
            if not all(geo_data_dict['elem_no'] == lc.state_data_dict['elem_no']):
                raise ValueError("element numbers in loading case '%s' and loading case '%s' are not identical. Check input files for consistency!" \
                        % (lc_list[0].name, lc.name))
                return False

        print('*** input files checked for consistency (OK) ***')
        return True


class LCCReaderInfoCAD(LCCReader):

    def read_state_data(self, f_name):

        file_name = os.path.join(self.data_dir, 'state_data',
                                 'Flaechenschnittgroessen',
                                 'Schwerpunkt', f_name)

        print('*** read state data from file: %s ***' % (file_name))

        input_arr = np.loadtxt(file_name)

        elem_no_idx, nx_idx, ny_idx, nxy_idx, mx_idx, my_idx, mxy_idx = list(range(0, 7))

        # element number:
        #
        elem_no = input_arr[:, [elem_no_idx]]
        
        # moments [kNm/m]
        #
        mx = input_arr[:, [mx_idx]]
        my = input_arr[:, [my_idx]]
        mxy = input_arr[:, [mxy_idx]]

        # normal forces [kN/m]:
        #
        nx = input_arr[:, [nx_idx]]
        ny = input_arr[:, [ny_idx]]
        nxy = input_arr[:, [nxy_idx]]

        file_name = os.path.join(self.data_dir, 'state_data',
                                 'Knotendeformationen',
                                 f_name)

        print('*** read deformation data from file: %s ***' % (file_name))

        input_arr = np.loadtxt(file_name)

        node_no_idx, ux_idx, uy_idx, uz_idx, phix_idx, phiy_idx, phiz_idx = list(range(0, 7))

        ux = input_arr[:, [ux_idx]]
        uy = input_arr[:, [uy_idx]]
        uz = input_arr[:, [uz_idx]]

        # nodal displacements [m] (arranged as array with 3 columns)
        #
        node_U = input_arr[:, [ux_idx, uy_idx, uz_idx]]

        #-----------
        # get the element displacement at the center of gravity
        # the calculation based on the nodal displacement needs the information stored in 'geo_data'
        #-----------

        # the path to the 'geo_data' files is specified specifically in the definition of 'read_geo_data'
        gd = self.read_geo_data( f_name )

        # get mapping from 'geo_data'
        #
        t_elem_node_map = gd['t_elem_node_map']
        q_elem_node_map = gd['q_elem_node_map']
        t_idx = gd['t_idx']
        q_idx = gd['q_idx']
        
        # average element displacements (unordered)
        #
        t_elem_node_U = node_U[ t_elem_node_map ]
        q_elem_node_U = node_U[ q_elem_node_map ]
        t_elem_U = np.average(t_elem_node_U, axis = 1)
        q_elem_U = np.average(q_elem_node_U, axis = 1)
    
        # average element displacements (ordered in ascending element number)
        #
        elem_U = np.zeros((len(t_elem_U)+len(q_elem_U), 3), dtype = 'float')
        elem_U[t_idx, :] = t_elem_U
        elem_U[q_idx, :] = q_elem_U
    
        # average element displacements stored in 1d column arrays
        #
        ux_elem = elem_U[:,0,None]
        uy_elem = elem_U[:,1,None]
        uz_elem = elem_U[:,2,None]

        return { 'elem_no' : elem_no,
                 'mx' : mx, 'my' : my, 'mxy' : mxy,
                 'nx' : nx, 'ny' : ny, 'nxy' : nxy,
                 'ux' : ux, 'uy' : uy, 'uz' : uz,
                 'node_U' : node_U,
                 #-------------------
                 'ux_elem' : ux_elem,
                 'uy_elem' : uy_elem,
                 'uz_elem' : uz_elem
               }

    # stress resultants to be multiplied within the LCC combinations
    #
    sr_columns = ['mx', 'my', 'mxy', 'nx', 'ny', 'nxy', 'ux_elem', 'uy_elem', 'uz_elem']

    def read_geo_data(self, file):
        '''to read the thickness file exported from InfoCAD
        using 'tab' as filed delimiter.
        '''
        geo_dir = os.path.join(self.data_dir, 'geo_data')
        node_file = os.path.join(geo_dir, 'Knotenkoordinaten.txt')
        elem_file = os.path.join(geo_dir, 'Elementbeschreibung.txt')
        thic_file = os.path.join(geo_dir, 'Querschnittswerte.txt')

        print('*** read geo data from node file: %s ***' % (node_file))
        print('*** read geo data from elem file: %s ***' % (elem_file))
        print('*** read geo data from thic files: %s ***' % (thic_file))

        node_arr = np.loadtxt(node_file)

        elem_line_arr = np.loadtxt(elem_file, usecols = (0, 1,), dtype = str)

        elem_no_arr, elem_type_arr = elem_line_arr[:, (0, 1)].T
        t_idx = np.argwhere(elem_type_arr == 'SH36')[:, 0]
        q_idx = np.argwhere(elem_type_arr == 'SH46')[:, 0]

        elem_file_ = open(elem_file, 'r')
        lines = elem_file_.readlines()
        line_arr = np.array(lines)
        elem_file_.close()

        t_line_arr = line_arr[t_idx]
        q_line_arr = line_arr[q_idx]

        t_str = StringIO(''.join(t_line_arr))
        q_str = StringIO(''.join(q_line_arr))

        t_elems = np.loadtxt(t_str, usecols = (0, 2, 3, 4, 5), dtype = int)
        q_elems = np.loadtxt(q_str, usecols = (0, 2, 3, 4, 5, 6), dtype = int)

        t_elem_node_map = t_elems[:, 1:-1] - 1
        q_elem_node_map = q_elems[:, 1:-1] - 1
        t_thickness_idx = t_elems[:, -1] - 1
        q_thickness_idx = q_elems[:, -1] - 1

        node_idx = np.array(node_arr[:, 0] - 1, dtype = 'int')
        node_X = node_arr[:, 1:][node_idx]

        t_elem_node_X = node_X[ t_elem_node_map ]
        q_elem_node_X = node_X[ q_elem_node_map ]

        t_elem_X = np.average(t_elem_node_X, axis = 1)
        q_elem_X = np.average(q_elem_node_X, axis = 1)

        n_X = np.zeros((len(line_arr), 3), dtype = 'float')

        n_X[t_idx, :] = t_elem_X
        n_X[q_idx, :] = q_elem_X

        X, Y, Z = n_X.T

        thic_file_ = open(thic_file, 'r')
        lines = thic_file_.readlines()
        line_arr = np.array(lines)
        thic_file_.close()

        idx_list = []
        d_list = []
        nr_lines = list(range(0, len(lines), 2))
        for line in line_arr[nr_lines]:
            nr, idx_str, id, d_str = line.split()
            idx = int(idx_str)
            d = float(d_str.split('=')[1])
            idx_list.append(idx)
            d_list.append(d)

        d_arr = np.array(d_list, dtype = 'f')
        d_arr = d_arr[np.array(idx_list, dtype = 'int_') - 1]

        thickness_arr = np.zeros((elem_line_arr.shape[0],), dtype = 'f')
        thickness_arr[t_idx] = d_arr[ t_thickness_idx ]
        thickness_arr[q_idx] = d_arr[ q_thickness_idx ]

        # convert strings entries to floats
        #
        elem_no_arr = array( elem_no_arr, dtype = float )

        # convert 1d-arrays to 2d-column arrays
        #
        elem_no_arr = elem_no_arr[:,np.newaxis]

        X = X[:,np.newaxis]
        Y = Y[:,np.newaxis]
        Z = Z[:,np.newaxis]
        thickness_arr = thickness_arr[:,np.newaxis]

        return  {'elem_no':elem_no_arr,
                 'X':X, 'Y':Y, 'Z':Z,
                 'node_X' : node_X,
                 't_elem_node_map' : t_elem_node_map,
                 'q_elem_node_map' : q_elem_node_map,
                 't_idx' : t_idx,
                 'q_idx' : q_idx,
                 'thickness':thickness_arr
                 }

    def plot_mesh(self, mlab, geo):

        points = geo['node_X'] 
        triangles = geo['t_elem_node_map']
        quads = geo['q_elem_node_map']
        #scalars = random.random(points.shape)

        # The TVTK dataset.
        qmesh = tvtk.PolyData(points = points, polys = quads)
        mlab.pipeline.surface(qmesh, representation = 'wireframe')

        # The TVTK dataset.
        tmesh = tvtk.PolyData(points = points, polys = triangles)
        mlab.pipeline.surface(tmesh, representation = 'wireframe')

    def plot_deformed_mesh(self, mlab, geo_data, state_data = {'node_U' : array([0., 0., 0.])}, warp_factor = 1.0):
        '''plot the deformed mesh based on the nodal displacement defined in 'state_data'
        '''
        points = geo_data['node_X'] 
        node_U = state_data['node_U']

        node_U_warped = node_U * warp_factor 
        points += node_U_warped

        triangles = geo_data['t_elem_node_map']
        quads = geo_data['q_elem_node_map']

        # The TVTK dataset.
        qmesh = tvtk.PolyData(points = points, polys = quads)
        mlab.pipeline.surface(qmesh, representation = 'wireframe')

        # The TVTK dataset.
        tmesh = tvtk.PolyData(points = points, polys = triangles)
        mlab.pipeline.surface(tmesh, representation = 'wireframe')

    def plot_sd(self, mlab, geo_data, sd_key, state_data = {'node_U' : array([0., 0., 0.])}, warp_factor = 1.0):
        '''plot the chosen state data defined by 'sd_key' at the center of gravity of the elements
        together with the element mesh; 'warp_factor' can be used to warp the deformation state in the plot.
        '''
        gd = geo_data
        sd = state_data

        # plot the deformed geometry (as mesh)
        #
        self.plot_deformed_mesh(mlab, gd, state_data = sd, warp_factor = warp_factor)

        # get mapping from 'geo_data'
        #
        t_elem_node_map = gd['t_elem_node_map']
        q_elem_node_map = gd['q_elem_node_map']
        t_idx = gd['t_idx']
        q_idx = gd['q_idx']
        
        # nodal displacement 
        #
        node_U = sd['node_U']

        # average element displacement (unordered)
        #
        t_elem_node_U = node_U[ t_elem_node_map ]
        q_elem_node_U = node_U[ q_elem_node_map ]
        t_elem_U = np.average(t_elem_node_U, axis = 1)
        q_elem_U = np.average(q_elem_node_U, axis = 1)
    
        # element displacement (ordered in ascending element number)
        #
        elem_U = np.zeros((len(t_elem_U)+len(q_elem_U), 3), dtype = 'float')
        elem_U[t_idx, :] = t_elem_U
        elem_U[q_idx, :] = q_elem_U
    
        # element coordinates of the undeformed shape 
        # (2d column arrays)
        #
        X = gd['X']
        Y = gd['Y']
        Z = gd['Z']
    
        # average element deformations 
        #
        ux_elem = elem_U[:,0,None]
        uy_elem = elem_U[:,1,None]
        uz_elem = elem_U[:,2,None]
    
        # element coordinates of the deformed state 
        # considering the specified warp factor 
        #
        X_def = X + ux_elem * warp_factor
        Y_def = Y + uy_elem * warp_factor
        Z_def = Z + uz_elem * warp_factor
        
        # plot state data in the deformed geometry  
        #
        mlab.points3d(X_def, Y_def, Z_def, sd[sd_key],
                      mode = "cube")

    def plot_col(self, mlab, plot_col, geo_data,  
                 state_data = {'ux_elem' : array([[0.], [0.], [0.]]),
                               'uy_elem' : array([[0.], [0.], [0.]]), 
                               'uz_elem' : array([[0.], [0.], [0.]])}, 
                 warp_factor = 1.0):
        '''
        plot the chosen plot_col array at the center of gravity of the elements;
        method is used by 'ls_table' to plot the selected plot variable
        ('warp_factor' can be used to warp the deformation state in the plot).
        '''
        gd = geo_data
        sd = state_data

        # element coordinates of the undeformed shape 
        # (2d column arrays)
        #
        X = gd['X']
        Y = gd['Y']
        Z = gd['Z']
    
        # average element deformations 
        #
        ux_elem = sd['ux_elem']
        uy_elem = sd['uy_elem']
        uz_elem = sd['uz_elem']
    
        # element coordinates of the deformed state 
        # considering the specified warp factor 
        #
        X_def = X + ux_elem * warp_factor
        Y_def = Y + uy_elem * warp_factor
        Z_def = Z + uz_elem * warp_factor
        
        X_def = X_def.flatten()
        Y_def = Y_def.flatten()
        Z_def = Z_def.flatten()

        # plot state data in the deformed geometry  
        #
        mlab.points3d(X_def, Y_def, Z_def, plot_col,
                      mode = "cube",
                      scale_mode = 'none',
                      scale_factor = 0.05)
    
    def check_for_consistency(self, lc_list, geo_data_dict):
        print('*** check for consistency ***')

        for lc in lc_list:
            # check internal LC-consitency: 
            # (compare elem_no of first LC with all other LC's in 'lc_list')
            #
            if not all(lc_list[0].state_data_dict['elem_no'] == lc.state_data_dict['elem_no']):
                raise ValueError("element numbers in loading case '%s' and loading case '%s' are not identical. Check input files for internal consistency!" \
                        % (self.lc_list[0].name, lc.name))
                return False

            # check external consistency:
            # (compare 'elem_no' in 'geo_data' and 'elem_no' in state data of all loading-cases 
            # input files (e.g. 'LC1.txt') defined in 'lc_list')
            #
            if not all(geo_data_dict['elem_no'] == lc.state_data_dict['elem_no']):
                raise ValueError("element numbers in loading case '%s' and loading case '%s' are not identical. Check input files for external consistency!" \
                        % (lc_list[0].name, lc.name))
                return False

        print('*** input files checked for consistency (OK) ***')
        return True

class LCCReaderInfoCADRxyz(LCCReader):

    def read_state_data(self, f_name):

        file_name = os.path.join(self.data_dir, 'state_data',\
                                 'Auflagerreaktionen',\
                                 f_name)

        print('*** read state data from file: %s ***' % (file_name))

        input_arr = np.loadtxt(file_name)

        node_no_idx, Rx_idx, Ry_idx, Rz_idx, Mx_idx, My_idx, Mz_idx = list(range(0, 7))

        # node number:
        #
        node_no = input_arr[:, [node_no_idx]]

        # forces [kN]:
        #
        Rx = input_arr[:, [Rx_idx]]
        Ry = input_arr[:, [Ry_idx]]
        Rz = input_arr[:, [Rz_idx]]

        # moments [kNm]
        #
        Mx = input_arr[:, [Mx_idx]]
        My = input_arr[:, [My_idx]]
        Mz = input_arr[:, [Mz_idx]]

        return { 'node_no' : node_no,
                 'Rx' : Rx, 'Ry' : Ry, 'Rz' : Rz,
                 'Mx' : Mx, 'My' : My, 'Mz' : Mz }

    def read_geo_data(self, f_name):
        '''read the thickness file exported from InfoCAD
        using 'tab' as filed delimiter.
        '''
        sd = self.read_state_data( f_name )
        node_no = sd['node_no']
        
        geo_dir = os.path.join(self.data_dir, 'geo_data')
        node_file = os.path.join(geo_dir, 'Knotenkoordinaten.txt')
        print('*** read support node file: %s ***' % (node_file))
        node_arr = np.loadtxt(node_file)

        idx_spprt_nodes = [ np.where( node_arr[:,0] == node_no[i] )[0] for i in range( node_no.shape[0] ) ]
        print('idx_spprt_nodes', idx_spprt_nodes)
        X_spprt = node_arr[ idx_spprt_nodes, 1 ]
        Y_spprt = node_arr[ idx_spprt_nodes, 2 ]
        Z_spprt = node_arr[ idx_spprt_nodes, 3 ]
        
        sd['X_spprt'] = X_spprt
        sd['Y_spprt'] = Y_spprt
        sd['Z_spprt'] = Z_spprt
        return sd

    
    def check_for_consistency(self, lc_list, geo_data_dict):
        pass



if __name__ == '__main__':

    from matresdev.db.simdb import \
        SimDB

    import os

    from etsproxy.mayavi import \
        mlab

    from lcc_table import LCCTableULS, LC

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    #---------------------------------------------
    # 2 shells: 
    # new geometry with new loading cases plus waterfilling
    #---------------------------------------------

    data_dir = os.path.join(simdb.simdb_dir,
                            'simdata',
                            'input_data_barrelshell',
                            '2cm'
#                            '3cm'
#                            '3-4cm'
                            )

    r = LCCReaderInfoCADRxyz(data_dir = data_dir)

    mlab.figure(figure = "SFB532Demo",
                bgcolor = (1.0, 1.0, 1.0),
                fgcolor = (0.0, 0.0, 0.0))

    # get 'geo_data' and 'state_data'
    #
#    gd = r.read_geo_data('')
    sd = r.read_state_data('LC1.txt')

    # plot the undeformed geometry (as mesh)
    #
#    r.plot_mesh(mlab, gd)

    # plot the deformed geometry (as mesh)
    #
#    r.plot_deformed_mesh(mlab, gd, sd, warp_factor = 1000.)

    # plot state data in the deformed geometry (as mesh)
    #
#    r.plot_sd(mlab, gd, 'mx', state_data = sd, warp_factor = 1000.)

    
#    # get mapping from 'geo_data'
#    #
#    t_elem_node_map = gd['t_elem_node_map']
#    q_elem_node_map = gd['q_elem_node_map']
#    t_idx = gd['t_idx']
#    q_idx = gd['q_idx']
##    print 't_elem_node_map',t_elem_node_map.shape
##    print 'q_elem_node_map',q_elem_node_map.shape
#
#    
#    # nodal displacement 
#    #
#    node_U = sd['node_U']
#    
#    t_elem_node_U = node_U[ t_elem_node_map ]
#    q_elem_node_U = node_U[ q_elem_node_map ]
#    
#    t_elem_U = np.average(t_elem_node_U, axis = 1)
#    q_elem_U = np.average(q_elem_node_U, axis = 1)
#
#    # element displacement (ordered)
#    #
#    elem_U = np.zeros((len(t_elem_U)+len(q_elem_U), 3), dtype = 'float')
#    
#    elem_U[t_idx, :] = t_elem_U
#    elem_U[q_idx, :] = q_elem_U
#
#    # plot the deformed geometry (as mesh)
#    #
#    warp_factor = 2000.
#    r.plot_deformed_mesh(mlab, gd, sd, warp_factor = warp_factor)
#
#    node_X = gd['node_X'] 
#
#    # element coordinates of the undeformed shape 
#    # (2d column arrays)
#    #
#    X = gd['X']
#    Y = gd['Y']
#    Z = gd['Z']
#
#    # element coordinates of the deformed shape  
#    # considering 'warp_factor' 
#    #
#    ux_elem = elem_U[:,0,None]
#    uy_elem = elem_U[:,1,None]
#    uz_elem = elem_U[:,2,None]
#
#    # deformed element coordinates with 
#    # increase deformation by warp factor 
#    #
#    X_def = X + ux_elem * warp_factor
#    Y_def = Y + uy_elem * warp_factor
#    Z_def = Z + uz_elem * warp_factor
#    
#    # plot state data in the deformed geometry  
#    #
#    points_pipeline = mlab.points3d(X_def, Y_def, Z_def, sd['mx'],
#                                    mode = "cube",
#                                    scale_factor = 0.2,
#                                    scale_mode = 'none')
##    src = points_pipeline.mlab_source
##    warp = src.vectors = np.hstack([ux_elem, Y_def, 2.0 + 2*Z + uz_elem])

#    mlab.show()

    #-------
    #@todo: use mayavi functionality using pipeline!
#    src = points_pipeline.mlab_source
#    warp = src.vectors = elem_U
#    mlab.show_pipeline()
#    glyph.glyph.glyph.vector_mode = 'use_normal'
#    engine.add_filter(<enthought.mayavi.filters.warp_vector.WarpVector object at 0x8dc35f0>, vtk_data_source2)
#    enthought.mayavi.modules.glyph.Glyph

