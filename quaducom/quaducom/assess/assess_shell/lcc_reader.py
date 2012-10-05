
from etsproxy.traits.api import \
    HasTraits, Directory, \
    Property, WeakRef

import os

import numpy as np

from StringIO import StringIO

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

    lcc_table = WeakRef

class LCCReaderRFEM(LCCReader):

    def read_state_data(self, f_name):

        '''to read the stb-stress resultants save the xls-worksheet
        to a csv-file using ';' as filed delimiter and ' ' (blank)
        as text delimiter.
        '''

        file_name = os.path.join(self.data_dir, f_name)

        print '*** read state data from file: %s ***' % (file_name)

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

    def read_geo_data(self, file_name):
        '''to read the stb - thickness save the xls - worksheet
        to a csv - file using ';' as filed delimiter and ' ' ( blank )
        as text delimiter.
        '''
        print '*** read geo data from file: %s ***' % (file_name)


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

        # coordinates [m]:
        # (NOTE: corrds are taken from the state data file of the first loading case) 
        #
        X = self.lcc_table.lc_list[0].state_data_orig['X']
        Y = self.lcc_table.lc_list[0].state_data_orig['Y']
        Z = self.lcc_table.lc_list[0].state_data_orig['Z']

        return  {'elem_no':elem_no,
                 'X':X, 'Y':Y, 'Z':Z,
                 'thickness':thickness }

class LCCReaderInfoCAD(LCCReader):

    def read_state_data(self, f_name):

        file_name = os.path.join(self.data_dir, 'state_data', f_name)

        print '*** read state data from file: %s ***' % (file_name)

        input_arr = np.loadtxt(file_name)

        elem_no_idx, nx_idx, ny_idx, nxy_idx, mx_idx, my_idx, mxy_idx = range(0, 7)

        # element number:
        #
        elem_no = input_arr[:, elem_no_idx]

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

        return { 'elem_no' : elem_no,
                 'mx' : mx, 'my' : my, 'mxy' : mxy,
                 'nx' : nx, 'ny' : ny, 'nxy' : nxy,
               }

    def read_geo_data(self):
        '''to read the stb - thickness save the xls - worksheet
        to a csv - file using ';' as filed delimiter and ' ' ( blank )
        as text delimiter.
        '''
        geo_dir = os.path.join(self.data_dir, 'geo_data')
        node_file = os.path.join(geo_dir, 'Knotenkoordinaten.txt')
        elem_file = os.path.join(geo_dir, 'Elementbeschreibung.txt')
#        thic_file = os.path.join(geo_dir, 'Querschnittswerte.txt')

        node_arr = np.loadtxt(node_file)

        elem_line_arr = np.loadtxt(elem_file, usecols = (0, 1,), dtype = str)

        elem_no_arr, elem_type_arr = elem_line_arr[:, (0, 1)].T
        t_idx = np.argwhere(elem_type_arr == 'SH36')[:, 0]
        q_idx = np.argwhere(elem_type_arr == 'SH46')[:, 0]

        elem_file_ = open(elem_file, 'r')
        lines = elem_file_.readlines()

        line_arr = np.array(lines)

        t_line_arr = line_arr[t_idx]
        q_line_arr = line_arr[q_idx]

        t_str = StringIO(''.join(t_line_arr))
        q_str = StringIO(''.join(q_line_arr))

        t_elems = np.loadtxt(t_str, usecols = (0, 2, 3, 4, 5), dtype = int)
        q_elems = np.loadtxt(q_str, usecols = (0, 2, 3, 4, 5, 6), dtype = int)

        t_elem_node_map = t_elems[:, 1:] - 1
        q_elem_node_map = q_elems[:, 1:] - 1

        node_idx = np.array(node_arr[:, 0] - 1, dtype = 'int')
        node_X = node_arr[:, 1:][node_idx]

        print 'shape', node_X.shape

        t_elem_node_X = node_X[ t_elem_node_map ]
        q_elem_node_X = node_X[ q_elem_node_map ]

        t_elem_X = np.average(t_elem_node_X, axis = 1)
        q_elem_X = np.average(q_elem_node_X, axis = 1)

        n_X = np.zeros((len(line_arr), 3), dtype = 'float')

        print t_idx.shape
        print n_X.shape
        print t_elem_X.shape

        n_X[t_idx, :] = t_elem_X
        n_X[q_idx, :] = q_elem_X

        X, Y, Z = n_X.T

        #thickness = thic_arr[elem_arr[:, 6]]

        return  {'elem_no':elem_no_arr,
                 'X':X, 'Y':Y, 'Z':Z,
                 #'thickness':thickness 
                 }

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
                            'simdata', 'input_data_barrelshell')

    r = LCCReaderInfoCAD(data_dir = data_dir)

    gd = r.read_geo_data()

    mlab.figure(figure = "SFB532Demo",
                 bgcolor = (1.0, 1.0, 1.0),
                 fgcolor = (0.0, 0.0, 0.0))

    mlab.points3d(gd['X'], gd['Y'], -gd['Z'], 0.5 * np.ones_like(gd['Z']),
                   #colormap = "YlOrBr",
                   mode = "cube",
                   scale_factor = 0.10)

    mlab.show()
