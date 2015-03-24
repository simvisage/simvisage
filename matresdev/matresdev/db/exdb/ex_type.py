# -------------------------------------------------------------------------
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
# Created on Feb 15, 2010 by: rch

# @todo: introduce the activation of filters - ironing, smoothing

from etsproxy.traits.api import \
    File, \
    Array, Str, Property, cached_property, \
    Dict, Bool, implements, Float

import ConfigParser

import string

from matresdev.db import SimDBClass

import os

from numpy import \
    loadtxt

from os.path import exists

from loadtxt_novalue import loadtxt_novalue

from string import split

from i_ex_type import \
    IExType

from matresdev.db.simdb import SFTPServer

import zipfile

from matresdev.db.simdb import SimDB
simdb = SimDB()


class ExType(SimDBClass):

    '''Read the data from the directory
    '''

    implements(IExType)

    data_file = File

    file_ext = Str('DAT')

    def validate(self):
        '''Validate the input data return the info whether or not
         the input is valid. This is the condition for processing
         of the derived data.
        '''
        return True

    # set a flag for the view to check whether derived data is available
    #
    derived_data_available = Bool(False)

    # specify inputs
    #
    key = Property(Str, trantient=True, depends_on='data_file')

    def _get_key(self):
        return split(os.path.basename(self.data_file), '.')[0]

    def _set_key(self, value):
        genkey = split(os.path.basename(self.data_file), '.')[0]
        if genkey != value:
            raise KeyError, 'key mismatch %s != %s' % (genkey, value)

    def __setstate__(self, state, kw={}):
        if 'key' in state:
            del state['key']
        super(SimDBClass, self).__setstate__(state, **kw)

    # indicate whether the test is suitable and prepared for
    # calibration.
    ready_for_calibration = Property(Bool)

    def _get_ready_for_calibration(self):
        # return False by default
        # the subclasses shall overload this
        # and define the rules
        return False

    # specify plot templates that can be chosen for viewing
    #
    plot_templates = Dict(transient=True)

    # define processing
    #
    processed_data_array = Array('float_', transient=True)

    def process_source_data(self):
        '''process the source data and assign
        attributes to the DAT-file channel names.
        '''
        print '*** process data ***'
        self._read_data_array()
        self.processed_data_array = self.data_array
        self._set_array_attribs()

    data_array = Array(float, transient=True)

    unit_list = Property(depends_on='data_file')

    def _get_unit_list(self):
        return self.names_and_units[1]

    factor_list = Property(depends_on='data_file')

    def _get_factor_list(self):
        return self.names_and_units[0]

    names_and_units = Property(depends_on='data_file')

    @cached_property
    def _get_names_and_units(self):
        ''' Extract the names and units of the measured data.
        The order of the names in the .DAT-file corresponds
        to the order of the .ASC-file.
        '''
        file_ = open(self.data_file, 'r')
        lines = file_.read().split()
        names = []
        units = []
        for i in range(len(lines)):
            if lines[i] == '#BEGINCHANNELHEADER':
                name = lines[i + 1].split(',')[1]
                unit = lines[i + 3].split(',')[1]
                names.append(name)
                units.append(unit)
        return names, units

    def _names_and_units_default(self):
        ''' Extract the names and units of the measured data.
        The order of the names in the .DAT-file corresponds
        to the order of the .ASC-file.
        '''
        file_ = open(self.data_file, 'r')
        lines = file_.read().split()
        names = []
        units = []
        for i in range(len(lines)):
            if lines[i] == '#BEGINCHANNELHEADER':
                name = lines[i + 1].split(',')[1]
                unit = lines[i + 3].split(',')[1]
                names.append(name)
                units.append(unit)
        print 'names, units', names, units
        return names, units

    def _set_array_attribs(self):
        '''Set the measured data as named attributes defining slices into
        the processed data array.
        '''
        for i, factor in enumerate(self.factor_list):
            self.add_trait(
                factor, Array(value=self.processed_data_array[:, i],
                              transient=True))

    # ------------------

    def _read_data_array(self):
        ''' Read the experiment data.
        '''
        if exists(self.data_file):

            print 'READ FILE'
            # change the file name dat with asc
            file_split = self.data_file.split('.')

            file_name = file_split[0] + '.csv'
            if not os.path.exists(file_name):

                file_name = file_split[0] + '.ASC'
                if not os.path.exists(file_name):
                    raise IOError, 'file %s does not exist' % file_name

            print 'file_name', file_name

            # try to use loadtxt to read data file
            try:
                _data_array = loadtxt(file_name,
                                      delimiter=';')

            # loadtxt returns an error if the data file contains
            # 'NOVALUE' entries. In this case use the special
            # method 'loadtxt_novalue'
            except ValueError:
                _data_array = loadtxt_novalue(file_name)

            self.data_array = _data_array

    data_dir = Property()
    '''Local directory path of the data file.
    '''

    def _get_data_dir(self):
        return os.path.dirname(self.data_file)

    relative_path = Property
    '''Relative path inside database structure - the path is same for experiment
    in both database structures (remote and local)
    '''

    def _get_relative_path(self):
        return self.data_dir.replace(simdb.simdb_dir, '')[1:]

    hook_up_file = Property
    '''File specifying the access to extended data.
    The cfg file is used to hook up arbitrary type
    of data stored anywhere that can be downloaded
    on demand to the local cache.
    '''

    def _get_hook_up_file(self):
        dir_path = os.path.dirname(self.data_file)
        file_name = os.path.basename(self.data_file)
        file_split = file_name.split('.')
        file_name = os.path.join(dir_path,
                                 file_split[0] + '.cfg')
        if not os.path.exists(file_name):
            file_name = ''
        return file_name

    aramis_start_offset = Property(Float, depends_on='data_file')
    '''Get time offset of aramis start specified in the hookup file.
    '''
    @cached_property
    def _get_aramis_start_offset(self):
        # hook_up an extended file if available.
        aramis_start_offset = 0.0
        if self.hook_up_file:
            config = ConfigParser.ConfigParser()
            config.read(self.hook_up_file)
            try:
                aramis_start_offset = config.get('aramis_data',
                                                 'aramis_start_offset')
            except ConfigParser.NoOptionError:
                pass
        return float(aramis_start_offset)

    aramis_files = Property(depends_on='data_file')
    '''Get the list of available aramis files specified in the hookup file.
    '''
    @cached_property
    def _get_aramis_files(self):
        # hook_up an extended file if available.
        aramis_files = []
        if self.hook_up_file:
            config = ConfigParser.ConfigParser()
            config.read(self.hook_up_file)
            aramis_files = config.get(
                'aramis_data', 'aramis_files').split(',\n')
        return aramis_files

    aramis_dict = Property(depends_on='data_file')
    '''Use the last two specifiers of the aramis file name
    as a key to access the proper file.
    '''
    @cached_property
    def _get_aramis_dict(self):
        # hook_up an extended file if available.
        af_dict = {}
        for af in self.aramis_files:
            fx, fy = af.split('-')[-2:]
            af_dict[fx + '-' + fy] = af
        return af_dict

    def download_aramis_file(self, arkey):
        af = self.aramis_dict[arkey]
        af_rel_dir = os.path.join(self.relative_path, 'aramis')
        af_local_dir = os.path.join(simdb.simdb_cache_dir, af_rel_dir)
        if not os.path.exists(af_local_dir):
            os.makedirs(af_local_dir)
        try:
            s = SFTPServer(simdb.server_username, '', simdb.server_host)
            if hasattr(s, 'sftp'):
                zip_filename = af + '.zip'
                zipfile_server = os.path.join(
                    simdb.simdb_cache_remote_dir, af_rel_dir, zip_filename)

                zipfile_server = string.replace(zipfile_server, '\\', '/')
                zipfile_local = os.path.join(af_local_dir, zip_filename)

                print 'downloading', zipfile_server
                print 'destination', zipfile_local

                s.download(zipfile_server, zipfile_local)
                s.sftp.stat(zipfile_server)
                s.close()
        except IOError, e:
            raise IOError(e)

    def uncompress_aramis_file(self, arkey):
        af = self.aramis_dict[arkey]
        af_rel_dir = os.path.join(self.relative_path, 'aramis')
        af_local_dir = os.path.join(simdb.simdb_cache_dir, af_rel_dir)
        zip_filename = af + '.zip'
        zipfile_local = os.path.join(af_local_dir, zip_filename)
        if not os.path.exists(zipfile_local):
            self.download_aramis_file(arkey)

        print 'uncompressing'
        zf = zipfile.ZipFile(zipfile_local, 'r')
        zf.extractall(af_local_dir)
        zf.close()

    def get_cached_aramis_file(self, arkey):
        '''For the specified aramis resolution key check if the file
        has already been downloaded.
        '''
        af = self.aramis_dict.get(arkey, None)
        if af is None:
            print 'Aramis data not available for resolution %s of the'\
                'test data\n%s' % (arkey, self.data_file)
            return None
        af_path = os.path.join(
            simdb.simdb_cache_dir, self.relative_path, 'aramis', af)

        if not os.path.exists(af_path):
            self.uncompress_aramis_file(arkey)

        return af_path
