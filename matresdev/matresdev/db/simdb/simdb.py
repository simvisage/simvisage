#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license. The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Mar 29, 2010 by: rch

from etsproxy.traits.api import \
    HasTraits, Property, Str

from os.path import \
    join

from etsproxy.util.home_directory import \
    get_home_directory

class SimDB(HasTraits):
    '''
Basic structure of the database directory.
Implements the relative paths for the three different
categories of data that are managed using the svn server.
- experimental data
- simulation data
- material data
Repository of raw data
- local data
'''
    home_dir = Property
    def _get_home_dir(self):
        return get_home_directory()

    simdb_dir = Property
    def _get_simdb_dir(self):
        return join(self.home_dir, 'simdb')

    simdb_cache_dir = Property
    def _get_simdb_cache_dir(self):
        return join(self.home_dir, '.simdb_cache')

    exdata_dir = Property
    def _get_exdata_dir(self):
        return join(self.simdb_dir, 'exdata')

    matdata_dir = Property
    def _get_matdata_dir(self):
        return join(self.simdb_dir, 'matdata')

    '''
-remote data
'''
    server_username = Str('simdb')

    server_host = Str('mordred.imb.rwth-aachen.de')

    simdb_cache_remote_dir = Str('/home/simdb/simdb/')