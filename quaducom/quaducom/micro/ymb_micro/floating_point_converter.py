'''
Created on Dec 13, 2010

@author: kelidas
'''
'''
 Replace all commas by points in all *.txt files in the given folder (data_dir).  
'''

from os.path import join
from promod.simdb import SimDB
from string import replace
import os

simdb = SimDB()
#data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'VET', 'raw_data' )
#data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'MAG', 'raw_data' )
data_dir = join( simdb.exdata_dir, 'trc', 'yarn_structure', 'TEST' )

def fcomma2fdot( x ): return replace( x, ',', '.' )

paths = os.listdir( data_dir )  # list of paths in that dir

for file in paths:
    print 'reading %s' % file 
    if file.endswith( '.txt' ):
        print os.path.join( data_dir, file )
        infile = open( os.path.join( data_dir, file ), 'r' )
        data = infile.read()
        #print data
        data = fcomma2fdot( data )    
        infile.close()
        outfile = open( os.path.join( data_dir, file ), 'w' )
        outfile.write( data )
        outfile.close()



