'''
Created on Jan 11, 2013

@author: alexander
'''

if __name__ == '__main__':

    from matresdev.db.simdb import \
        SimDB
    
    import numpy as np

    import os
    import string 

    # Access to the top level directory of the database
    #
    simdb = SimDB()

    dir_list = [
#                            '2cm-feines-Netz', 
#                            '2cm-Elastomer', 
                            '2cm-feines-Netz-EC1', 
#                            '2cm', 
#                            '2cm-Zugausfall-komplett', 
#                            '2cm-Zugausfall-vorne', 
#                            '2cm-Zugausfall-hinten', 
#                            '2cm-Zugausfall-links', 
#                            '2cm-Zugausfall-rechts',
#                            '2cm-Zugausfall-hinten-links',
#                            '2cm-Zugausfall-vorne-links',
#                            '2cm-Zugausfall-vorne-rechts',
#                            '2cm-Zugausfall-hinten-rechts',
#                            '2cm-Zugausfall-nie-hinten-links',
#                            '2cm-Zugausfall-nie-vorne-links',
#                            '2cm-Zugausfall-nie-vorne-rechts',
#                            '2cm-Zugausfall-nie-hinten-rechts',
#                            '2cm-Zugausfall-hinten-links-vorne-rechts',
#                            '2cm-Zugausfall-vorne-links-hinten-rechts'
               ]
   
    
    def fcomma2fdot( x ): 
        return string.replace( x, ',', '.' )

    for dir_name in dir_list:
        data_dir = os.path.join(simdb.simdb_dir,
                    'simdata', 
                    'input_data_barrelshell',
                    dir_name) 

        data_dir_ = os.path.join(data_dir, 'state_data', 'Auflagerreaktionen')
        print 'data_dir_', data_dir_
        filenames = os.listdir( data_dir_ )  # list of filenames in that dir
        print 'filenames', filenames
    
        for filename in filenames:
            print 'reading %s' % filename 
            if filename.endswith( '.txt' ):
                print os.path.join( data_dir_, filename )
                infile = open( os.path.join( data_dir_, filename ), 'r' )
                data = infile.read()
                #print data
                data = fcomma2fdot( data )    
                infile.close()
                outfile = open( os.path.join( data_dir_, filename ), 'w' )
                outfile.write( data )
                outfile.close()

        data_dir_ = os.path.join(data_dir, 'state_data', 'Flaechenschnittgroessen', 'Schwerpunkt')
        print 'data_dir_', data_dir_
        filenames = os.listdir( data_dir_ )  # list of filenames in that dir
        print 'filenames', filenames
    
        for filename in filenames:
            print 'reading %s' % filename 
            if filename.endswith( '.txt' ):
                print os.path.join( data_dir_, filename )
                infile = open( os.path.join( data_dir_, filename ), 'r' )
                data = infile.read()
                #print data
                data = fcomma2fdot( data )    
                infile.close()
                outfile = open( os.path.join( data_dir_, filename ), 'w' )
                outfile.write( data )
                outfile.close()
                
        data_dir_ = os.path.join(data_dir, 'state_data', 'Knotendeformationen')
        print 'data_dir_', data_dir_
        filenames = os.listdir( data_dir_ )  # list of filenames in that dir
        print 'filenames', filenames
    
        for filename in filenames:
            print 'reading %s' % filename 
            if filename.endswith( '.txt' ):
                print os.path.join( data_dir_, filename )
                infile = open( os.path.join( data_dir_, filename ), 'r' )
                data = infile.read()
                #print data
                data = fcomma2fdot( data )    
                infile.close()
                outfile = open( os.path.join( data_dir_, filename ), 'w' )
                outfile.write( data )
                outfile.close()

        data_dir_ = os.path.join(data_dir, 'geo_data')
        print 'data_dir_', data_dir_
        filenames = os.listdir( data_dir_ )  # list of filenames in that dir
        print 'filenames', filenames
    
        for filename in filenames:
            print 'reading %s' % filename 
            if filename.endswith( '.txt' ):
                print os.path.join( data_dir_, filename )
                infile = open( os.path.join( data_dir_, filename ), 'r' )
                data = infile.read()
                #print data
                data = fcomma2fdot( data )    
                infile.close()
                outfile = open( os.path.join( data_dir_, filename ), 'w' )
                outfile.write( data )
                outfile.close()
