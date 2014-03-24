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

    def fcomma2fdot(x):
        return string.replace(x, ',', '.')

#    data_dir = os.path.join(simdb.simdb_dir,
#                'simdata',
#                'input_data_barrelshell',
#                '2cm-feines-Netz-EC1')

    data_dir = os.path.join(simdb.simdb_dir,
                            'simdata',
                            'input_data_slabtest',
                            'slabtest_125x125x6cm_quarter_lineload'
                            )

    data_dir_ = os.path.join(data_dir, 'state_data', 'Auflagerreaktionen')
    print 'data_dir_', data_dir_
    filenames = os.listdir(data_dir_)  # list of filenames in that dir
    print 'filenames', filenames

    for filename in filenames:
        print 'reading %s' % filename
        if filename.endswith('.txt'):
            print os.path.join(data_dir_, filename)
            infile = open(os.path.join(data_dir_, filename), 'r')
            data = infile.read()
            # print data
            data = fcomma2fdot(data)
            infile.close()
            outfile = open(os.path.join(data_dir_, filename), 'w')
            outfile.write(data)
            outfile.close()

    data_dir_ = os.path.join(data_dir, 'state_data', 'Flaechenschnittgroessen', 'Schwerpunkt')
    print 'data_dir_', data_dir_
    filenames = os.listdir(data_dir_)  # list of filenames in that dir
    print 'filenames', filenames

    for filename in filenames:
        print 'reading %s' % filename
        if filename.endswith('.txt'):
            print os.path.join(data_dir_, filename)
            infile = open(os.path.join(data_dir_, filename), 'r')
            data = infile.read()
            # print data
            data = fcomma2fdot(data)
            infile.close()
            outfile = open(os.path.join(data_dir_, filename), 'w')
            outfile.write(data)
            outfile.close()

    data_dir_ = os.path.join(data_dir, 'state_data', 'Knotendeformationen')
    print 'data_dir_', data_dir_
    filenames = os.listdir(data_dir_)  # list of filenames in that dir
    print 'filenames', filenames

    for filename in filenames:
        print 'reading %s' % filename
        if filename.endswith('.txt'):
            print os.path.join(data_dir_, filename)
            infile = open(os.path.join(data_dir_, filename), 'r')
            data = infile.read()
            # print data
            data = fcomma2fdot(data)
            infile.close()
            outfile = open(os.path.join(data_dir_, filename), 'w')
            outfile.write(data)
            outfile.close()

    data_dir_ = os.path.join(data_dir, 'geo_data')
    print 'data_dir_', data_dir_
    filenames = os.listdir(data_dir_)  # list of filenames in that dir
    print 'filenames', filenames

    for filename in filenames:
        print 'reading %s' % filename
        if filename.endswith('.txt'):
            print os.path.join(data_dir_, filename)
            infile = open(os.path.join(data_dir_, filename), 'r')
            data = infile.read()
            # print data
            data = fcomma2fdot(data)
            infile.close()
            outfile = open(os.path.join(data_dir_, filename), 'w')
            outfile.write(data)
            outfile.close()
