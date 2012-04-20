'''
Created on Mar 16, 2012

The script demonstrates the functionality of the simple
object-oriented persistent storage with cross references.

Two classes with extensions for storing the data 
are defined first, the database is populated with
random values. Afterwards the evolution of the class
diagram is demonstrated.

@author: rch
'''

from matresdev.db import SimDBClass
from matresdev.db import SimDBClassExt
from enthought.traits.api import Int, Property, Str, Float

if __name__ == '__main__':

    #===========================================================================
    # Simple - persistent class
    #===========================================================================
    class MyDB(SimDBClass):
        v1 = Int(10, simdb = True)
        v2 = Int(20, simdb = True)

    MyDB.db = SimDBClassExt(klass = MyDB, verbose = 'io')
    MyDB.db.delete_instances()
    
    # populate the extension with program constants (non-editable)
    MyDB.db.constants = {'DB_1' : MyDB(v1 = 10, v2 = 40),
                         'DB_2' : MyDB(v1 = 10, v2 = 50),
                         'DB_3' : MyDB(v1 = 10, v2 = 60) }

    print 'MyDB.db.instances', MyDB.db.instances

    #===========================================================================
    # Class with a cross reference
    #===========================================================================
    class MyDB2(SimDBClass):
        my_key = Str(simdb = True)
        my_ref = Property
        def _get_my_ref(self):
            return MyDB.db[ self.my_key ]

    MyDB2.db = SimDBClassExt(klass = MyDB2, verbose = 'io')
    
    print 'MyDB2.db.delete_instances()'
    
    MyDB2.db.delete_instances()
    
    print 'MyDB2.db.instances', MyDB2.db.instances
    
    # Insert some instances into MyDB2 extension
    # assignment to the .db extension stores the
    # object within as a pickle file.
    #
    MyDB2.db['DB2_1'] = MyDB2(my_key = 'DB_2')
    
    # Another way how to store an object is to 
    # use the *save* method of the SimDB base class
    # Note that if DB2_2 already exists in the database
    # it will be rewritten here.
    DB_2 = MyDB2(key = 'DB2_2', my_key = 'DB_1')
    DB_2.insert()

    print MyDB2.db['DB2_2'].my_ref.key
    print MyDB2.db.instances

    # now the two instances are stored
    # add a new trait to all instances of MyDB2 with the
    
    MyDB2.add_class_trait('humidity', Float)
    
    for inst in MyDB2.db.inst_list:
        #inst.add_trait('humidity', Float())
        inst.save() 

    print 'humidity'
    print MyDB2.db['DB2_2'].humidity

#    print 'v2'
#    print MyDB2.db['MDBO1'].my_ref.v2
#    
#    print MyDB2.db['MDBO1'].my_ref.key
#    
#    MyDB2.db['MDBO1'].my_key = 'object_2'
#    
#    print MyDB2.db['MDBO1'].my_ref.key
#    
#    MyDB2.db['MDBO1'].save()
#    
    MyDB2.db.configure_traits()
