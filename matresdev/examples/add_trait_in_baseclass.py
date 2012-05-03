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
from enthought.traits.api import \
    Int, Property, Str, Float, Date

if __name__ == '__main__':

    #===========================================================================
    # Simple - persistent class
    #===========================================================================
    class ExpBase(SimDBClass):
        pass

    class ExpBT3Pt(ExpBase):
        length = Float(0.50, simdb = True)

    ExpBT3Pt.db = SimDBClassExt(klass = ExpBT3Pt, verbose = 'io')
    ExpBT3Pt.db.delete_instances()
    
    # Insert some instances into MyDB2 extension
    # assignment to the .db extension stores the
    # object within as a pickle file.
    #
    ExpBT3Pt.db['DB2_1'] = ExpBT3Pt(length = 0.5)

    # Another way how to store an object is to 
    # use the *save* method of the SimDB base class
    # Note that if DB2_2 already exists in the database
    # it will be rewritten here.
    DB_2 = ExpBT3Pt(key = 'DB2_2', length = 0.8)
    DB_2.insert()

    print ExpBT3Pt.db.instances

    ExpBT3Pt.db.configure_traits()
    
    # now the two instances are stored
    # add a new trait to all instances of MyDB2 with the
    
    ExpBase.add_class_trait('production_date', Date)
    ExpBase.add_class_trait('testing_date', Date)
    
    for num, inst in enumerate(ExpBT3Pt.db.inst_list):
        #inst.add_trait('humidity', Float())
        inst.save() 

    ExpBT3Pt.db.configure_traits()

    #===========================================================================
    # Redefine the class with the added trait
    #===========================================================================

    class ExpBase(SimDBClass):
        production_date = Date(simdb = True)
        testing_date = Date(simdb = True)

    class ExpBT3Pt(ExpBase):
        length = Float(0.50, simdb = True)
        
    ExpBT3Pt.db = SimDBClassExt(klass = ExpBT3Pt, verbose = 'io')
    
    ExpBT3Pt.db.configure_traits()            
