=================================
Material research and development
=================================

This package summarizes procedures included in the process
of material development accompanied with the development 
of numerical models. The following tasks are supported

Data management
===============

The included packages provide support for class extensions
associated with the traited python class.

matresdev.db
------------
Provides the basic functionality for the management of 
class instances within its database extension.

matresdev.db.exdb
-----------------
Persistent data storage for experimental results within a directory
tree. The subdirectories can be defined as class extensions
of python classes. Measured data can be inserted in form
of files produced by the measuring equipment. They are
automatically incorporated into the class extension defined
in the subdirectory. A class extension can be distributed
over several directories. This means that one can keep both
the hierarchical and chronological order of the experiments
within the file tree on the one hand, and one has the 
overall view over all instances of a certain experimental setup 
stored in the directory tree of the database store.

matresdev.db.matdb
------------------
Storage for *objective* material characteristics.
These characteristics are valid for a general material point.
