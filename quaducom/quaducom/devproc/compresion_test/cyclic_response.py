'''
Created on Nov 3, 2017

This module provides the basic infrastructure
for extracting the cyclic response from a
test

It identifies the cycles, and construct 
cycle-based indexing and filtering of the data

  * define the methods as cached properties
  * define a view to cyclic response
  
Storage issues
-----------------
  * pandas
  * binary data format - npy 

Processing issues
-----------------
  * binary criteria for loop identification
  * integral values - hysteretic loops
  * derivatives - stiffness / damage

Technical issues
----------------
  * data - reader - file selection 
  * database coupling
  
Plotting
--------
  * which views are useful - fatigue creep
  * upper-lower bound
  * selection of loops (regular, logarithmic, customized)
  * smoothing 
  


@author: rch
'''


import os

from traits.api import \
    HasStrictTraits, Array, Int, Property, cached_property

from matresdev.db.simdb import SimDB
import numpy as np


simdb = SimDB()


class CyclicResponse(HasStrictTraits):

    data = Array(np.float_, input=True)

    def get_index_at_max_force(self):
        f = data[:, 1]
        return np.argmax(f)

    argmax_force = Property(depends_on='+input')
    '''Index at the maximum force
    '''
    @cached_property
    def _get_argmax_force(self):
        f = data[:, 1]
        return np.argmin(f)

    u = Property(depends_on='+input')
    '''Displacement
    '''
    @cached_property
    def _get_u(self):
        return -self.data[:self.argmax_force, 2]

    f = Property(depends_on='+input')
    '''Force
    '''
    @cached_property
    def _get_f(self):
        return -self.data[:self.argmax_force, 1]

    t = Property(depends_on='+input')
    '''Force
    '''
    @cached_property
    def _get_t(self):
        return self.data[:self.argmax_force, 0]


if __name__ == '__main__':

    import pylab as p

    # define the path to the file
    test_file_path = os.path.join(simdb.exdata_dir,
                                  'compression_tests', 'cylinder_tests',
                                  )

    test_file = 'WinConFat_CT_120_1_5.csv'
    f_file = os.path.join(test_file_path, test_file)

    # read the data
    data = np.loadtxt(f_file, dtype=np.float_,
                      skiprows=2, delimiter=';')

    cr = CyclicResponse(data=data)

    ax1 = p.subplot(2, 3, 1)
    ax2 = p.subplot(2, 3, 2)
    ax3 = p.subplot(2, 3, 4)
    ax4 = p.subplot(2, 3, 5)
    ax5 = p.subplot(2, 3, 3)
    ax6 = p.subplot(2, 3, 6)

    ax1.plot(cr.t, cr.f)

    delta_arg2 = 20

    df = cr.f[2 * delta_arg2:] - cr.f[:-2 * delta_arg2]

    ddf = df[2 * delta_arg2:] - df[:-2 * delta_arg2]

    ax2.plot(cr.u[delta_arg2:-delta_arg2], df)

    ax3.plot(cr.t[delta_arg2:-delta_arg2], df)

    ax4.plot(cr.t[2 * delta_arg2:-2 * delta_arg2], ddf)

    df_threshold = 0.0
    ddf_threshold = 0.0
    Df = df[delta_arg2:-delta_arg2]
    print Df.shape
    print ddf.shape
    up_args_dd = np.where(((Df[1:] * Df[:-1] < df_threshold) *
                           ((ddf[1:] + ddf[:-1]) / 2.0 < ddf_threshold)))[0]
    up_args_d = up_args_dd + delta_arg2
    up_args = up_args_d + delta_arg2
    ax1.plot(cr.t[up_args], cr.f[up_args], 'go')
    ax2.plot(cr.u[delta_arg2:-delta_arg2][up_args_d], df[up_args_d], 'go')
    ax3.plot(cr.t[delta_arg2:-delta_arg2][up_args_d], df[up_args_d], 'go')
    ax4.plot(cr.t[up_args], ddf[up_args_dd], 'go')

    down_args_dd = np.where(((Df[1:] * Df[:-1] < df_threshold) *
                             ((ddf[1:] + ddf[:-1]) / 2.0 > ddf_threshold)))[0]
    down_args_d = down_args_dd + delta_arg2
    down_args = down_args_d + delta_arg2
    ax1.plot(cr.t[down_args], cr.f[down_args], 'ro')
    ax2.plot(cr.u[delta_arg2:-delta_arg2][down_args_d], df[down_args_d], 'ro')
    ax3.plot(cr.t[delta_arg2:-delta_arg2][down_args_d], df[down_args_d], 'ro')
    ax4.plot(cr.t[down_args], ddf[down_args_dd], 'ro')

    f_envelope = np.hstack([cr.f[:up_args[0]], cr.f[up_args[1:]]])
    u_envelope = np.hstack([cr.u[:up_args[0]], cr.u[up_args[1:]]])
    ax5.plot(u_envelope, f_envelope)
    ax5.plot(cr.u[down_args], cr.f[down_args], color='red')

    selected_cycles = np.array([2, 20, 50], dtype=np.int_)

    selected_cycles_start_args = up_args[selected_cycles]
    selected_cycles_end_args = up_args[selected_cycles + 1]

    for sc in selected_cycles:
        s_arg = up_args[sc]
        e_arg = up_args[sc + 1]
        ax6.plot(cr.u[s_arg:e_arg],
                 cr.f[s_arg:e_arg],
                 )

    p.show()
