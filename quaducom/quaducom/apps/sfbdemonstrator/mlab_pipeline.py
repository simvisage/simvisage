'''
Created on Mar 14, 2014

@author: rch
'''

import mayavi.tools as tools
import mayavi.mlab as mlab
import numpy as np
from matresdev.db import SimDB
from os.path import join

simdb = SimDB()

simdata_dir = join(simdb.simdata_dir, 'MRquarterDB')
f = mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(1200, 600))

for i in range(0, 15):

    in_fname = join(simdata_dir, 'max_omega_i0nodes_%d.vtk' % i)
    out_fname_top = '/tmp/kde-rch/demonstrator_damage_top_%02d.png' % i
    out_fname_bot = '/tmp/kde-rch/demonstrator_damage_bot_%02d.png' % i

    print 'fname', in_fname, out_fname_top
    pipelines = [tools.pipeline.open(in_fname) for i in range(4)]
    points_arr = [pipeline.outputs[0].points.to_array() for pipeline in pipelines ]

    points_arr[1][:, 0] *= -1.0
    points_arr[2][:, 1] *= -1.0
    points_arr[3][:, 0] *= -1.0
    points_arr[3][:, 1] *= -1.0

    for pipeline in pipelines:
        surf = mlab.pipeline.surface(pipeline, colormap='OrRd')
        lut_manager = surf.module_manager.scalar_lut_manager
        lut_manager.data_range = np.array([0, 1.0], dtype='f')

#    mlab.view(35.0, 65.0, figure=f)
    mlab.view(35.0, 65.0, distance=10.0, focalpoint=np.array([ 0.5, 0.0, 0.0]), figure=f)
    mlab.savefig(out_fname_top)

    mlab.view(35.0, 120.0, distance=10.0, focalpoint=np.array([ 0.5, 0.0, 1.85]), figure=f)
    mlab.savefig(out_fname_bot)

mlab.show()

