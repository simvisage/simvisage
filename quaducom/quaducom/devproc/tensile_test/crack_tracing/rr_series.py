'''
Created on Jan 28, 2015

'''

from matresdev.db.simdb import SimDB
simdb = SimDB()
from matresdev.db.exdb import ExRun
from aramis_cdt import AramisUI

import os
import numpy as np
import pylab as p
params = {'legend.fontsize': 10,
          # 'legend.linewidth': 2
          }
p.rcParams.update(params)

test_files = ['TTb-6c-2cm-0-TU-V%d.DAT' % i for i in [3]]
test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2013-12-02_TTb-6c-2cm-0-TU_Aramis2d_RR'
                              )
e_list_6c = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

test_files = ['TTb-4c-2cm-0-TU-V%d.DAT' % i for i in [3]]
test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2013-12-01_TTb-4c-2cm-0-TU_Aramis2d_RR'
                              )
e_list_4c = [ExRun(data_file=os.path.join(test_file_path, test_file))
             for test_file in test_files]

res_key = 'Xf15s1-Yf15s4'

# (integ_radius, dux_avg, ddd_ux_avg)
param_dict = {'TTb-4c-2cm-0-TU-V1.DAT': (22, -4e-4, 0.003),
              'TTb-4c-2cm-0-TU-V2.DAT': (20, -4e-4, 0.0025),
              'TTb-4c-2cm-0-TU-V3.DAT': (20, -4e-4, 0.0035),
              'TTb-4c-2cm-0-TU-V4.DAT': (22, -0.00028, 0.003),
              'TTb-4c-2cm-0-TU-V5.DAT': (20, -4e-4, 0.0035),
              'TTb-6c-2cm-0-TU-V1.DAT': (20, -2e-4, 0.002),
              'TTb-6c-2cm-0-TU-V2.DAT': (22, -2.5e-4, 0.0015),
              'TTb-6c-2cm-0-TU-V3.DAT': (20, -2e-4, 0.0025),
              'TTb-6c-2cm-0-TU-V4.DAT': (20, -2.5e-4, 0.0015),
              'TTb-6c-2cm-0-TU-V5.DAT': (20, -2.5e-4, 0.0015)
              }

e_list = e_list_6c


def plot_all():

    fig = p.figure()
    fig.subplots_adjust(
        left=0.03, right=0.97, bottom=0.04, top=0.96, wspace=0.25, hspace=0.2)

    for idx, e_run in enumerate(e_list):
        e = e_run.ex_type
        e.aramis_resolution_key = res_key

        print
        # print 'crack filter', e.crack_filter_avg
        params = param_dict.get(os.path.basename(e.data_file))
        e.aramis_cdt.aramis_data.integ_radius = params[0]
        e.aramis_cdt.ddd_ux_avg_threshold = params[1]
        e.aramis_cdt.ddd_ux_threshold = params[1]
        e.aramis_cdt.d_ux_avg_threshold = params[2]
        # e.aramis_cdt.crack_detection_step = e.maxF_idx_aramis
        # e.number_of_cracks_aramis
        # e.process_aramis_cdt_data()

        dux_t_cr = e.dux_t_cr
        eps_t_cr = e.eps_t_cr
<<<<<<< HEAD

        # evaluate the derivative of strain field
        twin = 5
        deps_t_cr = eps_t_cr[twin:, :] - eps_t_cr[:-twin,:]
=======
        
        # evaluate the derivative of strain field
        twin = 5
        deps_t_cr = eps_t_cr[twin:,:] - eps_t_cr[:-twin,:]
>>>>>>> branch 'master' of https://rosoba@github.com/simvisage/simvisage.git
        argmax_deps_t_cr = np.argmax(deps_t_cr, axis=0)
        cr_enum = np.arange(len(argmax_deps_t_cr))

        max_deps_t_cr = deps_t_cr[argmax_deps_t_cr, cr_enum]
        max_eps_t_cr = eps_t_cr[argmax_deps_t_cr + twin, cr_enum]
        max_t_cr = e.t_aramis_asc[argmax_deps_t_cr]

        axes = p.subplot(231)
        axes.plot(e.t_aramis_asc, dux_t_cr, color='darkblue', label='1')
        p.ylim((-0.05, 0.30))

        axes = p.subplot(232)
        axes.plot(e.t_aramis_asc, eps_t_cr, label='1')

        axes.plot(e.t_aramis_asc, e.min_eps_t_cr, color='red',
                  linewidth=3, label='1')

        axes.plot(max_t_cr, max_eps_t_cr, 'ro')
        
        p.ylim((-0.002, 0.02))

<<<<<<< HEAD
        p.ylim((-0.002, 0.02))
=======
>>>>>>> branch 'master' of https://rosoba@github.com/simvisage/simvisage.git

        axes = p.subplot(233)

        axes.plot(max_t_cr, max_deps_t_cr, 'ro')

        axes.plot(e.t_aramis_asc[:deps_t_cr.shape[0]], deps_t_cr,
                  linewidth=1, label='1')
        p.ylim((-0.0005, 0.005))

<<<<<<< HEAD
=======
        

>>>>>>> branch 'master' of https://rosoba@github.com/simvisage/simvisage.git
        axes = p.subplot(234)
        e._plot_sigc_eps_ironed(axes, color='darkblue', label='eps')
        
        axes = p.subplot(235)
        axes.plot(e.F_t_aramis_asc, eps_t_cr)
        p.ylim((-0.002, 0.02))
        
        axes = p.subplot(236)
        axes.plot(e.t_aramis_asc, e.F_t_aramis_asc)


        #print the crack positions
        m = e.aramis_cdt.crack_detect_mask_avg.copy()
        cr_idx = np.where(m)[0]
        x = e.aramis_field_data.x_arr_0[0,:]
        print 'position', [x[cr_idx]]
               
        # print the crack initiating force        
        print 'force', [e.F_t_aramis_asc[argmax_deps_t_cr]]

        axes = p.subplot(235)
        axes.plot(e.F_t_aramis_asc, eps_t_cr)
        p.ylim((-0.002, 0.02))

        axes = p.subplot(236)
        axes.plot(e.t_aramis_asc, e.F_t_aramis_asc)

        # print the crack positions
        m = e.aramis_cdt.crack_detect_mask_avg.copy()
        cr_idx = np.where(m)[0]
        x = e.aramis_field_data.x_arr_0[0, :]
        print 'position', [x[cr_idx]]

        # print the crack initiating force
        print 'force', [e.F_t_aramis_asc[argmax_deps_t_cr]]

        p.show()
        
        aui = AramisUI(aramis_info=e.aramis_info,
               aramis_data=e.aramis_field_data,
               aramis_cdt=e.aramis_cdt)
        aui.configure_traits()

if __name__ == '__main__':
    plot_all()
