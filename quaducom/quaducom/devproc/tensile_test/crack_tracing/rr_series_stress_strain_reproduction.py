'''
Trying to reproduce the stress-strain curves from RR Thesis
Created on Jan 28, 2015

'''

from exp_att_db import ExpATTDB, ExpTTDB
from matresdev.db.simdb import SimDB
from matresdev.db.exdb import ExRun
simdb = SimDB()
import os
import pylab as p
params = {'legend.fontsize': 10,
          #'legend.linewidth': 2
          }
p.rcParams.update(params)

test_files = ['TTb-4c-2cm-0-TU-V%d.DAT' % i for i in range(1, 6)]
test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2013-12-01_TTb-4c-2cm-0-TU_Aramis2d_RR'
                              )
e_list_4c = [ExRun(data_file=os.path.join(test_file_path, test_file),
                   delta_t_aramis=5)
             for test_file in test_files]

test_files = ['TTb-6c-2cm-0-TU-V%d.DAT' % i for i in range(1, 6)]
test_file_path = os.path.join(simdb.exdata_dir,
                              'tensile_tests', 'buttstrap_clamping',
                              '2013-12-02_TTb-6c-2cm-0-TU_Aramis2d_RR'
                              )
e_list_6c = [ExRun(data_file=os.path.join(test_file_path, test_file),
                   delta_t_aramis=5)
             for test_file in test_files]


def plot_all():

    fig = p.figure()
    fig.subplots_adjust(
        left=0.03, right=0.97, bottom=0.04, top=0.96, wspace=0.25, hspace=0.2)

    axes = p.subplot(221)

    for idx, ex_run in enumerate(e_list_4c):
        e = ex_run.ex_type
        e._plot_sigtex_eps_ironed(axes)

    axes = p.subplot(222)

    for idx, ex_run in enumerate(e_list_6c):
        e = ex_run.ex_type
        e._plot_sigtex_eps_ironed(axes)

    axes = p.subplot(223)

    for idx, ex_run in enumerate(e_list_4c):
        e = ex_run.ex_type
        print 'reinf', e.ccs.rho_c
        axes.plot(e.eps_ironed, e.sig_c_ironed, color='blue')
        axes.legend()

    max_eps = e.eps_ironed[-1]
    axes.plot([0, max_eps], [0, e.ccs.rho_c * e.ccs.E_tex * max_eps],
              linestyle='--')

    for idx, ex_run in enumerate(e_list_6c):
        e = ex_run.ex_type
        axes.plot(e.eps_ironed, e.sig_c_ironed, color='green')
        axes.legend()

    max_eps = e.eps_ironed[-1]
    axes.plot([0, max_eps], [0, e.ccs.rho_c * e.ccs.E_tex * max_eps],
              linestyle='--')


if __name__ == '__main__':
    plot_all()
    p.show()
