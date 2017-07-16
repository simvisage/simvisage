
import os.path

from mathkit.array.smoothing import smooth
from mathkit.mfn import MFnLineArray
from traits.api import File
from matresdev.db.simdb.simdb import simdb
from mats_calib_damage_fn import MATSCalibDamageFn
import numpy as np
import pylab as p


class MATSCalibDamageFnSigEps(MATSCalibDamageFn):

    test_file = File

    def _get_mfn_line_array_target(self):
        data = np.loadtxt(self.test_file)
        xdata, ydata = data.T
        xdata *= 0.001

        xdata = xdata[:np.argmax(ydata)]
        ydata = ydata[:np.argmax(ydata)]

        eps = smooth(xdata, 60, 'flat')
        sig = smooth(ydata, 60, 'flat')

        return MFnLineArray(xdata=eps, ydata=sig)

    def get_E_c(self):

        data = np.loadtxt(self.test_file)
        xdata, ydata = data.T
        xdata *= 0.001
        xdata = xdata[:np.argmax(ydata)]
        ydata = ydata[:np.argmax(ydata)]
        del_idx = np.arange(60)
        xdata = np.delete(xdata, del_idx)
        ydata = np.delete(ydata, del_idx)
        xdata = np.append(0.0, xdata)
        ydata = np.append(0.0, ydata)
        E_c = np.average((ydata[2:4] - ydata[0]) / (xdata[2:4] - xdata[0]))
        return E_c

    def _get_test_key(self):
        return 'girder_dresden'


file_names_800 = ['tt-dk1-800tex.txt',
                  'tt-dk2-800tex.txt',
                  'tt-dk3-800tex.txt',
                  'tt-dk4-800tex.txt']


file_names_3300 = ['tt-dk1-3300tex.txt',
                   'tt-dk2-3300tex.txt',
                   'tt-dk3-3300tex.txt',
                   'tt-dk4-3300tex.txt']


def get_test_files(test_file_names):
    return [os.path.join(simdb.exdata_dir,
                         'tensile_tests',
                         'buttstrap_clamping',
                         '2017-06-22-TTb-sig-eps-dresden-girder',
                         file_name)
            for file_name in test_file_names
            ]


def show_input(ax, test_file):
    cf = MATSCalibDamageFnSigEps(test_file=test_file, xtol=1e-3)

    data = np.loadtxt(test_file)
    xdata, ydata = data.T
    xdata *= 0.001
    xdata = xdata[:np.argmax(ydata)]
    ydata = ydata[:np.argmax(ydata)]
    del_idx = np.arange(60)
    xdata = np.delete(xdata, del_idx)
    ydata = np.delete(ydata, del_idx)
    xdata = np.append(0.0, xdata)
    ydata = np.append(0.0, ydata)

    eps = smooth(xdata, 60, 'flat')
    sig = smooth(ydata, 60, 'flat')

    E_c = np.average((ydata[2:4] - ydata[0]) / (xdata[2:4] - xdata[0]))
    print 'e_mod', E_c

    p.plot([0, 0.001], [0, E_c * 0.001], color='blue')

    p.plot(xdata, ydata)
    cf.mfn_line_array_target.mpl_plot(ax)


def show_test_files(ax, file_names):
    for test_file in get_test_files(file_names):
        show_input(ax, test_file)


def run(test_file):
    #-------------------------------------------------------------------------
    # Example using the mats2d_explore
    #-------------------------------------------------------------------------
    from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore

    from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import \
        MATS2DMicroplaneDamage
    from ibvpy.mats.matsXD.matsXD_cmdm import \
        PhiFnGeneral

    from ibvpy.api import RTraceGraph

    from ibvpy.mats.mats2D.mats2D_explorer_bcond import BCDofProportional
    from os.path import join

    ec = {
        # overload the default configuration
        'bcond_list': [BCDofProportional(max_strain=1.0, alpha_rad=0.0)],
        'rtrace_list': [
            RTraceGraph(name='stress - strain',
                        var_x='eps_app', idx_x=0,
                        var_y='sig_app', idx_y=0,
                        record_on='iteration'),
        ],
    }

    mats_eval = MATS2DMicroplaneDamage(
        n_mp=30,
        # mats_eval = MATS1DMicroplaneDamage(
        elastic_debug=False,
        stress_state='plane_stress',
        symmetrization='sum-type',
        model_version='compliance',
        phi_fn=PhiFnGeneral,
    )

#    print 'normals', mats_eval._MPN
#    print 'weights', mats_eval._MPW

    fitter = MATSCalibDamageFnSigEps(test_file=test_file,
                                     KMAX=300,
                                     tolerance=5e-4,  # 0.01,
                                     RESETMAX=0,
                                     dim=MATS2DExplore(
                                         mats_eval=mats_eval,
                                         explorer_config=ec,
                                     ),
                                     store_fitted_phi_fn=True,
                                     log=False
                                     )

    #------------------------------------------------------------------
    # specify the parameters used within the calibration
    #------------------------------------------------------------------
    #

    # smallest value for matrix E-modulus obtained from cylinder tests (d=150mm)
#        E_m = 18709.5

    # set 'nu'
    # @todo: check values stored in 'mat_db'
    #
    nu = 0.20

    n_steps = 200
    fitter.n_steps = n_steps

    fitter.format_ticks = True
    E_c = fitter.get_E_c()
    fitter.dim.mats_eval.E = E_c
    fitter.dim.mats_eval.nu = nu

    print 'n_steps = %g used for calibration' % n_steps

    max_eps = fitter.max_eps
    print 'max_eps = %g used for calibration' % max_eps

    #------------------------------------------------------------------
    # set 'param_key' of 'fitter' to store calibration params in the name
    #------------------------------------------------------------------
    #
    age = 28
#        param_key = '_age%g_Em%g_nu%g_nsteps%g' % (age, E_m, nu, n_steps)
#        param_key = '_age%g_Ec%g_nu%g_nsteps%g__smoothed' % (age, E_c, nu, n_steps, max_eps)
    param_key = '_age%g_Ec%g_nu%g_nsteps%g_smoothed' % (
        age, E_c, nu, n_steps)

    fitter.param_key = param_key
    print 'param_key = %s used in calibration name' % param_key

    #------------------------------------------------------------------
    # run fitting procedure
    #------------------------------------------------------------------
    #
    fitter.init()
    fitter.fit_response()
    fitter.fitted_phi_fn
    ax = p.subplot(122)
    print 'plotting'
    fitter.fitted_phi_fn.mpl_plot(ax)
    # p.show()

    return fitter.fitted_phi_fn.xdata, fitter.fitted_phi_fn.ydata, E_c


def calibrate(file_names):
    for test_file in get_test_files(file_names):
        xdata, ydata, E_c = run(test_file)
        results = np.c_[xdata, ydata]
        np.savetxt(test_file + 'phi_data', results)
        with open(test_file + 'E_c', 'w') as f:
            f.write(r'''E_c = %g''' % E_c)
    p.show()


def show_results(ax, file_names):
    for test_file in get_test_files(file_names):
        fname = os.path.basename(test_file)
        results = np.loadtxt(test_file + 'phi_data')
        xdata, ydata = results.T
        with open(test_file + 'E_c', 'r') as f:
            E_c_str = f.read()
#         xdata = smooth(xdata, 6, 'flat')
#         ydata = smooth(ydata, 6, 'flat')
        ax.plot(xdata, ydata, label='%s: %s' % (fname, E_c_str))
        ax.set_ylim(ymin=0.0)
    p.legend(loc=1)


if __name__ == '__main__':
    ax = p.subplot(221)
    show_test_files(ax, file_names_800)
    ax = p.subplot(222)
    show_results(ax, file_names_800)
    ax = p.subplot(223)
    show_test_files(ax, file_names_3300)
    ax = p.subplot(224)
    show_results(ax, file_names_3300)
    p.show()
