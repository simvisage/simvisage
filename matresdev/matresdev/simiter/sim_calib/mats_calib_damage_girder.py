
import os.path

from mathkit.array.smoothing import smooth
from mathkit.mfn import MFnLineArray

from matresdev.db.simdb.simdb import simdb
from mats_calib_damage_fn import MATSCalibDamageFn
import numpy as np
import pylab as p


test_file = os.path.join(simdb.exdata_dir,
                         'tensile_tests',
                         'buttstrap_clamping',
                         '2017-06-22-TTb-sig-eps-dresden-girder',
                         'tt-sig-eps-800tex-test.txt')


class MATSCalibDamageFnSigEps(MATSCalibDamageFn):

    def _get_mfn_line_array_target(self):
        data = np.loadtxt(test_file)
        xdata, ydata = data.T
        xdata *= 0.001
        eps = xdata  # smooth(xdata, 8, 'flat')
        sig = ydata  # smooth(ydata, 8, 'flat')

        return MFnLineArray(xdata=eps, ydata=sig)

    def _get_test_key(self):
        return 'girder_dresden'


cf = MATSCalibDamageFnSigEps(xtol=1e-3)
ax = p.subplot(121)

data = np.loadtxt(test_file)
xdata, ydata = data.T
xdata *= 0.001

E_c = np.average((ydata[1:3] - ydata[0]) / (xdata[1:3] - xdata[0]))
print 'e_mod', E_c

p.plot([0, 0.001], [0, E_c * 0.001], color='blue')

p.plot(xdata, ydata)
cf.mfn_line_array_target.mpl_plot(ax)


def run():
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

    fitter = MATSCalibDamageFnSigEps(KMAX=300,
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
    p.show()

    return

    #---------------------------
    # basic testing of fitter methods:
    #---------------------------

    # set to True for basic testing of the methods:
    basic_tests = False

    if basic_tests:
        fitter.run_through()
        #    fitter.tloop.rtrace_mngr.rtrace_bound_list[0].configure_traits()
        fitter.tloop.rtrace_mngr.rtrace_bound_list[0].redraw()
        last_strain_run_through = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.xdata[:]
        last_stress_run_through = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.ydata[:]
        print 'last strain (run-through) value', last_strain_run_through
        print 'last stress (run-through) value', last_stress_run_through

        fitter.tloop.reset()
        fitter.run_step_by_step()
        # fitter.tloop.rtrace_mngr.rtrace_bound_list[0].configure_traits()
        fitter.tloop.rtrace_mngr.rtrace_bound_list[0].redraw()
        last_strain_step_by_step = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.xdata[:]
        last_stress_step_by_step = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.ydata[:]
        print 'last stress (step-by-step) value', last_stress_step_by_step

        fitter.run_trial_step()
        fitter.run_trial_step()
        fitter.tloop.rtrace_mngr.rtrace_bound_list[0].redraw()
        strain_after_trial_steps = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.xdata[:]
        stress_after_trial_steps = fitter.tloop.rtrace_mngr.rtrace_bound_list[0].trace.ydata[:]
        print 'stress after trial', stress_after_trial_steps

        fitter.init()
        # fitter.mats2D_eval.configure_traits()
        lof = fitter.get_lack_of_fit(1.0)
        print '1', lof
        lof = fitter.get_lack_of_fit(0.9)
        print '2', lof

        # fitter.tloop.rtrace_mngr.configure_traits()
        fitter.run_trial_step()

    else:
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        ibvpy_app = IBVPyApp(ibv_resource=fitter)
        ibvpy_app.main()


if __name__ == '__main__':
    run()
