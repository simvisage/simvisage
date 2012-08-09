'''
Created on Sep 30, 2011

The example implements micromechanical responses of a crack opening
in a composite material where a single discrete crack is bridged
by:
a) single homogeneous filament
b) yarn with random filament properties
The diagram shows a comparison between the formulation with discrete
filament strength and length dependent strength with residual filament
stress due to pullout of the broken filament from the matrix.
BC - fixed displacement at both ends (Ll and Lr) of the filament.
@author: rostar
'''

from quaducom.micro.resp_func.cb_emtrx_clamped_fiber_stress import \
    CBEMClampedFiberStress
from quaducom.micro.resp_func.cb_emtrx_clamped_fiber_stress_residual import \
    CBEMClampedFiberStressResidual

from stats.spirrid.spirrid import SPIRRID
from stats.spirrid.rv import RV
from matplotlib import pyplot as plt
import numpy as np

if __name__ == '__main__':

    # filaments
    r = 0.00345
    V_f = 0.0103
    tau = RV('uniform', loc=0.1, scale=.3)
    E_f = 200e3
    E_m = 25e3
    l = RV('uniform', scale=10., loc=2.)
    theta = 0.0
    xi = RV('weibull_min', scale=0.02, shape=5)
    phi = 1.
    Ll = 20.
    Lr = 30.
    s0 = 0.0205
    m = 5.0
    Pf = RV('uniform', loc=0., scale=1.0)

    w = np.linspace(0, 1.3, 100)

    cb_emtrx = CBEMClampedFiberStress()
    s = SPIRRID(q=cb_emtrx,
         sampling_type='PGrid',
         evars=dict(w=w),
         tvars=dict(tau=tau, l=l, E_f=E_f, theta=theta, xi=xi, phi=phi,
                    E_m=E_m, r=r, V_f=V_f, Ll=Ll, Lr=Lr),
         n_int=30)

    cb_emtrx_r = CBEMClampedFiberStressResidual()
    s_r = SPIRRID(q=cb_emtrx_r,
         sampling_type='PGrid',
         evars=dict(w=w),
         tvars=dict(tau=tau, l=l, E_f=E_f, theta=theta, Pf=Pf, phi=phi,
                    E_m=E_m, r=r, V_f=V_f, Ll=Ll, Lr=Lr, s0=s0, m=m),
         n_int=30)

    mean_params = []
    for param in s.tvar_lst:
        if isinstance(param, RV):
            mean_params.append(param._distr.mean)
        else:
            mean_params.append(param)

    mean_params_r = []
    for param in s_r.tvar_lst:
        if isinstance(param, RV):
            mean_params_r.append(param._distr.mean)
        else:
            mean_params_r.append(param)
    fil = cb_emtrx(w, *mean_params)
    fil_r = cb_emtrx_r(w, *mean_params_r)

    def w_sigmaf():
        plt.plot(w, s.mu_q_arr, lw=2, color='blue',
                 ls='dashed', label='no res stress')
        plt.plot(w, fil, color='blue', lw=2,
                 ls='solid', label='filament')
        plt.plot(w, s_r.mu_q_arr, lw=2, color='red',
                 ls='dashed', label='residual stress')
        plt.plot(w, fil_r, color='red', lw=2,
                 ls='solid', label='filament')
        plt.title('comparison no res stress vs. res stress after breakage')
        plt.legend(loc='best')
        plt.show()

    w_sigmaf()
