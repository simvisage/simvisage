'''
Created on Jan 15, 2013

@author: rostar
'''
#---#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Jun 14, 2010 by: rch

from quaducom.micro.resp_func.CB_clamped_rand_xi import CBClampedRandXi
import numpy as np
from matplotlib import pyplot as plt
from spirrid import SPIRRID
from spirrid.rv import RV

if __name__ == '__main__':
    def CB_composite_stress(w, tau, E_f, V_f, r, m, sV0, Pf, n_int):
        cb = CBClampedRandXi()
        spirrid = SPIRRID(q=cb,
                    sampling_type='PGrid',
                    eps_vars=dict(w=w),
                    theta_vars=dict(tau=tau, E_f=E_f, V_f=V_f, r=r,
                               m=m, sV0=sV0, Pf=Pf),
                    n_int=n_int)
        if isinstance(r, RV):
            r_arr = np.linspace(r.ppf(0.001), r.ppf(0.999), 300)
            Er = np.trapz(r_arr ** 2 * r.pdf(r_arr), r_arr)
        else:
            Er = r ** 2
        sigma_c = spirrid.mu_q_arr / Er    
        plt.plot(w, sigma_c, lw=2, label='sigma c')
        #plt.ylim(0, 35)
        plt.xlabel('w [mm]')
        plt.ylabel('sigma_c [MPa]')
        plt.legend(loc='best')
        plt.show()
    w = np.linspace(0, .5, 300)
    tau = RV('weibull_min', shape=3., scale=.28)
    E_f = 170e3
    V_f = 0.01
    r = 0.003#RV('uniform', loc=0.001, scale=0.004)
    m = 5.
    # sV0=XXX corresponds to sL0=0.02 at L0=100 and r=0.002
    sV0 = 3e-3
    Pf = RV('uniform', loc=0., scale=1.0)
    n_int = 30
    #cb = CBResidual()
    #plt.plot(w, cb(w, 0.5, E_f, V_f, 0.001, m, sV0, 0.5) / 0.001 ** 2, label='r=0.001')
    #plt.plot(w, cb(w, 0.5, E_f, V_f, 0.002, m, sV0, 0.5) / 0.002 ** 2, label='r=0.002')
    #plt.plot(w, cb(w, 0.5, E_f, V_f, 0.003, m, sV0, 0.5) / 0.003 ** 2, label='r=0.003')
    #plt.legend()
    #plt.show()
    CB_composite_stress(w, tau, E_f, V_f, r, m, sV0, Pf, n_int)
