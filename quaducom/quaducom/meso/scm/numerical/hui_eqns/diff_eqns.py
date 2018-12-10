'''
Created on 17 Aug 2013

@author: Q
'''

from etsproxy.traits.api import HasTraits, Float, Property, cached_property, \
                                 Instance
import numpy as np
from scipy.integrate import cumtrapz
from scipy.integrate import ode
from quaducom.meso.scm.analytical.hui_pho_ibn_smi_model import SFC_Hui

class SCM_Fragmentation(HasTraits):

    hui_model = Property(Instance(SFC_Hui), depends_on='+params')
    @cached_property
    def _get_hui_model(self):
        return SFC_Hui(l0=self.l0, d=self.d, tau=self.tau,
                       sigma0=self.sigma0, rho=self.rho)
    
    l0 = Float(params=True)
    d = Float(params=True)
    tau = Float(params=True)
    sigma0 = Float(params=True)
    rho = Float(params=True)

    def dp_const_tau(self, si, p_i):
        x = self.x
        integ = np.trapz(p_i * (x >= si / 2.), x)
        integs = integ - np.hstack((0.0, cumtrapz(p_i * (x >= si / 2.), x)))
        integs = np.hstack((integs[x > si/2.], np.repeat(integs[-1], np.sum(x < si/2.))))
        dp2 = 2 *self.hui_model.h(si) * integs * (x >= si / 2.) * (x < si)
        dp3 = (2* self.hui_model.h(si) * integs - (x - si) * p_i * self.hui_model.h(si)) * (x >= si)  
        return dp2 + dp3

    def h_TRC(self, s, x):
        return self.rho * s ** (self.rho - 1.) * (x > s / 2.)

    def dp_TRC(self, si, p_i):
        x = self.x
        integ = np.trapz(p_i * (x >= si / 2.), x)
        integs = integ - np.hstack((0.0, cumtrapz(p_i * (x >= si / 2.), x)))
        integs = np.hstack((integs[x > si/2.], np.repeat(integs[-1], np.sum(x < si/2.))))
        dp = 2* self.h_TRC(si, x) * integs - (x - si) * p_i * self.h_TRC(si, x - si/2.)
        return dp

    def p_num_euler(self, s):
        p_i = self.hui_model.p_x(1., self.x)
        s_arr = np.linspace(1., s, 500)
        ds = s_arr[1] - s_arr[0]
        for si in s_arr:
            p_i += self.dp_const_tau(si, p_i) * ds
        return p_i
    
    def p_num_scipy(self, s):
        p_0 = self.hui_model.p_x(1., self.x)
        r = ode(self.dp_TRC).set_integrator('dopri5')
        r.set_initial_value(p_0, 1.0)
        t1 = s
        dt = 0.1
        while r.successful() and r.t < t1:
            r.integrate(r.t+dt)
        return r.y
        

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    import time
    sfc = SCM_Fragmentation(l0=1., d=0.007, tau=0.1, sigma0=2200., rho=5.0)
    
    for rho in np.array([10.]):
        sfc.rho = rho
        sfc.x = np.linspace(0.01, 7., 300)
        #print np.trapz(x*pdf_x_num, x)
        #cdf_x = np.hstack((0., cumtrapz(pdf_x, x)))
        plt.plot(sfc.x, sfc.hui_model.p_x(1., sfc.x), label='analytical 1.0')
        plt.plot(sfc.x, sfc.hui_model.p_x(1.5, sfc.x), label='analytical 1.5')
        t = time.clock()
        euler = sfc.p_num_euler(1.5)
        print('euler', time.clock() - t, 's')
        t = time.clock()
        rk = sfc.p_num_scipy(1.5)
        print('RK', time.clock() - t, 's')
        plt.plot(sfc.x, euler, label='euler')
        plt.plot(sfc.x, rk, label='ru-ku')
    plt.ylim(0,4.0)
    plt.legend()
    plt.show()
    