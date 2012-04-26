'''
Created on Apr 20, 2012
stress - strain relation for composites according to Ahn and Curtin 1996
@author: rostar
'''
from material import Material
from enthought.traits.api import Property, cached_property
from scipy.stats import weibull_min
import numpy as np
from math import e
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from scipy.optimize import brentq

class Curtin(Material):
    
    L = Property()
    def _get_S(self):
        '''distribution of the short fragments'''
        return weibull_min(8, scale = 10.,loc = 1.)
    
    M = Property()
    def _get_M(self):
        '''distribution of the medium fragments'''
        return weibull_min(10, scale = 15.,loc = 3.)
    
    S = Property()
    def _get_L(self):
        '''distribution of the long fragments'''
        return weibull_min(12, scale = 20.,loc = 5.)

    alpha = Property()
    def _get_alpha(self):
        return self.Vm*self.Em/self.Vf/self.Ec
    
    delta = Property()
    def _get_delta(self):
        return self.r/2./self.tau * self.alpha*self.sigma
        
    sigma_fL = Property()
    def _get_sigma_fL(self):
        return self.delta * self.alpha * self.sigma * self.L.moment(1) + \
                self.sigma * self.Ef/self.Ec * self.L.moment(2)

    sigma_fM = Property()
    def _get_sigma_fM(self):
        return self.sigma/self.Vf*self.M.moment(1) - \
                self.alpha/4./self.delta*self.sigma*self.M.moment(2)

    sigma_fS = Property()
    def _get_sigma_fS(self):
        return self.sigma/self.Vf*self.S.moment(1) - \
                self.alpha/4./self.delta*self.sigma*self.S.moment(2)

    eps = Property()
    def _get_eps(self):
        return 1./self.Lc/self.Ef*(self.sigma_fL+self.sigma_fM+self.sigma_fS)

# auxiliary function
    def ints(self,t):
        s = np.linspace(1e-10,t,500)
        return np.trapz((1.-e**(-s))/s,s)
       
    def intt(self,psi):
        t = np.linspace(1e-10,psi,500)
        vect_ints = np.vectorize(self.ints)
#        print e**(-2*vect_ints(t))
#        print 'integ value', np.trapz(e**(-2*vect_ints(t)), t)
#        plt.plot(t, e**(-2*vect_ints(t)))
#        plt.show()
        return np.trapz(e**(-2*vect_ints(t)), t)

    psi_line = Property()
    @cached_property
    def _get_psi_line(self):
        eta_psi_vect = np.vectorize(self.intt)
        psi = np.linspace(1e-3,4.99,500)
        eta_psi_func = eta_psi_vect(psi)
        eta2 = np.linspace(0.7476-e**(-2.*0.577216)/5,0.74759,500)
        psi2 = e**(-2.*0.577216)/(0.7476-eta2)
        return MFnLineArray(xdata = np.hstack((eta_psi_func,eta2)),
                            ydata = np.hstack((psi,psi2)))
        #return MFnLineArray(xdata = eta_psi_func, ydata = psi)  
    def psi(self, eta):
        return self.psi_line.get_values(eta, k = 1)
    
    def diff_psi(self,eta):
        return self.psi_line.get_diffs(eta, k = 1)

    def eta(self, psi):
        psi_line = MFnLineArray(xdata = self.psi_line.ydata,
                                ydata = self.psi_line.xdata)
        return psi_line.get_values(psi, k = 1)

# distribution of long segments
    
    def P(self,x,eta):
        eta_linsp = np.linspace(1e-10,eta,200)
        M = 2.*np.trapz(self.psi(eta_linsp)*
            e**(-(x[:,np.newaxis]/self.delta - 1.)*self.psi(eta_linsp)), eta_linsp)
        
        L = self.psi(eta)**2/self.diff_psi(eta) * e**(-(x/self.delta - 2.)*self.psi(eta))
        
        return M.flatten()/eta/self.delta, L/eta/self.delta

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    sigma = 5.0
    eta = np.linspace(0.001,0.74, 500)
    psi = np.linspace(1e-10, 20., 200)
    c = Curtin(sigma = sigma)
    x1 = np.linspace(c.delta,2*c.delta,300)
    x2 = np.linspace(2*c.delta,120,300)
    
#    plt.plot(eta, c.psi(eta), label = 'psi')
#    plt.plot(eta, e**(-2.*0.577216)/(0.7476-eta), label = 'psi_1')
#    plt.plot(eta, 1./(1-eta/0.7476)**0.7476, label = 'psi_2')
#    plt.plot(eta, eta + eta**2 + 7./6.*eta**3 + 13./9.* eta**4, label = 'psi_expansion')
#    plt.plot(eta, c.diff_psi(eta), label = 'diff psi')
#    plt.plot(psi, c.eta(psi), label = 'eta')
    plt.plot(x1,c.P(x1,0.2)[0], label = 'M fragments')
    plt.plot(x2,c.P(x2,0.2)[1], label = 'L fragments')
    print 'integal of P(x,eta) =', np.trapz(np.hstack((c.P(x1,0.5)[0],c.P(x2,0.5)[1])),np.hstack((x1,x2)))
    plt.legend(loc = 'best')
    plt.show() 
    
    
    
    
    
    
    
    
    
    
    
    