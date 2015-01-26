'''
Created on Jan 21, 2015

@author: rostislavrypl
'''

from traits.api import HasTraits, List, Instance, Float, Array, Property, cached_property, implements
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import Reinforcement
import numpy as np
from spirrid import SPIRRID, RV
from spirrid.i_rf import IRF
from spirrid.rf import RF
from matplotlib import pyplot as plt


class CBShortFiber(RF):
    '''
    Micromechanical response of a short fiber bridging a crack
    '''

    implements(IRF)

    xi = Float(distr=['weibull_min', 'uniform'])
    E_f = Float(distr=['uniform', 'norm'])
    r = Float(distr=['uniform', 'norm'])
    le = Float(distr=['uniform'])
    tau = Float(distr=['norm', 'uniform', 'weibull_min'])
    snub = Float(distr=['uniform', 'norm'])
    phi = Float(distr=['sin2x', 'uniform'])
    w = Float
    C_code = ''

    def __call__(self, w, tau, r, E_f, le, phi, snub, xi):
        T = 2. * tau / r
        # debonding stage
        ef0_deb = np.sqrt(T * w / E_f)
        # crack opening at which debonding is finished
        w0 = le ** 2 * T / E_f
        # pulling out stage - the fiber is pulled out from the
        # side with the shorter embedded length only
        ef0_pull = (le + w0 - w) * T / E_f
        ef0 = (ef0_deb * (w < w0) + ef0_pull * (w > w0)) * np.exp(phi*snub)
        # include breaking strain
        ef0 = ef0 * (ef0 < xi) * (ef0 > 0.0)
        return ef0
    
class CBShortFiberSP(CBShortFiber):
    '''
    stress profile for a crack bridged by a short fiber
    '''

    x = Float(distr=['uniform'])  

    C_code = ''

    def __call__(self, w, x, tau, r, E_f, le, phi, snub, xi):
        epsf0 = super(CBShortFiberSP, self).__call__(w, tau, r, E_f, le, phi, snub, xi)
        T = 2. * tau / r
        epsf_x = epsf0 / np.exp(snub*phi) - np.abs(x) * T / E_f
        epsf_x = epsf_x * (epsf_x > 0.0) * np.exp(snub * phi)
        return epsf_x * E_f


class CrackBridgeShortFibers(HasTraits):
    
    short_reinf_lst = List(Instance(Reinforcement))
    w = Float
    E_m = Float
    Ll = Float
    Lr = Float
    
    x_arr_lst = Property(Array, depends_on='short_reinf_lst,w,E_m')
    @cached_property
    def _get_x_arr_lst(self):
        x_arr_lst = []
        for reinf in self.short_reinf_lst: 
            if self.Ll < reinf.Lf/2.:
                x_left = np.linspace(-self.Ll, 0.0, 100)
            else:
                x_left = np.hstack((-self.Ll, np.linspace(-reinf.Lf/2., 0.0, 99)))
            if self.Lr < reinf.Lf/2.:
                x_right = np.linspace(1e-10, self.Lr, 100)
            else:
                x_right = np.hstack((np.linspace(1e-10, reinf.Lf/2., 99), self.Lr))
            x_arr_lst.append(np.hstack((x_left, x_right)))
        return x_arr_lst
    
    x_arr = Property(Array, depends_on='short_reinf_lst,w,E_m')
    @cached_property
    def _get_x_arr(self):
        x_arr = np.array([])
        for x in self.x_arr_lst:
            x_arr = np.hstack((x_arr, x))
        return x_arr

    epsf0_arr = Property(Array, depends_on='short_reinf_lst,w,E_m')
    @cached_property
    def _get_epsf0_arr(self):
        epsf0_arr = np.zeros(len(self.short_reinf_lst))
        for i, reinf in enumerate(self.short_reinf_lst):
            cb = CBShortFiber()
            spirrid = SPIRRID(q=cb,
                              sampling_type='PGrid',
                              eps_vars=dict(w=np.array([self.w])),
                              theta_vars=dict(tau=reinf.tau,
                                              E_f=reinf.E_f,
                                              r=reinf.r,
                                              xi=reinf.xi,
                                              snub=reinf.snub,
                                              le=reinf.le,
                                              phi=reinf.phi),
                              n_int=reinf.n_int)
            epsf0_arr[i] = spirrid.mu_q_arr        
        return epsf0_arr

    epsm_arr = Property(Array, depends_on='short_reinf_lst,w,E_m')
    @cached_property
    def _get_epsm_arr(self):
        epsm_x_arr_lst = []
        for i, reinf in enumerate(self.short_reinf_lst):
            cb = CBShortFiberSP()
            spirrid = SPIRRID(q=cb,
                              sampling_type='PGrid',
                              eps_vars=dict(w=np.array([self.w]),
                                            x=self.x_arr_lst[i]),
                              theta_vars=dict(tau=reinf.tau,
                                              E_f=reinf.E_f,
                                              r=reinf.r,
                                              xi=reinf.xi,
                                              snub=reinf.snub,
                                              le=reinf.le,
                                              phi=reinf.phi),
                              n_int=reinf.n_int)
            sigf_x_i = spirrid.mu_q_arr
            Ff_x_i = sigf_x_i * reinf.V_f
            Fmax_i = np.max(Ff_x_i)
            epsm_x_i = (Fmax_i - Ff_x_i) / (1. - reinf.V_f) / self.E_m 
            epsm_x_arr_lst.append(epsm_x_i.flatten())
        return epsm_x_arr_lst[0]
        
    sorted_V_f = Property(depends_on='short_reinf_lst')
    @cached_property
    def _get_sorted_V_f(self):
        sorted_V_f = np.array([])
        for reinf in self.short_reinf_lst:
            sorted_V_f = np.hstack((sorted_V_f,reinf.V_f))
        return sorted_V_f
    
    sorted_E_f = Property(depends_on='short_reinf_lst')
    @cached_property
    def _get_sorted_E_f(self):
        sorted_E_f = np.array([])
        for reinf in self.short_reinf_lst:
            sorted_E_f = np.hstack((sorted_E_f,reinf.E_f))
        return sorted_E_f        
            
if __name__ == '__main__':
    pass
    
