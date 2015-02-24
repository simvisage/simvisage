'''
Created on Jan 21, 2015

@author: rostislavrypl
'''

from traits.api import HasTraits, List, Instance, Float, Array, Property, cached_property, implements
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import Reinforcement
import numpy as np
from spirrid import SPIRRID
from spirrid.i_rf import IRF
from spirrid.rf import RF
from scipy.optimize import fminbound
from scipy.interpolate.interpolate import interp2d
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
    E_c = Float
    E_m = Float
        
    x_arr = Property(Array, depends_on='short_reinf_lst+')
    @cached_property
    def _get_x_arr(self):
        Lf_lst = []
        for reinf in self.short_reinf_lst:
            Lf_lst.append(reinf.Lf)
        max_Lf = np.max(np.array(Lf_lst))
        # !!! an even number has to be set as step for the zero position to be in the linspace !!!
        x_arr = np.linspace(-max_Lf/2., max_Lf/2.,61)
        return x_arr

    spirrid_lst = Property(List(Instance(SPIRRID)), depends_on='short_reinf_lst+')
    @cached_property
    def _get_spirrid_lst(self):
        spirrid_epsm_list = []
        for reinf in self.short_reinf_lst:
            cb = CBShortFiberSP()
            spirrid = SPIRRID(q=cb,
                              sampling_type='LHS',
                              theta_vars=dict(tau=reinf.tau,
                                              E_f=reinf.E_f,
                                              r=reinf.r,
                                              xi=reinf.xi,
                                              snub=reinf.snub,
                                              le=reinf.le,
                                              phi=reinf.phi),
                              n_int=reinf.n_int)
            spirrid_epsm_list.append(spirrid)
        return spirrid_epsm_list

    spirrid_evaluation_cached = Property(Array, depends_on='short_reinf_lst+')
    @cached_property
    def _get_spirrid_evaluation_cached(self):
        interpolators_lst = []
        for i, spirr in enumerate(self.spirrid_lst):
            Lfi = self.short_reinf_lst[i].Lf
            def minfunc_short_fibers(w):
                spirr.eps_vars=dict(w=np.array([w]),
                                    x=np.array([0.0]))
                return -spirr.mu_q_arr.flatten()
            w_maxi = fminbound(minfunc_short_fibers, 0.0, Lfi/3., maxfun=20, disp=0)
            w_arri = np.hstack((np.linspace(0.0, w_maxi, 15),np.linspace(w_maxi + 1e-10, Lfi/2., 15)))
            spirr.eps_vars = dict(w=w_arri,
                                  x=self.x_arr)
            interpolators_lst.append(interp2d(self.x_arr, w_arri, spirr.mu_q_arr, fill_value=0.0))
        return interpolators_lst

    epsm_arr = Property(Array, depends_on='short_reinf_lst+,w,E_m')
    @cached_property
    def _get_epsm_arr(self):
        epsm_x_arr = np.zeros(len(self.x_arr))
        for i, interpolator in enumerate(self.spirrid_evaluation_cached):
            sigf_x_i = interpolator(self.x_arr, self.w)
            Ff_x_i = sigf_x_i * self.sorted_V_f[i]
            Fmax_i = np.max(Ff_x_i)
            epsm_x_i = (Fmax_i - Ff_x_i) / self.E_c
            epsm_x_arr += epsm_x_i.flatten()      
        return epsm_x_arr       

    epsf0_arr = Property(Array, depends_on='short_reinf_lst+,w')
    @cached_property
    def _get_epsf0_arr(self):
        return np.array([interpolator(0.,self.w) / self.sorted_E_f[i] for i, interpolator
                         in enumerate(self.spirrid_evaluation_cached)]).flatten()

    sorted_V_f = Property(depends_on='short_reinf_lst+')
    @cached_property
    def _get_sorted_V_f(self):
        return np.array([reinf.V_f for reinf in self.short_reinf_lst])
   
    sorted_E_f = Property(depends_on='short_reinf_lst+')
    @cached_property
    def _get_sorted_E_f(self):
        return np.array([reinf.E_f for reinf in self.short_reinf_lst])        
            
if __name__ == '__main__':
    pass
    
