'''
Created on Sep 20, 2012

The CompositeCrackBridge class has a method for evaluating fibers and matrix
strain in the vicinity of a crack bridge.
Fiber diameter and bond coefficient can be set as random variables.
Reinforcement types can be combined by creating a list of Reinforcement
instances and defining it as the reinforcement_lst Trait in the
CompositeCrackBridge class.
The evaluation is array based.

@author: rostar
'''
import numpy as np
from etsproxy.traits.api import HasTraits, cached_property, \
    Float, Property, Instance, List, Array, Tuple
from hom_CB_cont_fibers import CrackBridgeContFibers
from hom_CB_short_fibers import CrackBridgeShortFibers
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import Reinforcement, ContinuousFibers, ShortFibers

class CompositeCrackBridge(HasTraits):

    reinforcement_lst = List(Instance(Reinforcement))
    w = Float
    E_m = Float
    Ll = Float
    Lr = Float

    V_f_tot = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_V_f_tot(self):
        V_f_tot = 0.0
        for reinf in self.reinforcement_lst:
            V_f_tot += reinf.V_f
        return V_f_tot

    E_c = Property(depends_on='reinforcement_lst+')
    @cached_property
    def _get_E_c(self):
        E_fibers = 0.0
        for reinf in self.reinforcement_lst:
            E_fibers += reinf.V_f * reinf.E_f
        E_c = self.E_m * (1. - self.V_f_tot) + E_fibers
        return E_c * (1. + 1e-15)

    sorted_reinf_lst = Property(Tuple(List,List), depends_on='reinforcement_lst')
    @cached_property
    def _get_sorted_reinf_lst(self):
        cont_reinf_lst = []
        short_reinf_lst = []
        for reinf in self.reinforcement_lst:
            if reinf.__class__ == ContinuousFibers:
                cont_reinf_lst.append(reinf)
            elif reinf.__class__ == ShortFibers:
                short_reinf_lst.append(reinf)
        return cont_reinf_lst, short_reinf_lst
    
    cont_fibers = Property(Instance(CrackBridgeContFibers), depends_on = 'w,E_m,Ll,Lr,reinforcement_lst+')
    @cached_property
    def _get_cont_fibers(self):
        cbcf = CrackBridgeContFibers(w=self.w,
                                     Ll=self.Ll,
                                     Lr=self.Lr,
                                     E_m=self.E_m,
                                     cont_reinf_lst=self.sorted_reinf_lst[0])
        return cbcf

    short_fibers = Property(Instance(CrackBridgeShortFibers), depends_on = 'w,E_m,Ll,Lr,reinforcement_lst+')
    @cached_property
    def _get_short_fibers(self):
        cbsf = CrackBridgeShortFibers(w=self.w,
                                     Ll=self.Ll,
                                     Lr=self.Lr,
                                     E_m=self.E_m,
                                     short_reinf_lst=self.sorted_reinf_lst[1])
        return cbsf
    
    _x_arr = Property(Array, depends_on = 'w,E_m,Ll,Lr,reinforcement_lst+')
    @cached_property
    def _get__x_arr(self):
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) != 0:
            added_x = np.hstack((self.cont_fibers.x_arr, self.short_fibers.x_arr))
            return np.sort(added_x)
        elif len(self.sorted_reinf_lst[0]) != 0:
            return self.cont_fibers.x_arr
        elif len(self.sorted_reinf_lst[1]) != 0:
            return self.short_fibers.x_arr
        
    _epsm_arr = Property(Array, depends_on = 'w,E_m,Ll,Lr,reinforcement_lst+')
    @cached_property
    def _get__epsm_arr(self):
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) != 0:
            epsm_cont_interp = MFnLineArray(xdata=self.cont_fibers.x_arr, ydata=self.cont_fibers.epsm_arr)
            epsm_short_interp = MFnLineArray(xdata=self.short_fibers.x_arr, ydata=self.short_fibers.epsm_arr)
            added_epsm_cont = self.cont_fibers.epsm_arr + epsm_short_interp.get_values(self.cont_fibers.x_arr) 
            added_epsm_short = self.short_fibers.epsm_arr + epsm_cont_interp.get_values(self.short_fibers.x_arr) 
            sorted_idx = np.argsort(np.hstack((self.cont_fibers.x_arr, self.short_fibers.x_arr)))
            return np.hstack((added_epsm_cont, added_epsm_short))[sorted_idx]
        elif len(self.sorted_reinf_lst[0]) != 0:
            return self.cont_fibers.epsm_arr
        elif len(self.sorted_reinf_lst[1]) != 0:
            return self.short_fibers.epsm_arr
        
    _epsf0_arr = Property(Array, depends_on = 'w,E_m,Ll,Lr,reinforcement_lst+')
    @cached_property
    def _get__epsf0_arr(self):
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) != 0:
            epsf0_cont = self.cont_fibers.epsf0_arr
            epsf0_short = self.short_fibers.epsf0_arr
        elif len(self.sorted_reinf_lst[0]) != 0:
            epsf0_cont = self.cont_fibers.epsf0_arr
            epsf0_short = np.array([])
        elif len(self.sorted_reinf_lst[1]) != 0:
            epsf0_cont = np.array([])
            epsf0_short = self.short_fibers.epsf0_arr
        return epsf0_cont, epsf0_short

    _epsf0_arr_cont = Property(Array, depends_on = 'w,E_m,Ll,Lr,reinforcement_lst+')
    @cached_property
    def _get__epsf0_arr_cont(self):
        return self._epsf0_arr[0]
    
    _epsf0_arr_short = Property(Array, depends_on = 'w,E_m,Ll,Lr,reinforcement_lst+')
    @cached_property
    def _get__epsf0_arr_short(self):
        return self._epsf0_arr[1]

    sigma_c = Property(depends_on = 'w,E_m,Ll,Lr,reinforcement_lst+')
    @cached_property
    def _get_sigma_c(self):
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) != 0:
            sigma_c_cont = np.sum(self._epsf0_arr_cont * self.cont_fibers.sorted_stats_weights *
                              self.cont_fibers.sorted_V_f * self.cont_fibers.sorted_nu_r *
                              self.cont_fibers.sorted_E_f * (1. - self.cont_fibers.damage))
            sigma_c_short = np.sum(self._epsf0_arr_short * self.short_fibers.sorted_V_f *
                               self.short_fibers.sorted_E_f)
        elif len(self.sorted_reinf_lst[0]) != 0:
            sigma_c_cont = np.sum(self._epsf0_arr_cont * self.cont_fibers.sorted_stats_weights *
                              self.cont_fibers.sorted_V_f * self.cont_fibers.sorted_nu_r *
                              self.cont_fibers.sorted_E_f * (1. - self.cont_fibers.damage))
            sigma_c_short = 0.0
        elif len(self.sorted_reinf_lst[1]) != 0:
            sigma_c_cont = 0.0
            sigma_c_short = np.sum(self._epsf0_arr_short * self.short_fibers.sorted_V_f *
                               self.short_fibers.sorted_E_f)
        return sigma_c_cont + sigma_c_short
        
    
if __name__ == '__main__':
    
    from matplotlib import pyplot as plt
    from stats.pdistrib.weibull_fibers_composite_distr import fibers_MC
    from spirrid.rv import RV

    tau_scale = 0.559582502104
    tau_shape = 0.120746270203
    tau_loc = 0.00
    xi_shape = 9.0
    xi_scale = 0.0075

    reinf1 = ContinuousFibers(r=3.5e-3,
                              tau=RV('gamma', loc=tau_loc, scale=tau_scale, shape=tau_shape),
                              V_f=0.1,
                              E_f=200e3,
                              xi=fibers_MC(m=xi_shape, sV0=xi_scale),
                              label='carbon',
                              n_int=100)
    
    reinf2 = ShortFibers(bond_law = 'plastic',
                         r=3.5e-3,
                        tau=.1,
                        V_f=0.01,
                        E_f=200e3,
                        xi=RV('weibull_min', scale=20., shape=5.),
                        snub=0.5,
                        phi=RV('sin2x', scale=1.0, shape=0.0),
                        Lf=14.,
                        label='carbon',
                        n_int=100)

    ccb = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf1, reinf2],
                                 Ll=500.,
                                 Lr=500.,
                                 w=.1)
    
    epsf0_combined = np.hstack((ccb._epsf0_arr[0], ccb._epsf0_arr[1]))
    plt.plot(np.zeros_like(epsf0_combined), epsf0_combined, 'ro', label='maximum')
    plt.plot(ccb._x_arr, ccb._epsm_arr, lw=2, color='blue', label='matrix')
    plt.legend(loc='best')
    plt.ylabel('matrix and fiber strain [-]')
    plt.ylabel('long. position [mm]')
    plt.show()
