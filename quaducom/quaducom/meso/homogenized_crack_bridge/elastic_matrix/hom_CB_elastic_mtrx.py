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
    Float, Property, Instance, List, Array, Tuple, Bool
from hom_CB_cont_fibers import CrackBridgeContFibers
from hom_CB_short_fibers import CrackBridgeShortFibers
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import Reinforcement, ContinuousFibers, ShortFibers
from traits.has_traits import on_trait_change

class CompositeCrackBridge(HasTraits):

    reinforcement_lst = List(Instance(Reinforcement))
    w = Float
    E_m = Float
    Ll = Float
    Lr = Float
    ft = Float
    Gf = Float(1.0)
    w_unld = Float(0.0)
    damage_switch = Bool(True)
    
    @on_trait_change('damage_switch')
    def switch_damage(self):
        self.cont_fibers_instance.damage_switch = self.damage_switch
        self.short_fibers_instance.damage_switch = self.damage_switch
    
    epsm_softening = Property(depends_on='w,ft,Gf,E_m')
    @cached_property
    def _get_epsm_softening(self):
        if self.w >= self.w_unld:
            if self.damage_switch == True:
                self.w_unld += self.w - self.w_unld 
            return self.ft * np.exp(-self.ft/self.Gf * self.w) / self.E_m
        else:
            return self.w / self.w_unld * self.ft * np.exp(-self.ft/self.Gf * self.w_unld) / self.E_m

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
    
    cont_fibers_instance = Instance(CrackBridgeContFibers)
    def _cont_fibers_instance_default(self):
        return CrackBridgeContFibers()

    cont_fibers = Property(Instance(CrackBridgeContFibers), depends_on = 'reinforcement_lst+,Ll,Lr,E_m,w')
    @cached_property
    def _get_cont_fibers(self):
        cbcf = self.cont_fibers_instance
        cbcf.w=self.w
        cbcf.Ll=self.Ll
        cbcf.Lr=self.Lr
        cbcf.E_m=self.E_m
        cbcf.E_c=self.E_c
        cbcf.w_unld=self.w_unld
        cbcf.cont_reinf_lst=self.sorted_reinf_lst[0]
        cbcf.epsm_softening=self.epsm_softening
        return cbcf
    
    short_fibers_instance = Instance(CrackBridgeShortFibers)
    def _short_fibers_instance_default(self):
        return CrackBridgeShortFibers()
    
    short_fibers = Property(Instance(CrackBridgeShortFibers), depends_on = 'reinforcement_lst+,E_m,w')
    @cached_property
    def _get_short_fibers(self):
        cbsf = self.short_fibers_instance
        cbsf.w = self.w
        cbsf.E_m = self.E_m
        cbsf.E_c = self.E_c
        cbsf.short_reinf_lst = self.sorted_reinf_lst[1]
        cbsf.epsm_softening = self.epsm_softening
        return cbsf
    
    _x_arr = Property(Array, depends_on = 'w,E_m,Ll,Lr,reinforcement_lst+')
    @cached_property
    def _get__x_arr(self):
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) != 0:
            added_x = np.hstack((self.cont_fibers.x_arr, self.short_fibers.x_arr))
            sorted_unique_x = np.unique(added_x)
            return sorted_unique_x
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
            sorted_unique_idx = np.unique(np.hstack((self.cont_fibers.x_arr, self.short_fibers.x_arr)), return_index=True)[1]
            return np.hstack((added_epsm_cont, added_epsm_short))[sorted_unique_idx] - self.epsm_softening
        elif len(self.sorted_reinf_lst[0]) != 0:
            return self.cont_fibers.epsm_arr
        elif len(self.sorted_reinf_lst[1]) != 0:
            self.short_fibers.w = self.w
            return self.short_fibers.epsm_arr

    _epsf_arr = Property(Array, depends_on = 'w,E_m,Ll,Lr,reinforcement_lst+')
    @cached_property
    def _get__epsf_arr(self):
        ''' only for continuous reinforcement '''
        if len(self.sorted_reinf_lst[0]) != 0 and len(self.sorted_reinf_lst[1]) == 0:
            self.cont_fibers.w = self.w
            return self.cont_fibers.epsf_arr
        else:
            raise ValueError('epsf can only be computed for continuous fibers')
        
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
        return sigma_c_cont + sigma_c_short + self.epsm_softening * self.E_m * (1.-self.V_f_tot)
        
    
if __name__ == '__main__':
    from matplotlib import pyplot as plt
    from stats.pdistrib.weibull_fibers_composite_distr import fibers_MC
    from spirrid.rv import RV

    tau_scale =1.53419049
    tau_shape = 0.90615
    tau_loc = 0.00
    xi_shape = 8.6
    xi_scale = 1.0114
    ft = 6.0
    Gf = 3.0

    reinf1 = ContinuousFibers(r=3.5e-3,
                              tau=RV('weibull_min', loc=0.01, scale=.1, shape=2.),
                              V_f=0.005,
                              E_f=200e3,
                              xi=fibers_MC(m=7., sV0=0.005),
                              label='carbon',
                              n_int=100)

    reinf_cont = ContinuousFibers(r=3.5e-3,
                          tau=RV('gamma', loc=tau_loc, scale=tau_scale, shape=tau_shape),
                          V_f=0.001,
                          E_f=181e3,
                          xi=fibers_MC(m=xi_shape, sV0=xi_scale),
                          label='carbon',
                          n_int=20)
    
    reinf_short = ShortFibers(bond_law = 'plastic',
                         r=3.5e-3,
                        tau=2.1,
                        V_f=0.03,
                        E_f=200e3,
                        xi=20.,
                        snub=0.5,
                        phi=RV('sin2x', scale=1.0, shape=0.0),
                        Lf=5.,
                        label='carbon',
                        n_int=201)

    ccb = CompositeCrackBridge(E_m=25e3,
                               ft=ft,
                               Gf=Gf,
                               reinforcement_lst=[reinf_cont],
                               Ll=50.,
                               Lr=50.,
                               w=.1)
    
    for i, depsf in enumerate(ccb.cont_fibers.sorted_depsf):
        epsf_x = np.maximum(ccb._epsf0_arr[0][i] - depsf * np.abs(ccb._x_arr), ccb._epsm_arr)
#             if i == 0:
#                 plt.plot(ccb._x_arr, epsf_x, color='blue', label='fibers')
#             else:
        #plt.plot(ccb._x_arr, epsf_x, color='black', alpha=1-0.5*ccb.cont_fibers.damage[i])
    
    
    epsf0_combined = np.hstack((ccb._epsf0_arr[0], ccb._epsf0_arr[1]))
    #plt.plot(np.zeros_like(epsf0_combined), epsf0_combined, 'ro', label='maximum')
    plt.plot(ccb._x_arr, ccb._epsm_arr, lw=2, color='blue', label='matrix')
    #plt.plot(ccb._x_arr, ccb._epsf_arr, lw=2, color='red', label='mean fiber')
    plt.legend(loc='best')
    plt.ylabel('matrix and fiber strain [-]')
#     plt.ylabel('long. position [mm]')
    plt.show()
