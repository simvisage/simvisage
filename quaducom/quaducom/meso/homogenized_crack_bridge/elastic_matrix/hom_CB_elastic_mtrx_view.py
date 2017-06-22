'''
Created on 16.11.2012

@author: Q
'''
from etsproxy.traits.ui.api import ModelView
from etsproxy.traits.api import Instance, Property, cached_property, Array
import numpy as np
from spirrid.rv import RV
from scipy.optimize import brentq, fminbound
from scipy.integrate import cumtrapz
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge


class CompositeCrackBridgeView(ModelView):

    model = Instance(CompositeCrackBridge)
    results = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_results(self):
        if self.model.w <= 0.0:
            self.model.w = 1e-15
        sigma_c = self.model.sigma_c
#         Kf_broken = np.sum(self.model.cont_fibers.sorted_V_f * self.model.cont_fibers.sorted_nu_r *
#                            self.model.cont_fibers.sorted_stats_weights * self.model.cont_fibers.sorted_E_f *
#                            self.model.cont_fibers.damage)
        if self.model.Ll > self.model.Lr:
            return -self.model._x_arr[::-1], self.model._epsm_arr[::-1], sigma_c, self.model._epsf_arr[::-1]
        else:
            return self.model._x_arr, self.model._epsm_arr, sigma_c, self.model._epsf_arr

    x_arr = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_x_arr(self):
        return self.results[0]
 
    epsm_arr = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_epsm_arr(self):
        return self.results[1]

    epsf_arr = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_epsf_arr(self):
        return self.results[3]

    sigma_c = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_sigma_c(self):
        return self.results[2]

    def sigma_c_arr(self, w_arr, u=False, damage=False):
        sigma_c_lst = []
        u_lst = []
        damage_lst = []
        for i, w in enumerate(w_arr):
            self.model.w = w
            sigma_c_lst.append(self.sigma_c)
            if u == True:
                u_lst.append(self.u_evaluated)
            if damage == True:
                damage_lst.append(np.sum(self.model.cont_fibers.damage *
                                         self.model.cont_fibers.sorted_stats_weights *
                                         self.model.cont_fibers.sorted_nu_r))
        if u == True or damage == True:
            return np.array(sigma_c_lst), np.array(u_lst), np.array(damage_lst)
        else:
            return np.array(sigma_c_lst)

    def secant_K(self, w_arr):
        secant_K_lst = []
        for w_i in w_arr:
            self.model.w = w_i
            secant_K_lst.append(self.model.secant_K)
        return np.array(secant_K_lst)

    u_evaluated = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_u_evaluated(self):
        return self.model.w + np.trapz(self.epsm_arr, self.x_arr)

    sigma_c_max = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_sigma_c_max(self):
        def minfunc_sigma(w):
            self.model.w = w
            stiffness_loss = np.sum(self.model.cont_fibers.Kf * self.model.cont_fibers.damage) / np.sum(self.model.cont_fibers.Kf)
            if stiffness_loss > 0.90:
                return 1. + w
            #plt.plot(w, self.sigma_c, 'ro')
            return -self.sigma_c
        def residuum_stiffness(w):
            self.model.w = w
            stiffness_loss = np.sum(self.model.Kf * self.model.damage) / np.sum(self.model.Kf)
            if stiffness_loss > 0.90:
                return 1. + w
            if stiffness_loss < 0.65 and stiffness_loss > 0.45:
                residuum = 0.0
            else:
                residuum = stiffness_loss - 0.5
            return residuum

        if len(self.model.sorted_reinf_lst[0]) == 0:
            # there are only short fibers
            def minfunc_short_fibers(w):
                self.model.w = w
                return -self.sigma_c
            w_max = fminbound(minfunc_short_fibers, 0.0, 3.0, maxfun=10, disp=0)
            return self.sigma_c, w_max
        else:
            # continuous or mixed fibers
            try:
                w_max = brentq(residuum_stiffness, 0.0, min(0.1 * (self.model.Ll + self.model.Lr), 20.))
            except:
                w_max = 0.03 * (self.model.Ll + self.model.Lr)
            w_points = np.linspace(0, w_max, len(self.model.reinforcement_lst) + 1)
            w_maxima = []
            sigma_maxima = []
            for i, w in enumerate(w_points[1:]):
                w_maxima.append(fminbound(minfunc_sigma, w_points[i], w_points[i + 1], maxfun=10, disp=0))
                sigma_maxima.append(self.sigma_c)
            return sigma_maxima[np.argmax(np.array(sigma_maxima))], w_maxima[np.argmax(np.array(sigma_maxima))]

    def apply_load(self, sigma):
        if sigma > self.sigma_c_max[0]:
            raise ValueError('applied load ', sigma , 'MPa is larger than composite strength ', self.sigma_c_max[0], 'MPa')
        else:
            def residuum(w):
                self.model.w = float(w)
                return sigma - self.sigma_c
            brentq(residuum, 0.0, min(self.sigma_c_max[1], 20.))

    def sigma_f_lst(self, w_arr):
        sigma_f_arr = np.zeros(len(w_arr) * 
                               len(self.model.reinforcement_lst)).reshape(len(w_arr),
                                len(self.model.reinforcement_lst))
        masks = [((self.model.sorted_xi == reinf.xi) * 
                          (self.model.sorted_E_f == reinf.E_f) * 
                          (self.model.sorted_V_f == reinf.V_f))
                 for reinf in self.model.reinforcement_lst]
        for i, w in enumerate(w_arr):
            if w == 0.0:
                self.model.w = 1e-15
            else:
                self.model.w = w
            self.model.damage
            for j, reinf in enumerate(self.model.reinforcement_lst):
                sigma_fi = np.sum(self.model._epsf0_arr * self.model.sorted_stats_weights * self.model.sorted_nu_r * 
                              self.model.sorted_E_f * (1. - self.model.damage) * masks[j])
                sigma_f_arr[i, j] = sigma_fi
        return sigma_f_arr

    Welm = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_Welm(self):
        Km = self.results[4]
        bonded_l = self.epsm_arr[0] ** 2 * Km * (self.model.Ll - np.abs(self.x_arr[0]))
        bonded_r = self.epsm_arr[-1] ** 2 * Km * (self.model.Lr - np.abs(self.x_arr[-1]))
        return 0.5 * (np.trapz(self.epsm_arr ** 2 * Km, self.x_arr) + bonded_l + bonded_r) 

    Welf = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_Welf(self):
        Kf = self.model.E_c - self.results[4]
        bonded_l = self.mu_epsf_arr[0] ** 2 * Kf * (self.model.Ll - np.abs(self.x_arr[0]))
        bonded_r = self.mu_epsf_arr[-1] ** 2 * Kf * (self.model.Lr - np.abs(self.x_arr[-1]))
        return 0.5 * (np.trapz(self.mu_epsf_arr ** 2 * Kf, self.x_arr) + bonded_l + bonded_r)

    W_el_tot = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_W_el_tot(self):
        '''total elastic energy stored in the specimen'''
        return self.Welf + self.Welm

    W_inel_tot = Property(depends_on='model.E_m, model.w, model.Ll, model.Lr, model.reinforcement_lst+')
    @cached_property
    def _get_W_inel_tot(self):
        '''total inelastic energy dissipated during loading up to w'''
        return self.U - self.W_el_tot

    U_line = Property(depends_on='model.E_m, model.Ll, model.Lr, model.reinforcement_lst+, w_arr_energy')
    @cached_property
    def _get_U_line(self):
        '''work done by external force - mfn_line'''
        w_arr = self.w_arr_energy
        u_lst = []
        F_lst = []
        for w in w_arr:
            self.model.w = w
            u_lst.append(self.u_evaluated)
            F_lst.append(self.sigma_c)
        u_arr = np.array(u_lst)
        F_arr = np.array(F_lst)
        U_line = MFnLineArray(xdata=w_arr, ydata=np.hstack((0, cumtrapz(F_arr, u_arr))))
        return U_line

    U = Property(depends_on='model.E_m, model.Ll, model.Lr, model.reinforcement_lst+, model.w')
    @cached_property
    def _get_U(self):
        '''work done by external force U(w)'''
        return self.U_line.get_values(self.model.w)

    w_arr_energy = Array

    def get_sigma_m_x_input(self, sigma):
        self.apply_load(sigma)
        line = MFnLineArray(xdata=self.x_arr,
                            ydata=self.epsm_arr)
        return line.get_values(self.x_input)

if __name__ == '__main__':

    from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import ContinuousFibers, ShortFibers
    from stats.pdistrib.weibull_fibers_composite_distr import fibers_MC
    import matplotlib.pyplot as plt

    tau_scale =1.53419049
    tau_shape = 1.
    tau_loc = 0.00
    xi_shape = 8.6
    xi_scale = .0114

    reinf1 = ContinuousFibers(r=3.5e-3,
                              tau=RV('weibull_min', loc=0.01, scale=.1, shape=2.),
                              V_f=0.005,
                              E_f=200e3,
                              xi=fibers_MC(m=7., sV0=0.005),
                              label='carbon',
                              n_int=100)

    reinf_cont = ContinuousFibers(r=3.5e-3,
                          tau=RV('gamma', loc=tau_loc, scale=tau_scale, shape=tau_shape),
                          V_f=0.00001,
                          E_f=181e3,
                          xi=fibers_MC(m=xi_shape, sV0=xi_scale),
                          label='carbon',
                          n_int=500)

    reinf_short = ShortFibers(bond_law = 'plastic',
                    r=.2,
                    tau=1.,
                    V_f=0.01,
                    E_f=200e3,
                    xi=10.,
                    snub=0.5,
                    phi=RV('sin2x', scale=1.0, shape=0.0),
                    Lf=20.,
                    label='short steel fibers',
                    n_int=50)

    model = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf_cont],
                                 Ll=10.,
                                 Lr=10.,
                                 ft=3.0,
                                 Gf=.3
                                 )

    ccb_view = CompositeCrackBridgeView(model=model)

    def profile(w):
        ccb_view.model.w = w
        plt.plot(ccb_view.x_arr, ccb_view.epsm_arr, lw=2, label='Ll=' + str(ccb_view.model.Ll))
        plt.plot(ccb_view.x_arr, ccb_view.epsf_arr, color='red', lw=2)
        plt.xlabel('position [mm]')
        plt.ylabel('strain')

    def sigma_c_w(w_arr):
        sigma_c_arr, u_arr, damage_arr = ccb_view.sigma_c_arr(w_arr, damage=False, u=True)
        plt.plot(w_arr, sigma_c_arr, lw=2, color='black', label='w-sigma')
        #plt.plot(w_arr, ccb_view.secant_K(w_arr), lw=2, color='blue', label='K secant')
        #plt.plot(w_arr, damage_arr, lw=2, color='red', label='damage')
        #plt.plot(u_arr, sigma_c_arr, lw=2, label='u-sigma')
        #plt.plot(ccb_view.sigma_c_max[1], ccb_view.sigma_c_max[0], 'bo')
        plt.xlabel('w,u [mm]')
        plt.ylabel('$\sigma_c$ [MPa]')
        plt.legend(loc='best')

    def sigma_f(w_arr):
        sf_arr = ccb_view.sigma_f_lst(w_arr)
        for i, reinf in enumerate(ccb_view.model.reinforcement_lst):
            plt.plot(w_arr, sf_arr[:, i], label=reinf.label)

    def energy(w_arr):
        ccb_view.w_arr_energy = w_arr
        Welm = []
        Welf = []
        Wel_tot = []
        U = []
        Winel = []
        u = []
        ccb_view.U_line
        for w in w_arr:
            ccb_view.model.w = w
            Wel_tot.append(ccb_view.W_el_tot)
            Welm.append(ccb_view.Welm)
            Welf.append(ccb_view.Welf)
            U.append(ccb_view.U)
            Winel.append(ccb_view.W_inel_tot)
            u.append(ccb_view.u_evaluated)
        # plt.plot( w_arr, Welm, lw = 2, label = 'Welm' )
        # plt.plot( w_arr, Welf, lw = 2, label = 'Welf' )
        plt.plot(w_arr, Wel_tot, lw=2, color='black', label='elastic strain energy')
        plt.plot(w_arr, Winel, lw=2, ls='dashed', color='black', label='inelastic energy')
        plt.plot(w_arr, U, lw=3, color='red', label='work of external force')
        plt.xlabel('w [mm]')
        plt.ylabel('W')
        plt.ylim(0.0)
        plt.legend(loc='best')

    # TODO: check energy for combined reinf
    # energy(np.linspace(.0, .15, 100))
    #sigma_c = np.linspace(1., 7., 7)
    #profile(0.01)
    w_ld = np.linspace(0.0, 0.1, 50)
    w_unl1 = np.linspace(0.1, 0.0, 50)
    w_rel1 = np.linspace(0.0, 0.2, 50)
    w_unl2 = np.linspace(0.2, 0.0, 50)
    w_rel2 = np.linspace(0.0, 0.5, 50)
    w_unl3 = np.linspace(0.4, 0.0, 50)
    w_rel3 = np.linspace(0.0, 0.6, 50)

    w = np.hstack((w_ld, w_unl1, w_rel1, w_unl2, w_rel2))
    #w = np.linspace(0,0.6,100)
    sigma_c_w(w)
    # energy(w)
    # bundle at 20 mm
    # sigma_bundle = 70e3*w/20.*np.exp(-(w/20./0.03)**5.)
    # plt.plot(w,sigma_bundle)
    # plt.plot(ccb_view.sigma_c_max[1], ccb_view.sigma_c_max[0], 'bo')
    # sigma_f(np.linspace(.0, .16, 50))
    plt.legend(loc='best')
    plt.show()
