'''
Created on Aug 17, 2011

@author: rostar
'''
import numpy as np
from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Range, Bool, on_trait_change, Int, Array, Tuple, \
    List
from enthought.traits.ui.api import  VGroup, View, Item
from math import pi as Pi
from quaducom.resp_func.cb_short_fiber import CBShortFiber
from scipy.stats import  norm
from stats.pdistrib.sin2x_distr import sin2x
from scipy.stats import weibull_min
from scipy.interpolate import interp1d
from math import e


def H(x):
    return np.sign(np.sign(x) + 1.)

class ShortFibersMonteCarlo(HasTraits):

    # applied force
    P = Float(modified = True) # [N]

    # closest crack from left
    Ll = Float(modified = True) # [mm]

    # closest crack from right
    Lr = Float(modified = True) # [mm]

    Er = Float(200e3, auto_set = False, enter_set = True,
               desc = 'steel modulus of elasticity [N/mm^2]', modified = True, param = True)

    Em = Float(30e3, auto_set = False, enter_set = True,
               desc = 'matrix modulus of elasticity [N/mm^2]', modified = True, param = True)

    tau = Float(2.3, auto_set = False, enter_set = True,
               desc = 'sheer force per unit area [N/mm^2]', modified = True, param = True)


    height = Float (30., desc = 'total specimen height [mm]',
                     auto_set = False, enter_set = True, modified = True)
    width = Float(30. , desc = 'total specimen width [mm]',
                   auto_set = False, enter_set = True, modified = True)

    Vf = Float(3., auto_set = False, enter_set = True,
             desc = 'volume fraction of steel fibers [%]', modified = True, param = True)

    r = Float(0.075, auto_set = False, enter_set = True,
             desc = 'fiber radius[mm]', modified = True, param = True)

    f = Float(1., auto_set = False, enter_set = True, # [mm]
                 desc = 'snubbing coefficient', modified = True)
    discr = Float(400.)

    fiber_length = Float(9., desc = 'in mm')
    length = Float(400.)

    with_tau_distr = Bool(True , desc = 'enables Monte-Carlo simulation on every fiber"s tau', modified = True, param = True)
    weibull_tau_shape = Float(8.50666480278)
    weibull_tau_scale = Float(2.49019483074)
    with_f_distr = Bool(False , desc = 'enables Monte-Carlo simulation on every fiber"s snubbing coefficent')
    specimen_broken = Bool(False)

    Ac = Property(depends_on = 'height,width')
    @cached_property
    def _get_Ac(self):
        return self.height * self.width


    Ar = Property(depends_on = 'Vf')
    @cached_property
    def _get_Ar(self):
        return self.Ac * self.Vf / 100.

    Am = Property(depends_on = 'Ac, Vf')
    @cached_property
    def _get_Am(self):
        return self.Ac - self.Ar

    Kr = Property(depends_on = 'Vf, Er')
    @cached_property
    def _get_Kr(self):
        return self.Ar * self.Er

    Km = Property(depends_on = 'Ac, Vf, Em')
    @cached_property
    def _get_Km(self):
        return self.Am * self.Em

    Kc = Property(depends_on = 'Ac, Vf, Er, Em')
    @cached_property
    def _get_Kc(self):
        return self.Kr + self.Km

    #################################################################################


    specimen_broken = Bool(False)
    cbs_pulloutlist = List
    weakest_cb = Float(1e15)
    fiber_data_tuple = Tuple(desc = 'Tuple including fiber data no, le and phi arrays')
    fiber_data_tuple = [[], [], [], [], []]
    crack_list = []




    w = Property (Array , depends_on = 'fiber_length,tau,Er,r')
    @cached_property
    def _get_w(self):
        #defines x_array for pullout data 
        return np.linspace(0, 1.05 * self.fiber_length ** 2. * 1.5 * self.tau / (self.Er * self.r * 2), self.discr)


    all_resp_list = List
    invers_list = List
    def _get_P_w(self, crack_pos):
        #to be worked out
        return self.invers_list[self.crack_list.index(crack_pos)](self.P)



    def cbs_allocation(self, crack_x):
        if crack_x in self.crack_list:
            pass
        else:
            self.crack_list.append(crack_x)
            Cbsf = CBShortFiber()
            self.give_cb_attr()
            le_array, phi_array, tau_array = self.fiber_data_tuple[2][self.crack_list.index(crack_x)] , self.fiber_data_tuple[1][self.crack_list.index(crack_x)], self.fiber_data_tuple[4][self.crack_list.index(crack_x)]
            le_array = le_array.reshape(len(le_array)  , 1)
            phi_array = phi_array.reshape(len(phi_array), 1)
            tau_array = tau_array.reshape(len(tau_array), 1)
            w = self.w.reshape(1, self.discr)
            all_resp = Cbsf(w, tau_array, self.fiber_length, 2 * self.r, self.Er, le_array, phi_array, self.f, 0, 0, 1e15)
            resp = np.sum(all_resp, axis = 0)
            resp = resp[0:np.argmax(resp) + 2  ]
            if resp.max() <= self.weakest_cb:
                self.weakest_cb = resp.max()
            #Lists with results
            function = interp1d(resp, self.w[0:np.argmax(resp) + 2 ])
            self.invers_list.append(function)

    def weakest_cb_check(self):
            print 'weakest cb', self.weakest_cb, 'current force', self.P
            if self.weakest_cb <= self.P:
                print 'SPECIMEN BROKEN'
                #PROGR SHOULD STOP
                self.specimen_broken = True


    diff_list = List
    active_fibers_list = List
    unsorted_crack_list = List
    ar_list = List
    def get_eps_x_reinf (self, crack_x):
        self.cbs_allocation(abs(crack_x[0]))
        self.weakest_cb_check()
        if  crack_x[0]  not in self.unsorted_crack_list:
            self.unsorted_crack_list.append(crack_x[0])
            le_array, phi_array, f_array , tau_array = self.fiber_data_tuple[2][self.crack_list.index(abs(crack_x[0]))] , self.fiber_data_tuple[1][self.crack_list.index(abs(crack_x[0]))], self.fiber_data_tuple[3][self.crack_list.index(abs(crack_x[0]))], self.fiber_data_tuple[4][self.crack_list.index(abs(crack_x[0]))]
            diff = self.sum_of_fiberforces(crack_x, le_array, phi_array , f_array, tau_array)
            self.diff_list.append(diff)
        #print len( self.unsorted_crack_list ), len( self.diff_list )
        index_lists = self.unsorted_crack_list.index(crack_x[0])
        F = self.P - self.diff_list[index_lists]
        #Ar_fibers = self.active_fibers_list[index_lists] * Pi * self.r **2
        Ar_fibers = self.ar_list[index_lists]
        eps = F / self.Er / Ar_fibers
        print 'max', max(self.diff_list[index_lists])

        '''
        plt.plot( crack_x, F, 'k' )
        plt.show()
        plt.plot( crack_x, eps, 'k' )
        plt.show()
        plt.plot( crack_x, self.active_fibers_list[index_lists] , 'k', label = 'aktive Fasern der Rissbruecke' )
        plt.xlim( -10, 10 )
        plt.legend()
        plt.show()
        '''
        return eps * H(eps)

    def get_sigma_x_reinf(self, x):
        '''
        evaluation of strain profile in the vicinity of a crack bridge
        '''
        return self.get_eps_x_reinf(x) * self.Er


    def sum_of_fiberforces(self, crack_x , le_array , phi_array , f_array  , tau_array):
        #creates a 'True' , 'False' matrix. Every true entry is a active fiber. returns force in fibers in array
        snubbing_fact = e ** (self.f * phi_array)
        crack_x = crack_x.reshape(len(crack_x), 1)

        active_fibers_matrix = np.abs(crack_x) <= le_array
        phi_matrix = active_fibers_matrix * phi_array
        rv = self.r / np.sin(Pi / 2. - phi_matrix)
        mask1 = rv <= self.fiber_length / 2.
        mask2 = rv >= self.fiber_length / 2.
        rv *= mask1
        rv += mask2 * self.fiber_length / 2.
        self.ar_list.append(np.sum(rv * self.r * Pi, axis = 1))
        active_fibers_array = np.sum (active_fibers_matrix, axis = 1)
        mask = (active_fibers_array <= np.ones(len(active_fibers_array)))
        active_fibers_array += mask
        self.active_fibers_list.append(active_fibers_array)
        T = tau_array * self.r * Pi * 2.
        product_matrix = active_fibers_matrix * snubbing_fact * (abs(crack_x) - le_array) * T
        result_arr = np.sum (product_matrix, axis = 1)
        result_arr += np.min(result_arr) * -1
        return result_arr



    distr_of_fibers = Property (Tuple)
    @cached_property
    def _get_distr_of_fibers(self):
        #Giving mean,stdev of Fibers
        specimen_volume = self.length * self.width * self.height
        no_of_fibers_in_specimen = (specimen_volume * self.Vf / 100) / (Pi * self.r ** 2 * self.fiber_length)
        prob_crackbridging_fiber = .5 * self.fiber_length / self.length
        mean = prob_crackbridging_fiber * no_of_fibers_in_specimen
        stdev = (prob_crackbridging_fiber * no_of_fibers_in_specimen * (1 - prob_crackbridging_fiber)) ** 0.5
        return [mean, stdev]



    def give_cb_attr(self):
        #gives every crack bridge attributes
        no_of_fibers = norm(self.distr_of_fibers[0], self.distr_of_fibers[1]).ppf(np.random.rand(1))
        phi_array = sin2x._ppf(np.random.rand(no_of_fibers))
        le_array = (np.random.rand(no_of_fibers) * self.fiber_length / 2.)
        if self.with_tau_distr:
            tau_array = abs(weibull_min(self.weibull_tau_shape, self.weibull_tau_scale).ppf(np.random.rand(no_of_fibers)))

        else:
            tau_array = self.tau


        if self.with_f_distr:
            f_array = abs(norm(self.mean_f, self.stdev_f).ppf(np.random.rand(no_of_fibers)))

        else: f_array = self.f
        self.fiber_data_tuple[0].append(int(no_of_fibers))
        self.fiber_data_tuple[1].append(phi_array)
        self.fiber_data_tuple[2].append(le_array)
        self.fiber_data_tuple[3].append(f_array)
        self.fiber_data_tuple[4].append(tau_array)
        #return phi_array, le_array , f_array, tau_array



    traits_view = View(
                       VGroup(
                           Item('r', label = 'fiber radius', resizable = False, springy = True),
                           Item('Er', resizable = False, springy = False),
                           Item('Em', resizable = False, springy = False),
                           Item('tau', resizable = False, springy = False),
                           Item('height', label = 'specimen height', resizable = False, springy = False),
                           Item('width', label = 'specimen width', resizable = False, springy = False),
                           Item('f', label = 'snubbing coefficient f', resizable = False, springy = False),
                           Item('fiber_length', resizable = False, springy = False),
                           Item('with_tau_distr', label = 'with tau gauss distr', springy = False),
                           Item('weibull_tau_shape', label = 'mean tau', springy = False),
                           Item('weibull_tau_scale', label = 'Stdev tau', springy = False),
                           Item('with_f_distr', label = 'with f gauss distr', springy = False),
                           Item('mean_f', label = 'mean f', springy = False),
                           Item('stdev_f', label = 'Stdev f', springy = False),

                           springy = True,
                           label = 'CB parameters',
                           dock = 'tab',
                           id = 'cb.steel_bar.params',
                        ),
                            id = 'cb.steel_bar',
                            dock = 'fixed',
                            scrollable = True,
                            resizable = True,
                            height = 0.8, width = 0.8
                                   )
if __name__ == '__main__':
    from matplotlib import pyplot as plt
    sb = ShortFibersMonteCarlo(Ll = 60., tau = 3.5 , P = 4000)
    #sb.configure_traits()
    x = np.linspace(0, 20, 200)
    #print x[35]
    x -= x[70]

    #print x
    #crack_x = x -
    #x = np.linspace( -10 , 10 , 100 )
#    eps = sb.get_eps_x_reinf( x )
#
#    plt.plot( x, eps, lw = 2, color = 'black' )
#    plt.show()

    import enthought.mayavi.mlab as m
    from stats.spirrid import orthogonalize

    resolution = 200
    profiles_list = []
    forces_arr = np.linspace(0, 6200, resolution)
    for i, p in enumerate(forces_arr):
        sb.P = p
        profiles_list.append(sb.get_eps_x_reinf(x))

    profiles_arr = np.array(profiles_list)

    param_arr = orthogonalize([x, forces_arr])
    norm_param_arr = [ e / np.max(np.fabs(e)) for e in param_arr ]

    norm_profiles_arr = profiles_arr / np.max(np.fabs(profiles_arr))
    m.surf(norm_param_arr[0], norm_param_arr[1], norm_profiles_arr)
    m.show()
