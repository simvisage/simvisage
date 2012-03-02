'''
Created on 19.04.2011

Model of a composite specimen loaded in uniaxial tension. Matrix and reinforcement, respectively,
are assumed to be 1D homogeneous coupled by a (linear or nonlinear) bond law.
The overall response is evaluated as:
1) the integrated strain in the reinforcement over the specimen length
2) the sum of individual displacements localized in crack bridges
   plus integrated strain in the uncracked matrix

For the evaluation the CSMModel class requires an instance of a crack bridge model which includes:
1) the method get_eps_x_reinf(P, Ll, Lr, material params) and returning stress/strain profile
   in the reinforcement over the specimen length,
2) the method get_P_w(P, Ll, Lr, material params) and returning force-crack width relation
   for a crack bridge. 

@author: rrypl
'''

from enthought.traits.api import \
    HasTraits, Instance, on_trait_change, Int, Array, Tuple, List

from enthought.traits.api import HasTraits, Float, Property, \
                                cached_property, Range, Button
from enthought.traits.ui.api import View, Item, Tabbed, VGroup, \
                                VSplit, Group
from enthought.traits.ui.menu import OKButton

import numpy as np
import scipy as sp

from quaducom.ctt import ICBM
from stats.misc.random_field.gauss_1D import GaussRandomField
from quaducom.ctt.homogenized_crack_bridges.steel_bar import SteelBar
from quaducom.ctt.homogenized_crack_bridges.multifilament_yarn import MultifilamentYarn
from quaducom.ctt.homogenized_crack_bridges.short_fibers_monte_carlo import ShortFibersMonteCarlo

class SCMModel(HasTraits):
    '''
    class evaluating the l-d digram of a composite specimen loaded in tension;
    see description above for details
    '''

    # instance of the crack bridge class
    # @todo (crime): Introduce ICBModel interface and use it here.
    cb = Instance(ICBM)
    def _cb_default(self):
        return SteelBar(length = self.length)
    #    return MultifilamentYarn( length = self.length )
    #    return ShortFibersMonteCarlo( length = self.length )

    # @todo (crime): Introduce interface class for random field.
    gauss_rf = Instance(HasTraits)
    def _gauss_rf_default(self):
        return GaussRandomField()

    length = Float(1000., auto_set = False, enter_set = True,
                 desc = 'total specimen length [mm]', modified = True)

    nx = Int(500, auto_set = False, enter_set = True,
                 desc = 'number of length discretization points', modified = True)

    applied_force = Range(low = 1e-10, high = 10000.0, value = 0.2,
                         auto_set = False, enter_set = True,
                 desc = 'applied force [N]', ctrl_param = True)

    # RANDOM MATRIX PROPERTIES
    lacor = Float(10)
    mean = Float(1)
    stdev = Float(.2)

    current_stress = Float

    cracks = List

    x_arr = Property(Array, depends_on = 'length, nx')
    @cached_property
    def _get_x_arr(self):
        '''discretizes the specimen length'''
        return np.linspace(0, self.length, self.nx)

    sigma_m_ff = Property(depends_on = '+ctrl_param, +modified, cb.+params')
    @cached_property
    def _get_sigma_m_ff(self):
        '''stress in the matrix in an uncracked composite'''
        return self.applied_force * self.cb.Km / self.cb.Kc / self.cb.Am

    sigma_mx_ff = Property(depends_on = '+ctrl_param, +modified, cb.+params')
    @cached_property
    def _get_sigma_mx_ff(self):
        '''stress in the matrix along an uncracked composite - far field matrix stress'''
        return np.ones(len(self.x_arr)) * self.sigma_m_ff

    sigma_cx_ff = Property(depends_on = '+ctrl_param, +modified, cb.+params')
    @cached_property
    def _get_sigma_cx_ff(self):
        '''composite stress along an uncracked specimen - far field composite stress'''
        return self.sigma_mx_ff * self.cb.Kc / self.cb.Km

    random_field = Property(Array, depends_on = 'gauss_rf.+modified')
    #@cached_property
    def _get_random_field(self):
        '''generates an array of random matrix strength'''
        # where should the matrix parameters be stored???
        self.gauss_rf.lacor = self.lacor
        self.gauss_rf.xgrid = self.x_arr
        self.gauss_rf.nsim = 1
        self.gauss_rf.mean = self.mean
        self.gauss_rf.stdev = self.stdev
        rf = self.gauss_rf.random_field
        #checks all values of the array for negative values
        if rf.any() <= 0:
            print 'random field created matrix strength smaller than 0'
        return rf

    def x_sig(self, crack, Ll, Lr):
        '''
        calls the crack bridge instance cb and uses its method get_eps_x to
        evaluate the strain (stress) profile in the matrix affected by a single crack
        '''
        cb = self.cb
        cb.P = self.applied_force
        cb.Ll = Ll
        cb.Lr = Lr
        sigma_f_x = cb.get_sigma_x_reinf(self.x_arr - crack)
        # if the crack stress is too low, the reinforcement stiffness too high 
        # and the x grid too coarse, strain localization may not be visible
        # which could cause division by zero, therefore the following if condition
        if max(sigma_f_x) == 0:
            single_crack_profile_m_x = self.sigma_mx_ff
        else:
            single_crack_profile_m_x = self.sigma_m_ff - (self.sigma_m_ff / max(sigma_f_x)) * sigma_f_x
        return single_crack_profile_m_x

    sigma_m_x = Property(depends_on = 'applied_force, cracks')
    @cached_property
    def _get_sigma_m_x(self):
        '''creates the matrix stress array given the current load and array of crack positions'''
        # list of crack and symmetry points positions 
        sym = [0.0, self.length]
        sigma_m_x = np.copy(self.sigma_mx_ff)

        # creates the symmetry points at half the distance between cracks
        # boundaries are symmetry points too - equals homogeneous introduction
        # of loading across the composite cross section
        # 2nd option = introduction of the load through the reinforcement
        # this would equal a crack at the boundaries (TODO as option) 
        sym += list(np.diff(self.cracks) / 2. + self.cracks[:-1])
        sym.sort()

        # creates the stress profile in the matrix from
        # stress profiles of individual cracks   
        for i, crack in enumerate(self.cracks):
            # finds the symmetry point left and right from every crack
            Ll = crack - sym[i]
            Lr = sym[i + 1] - crack
            # evaluates single crack stress profile
            single_cr_profile = self.x_sig(crack, Ll, Lr)
            mask = sigma_m_x <= single_cr_profile
            sigma_m_x = sigma_m_x * mask + single_cr_profile * (mask == False)
        return sigma_m_x

    eps_m_x = Property(depends_on = 'applied_force')
    @cached_property
    def _get_eps_m_x(self):
        '''evaluates the matrix strain array given the current load and array of crack positions'''
        return self.sigma_m_x / self.cb.Em

    eps_r_x = Property(depends_on = 'applied_force')
    @cached_property
    def _get_eps_r_x(self):
        '''evaluates the reinforcement strain array given the current load and array of crack positions'''
        return (self.applied_force - self.sigma_m_x * self.cb.Am) / self.cb.Kr

    def crack_development(self):
        '''
        finds the points, where the matrix tensile strength
        is lower than the applied load and adds new cracks
        '''
        while sum(self.sigma_m_x >= self.random_field) > 0.:
            cr_pos = np.argmax(self.sigma_m_x - self.random_field)
            self.cracks.append(self.x_arr[cr_pos])
            self.cracks.sort()

    def load_step(self):
        '''launches the computation'''

        # loading
        if self.sigma_m_ff > self.current_stress:
            self.current_stress = self.sigma_m_ff
            self.crack_development()

        # unloading           
        else:
            self.current_stress <= self.sigma_m_ff

        self.sig_eps[0].append(np.mean(self.eps_r_x))
        self.sig_eps[1].append(self.applied_force / self.cb.Ac)

    @on_trait_change('+modified, cb.+params')
    def reset_history(self):
        ''' if some params are changed, the sigma - eps history is not valid any more and is reseted'''
        self.no_of_cracks = ([0.], [0.])
        self.sig_eps = ([0.], [0.])
        self.cracks = []
        self.current_stress = 0.
        self.applied_force = 1e-10


    sig_eps = Tuple(List([0.0]), List([0.0]))


    no_of_cracks = Tuple(List, List)
    def get_no_of_cracks(self):
        self.no_of_cracks[0].append(self.current_stress)
        self.no_of_cracks[1].append(len(self.cracks))


    def launch(self, max, points):
        self.reset_history()
        for F in np.linspace(0.01, max, points):
            try:
                if self.cb.weakest_cb <= F:
                    print 'SPECIMEN IS BROKEN'
                    break
            except AttributeError:
                    pass
            self.applied_force = F
    ########################################################################################
    '''Axel's additions'''
    ########################################################################################

#    uncracked_length = Property( Float )
#    def _get_uncracked_length( self ):
#        #aspect ratio 
#        return self.length * np.sum( self.sigma_m_x == self.sigma_mx_ff ) / self.nx
#
#    cbs_sum_strain = List
#    P_list = List
#    cbs_responses = Property( Float )
#    def _get_cbs_responses( self ):
#        #sums cbs_responses and uncracked composite strain up
#        U = self.applied_force * self.uncracked_length / self.cb.Kc
#        if len( self.cracks ) > 0:
#            for crack in self.cracks :
#               try:
#                   U += self.cb.P_w
#               except:
#                   U += self.cb._get_P_w( crack )
#        self.cbs_sum_strain.append( U / self.length )
#        self.P_list.append( self.applied_force / self.cb.Ac )


    traits_view = View(
                           VGroup(
                               Item('length', resizable = False),
                               Item('nx', resizable = False),
                               label = 'CCT parameters',
                               dock = 'tab',
                               id = 'scm.model.numerical.params',
                                    ),
                            VGroup(
                               Item('@cb', resizable = False, show_label = False),
                               label = 'crack bridge parameters',
                               dock = 'tab',
                               id = 'scm.model.numerical.cb',
                                    ),

                        id = 'scm.model.numerical',
                        dock = 'tab',
                        scrollable = False,
                        resizable = True,
                        buttons = [OKButton],
                        height = 0.8, width = 0.8
                        )

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    import enthought.mayavi.mlab as m
    from stats.spirrid import orthogonalize

    P = np.linspace(1, 8000, 500)

    cs = []
    scm = SCMModel(mean = 6.0, stdev = .6, lacor = 1.)
    ni = 0

    x = scm.x_arr
    profiles = []

    for p in P:
        print p
        scm.applied_force = p
        scm.load_step()
        profiles.append(scm.eps_r_x)

    e_arr = orthogonalize([x, P])
    n_e_arr = [ e / np.max(np.fabs(e)) for e in e_arr ]

    mu_q_arr = np.array(profiles)

    n_mu_q_arr = mu_q_arr / np.max(np.fabs(mu_q_arr))

    print n_mu_q_arr.shape
    print n_e_arr[0].shape, n_e_arr[1].shape

    m.surf(n_e_arr[0], n_e_arr[1], n_mu_q_arr[24:])

    m.show()

#    for sim in range(5):
#        ni += 1
#        print ni
#        for p in P:
#            print p
#            scm.applied_force = p
#            scm.load_step()
#        plt.plot(scm.sig_eps[0], scm.sig_eps[1], color = 'k',)
#        cs += list(np.diff(scm.cracks))
#        scm.reset_history()
#        scm.gauss_rf.reevaluate = True
#
#    plt.show()
#    plt.hist(cs, bins = 50, normed = True)
#    plt.show()

