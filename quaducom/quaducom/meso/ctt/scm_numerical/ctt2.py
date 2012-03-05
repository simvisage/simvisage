#-------------------------------------------------------------------------------
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
# Created on Oct 21, 2011 by: rch

from enthought.traits.api import \
    HasTraits, Float, List, Range, \
    Callable, Instance, Trait, Int, \
    cached_property, Property, \
    implements, Interface, Array, WeakRef, \
    DelegatesTo, on_trait_change

from stats.spirrid.spirrid import SPIRRID, FunctionRandomization, MonteCarlo
from stats.spirrid.rv import RV
from interpolated_spirrid import InterpolatedSPIRRID
from stats.misc.random_field.gauss_1D import GaussRandomField
from operator import attrgetter
import pickle



import numpy as np

class CBRandomSample(HasTraits):

    randomization = Instance(FunctionRandomization)

    sampling = Property
    @cached_property
    def _get_sampling(self):
        return MonteCarlo(randomization = self.randomization)

    def __call__(self, e):
        q = self.randomization.q
        sig = q(e, *self.sampling.theta)
        return np.sum(sig)

class CB(HasTraits):

    get_force_x_reinf = Callable
    randomization = Instance(FunctionRandomization)

    Em = Property
    @cached_property
    def _get_Em(self):
        return self.randomization.tvars['E_m']

    Am = Property
    @cached_property
    def _get_Am(self):
        return self.randomization.tvars['A_m']

    Er = Property
    @cached_property
    def _get_Er(self):
        return self.randomization.tvars['E_r']

    Ar = Property
    @cached_property
    def _get_Ar(self):
        return self.randomization.tvars['A_r']

    Nf = Property
    @cached_property
    def _get_Nf(self):
        return self.randomization.tvars['Nf']

    Kr = Property
    @cached_property
    def _get_Kr(self):
        return self.Ar * self.Er * self.Nf

    Km = Property
    @cached_property
    def _get_Km(self):
        return self.Em * self.Am

    position = Float
    Ll = Float
    Lr = Float
    x = Array
    P = Array
    sigma_ff = Float

    def get_eps_x_reinf(self):
        '''
        evaluation of strain profile in the vicinity of a crack bridge
        '''
        return self.get_force_x_reinf(self.P, self.x, self.Ll, self.Lr)[0] / self.Kr

    def get_sigma_x_reinf(self):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        return self.get_force_x_reinf(self.P, self.x, self.Ll, self.Lr)[0] / self.Ar / self.Nf

    def get_eps_x_matrix(self):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        P = self.P
        Pff = np.ones(len(self.x))[:, np.newaxis] * P[np.newaxis, :]
        return (Pff - self.get_force_x_reinf(P, self.x, self.Ll, self.Lr)[0]) / self.Km

    def get_sigma_x_matrix(self):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        return self.get_eps_x_matrix() * self.Em

    def get_crack_width(self):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        return self.get_force_x_reinf(self.P, self.x, self.Ll, self.Lr)[1]


class ICBMFactory(Interface):

    def new_cb_model(self):
        pass

class CBMFactory(HasTraits):

    implements(ICBMFactory)

    randomization = Instance(FunctionRandomization)

class CBMeanFactory(CBMFactory):

    # homogenization model
    spirrid = Property(Instance(SPIRRID))
    @cached_property
    def _get_spirrid(self):
        args = self.randomization.trait_get(['q', 'evars', 'tvars'])
        return SPIRRID(**args)

    interpolated_spirrid = Property()
    @cached_property
    def _get_interpolated_spirrid(self):
        return InterpolatedSPIRRID(spirrid = self.spirrid)

    #===========================================================================
    # Construct new crack bridge (the mean response is shared by all instances)
    #===========================================================================
    def new_cb(self):
        return CB(randomization = self.randomization,
                  get_force_x_reinf = self.interpolated_spirrid
                  )

class CBRandomFactory(CBMFactory):

    #===========================================================================
    # Construct new crack bridge (each crack bridge has its own random realization)
    #===========================================================================
    def new_cb(self):
        sample = CBRandomSample(randomization = self.randomization)
        return CB(get_force_x_reinf = sample,
                  randomization = self.randomization)

class CTT(HasTraits):

    cb_randomization = Instance(FunctionRandomization)

    cb_type = Trait('mean', dict(mean = CBMeanFactory,
                                 random = CBRandomFactory))

    cb_factory = Property(depends_on = 'cb_type')
    @cached_property
    def _get_cb_factory(self):
        return self.cb_type_(randomization = self.cb_randomization)

    cb_list = List

    length = Float(desc = 'composite specimen length')

    nx = Int(desc = 'number of discretization points')

    applied_force = Property(depends_on = '+force')
    def _get_applied_force(self):
        return np.linspace(self.force_min, self.force_max, self.n_force)

    force_min = Float(force = True)
    force_max = Float(force = True)
    n_force = Int(force = True)

    x_arr = Property(Array, depends_on = 'length, nx')
    @cached_property
    def _get_x_arr(self):
        '''discretizes the specimen length'''
        return np.linspace(0., self.length, self.nx)

    sigma_mxP_ff = Property(depends_on = '+force, cb_randomization')
    @cached_property
    def _get_sigma_mxP_ff(self):
        '''2D array of stress in the matrix along an uncracked composite X force array'''
        Em = self.cb_randomization.tvars['E_m']
        Am = self.cb_randomization.tvars['A_m']
        Er = self.cb_randomization.tvars['E_r']
        Ar = self.cb_randomization.tvars['A_r']
        Nf = self.cb_randomization.tvars['Nf']
        sig_ff = self.applied_force[:, np.newaxis] * Em / (Em * Am + Er * Ar * Nf)
        return np.ones(len(self.x_arr)) * sig_ff

    gauss_rf = Instance(GaussRandomField)
    matrix_strength = Property(depends_on = 'gauss_rf.+modified')
    @cached_property
    def _get_matrix_strength(self):
        rf = self.gauss_rf.random_field
        #rf[-self.nx / 5:] *= 3.
        #rf[:self.nx / 5] *= 3.
        return np.ones_like(self.applied_force[:, np.newaxis]) * rf


    def sort_cbs(self, P):
        '''sorts the CBs by position and adjusts the boundary conditions'''
        # sort the CBs
        self.cb_list = sorted(self.cb_list, key = attrgetter('position'))
        # specify the boundaries at the ends (0 and self.length) of the specimen 
        self.cb_list[0].Ll = self.cb_list[0].position
        self.cb_list[-1].Lr = self.length - self.cb_list[-1].position

        # specify the boundaries between the cracks
        for i, cb in enumerate(self.cb_list[:-1]):
            self.cb_list[i].Lr = (self.cb_list[i + 1].position - cb.position) / 2.
        for i, cb in enumerate(self.cb_list[1:]):
            self.cb_list[i + 1].Ll = (cb.position - self.cb_list[i].position) / 2.

        # specify the x range within the specimen length for every crack
        for i, cb in enumerate(self.cb_list):
            mask1 = self.x_arr >= (cb.position - cb.Ll)
            if i == 0:
                mask1[0] = True
            mask2 = self.x_arr <= (cb.position + cb.Lr)
            cb.x = self.x_arr[mask1 * mask2] - cb.position

        for idx, cb in enumerate(self.cb_list):
            #print 'IDX', idx, 'of', len(self.cb_list) - 1, 'position', cb.position, 'force', P[0]
            self.cb_list[idx].P = P
            crack_position_idx = np.argwhere(self.x_arr == cb.position)
            left = crack_position_idx - np.count_nonzero(cb.x < 0.)
            right = crack_position_idx + np.count_nonzero(cb.x > 0.) + 1
            self.sigma_m[-len(P):, left:right] = cb.get_sigma_x_matrix().T

    x_area = Property(depends_on = 'length, nx')
    def _get_x_area(self):
        return  np.ones_like(self.applied_force)[:, np.newaxis] * self.x_arr[np.newaxis, :]

    def evaluate(self):
        '''seek for the minimum strength redundancy to find the position
        of the next crack'''

        while np.any(self.sigma_m > self.matrix_strength):

            mask = self.sigma_m < self.matrix_strength
            idx = np.argmin(mask.sum(0))
            f = mask.sum(0)[idx]
            crack_position = self.x_arr[idx]
            P = self.applied_force[f:]
            cbf = self.cb_factory
            new_cb = cbf.new_cb()
            new_cb.position = float(crack_position)
            self.crack_positions.append(float(crack_position))
            self.crack_positions.sort()
            new_cb.P = P
            self.cb_list.append(new_cb)
            self.sort_cbs(P)
            if self.last_pos == crack_position:
                print 'FAIL'
                break
            self.last_pos = crack_position

    last_pos = Float
    crack_positions = List

    sigma_m = Array
    def _sigma_m_default(self):
        return self.sigma_mxP_ff

    eps_m = Property(Array, depends_on = 'cb_list')
    @cached_property
    def _get_eps_m(self):
        return self.sigma_m / self.cb_randomization.q.E_m

    sigma_r = Property(Array, depends_on = 'cb_list')
    @cached_property
    def _get_sigma_r(self):
        return (self.applied_force[:, np.newaxis] - self.sigma_m * self.cb_randomization.q.A_m) / self.cb_randomization.q.A_r

    eps_r = Property(Array, depends_on = 'cb_list')
    @cached_property
    def _get_eps_r(self):
        return self.sigma_r / self.cb_randomization.q.E_r

    eps_sigma = Property(depends_on = '+force, length, nx, gauss_rf.+modified')
    def _get_eps_sigma(self):

        eps = np.trapz(self.eps_r, self.x_area, axis = 1) / self.length
        sigma = self.applied_force / (self.cb_randomization.q.A_m + self.cb_randomization.q.A_r)
        return eps, sigma

if __name__ == '__main__':

    from quaducom.resp_func.cb_emtrx_clamped_fiber import \
        CBEMClampedFiberSP
    #import enthought.mayavi.mlab as m
    from stats.spirrid import make_ogrid as orthogonalize
    from matplotlib import pyplot as plt

    length = 1000.
    nx = 1000
    gauss_rf = GaussRandomField(lacor = 5.,
                                    xgrid = np.linspace(0., length, nx),
                                    nsim = 1,
                                    mean = 4.,
                                    stdev = .7,
                                    non_negative_check = True
                               )

    rf = CBEMClampedFiberSP()
    rand = FunctionRandomization(q = rf,
         evars = dict(w = np.linspace(0.0, .25, 40),
                       x = np.linspace(-50., 50., 60),
                       Ll = np.linspace(1., 70., 5),
                       Lr = np.linspace(1., 70., 5),
                        ),
         tvars = dict(tau = .3, #RV('uniform', 0.7, 1.0),
                       l = RV('uniform', .0, 20.0),
                       A_r = 5.31e-4,
                       E_r = 72e3,
                       theta = 0.01, #RV('uniform', 0.0, .05),
                       xi = 0.5, #RV( 'weibull_min', scale = 0.017, shape = 5, n_int = 10 ),
                       phi = 1.0,
                       E_m = 30e3,
                       A_m = 50.,
                       Nf = 1700.,
                        ),
         n_int = 10)


#    t = .1
#    Af = 5.31e-4
#    Ef = 72e3
#    Am = 50. / 1700
#    Em = 30e3
#    l = 10.
#    theta = 0.01
#    xi = 0.0179
#    phi = 1.
#    Ll = 50.
#    Lr = 20.
#    Nf = 1.



    ctt = CTT(length = length,
              nx = nx,
              gauss_rf = gauss_rf,
              cb_randomization = rand,
              cb_type = 'mean',
              force_min = 0.1,
              force_max = 200,
              n_force = 3000
              )


    ctt.evaluate()

    def plot():
#        e_arr = orthogonalize([ctt.applied_force, ctt.x_arr])
#        n_e_arr = [ e / np.max(np.fabs(e)) for e in e_arr ]
#
#        scalar1 = ctt.matrix_strength
#        scalar2 = ctt.sigma_m
#    #    scalar3 = ctt.x_area
#
#        n_scalar1 = scalar1 / np.max(np.fabs(scalar1))
#        n_scalar2 = scalar2 / np.max(np.fabs(scalar1))
#    #    n_scalar3 = scalar3 / np.max(np.fabs(scalar3))
#
#        #m.surf(n_e_arr[0], n_e_arr[1], n_scalar1)
#        #m.surf(n_e_arr[0], n_e_arr[1], n_scalar2)
#    #    m.surf(n_e_arr[0], n_e_arr[1], n_scalar3)
#
#        eps, sigma = ctt.eps_sigma
        w_list = []
        for i, crack in enumerate(ctt.cb_list):
            w_list.append(crack.get_crack_width()[-1])

        plt.hist(w_list, bins = 20, normed = True)
        #plt.plot(eps, sigma)
        plt.title('Histogram der Rissbreiten')
        plt.ylabel('PMF')
        plt.xlabel('crack width [mm]')
        plt.legend(loc = 'best')
        plt.show()
        #m.show()

    plot()


    #for crack in range(20):
    #ctt.evaluate()

#    file1 = open('sigma_m.pickle', 'w')
#    file2 = open('eps_m.pickle', 'w')
#    file3 = open('eps_r.pickle', 'w')
#    file4 = open('ld.pickle', 'w')
#    file5 = open('matr_strength.pickle', 'w')
#    file6 = open('no_crack_sigma.pickle', 'w')
#    pickle.dump(ctt.sigma_m, file1)
#    pickle.dump(ctt.eps_m, file2)
#    pickle.dump(ctt.eps_r, file3)
#    pickle.dump(ctt.eps_sigma, file4)
#    pickle.dump(ctt.matrix_strength, file5)
#    pickle.dump(ctt.sigma_mxP_ff, file6)
#    file1.close()
#    file2.close()
#    file3.close()
#    file4.close()
#    file5.close()
#    file6.close()


