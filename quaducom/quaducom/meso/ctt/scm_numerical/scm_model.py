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

from etsproxy.traits.api import \
    HasTraits, Instance, Int, Array, List, Callable, Interface, \
    implements, Trait, cached_property, Property, Float
from stats.spirrid.spirrid import SPIRRID, FunctionRandomization, MonteCarlo
from interpolated_spirrid import InterpolatedSPIRRID, RangeAdaption
from stats.misc.random_field.random_field_1D import RandomField
from operator import attrgetter
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
    '''crack bridge class - includes informations about position, stress range and evaluates
    the stress and strain profiles in the composite components'''
    get_sigma_f_x_reinf = Callable
    randomization = Instance(FunctionRandomization)

    E_m = Property
    @cached_property
    def _get_Em(self):
        return self.randomization.tvars['E_m']

    V_f = Property
    @cached_property
    def _get_V_f(self):
        return self.randomization.tvars['V_f']

    E_f = Property
    @cached_property
    def _get_E_f(self):
        return self.randomization.tvars['E_f']

    position = Float
    Ll = Float
    Lr = Float
    x = Array
    load_sigma_c = Array

    def get_sigma_x_reinf(self):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        return self.get_sigma_f_x_reinf(self.load_sigma_c, self.x, self.Ll, self.Lr)

    def get_eps_x_reinf(self):
        '''
        evaluation of strain profile in the vicinity of a crack bridge
        '''
        return self.get_sigma_x_reinf() / self.E_f

    def get_eps_x_matrix(self):
        '''
        evaluation of strain profile in the vicinity of a crack bridge
        '''
        return self.get_sigma_x_matrix()/self.E_m

    def get_sigma_x_matrix(self):
        '''
        evaluation of stress profile in the vicinity of a crack bridge
        '''
        load_sigma_c = self.load_sigma_c
        sigma_c_ff = np.ones(len(self.x))[:, np.newaxis] * load_sigma_c[np.newaxis, :]
        sigma_x_matrix = (sigma_c_ff - self.get_sigma_x_reinf() * self.V_f)/(1.-self.V_f)
        return sigma_x_matrix

class ICBMFactory(Interface):

    def new_cb_model(self):
        pass

class CBMFactory(HasTraits):

    implements(ICBMFactory)

    randomization = Instance(FunctionRandomization)

class CBMeanFactory(CBMFactory):

    load_sigma_c_max = Float
    load_n_sigma_c = Int
    n_w = Int
    n_x = Int
    n_BC = Int
    # homogenization model
    spirrid = Property(Instance(SPIRRID))
    @cached_property
    def _get_spirrid(self):
        args = self.randomization.trait_get(['q', 'evars', 'tvars'])
        return SPIRRID(**args)
    
    adaption = Instance(RangeAdaption)
    def _adaption_default(self):
        return RangeAdaption(spirrid = self.spirrid,
                             load_sigma_c_max = self.load_sigma_c_max,
                             load_n_sigma_c = self.load_n_sigma_c,
                             n_w = self.n_w,
                             n_x = self.n_x,
                             n_BC = self.n_BC)
        
    interpolated_spirrid = Property()
    @cached_property
    def _get_interpolated_spirrid(self):
        return InterpolatedSPIRRID(spirrid = self.spirrid,
                                   adaption = self.adaption)

    #===========================================================================
    # Construct new crack bridge (the mean response is shared by all instances)
    #===========================================================================
    def new_cb(self):
        return CB(randomization = self.randomization,
                  get_sigma_f_x_reinf = self.interpolated_spirrid
                  )

class CBRandomFactory(CBMFactory):

    #===========================================================================
    # Construct new crack bridge (each crack bridge has its own random realization)
    #===========================================================================
    def new_cb(self):
        sample = CBRandomSample(randomization = self.randomization)
        return CB(get_sigma_f_x_reinf = sample,
                  randomization = self.randomization)

from scipy.interpolate import RectBivariateSpline

class SCM(HasTraits):
    '''Stochastic Cracking Model - compares matrix strength and stress, inserts new CS instances
    at positions, where the matrix strength is lower than the stress; evaluates stress-strain diagram
    by integrating the strain profile along the composite'''
    cb_randomization = Instance(FunctionRandomization)
    
    cb_type = Trait('mean', dict(mean = CBMeanFactory,
                                 random = CBRandomFactory))

    cb_factory = Property(depends_on = 'cb_type')
    @cached_property
    def _get_cb_factory(self):
        return self.cb_type_(randomization = self.cb_randomization,
                             load_sigma_c_max = self.load_sigma_c_max,
                             load_n_sigma_c = self.load_n_sigma_c,
                             n_w = self.n_w,
                             n_x = self.n_x,
                             n_BC = self.n_BC
                             )

    n_w = Int
    n_x = Int
    n_BC = Int
    cb_list = List
    length = Float(desc = 'composite specimen length')
    nx = Int(desc = 'number of discretization points')

    load_sigma_c = Property(depends_on = '+load')
    def _get_load_sigma_c(self):
        #applied external load in terms of composite stress
        return np.linspace(self.load_sigma_c_min, self.load_sigma_c_max, self.load_n_sigma_c)

    load_sigma_c_min = Float(load = True)
    load_sigma_c_max = Float(load = True)
    load_n_sigma_c = Int(load = True)

    x_arr = Property(Array, depends_on = 'length, nx')
    @cached_property
    def _get_x_arr(self):
        # discretizes the specimen length
        return np.linspace(0., self.length, self.nx)

    sigma_mx_load_ff = Property(depends_on = '+load, +cb_randomization.tvars')
    @cached_property
    def _get_sigma_mx_load_ff(self):
        # 2D array of stress in the matrix along an uncracked composite of shape (load_sigma_c, x_arr)
        Em = self.cb_randomization.tvars['E_m']
        Ef = self.cb_randomization.tvars['E_f']
        Vf = self.cb_randomization.tvars['V_f']
        Ec = Ef*Vf + Em*(1.-Vf)
        sigma_m_ff = self.load_sigma_c[:, np.newaxis] * Em / Ec
        sigma_mx_load_ff = np.ones(len(self.x_arr)) * sigma_m_ff
        return sigma_mx_load_ff

    random_field = Instance(RandomField)
    matrix_strength = Property(depends_on = 'random_field.+modified')
    @cached_property
    def _get_matrix_strength(self):
        # evaluates a random field realization and creates a 2D spline reprezentation
        rf = self.random_field.random_field
        field_2D = np.ones_like(self.load_sigma_c[:, np.newaxis]) * rf
        rf_spline = RectBivariateSpline(self.load_sigma_c, self.random_field.xgrid, field_2D)
        return rf_spline(self.load_sigma_c, self.x_arr)

    def sort_cbs(self):
        # sorts the CBs by position and adjusts the boundary conditions
        # sort the CBs
        crack_position = self.cb_list[-1].position
        self.cb_list = sorted(self.cb_list, key = attrgetter('position'))
        #find idx of the new crack
        for i, crack in enumerate(self.cb_list):
            if crack.position == crack_position:
                idx = i
        # specify the boundaries
        if idx != 0:
            # there is a crack at the left hand side
            cbl = self.cb_list[idx - 1]
            cb = self.cb_list[idx]
            cbl.Lr = (cb.position - cbl.position) / 2. 
            cb.Ll = cbl.Lr
        else:
            # the new crack is the first from the left hand side
            self.cb_list[idx].Ll = self.cb_list[idx].position

        if idx != len(self.cb_list) - 1:
            # there is a crack at the right hand side
            cb, cbr = self.cb_list[idx], self.cb_list[idx+1]
            cbr.Ll = (cbr.position - cb.position) / 2.
            cb.Lr = cbr.Ll
        else:
            # the new crack is the first from the right hand side
            self.cb_list[idx].Lr = self.length - self.cb_list[idx].position

        # specify the x range and stress profile for the new crack and its neighbors
        idxs = [idx-1,idx,idx+1]
        if idx == 0:
            idxs.remove(-1)
        if idx == len(self.cb_list) - 1:
            idxs.remove(len(self.cb_list))
        for idx in idxs:
            mask1 = self.x_arr >= (self.cb_list[idx].position - self.cb_list[idx].Ll)
            if idx == 0:
                mask1[0] = True
            mask2 = self.x_arr <= (self.cb_list[idx].position + self.cb_list[idx].Lr)
            self.cb_list[idx].x = self.x_arr[mask1 * mask2] - self.cb_list[idx].position
            crack_position_idx = np.argwhere(self.x_arr == self.cb_list[idx].position)
            left = crack_position_idx - len(np.nonzero(self.cb_list[idx].x < 0.)[0])
            right = crack_position_idx + len(np.nonzero(self.cb_list[idx].x > 0.)[0]) + 1
            self.sigma_m[-len(self.cb_list[idx].load_sigma_c):, left:right] = self.cb_list[idx].get_sigma_x_matrix().T

    def evaluate(self):
        #seek for the minimum strength redundancy to find the position of the next crack
        while np.any(self.sigma_m > self.matrix_strength):
            mask = self.sigma_m < self.matrix_strength
            idx = np.argmin(mask.sum(0))
            f = mask.sum(0)[idx]
            crack_position = self.x_arr[idx]
            load_sigma_c = self.load_sigma_c[f:]
            cbf = self.cb_factory
            new_cb = cbf.new_cb()
            new_cb.position = float(crack_position)
            new_cb.load_sigma_c = load_sigma_c
            self.cb_list.append(new_cb)
            self.sort_cbs()
#            e_arr = orthogonalize([np.arange(len(self.x_arr)), np.arange(len(self.load_sigma_c))])
#            m.surf(e_arr[0], e_arr[1], self.sigma_m)
#            m.surf(e_arr[0], e_arr[1], self.matrix_strength)
#            m.show()
            if self.last_pos == crack_position:
                print 'FAIL - refine the Ll, Lr, w or x ranges'
                break
            self.last_pos = crack_position
    
    last_pos = Float

    sigma_m = Array
    def _sigma_m_default(self):
        return self.sigma_mx_load_ff
