'''
Created on 12 May 2013

@author: Q
'''

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
    HasTraits, Instance, Int, Array, List, \
    cached_property, Property, Float
from quaducom.meso.scm.numerical.interdependent_fibers.representative_CB import RepresentativeCB
from stats.misc.random_field.random_field_1D import RandomField
from operator import attrgetter
import numpy as np
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from scipy.optimize import brentq, newton
import copy
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge
from spirrid.rv import RV
from math import pi
import time as t
import matplotlib.pyplot as plt


class CB(HasTraits):
    '''crack bridge class - includes informations about position,
    stress range and evaluates the stress and strain profiles in the
    composite components'''
    representative_cb= Instance(RepresentativeCB)
    position = Float
    Ll = Float
    Lr = Float
    x = Array
    crack_load_sigma_c = Float

    max_sigma_c = Property(Float, depends_on='Ll, Lr')
    @cached_property
    def _get_max_sigma_c(self):
        '''
        evaluates the strength of the particular crack bridge
        '''
        return self.representative_cb.interpolate_max_sigma_c(self.Ll, self.Lr)

    def get_w(self, load):
        '''
        evaluates the crack width at given load
        '''
        return self.representative_cb.interpolate_w(self.Ll, self.Lr, load)
            
    def get_epsm_x(self, load):
        '''
        evaluates the matrix strain profile at given load
        '''
        if load > self.max_sigma_c:
            return np.nan * np.ones_like(self.x)
        else:
            return self.representative_cb.interpolate_epsm(self.Ll, self.Lr, load, self.x)

    def get_epsf_x(self, load):
        '''
        evaluates the mean fiber strain profile at given load
        '''
        if load > self.max_sigma_c:
            return np.nan * np.ones_like(self.x)
        else:
            return self.representative_cb.interpolate_epsf(self.Ll, self.Lr, load, self.x)
            

class SCM(HasTraits):
    '''Stochastic Cracking Model - compares matrix strength and stress,
    inserts new CS instances at positions, where the matrix strength
    is lower than the stress; evaluates stress-strain diagram
    by integrating the strain profile along the composite'''

    length = Float(desc='composite specimen length')
    nx = Int(desc='# of discretization points for the whole specimen')
    CB_model = Instance(CompositeCrackBridge)
    load_sigma_c_arr = Array
    n_w_CB = Int(100, desc='# of discretization points for w')
    n_BC_CB = Int(15, desc='# of discretization points for boundary conditions')

    representative_cb = Instance(RepresentativeCB)
    def _representative_cb_default(self):
        return RepresentativeCB(CB_model=self.CB_model,
                            load_sigma_c_arr=self.load_sigma_c_arr,
                            length=self.length, n_w=self.n_w_CB, n_BC=self.n_BC_CB
                            )

    # list of composite stress at which cracks occur
    cracking_stress_lst = List
    # list of cracks positions
    crack_positions_lst = List
    
    x_arr = Property(Array, depends_on='length, nx')
    @cached_property
    def _get_x_arr(self):
        # discretizes the specimen length
        return np.linspace(0., self.length, self.nx)

    random_field = Instance(RandomField)

    matrix_strength = Property(depends_on='random_field.+modified')
    @cached_property
    def _get_matrix_strength(self):
        # evaluates a random field
        # realization and creates a spline reprezentation
        rf = self.random_field.random_field
        rf_spline = MFnLineArray(xdata=self.random_field.xgrid, ydata=rf)
        return rf_spline.get_values(self.x_arr)

    def get_cbs_lst(self, positions, cracking_stresses):
        # sorts the CBs by position and adjusts the boundary conditions
        # sort the CBs
        positions = np.array(positions)
        cracking_stresses = np.array(cracking_stresses)
        argsort = np.argsort(positions)
        sorted_positions = positions[argsort]
        sorted_cracking_stresses = cracking_stresses[argsort]
        
        cb_lst = []
        for i, position_i in enumerate(list(sorted_positions)):
            cb_i = CB(position=position_i,
                      representative_cb=self.representative_cb,
                    crack_load_sigma_c=sorted_cracking_stresses[i])
            # specify the boundaries
            if i == 0:
                # the leftmost crack
                cb_i.Ll = cb_i.position
            else:
                # there is a crack at the left hand side
                cb_i.Ll = (cb_i.position - cb_lst[-1].position)/2.
                cb_lst[-1].Lr = cb_i.Ll
            if i == len(positions) - 1:
                # the rightmost crack
                cb_i.Lr = self.length - cb_i.position
            cb_lst.append(cb_i)
            
        for cb_i in cb_lst: 
            # specify the x range for cracks
            mask_right = self.x_arr >= (cb_i.position - cb_i.Ll)
            mask_left = self.x_arr <= (cb_i.position + cb_i.Lr)
            cb_i.x = self.x_arr[mask_left * mask_right] - cb_i.position           
        return cb_lst

    def get_current_cracking_state(self, load):
        '''Creates the list of CB objects that have been formed up to the current load
        '''
        idx_load_level = np.sum(np.array(self.cracking_stress_lst) <= load)
        cb_lst = self.get_cbs_lst(self.crack_positions_lst[:idx_load_level],
                                      self.cracking_stress_lst[:idx_load_level])
        return cb_lst
    
    def get_current_strnegth(self, load):
        if len(self.crack_positions_lst) is not 0:
            cb_lst = self.get_current_cracking_state(load)
            strengths = np.array([cb_i.max_sigma_c for cb_i in cb_lst])
            return np.min(strengths)
        else:
            return np.inf

    def sigma_m(self, load):
        Em = self.CB_model.E_m
        Ec = self.CB_model.E_c
        sigma_m = load * Em / Ec * np.ones(len(self.x_arr))
        if len(self.crack_positions_lst) != 0:
            cb_lst = self.get_current_cracking_state(load)
            for cb_i in cb_lst:
                crack_position_idx = np.argwhere(self.x_arr == cb_i.position)
                idx_l = crack_position_idx - len(np.nonzero(cb_i.x < 0.)[0])
                idx_r = crack_position_idx + len(np.nonzero(cb_i.x > 0.)[0]) + 1
                sigma_m[idx_l:idx_r] = cb_i.get_epsm_x(float(load)) * Em
        return sigma_m
    
    def sigma_m_given_crack_lst(self, load):
        Em = self.CB_model.E_m
        Ec = self.CB_model.E_c
        sigma_m = load * Em / Ec * np.ones(len(self.x_arr))
        if len(self.crack_positions_lst) != 0:
            cb_lst = self.get_cbs_lst(self.crack_positions_lst, self.cracking_stress_lst)
            for cb_i in cb_lst:
                crack_position_idx = np.argwhere(self.x_arr == cb_i.position)
                idx_l = crack_position_idx - len(np.nonzero(cb_i.x < 0.)[0])
                idx_r = crack_position_idx + len(np.nonzero(cb_i.x > 0.)[0]) + 1
                sigma_m[idx_l:idx_r] = cb_i.get_epsm_x(float(load)) * Em
        return sigma_m
    

    def residuum(self, q):
        '''Callback method for the identification of the
        next emerging crack calculated as the difference between
        the current matrix stress and strength. See the scipy newton call below.
        '''
        residuum = np.min(self.matrix_strength - self.sigma_m_given_crack_lst(q))
        return residuum
    
    def evaluate(self):
        # seek for the minimum strength redundancy to find the position
        # of the next crack
        last_pos = pi
        sigc_min = 0.0
        while True:
            try:
                s = t.clock()
                sigc_min = newton(self.residuum, sigc_min)
                try:
                    sigc_min = brentq(self.residuum, 0.0, sigc_min - 1e-10)
                    print 'another root found!!!'
                except:
                    pass
                print 'evaluation of the matrix crack #'+str(len(self.cracking_stress_lst) + 1), t.clock() - s, 's'
            except:
                print 'composite saturated'
                break
            print 'current strength = ', self.get_current_strnegth(sigc_min)
            crack_position = self.x_arr[np.argmin(self.matrix_strength - 
                                                  self.sigma_m(sigc_min))]
            self.crack_positions_lst.append(crack_position)
            self.cracking_stress_lst.append(sigc_min - 1e-10)
            #plt.plot(self.x_arr, self.sigma_m(sigc_min) / self.CB_model.E_m, color='blue', lw=2)
            #plt.plot(self.x_arr, self.matrix_strength / self.CB_model.E_m, color='black', lw=2)
            #plt.show()
            if float(crack_position) == last_pos:
                print last_pos
                raise ValueError('''got stuck in loop,
                try to adapt x, w, BC ranges''')
            last_pos = float(crack_position)
                

if __name__ == '__main__':
    from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import ContinuousFibers, ShortFibers
    from stats.pdistrib.weibull_fibers_composite_distr import fibers_MC
    from matplotlib import pyplot as plt
    length = 100.
    nx = 500
    random_field = RandomField(seed=True,
                               lacor=1.,
                               length=length,
                               nx=500,
                               nsim=1,
                               loc=.0,
                               shape=15.,
                               scale=4.0,
                               distr_type='Weibull'
                               )

    reinf1 = ContinuousFibers(r=3.5e-3,
                              tau=RV('weibull_min', loc=0.01, scale=.01, shape=2.),
                              V_f=0.05,
                              E_f=200e3,
                              xi=fibers_MC(m=7., sV0=0.01),
                              label='carbon',
                              n_int=100)
    
    reinf2 = ContinuousFibers(r=3.5e-3,
                              tau=RV('weibull_min', loc=0.01, scale=.02, shape=2.),
                              V_f=0.008,
                              E_f=200e3,
                              xi=fibers_MC(m=7., sV0=0.011),
                              label='carbon',
                              n_int=100)   

    reinf_short = ShortFibers(bond_law = 'plastic',
                        r=.2,
                        tau=1.5,
                        V_f=0.01,
                        E_f=200e3,
                        xi=10.,
                        snub=0.5,
                        phi=RV('sin2x', scale=1.0, shape=0.0),
                        Lf=20.,
                        label='carbon',
                        n_int=50)

    CB_model = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf1],
                                 )

    scm = SCM(length=length,
              nx=nx,
              random_field=random_field,
              CB_model=CB_model,
              load_sigma_c_arr=np.linspace(0.01, 20., 100),
              n_BC_CB=3
              )

    scm.evaluate()
