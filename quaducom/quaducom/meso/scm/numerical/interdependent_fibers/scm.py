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
from stats.misc.random_field.random_field_1D import RandomField
import numpy as np
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from scipy.optimize import brentq, newton, root
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge
from spirrid.rv import RV
from math import pi
import time as t
import matplotlib.pyplot as plt


class CB(HasTraits):
    '''crack bridge class - includes informations about position,
    stress range and evaluates the stress and strain profiles in the
    composite components'''
    position = Float
    x = Array
    cracking_stress = Float
    
    CB_model = Instance(CompositeCrackBridge)
    
    def get_sigma_c(self, w):
        self.CB_model.w = w
        return self.CB_model.sigma_c
            
    def get_epsm_x(self, w):
        '''
        evaluates the matrix strain profile at given load
        '''
        self.CB_model.w = w
        if self.CB_model.Ll > self.CB_model.Lr:
            epsm_interp = MFnLineArray(xdata=-self.CB_model._x_arr[::-1], ydata=self.CB_model._epsm_arr[::-1])
        else:
            epsm_interp = MFnLineArray(xdata=self.CB_model._x_arr, ydata=self.CB_model._epsm_arr)
        return epsm_interp.get_values(self.x)

    def get_epsf_x(self, w):
        '''
        evaluates the mean fiber strain profile at given load
        '''
        self.CB_model.w = w
        if self.CB_model.Ll > self.CB_model.Lr:
            epsf_interp = MFnLineArray(xdata=-self.CB_model._x_arr[::-1], ydata=self.CB_model._epsf_arr[::-1])
        else:
            epsf_interp = MFnLineArray(xdata=self.CB_model._x_arr, ydata=self.CB_model._epsf_arr)
        return epsf_interp.get_values(self.x)


class SCM(HasTraits):
    '''Stochastic Cracking Model - compares matrix strength and stress,
    inserts new CS instances at positions, where the matrix strength
    is lower than the stress; evaluates stress-strain diagram
    by integrating the strain profile along the composite'''

    length = Float(desc='composite specimen length')
    nx = Int(desc='# of discretization points for the whole specimen')
    CB_model = Instance(CompositeCrackBridge)
    load_sigma_c_arr = Array

    # list of composite stress at which cracks occur
    cracking_stresses_lst = List
    # list of cracks positions
    crack_positions_lst = List
    # list of cracking strains
    cracking_strains_lst = List
    
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

    CB_objects_lst = List

    current_CB_lst = Property(List, depends_on='crack_positions_lst')
    @cached_property
    def _get_current_CB_lst(self):
        # sorts the CBs by position and adjusts the boundary conditions
        self.CB_objects_lst.sort(key=lambda x: x.position)
        cb_lst = []
        for i, cb_i in enumerate(self.CB_objects_lst):
            # specify the boundaries
            if i == 0:
                # the leftmost crack
                cb_i.CB_model.Ll = cb_i.position
            else:
                # there is a crack at the left hand side
                cb_i.CB_model.Ll = (cb_i.position - cb_lst[-1].position)/2.
                cb_lst[-1].CB_model.Lr = cb_i.CB_model.Ll
            if i == len(self.CB_objects_lst) - 1:
                # the rightmost crack
                cb_i.CB_model.Lr = self.length - cb_i.position
            cb_lst.append(cb_i)
            
        for cb_i in cb_lst: 
            # specify the x range for cracks
            mask_right = self.x_arr >= (cb_i.position - cb_i.CB_model.Ll)
            mask_left = self.x_arr <= (cb_i.position + cb_i.CB_model.Lr)
            cb_i.x = self.x_arr[mask_left * mask_right] - cb_i.position
        self.CB_objects_lst = cb_lst           
        return cb_lst

    def get_current_cracking_state(self, load):
        '''Creates the list of CB objects that have been formed up to the current load
        '''
        idx_load_level = np.sum(np.array(self.cracking_stresses_lst) <= load)
        cb_lst = self.get_cbs_lst(self.crack_positions_lst[:idx_load_level],
                                      self.cracking_stresses_lst[:idx_load_level])
        return cb_lst
    
    def get_current_strnegth(self, load):
        if len(self.crack_positions_lst) is not 0:
            cb_lst = self.get_current_cracking_state(load)
            strengths = np.array([cb_i.max_sigma_c for cb_i in cb_lst])
            return np.min(strengths)
        else:
            return np.inf
    
    def sigma_m(self, W=None, w_lst=None):
        '''evaluates sigma_m along the composite
        the crack widths can be provided directly or computed given
        the sum of crack widths (eval_w_vect(W) call)'''
        if w_lst is None:
            w_lst = self.eval_w_vect(W)
        sigma_m = np.ones_like(self.matrix_strength)
        for i, cb_i in enumerate(self.current_CB_lst):
            crack_position_idx = np.argwhere(self.x_arr == cb_i.position)
            idx_l = crack_position_idx - len(np.nonzero(cb_i.x < 0.)[0])
            idx_r = crack_position_idx + len(np.nonzero(cb_i.x > 0.)[0]) + 1
            sigma_m[idx_l:idx_r] = cb_i.get_epsm_x(w_lst[i]) * self.CB_model.E_m
        return sigma_m
    
    def find_next_crack(self):
        '''Finds the crack openings given W = sum(w_i) and crack list.
        an iterative procedure of solving non-linear equations'''
        if len(self.crack_positions_lst) == 0:
            #first crack
            position = self.x_arr[np.argmin(self.matrix_strength)]
            sigmac = np.min(self.matrix_strength) / self.CB_model.E_m * self.CB_model.E_c
            epsc = sigmac / self.CB_model.E_c
            return epsc, sigmac, position, 0.002
        else:
            init_guess = self.cracking_strains_lst[-1] * self.length
            W = newton(self.residuum, init_guess)
            position = self.x_arr[np.argmin(self.matrix_strength - self.sigma_m(W))]
            w_lst = self.eval_w_vect(W)
            sigmac = self.CB_objects_lst[0].get_sigma_c(w_lst[0])
            epsc = (np.sum(np.array(w_lst)) + np.trapz(self.sigma_m(w_lst=w_lst),self.x_arr) / self.CB_model.E_m) / self.length
            return epsc, sigmac, position, W

    def residuum(self, W):
        '''Callback method for the identification of the
        next emerging crack calculated as the difference between
        the current matrix stress and strength. See the scipy newton call above.
        '''
        if W < 0.0:
            mu_sigmamu = np.mean(self.matrix_strength)
            residuum = mu_sigmamu - mu_sigmamu * W
        else: 
            residuum = np.min(self.matrix_strength - self.sigma_m(W))
        return residuum
    
    def eval_w_vect(self, W):
        if len(self.current_CB_lst) == 1:
            return np.array([W])
        else:
            for i in range(len(self.current_CB_lst)):
                self.current_CB_lst[i].damage_switch = False
            CB_arr = np.array(self.current_CB_lst)
            def scalar_func(f, w):
                return f.get_sigma_c(w)       
            vect_CB = np.vectorize(scalar_func)
            def min_func(w_arr):
                w_arr = np.abs(w_arr)
                if W - np.sum(w_arr) < 0.0:
                    return np.repeat(np.mean(w_arr),len(w_arr)) + 1.0
                w_arr_stacked = np.hstack((w_arr,W-np.sum(w_arr)))
                CB_sigmac_arr = vect_CB(CB_arr, w_arr_stacked)
                residuum = CB_sigmac_arr[:-1] - CB_sigmac_arr[1:]
                return residuum
            opt_result = root(min_func, np.repeat(W/float(len(self.current_CB_lst)),len(self.current_CB_lst)-1),
                              method='krylov', options={'maxiter':10})
            equil_w_arr = np.hstack((opt_result.x, W-np.sum(opt_result.x)))
            for i in range(len(self.current_CB_lst)):
                self.current_CB_lst[i].damage_switch = True
            return equil_w_arr
        
    
    def evaluate(self):
        # seek for the minimum strength redundancy to find the position
        # of the next crack
        while True:
            try:
                s = t.clock()
                epsc, sigmac, position, W = self.find_next_crack()
                print 'evaluation of the matrix crack #'+str(len(self.cracking_stresses_lst) + 1), t.clock() - s, 's'
                if len(self.cracking_strains_lst) > 1000:
                    plt.plot(self.x_arr, self.sigma_m(W) / self.CB_model.E_m, color='blue', lw=2)
                    plt.plot(self.x_arr, self.matrix_strength / self.CB_model.E_m, color='black', lw=2)
                    plt.show()
                self.crack_positions_lst.append(position)
                self.cracking_stresses_lst.append(sigmac)
                self.cracking_strains_lst.append(epsc)
                #append the last crack object
                self.CB_objects_lst.append(CB(position=self.crack_positions_lst[-1],
                                              cracking_stress=self.cracking_stresses_lst[-1],
                                              CB_model=CompositeCrackBridge(reinforcement_lst=self.CB_model.reinforcement_lst,
                                                                            E_m=self.CB_model.E_m,
                                                                            Gf=self.CB_model.Gf,
                                                                            ft=float(self.matrix_strength[np.argwhere(self.crack_positions_lst[-1] == self.x_arr)]) * 0.99
                                                                            )
                                              )
                                           )
            except:
                #plt.plot(np.hstack((0.0,self.cracking_strains_lst)), np.hstack((0.0,self.cracking_stresses_lst)), color='blue', lw=2)
                #plt.show()
                print 'composite saturated'
                break

    
if __name__ == '__main__':
    from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import ContinuousFibers, ShortFibers
    from stats.pdistrib.weibull_fibers_composite_distr import fibers_MC
    from matplotlib import pyplot as plt
    length = 250.
    nx = 5000
    random_field = RandomField(seed=True,
                               lacor=1.,
                               length=length,
                               nx=500,
                               nsim=1,
                               loc=.0,
                               shape=45.,
                               scale=3.6,
                               distr_type='Weibull'
                               )

    reinf1 = ContinuousFibers(r=3.5e-3,
                              tau=RV('gamma', loc=1e-6, shape=0.1040, scale=0.6554),
                              V_f=0.01,
                              E_f=182e3,
                              xi=fibers_MC(m=7.1, sV0=0.0069),
                              label='carbon',
                              n_int=500)

    CB_model = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf1],
                                 ft=1.0,
                                 Gf=0.15
                                 )

    scm = SCM(length=length,
              nx=nx,
              random_field=random_field,
              CB_model=CB_model
              )

    scm.evaluate()
