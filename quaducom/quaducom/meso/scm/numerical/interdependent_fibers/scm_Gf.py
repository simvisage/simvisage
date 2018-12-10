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
    pre_cracking_W = Float
    
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

    # crack objects
    CB_objects_lst = List
    # current prescribed sum of crack openings
    W = Float
    
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
        # realization and creates a spline representation
        rf = self.random_field.random_field
        rf_spline = MFnLineArray(xdata=self.random_field.xgrid, ydata=rf)
        return rf_spline.get_values(self.x_arr)

    sorted_CB_lst = Property(List, depends_on='CB_objects_lst')
    @cached_property
    def _get_sorted_CB_lst(self):
        ''' sorts the CBs by position and adjusts the boundary conditions '''
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
        return cb_lst
    
    w_lst = Property(List, depends_on='W,CB_objects_lst')
    @cached_property
    def _get_w_lst(self):
        '''evaluates the crack openings'''
        cb_lst = self.sorted_CB_lst
        for cb_i in cb_lst:
            cb_i.damage_switch = False
        CB_arr = np.array(cb_lst)
        def scalar_sigma_c(w, f):
            return f.get_sigma_c(w)      
        vect_sigma_c = np.vectorize(scalar_sigma_c)
        
        res = []
        def min_func(w_arr):
            w_arr = np.abs(w_arr)
            CB_sigmac_arr = vect_sigma_c(w_arr, CB_arr)
            residuum = np.hstack((CB_sigmac_arr[:-1] - CB_sigmac_arr[1:], self.W-np.sum(w_arr)))
            res.append(CB_sigmac_arr)
            return residuum
        opt_result = root(min_func, np.repeat(self.W/float(len(cb_lst)),len(cb_lst)),
                          method='krylov', options={'maxiter':200})
#         res = np.array(res)
#         if len(CB_arr) > 2:
#             for i in np.arange(res.shape[1] - 1):
#                 plt.plot(res[:,i])
#             plt.show()
        for cb_i in cb_lst:
            cb_i.damage_switch = True 
        return np.abs(opt_result.x)
    
    eps_m = Property(Array, depends_on='W,CB_objects_lst')
    @cached_property
    def _get_eps_m(self):
        '''evaluates sigma_m along the composite
        the crack widths can be provided directly or computed given
        the sum of crack widths (eval_w_vect(W) call)'''
        eps_m = np.ones_like(self.matrix_strength)
        for i, cb_i in enumerate(self.sorted_CB_lst):
            crack_position_idx = np.argwhere(self.x_arr == cb_i.position)
            idx_l = crack_position_idx - len(np.nonzero(cb_i.x < 0.)[0])
            idx_r = crack_position_idx + len(np.nonzero(cb_i.x > 0.)[0]) + 1
            eps_m[idx_l:idx_r] = cb_i.get_epsm_x(self.w_lst[i])
        return eps_m
    
    
    eps_f = Property(Array, depends_on='W,CB_objects_lst')
    @cached_property
    def _get_eps_f(self):
        '''evaluates sigma_m along the composite
        the crack widths can be provided directly or computed given
        the sum of crack widths (eval_w_vect(W) call)'''
        eps_f = np.ones_like(self.matrix_strength)
        for i, cb_i in enumerate(self.sorted_CB_lst):
            crack_position_idx = np.argwhere(self.x_arr == cb_i.position)
            idx_l = crack_position_idx - len(np.nonzero(cb_i.x < 0.)[0])
            idx_r = crack_position_idx + len(np.nonzero(cb_i.x > 0.)[0]) + 1
            eps_f[idx_l:idx_r] = cb_i.get_epsf_x(self.w_lst[i])
        return eps_f
    
    def find_next_crack(self):
        '''Finds the crack openings given W = sum(w_i) and crack list.
        an iterative procedure of solving non-linear equations'''
        def residuum(W):
            '''Callback method for the identification of the
            next emerging crack calculated as the difference between
            the current matrix stress and strength. See the scipy newton call above.
            '''
            self.W = W
            if self.W <= 0.0 or self.W > 0.3 * len(self.CB_objects_lst):
                min_strength = np.min(self.matrix_strength)
                residuum = min_strength - self.CB_model.E_c * self.W / self.length
            else:
                residuum = np.min(self.matrix_strength - self.eps_m * self.CB_model.E_m)
            return residuum

        W_pre_crack = newton(residuum, self.CB_objects_lst[-1].pre_cracking_W)
        position_new_crack = self.x_arr[np.argmin(self.matrix_strength - self.eps_m * self.CB_model.E_m)]         
        sigmac_pre_crack = self.sorted_CB_lst[0].get_sigma_c(self.w_lst[0])
        if len(self.CB_objects_lst) > 2:
            plt.plot(self.x_arr, self.eps_f, color='red', lw=2)
            plt.plot(self.x_arr, self.eps_m, color='blue', lw=2)
            plt.plot(self.x_arr, self.matrix_strength / self.CB_model.E_m, color='black', lw=2)
            plt.title('lm')
            plt.show()
#         for i, w in enumerate(self.w_lst):
#             print self.sorted_CB_lst[i].get_sigma_c(w)
        #epsc = (np.sum(np.array(w_lst)) + np.trapz(self.sigma_m(w_lst=w_lst),self.x_arr) / self.CB_model.E_m) / self.length
        return sigmac_pre_crack, position_new_crack, W_pre_crack, self.w_lst
    
    def evaluate(self):
        '''evaluates the cracking history and saves it in a list'''
        if len(self.CB_objects_lst) == 0:
            # first crack
            position_first_crack = self.x_arr[np.argmin(self.matrix_strength)]
            sigmac_pre_crack = np.min(self.matrix_strength) / self.CB_model.E_m * self.CB_model.E_c
            pre_cracking_W = sigmac_pre_crack  / self.CB_model.E_c * self.length
            first_CB = CB(position=position_first_crack,
                          cracking_stress=sigmac_pre_crack,
                          pre_cracking_W=pre_cracking_W,
                          CB_model=CompositeCrackBridge(reinforcement_lst=self.CB_model.reinforcement_lst,
                                                        E_m=self.CB_model.E_m,
                                                        Gf=self.CB_model.Gf,
                                                        ft=float(self.matrix_strength[np.argwhere(position_first_crack == self.x_arr)]) * 0.99
                                                        )
                          )
            self.CB_objects_lst.append(first_CB)
            
        while True:
            try:
                s = t.clock()
                sigmac_pre_crack, position_new_crack, pre_cracking_W, pre_cracking_w_lst = self.find_next_crack()
                cb_lst = []
                for i, cb_i in enumerate(self.sorted_CB_lst):
                    if cb_i.CB_model.w_unld < pre_cracking_w_lst[i]:
                        cb_i.CB_model.w_unld = pre_cracking_w_lst[i]
                    cb_lst.append(cb_i)
                self.CB_objects_lst = cb_lst
    
                new_CB = CB(position=position_new_crack,
                            cracking_stress=sigmac_pre_crack,
                            pre_cracking_W=pre_cracking_W,
                            CB_model=CompositeCrackBridge(reinforcement_lst=self.CB_model.reinforcement_lst,
                                                          E_m=self.CB_model.E_m,
                                                          Gf=self.CB_model.Gf,
                                                          ft=float(self.matrix_strength[np.argwhere(position_new_crack == self.x_arr)]) * 0.99
                                                        )
                            )
                
                self.CB_objects_lst.append(new_CB)
                print('crack #'+str(len(self.CB_objects_lst) + 1), 'evaluated, time: ', t.clock() - s, 's, Gf = ', self.CB_model.Gf)

            except:
#                 #plt.plot(np.hstack((0.0,self.cracking_W_lst)), np.hstack((0.0,self.cracking_stresses_lst)), color='blue', lw=2)
#                 #plt.show()
                print('composite saturated')
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
                               shape=12.,
                               scale=5.2,
                               distr_type='Weibull'
                               )

    reinf1 = ContinuousFibers(r=3.5e-3,
                              tau=RV('gamma', loc=1e-6, shape=0.1040, scale=0.6554),
                              V_f=0.01,
                              E_f=182e3,
                              xi=fibers_MC(m=7.1, sV0=0.0069),
                              label='carbon',
                              n_int=500)

    cracks = []
    Gf=[0.2, 0.05, 0.1, 0.15, 0.2]
    for Gfi in Gf:
        CB_model = CompositeCrackBridge(E_m=25e3,
                                     reinforcement_lst=[reinf1],
                                     ft=1.0,
                                     Gf=Gfi
                                     )
    
        scm = SCM(length=length,
                  nx=nx,
                  random_field=random_field,
                  CB_model=CB_model
                  )
        scm.evaluate()
        cracks.append(len(scm.CB_objects_lst))
    print(Gf)
    print(cracks)
