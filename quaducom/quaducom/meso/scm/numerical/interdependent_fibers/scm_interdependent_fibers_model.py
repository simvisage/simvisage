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
from scipy.interpolate import interp1d
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
from scipy.optimize import brentq, newton
import copy
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import ContinuousFibers
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers
from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx import CompositeCrackBridge
from spirrid.rv import RV
from matplotlib import pyplot as plt
from math import pi
import time as t


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

    def get_w_sigma_c(self, load):
        '''
        evaluates the crack width at a particular load
        '''
        if load > self.max_sigma_c:
            print 'Warning: applied load', load, 'MPa higher then strength ', self.max_sigma_c, 'MPa'
            return np.nan * np.ones_like(self.x)
        else:
            return self.representative_cb.interpolate_w(self.Ll, self.Lr, load)
            
    def get_epsm_x_w(self, load):
        '''
        evaluates the matrix strain profile
        '''
        if load > self.max_sigma_c:
            print 'Warning: applied load', load, 'MPa higher then strength ', self.max_sigma_c, 'MPa'
            return np.nan * np.ones_like(self.x)
        else:
            return self.representative_cb.interpolate_epsm(self.Ll, self.Lr, load, self.x)
            

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
    n_BC_CB = Int(10, desc='# of discretization points for boundary conditions')
    n_x_CB = Int(200, desc='# of discretization points for the long. axis of a single CB')

    representative_cb = Instance(RepresentativeCB)
    def _representative_cb_default(self):
        return RepresentativeCB(CB_model=self.CB_model,
                            load_sigma_c_arr=self.load_sigma_c_arr,
                            length=self.length, n_w=self.n_w_CB, n_BC=self.n_BC_CB, n_x=self.n_x_CB
                            )

    sigma_c_crack = List
    # @todo: rr comment - what's the difference between sigma_c_crack and cracks_list 
    cracks_list = List

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

    def sort_cbs(self):
        # sorts the CBs by position and adjusts the boundary conditions
        # sort the CBs
        cb_list = self.cracks_list[-1]
        crack_position = cb_list[-1].position
        cb_list = sorted(cb_list, key=attrgetter('position'))
        # find idx of the new crack
        for i, crack in enumerate(cb_list):
            if crack.position == crack_position:
                idx = i
        # specify the boundaries
        if idx != 0:
            # there is a crack at the left hand side
            cbl = cb_list[idx - 1]
            cb = cb_list[idx]
            cbl.Lr = (cb.position - cbl.position) / 2.
            cb.Ll = cbl.Lr
        else:
            # the new crack is the first from the left hand side
            cb_list[idx].Ll = cb_list[idx].position

        if idx != len(cb_list) - 1:
            # there is a crack at the right hand side
            cb, cbr = cb_list[idx], cb_list[idx + 1]
            cbr.Ll = (cbr.position - cb.position) / 2.
            cb.Lr = cbr.Ll
        else:
            # the new crack is the first from the right hand side
            cb_list[idx].Lr = self.length - cb_list[idx].position

        # specify the x range and stress profile for
        # the new crack and its neighbors
        idxs = [idx - 1, idx, idx + 1]
        if idx == 0:
            idxs.remove(-1)
        if idx == len(cb_list) - 1:
            idxs.remove(len(cb_list))
        for idx in idxs:
            mask1 = self.x_arr >= (cb_list[idx].position - cb_list[idx].Ll)
            if idx == 0:
                mask1[0] = True
            mask2 = self.x_arr <= (cb_list[idx].position + cb_list[idx].Lr)
            cb_list[idx].x = self.x_arr[mask1 * mask2] - cb_list[idx].position
        self.cracks_list[-1] = cb_list

    def cb_list(self, load):
        '''Get the number of cracks initiated below the supplied load.
        @todo: rename cracks_list with cracking_state
        '''
        if len(self.cracks_list) is not 0:
            idx = np.sum(np.array(self.sigma_c_crack) < load) - 1
            if idx == -1:
                return [None]
            else:
                return self.cracks_list[idx]
        else:
            return [None]

    def sigma_m(self, load):
        Em = self.CB_model.E_m
        Ec = self.CB_model.E_c
        sigma_m = load * Em / Ec * np.ones(len(self.x_arr))
        cb_load = self.cb_list(load)
        if cb_load[0] is not None:
            for cb in cb_load:
                crack_position_idx = np.argwhere(self.x_arr == cb.position)
                left = crack_position_idx - len(np.nonzero(cb.x < 0.)[0])
                right = crack_position_idx + len(np.nonzero(cb.x > 0.)[0]) + 1
                sigma_m[left:right] = cb.get_epsm_x_w(float(load)) * Em
        return sigma_m

    def residuum(self, q):
        '''Callback method for the identification of the
        next emerging crack calculated as the difference between
        the current matrix stress and strength. See the np.brentq call below.
        '''
        residuum = np.min(self.matrix_strength - self.sigma_m(q))
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
                print 'evaluation of the matrix crack #'+str(len(self.sigma_c_crack)), t.clock() - s, 's'
            except:
                print 'composite saturated'
                break
            crack_position = self.x_arr[np.argmin(self.matrix_strength - 
                                                  self.sigma_m(sigc_min))]
            new_cb = CB(position=float(crack_position),
                     crack_load_sigma_c=sigc_min - self.load_sigma_c_arr[-1] / 1000.,
                     representative_cb=self.representative_cb)
            self.sigma_c_crack.append(sigc_min - 1e-10)
            if len(self.cracks_list) is not 0:
                self.cracks_list.append(copy.copy(self.cracks_list[-1])
                                        + [new_cb])
            else:
                self.cracks_list.append([new_cb])
            
            self.sort_cbs()
            #plt.plot(self.x_arr, self.sigma_m(sigc_min) / self.CB_model.E_m, color='blue', lw=2)
            #plt.plot(self.x_arr, self.matrix_strength / self.CB_model.E_m, color='black', lw=2)
            #plt.show()
            if float(crack_position) == last_pos:
                print last_pos
                raise ValueError('''got stuck in loop,
                try to adapt x, w, BC ranges''')
            last_pos = float(crack_position)
                

if __name__ == '__main__':
    from quaducom.meso.homogenized_crack_bridge.elastic_matrix.reinforcement import ContinuousFibers
    from stats.pdistrib.weibull_fibers_composite_distr import fibers_MC
    length = 500.
    nx = 2000
    random_field = RandomField(seed=False,
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
                              tau=RV('weibull_min', loc=0.01, scale=.02, shape=2.),
                              V_f=0.03,
                              E_f=200e3,
                              xi=fibers_MC(m=7., sV0=0.01),
                              label='carbon',
                              n_int=500)

    reinf2 = ContinuousFibers(r=3.5e-3,
                              tau=RV('weibull_min', loc=0.01, scale=.1, shape=2.),
                              V_f=0.005,
                              E_f=200e3,
                              xi=fibers_MC(m=7., sV0=0.005),
                              label='carbon',
                              n_int=500)

    CB_model = CompositeCrackBridge(E_m=25e3,
                                 reinforcement_lst=[reinf1],
                                 )

    scm = SCM(length=length,
              nx=nx,
              random_field=random_field,
              CB_model=CB_model,
              load_sigma_c_arr=np.linspace(0.01, 20., 100),
              )

    scm.evaluate()
