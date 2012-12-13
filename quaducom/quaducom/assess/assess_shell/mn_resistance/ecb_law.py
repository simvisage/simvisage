'''
Created on Aug 23, 2012

@author: rch
'''

from etsproxy.traits.api import \
    Float, Property, \
    cached_property, Str, Array, \
    Int, List

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, ModelView, VGroup, HGroup, RangeEditor, InstanceEditor

from constitutive_law import CLBase

import numpy as np

from math import exp, log

from mathkit.mfn import MFnLineArray

class ECBLBase(CLBase):
    '''Base class for Effective Crack Bridge Laws.'''

    u0 = List([0.0, 0.0])

class ECBLLinear(ECBLBase):
    '''Effective crack bridge Law with linear elastic response.'''

    eps_tex_u = Float(0.01, enter_set = True, auto_set = False, input = True)
    E_tex = Float(80000, enter_set = True, auto_set = False, input = True)
    u0 = List([ 0.01, 80000. ], enter_set = True, auto_set = False)

    sig_tex_u = Property(depends_on = '+input')
    @cached_property
    def _get_sig_tex_u(self):
        return self.E_tex * self.eps_tex_u

    cnames = ['eps_tex_u', 'E_tex']

    eps_arr = Property(depends_on = '+input')
    @cached_property
    def _get_eps_arr(self):
        return np.array([ 0., self.eps_tex_u])

    sig_arr = Property(depends_on = '+input')
    @cached_property
    def _get_sig_arr(self):
        # with limit for eps_tex
        #
        return self.E_tex * self.eps_arr

class ECBLFBM(ECBLBase):
    '''Effective crack bridge Law based on fiber-bundle-model.'''

    sig_tex_u = Float(1216, input = True)
    eps_tex_u = Float(0.014, input = True)
    m = Float(0.5, input = True)

    cnames = ['eps_tex_u', 'm']

    u0 = List([0.01266923, 0.5 ])

    eps_arr = Property(depends_on = '+input')
    @cached_property
    def _get_eps_arr(self):
        return np.linspace(0, self.eps_tex_u, num = 100.)

    sig_arr = Property(depends_on = '+input')
    @cached_property
    def _get_sig_arr(self):
        sig_tex_u = self.sig_tex_u
        eps_tex_u = self.eps_tex_u
        eps_arr = self.eps_arr
        m = self.m
        return (sig_tex_u / eps_tex_u /
                exp(-pow(exp(-log(m) / self.m), 1.0 * m)) *
                eps_arr * np.exp(-np.power(eps_arr / eps_tex_u * exp(-log(m) / m), 1.0 * m)))

class ECBLCubic(ECBLBase):
    '''Effective crack bridge Law using a cubic polynomial.'''

    sig_tex_u = Float(1250, input = True)
    eps_tex_u = Float(0.016, input = True)
    var_a = Float(-5e+6, input = True)

    cnames = ['eps_tex_u', 'var_a']

    u0 = List([ 0.016, -5000000. ])

    eps_arr = Property(depends_on = '+input')
    @cached_property
    def _get_eps_arr(self):
        return np.linspace(0, self.eps_tex_u, num = 100.)

    sig_arr = Property(depends_on = '+input')
    @cached_property
    def _get_sig_arr(self):
        # for horizontal tangent at eps_tex_u
        sig_tex_u, var_a, eps_tex_u = self.sig_tex_u, self.var_a, self.eps_tex_u
        eps_arr = self.eps_arr
        var_b = -(sig_tex_u + 2. * var_a * eps_tex_u ** 3.) / eps_tex_u ** 2.
        var_c = -3. * var_a * eps_tex_u ** 2. - 2. * var_b * eps_tex_u
        sig_arr = var_a * eps_arr ** 3. + var_b * eps_arr ** 2. + var_c * eps_arr
        return sig_arr

class ECBLBilinear(ECBLBase):
    '''Effective crack bridge Law using a cubic polynomial.'''

    sig_tex_u = Float(1250, input = True)
    eps_tex_u = Float(0.014, input = True)
    var_a = Float(50000, input = True)
    eps_el_fraction = Float(0.0001, input = True)

    cnames = ['eps_tex_u', 'var_a']

    u0 = List([ 0.014, 50000. ])

    eps_arr = Property(depends_on = '+input')
    @cached_property
    def _get_tex_arr(self):
        return np.hstack([0., self.eps_el_fraction * self.eps_tex_u, self.eps_tex_u ])

    sig_arr = Property(depends_on = '+input')
    @cached_property
    def _get_sig_arr(self):
        return np.hstack([0., self.var_a * self.sig_tex_u, self.sig_tex_u])

class ECBLPiecewiseLinearOld(ECBLBase):
    '''Effective crack bridge Law using a piecewise linear function.'''

    sig_level_arr = Array(float, input = True)
    def _sig_level_arr_default(self):
        E_eff = self.sig_u_trial / self.eps_tex_u
        return E_eff * self.eps_arr

    sig_tex_u = Float(1000.0, input = True)
    eps_tex_u = Float(0.01, input = True)
    n_eps = Int(2, input = True)

    eps_fraction_arr = Property(depends_on = 'n_eps')
    @cached_property
    def _get_eps_fraction_arr(self):
        return np.linspace(0.0, 1.0, self.n_eps)

    cnames = ['eps_tex_u', 'sig_level_arr']

    sig_u_trial = Float(1000)

    u0 = Property(depends_on = 'eps_tex_u, sig_u_trial, eps_fraction_list')
    @cached_property
    def _get_u0(self):
        return self.sig_level_arr[1:]

    eps_arr = Property(depends_on = '+input')
    @cached_property
    def _get_eps_arr(self):
        return self.eps_fraction_arr * self.eps_tex_u

    sig_arr = Property(depends_on = '+input')
    @cached_property
    def _get_sig_arr(self):
        return self.sig_level_arr

class ECBLPiecewiseLinear(ECBLBase):
    '''Effective crack bridge Law using a piecewise linear function.'''

    sig_level_arr = Array(float, input = True)
    def _sig_level_arr_default(self):
        E_eff = self.sig_tex_u / self.eps_tex_u
        return E_eff * self.eps_arr

    sig_tex_u = Float(800.0, input = True)
    eps_tex_u = Float(0.01, input = True)
    n_eps = Int(2, input = True)

    eps_fraction_arr = Property(depends_on = 'n_eps')
    @cached_property
    def _get_eps_fraction_arr(self):
        return np.linspace(0.0, 1.0, self.n_eps)

    cnames = ['eps_tex_u', 'sig_level_arr']

    u0 = Property(depends_on = 'eps_tex_u, sig_tex_u, eps_fraction_list')
    @cached_property
    def _get_u0(self):
        return self.sig_level_arr[1:]

    eps_arr = Array(float)
    def _eps_arr_default(self):
        return self.eps_fraction_arr * self.eps_tex_u

    sig_arr = Array(float)
    def _sig_arr_default(self):
        return self.sig_level_arr

    def set_sig_eps_arr(self, eps_arr, sig_arr):
        self.eps_arr = eps_arr
        self.sig_arr = sig_arr

if __name__ == '__main__':
    from constitutive_law import ConstitutiveLawModelView
    #ecbl = ECBLFBM()
    ecbl = ECBLPiecewiseLinear()
    ew = ConstitutiveLawModelView(model = ecbl)
    ew.configure_traits()
