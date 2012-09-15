'''
Created on Aug 23, 2012

@author: rch
'''

from etsproxy.traits.api import \
    HasStrictTraits, List, Float, Property, cached_property, Str
    
import numpy as np

from math import exp, log

from mathkit.mfn import MFnLineArray
    
class ECBLBase(HasStrictTraits):
    '''Base class for Effective Crack Bridge Laws.'''
    
    u0 = List([0.0, 0.0])
    
    # names of calibrated parameters
    cnames = List(Str)
    
    def set_cparams(self, *args):
        for name, value in zip(self.cnames, args):
            setattr(self, name, value)

    def __call__(self, *args):
        self.set_cparams(*args)
        return self.eps_tex_arr, self.sig_tex_arr

    arr = Property()
    def _get_arr(self):
        return self.eps_tex_arr, self.sig_tex_arr
    
    mfn = Property()
    def _get_mfn(self):
        return MFnLineArray(xdata = self.eps_tex_arr,
                            ydata = self.sig_tex_arr)

    mfn_vct = Property()
    def _get_mfn_vct(self):
        return np.vectorize(self.mfn.get_value)

class ECBLLinear(ECBLBase):
    '''Effective crack bridge Law with linear elastic response.'''
    
    eps_tex_u = Float(0.01, input = True)
    E_tex = Float(80000, input = True)
    sig_tex_u = Float
    u0 = List([ 0.01, 80000. ])

    cnames = ['eps_tex_u', 'E_tex']

    eps_tex_arr = Property(depends_on = '+input')
    @cached_property
    def _get_eps_tex_arr(self):
        return np.array([ 0., self.eps_tex_u])
        
    sig_tex_arr = Property(depends_on = '+input')
    @cached_property
    def _get_sig_tex_arr(self):
        # with limit for eps_tex
        #
        return self.E_tex * self.eps_tex_arr

class ECBLFBM(ECBLBase):
    '''Effective crack bridge Law based on fiber-bundle-model.'''
    
    sig_tex_u = Float(1250, input = True)
    eps_tex_u = Float(0.014, input = True)
    m = Float(0.5, input = True)

    cnames = ['eps_tex_u', 'm']
        
    u0 = List([0.014, 0.5 ])
        
    eps_tex_arr = Property(depends_on = '+input')
    @cached_property
    def _get_eps_tex_arr(self):
        return np.linspace(0, self.eps_tex_u, num = 100.)
    
    sig_tex_arr = Property(depends_on = '+input')
    @cached_property
    def _get_sig_tex_arr(self):
        sig_tex_u = self.sig_tex_u
        eps_tex_u = self.eps_tex_u
        eps_tex_arr = self.eps_tex_arr
        m = self.m
        return (sig_tex_u / eps_tex_u / 
                exp(-pow(exp(-log(m) / self.m), 1.0 * m)) * 
                eps_tex_arr * np.exp(-np.power(eps_tex_arr / eps_tex_u * exp(-log(m) / m), 1.0 * m)))            
            
class ECBLCubic(ECBLBase):
    '''Effective crack bridge Law using a cubic polynomial.'''
    
    sig_tex_u = Float(1250, input = True)
    eps_tex_u = Float(0.016, input = True)
    var_a = Float(-5e+6, input = True)

    cnames = ['eps_tex_u', 'var_a']

    u0 = List([ 0.016, -5000000. ])

    eps_tex_arr = Property(depends_on = '+input')
    @cached_property
    def _get_eps_tex_arr(self):
        return np.linspace(0, self.eps_tex_u, num = 100.)
                                       
    sig_tex_arr = Property(depends_on = '+input')
    @cached_property
    def _get_sig_tex_arr(self):
        # for horizontal tangent at eps_tex_u
        sig_tex_u, var_a, eps_tex_u = self.sig_tex_u, self.var_a, self.eps_tex_u
        eps_tex_arr = self.eps_tex_arr
        var_b = -(sig_tex_u + 2. * var_a * eps_tex_u ** 3.) / eps_tex_u ** 2. 
        var_c = -3. * var_a * eps_tex_u ** 2. - 2. * var_b * eps_tex_u 
        sig_tex_arr = var_a * eps_tex_arr ** 3. + var_b * eps_tex_arr ** 2. + var_c * eps_tex_arr
        return sig_tex_arr         

class ECBLBilinear(ECBLBase):
    '''Effective crack bridge Law using a cubic polynomial.'''
    
    sig_tex_u = Float(1250, input = True)
    eps_tex_u = Float(0.014, input = True)
    var_a = Float(50000, input = True)
    eps_el_fraction = Float(0.0001, input = True)

    cnames = ['eps_tex_u', 'var_a']

    u0 = List([ 0.014, 50000. ])
    
    eps_tex_arr = Property(depends_on = '+input')
    @cached_property
    def _get_eps_tex_arr(self):
        return np.hstack([0., self.eps_el_fraction * self.eps_tex_u, self.eps_tex_u ])

    sig_tex_arr = Property(depends_on = '+input')
    @cached_property
    def _get_sig_tex_arr(self):
        return np.hstack([0., self.var_a * self.sig_tex_u, self.sig_tex_u])

if __name__ == '__main__':
    ecbl = ECBLFBM(sig_tex_u = 2e+3, m = 2)
    import pylab as p
    print ecbl(0.04, 2)[0]
    print ecbl(0.04, 2)[1]
    
    p.plot(*ecbl(0.04, 2))
    p.show()
    
