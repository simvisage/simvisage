'''
Created on Aug 23, 2012

@author: rch
'''

from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, WeakRef
    
import numpy as np

from mathkit.mfn import MFnLineArray

class CCLawBase(HasStrictTraits):
    '''Base class for concrete constitutive laws.'''
    # characteristic compressive stress [MPa]
    #
    f_ck = Float(60.0, input = True)
    eps_c_u = Float(0.0033, input = True)

    cs = WeakRef
    
    def __init__(self, *args, **kw):
        super(HasStrictTraits, self).__init__(*args, **kw)
        self.on_trait_change(self._notify_cs, '+input')

    def _notify_cs(self):
        if self.cs:
            self.cs.cc_modified = True

    mfn_vct = Property()
    def _get_mfn_vct(self):
        return np.vectorize(self.mfn.get_value)
            
class CCLawBlock(CCLawBase):
    '''Effective crack bridge Law with linear elastic response.'''
    
    #-----------------------------
    # 
    # for simplified constant stress-strain-diagram of the concrete (EC2)
    #-----------------------------
    mfn = Property(depends_on = '+input')
    @cached_property
    def _get_mfn(self):
        '''simplified constant stress-strain-diagram of the concrete (EC2)
        '''
        #(for standard concrete)
        f_ck = self.f_ck + 8. 
        if f_ck <= 50:
            lamda = 0.8
            eta = 1.0  
            eps_cu3 = self.eps_c_u
        # (for high strength concrete)
        #
        else:
            eta = 1.0 - (f_ck / 50.) / 200.
        # factor [-] to calculate the height of the compressive zone  
            lamda = 0.8 - (f_ck - 50.) / 400.
            eps_cu3 = (2.6 + 35. * ((90. - f_ck) / 100) ** 4.) / 1000. 
            
        xdata = np.hstack([0., (1. - lamda) * eps_cu3 - 0.00001, (1 - lamda) * eps_cu3, eps_cu3 ]) 
        ydata = np.hstack([0., 0., eta * (f_ck), eta * (f_ck), ])
       
        return MFnLineArray(xdata = xdata, ydata = ydata)       
    
class CCLawLinear(CCLawBase):
    '''Effective crack bridge Law based on fiber-bundle-model.'''
    #-----------------------------
    # for bilinear stress-strain-diagram of the concrete (EC2)
    #-----------------------------
    mfn = Property(depends_on = '+input')
    @cached_property
    def _get_mfn(self):
        '''bilinear stress-strain-diagram of the concrete
        '''
        #(for standard concrete)
        f_ck = self.f_ck + 8.  
        if f_ck <= 50.:
            eps_c3 = 0.00175
            eps_cu3 = self.eps_c_u
        #(for high strength concrete)
        else :
            eps_c3 = (1.75 + 0.55 * (f_ck - 50.) / 40.) / 1000.
            eps_cu3 = (2.6 + 35 * (90 - f_ck) ** 4.) / 1000.      
        # concrete law with limit for eps_c
       
        xdata = np.hstack([0., eps_c3, eps_cu3])
        ydata = np.hstack([0., (f_ck), (f_ck)])
        
        return MFnLineArray(xdata = xdata, ydata = ydata)
        
class CCLawQuadratic(CCLawBase):
    '''Effective crack bridge Law using a cubic polynomial.'''
    
    mfn = Property(depends_on = '+input')
    @cached_property
    def _get_mfn(self):
        '''quadratic stress-strain-diagram of the concrete
        '''
        # (for all concretes up to f_cm=88 N/mm2) #max epislon_c1u
        f_cm = self.f_ck + 8 
        E_tan = 22. * (f_cm / 10) ** 0.3 * 1000.
        eps_c1 = min(0.7 * f_cm ** 0.31, 2.8) / 1000. #EC2
        # @todo: with constant value this yields negative values for strains close to 'eps_c1u'
#        eps_c1 = 0.0022 #Brockmann
        E_sec = f_cm / eps_c1   
       
        if self.f_ck <= 50.:
            eps_c1u = self.eps_c_u
            eps_arr = np.linspace(0., eps_c1u, num = 100.)
        
        elif self.f_ck > 50.:   
            eps_c1u = (2.8 + 27. * (((98. - f_cm) / 100.) ** 4.)) / 1000.
            eps_arr = np.linspace(0., eps_c1u, num = 100.) 

        k = E_tan / E_sec
        sig_c_arr = ((k * eps_arr / eps_c1 - (eps_arr / eps_c1) ** 2.) / (1. + (k - 2.) * eps_arr / eps_c1)) * f_cm
        
        xdata = eps_arr  
        ydata = sig_c_arr
        
        return MFnLineArray(xdata = xdata, ydata = ydata)
    
    
if __name__ == '__main__':
    cc_law = CCLawQuadratic()
    
    print cc_law.mfn.get_value(-0.00000001)
    import pylab as p
    cc_law.mfn.plot(p)
    p.show()
