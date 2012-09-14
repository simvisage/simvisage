'''
Created on Sep 4, 2012

@author: rch
'''
from etsproxy.traits.api import \
    HasTraits, Float, Property, cached_property, Int, Array, \
    Trait, Event, on_trait_change, Instance
    
from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit
    
from ecb_law import \
    ECBLLinear, ECBLFBM, ECBLCubic, ECBLBilinear

from cc_law import \
    CCLawBlock, CCLawLinear, CCLawQuadratic

import numpy as np

import math

class ECBCrossSection(HasTraits):

    #---------------------------------------------------------------
    # Cross section characteristics needed for tensile specimens 
    #---------------------------------------------------------------
    
    # thickness of reinforced cross section
    #
    thickness = Float(0.06, geo_input = True)

    # total number of reinforcement layers [-]
    # 
    n_layers = Int(12, geo_input = True)

    #---------------------------------------------------------------
    # Cross section characteristics needed for bending specimens 
    #---------------------------------------------------------------

    # width of the cross section [m]
    #
    width = Float(0.20, geo_input = True)
    
    # number of rovings in 0-direction of one composite 
    # layer of the bending test [-]:
    #
    n_rovings = Int(23, geo_input = True)
    
    # cross section of one roving [mm**2]:
    #
    A_roving = Float(0.461, geo_input = True)

    #===========================================================================
    # Distribution of reinforcement
    #===========================================================================

    # spacing between the layers [m]
    #
    s_tex_z = Property(depends_on = '+geo_input')
    @cached_property
    def _get_s_tex_z(self):
        return self.thickness / (self.n_layers + 1)

    # distance from the top of each reinforcement layer [m]:
    #
    z_ti_arr = Property(depends_on = '+geo_input')
    @cached_property
    def _get_z_ti_arr(self):
        return np.array([ self.thickness - (i + 1) * self.s_tex_z for i in range(self.n_layers) ],
                      dtype = float)

    # distance of reinforcement layers from the bottom 
    #
    zz_ti_arr = Property
    def _get_zz_ti_arr(self):
        return self.thickness - self.z_ti_arr
    
    # number of subdivisions of the compressive zone
    #
    n_cj = Float(20., sig_c_modified = True)

    zeta_cj_arr = Property(Array)
    @cached_property
    def _get_zeta_cj_arr(self):
        '''subdivide the compression zone 'x' in 'n_cj' sub-areas;
        'zeta_cj_arr' giving the fraction of each distance of the sub-area from the upper surface 
        with respect to the compressive zone 'x'
        '''
        return np.arange(self.n_cj) / self.n_cj + 1. / (2. * self.n_cj)

    #===========================================================================
    # Strain state
    #===========================================================================

    eps_up = Float(0.0033, eps_input = True)
    eps_lo = Float(-0.0172, eps_input = True)
        
    #===========================================================================
    # Discretization conform to the tex layers
    #===========================================================================

    eps_i_arr = Property(depends_on = '+eps_input,+geo_input')
    @cached_property
    def _get_eps_i_arr(self):
        '''CALIBRATION: derive the unknown constitutive law of the layer
        (effective crack bridge law)
        '''
        # ------------------------------------------------------------------------                
        # geometric params independent from the value for 'eps_t'
        # ------------------------------------------------------------------------                
        thickness = self.thickness
        z_ti_arr = self.z_ti_arr
        
        # ------------------------------------------------------------------------                
        # derived params depending on value for 'eps_t'
        # ------------------------------------------------------------------------                

        # heights of the compressive zone:
        #
        x = abs(self.eps_up) / (abs(self.eps_lo) + abs(self.eps_lo)) * thickness

        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = math.fabs(self.eps_lo) / (thickness - x) * (z_ti_arr - x)

        # use a ramp function to consider only negative strains
        #eps_ci_arr = -eps_i_arr[ np.where(eps_i_arr <= 0.0)]
        eps_ci_arr = np.fabs(-np.fabs(eps_i_arr) + eps_i_arr) / 2.0

        # use a ramp function to consider only positive strains
        #eps_ti_arr = eps_i_arr[ np.where(eps_i_arr <= 0.0)]
        eps_ti_arr = (np.fabs(eps_i_arr) + eps_i_arr) / 2.0     
        
        return x, eps_ti_arr, eps_ci_arr

    x = Property
    def _get_x(self):
        return self.eps_i_arr[0]
    
    eps_ti_arr = Property
    def _get_eps_ti_arr(self):
        return self.eps_i_arr[1]
    
    eps_ci_arr = Property
    def _get_eps_ci_arr(self):
        return self.eps_i_arr[2]

    eps_cj_arr = Property(depends_on = '+eps_input,+geo_input')    
    @cached_property
    def _get_eps_cj_arr(self):
        '''get compressive strain at each integration layer of the compressive zone [-]:
        for 'stress_case' flexion
        '''
        # for calibration us measured compressive strain
        # @todo: use mapped traits instead
        #
        eps_cj_arr = (1. - self.zeta_cj_arr) * abs(self.eps_up)
        return eps_cj_arr

    z_cj_arr = Property(depends_on = '+eps_input,+geo_input')
    @cached_property
    def _get_z_cj_arr(self):
        '''Get the discretizaton of the  compressive zone
        '''
        return self.x * self.zeta_cj_arr


    #===========================================================================
    # Compressive concrete constitutive law
    #===========================================================================
    cc_law_type = Trait('constant', dict(constant = CCLawBlock,
                                         linear = CCLawLinear,
                                         quadratic = CCLawQuadratic),
                      config_modified = True)
    
    cc_law = Property(depends_on = 'cc_law_type')
    @cached_property
    def _get_cc_law(self):
        '''Construct the compressive concrete law'''
        return self.cc_law_type_(f_ck = self.f_ck)

    sig_c_mfn = Property(depends_on = '+sig_c_modified')
    def _get_sig_c_mfn(self):
        return self.cc_law.mfn

    sig_c_mfn_vct = Property(depends_on = '+sig_c_modified')
    @cached_property
    def _get_sig_c_mfn_vct(self):
        return np.vectorize(self.sig_c_mfn.get_value)

    #===========================================================================
    # Effective crack bridge law
    #===========================================================================
    ecbl_type = Trait('fbm', dict(fbm = ECBLFBM,
                                  cubic = ECBLCubic,
                                  linear = ECBLLinear,
                                  bilinear = ECBLBilinear),
                      config_modified = True)
    
    ecbl = Property(depends_on = 'ecbl_type')
    @cached_property
    def _get_ecbl(self):
        return self.ecbl_type_(sig_tex_u = self.sig_tex_u)

    ecbl_modified = Event
    @on_trait_change('ecbl.+input')
    def _set_ecbl_modified(self):
        self.ecbl_modified = True
    
    ecbl_tex_mfn = Property(depends_on = '+config_modified, ecbl_modified')
    @cached_property
    def _get_ecbl_tex_arr(self):
        '''Get the arrays sigma - epsilon defining the crack bridge law.
        '''
        return self.ecbl.arr
    
    ecbl_tex_mfn = Property(depends_on = '+config_modified, ecbl_modified')
    @cached_property
    def _get_ecbl_tex_mfn(self):
        '''Get the callable function for effective crack brige law.
        '''
        return self.ecbl.mfn

    ecbl_tex_mfn_vct = Property(depends_on = '+config_modified, ecbl_modified')
    @cached_property
    def _get_ecbl_tex_mfn_vct(self):
        '''Get the callable function for effective crack brige law.
        '''
        return self.ecbl.mfn_vct

    u0 = Property(Array(float), depends_on = '+config_modified, ecbl_modified')
    @cached_property
    def _get_u0(self):
        return self.ecbl.u0

    sig_ti_arr = Property
    def _get_sig_ti_arr(self):
        '''force at the height of each reinforcement layer [kN]:
        '''
        return self.ecbl.mfn_vct(self.cs.eps_ti_arr)
    
    #===========================================================================
    # Layer conform discretization of the tensile zone
    #===========================================================================
    f_ti_arr = Property
    def _get_f_ti_arr(self):
        '''force at the height of each reinforcement layer [kN]:
        '''
        # tensile force of one reinforced composite layer [kN]:
        #
        sig_ti_arr = self.sig_ti_arr
        n_rovings = self.cs.n_rovings
        A_roving = self.cs.A_roving
        return sig_ti_arr * n_rovings * A_roving / 1000. 

    #===========================================================================
    # Discretization of the compressive zone - tex layer conform
    #===========================================================================
    f_ci_arr = Property
    def _get_f_ci_arr(self):
        '''Compressive stress in the compresive zone 'x' for each layer i.
        '''
        sig_ci_arr = self.sig_c_mfn_vct()
        return sig_ci_arr * self.cs.width * self.cs.s_tex_z * 1000.

    #===========================================================================
    # Fine discretization of the compressive zone
    #===========================================================================
    sig_cj_arr = Property
    def _get_sig_cj_arr(self):
        return self.sig_c_mfn_vct(self.cs.eps_cj_arr)

    f_cj_arr = Property
    def _get_f_cj_arr(self):

        sig_cj_arr = self.sig_cj_arr
        x = self.cs.x
        return sig_cj_arr * self.cs.width * x / self.cs.n_cj * 1000. 

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor = 'white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    @on_trait_change('+geo_input,+input')
    def _replot(self):
        print 'PLOTTING'
        axes = self.figure.gca()
        axes.clear()

        axes.plot(self.eps_ti_arr, self.z_ti_arr)
        self.data_changed = True


    view = View(HSplit(Group(
                Group(Item('thickness'),
                      Item('width'),
                      Item('n_layers'),
                      Item('n_rovings'),
                      Item('A_roving'),
                      label = 'Geometry',
                      ),
                Group(Item('eps_lo'),
                      Item('eps_up'),
                      label = 'Strain'
                      ),
                Group(Item('s_tex_z', label = 'vertical spacing', style = 'readonly'),
                      label = 'Layout',
                      ),
                             ),
                Group(Item('figure', editor = MPLFigureEditor(),
                           resizable = True, show_label = False),
                      id = 'simexdb.plot_sheet',
                      label = 'plot sheet',
                      dock = 'tab',
                      ),
                       ),
                resizable = True,
                buttons = ['OK', 'Cancel'])
    
if __name__ == '__main__':
    ecs = ECBCrossSection()
    ecs.configure_traits()
