'''
Created on Sep 4, 2012

@author: rch
'''
from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, Array, \
    Trait, Event, on_trait_change, Instance
    
from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit
    
from ecb_law import \
    ECBLBase, ECBLLinear, ECBLFBM, ECBLCubic, ECBLBilinear

from cc_law import \
    CCLawBase, CCLawBlock, CCLawLinear, CCLawQuadratic

import numpy as np

class ECBCrossSection(HasStrictTraits):

    #---------------------------------------------------------------
    # Cross section characteristics needed for tensile specimens 
    #---------------------------------------------------------------
    
    # thickness of reinforced cross section
    #
    thickness = Float(0.06, auto_set = False, enter_set = True, geo_input = True)

    # total number of reinforcement layers [-]
    # 
    n_layers = Int(12, auto_set = False, enter_set = True, geo_input = True)

    #---------------------------------------------------------------
    # Cross section characteristics needed for bending specimens 
    #---------------------------------------------------------------

    # width of the cross section [m]
    #
    width = Float(0.20, auto_set = False, enter_set = True, geo_input = True)
    
    # number of rovings in 0-direction of one composite 
    # layer of the bending test [-]:
    #
    n_rovings = Int(23, auto_set = False, enter_set = True, geo_input = True)
    
    # cross section of one roving [mm**2]:
    #
    A_roving = Float(0.461, auto_set = False, enter_set = True, geo_input = True)

    #===========================================================================
    # material properties 
    #===========================================================================

    f_ck = Float(55.7, auto_set = False, enter_set = True,
                 cc_input = True)
    # ultimate textile stress measured in the tensile test [MPa]
    #
    sig_tex_u = Float(1216., auto_set = False, enter_set = True,
                      tt_input = True)
    
    # compressive strain at the top at rupture [-]
    # (positive value is used)
    # value measured at the experiment (with DMS)
    # for 0: about 3.3
    # for 90: about 3.0
    # ultimate strain theoretically (Brockman): about 4.5
    # NOTE: strain was meassured at a distance of 5 cm
    #
    eps_cu = Float(0.0033, auto_set = False, enter_set = True,
                   cc_input = True) # float value corresponds to 3 promile


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
    n_cj = Int(200, auto_set = False, enter_set = True,
                 cc_input = True, eps_input = True)

    #===========================================================================
    # Strain state
    #===========================================================================

    eps_up = Float(0.0033, auto_set = False, enter_set = True, eps_input = True)
    eps_lo = Float(-0.0172, auto_set = False, enter_set = True, eps_input = True)
        
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
        
        # strain at the height of each reinforcement layer [-]:
        #
        return self.eps_up + (self.eps_lo - self.eps_up) * self.z_ti_arr / thickness

    x = Property(depends_on = '+eps_input,+geo_input')
    @cached_property
    def _get_x(self):
        # heights of the compressive zone:
        #
        return (abs(self.eps_up) / (abs(self.eps_up - self.eps_lo)) * 
                 self.thickness)
    
    eps_ti_arr = Property(depends_on = '+eps_input,+geo_input')
    @cached_property
    def _get_eps_ti_arr(self):
        return (-np.fabs(self.eps_i_arr) + self.eps_i_arr) / 2.0     
    
    eps_ci_arr = Property(depends_on = '+eps_input,+geo_input')
    @cached_property
    def _get_eps_ci_arr(self):
        return (np.fabs(self.eps_i_arr) + self.eps_i_arr) / 2.0  

    eps_cj_arr = Property(depends_on = '+eps_input,+geo_input,n_cj')    
    @cached_property
    def _get_eps_cj_arr(self):
        '''get compressive strain at each integration layer of the compressive zone [-]:
        for 'stress_case' flexion
        '''
        # for calibration us measured compressive strain
        # @todo: use mapped traits instead
        #
        eps_j_arr = (self.eps_up + (self.eps_lo - self.eps_up) * self.z_cj_arr / 
                     self.thickness)        
        return (np.fabs(eps_j_arr) + eps_j_arr) / 2.0  

    z_cj_arr = Property(depends_on = '+eps_input,+geo_input,n_cj')
    @cached_property
    def _get_z_cj_arr(self):
        '''Get the discretizaton of the  compressive zone
        '''
        return np.linspace(0, self.thickness, self.n_cj)

    # distance of reinforcement layers from the bottom 
    #
    zz_cj_arr = Property(depends_on = '+eps_input,+geo_input,n_cj')
    @cached_property
    def _get_zz_cj_arr(self):
        return self.thickness - self.z_cj_arr

    #===========================================================================
    # Compressive concrete constitutive law
    #===========================================================================
    cc_law_type = Trait('constant', dict(constant = CCLawBlock,
                                         linear = CCLawLinear,
                                         quadratic = CCLawQuadratic),
                        cc_input = True)
    
    cc_law = Property(Instance(CCLawBase), depends_on = 'cc_law_type')
    @cached_property
    def _get_cc_law(self):
        '''Construct the compressive concrete law'''
        return self.cc_law_type_(f_ck = self.f_ck)

    sig_c_mfn = Property(depends_on = '+cc_input')
    def _get_sig_c_mfn(self):
        return self.cc_law.mfn

    sig_c_mfn_vct = Property(depends_on = '+cc_input')
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
                      tt_input = True)
    
    ecbl = Property(Instance(ECBLBase), depends_on = 'ecbl_type, +tt_input')
    @cached_property
    def _get_ecbl(self):
        return self.ecbl_type_(sig_tex_u = self.sig_tex_u)

    ecbl_modified = Event
    @on_trait_change('ecbl.+input')
    def _set_ecbl_modified(self):
        self.ecbl_modified = True
    
    ecbl_tex_mfn = Property(depends_on = '+tt_input, ecbl_modified')
    @cached_property
    def _get_ecbl_tex_arr(self):
        '''Get the arrays sigma - epsilon defining the crack bridge law.
        '''
        return self.ecbl.arr
    
    ecbl_tex_mfn = Property(depends_on = '+tt_input, ecbl_modified')
    @cached_property
    def _get_ecbl_tex_mfn(self):
        '''Get the callable function for effective crack brige law.
        '''
        return self.ecbl.mfn

    ecbl_tex_mfn_vct = Property(depends_on = '+tt_input, ecbl_modified')
    @cached_property
    def _get_ecbl_tex_mfn_vct(self):
        '''Get the callable function for effective crack brige law.
        '''
        return self.ecbl.mfn_vct

    sig_ti_arr = Property
    def _get_sig_ti_arr(self):
        '''force at the height of each reinforcement layer [kN]:
        '''
        return self.ecbl.mfn_vct(self.eps_ti_arr)
    
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
        n_rovings = self.n_rovings
        A_roving = self.A_roving
        return sig_ti_arr * n_rovings * A_roving / 1000. 

    #===========================================================================
    # Discretization of the compressive zone - tex layer conform
    #===========================================================================
    f_ci_arr = Property
    def _get_f_ci_arr(self):
        '''Compressive stress in the compresive zone 'x' for each layer i.
        '''
        sig_ci_arr = self.sig_c_mfn_vct()
        return sig_ci_arr * self.width * self.s_tex_z * 1000.

    #===========================================================================
    # Fine discretization of the compressive zone
    #===========================================================================
    sig_cj_arr = Property
    def _get_sig_cj_arr(self):
        return self.sig_c_mfn_vct(self.eps_cj_arr)

    f_cj_arr = Property
    def _get_f_cj_arr(self):

        sig_cj_arr = self.sig_cj_arr
        return sig_cj_arr * self.width * self.thickness / self.n_cj * 1000. 

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor = 'white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    data_changed = Event
    
    @on_trait_change('+geo_input,+eps_input,+tt_input,+cc_input')
    def replot(self):
        
        d = self.thickness - self.x

        self.figure.clear()
        ax = self.figure.add_subplot(1, 2, 1)        
        
        #ax = self.figure.gca()
        
        # eps ti
        ax.plot([self.eps_lo, self.eps_up], [0, self.thickness], color = 'black')
        ax.hlines(self.zz_ti_arr, [0], self.eps_ti_arr, lw = 4, color = 'red')

        # eps cj
        ec = np.hstack([self.eps_cj_arr] + [0, 0])
        zz = np.hstack([self.zz_cj_arr] + [0, self.thickness]) 
        ax.fill(ec, zz, color = 'blue')

        # reinforcement layers
        eps_range = np.array([min(0.0, self.eps_lo),
                              max(0.0, self.eps_up)], dtype = 'float')
        z_ti_arr = np.ones_like(eps_range)[:, None] * self.z_ti_arr[None, :]
        ax.plot(eps_range, z_ti_arr, 'k--', color = 'black')
        
        # neutral axis
        ax.plot(eps_range, [d, d], 'k--', color = 'green', lw = 2)
        
        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax = self.figure.add_subplot(1, 2, 2)        

        # f ti
        ax.hlines(self.zz_ti_arr, [0], self.f_ti_arr, lw = 4, color = 'red')

        # f cj
        f_c = self.width * np.hstack([self.sig_cj_arr] + [0, 0])
        zz = np.hstack([self.zz_cj_arr] + [0, self.thickness]) 
        ax.fill(f_c, zz, color = 'blue')

        f_range = [np.min(self.f_ti_arr), np.max(f_c)]        
        # neutral axis
        ax.plot(f_range, [d, d], 'k--', color = 'green', lw = 2)
        
        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
        self.data_changed = True

    view = View(HSplit(Group(
                Group(Item('thickness'),
                      Item('width'),
                      Item('n_layers'),
                      Item('n_rovings'),
                      Item('A_roving'),
                      label = 'Geometry',
                      ),
                Group(Item('eps_up', label = 'Upper strain'),
                      Item('eps_lo', label = 'Lower strain'),
                      label = 'Strain'
                      ),
                Group(Item('f_ck', label = 'Compressive strength'),
                      Item('n_cj', label = 'Discretization'),
                      Item('cc_law_type', show_label = False),
                      Item('@cc_law', show_label = False),
                      label = 'Concrete'
                      ),
                Group(Item('ecbl_type', show_label = False),
                      Item('ecbl', show_label = False),
                      label = 'Reinforcement'
                      ),
                Group(Item('s_tex_z', label = 'vertical spacing', style = 'readonly'),
                      label = 'Layout',
                      ),
                             scrollable = True,
                             ),
                Group(Item('figure', editor = MPLFigureEditor(),
                           resizable = True, show_label = False),
                      id = 'simexdb.plot_sheet',
                      label = 'plot sheet',
                      dock = 'tab',
                      ),
                       ),
                width = 0.8,
                height = 0.5,
                resizable = True,
                buttons = ['OK', 'Cancel'])
    
if __name__ == '__main__':
    ecs = ECBCrossSection(# 7d: f_ck,cube = 62 MPa; f_ck,cyl = 62/1.2=52
                           # 9d: f_ck,cube = 66.8 MPa; f_ck,cyl = 55,7
                           f_ck = 55.7,

                           # measured strain at bending test rupture (0-dir)
                           #
                           eps_cu = 3.3 / 1000.,
                           ecbl_type = 'fbm',
                           cc_law_type = 'quadratic'
                           )
    ecs.replot()
    ecs.configure_traits()
