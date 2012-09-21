'''
Created on Sep 4, 2012

@todo: define views for ecbl and cc_law
@todo: introduce the dock feature for the views
@todo: classify the state changes and provide examples.
@todo: perform the calibration using ecbl_calib


@author: rch
'''
from etsproxy.traits.api import \
    HasStrictTraits, Float, Property, cached_property, Int, \
    Trait, Event, on_trait_change, Instance, Button
    
from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, VGroup, HGroup
    
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
    # State management
    #===========================================================================
    modified = Event
    @on_trait_change('+geo_input, cc_modified,+tt_input, +eps_input, tt_modified')
    def _set_cs_modified(self):
        self.modified = True
    
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
    n_cj = Int(20, auto_set = False, enter_set = True,
                 cc_input = True, eps_input = True)

    #===========================================================================
    # Strain state
    #===========================================================================

    eps_up = Float(-0.0033, auto_set = False, enter_set = True, eps_input = True)
    eps_lo = Float(0.0140, auto_set = False, enter_set = True, eps_input = True)
        
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
        return (np.fabs(self.eps_i_arr) + self.eps_i_arr) / 2.0  
    
    def convert_eps_tex_u_2_lo(self, eps_tex_u):
        eps_up = self.eps_up
        return eps_up + (eps_tex_u - eps_up) / self.z_ti_arr[0] * self.thickness
    
    def convert_eps_lo_2_tex_u(self, eps_lo):
        eps_up = self.eps_up
        return (eps_up + (eps_lo - eps_up) / self.thickness * self.z_ti_arr[0])
    
    eps_ci_arr = Property(depends_on = '+eps_input,+geo_input')
    @cached_property
    def _get_eps_ci_arr(self):
        return (-np.fabs(self.eps_i_arr) + self.eps_i_arr) / 2.0     

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
        return (-np.fabs(eps_j_arr) + eps_j_arr) / 2.0  

    z_cj_arr = Property(depends_on = '+eps_input,+geo_input,n_cj')
    @cached_property
    def _get_z_cj_arr(self):
        '''Get the discretizaton of the  compressive zone
        '''
        return np.linspace(0, self.x, self.n_cj)

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
    
    cc_law = Property(Instance(CCLawBase), depends_on = '+cc_input')
    @cached_property
    def _get_cc_law(self):
        '''Construct the compressive concrete law'''
        return self.cc_law_type_(f_ck = self.f_ck, cs = self)

    cc_modified = Event

    #===========================================================================
    # Fine discretization of the compressive zone
    #===========================================================================
    sig_cj_arr = Property(depends_on = '+eps_input, +cc_input, cc_modified')
    @cached_property
    def _get_sig_cj_arr(self):
        return -self.cc_law.mfn_vct(-self.eps_cj_arr)

    #===========================================================================
    # Effective crack bridge law
    #===========================================================================
    ecbl_type = Trait('fbm', dict(fbm = ECBLFBM,
                                  cubic = ECBLCubic,
                                  linear = ECBLLinear,
                                  bilinear = ECBLBilinear),
                      tt_input = True)
    
    ecbl = Property(Instance(ECBLBase), depends_on = '+tt_input')
    @cached_property
    def _get_ecbl(self):
        return self.ecbl_type_(sig_tex_u = self.sig_tex_u, cs = self)

    tt_modified = Event
    
    sig_ti_arr = Property(depends_on = '+eps_input, +tt_input, tt_modified')
    @cached_property    
    def _get_sig_ti_arr(self):
        '''force at the height of each reinforcement layer [kN]:
        '''
        return self.ecbl.mfn_vct(self.eps_ti_arr)
    
    #===========================================================================
    # Layer conform discretization of the tensile zone
    #===========================================================================
    f_ti_arr = Property(depends_on = '+eps_input, +tt_input, tt_modified')
    @cached_property
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
    # Stress resultants
    #===========================================================================

    N = Property(depends_on = '+geo_input,+eps_input,+tt_input,+cc_input, cc_modified, tt_modified')
    @cached_property
    def _get_N(self):
        '''Get the resulting normal force.
        '''
        N_tk = sum(self.f_ti_arr)
        N_ck = np.trapz(self.sig_cj_arr * self.width, self.z_cj_arr) * 1000.0
        
        N_internal = N_ck + N_tk
        return N_internal 

    M = Property(depends_on = '+geo_input,+eps_input,+tt_input, +cc_input, cc_modified, tt_modified')
    @cached_property
    def _get_M(self):
        M_tk = np.dot(self.f_ti_arr, self.z_ti_arr)
        M_ck = np.trapz(self.sig_cj_arr * self.width * self.z_cj_arr, self.z_cj_arr) * 1000.0

        M_internal_ = M_tk + M_ck

        # moment evaluated with respect to the center line
        #
        M_internal = M_internal_ - self.N * self.thickness / 2.
        
        # return the internal stress resultants

        return M_internal

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor = 'white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    data_changed = Event
    
    replot = Button
    def _replot_fired(self):
        
        d = self.thickness - self.x

        self.figure.clear()
        ax = self.figure.add_subplot(1, 2, 1)        
        
        #ax = self.figure.gca()
        
        # eps ti
        ax.plot([self.eps_lo, self.eps_up], [0, self.thickness], color = 'black')
        ax.hlines(self.zz_ti_arr, [0], self.eps_ti_arr, lw = 4, color = 'red')

        # eps cj
        ec = np.hstack([self.eps_cj_arr] + [0, 0])
        zz = np.hstack([self.zz_cj_arr] + [0, self.thickness ]) 
        ax.fill(ec, zz, color = 'blue')

        # reinforcement layers
        eps_range = np.array([max(0.0, self.eps_lo),
                              min(0.0, self.eps_up)], dtype = 'float')
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
        zz = np.hstack([self.zz_cj_arr] + [0, self.thickness ]) 
        ax.fill(f_c, zz, color = 'blue')

        f_range = [np.max(self.f_ti_arr), np.min(f_c)]        
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
                HGroup(
                Group(Item('thickness', springy = True),
                      Item('width'),
                      Item('n_layers'),
                      Item('n_rovings'),
                      Item('A_roving'),
                      label = 'Geometry',
                      springy = True
                      ),
                Group(Item('eps_up', label = 'Upper strain', springy = True),
                      Item('eps_lo', label = 'Lower strain'),
                      label = 'Strain',
                      springy = True
                      ),
                springy = True,
                ),
                HGroup(
                Group(HGroup(
                      Item('cc_law_type', show_label = False, springy = True),
                      Item('cc_law', label = 'Edit', show_label = False, springy = True),
                      springy = True
                      ),
                      Item('f_ck', label = 'Compressive strength'),
                      Item('n_cj', label = 'Discretization'),
                      label = 'Concrete',
                      springy = True
                      ),
                Group(HGroup(
                      Item('ecbl_type', show_label = False, springy = True),
                      Item('ecbl', label = 'Edit', show_label = False, springy = True),
                      springy = True,
                      ),
                      label = 'Reinforcement',
                      springy = True
                      ),
                springy = True,
                ),
                Group(Item('s_tex_z', label = 'vertical spacing', style = 'readonly'),
                      label = 'Layout',
                      ),
                Group(
                HGroup(Item('M', springy = True, style = 'readonly'),
                       Item('N', springy = True, style = 'readonly'),
                       ),
                       label = 'Stress resultants'
                       ),
                scrollable = True,
                             ),
                Group(Item('replot', show_label = False),
                      Item('figure', editor = MPLFigureEditor(),
                           resizable = True, show_label = False),
                      id = 'simexdb.plot_sheet',
                      label = 'plot sheet',
                      dock = 'tab',
                      ),
                       ),
                width = 0.8,
                height = 0.8,
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

    print 'initial'
    ecs.eps_up = -0.0033
    ecs.eps_lo = 0.014

    print 'z_cj'
    print ecs.z_cj_arr
    print ecs.eps_cj_arr
    print 'x', ecs.x
    
    eps_tex_u = ecs.convert_eps_lo_2_tex_u(0.014)
    print 'eps_tex_u', eps_tex_u
    eps_lo = ecs.convert_eps_tex_u_2_lo(eps_tex_u)
    print 'eps_lo', eps_lo

    print 'M', ecs.M
    print 'N', ecs.N
    
    ecs.configure_traits()
