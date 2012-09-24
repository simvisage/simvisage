'''
Created on Sep 4, 2012

@author: rch
'''
from etsproxy.traits.api import \
    HasTraits, Int, Instance, Property, cached_property, DelegatesTo, \
    Event, Button

from util.traits.editors.mpl_figure_editor import \
    MPLFigureEditor

from matplotlib.figure import \
    Figure

from etsproxy.traits.ui.api import \
    View, Item, Group, HSplit, VGroup, HGroup
    
from ecb_law_calib import \
    ECBLCalib
    
import numpy as np
    
class ECBLMNDiagram(HasTraits):

    # calibrator supplying the effective material law
    calib = Instance(ECBLCalib)
    def _calib_default(self):
        return ECBLCalib()
    
    # cross section
    cs = DelegatesTo('calib')

    calibrated_ecbl = Property
    def _get_calibrated_ecbl(self):
        return self.calib.calibrated_ecbl
    
    eps_cu = Property()
    def _get_eps_cu(self):
        return -self.cs.cc_law.eps_c_u
    
    eps_tu = Property()
    def _get_eps_tu(self):
        return self.calibrated_ecbl.eps_tex_u

    n_eps = Int(5, auto_set = False, enter_set = True)
    eps_range = Property
    def _get_eps_range(self):
        eps_c_space = np.linspace(self.eps_cu, 0, self.n_eps)
        eps_t_space = np.linspace(0, self.eps_tu, self.n_eps)

        eps_cc = self.eps_cu * np.ones_like(eps_c_space) 
        eps_ct = self.eps_cu * np.ones_like(eps_t_space) 
        eps_tc = self.eps_tu * np.ones_like(eps_c_space) 
        eps_tt = self.eps_tu * np.ones_like(eps_t_space) 
        
        eps1 = np.vstack([eps_c_space, eps_cc])
        eps2 = np.vstack([eps_t_space, eps_ct])
        eps3 = np.vstack([eps_tc, eps_c_space])
        eps4 = np.vstack([eps_tt, eps_t_space])
        
        return np.hstack([eps1, eps2, eps3, eps4])
        
        np.hstack(eps_c_space,)

    def _get_MN_fn(self, eps_lo, eps_up):
        print 'eps_lo', eps_lo, 'eps_up', eps_up
        self.cs.set(eps_lo = eps_lo,
                    eps_up = eps_up)
        return self.cs.M, self.cs.N
    
    MN_vct = Property(depends_on = 'calib.+config_changed')
    def _get_MN_vct(self):
        return np.vectorize(self._get_MN_fn)
        
    MN_arr = Property(depends_on = 'calib.+config_changed')
    @cached_property
    def _get_MN_arr(self):
        return self.MN_vct(self.eps_range[0, :], self.eps_range[1, :])

    figure = Instance(Figure)
    def _figure_default(self):
        figure = Figure(facecolor = 'white')
        figure.add_axes([0.08, 0.13, 0.85, 0.74])
        return figure

    data_changed = Event
    
    replot = Button
    def _replot_fired(self):
        
        self.figure.clear()
        ax = self.figure.add_subplot(1, 2, 1)        

        ax.plot(-self.eps_range, [0, 0.06], color = 'black')

        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax = self.figure.add_subplot(1, 2, 2)        

        ax.plot(self.MN_arr[0], -self.MN_arr[1], color = 'blue')
        
        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.grid(b = None, which = 'major')
        
        self.data_changed = True

    view = View(HSplit(Group(
                HGroup(
                Group(Item('n_eps', springy = True),
                      label = 'Discretization',
                      springy = True
                      ),
                springy = True,
                ),
                HGroup(
                Group(HGroup(
                      Item('cs', label = 'Cross section', show_label = False, springy = True),
                      springy = True,
                      ),
                      label = 'Cross sectoin',
                      springy = True
                      ),
                springy = True,
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
                height = 0.3,
                resizable = True,
                buttons = ['OK', 'Cancel'])

if __name__ == '__main__':
    c = ECBLCalib(
                  Mu = 3.49,
                  width = 0.20,
                  n_rovings = 23,
                  ecbl_type = 'fbm',
                  cc_law_type = 'linear'             #eps_tu 0.0137279096658                              
                  )
                     
    mn = ECBLMNDiagram(calib = c,
                       n_eps = 30,
                      )
    
    mn.configure_traits()

