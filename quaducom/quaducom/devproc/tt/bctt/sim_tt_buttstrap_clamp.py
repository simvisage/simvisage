'''
Optimization of tensile test buttstrap clamping
'''

from etsproxy.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool, HasTraits, Enum, Dict, Str

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel, FEDomain, FERefinementGrid

from ibvpy.mats.mats2D import \
    MATS2DElastic

from ibvpy.fets.fets1D5.fets1D52l4uLRH import \
    MATS1DElastic, MATS1DPlastic, MATS1D5Bond

from ibvpy.fets.fets1D5.fets1D52l6uLRH import \
    FETS1D52L6ULRH

from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets2D import \
    FETS2D4Q8U

from ibvpy.mesh.fe_grid import \
    FEGrid

import math

from ibvpy.mesh.fe_spring_array import \
     FESpringArray

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

import numpy as np

# Interpolation
from scipy.interpolate import Rbf

from math import \
    sqrt, asin, acos, pi as Pi

from matplotlib.pyplot import \
    bar, show, axhline, ion, ioff, xlabel, ylabel, title, figure, savefig

from matresdev.simiter.sim_pstudy import\
    ISimModel, SimOut, SimPStudy, SimArray, SimArrayView

from mathkit.geo.transform.square2circle import GeoSquare2Circle

class ButtstrapClamping(IBVModel):

    implements(ISimModel)

    #===========================================================================
    # Speciment geometry and discretization parameters
    #===========================================================================
    specimen_ne_x = Int(2, input = True,
                        label = 'number of elements in x-direction')
    specimen_ne_y = Int(2, input = True,
                        label = 'number of elements in y-direction')
    specimen_thickness = Float(0.01, input = True,
                           label = 'thickness of the tensile specimen')
    specimen_length = Float(0.6, input = True,
                        label = 'length of the tensile specimen')

    specimen_width = Float(0.1, input = True,
                        label = 'length of the tensile specimen')

    #===========================================================================
    # Buttstrap geometry and discretization parameters
    #===========================================================================
    buttstrap_length = Int(0.2, input = True,
                           label = 'length of the buttstrap')
    buttstrap_ne_y = Int(4, input = True,
                          label = 'number of elements in buttstrap in y-direction')
    buttstrap_ne_x = Int(7, input = True,
                          label = 'number of elements in buttstrap in x-direction')

    buttstrap_max_thickness = Int(0.04, input = True,
                                  label = 'maximum thickness of the buttstrap')
    buttstrap_min_thickness = Int(0.004, input = True,
                                  label = 'mimnimum thickness of the buttstrsp')

    #===========================================================================
    # Elastomer geometry
    #===========================================================================
    elastomer_thickness = Float(0.002, input = True,
                                label = 'thickness of the elastomer')

    elastomer_ne_y = Int(1, input = True,
                          label = 'number of elements in buttstrap in y-direction')
    elastomer_ne_x = Property
    def _get_elastomer_ne_x(self):
        return self.buttstrap_ne_x

    #===========================================================================
    # Friction interface geometry
    #===========================================================================
    friction_thickness = Float(0.005, input = True,
                                label = 'thickness of the elastomer')

    friction_ne_y = Int(1, input = True,
                          label = 'number of elements in buttstrap in y-direction')
    friction_ne_x = Property
    def _get_friction_ne_x(self):
        return self.buttstrap_ne_x

    vtk_r = Float(0.9, input = True)

    #===========================================================================
    # Element types
    #===========================================================================

    E_specimen = Int(31000.0, input = True,
                     label = 'E Modulus of the specimen [MPa]')

    nu_specimen = Int(0.2, input = True,
                      label = 'E Modulus of the specimen [-]')

    specimen_fets = Instance((FETSEval), depends_on = '+ps_levels, +input')
    def _specimen_fets_default(self):
        # fe_quad_serendipity_roof is defined in base class 
        # connected with material properties of the roof
        #
        E_eff = self.E_specimen * self.specimen_width
        mats = MATS2DElastic(E = E_eff, nu = self.nu_specimen,
                             stress_state = 'plane_stress')
        fets = FETS2D4Q8U(mats_eval = mats)
        fets.vtk_r *= self.vtk_r
        return fets

    E_buttstrap = Int(210000.0, input = True,
                     label = 'E Modulus of the buttstrap [MPa]')

    nu_buttstrap = Int(0.2, input = True,
                       label = 'Poisson ratio of the buttstrap [-]')

    buttstrap_fets = Instance((FETSEval), depends_on = '+ps_levels, +input')
    def _buttstrap_fets_default(self):
        # fe_quad_serendipity_roof is defined in base class 
        # connected with material properties of the roof
        #
        E_eff = self.E_buttstrap * self.specimen_width
        mats = MATS2DElastic(E = E_eff, nu = self.nu_buttstrap,
                             stress_state = 'plane_stress')
        fets = FETS2D4Q8U(mats_eval = mats)
        fets.vtk_r *= self.vtk_r
        return fets

    E_elastomer = Int(100.0, input = True,
                     label = 'E Modulus of the elastomer [MPa]')

    nu_elastomer = Int(0.3, input = True,
                       label = 'Poisson ratio of the elastomer [-]')

    elastomer_fets = Instance((FETSEval), depends_on = '+ps_levels, +input')
    def _elastomer_fets_default(self):
        # fe_quad_serendipity_roof is defined in base class 
        # connected with material properties of the roof
        #
        E_eff = self.E_buttstrap * self.specimen_width
        mats = MATS2DElastic(E = E_eff, nu = self.nu_elastomer,
                             stress_state = 'plane_stress')
        fets = FETS2D4Q8U(mats_eval = mats)
        fets.vtk_r *= self.vtk_r
        return fets

    g = Float(400000, input = True,
            label = 'Shear stiffness between elastomer and specimen')

    t_max = Float(50, input = True,
            label = 'Frictional stress between elastomer and specimen')

    friction_fets = Instance((FETSEval), depends_on = '+ps_levels, +input')
    def _friction_fets_default(self):
        width = self.specimen_width
        G = self.g * width
        T_max = self.t_max * width

        ifopen_law = MFnLineArray(ydata = [ -1.0e+1, 0., 1.0e+10 ],
                                  xdata = [ -1., 0., 1.])
        mats_ifopen = MATS1DElastic(stress_strain_curve = ifopen_law)

        mats = MATS1D5Bond(mats_phase1 = MATS1DElastic(E = 0),
                           mats_phase2 = MATS1DElastic(E = 0),
                           mats_ifslip = MATS1DPlastic(E = G,
                                                       sigma_y = T_max,
                                                       K_bar = 0.,
                                                       H_bar = 0.),
                           mats_ifopen = mats_ifopen
                           )

        fets = FETS1D52L6ULRH(mats_eval = mats) #quadratic
        fets.vtk_r *= self.vtk_r
        return fets

    #===========================================================================
    # Geometry transformation
    #===========================================================================
    specimen_cl_geo = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_specimen_cl_geo(self):

        def geo_transform(points):
            points[:, 0] *= self.buttstrap_length
            points[:, 1] *= self.specimen_thickness
            return points

        return geo_transform

    specimen_fl_geo = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_specimen_fl_geo(self):

        def geo_transform(points):
            points[:, 0] *= self.specimen_length
            points[:, 0] += self.buttstrap_length
            points[:, 1] *= self.specimen_thickness
            return points

        return geo_transform

    buttstrap_geo = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_buttstrap_geo(self):

        def geo_transform(points):
            x, y = points.T
            x *= self.buttstrap_length
            h_max = self.buttstrap_max_thickness
            h_min = self.buttstrap_min_thickness
            y *= (h_max - (h_max - h_min) / self.buttstrap_length * x)
            y_offset = (self.specimen_thickness +
                        self.friction_thickness +
                        self.elastomer_thickness)
            points[:, 0], points[:, 1] = x, y + y_offset
            return points

        return geo_transform

    buttstrap_clamp_geo = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_buttstrap_clamp_geo(self):

        def geo_transform(points):
            x, y = points.T
            x *= self.buttstrap_length
            x -= self.buttstrap_length
            y *= self.buttstrap_max_thickness
            y_offset = (self.specimen_thickness +
                        self.friction_thickness +
                        self.elastomer_thickness)
            points[:, 0], points[:, 1] = x, y + y_offset
            return points

        return geo_transform

    elastomer_geo = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_elastomer_geo(self):

        def geo_transform(points):
            x, y = points.T
            x *= self.buttstrap_length
            h = self.elastomer_thickness
            y *= h
            y_offset = self.specimen_thickness + self.friction_thickness
            points[:, 0], points[:, 1] = x, y + y_offset
            return points

        return geo_transform

    friction_geo = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_friction_geo(self):

        def geo_transform(points):
            x, y = points.T
            x *= self.buttstrap_length
            h = self.friction_thickness
            y *= h
            y_offset = self.specimen_thickness
            points[:, 0], points[:, 1] = x, y + y_offset
            return points

        return geo_transform

    #===========================================================================
    # FE-grids 
    #===========================================================================

    fe_domain = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_domain(self):
        return FEDomain()

    specimen_fl_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_specimen_fl_fe_level(self):
        return  FERefinementGrid(name = 'specimen free level',
                                 fets_eval = self.specimen_fets,
                                 domain = self.fe_domain)

    specimen_cl_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_specimen_cl_fe_level(self):
        return  FERefinementGrid(name = 'specimen clamped level',
                                 fets_eval = self.specimen_fets,
                                 domain = self.fe_domain)

    buttstrap_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_buttstrap_fe_level(self):
        return FERefinementGrid(name = 'buttstrap level',
                                       fets_eval = self.buttstrap_fets,
                                       domain = self.fe_domain)

    buttstrap_clamp_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_buttstrap_clamp_fe_level(self):
        return FERefinementGrid(name = 'buttstrap clamp level',
                                       fets_eval = self.buttstrap_fets,
                                       domain = self.fe_domain)

    elastomer_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_elastomer_fe_level(self):
        return FERefinementGrid(name = 'elastomer level',
                                       fets_eval = self.elastomer_fets,
                                       domain = self.fe_domain)

    friction_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_friction_fe_level(self):
        return FERefinementGrid(name = 'friction level',
                                       fets_eval = self.friction_fets,
                                       domain = self.fe_domain)

    specimen_fl_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_specimen_fl_fe_grid(self):
        fe_grid = FEGrid(coord_min = (0, 0),
                         coord_max = (1, 1),
                         level = self.specimen_fl_fe_level,
                         geo_transform = self.specimen_fl_geo,
                         shape = (self.specimen_ne_x, self.specimen_ne_y),
                         fets_eval = self.specimen_fets)
        return fe_grid

    specimen_cl_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_specimen_cl_fe_grid(self):
        fe_grid = FEGrid(coord_min = (0, 0),
                         coord_max = (1, 1),
                         level = self.specimen_cl_fe_level,
                         geo_transform = self.specimen_cl_geo,
                         shape = (self.buttstrap_ne_x, self.specimen_ne_y),
                         fets_eval = self.specimen_fets)
        return fe_grid

    buttstrap_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_buttstrap_fe_grid(self):
        fe_grid = FEGrid(coord_min = (0, 0),
                         coord_max = (1, 1),
                         level = self.buttstrap_fe_level,
                         geo_transform = self.buttstrap_geo,
                         shape = (self.buttstrap_ne_x, self.buttstrap_ne_y),
                         fets_eval = self.buttstrap_fets)
        return fe_grid

    buttstrap_clamp_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_buttstrap_clamp_fe_grid(self):
        fe_grid = FEGrid(coord_min = (0, 0),
                         coord_max = (1, 1),
                         level = self.buttstrap_clamp_fe_level,
                         geo_transform = self.buttstrap_clamp_geo,
                         shape = (1, self.buttstrap_ne_y),
                         fets_eval = self.buttstrap_fets)
        return fe_grid

    elastomer_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_elastomer_fe_grid(self):
        fe_grid = FEGrid(coord_min = (0, 0),
                         coord_max = (1, 1),
                         level = self.elastomer_fe_level,
                         geo_transform = self.elastomer_geo,
                         shape = (self.elastomer_ne_x, self.elastomer_ne_y),
                         fets_eval = self.elastomer_fets)
        return fe_grid

    friction_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_friction_fe_grid(self):
        fe_grid = FEGrid(coord_min = (0, 0),
                         coord_max = (1, 1),
                         level = self.friction_fe_level,
                         geo_transform = self.friction_geo,
                         shape = (self.friction_ne_x, self.friction_ne_y),
                         fets_eval = self.friction_fets)
        return fe_grid

    #===========================================================================
    # Boundary conditions
    #===========================================================================
    bc_list = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_list(self):

        sp_fl_grid = self.specimen_fl_fe_grid
        sp_cl_grid = self.specimen_cl_fe_grid
        bs_grid = self.buttstrap_fe_grid
        bs_clamp_grid = self.buttstrap_clamp_fe_grid
        el_grid = self.elastomer_fe_grid
        fl_grid = self.friction_fe_grid

        link_sp_fl_cl = BCDofGroup(var = 'u', value = 0., dims = [0, 1],
                                    get_dof_method = sp_cl_grid.get_right_dofs,
                                    get_link_dof_method = sp_fl_grid.get_left_dofs,
                                    link_coeffs = [1.])

        link_fl_sp = BCDofGroup(var = 'u', value = 0., dims = [0, 1],
                                    get_dof_method = fl_grid.get_bottom_dofs,
                                    get_link_dof_method = sp_cl_grid.get_top_dofs,
                                    link_coeffs = [1.])

        link_el_fl = BCDofGroup(var = 'u', value = 0., dims = [0, 1],
                                    get_dof_method = el_grid.get_bottom_dofs,
                                    get_link_dof_method = fl_grid.get_top_dofs,
                                    link_coeffs = [1.])

        link_bs_el = BCDofGroup(var = 'u', value = 0., dims = [0, 1],
                                    get_dof_method = bs_grid.get_bottom_dofs,
                                    get_link_dof_method = el_grid.get_top_dofs,
                                    link_coeffs = [1.])

        link_bs_clamp = BCDofGroup(var = 'u', value = 0., dims = [0, 1],
                                    get_dof_method = bs_clamp_grid.get_right_dofs,
                                    get_link_dof_method = bs_grid.get_left_dofs,
                                    link_coeffs = [1.])

        symx_bc = BCSlice(var = 'u',
                          slice = sp_fl_grid[-1, :, -1, :],
                          dims = [0],
                          value = 0)

        symy_fl_bc = BCSlice(var = 'u',
                              slice = sp_fl_grid[:, 0, :, 0],
                              dims = [1],
                              value = 0)

        symy_cl_bc = BCSlice(var = 'u',
                              slice = sp_cl_grid[:, 0, :, 0],
                              dims = [1],
                              value = 0)

        cntl_bc = BCSlice(var = 'u',
                          slice = bs_clamp_grid[0, :, 0, :],
                          dims = [0],
                          value = -0.001)

        return [symx_bc, symy_fl_bc, symy_cl_bc,
                cntl_bc, link_sp_fl_cl, link_bs_clamp,
                link_fl_sp, link_el_fl, link_bs_el]

    #----------------------------------------------------
    # ps_study
    #----------------------------------------------------
    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''

        # DISPLACEMENT
        #
        U = self.tloop.eval()

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [  SimOut(name = '$u_z$', unit = '[mm]'), ]

    #----------------------------------------------------
    # response tracer
    #----------------------------------------------------

    rtrace_list = List
    def _rtrace_list_default(self):
        return [  self.max_princ_stress, self.sig_app, self.u]

    max_princ_stress = Instance(RTraceDomainListField)
    def _max_princ_stress_default(self):
        return RTraceDomainListField(name = 'max principle stress' , idx = 0,
                                      var = 'max_principle_sig', warp = True,
#                                      position = 'int_pnts',
                                      record_on = 'update',)

    sig_app = Property(Instance(RTraceDomainListField), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_sig_app(self):
        return RTraceDomainListField(name = 'sig_app' ,
                                      position = 'int_pnts',
                                      var = 'sig_app',
                                      record_on = 'update',)

    u = Property(Instance(RTraceDomainListField), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_u(self):
        return RTraceDomainListField(name = 'displacement' ,
                                      var = 'u', warp = True,
                                      record_on = 'update',)

    #----------------------------------------------------
    # time loop
    #----------------------------------------------------

    tline = Instance(TLine)
    def _tline_default(self):
        return TLine(min = 0.0, step = 1.0, max = 1.0)

    tloop = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_tloop(self):

        ts = TS(sdomain = self.fe_domain,
                 dof_resultants = True,
                 bcond_list = self.bc_list,
                 rtrace_list = self.rtrace_list

               )

        # Add the time-loop control
        #
        tloop = TLoop(tstepper = ts,
                      tolerance = 1e-4,
                      tline = self.tline)

        return tloop

if __name__ == '__main__':

    sim_model = ButtstrapClamping()

    do = 'ui'

    if do == 'ui':
        print '*** ui ***'

        # evaluation
        #
        sim_model.peval()

        # visualisation
        # 
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource = sim_model)
        app.main()
