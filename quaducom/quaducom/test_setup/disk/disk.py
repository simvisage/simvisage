'''
MUSHROOF - For double symmetric loading cases, one Roof is sufficient 

TODO: @ Andreas
     - split dead load cases, include additional dead load case
'''

from enthought.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool, HasTraits, Enum, Dict, Str

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel

from ibvpy.mats.mats2D import \
    MATS2DElastic

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

from simiter.sim_pstudy import\
    ISimModel, SimOut, SimPStudy, SimArray, SimArrayView

from mathkit.geo.transform.square2circle import GeoSquare2Circle

class DiskTestSetup(IBVModel):

    implements(ISimModel)

    material_density_roof = Float(-22.4e-3)    # [MN/m^3]

    #----------------------------------------------------
    # elements
    #----------------------------------------------------

    n_elems = Int(50, input = True)
    thickness = Float(1.0, input = True)

    vtk_r = Float(0.9, input = True)

    # fets used for roof
    #
    fets_disk = Instance((FETSEval), depends_on = '+ps_levels, +input')
    def _fets_disk_default(self):
        # fe_quad_serendipity_roof is defined in base class 
        # connected with material properties of the roof
        #
        mats = MATS2DElastic(E = 3.1e5, nu = 0.0, stress_state = 'plane_stress')
        fets = FETS2D4Q8U(mats_eval = mats)
        fets.vtk_r *= self.vtk_r
        return fets

    geo_disk = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_geo_disk(self):
        # list of local origins of each roof defined in global coordinates
        # (for MRDetail two roofs are considered due to symmetry)
        #
        gt = GeoSquare2Circle(circle_center = [0.0, 0.0],
                              circle_radius = 1.0,
                              square_edge = 4.0)
        return gt

    #----------------------------------------------------
    # fe-grids 
    #----------------------------------------------------

    fe_disk_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_fe_disk_grid(self):
        fe_grid = FEGrid(coord_min = (-1, -1),
                          coord_max = (1, 1),
                          geo_transform = self.geo_disk,
                          shape = (self.n_elems, self.n_elems),
                          fets_eval = self.fets_disk)
        return fe_grid

    bc_fixed = Property(depends_on = '+ps_levels, +input')
    @cached_property
    def _get_bc_fixed(self):
        fe_disk_grid = self.fe_disk_grid

        extension = 0.1
        alpha_45 = math.pi / 4.0

        bc00 = BCDof(var = 'u',
                     dof = fe_disk_grid[0, 0, 0, 0].dofs[0, 0, 1],
                     value = -extension * math.cos(alpha_45))

        bc00_link = BCDof(var = 'u',
                     dof = fe_disk_grid[0, 0, 0, 0].dofs[0, 0, 0],
                     link_dofs = [fe_disk_grid[0, 0, 0, 0].dofs[0, 0, 1]],
                     link_coeffs = [1.0],
                     value = 0)

        bc10 = BCDof(var = 'u',
                     dof = fe_disk_grid[-1, 0, -1, 0].dofs[0, 0, 1],
                     value = -extension * math.cos(alpha_45))

        bc10_link = BCDof(var = 'u',
                     dof = fe_disk_grid[-1, 0, -1, 0].dofs[0, 0, 0],
                     link_dofs = [fe_disk_grid[-1, 0, -1, 0].dofs[0, 0, 1]],
                     link_coeffs = [-1.0],
                     value = 0)

        bc11 = BCDof(var = 'u',
                     dof = fe_disk_grid[-1, -1, -1, -1].dofs[0, 0, 1],
                     value = extension * math.cos(alpha_45))

        bc11_link = BCDof(var = 'u',
                     dof = fe_disk_grid[-1, -1, -1, -1].dofs[0, 0, 0],
                     link_dofs = [fe_disk_grid[-1, -1, -1, -1].dofs[0, 0, 1]],
                     link_coeffs = [1.0],
                     value = 0)

        bc01 = BCDof(var = 'u',
                     dof = fe_disk_grid[0, -1, 0, -1].dofs[0, 0, 1],
                     value = extension * math.cos(alpha_45))

        bc01_link = BCDof(var = 'u',
                     dof = fe_disk_grid[0, -1, 0, -1].dofs[0, 0, 0],
                     link_dofs = [fe_disk_grid[0, -1, 0, -1].dofs[0, 0, 1]],
                     link_coeffs = [-1.0],
                     value = 0)

        n_xy = self.n_elems / 2

        bc_bottom = BCDof(var = 'u',
                     dof = fe_disk_grid[n_xy, 0, 0, 0].dofs[0, 0, 1],
                     value = -extension)
        bc_right = BCDof(var = 'u',
                     dof = fe_disk_grid[-1, n_xy, -1, 0].dofs[0, 0, 0],
                     value = extension)
        bc_top = BCDof(var = 'u',
                     dof = fe_disk_grid[n_xy, -1, 0, -1].dofs[0, 0, 1],
                     value = extension)
        bc_left = BCDof(var = 'u',
                     dof = fe_disk_grid[0, n_xy, 0, 0].dofs[0, 0, 0],
                     value = -extension)

        return [bc_left, bc_right, bc_top, bc_bottom,
                bc00, bc00_link, bc10, bc10_link, bc11, bc11_link, bc01, bc01_link]

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

        roof = self.fe_disk_grid

        # self.bc_roof_deadweight + \
        bc_list = self.bc_fixed

        ts = TS(sdomain = [ roof ],
                 dof_resultants = True,
                 bcond_list = bc_list,
                 rtrace_list = self.rtrace_list

               )

        # Add the time-loop control
        #
        tloop = TLoop(tstepper = ts,
                      tolerance = 1e-4,
                      tline = self.tline)

        return tloop

if __name__ == '__main__':

    sim_model = DiskTestSetup(n_elems = 40)

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
