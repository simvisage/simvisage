from etsproxy.traits.api import \
    Array, Bool, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, \
    Callable, List, TraitDict, Any, Range, \
    Delegate, Event, on_trait_change, Button, \
    Interface, implements, Property, cached_property

from ibvpy.api import \
    TStepper as TS, TLoop, TLine, \
    IBVModel, DOTSEval, \
    RTraceGraph, RTraceDomainListField, \
    BCDof, BCDofGroup, BCSlice, FERefinementGrid, FEDomain

from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U
from ibvpy.fets.fets3D.fets3D8h32u import \
    FETS3D8H32U
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U
from ibvpy.fets.fets2D5.fets2D58h32u import \
    FETS2D58H32U

from ibvpy.mesh.fe_grid import \
    FEGrid

from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
    MATS2D5MicroplaneDamage

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
    IPhiFn, PhiFnGeneralExtended, \
    PhiFnGeneral, PhiFnStrainHardening, PhiFnStrainHardeningLinear, \
    PhiFnGeneralExtendedExp, PhiFnStrainSoftening

from mathkit.geo.geo_ndgrid import \
    GeoNDGrid

from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import \
    MFnNDGrid, GridPoint

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

from numpy import \
    sin, cos, c_, arange, hstack, array, max, frompyfunc, linspace

import numpy as np

from time import time
from os.path import join

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from matresdev.simiter.sim_pstudy import \
    SimPStudy, SimOut, ISimModel

from matresdev.db.exdb.ex_run_view import \
    ExRunView

from matresdev.db.matdb.trc.ccs_unit_cell import \
    CCSUnitCell, DamageFunctionEntry

from matresdev.db.simdb import \
    SimDB

from matresdev.db.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

simdb = SimDB()

from geo_st import \
    GeoST

from geo_lip import \
    GeoLIP

from geo_supprt import \
    GeoSUPPRT

from pickle import dump, load

from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_name('Script MT')
font.set_family('serif')
font.set_style('normal')
font.set_size('large')
font.set_variant('normal')
font.set_weight('medium')

def format_plot(axes, xlim=None, ylim=None, xlabel='', ylabel=''):
    '''format 2d-plot black and with with times legends
    '''
    #-------------------------------------------------------------------
    # configure the style of the font to be used for labels and ticks
    #-------------------------------------------------------------------
    #
    from matplotlib.font_manager import FontProperties
    font = FontProperties()
    font.set_name('Script MT')
    font.set_family('serif')
    font.set_style('normal')
#    font.set_size('small')
    font.set_size('large')
    font.set_variant('normal')
    font.set_weight('medium')

    if xlim != None and ylim != None:
        axes.axis([0, xlim, 0., ylim], fontproperties=font)

    # format ticks for plot
    #
    locs, labels = axes.xticks()
    axes.xticks(locs, map(lambda x: "%.0f" % x, locs), fontproperties=font)
    axes.xlabel(xlabel, fontproperties=font)

    locs, labels = axes.yticks()
    axes.yticks(locs, map(lambda x: "%.0f" % x, locs), fontproperties=font)
    axes.ylabel(ylabel, fontproperties=font)


class SimST(IBVModel):
    '''Slab test prepared for parametric study. Due to symmetry only a quarter of the slab is simulated.
    @params:
    'geo_st_flag': specify weather a geometric transformation of the mesh of the slab is to be performed
                   in order to capture the round load introduction area

    'elstrm_flag': specify weather a elastomer is modeled for the load introduction

    'supprt_flag': specify weather the support is modeled as separate fe-grid allowing for realistic
                   kinematics at the support (rotation/translation in the plane between slab and support
    '''
    i = 0
    input_change = Event
    @on_trait_change('+input,ccs_unit_cell.input_change')
    def _set_input_change(self):
        print 'xxx input_change Event', self.i
        self.i += 1
        self.input_change = True

    implements(ISimModel)

    #-----------------
    # discretization:
    #-----------------

    # slab discretization in x,y-direction (quarter of the slab):
    #
    shape_xy = Int(14, input=True,
                      ps_levels=(8, 12, 3))

    # discretization of the load introduction plate (region 'R')
    #
    shape_R = Int(3, input=True)

    # slab discretization in z-direction:
    #
    shape_z = Int(3, input=True,
                      ps_levels=(2, 3, 3))

    # support discretization in xy-direction:
    # NOTE: chose '2' for 2x2-grid or '4' for 4x4-grid
    #
    shape_supprt_xy = Int(2, input=True,
                      ps_levels=(2, 4, 2))

    #-----------------
    # geometry:
    #-----------------

    # get an element of the support mesh adjacent to the center node, i.e. the fixed node
    # NOTE: property only distinguish weather support mesh is discretized with 2x2 or 4x4 grid
    #
    idx_supprt_center_elem = Property(Float, depends_on='shape_supprt_xy')
    @cached_property
    def _get_idx_supprt_center_elem(self):
        return self.shape_supprt_xy / 2

    # edge length of the quadratic plate (entire length without symmetry)
    #
    length = Float(0.80, input=True)
    thickness = Float(0.01714, input=True)
    thickness_plain_concrete = Float(0.00286, input=True)

    elem_length_xy = Property(Float, depends_on='shape_xy, length')
    @cached_property
    def _get_elem_length_xy(self):
        return self.length / self.shape_xy

    # properties of the load introduction plate
    # NOTE: only available when geometric transformation for slab test is used
    #
    radius_plate = Float(0.095, input=True)
    thickness_plate = Float(0.02, input=True)

    # thickness of the tappered support
    #
    thickness_supprt = Float(0.02, input=True)

    #-----------------
    # specify the refinement of the idealization
    #-----------------

    # specify weather rounded load introduction plate is modeled. If set to 'False' a complete square element mesh is used.
    #
    plain_concrete_flag = False

    # specify weather rounded load introduction plate is modeled. If set to 'False' a complete square element mesh is used.
    #
    geo_st_flag = True

    # specify weather elastomer is to be modeled for load introduction
    #
    elstmr_flag = False

    # specify weather steel support is to be modeled
    #
    supprt_flag = True

    #-----------------
    # 'geo_transform'
    #-----------------

    # geometry transformation for slab with rounded load introduction plate
    #
    geo_st = Property(Instance(GeoST), depends_on='+ps_levels, +input')
    @cached_property
    def _get_geo_st(self):
        return GeoST(shape_xy=self.shape_xy,
                     shape_R=self.shape_R,
                     shape_z=self.shape_z,
                     radius_plate=self.radius_plate,
                     length_quarter=self.length / 2.,
                     thickness=self.thickness,
                     zoffset=self.thickness_plain_concrete)

    # geometry transformation for slab with rounded load introduction plate
    # (for unreinforced concrete layer corresponding to 's_tex,z')
    #
    geo_st_plain_concrete = Property(Instance(GeoST), depends_on='+ps_levels, +input')
    @cached_property
    def _get_geo_st_plain_concrete(self):
        return GeoST(shape_xy=self.shape_xy,
                     shape_R=self.shape_R,
                     shape_z=1,
                     radius_plate=self.radius_plate,
                     length_quarter=self.length / 2.,
                     thickness=self.thickness_plain_concrete,
                     zoffset=0.)

    # geometry transformation rounded load introduction plate (if modeled)
    #
    geo_elstmr = Property(Instance(GeoLIP), depends_on='+ps_levels, +input')
    @cached_property
    def _get_geo_elstmr(self):
        return GeoLIP(radius_plate=self.radius_plate,
                      thickness_plate=self.thickness_plate,
                      xyoffset=self.length / 2. - self.radius_plate,
                      zoffset=self.thickness)

    # geometry transformation for four sided tappered support with reduced stiffness towards the edges (if modeled)
    #
    geo_supprt = Property(Instance(GeoSUPPRT), depends_on='+ps_levels, +input')
    @cached_property
    def _get_geo_supprt(self):
        # element length within the range of the slab without the area of the
        # load introduction plate
        #
        elem_size = (self.length / 2.0 - self.radius_plate) / (self.shape_xy - self.shape_R)
        width_supprt = self.shape_supprt_xy * elem_size
        print 'width_supprt = ', width_supprt
        return GeoSUPPRT(thickness_supprt=self.thickness_supprt,
                         width_supprt=width_supprt,
                         xyoffset=0.,
                         zoffset=-self.thickness_supprt)

    #----------------------------------------------------------------------------------
    # mats
    #----------------------------------------------------------------------------------

    # age of the specimen at the time of testing
    # determines the E-modulus and nu based on the time dependent function stored
    # in 'CCSUniteCell' if params are not explicitly set
    #
    age = Int(28, auto_set=False, enter_set=True, input=True)

    # composite E-modulus / Poisson's ratio
    # NOTE 1: value are retrieved from the database
    #         same values are used for calibration (as taken from tensile test) and for slab test
    # NOTE 2: alternatively the same phi-function could be used independent of age. This would assume
    #         an afine/proportional damage evolution for different ages, i.e. a different E would be
    #         used in the simulation of the slab test and within the calibration procedure.

    # time stepping params
    #
    tstep = Float(0.05, auto_set=False, enter_set=True, input=True)
    tmax = Float(1.0, auto_set=False, enter_set=True, input=True)
    tolerance = Float(0.001, auto_set=False, enter_set=True, input=True)

    # specify type of 'linalg.norm'
    # default value 'None' sets norm to 2-norm,
    # i.e "norm = sqrt(sum(x_i**2))
    #
    # set 'ord=np.inf' to switch norm to
    # "norm = max(abs(x_i))"
    #
    ord = Enum(None, np.inf)

    n_mp = Int(30., auto_set=False, enter_set=True, input=True)

#    E_specmn_mats = Float(input=True)

    # @todo: for mats_eval the information of the unit cell should be used
    # in order to use the same number of microplanes and model version etc...
    #
    specmn_mats = Property(Instance(MATS2D5MicroplaneDamage),
                          depends_on='input_change')
    @cached_property
    def _get_specmn_mats(self):
#        import pylab as p
#        self.phi_fn.mfn.mpl_plot(p)
#        p.show()
        print 'self.n_mp', self.n_mp
        specmn_mats = MATS2D5MicroplaneDamage(
                                E=self.E_c,
#                                E=self.E_specmn_mats,
#                                 E=self.E_m,  # relevant for compressive behavior/used for calibration of phi_fn
#                                 E=self.E_m,  # #@todo: only available if this parameter is defined by the class
                                nu=self.nu,
                                n_mp=self.n_mp,
                                symmetrization='sum-type',
                                model_version='compliance',
                                phi_fn=self.phi_fn)
        return specmn_mats

    plain_concrete_mats = Property(Instance(MATS2D5MicroplaneDamage),
                          depends_on='input_change')
    @cached_property
    def _get_plain_concrete_mats(self):
#        import pylab as p
#        self.phi_fn.mfn.mpl_plot(p)
#        p.show()
        print 'self.n_mp', self.n_mp
        specmn_mats = MATS2D5MicroplaneDamage(
                                E=self.E_c,
#                                 E=self.E_m,  # relevant for compressive behavior/used for calibration of phi_fn
#                                 E=self.E_m,  # #@todo: only available if this parameter is defined by the class
                                nu=self.nu,
                                n_mp=self.n_mp,
                                symmetrization='sum-type',
                                model_version='compliance',
                                phi_fn=PhiFnStrainSoftening)
        return specmn_mats

    elstmr_mats = Property(Instance(MATS3DElastic),
                                    depends_on='input_change')
    @cached_property
    def _get_elstmr_mats(self):
#        max_f = 0.020 # MN
#        max_eps = 1.0 # [-]
#        area = self.elem_length_xy ** 2
#        sig = max_f / area
#        E_elast = sig / max_eps
        E_elast = self.E_c / 10.
        print 'effective elastomer E_modulus', E_elast
        return MATS3DElastic(E=E_elast,
                             nu=0.2)
#                             nu = 0.4)

    # E-modulus and nu of steel support
    E_s = Float(210000., auto_set=False, enter_set=True, input=True)
    nu_s = Float(0.20, auto_set=False, enter_set=True, input=True)

    supprt_mats = Property(Instance(MATS3DElastic),
                                    depends_on='input_change')
    @cached_property
    def _get_supprt_mats(self):
        return MATS3DElastic(E=self.E_s,
                             nu=self.nu_s)

    #-----------------
    # fets:
    #-----------------

    # use quadratic serendipity elements
    # NOTE: 2D5 elements behave linear elastic in out of plane direction!
    #
    specmn_fets = Property(Instance(FETSEval),
                           depends_on='input_change')
    @cached_property
    def _get_specmn_fets(self):
        # NOTE: chose number of notes correspondingly to 'supprt_fets'!
        # quadratic serendipity elements:
#         fets = FETS2D58H20U(mats_eval=self.specmn_mats)
        # cubic serendipity elements:
        fets = FETS2D58H32U(mats_eval=self.specmn_mats)
        fets.vtk_r *= self.vtk_r
        return fets

    plain_concrete_fets = Property(Instance(FETSEval),
                           depends_on='input_change')
    @cached_property
    def _get_plain_concrete_fets(self):
        fets = FETS2D58H32U(mats_eval=self.plain_concrete_mats)
        fets.vtk_r *= self.vtk_r
        return fets

    elstmr_fets = Property(Instance(FETSEval),
                           depends_on='input_change')
    @cached_property
    def _get_elstmr_fets(self):
        fets = FETS2D58H32U(mats_eval=self.elstmr_mats)
        fets.vtk_r *= self.vtk_r
        return fets

    supprt_fets = Property(Instance(FETSEval),
                           depends_on='input_change')
    @cached_property
    def _get_supprt_fets(self):
        # linear-elastic behavior quadratic serendipity elements
        # quadratic serendipity elements:
#         fets = FETS3D8H20U(mats_eval=self.supprt_mats)
        # cubic serendipity elements:
        fets = FETS3D8H32U(mats_eval=self.supprt_mats)
        fets.vtk_r *= self.vtk_r
        return fets

    #-----------------
    # fe_grid:
    #-----------------

    # specify element shrink factor in plot of fe-model
    #
    vtk_r = Float(0.95)

    fe_domain = Property(depends_on='+ps_levels, +input')
    @cached_property
    def _get_fe_domain(self):
        return FEDomain()

    specmn_fe_level = Property(depends_on='+ps_levels, +input')
    @cached_property
    def _get_specmn_fe_level(self):
        return FERefinementGrid(name='specimen patch',
                                fets_eval=self.specmn_fets,
                                domain=self.fe_domain)

    plain_concrete_fe_level = Property(depends_on='+ps_levels, +input')
    @cached_property
    def _get_plain_concrete_fe_level(self):
        return FERefinementGrid(name='plain_concrete patch',
                                fets_eval=self.plain_concrete_fets,
                                domain=self.fe_domain)

    specmn_fe_grid = Property(Instance(FEGrid), depends_on='+ps_levels, +input')
    @cached_property
    def _get_specmn_fe_grid(self):
        # if 'geo_st' is used for modeling of round load introduction area
        #
        if self.geo_st_flag:
            fe_grid = FEGrid(coord_min=(0.0, 0.0, 0.0),
                             coord_max=(1.0, 1.0, 1.0),
                             shape=(self.shape_xy, self.shape_xy, self.shape_z),
                             level=self.specmn_fe_level,
                             geo_transform=self.geo_st,
                             fets_eval=self.specmn_fets)
        else:
            # total effective length = 1.25m - 2*0.05m (from center supports to center support)
            width_supprt = 0.1
            length_quarter_eff = self.length / 2. - width_supprt / 2.
            length_quarter_eff = self.length / 2.
            fe_grid = FEGrid(coord_min=(0.0, 0.0, 0.0),
                             coord_max=(length_quarter_eff, length_quarter_eff, self.thickness),
                             shape=(self.shape_xy, self.shape_xy, self.shape_z),
                             level=self.specmn_fe_level,
                             fets_eval=self.specmn_fets)
        return fe_grid

    plain_concrete_fe_grid = Property(Instance(FEGrid), depends_on='+ps_levels, +input')
    @cached_property
    def _get_plain_concrete_fe_grid(self):
        fe_grid = FEGrid(coord_min=(0.0, 0.0, 0.0),
                         coord_max=(1.0, 1.0, 1.0),
                         shape=(self.shape_xy, self.shape_xy, 1),
                         level=self.plain_concrete_fe_level,
                         geo_transform=self.geo_st_plain_concrete,
                         fets_eval=self.plain_concrete_fets)
        return fe_grid

    elstmr_fe_level = Property(depends_on='+ps_levels, +input')
    @cached_property
    def _get_elstmr_fe_level(self):
        return  FERefinementGrid(name='elastomer patch',
                                 fets_eval=self.elstmr_fets,
                                 domain=self.fe_domain)

    elstmr_fe_grid = Property(Instance(FEGrid), depends_on='+ps_levels, +input')
    @cached_property
    def _get_elstmr_fe_grid(self):
        # if 'geo_st' is used for modeling of round load introduction area
        #
        if self.geo_st_flag:
            fe_grid = FEGrid(coord_min=(0.0, 0.0, 0.0),
                             coord_max=(1.0, 1.0, 1.0),
                             level=self.elstmr_fe_level,
                             shape=(self.shape_R, self.shape_R, 1),
                             geo_transform=self.geo_elstmr,
                             fets_eval=self.elstmr_fets)
        # if 'geo_st' is not used a regular rectangular grid is used as mesh for the slab
        #
        else:
            elstmr_min = self.length / 2. - self.elem_length
            elstmr_max = self.length / 2.
            fe_grid = FEGrid(coord_min=(elstmr_min, elstmr_min, self.thickness),
                             coord_max=(elstmr_max, elstmr_max, self.thickness + self.thickness_elstmr),
                             level=self.elstmr_fe_level,
                             shape=(self.shape_R, self.shape_R, 1),
                             fets_eval=self.elstmr_fets)
        return fe_grid

    if supprt_flag:
        supprt_fe_level = Property(depends_on='+ps_levels, +input')
        @cached_property
        def _get_supprt_fe_level(self):
            return  FERefinementGrid(name='support patch',
                                     fets_eval=self.supprt_fets,
                                     domain=self.fe_domain)

        supprt_fe_grid = Property(Instance(FEGrid), depends_on='+ps_levels, +input')
        @cached_property
        def _get_supprt_fe_grid(self):
            return FEGrid(coord_min=(0, 0, 0),
                          coord_max=(1, 1, 1),
                          level=self.supprt_fe_level,
                          # use shape (2,2) in order to place support in the center of the steel support
                          # corresponding to 4 elements of the slab mesh
                          #
                          shape=(self.shape_supprt_xy, self.shape_supprt_xy, 1),
                          geo_transform=self.geo_supprt,
                          fets_eval=self.supprt_fets)

    #--------------------------------------------------------------
    # tloop
    #--------------------------------------------------------------

    w_max = Float(-0.030, input=True)  # [m]

    tloop = Property(depends_on='input_change')
    @cached_property
    def _get_tloop(self):

        specmn = self.specmn_fe_grid
        print 'XXX specmn_fe_grid.n_dofs', specmn.n_dofs
        print 'XXX specmn_fe_grid.dof_grid', specmn.dof_grid.dofs.shape

        if self.elstmr_flag:
            elstmr = self.elstmr_fe_grid

        if self.supprt_flag:
            supprt = self.supprt_fe_grid

        if self.plain_concrete_flag:
            plain_concrete = self.plain_concrete_fe_grid

        self.fe_domain.n_dofs

        #--------------------------------------------------------------
        # boundary conditions for the symmetry planes and the support
        #--------------------------------------------------------------

        # symmetry planes of the specimen
        #
        bc_symplane_yz = BCSlice(var='u', value=0., dims=[0], slice=specmn[-1, :, :, -1, :, :])
        bc_symplane_xz = BCSlice(var='u', value=0., dims=[1], slice=specmn[:, -1, :, :, -1, :])

        if self.plain_concrete_flag:
            # symmetry planes of the plain concrete layer
            #
            bc_symplane_yz_plain_concrete = BCSlice(var='u', value=0., dims=[0], slice=plain_concrete[-1, :, :, -1, :, :])
            bc_symplane_xz_plain_concrete = BCSlice(var='u', value=0., dims=[1], slice=plain_concrete[:, -1, :, :, -1, :])

        # tappered support with fixed node in the corner
        #
#        if self.supprt_flag:
#            print 'support with edge rotation'
#            bc_supprt_symplane_yz = BCSlice(var = 'u', value = 0., dims = [0], slice = supprt[0, :, :, 0, :, :])
#            bc_supprt_symplane_xz = BCSlice(var = 'u', value = 0., dims = [1], slice = supprt[:, 0, :, :, 0, :])
#            support_000 = BCSlice(var = 'u', value = 0., dims = [2], slice = supprt[0, 0, 0, 0, 0, 0])
#            link_sp_sprt_z = BCSlice(var = 'u', value = 0., dims = [2],
#                                   slice = specmn[0, 0, 0, :, :, 0],
#                                   link_slice = supprt[0, 0, -1, :, :, -1],
#                                   link_dims = [2],
#                                   link_coeffs = [1.])

        # four sided tappered support with fixed node in the center
        # use two element rows for the size of the steel plate rotation support
        #
        if self.supprt_flag:
            #---------------------------
            # fixed support with sliding connection to the plate
            #---------------------------

            print 'four sided tappered support with center rotation'
            # allow movement of the support's outer corner node only in the 45-degree direction
            # by linking the corresponding degrees of freedom
            #
            sdofs = supprt[0, 0, 0, 0, 0, 0].dofs
            s_dof_x, s_dof_y = sdofs[0, 0, (0, 1)].flatten()
            sym_support_link = BCDof(var='u', value=0., dof=s_dof_x,
                                      link_dofs=[s_dof_y],
                                      link_coeffs=[1.0])

            if self.plain_concrete_flag:
                # link specimen bottom nodes with plain concrete layer for all directions
                #
                link_sp_pc_xyz = BCSlice(var='u', value=0., dims=[2],
                                       slice=specmn[:, :, 1, :, :, 1],
                                       link_slice=plain_concrete[:, :, -1, :, :, -1],
                                       link_dims=[0, 1, 2],
                                       link_coeffs=[1.])

                # link plain concrete layer bottom nodes in the slab corner with the top nodes of the tappered support in z-direction
                #
                link_pc_sprt_z = BCSlice(var='u', value=0., dims=[2],
                                       slice=plain_concrete[0:self.shape_supprt_xy, 0:self.shape_supprt_xy, 0, :, :, 0],
                                       link_slice=supprt[:, :, -1, :, :, -1],
                                       link_dims=[2],
                                       link_coeffs=[1.])

                # fix bottom-middle node of the tappered support
                #
                support_000 = BCSlice(var='u', value=0., dims=[0, 1, 2],
                                      slice=supprt[self.idx_supprt_center_elem, self.idx_supprt_center_elem, 0, 0, 0, 0])

                supprt_dof_z = supprt[self.idx_supprt_center_elem, self.idx_supprt_center_elem, 0, 0, 0, 0].dofs[0, 0, 2]
                print 'supprt_dof_z used for internal force tracing ', supprt_dof_z

            else:
                # link specimen bottom nodes in the slab corner with the top nodes of the tappered support in z-direction
                #
                link_sp_sprt_z = BCSlice(var='u', value=0., dims=[2],
                                       slice=specmn[0:self.shape_supprt_xy, 0:self.shape_supprt_xy, 0, :, :, 0],
                                       link_slice=supprt[:, :, -1, :, :, -1],
                                       link_dims=[2],
                                       link_coeffs=[1.])

                # fix bottom-middle node of the tappered support
                #
                support_000 = BCSlice(var='u', value=0., dims=[0, 1, 2],
                                      slice=supprt[self.idx_supprt_center_elem, self.idx_supprt_center_elem, 0, 0, 0, 0])

                supprt_dof_z = supprt[self.idx_supprt_center_elem, self.idx_supprt_center_elem, 0, 0, 0, 0].dofs[0, 0, 2]
                print 'supprt_dof_z used for internal force tracing ', supprt_dof_z

        else:
            # Var_1: place support at the corner / reduce the edge length of the slab in the model to effective values,
            # i.e. l_ = length_quarter - width_supprt / 2
            #
#            support_000 = BCSlice(var = 'u', value = 0., dims = [2], slice = specmn[0, 0, 0, 0, 0, 0])
#            supprt_dof_z = specmn[0, 0, 0, 0, 0, 0].dofs[0, 0, 2]
#            print 'supprt_dof_z used for internal force tracing ', supprt_dof_z

            # Var_2: place single support at the center of the support position (position defined by 'idx_supprt_center_elem'
            # depending on the mesh fines; NOTE: mesh fines and index must correspond approximately to the correct support width
            #
            support_000 = BCSlice(var='u', value=0., dims=[2], slice=specmn[self.idx_supprt_center_elem, self.idx_supprt_center_elem, 0, 0, 0, 0])
            supprt_dof_z = specmn[self.idx_supprt_center_elem, self.idx_supprt_center_elem, 0, 0, 0, 0].dofs[0, 0, 2]
            print 'supprt_dof_z used for internal force tracing ', supprt_dof_z


        #--------------------------------------------------------------
        # boundary conditions for loading
        #--------------------------------------------------------------

        # center displacement
        #
        w_max = self.w_max
#        w_max = Float( -0.020, auto_set = False, enter_set = True, input = True) # [m]

        # if elastomer is modeled
        #
        if self.elstmr_flag:
            bc_el_symplane_xz = BCSlice(var='u', value=0., dims=[1],
                                        slice=elstmr[:, -1, :, :, -1, :])
            bc_el_symplane_yz = BCSlice(var='u', value=0., dims=[0],
                                        slice=elstmr[-1, :, :, -1, :, :])
            if self.geo_st_flag:
                # link all elements of the load introduction area with the corresponding elements of the elastomer
                #
                link_el_sp = BCSlice(var='u', value=0., dims=[2],
                                     slice=elstmr[:, :, 0, :, :, 0],
                                     link_slice=specmn[-self.shape_R :, -self.shape_R :, -1, :, :, -1],
                                     link_dims=[2],
                                     link_coeffs=[1.])
            else:
                # link the corner element of the slab mesh with the elastomer (one element)
                #
                link_el_sp = BCSlice(var='u', value=0., dims=[2],
                                     slice=elstmr[0, 0, 0, :, :, 0],
                                     link_slice=specmn[-1, -1, -1, :, :, -1],
                                     link_dims=[2],
                                     link_coeffs=[1.])

            #--------------------------------------------------------------
            # var 'elstmr': nodal displacement applied at top of the load introduction plate
            #--------------------------------------------------------------
#            time_function = MFnLineArray(xdata = [0.0, 0.2, 0.4, 1.0], ydata = [0.0, 0.2, 0.75, 1.0])

            # displacement applied at the top of the elastomer
            #
            bc_w = [BCSlice(var='u', value=w_max, dims=[2],  # time_function = time_function.get_value,
                           slice=elstmr[:, :, -1, :, :, -1])]
            # dofs for RTrace f-w-diagram
            #
            load_dofs_z = elstmr[:, :, -1, :, :, -1].dofs[:, :, 2].flatten()
            print 'load_dofs_z', load_dofs_z
            print 'load_dofs_z.shape', load_dofs_z.shape

        # if no elastomer is modeled
        #
        else:

            # if geometry transformation for the slab is used
            #
            if self.geo_st_flag:
                # z-displacement applied along the edge of the load introduction plate
                #
                nodes_line_Rx = specmn[-self.shape_R:, -self.shape_R, -1, :, 0, -1]
                nodes_line_Ry = specmn[-self.shape_R, -self.shape_R:, -1, 0, :, -1]
                bc_w_Rxline = BCSlice(var='u', value=w_max, dims=[2], slice=nodes_line_Rx)
                bc_w_Ryline = BCSlice(var='u', value=w_max, dims=[2], slice=nodes_line_Ry)
                bc_w = [ bc_w_Rxline, bc_w_Ryline ]
                # dofs for RTrace f-w-diagram
                #
                dofs_line_Rx = nodes_line_Rx.dofs[:, :, 2].flatten()
                dofs_line_Ry = nodes_line_Ry.dofs[:, :, 2].flatten()
                # use 'np.unique' in order to avoid that the same 'dof' is listed twice in the
                # array as it appears in two adjacent elements of the same slice or in both x- and y-line slice
                load_dofs_z = np.unique(np.hstack([dofs_line_Rx, dofs_line_Ry]))
                print 'dofs_line_Rx', dofs_line_Rx
                print 'dofs_line_Ry', dofs_line_Ry
                print 'load_dofs_z', load_dofs_z

            # if no geometry transformation for the slab is used (regular rectangular mesh only)
            #
            else:
                # ## var 1: single nodal displacement applied at the center top node (one node only)
                #
                bc_w = [BCSlice(var='u', value=w_max, dims=[2],
                                          slice=specmn[-1, -1, -1, -1, -1, -1])]  # only top node
                load_dofs_z = specmn[-1, -1, -1, -1, -1, -1].dofs[0, 0, 2]
                print 'load_dofs_z', load_dofs_z
                print 'load_dofs_z.shape', load_dofs_z.shape
                load_dofs_z = array([specmn[-1, -1, -1, -1, -1, -1].dofs[0, 0, 2]])
                print 'load_dofs_z', load_dofs_z
                print 'load_dofs_z.shape', load_dofs_z.shape

#                ### var 2 displacement applied at all center nodes (all nodes in the symmetry axis)
#                #
#                bc_w = [BCSlice(var = 'u', value = w_max, dims = [2],
#                                          slice = specmn[-1, -1, :, -1, -1, :])] # all nodes along z-axis
#                load_dofs_z = specmn[-1, -1, :, -1, -1, :].dofs[0, 0, 2]
#                print 'load_dofs_z', load_dofs_z

                # ## var 3:
                #
                # displacement applied along the edge of the 2nd element row (from the middle) for
                # shape_xy = 10 and an edge length of L_quarter = 0.575m the distance from the center equals 0.23m
#                nodes_line_x = specmn[-2:, -2, -1, :, 0, -1]
#                nodes_line_y = specmn[-2, -2:, -1, 0, :, -1]
#                bc_w_xline = BCSlice(var = 'u', value = w_max, dims = [2], slice = nodes_line_x)
#                bc_w_yline = BCSlice(var = 'u', value = w_max, dims = [2], slice = nodes_line_y)
#                bc_w = [ bc_w_xline, bc_w_yline ]
#                # dofs for RTrace f-w-diagram
#                #
#                dofs_line_x = nodes_line_x.dofs[:, :, 2].flatten()
#                dofs_line_y = nodes_line_y.dofs[:, :, 2].flatten()
#                load_dofs_z = np.unique( np.hstack([dofs_line_x, dofs_line_y]) )

        #--------------------------------------------------------------
        # force-displacement-diagram
        #--------------------------------------------------------------

        center_top_dof_z = specmn[-1, -1, -1, -1, -1, -1].dofs[0, 0, 2]
        print 'center_top_dof used for displacement tracing: ', center_top_dof_z

        self.f_w_diagram_center = RTraceGraph(name='displacement (center) - force',
                                       var_x='U_k'  , idx_x=center_top_dof_z,
                                       # elastomer load
                                       var_y='F_int', idx_y_arr=load_dofs_z,
 #                                       var_y = 'F_int', idx_y = center_dof,
                                       record_on='update',
                                       transform_x='-x * 1000',  # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the slab (2 symmetry axis):
                                       #
                                       transform_y='-4 * 1000. * y')

        # trace support reaction
        #
        self.f_w_diagram_supprt = RTraceGraph(name='displacement (supprt) - force',
                                       var_x='U_k'  , idx_x=center_top_dof_z,
                                       # elastomer load
                                       var_y='F_int', idx_y=supprt_dof_z,
                                       record_on='update',
                                       transform_x='-x * 1000',  # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the slab (2 symmetry axis):
                                       #
                                       transform_y='4 * 1000. * y')

        #------------------------
        # F-w-diagram edge and edge/center
        #------------------------
#        # get center-edge dof at the upper side of the slab
#        #
#        # if shape_xy is an even number:
#        if self.shape_xy % 2 == 0:
#            x_idx_center_edge = (self.shape_xy / 2) - 1
#            center_edge_top_dofs = specmn[x_idx_center_edge, -1, -1, -1, -1, -1].dofs
#        # if shape_xy is an odd number:
#        else:
#            # get the midside node of the center-edge-element
#            x_idx_center_edge = (self.shape_xy - 1) / 2
#            # valid only for quadratic elements
#            center_edge_top_dofs = specmn[x_idx_center_edge, -1, -1, 1, -1, -1].dofs

#        center_edge_dof = center_edge_top_dofs[0, 0, 2]
#        print 'center_edge_dof', center_edge_dof
#        self.f_w_diagram_center_edge = RTraceGraph(name = 'displacement (center/edge) - reaction 2',
#                                                   var_x = 'U_k'  , idx_x = center_edge_dof,
#                                                   var_y = 'F_int', idx_y_arr = load_dofs_z,
#                                                   record_on = 'update',
#                                                   transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
#                                                   transform_y = '-4 * 1000 * y')
#
#        edge_dof = specmn[ 0, -1, -1, 0, -1, -1].dofs[0, 0, 2]
#        print 'edge_dof', edge_dof
#        self.f_w_diagram_edge = RTraceGraph(name = 'displacement (edge) - reaction 2',
#                                            var_x = 'U_k'  , idx_x = edge_dof,
#                                            var_y = 'F_int', idx_y_arr = load_dofs_z,
#                                            record_on = 'update',
#                                            transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
#                                            transform_y = '-4 * 1000 * y')

        #--------------------------------------------------------------
        # boundary conditions
        #--------------------------------------------------------------

        bcond_list = [ bc_symplane_yz, bc_symplane_xz, support_000 ] + bc_w

        if self.supprt_flag:
            if self.plain_concrete_flag:
                bcond_list += [ link_sp_pc_xyz, link_pc_sprt_z, sym_support_link, bc_symplane_yz_plain_concrete, bc_symplane_xz_plain_concrete ]
            else:
                bcond_list += [ link_sp_sprt_z, sym_support_link ]

        if self.elstmr_flag:
            bcond_list += [ bc_el_symplane_yz, bc_el_symplane_xz, link_el_sp ]

        print 'bcond_list', bcond_list
        print 'bcond_list', len(bcond_list)
        print 'bcond_list', np.product(bcond_list[1].dofs.shape)
        print 'bcond_list', [np.product(bcond_list[0].dofs.shape),
                             np.product(bcond_list[1].dofs.shape),
                             np.product(bcond_list[2].dofs.shape),
                             np.product(bcond_list[3].dofs.shape),
                             np.product(bcond_list[4].dofs.shape),
                             np.product(bcond_list[5].dofs.shape)]

        #--------------------------------------------------------------
        # ts
        #--------------------------------------------------------------

        ts = TS(
                sdomain=self.fe_domain,
                bcond_list=bcond_list,
                rtrace_list=[
                              # @todo: check weather tracing yields same results:
                              #
                              self.f_w_diagram_center,
                              self.f_w_diagram_supprt,
#                             self.f_w_diagram_center_edge,
#                             self.f_w_diagram_edge,

                              RTraceDomainListField(name='Displacement' ,
                                                    var='u', idx=0, warp=True),
#                             RTraceDomainListField(name = 'Stress' ,
#                                            var = 'sig_app', idx = 0, warp = True,
#                                            record_on = 'update'),

#                             RTraceDomainListField(name = 'Strain' ,
#                                        var = 'eps_app', idx = 0, warp = True,
#                                        record_on = 'update'),
                             RTraceDomainListField(name='Damage' ,
                                        var='omega_mtx', idx=0, warp=True,
                                        record_on='update'),

                             # max[omega_i for i = 1..Nmp]
                             RTraceDomainListField(name='max_omega_i', warp=True,
                                        var='max_omega_i', idx=0,
                                        record_on='update'),
#                             RTraceDomainListField(name='max_omega_i2', warp=True,
#                                        var='max_omega_i2', idx=0,
#                                        record_on='update'),

                             # value in phi_mtx for 1st principle damage direction
                             RTraceDomainListField(name='phi_pdc', warp=True,
                                        var='phi_pdc', idx=0,
                                        record_on='update'),
#                             RTraceDomainListField(name='phi_pdc2', warp=True,
#                                        var='phi_pdc2', idx=0,
#                                        record_on='update'),


#                             RTraceDomainListField(name = 'IStress' ,
#                                            position = 'int_pnts',
#                                            var = 'sig_app', idx = 0,
#                                            record_on = 'update'),
#                             RTraceDomainListField(name = 'IStrain' ,
#                                            position = 'int_pnts',
#                                            var = 'eps_app', idx = 0,
#                                            record_on = 'update'),
                              ]
                )

        # Add the time-loop control
        #
        tloop = TLoop(tstepper=ts,
#                      KMAX = 300,
                      KMAX=100,
                      RESETMAX=0,
                      tolerance=self.tolerance,
                      tline=TLine(min=0.0, step=self.tstep, max=self.tmax),
                      ord=np.inf
                      )

        return tloop

    #--------------------------------------------------------------
    # prepare pstudy
    #--------------------------------------------------------------

    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()

        self.f_w_diagram_center.refresh()
        F_max = max(self.f_w_diagram_center.trace.ydata)

        u_center_top_z = U[ self.center_top_dofs ][0, 0, 2]
        return array([ u_center_top_z, F_max ],
                       dtype='float_')

    def get_sim_outputs(self):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut(name='u_center_top_z', unit='m'),
                 SimOut(name='F_max', unit='kN') ]


class SimSTDB(SimST):

    # extend the failure strain in 'PhiFnGeneralExtended' by this factor
    # before the integrety function drops down to zero
    #
    factor_eps_fail = Float(1.0, input=True,
                             ps_levels=(1.0, 1.2, 3))

    # specify exponential softening in 'PhiFnGeneralExtendedExp'
    #
    Efp_frac = Float(0.2, input=True)


    #-----------------
    # composite cross section unit cell:
    #-----------------
    #
    ccs_unit_cell_key = Enum('FIL-10-09_2D-05-11_0.00462_all0',
                              CCSUnitCell.db.keys(),
                              simdb=True, input=True,
                              auto_set=False, enter_set=True)

    ccs_unit_cell_ref = Property(Instance(SimDBClass),
                                 depends_on='ccs_unit_cell_key')
    @cached_property
    def _get_ccs_unit_cell_ref(self):
        return CCSUnitCell.db[ self.ccs_unit_cell_key ]

    #-----------------
    # damage function:
    #-----------------
    #
    material_model = Str(input=True)
    def _material_model_default(self):
        # return the material model key of the first DamageFunctionEntry
        # This is necessary to avoid an ValueError at setup
        return self.ccs_unit_cell_ref.damage_function_list[0].material_model

    calibration_test = Str(input=True)
    def _calibration_test_default(self):
        # return the material model key of the first DamageFunctionEntry
        # This is necessary to avoid an ValueError at setup
        return self.ccs_unit_cell_ref.damage_function_list[0].calibration_test

    damage_function = Property(Instance(MFnLineArray),
                                depends_on='input_change')
    @cached_property
    def _get_damage_function(self):
#        print'XXX', [self.ccs_unit_cell_ref.damage_function_list[i].calibration_test for i in range(14)]
        return self.ccs_unit_cell_ref.get_param(self.material_model, self.calibration_test)

    #-----------------
    # phi function General:
    #-----------------
    #
    # definition of 'phi_fn' in order to run scripting
    # @todo: how can the implementation be changed in order to allow for parameter passing, e.g. 'factor_eps_fail', etc.
    #
    phi_fn_class = Enum(PhiFnGeneral, [PhiFnGeneral, PhiFnGeneralExtended, PhiFnGeneralExtendedExp], input=True,)
    phi_fn = Property(Instance(phi_fn_class), depends_on='phi_fn_class, input_change, +ps_levels')
    @cached_property
    def _get_phi_fn(self):
        return self.phi_fn_class(mfn=self.damage_function)

    # @todo: is mapped traist a better alternative?
#    phi_fn_class = Trait('PhiFnGeneral', { 'PhiFnGeneral' : PhiFnGeneral( mfn = self.damage_function ),
#                                           'PhiFnGeneralExtended : PhiFnGeneralExtended( mfn = self.damage_function, factor_eps_fail = self.factor_eps_fail ),
#                                           'PhiFnGeneralExtendedExp : PhiFnGeneralExtendedExp( mfn = self.damage_function, Efp_frac = self.Efp_frac ) },
#                                           input = True,)

    #-----------------
    # phi function General:
    #-----------------
#    phi_fn = Property(Instance(PhiFnGeneral),
#                       depends_on = 'input_change,+ps_levels')
#    @cached_property
#    def _get_phi_fn(self):
#        return PhiFnGeneral( mfn = self.damage_function )
#
#    #-----------------
#    # phi function General Extended:
#    #-----------------
#    #
#    phi_fn_ext = Property(Instance(PhiFnGeneralExtended),
#                       depends_on = 'input_change,+ps_levels')
#    @cached_property
#    def _get_phi_fn_ext(self):
#        return PhiFnGeneralExtended( mfn = self.damage_function,
#                                     factor_eps_fail = self.factor_eps_fail )
#
#    #-----------------
#    # phi function General Extended Exp:
#    #-----------------
#    #
#    phi_fn_exp = Property(Instance(PhiFnGeneralExtendedExp),
#                       depends_on = 'input_change,+ps_levels')
#    @cached_property
#    def _get_phi_fn_exp(self):
#        return PhiFnGeneralExtendedExp(mfn = self.damage_function,
#                                       Efp_frac = self.Efp_frac)

    #----------------------------------------------------------------------------------
    # material properties stored in 'ccs'
    #----------------------------------------------------------------------------------

    # concrete matrix E-modulus (taken from 'ccs_unit_cell')
    #
    E_m = Property(Float, depends_on='input_change')
    @cached_property
    def _get_E_m(self):
        E_m = self.ccs_unit_cell_ref.get_E_m_time(self.age)
        print 'E_m (from ccs)', E_m
#        E_m = 29100.
#        print 'E_m set explicitly to ', E_m
        return E_m

    # composite E-modulus (taken from 'ccs_unit_cell')
    #
    E_c = Property(Float, depends_on='input_change')
    @cached_property
    def _get_E_c(self):
        E_c = self.ccs_unit_cell_ref.get_E_c_time(self.age)
        print 'E_c (from ccs)', E_c
        E_c = 18709.5
        print 'E_c set explicitly to 18709.5'
#        print 'E_c set explicitly to 29100.'
#        E_c = 29100.
        return E_c

    # Poisson's ratio
    #
    nu = Property(Float, depends_on='input_change')
    @cached_property
    def _get_nu(self):
        nu = self.ccs_unit_cell_ref.nu
        print 'nu (from ccs)', nu
        # set nu explicitly corresponding to settings in 'mats_calib_damage_fn'
        #
        print 'nu set explicitly to 0.20'
        nu = 0.2
        return nu


if __name__ == '__main__':

# paramters for AG-glas slab tests (3cm) with tricot-binding:
#                                ccs_unit_cell_key = 'FIL-10-09_2D-02-06a_0.00273_90_0',
#                                calibration_test = 'TT11-10a-average',
#                                calibration_test = 'TT11-10a-V2',
#                                age = 28 )

    # s_tex,z = 4.62 mm corresponds to: n_tex = 12; h = 60 mm
    #
    sim_model = SimSTDB(

#                        ccs_unit_cell_key = 'FIL-10-09_2D-02-06a_0.00273_90_0',
#                        calibration_test = 'TT11-10a-average',

                        ccs_unit_cell_key='barrelshell_2D-05-11_0.00286_all0',
                        calibration_test='TTb-6c-2cm-0-TU-V1_bs5_a23d-nu02',
                        thickness=0.02,
                        length=0.80,

#                        ccs_unit_cell_key = 'FIL-10-09_2D-05-11_0.00462_all0',

                         # calibration_test = 'TT-12c-6cm-TU-SH1F-V1',
#                         calibration_test = 'TT-12c-6cm-0-TU-SH2F-V2',

                         # calibration for: age = 26d; E_m = 28400 MPa; nu = 0.25
                         #
#                         calibration_test = 'TT-12c-6cm-0-TU-SH2F-V3',

                         # calibration for: age = 23d; E_m = 28400 MPa; nu = 0.25
                         #
#                         calibration_test = 'TT-12c-6cm-0-TU-SH2F-V3_a23d-nu02',
                         # age of the slab at the time of testing
                         # NOTE: that the same phi-function is used independent of age. This assumes a
                         # an afine/proportional damage evolution for different ages. get_E_m_time(self.age)
                         #
                         age=23,
#                         shape_xy = 10,
#                         shape_z = 3,
                         tstep=0.05,
                         tmax=0.4
#                         tstep = 0.05
#                         tstep = 0.025 # tmax / 34 steps with tmax = 0.85 corresponding to 1 mm / tstep
                         )

    #-----------------------------------------
    # ## Var.3 ### ST-12c-6cm; L = 1,25m; t = 60 mm
    #-----------------------------------------
    #
    sim_model = SimSTDB(

                        thickness=0.06,
                        length=1.25,

                        ccs_unit_cell_key='FIL-10-09_2D-05-11_0.00462_all0',

                         # calibration for: age = 23d; E_m = 28400 MPa; nu = 0.25
                         #
                         calibration_test='TT-12c-6cm-0-TU-SH2F-V3_a23d-nu02',
                         # age of the slab at the time of testing
                         age=23,
                         # NOTE: that the same phi-function is used independent of age. This assumes a
                         # an afine/proportional damage evolution for different ages.
                         #
                         shape_xy=28,
                         shape_R=5,
                         shape_z=3,
                         #
                         tstep=0.05,
                         tmax=1.0,
                         tolerance=0.0001
                         )

    #------------------------------
    # do
    #------------------------------
#    do = 'show_phi_fn'
    do = 'ui'
#    do = 'pstudy'
#    do = 'validation'
#    do = 'show_last_results'

    #------------------------------
    # ui
    #------------------------------
    if do == 'ui':
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource=sim_model)
#        sim_model.tloop.eval()
        app.main()

    #------------------------------
    # validation
    #------------------------------
    if do == 'validation':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = 'pickle_files'
        png_path = 'png_files'

        p.figure(facecolor='white')  # white background

        #--------------------
        # simulation
        #--------------------
        sim_model.tloop.eval()

        pickle_path = 'pickle_files'
        ccs_unit_cell_key = sim_model.ccs_unit_cell_key
        calibration_test = sim_model.calibration_test
        length = sim_model.length
        thickness = sim_model.thickness
        shape_xy = sim_model.shape_xy
        E_m = sim_model.E_m
        nu = sim_model.nu
        tolerance = sim_model.tolerance
        param_key = ccs_unit_cell_key + '_' + calibration_test + 'L=%g h=% sxy= %g Em=%g nu=%g tol=%g' % (length, thickness, shape_xy, E_m, nu, tolerance)

        # f-w-diagram_center
        #
        sim_model.f_w_diagram_center.refresh()
        file_name = 'f_w_diagram_c_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'w')
        dump(sim_model.f_w_diagram_center.trace, file)
        file.close()
        sim_model.f_w_diagram_center.trace.mpl_plot(p, color='red')

        # f-w-diagram_supprt
        #
        sim_model.f_w_diagram_supprt.refresh()
        file_name = 'f_w_diagram_supprt_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'w')
        dump(sim_model.f_w_diagram_supprt.trace, file)
        file.close()
        sim_model.f_w_diagram_supprt.trace.mpl_plot(p, color='blue')

#        # f-w-diagram_center-edge
#        #
#        sim_model.f_w_diagram_center.refresh()
#        file_name = 'f_w_diagram_ce_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'w')
#        dump(sim_model.f_w_diagram_center_edge.trace, file)
#        file.close()
#        sim_model.f_w_diagram_center_edge.trace.mpl_plot(p, color = 'red')
#
#        # f-w-diagram_edge
#        #
#        sim_model.f_w_diagram_center.refresh()
#        file_name = 'f_w_diagram_e_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'w')
#        dump(sim_model.f_w_diagram_edge.trace, file)
#        file.close()
#        sim_model.f_w_diagram_edge.trace.mpl_plot(p, color = 'red')

        #--------------------
        # experiments
        #--------------------

        # PT-12c-6cm-TU
#        path = join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU')
#        tests = [ 'ST-12c-6cm-u-TU.DAT' ]

        # PT-10a
#        path10 = join( simdb.exdata_dir, 'slab_tests', '2010-03-08_ST-10g-3cm-a-FR_TRC10', 'ST-10g-3cm-a-FR-TRC10.DAT' )
#        path11 = join( simdb.exdata_dir, 'slab_tests', '2010-03-09_ST-10g-3cm-a-FR_TRC11', 'ST-10g-3cm-a-FR-TRC11.DAT')
#        path12 = join( simdb.exdata_dir, 'slab_tests', '2010-03-10_ST-10g-3cm-a-FR_TRC12', 'ST-10g-3cm-a-FR-TRC12.DAT' )
#        tests = [path10, path11, path12]
#        for ex_path in tests:
#            ex_run = ExRun(ex_path)
#            ex_run.ex_type._plot_force_center_deflection_smoothed(p)

        # ST-6c-2cm-TU_bs2
        #
        ex_path = join(simdb.exdata_dir, 'slab_tests', '2013-07-10_ST-6c-2cm-TU_bs2', 'ST-6c-2cm-TU_bs2.DAT')
        ex_run = ExRun(ex_path)
        ex_run.ex_type._plot_force_center_deflection(p)

        # plot sim curve as time new roman within the predefined limits
        #
#        format_plot(p, xlim = 34, ylim = 54, xlabel = 'displacement [mm]', ylabel = 'force [kN]')
        png_file_path = join(png_path, param_key + '.png')
        p.title(param_key)
        p.savefig(png_file_path, dpi=600.)
        p.show()

    #------------------------------
    # show last results
    #------------------------------

    if do == 'show_last_results':
        from matresdev.db.exdb.ex_run import ExRun
        import pylab as p

        #------------------
        # simulation
        #------------------

        pickle_path = 'pickle_files'
#        param_key = 'SH2F-V3_nelems14-3-2_w_line_Ec28400'
#        param_key = 'SH2F-V3_nelems14-3-2_line_w_Em23d_a23d-nu02_ts0-0025'
#        param_key = 'SH2F-V3_nelems10-10-2_node_w_node_supprt'
        param_key = 'FIL-10-09_2D-02-06a_0.00273_90_0__TT11-10a-average'

        # f-w-diagram_supprt
        #
        file_name = 'f_w_diagram_supprt_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color='blue')

        # f-w-diagram_center
        #
        file_name = 'f_w_diagram_c_' + param_key + '.pickle'
        pickle_file_path = join(pickle_path, file_name)
        file = open(pickle_file_path, 'r')
        trace = load(file)
        p.plot(trace.xdata, trace.ydata, color='red')

#        # f-w-diagram_center-edge
#        #
#        file_name = 'f_w_diagram_ce_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'r')
#        trace = load(file)
#        p.plot(trace.xdata, trace.ydata, color = 'red')
#
#        # f-w-diagram_edge
#        #
#        file_name = 'f_w_diagram_e_' + param_key + '.pickle'
#        pickle_file_path = join(pickle_path, file_name)
#        file = open(pickle_file_path, 'r')
#        trace = load(file)
#        p.plot(trace.xdata, trace.ydata, color = 'red')

        #------------------
        # experiments
        #------------------

#        path = join(simdb.exdata_dir, 'slab_tests', '2011-12-15_ST-12c-6cm-u-TU')
#        tests = [ 'ST-12c-6cm-u-TU.DAT' ]
# #        path = join( simdb.exdata_dir, 'plate_tests', 'PT-10a' )
# #        tests = [ 'PT10-10a.DAT', 'PT11-10a.DAT' , 'PT12-10a.DAT' ]
# #        tests = [ 'PT10-10a.DAT' ]
#
#        for t in tests:
#            ex_path = join(path, t)
#            ex_run = ExRun(ex_path)
#            ex_run.ex_type._plot_force_deflection_avg_interpolated(p)
# #            ex_run.ex_type._plot_force_deflection_avg( p )

        # ST-6c-2cm-TU_bs2
        #
        ex_path = join(simdb.exdata_dir, 'slab_tests', '2013-07-10_ST-6c-2cm-TU_bs2', 'ST-6c-2cm-TU_bs2.DAT')
        ex_run = ExRun(ex_path)
        ex_run.ex_type._plot_force_center_deflection(p)
#            ex_run.ex_type._plot_force_center_deflection_smoothed(p)

        # plot sim curve as time new roman within the predefined limits
        #
#        format_plot(p, xlim = 34, ylim = 54, xlabel = 'displacement [mm]', ylabel = 'force [kN]')

        p.show()


    if do == 'pstudy':
        sim_ps = SimPStudy(sim_model=sim_model)
        sim_ps.configure_traits()








        # ## shape ###:

        # plate elem(12,12,4); w=0.08; step = 0.05; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-05_KMAX20_tol0-001_nelems12-12-4.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'blue' )

#        # plate elem(8,8,4); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_nelems8-8-4.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )

#        # plate elem(10,10,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_nelems10-10-2.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )


        # ## eps_factor ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001, factor_eps_fail = 1.4
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_epsfactor_1-4.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001, factor_eps_fail = 1.2
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_epsfactor_1-2.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open(pickle_file_path, 'r')
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )


        # ## average - V1 (c) ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )


        # ## c ce e - average ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_ce_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_e_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )


        # ## c ce e - V1 ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_ce_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_e_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )


        # ## step ###:

#        # plate elem(8,8,2); w=0.08; step = 0.02; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_V1_step0-02_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )


        # ## tolerance ###:

#        # plate elem(8,8,2); w=0.08; step = 0.02; KMAX = 50, tol = 0.0005
#        file_name = 'f_w_diagram_V1_step0-02_KMAX50_tol0-0005.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.02; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_V1_step0-02_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'yellow' )


#    print 'sim_model.ccs_unit_cell_key', sim_model.ccs_unit_cell_key
#    print 'sim_model.damage_function', sim_model.damage_function

#    import pylab as p
#    mfn = sim_model.damage_function
#    eps_last = sim_model.damage_function.xdata[-1]
#    print 'eps_last', eps_last
#    mfn.mpl_plot(p)
#    sim_model.phi_fn.mfn.mpl_plot(p)
#    phi_fn_vect = np.vectorize( sim_model.phi_fn.get_value )
#    xdata = np.linspace(0., eps_last*1.1, 400)
#    ydata = phi_fn_vect( xdata )
#    p.plot(xdata, ydata)
#    p.show()
