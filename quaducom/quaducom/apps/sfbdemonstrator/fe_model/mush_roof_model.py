'''
MUSHROOF Basis Class

TODO:  @ Andreas
        - depends_on, to many calculations of geometric transformation
        - update of internal strain bool is not possible default value always set in advance
'''

from etsproxy.traits.api import \
    Array, implements, Float, Property, cached_property, Instance, \
    Int, List, Bool, Dict, Enum

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel

#material model
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

#elements
from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U
from ibvpy.mesh.fe_grid import \
    FEGrid

from numpy import \
    diag, size, arange, append, sum, argsort, zeros_like, array, shape

from matplotlib.pyplot import \
    bar, show, axhline, ion, ioff, xlabel, ylabel, title, figure, savefig, ylim

import csv

from simiter.sim_pstudy import\
    ISimModel, SimOut, SimPStudy, SimArray, SimArrayView

# geo transform
#
from geo_column import GEOColumn
from hp_shell import HPShell
from hp_shell_conn_detail import HPShell


class MushRoofModel(IBVModel):
    '''Basis Class for Mushroof models (mr_quarter, mr_one, mr_one_free, mr_two, mr_four)
    '''

    #===========================================================================
    # initial strain
    #===========================================================================
    initial_strain_roof = False
    initial_strain_col = False

    alpha = Float(1.3e-5)
    t_up = Float(-100.)
    t_lo = Float(-100.)
    def temperature_strain_z(self, X_pnt, x_pnt):
        alpha, t_up, t_lo = self.alpha, self.t_up, self.t_lo
        delta_t = t_lo + (t_up - t_lo) * x_pnt[2]
        epsilon_0 = alpha * delta_t
        # return the initial volumetric strain tensor with n_dims  
        return diag([ epsilon_0 for i in range(3) ])

    #===========================================================================
    # material model
    #===========================================================================

    E_roof = Float(28700) # [MN/m^2]
    E_column = Float(32800) # [MN/m^2] E_cm for concrete C45/55 
    E_plate = Float(210000) # [MN/m^2] steel plate
    nu = Float(0.2) # [-]

    mats_roof = Property(Instance(MATS3DElastic), depends_on='+input')
    @cached_property
    def _get_mats_roof(self):
        if self.initial_strain_roof == True:
            return MATS3DElastic(E=self.E_roof, nu=self.nu, initial_strain=self.temperature_strain_z)
        else:
            return MATS3DElastic(E=self.E_roof, nu=self.nu)

    mats_column = Property(Instance(MATS3DElastic), depends_on='+input')
    @cached_property
    def _get_mats_column(self):
        if self.initial_strain_col == True:
            return MATS3DElastic(E=self.E_column, nu=self.nu, initial_strain=self.temperature_strain_z)
        else:
            return MATS3DElastic(E=self.E_column, nu=self.nu)

    mats_plate = Property(Instance(MATS3DElastic), depends_on='+input')
    @cached_property
    def _get_mats_plate(self):
        return MATS3DElastic(E=self.E_plate, nu=self.nu)

    #===========================================================================
    # finite elements
    #===========================================================================

    fe_linear_roof = Property(Instance(FETSEval, transient=True), depends_on='+input')
    def _get_fe_linear_roof(self):
        return FETS3D8H(mats_eval=self.mats_roof)

    fe_quad_serendipity_roof = Property(Instance(FETSEval, transient=True), depends_on='+input')
    def _get_fe_quad_serendipity_roof(self):
        return FETS3D8H20U(mats_eval=self.mats_roof)

    fe_quad_serendipity_column = Property(Instance(FETSEval, transient=True), depends_on='+input')
    def _get_fe_quad_serendipity_column(self):
        return FETS3D8H20U(mats_eval=self.mats_column)

    fe_linear_plate = Property(Instance(FETSEval, transient=True), depends_on='+input')
    def _get_fe_linear_plate(self):
        return FETS3D8H(mats_eval=self.mats_plate)

    fe_quad_serendipity_plate = Property(Instance(FETSEval, transient=True), depends_on='+input')
    def _get_fe_quad_serendipity_plate(self):
        return FETS3D8H20U(mats_eval=self.mats_plate)

    #===========================================================================
    # geometric dimensions
    #===========================================================================

    # dimensions of one quarter of the shell structure [m]
    #
    length_xy_quarter = Float(3.5, input=True)#, ps_levels = [4, 16, 5] )
    length_z = Float(0.927, input=True)#, ps_levels = [1, 2, 1] )

    # shell thickness is used only by option 'const_reinf_layer_elem'
    #
    t_shell = Float(0.06, input=True)

    # dimensions of the steel plate
    #
    t_plate = Float(0.03, unit='m', input=True)

    # dimensions of the column
    # NOTE: width of column at top is used also by option 'shift_elem' of 'HPShell'
    #
    width_top_col = Float(0.45, unit='m', input=True)
    width_bottom_col = Float(0.35, unit='m', input=True)

    # column length (from lowest point of the shell to the upper edge of the foundation:
    #
    h_col = Float(3.60, unit='m', input=True)

#    r_pipe = Float( 0.1, unit = 'm', input = True )

    scalefactor_delta_h = Float(1.00, input=True) # [-]
    scalefactor_length_xy = Float(1.00, input=True) # [-] 


    #-----------------------------------------------------------------
    # specify the relation of the total structure (in the 'mushroof'-model)
    # with respect to a quarter of one shell defined in 'HPShell'
    #-----------------------------------------------------------------

    # choose model and discretization for roof shell
    #
    mushroof_part = Enum('one', 'quarter', 'four', input=True)

    def _mushroof_part_default(self):
        return 'one'

    # @todo: add comment!
    # this is no input!
    #
    X0 = List([0, 0, 0], input=True)

    #-----------------------------------------------------------------
    # discretization
    #-----------------------------------------------------------------

    # shell discretization:
    #
    n_elems_xy_quarter = Int(10, input=True)
    n_elems_z = Int(1, input=True)

    # @todo: remove "+input", use mapped traits instead!
    #
    n_elems_xy_dict = Property(Dict, depends_on='+ps_levels, +input')
    def _get_n_elems_xy_dict(self):
        return  {'quarter': self.n_elems_xy_quarter,
                 'one'    : self.n_elems_xy_quarter * 2,
                 #@todo: include "scale_size" parameter used by HPShell!
                 'detail' : self.n_elems_xy_quarter * 2  }

    n_elems_xy = Property(Int , depends_on='+ps_levels, +input')
    def _get_n_elems_xy(self):
        # @todo: mushroff_part == "four" is not supported! (see HPShell)
        return int(self.n_elems_xy_dict[self.mushroof_part])

    # column discretization:
    #
    n_elems_col_xy = Int(2, input=True)#, ps_levels = [5, 20, 3 ] )

    #-----------------------------------------------------------------
    # option 'shift_elem'
    #-----------------------------------------------------------------

    # optional parameters
    #
    shift_elems = Bool(True, input=True)

    # set fixed points for shell discretization
    # NOTE: method is overwritten in subclass, e.g. 'MRtwo'
    #
    shift_array = Array
#    def _get_shift_array( self ):
#        return array( [[self.width_top_col / 2 ** 0.5,
#                        self.width_top_col / 2 ** 0.5,
#                        self.n_elems_col_xy / 2], ] )


    #-----------------------------------------------------------------
    # option 'const_reinf_layer_elem'
    #-----------------------------------------------------------------

    # element thickness defined (e.g to 3 cm) for bottom and top layer of the roof
    # needed to simulate reinforced area of the shell for non-linear simulation
    # @todo: remove "+input"
    #
    const_reinf_layer_elem = Bool (False, input=True)


    #-----------------------------------------------------------------
    # geometric transformations
    #-----------------------------------------------------------------

    #@todo: remove! geo_transforme is performed in the subclass 'MRTwo'


    #===========================================================================
    # evaluation
    #===========================================================================

    tline = Instance(TLine)
    def _tline_default(self):
        return TLine(min=0.0, step=1.0, max=1.0)

    max_princ_stress = Instance(RTraceDomainListField)
    def _max_princ_stress_default(self):
        return RTraceDomainListField(name='max principle stress' , idx=0,
                                      var='max_principle_sig', warp=True,
#                                      position = 'int_pnts',
                                      record_on='update',)

    dof = Int(1)
    f_dof = Property (Instance(RTraceGraph), depends_on='+ps_levels, +input')
    def _get_f_dof(self):
        return RTraceGraph(name='Fi,right over u_right (iteration)' ,
                           var_y='F_int', idx_y=self.dof,
                           var_x='U_k', idx_x=self.dof + 2,
                           record_on='update')

    sig_app = Property(Instance(RTraceDomainListField), depends_on='+ps_levels, +input')
    @cached_property
    def _get_sig_app(self):
        return RTraceDomainListField(name='sig_app' ,
#                                      position = 'int_pnts',
                                      var='sig_app',
                                      record_on='update',)

    u = Property(Instance(RTraceDomainListField), depends_on='+ps_levels, +input')
    @cached_property
    def _get_u(self):
        return RTraceDomainListField(name='displacement' ,
                                      var='u', warp=True,
                                      record_on='update',)

    damage = Property(Instance(RTraceDomainListField), depends_on='+ps_levels, +input')
    @cached_property
    def _get_damage(self):
        return RTraceDomainListField(name='Damage' ,
                                      var='omega_mtx', idx=0, warp=True,
                                      record_on='update')

    # sorting force of slice by number of dofs internal forces
    #
    def sort_by_dofs(self, dofs, unsorted):
        """
        for dof slices, the slice for one edge, results in a more
        dimensional array of form (a,b,c)
        where    a = elements connected to slice
                 b = nodes of each element
                 c = degree of freedom at each node
        however the degree of freedoms are no more ordered in global order,
        but on the local element basis. therefore sum([:,:-1,:])
        for the hinge force is not possible, sort_dofs solves the
        problem by sorting the values in respect to the dof number,
        which is ordered globally.
        """
        if size(unsorted) != size(dofs):
            raise "--- array must have same size ---"
        #order of dofs
        order = argsort(dofs, axis=1)[:, :, 0]
        sorted = zeros_like(unsorted)
        for i, elem in enumerate(order):
            sorted[i] = unsorted[i][elem]
        return sorted

