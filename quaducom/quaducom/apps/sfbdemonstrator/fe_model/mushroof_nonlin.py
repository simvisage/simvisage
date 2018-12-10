#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Jul 8, 2010 by: rch

from etsproxy.traits.api import \
    Enum, Property, Instance, cached_property, Str, \
    Int, Float, List

from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
    MATS2D5MicroplaneDamage

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
    PhiFnGeneralExtended, \
    PhiFnGeneral, PhiFnStrainHardening, PhiFnStrainHardeningLinear

from quaducom.apps.sfbdemonstrator.fe_model.mush_roof_model import \
    MushRoofModel
    
from quaducom.apps.sfbdemonstrator.fe_model.mr_quarter import \
    MRquarter

from matresdev.db.matdb.trc.ccs_unit_cell import \
    CCSUnitCell, DamageFunctionEntry

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

from matresdev.db.simdb import \
    SimDB

from matresdev.db.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

from ibvpy.api import \
    TLine

from ibvpy.rtrace.rt_domain_list_field import \
    RTraceDomainListField

simdb = SimDB()

from pickle import dump, load

class MushRoofModelNonLin( MRquarter ):
    '''Overload the nonlinear model.
    '''

    #-----------------
    # composite cross section unit cell:
    #-----------------
    #
    ccs_unit_cell_key = Enum( list(CCSUnitCell.db.keys()),
                              simdb = True, input = True,
                              auto_set = False, enter_set = True )

    ccs_unit_cell_ref = Property( Instance( SimDBClass ),
                                  depends_on = 'ccs_unit_cell_key' )
    @cached_property
    def _get_ccs_unit_cell_ref( self ):
        return CCSUnitCell.db[ self.ccs_unit_cell_key ]

    # vary the failure strain in PhiFnGeneralExtended:
    factor_eps_fail = Float( 1.0, input = True,
                             ps_levels = ( 1.0, 1.2, 3 ) )

    #-----------------
    # damage function:
    #-----------------
    #
    material_model = Str( input = True )
    def _material_model_default( self ):
        # return the material model key of the first DamageFunctionEntry
        # This is necessary to avoid an ValueError at setup  
        return self.ccs_unit_cell_ref.damage_function_list[0].material_model

    calibration_test = Str( input = True
                            )
    def _calibration_test_default( self ):
        # return the material model key of the first DamageFunctionEntry
        # This is necessary to avoid an ValueError at setup  
        return self.ccs_unit_cell_ref.damage_function_list[0].calibration_test

    damage_function = Property( Instance( MFnLineArray ),
                                depends_on = 'input_change' )
    @cached_property
    def _get_damage_function( self ):
        return self.ccs_unit_cell_ref.get_param( self.material_model, self.calibration_test )

    #-----------------
    # phi function extended:
    #-----------------
    #
    phi_fn = Property( Instance( PhiFnGeneralExtended ),
                       depends_on = 'input_change,+ps_levels' )
    @cached_property
    def _get_phi_fn( self ):
        return PhiFnGeneralExtended( mfn = self.damage_function,
                                     factor_eps_fail = self.factor_eps_fail )

    #----------------------------------------------------------------------------------
    # mats_eval
    #----------------------------------------------------------------------------------

    # age of the plate at the time of testing
    # NOTE: that the same phi-function is used independent of age. This assumes a 
    # an afine/proportional damage evolution for different ages. 
    #
    age = Int( 28, #input = True
                )

    # composite E-modulus 
    #
    E_c = Property( Float, depends_on = 'input_change' )
    @cached_property
    def _get_E_c( self ):
        return self.ccs_unit_cell_ref.get_E_c_time( self.age )

    # Poisson's ratio 
    #
    nu = Property( Float, depends_on = 'input_change' )
    @cached_property
    def _get_nu( self ):
        return self.ccs_unit_cell_ref.nu

    # @todo: for mats_eval the information of the unit cell should be used
    # in order to use the same number of microplanes and model version etc...
    #
    mats = Property( Instance( MATS2D5MicroplaneDamage ),
                          depends_on = 'input_change' )
    @cached_property
    def _get_mats( self ):
        mats = MATS2D5MicroplaneDamage( 
                                E = self.E_c,
                                nu = self.nu,
                                n_mp = 30,
                                symmetrization = 'sum-type',
                                model_version = 'compliance',
                                phi_fn = self.phi_fn )

        return mats

    tline = Instance( TLine )
    def _tline_default( self ):
        return TLine( min = 0.0, step = 1.0, max = 8.0 )

    rtrace_list = List
    def _rtrace_list_default( self ):
        return [
#                                 RTraceDomainListField( name = 'Displacement' ,
#                                                var = 'u', idx = 0, warp = True ),
#                                 RTraceDomainListField( name = 'Stress' ,
#                                                var = 'sig_app', idx = 0, warp = True,
#                                                record_on = 'update', ),
#                    self.max_princ_stress,
                             RTraceDomainListField( name = 'Damage' ,
                                        var = 'omega_mtx', idx = 0, warp = True,
                                        record_on = 'update' ),
                    ]

if __name__ == '__main__':

    sim_model = MushRoofModelNonLin( n_dofs_xy = 4, shape_z = 1,
                                     ccs_unit_cell_key = 'FIL-10-09_2D-02-06a_0.00273_90_0',
                                     calibration_test = 'TT11-10a-average',
                                     age = 28 )

    do = 'ui'

    if do == 'ui':
        sim_model.tloop.eval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp( ibv_resource = sim_model )
        app.main()

    if do == 'eval':
        sim_model.tloop.eval()
        import pylab as p
        sim_model.f_w_diagram_center.refresh()
        file_path = 'f_w_diagram.pickle'
        file = open( pickle_file_path, 'w' )
        dump( sim_model.f_w_diagram.trace, file )
        file.close()
        sim_model.f_w_diagram.trace.mpl_plot( p, color = 'red' )
        p.show()

    elif do == 'show_last':
        from promod.exdb.ex_run import ExRun
        import pylab as p
        pickle_path = 'pickle_files'
        file_path = 'f_w_diagram.pickle'
        file = open( pickle_file_path, 'r' )
        trace = load( file )
        p.plot( trace.xdata, trace.ydata, color = 'blue' )


