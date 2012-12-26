'''
Created on Dec 21, 2012

@author: rch
'''

if __name__ == '__main__':
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, BCDofGroup, IBVPSolve as IS, DOTSEval

    from ibvpy.fets.fets3D import FETS3D8H16U

    from mathkit.mfn import MFnLineArray

    #from lib.mats.mats2D.mats_cmdm2D.mats_mdm2d import MACMDM
    from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
    from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import *
    from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic, MATS3DElasticNL

    E = 28e+3

    fets_eval = FETS3D8H16U(mats_eval = MATS3DElasticNL(stress_strain_curve = MFnLineArray(ydata = [ E, 0., 0.1 ],
                                                                                          xdata = [ -1, 0., 1.])),)

    from ibvpy.mesh.fe_grid import FEGrid

    # Discretization
    domain = FEGrid(coord_max = (3., 3., 3.),
                           shape = (3, 3, 3),
                           fets_eval = fets_eval)

    # Put the tseval (time-stepper) into the spatial context of the
    # discretization and specify the response tracers to evaluate there.
    right_dof = 2
    ts = TS(
            sdomain = domain,
             # conversion to list (square brackets) is only necessary for slicing of 
             # single dofs, e.g "get_left_dofs()[0,1]" which elsewise retuns an integer only
             bcond_list = [ BCDofGroup(var = 'u', value = 0., dims = [0],
                                        get_dof_method = domain.get_left_dofs),
                        BCDofGroup(var = 'u', value = 0., dims = [1, 2],
                                  get_dof_method = domain.get_bottom_left_dofs),
                        BCDofGroup(var = 'u', value = -0.002, dims = [0],
                                  get_dof_method = domain.get_right_dofs) ],
             rtrace_list = [
#                        RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
#                                  var_y = 'F_int', idx_y = right_dof,
#                                  var_x = 'U_k', idx_x = right_dof,
#                                  record_on = 'update'),
#                        RTraceDomainListField(name = 'Deformation' ,
#                                       var = 'eps', idx = 0,
#                                       record_on = 'update'),
#                         RTraceDomainListField(name = 'Displacement_ip' ,
#                                        var = 'u', idx = 0,
#                                        position = 'int_pnts'),
                          RTraceDomainListField(name = 'Displacement' ,
                                        var = 'u', idx = 0),
#                         RTraceDomainListField(name = 'Stress' ,
#                                        var = 'sig', idx = 0,
#                                        record_on = 'update'),
#                        RTraceDomainListField(name = 'N0' ,
#                                       var = 'N_mtx', idx = 0,
#                                       record_on = 'update')
                        ]
            )

    # Add the time-loop control
    #
    tloop = TLoop(tstepper = ts,
         tline = TLine(min = 0.0, step = 0.5, max = 1.0))

    tloop.eval()

    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource = tloop)
    app.main()
