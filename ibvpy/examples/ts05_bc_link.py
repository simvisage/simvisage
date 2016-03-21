
from traits.api import \
    Int, implements, Array

from ibvpy.fets.fets_eval import IFETSEval, FETSEval

from numpy import array, dot

from scipy.linalg import \
     inv

from ts04_fets1D2l import FETS1D2L

def __demo__():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, BCSlice, FEDomain, FERefinementGrid
    from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic

    fets_eval = FETS1D2L(mats_eval = MATS1DElastic(E = 10.))
    from ibvpy.mesh.fe_grid import FEGrid

    fe_domain = FEDomain()
    
    r1 = FERefinementGrid(fets_eval = fets_eval,
                          domain = fe_domain)
    r2 = FERefinementGrid(fets_eval = fets_eval,
                          domain = fe_domain)
    # Discretization
    domain1 = FEGrid(coord_max = (3.,),
                     shape = (3,),
                     fets_eval = fets_eval,
                     level = r1)

    domain2 = FEGrid(coord_min = (3.,),
                     coord_max = (6.,),
                     shape = (3,),
                     fets_eval = fets_eval,
                     level = r2)

    ts = TS(dof_resultants = True,
            sdomain = fe_domain,
            bcond_list = [BCSlice(var = 'u', dims = [0], value = 0, slice = domain1[0, 0]),
                          BCSlice(var = 'u', dims = [0], value = 0, slice = domain1[-1, -1],
                                  link_slice = domain2[0, 0], link_coeffs = [1.]),
                          BCSlice(var = 'f', dims = [0], value = 1, slice = domain2[-1, -1])
                          ],
            rtrace_list = [RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                                        var_y = 'F_int', idx_y = 0,
                                        var_x = 'U_k', idx_x = 1),
                           RTraceDomainListField(name = 'Stress' ,
                                                 var = 'sig_app', idx = 0),
                           RTraceDomainListField(name = 'Displacement' ,
                                                 var = 'u', idx = 0,
                                                 warp = True),
                           RTraceDomainListField(name = 'N0' ,
                                                 var = 'N_mtx', idx = 0,
                                                 record_on = 'update')
                           ]
            )

    # Add the time-loop control
    tloop = TLoop(tstepper = ts,
                  tline = TLine(min = 0.0, step = 0.5, max = 1.0))

    print '---- result ----'
    print tloop.eval()
    print ts.F_int
    print ts.rtrace_list[0].trace.ydata

    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource = tloop)
    app.main()

if __name__ == '__main__':
    __demo__()
