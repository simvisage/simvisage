
from traits.api import \
    Int, implements, Array

from ibvpy.fets.fets_eval import IFETSEval, FETSEval

from numpy import array, dot

from scipy.linalg import \
     inv

from ibvpy.fets.fets2D import FETS2D4Q, FETS2D4Q8U, FETS2D4Q16U

def __demo__():
    from ibvpy.api import \
        TStepper as TS, RTraceDomainListField, TLoop, \
        TLine, BCSlice, FEDomain, FERefinementGrid, FEGrid
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic

    fets_eval = FETS2D4Q8U(mats_eval=MATS2DElastic())

    fe_domain = FEDomain()

    r1 = FERefinementGrid(fets_eval=fets_eval,
                          domain=fe_domain)
    r2 = FERefinementGrid(fets_eval=fets_eval,
                          domain=fe_domain)
    # Discretization
    domain1 = FEGrid(coord_max=(3., 3.),
                     shape=(10, 4),
                     fets_eval=fets_eval,
                     level=r1)

    domain2 = FEGrid(coord_min=(3., 0.),
                     coord_max=(6., 3),
                     shape=(10, 4),
                     fets_eval=fets_eval,
                     level=r2)

    ts = TS(dof_resultants=True,
            sdomain=[domain1, domain2],  # fe_domain,
            bcond_list=[
                          # Fix the left edge of domain1
                          BCSlice(var='u', dims=[0, 1], value=0, slice=domain1[0, :, 0, :]),
                          # Link the right edge of domain1 with the left edge of domain2
                          #
                          # note that following arrays must have the same lengths:
                          # slice and link_slice
                          # dims, link_dims and link_coeffs must have the same lengths

                          # VAR-1:
                          # linking along the complete line between 'domain1' and 'domain2'
                          # all nodes along the y-axis
                          # (used linking of more nodes at once in 'BCSlice')
                          #
                          BCSlice(var='u', dims=[0, 1], value=0.0,
                                  slice=domain1[-1, :, -1, :],
                                  link_slice=domain2[0, :, 0, :],
                                  link_dims=[0, 1],
                                  link_coeffs=[1., 1.]),

                          # VAR-2:
                          # linking along individual points between 'domain1' and 'domain2'
                          # (used linking of single nodes in 'BCSlice')
                          #
#                          BCSlice(var='u', dims=[0, 1], value=0.0,
#                                  slice=domain1[-1, -1, -1, -1],
#                                  link_slice=domain2[0, -1, 0, -1],
#                                  link_dims=[0, 1],
#                                  link_coeffs=[1., 1.]),
#                          BCSlice(var='u', dims=[0, 1], value=0.0,
#                                  slice=domain1[-1, 0, -1, 0],
#                                  link_slice=domain2[0, 0, 0, 0],
#                                  link_dims=[0, 1],
#                                  link_coeffs=[1., 1.]),

                          # Load the right edge of domain2
                          BCSlice(var='f', dims=[0], value=1, slice=domain2[-1, :, -1, :])
                          ],
            rtrace_list=[RTraceDomainListField(name='Stress' ,
                                                 var='sig_app', idx=0),
                           RTraceDomainListField(name='Displacement' ,
                                                 var='u', idx=0,
                                                 warp=True),
                           ]
            )

    # Add the time-loop control
    tloop = TLoop(tstepper=ts, debug=False,
                  tline=TLine(min=0.0, step=1.0, max=1.0))

    print '---- result ----'
    print tloop.eval()

    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource=tloop)
    app.main()

if __name__ == '__main__':
    __demo__()
