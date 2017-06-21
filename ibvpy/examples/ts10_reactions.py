'''
Created on 14.03.2011

@author: alexander

example for resulting moments derived from internal forces and coordinates
'''

if __name__ == '__main__':

    from ibvpy.fets.fets3D import FETS3D8H

    from numpy import sum, argsort, zeros_like, size

    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, BCSlice, IBVPSolve as IS, DOTSEval

    from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic

    fets_eval = FETS3D8H( mats_eval = MATS3DElastic( nu = 0.25 ) )

    from ibvpy.mesh.fe_grid import FEGrid

    # Discretization
    domain = FEGrid( coord_max = ( 1., 1., 1. ),
                           shape = ( 3, 3, 3 ),
                           fets_eval = fets_eval )

    support_z = BCSlice( var = 'u', value = 0., dims = [2],
                                  slice = domain[:, :, 0, :, :, 0] )
    support_y = BCSlice( var = 'u', value = 0., dims = [1],
                                  slice = domain[[0, -1], 0, 0, [0, -1], 0, 0] )
    support_xyz = BCSlice( var = 'u', value = 0., dims = [0, 1, 2],
                                  slice = domain[0, 0, 0, 0, 0, 0] )

    vertical_load = BCSlice( var = 'f', value = -1.00, dims = [2],
                                  slice = domain[:, :, -1, :, :, -1] )
#    horizontal_load = BCSlice( var = 'f', value = 1.00, dims = [0, 1],
#                                  slice = domain[-1, -1, -1, -1, -1, -1] )

    ts = TS( 
            sdomain = domain,
             bcond_list = [ support_z,
                            support_y,
                            support_xyz,
                           vertical_load
#                           horizontal_load,
                           ],
             rtrace_list = [
                        RTraceDomainListField( name = 'Deformation' ,
                                       var = 'eps_app', idx = 0,
                                       record_on = 'update' ),
                         RTraceDomainListField( name = 'Displacement' ,
                                        var = 'u', idx = 0 ),
                        ]
            )


    # Add the time-loop control
    tloop = TLoop( tstepper = ts,
         tline = TLine( min = 0.0, step = 0.5, max = 1.0 ) )

    tloop.eval()

    dofs = support_z.slice.dofs
    print "dofs", dofs

#TODO: check sorting,
#select corner nodes only one time

#    # corner nodes are doubled therefore undo slices
#    #
#    dofs_r0_right = append( dofs_r0_right[:, :-1, :],
#                            dofs_r0_right[-1, -1, :] ).reshape( -1, 3 )

    # sorting force of slice by number of dofs internal forces
    #
    def sort_by_dofs( dofs, unsorted ):
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
        if size( unsorted ) != size( dofs ):
            raise "--- array must have same size ---"
        #order of dofs
        order = argsort( dofs, axis = 1 )[:, :, 0]
        sorted = zeros_like( unsorted )
        for i, elem in enumerate( order ):
            sorted[i] = unsorted[i][elem]
        return sorted


    F_z = tloop.tstepper.F_int[ support_z.slice.dofs ][:, :, -1]
    X = support_z.slice.geo_X[:, :, 0]
    Y = support_z.slice.geo_X[:, :, 1]

    Xs = ( max( X.flatten() ) - min( X.flatten() ) ) / 2.0
    Ys = ( max( Y.flatten() ) - min( Y.flatten() ) ) / 2.0

    N = sum( F_z )
    Mx = sum( F_z * ( X - Xs ) )
    My = sum( F_z * ( Y - Ys ) )

    print 'N', N
    print 'Mx', Mx
    print 'My', My


    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
#    from ibvpy.plugins.ibvpy_app import IBVPyApp
#    app = IBVPyApp( ibv_resource = tloop )
#    app.main()

