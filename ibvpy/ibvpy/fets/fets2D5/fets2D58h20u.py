from etsproxy.traits.api import \
    Instance, Int, Property, Array, cached_property, Instance

from numpy import \
     zeros, dot, hstack, vstack, identity, cross, transpose, tensordot, outer

from scipy.linalg import \
     inv, norm

from ibvpy.fets.fets_eval import FETSEval
from ibvpy.fets.fets3D.fets3D8h20u import FETS3D8H20U

# import the 3D elastic matrix in order
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.mats.mats2D.mats2D_tensor import \
    map2d_eps_mtx_to_eng, map2d_sig_eng_to_mtx

from ibvpy.mats.mats3D.mats3D_tensor import \
    map3d_eps_eng_to_mtx, map3d_sig_eng_to_mtx, map3d_sig_mtx_to_eng, map3d_eps_mtx_to_eng, map3d_tns2_to_tns4, map3d_tns4_to_tns2, map3d_eps_mtx_to_eng

from coordinate_transformations import \
    get_Trans_strain_from_glb_to_loc, get_Trans_stress_from_loc_to_glb

# from i_fets_eval import IFETSEval

#-----------------------------------------------------------------------------------
# FETS2D58H20U -  degenerated subparametric 3D volume element using a 2D material model
#                 quadratic (serendipity)
#-----------------------------------------------------------------------------------

class FETS2D58H20U(FETS3D8H20U):

    # @todo: the calculation of the direction cosines
    #        should be cached (similar to J_mtx and B_mtx)

    def get_dircos_mtx(self, sctx):
        '''returns the direction cosine at the integration point.
        The three column vectors v1, v2, v3 of the direction cosine
        matrix build an orthonormal system and can be used to rotate
        the 3D strains into the local coordinate system.
        The vector v3 stands perpendicular to the plane defined by the
        curvilinear coordinates r and s. Note that for the orientation
        of the local in-plane coordinates different choices exist.'''

        # get the Jacobian matrix
        r_pnt = sctx.r_pnt
        X_mtx = sctx.X
        J_mtx = self.get_J_mtx(r_pnt, X_mtx)

        # the first two rows of the Jacobi matrix are tangential
        # to the mid-surface at the integration point.
        # (cf. Zienkiewicz Vol.2, Eq.(15.19))
        j1 = J_mtx[0, :]
        j2 = J_mtx[1, :]

        # v3 stands perdendicular to the surface
        v3 = cross(j1, j2)

        # One possible choice for the direction of the local x-axis
        # is to define it such that it coincides with the orientation
        # of the first curvilinear coordinate r (cf. FE-Skript)
        v1 = j1

        # The local y-axis is obtained as the cross product of the
        # local x and z direction
        v2 = cross(v3, v1)

        # normalize the vectors to obtain an ortho-normal basis
        v1_ = v1 / norm(v1)
        v2_ = v2 / norm(v2)
        v3_ = v3 / norm(v3)

        return transpose(vstack((v1_, v2_, v3_)))

    '''get_Trans_strain_from_glb_to_loc(dircos_mtx): NOTE: imported as global method
            Transformation matrix in 3d shape = (6,6).
            multiply the strain vector (engineering notation) with
            the returned matrix in order to obtain the strain vector
            (engineering notation) in the local coordinate system.
            The ordering of the strain vector in engineering notation
            corresponds to the ordering
            [eps11,eps22,eps33,eps23,eps13,eps12].
            The parameter 'dircos' is the matrix of the direction
            cosine describing the relation between the global and
            local coordinate system.
            (cf. Zienkiewicz Eq.(6.24))'''

    '''get_Trans_stress_from_loc_to_glb(dircos_mtx): NOTE: imported as global method
            Transformation matrix in 3d shape = (6,6)
            multiply the local stress vector (engineering notation) with
            the returned matrix in order to obtain the stress vector
            (engineering notation) in the global coordinate system.
            The ordering of the stress vector in engineering notation
            corresponds to the ordering [sig11,sig22,sig33,sig23,sig13,sig12].
            The parameter 'dircos' is the matrix of the direction
            cosine describing the relation between the global and
            local coordinate system.
            (cf. Zienkiewicz Eq.(6.23) reformulated to sigma instead of sigma' by twice multiplication with
            transformation matrix T (Note: T is orthonormal, i.e. T^(-1) = T^t)'''

    def get_eps_eng(self, sctx, u):

        eps_eng3D = super(FETS2D58H20U, self).get_eps_eng(sctx, u)

        # get the direction cosine at the ip:
        dircos_mtx = self.get_dircos_mtx(sctx)

        # get the Transformation matrices for switching between the
        # local and global coordinate system in 3D case:
        T_strain_to_loc_mtx = get_Trans_strain_from_glb_to_loc(dircos_mtx)

        # transform strain from global to local coordinates
        #
        eps_eng3D_ = dot(T_strain_to_loc_mtx, eps_eng3D)

        return eps_eng3D_

    def get_mtrl_corr_pred(self, sctx, eps_eng3D, d_eps_eng3D, tn, tn1, eps_avg=None):
        '''
        Overload the 3D material-corrector-predictor in
        order to transform the strain and stress tensor
        into the in-plane so that oriented material model formulation
        can be applied.

        This is performed
        employing the direction cosines at each integration
        point (ip). The direction cosines establish the link
        between the global and local coordinate system (the
        reference to the local coordinate system is indicated
        by the underline character at the end of a variable
        name).

        After employing the material model the returned
        stiffness and stresses are transformed back to the global
        coordinate system.
        '''

        # get the direction cosine at the ip:
        dircos_mtx = self.get_dircos_mtx(sctx)

        # get the Transformation matrices for switching between the
        # local and global coordinate system in 3D case:
        T_stress_from_loc_mtx = get_Trans_stress_from_loc_to_glb(dircos_mtx)
        T_strain_to_loc_mtx = get_Trans_strain_from_glb_to_loc(dircos_mtx)

        # transform strain from global to local coordinates
        #
        eps_eng3D_ = dot(T_strain_to_loc_mtx, eps_eng3D)
        d_eps_eng3D_ = dot(T_strain_to_loc_mtx, d_eps_eng3D)

        # evaluate the mtrl-model
        #
        sig_eng3D_, D_mtx3D_ = self.mats_eval.get_corr_pred(sctx, eps_eng3D_, d_eps_eng3D_,
                                                             tn, tn1)

        #------------------------------------------------------------------------
        # get predictor (D_mtx)
        #------------------------------------------------------------------------

        # transform D_mtx3D from local to global coordinates
        #
        D_mtx3D = dot(T_stress_from_loc_mtx, dot(D_mtx3D_, T_strain_to_loc_mtx))

        #------------------------------------------------------------------------
        # get corrector (sig_eng)
        #------------------------------------------------------------------------

        # transform sig_mtx3D from local to global coordinates
        #
        sig_eng3D = dot(T_stress_from_loc_mtx, sig_eng3D_)

        return sig_eng3D, D_mtx3D
