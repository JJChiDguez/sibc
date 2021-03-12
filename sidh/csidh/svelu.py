from sidh.common import attrdict
from sidh.csidh.poly_mul import Poly_mul
from sidh.csidh.poly_redc import Poly_redc
from sidh.math import hamming_weight, bitlength, isequal
from sidh.constants import ijk_data, parameters

import numpy
from sympy import symbols, floor, sqrt, sign

def Svelu(curve, tuned, multievaluation):
    fp = curve.fp
    poly_mul = Poly_mul(curve)
    poly_redc = Poly_redc(poly_mul)
    global_L = curve.L
    prime = curve.prime
    SCALED_REMAINDER_TREE = multievaluation
    n = parameters['csidh'][prime]['n']


    cEVAL = lambda l: numpy.array([2.0 * (l - 1.0), 2.0, (l + 1.0)])
    cISOG = lambda l: numpy.array(
        [
            (3.0 * l + 2.0 * hamming_weight(l) - 9.0 + isequal[l == 3] * 4.0),
            (l + 2.0 * bitlength(l) + 1.0 + isequal[l == 3] * 2.0),
            (3.0 * l - 7.0 + isequal[l == 3] * 6.0),
        ]
    )

    C_xEVAL = list(
        map(cEVAL, global_L)
    )  # list of the costs of each degree-l isogeny evaluation
    C_xISOG = list(
        map(cISOG, global_L)
    )  # list of the costs of each degree-l isogeny construction

    # Global variables to be used in KPs, xISOG, and xEVAL

    # Here, J is a set of cardinality sJ
    J = None
    sJ = None

    # Here, ptree_I corresponds with the product tree determined by I, and I is a set of cardinality sJ
    ptree_hI = None
    sI = (None,)

    # Here, K is a set of cardinality sK
    K = None
    sK = None

    # An extra nonlocal variable which is used in xISOG and xEVAL
    XZJ4 = None


    # Next functions is used for setting the cardinalities sI, sJ, and sK
    def set_parameters_velu(b, c, i):

        nonlocal sJ
        nonlocal sI
        nonlocal sK

        assert b <= c

        # At this step, everythin is correct
        sJ = b
        sI = c
        d = ((global_L[i] - 2 - 4 * b * c - 1) // 2) + 1
        assert d >= 0
        sK = d
        return None


    def print_parameters_velu():

        print("| sI: %3d, sJ: %3d, sK: %3d |" % (sI, sJ, sK), end="")
        return None


    # KPs computes x([i]P), x([j]P), and x([k]P) for each i in I, j in J, and j in J.
    # I, J, and K are defined according to Examples 4.7 and 4.12 of https://eprint.iacr.org/2020/341
    def KPs(P, A, i):

        # Global variables to be used
        nonlocal J
        nonlocal sJ
        nonlocal ptree_hI
        nonlocal sI
        nonlocal K
        nonlocal sK

        # This functions computes all the independent data from the input of algorithm 2 of https://eprint.iacr.org/2020/341
        if sI == 0:

            # Case global_L[i] = 3 is super very special case (nothing to do)
            J = []
            ptree_hI = None
            K = [list(P)]
            # J, b, ptree_hI, c, K, d
            assert sJ == 0 and sI == 0 and sK == 1
            return None

        # We need to ensure sI is greater than or equal sJ
        assert sI >= sJ

        # Additionally, sK should be greater than or equal to zero. If that is not the case, then sJ and sI were badly chosen
        assert sK >= 0

        if sI == 1:
            # Branch corresponds with global_L[i] = 5 and global_L[i] = 7
            # Recall sJ > 0, then sJ = 1
            assert sJ == 1
            P2 = curve.xDBL(P, A)

            J = [list(P)]

            I = [list(P2)]
            hI = [
                list([fp.fp_sub(0, P2[0]), P2[1]])
            ]  # we only need to negate x-coordinate of each point
            ptree_hI = poly_mul.product_tree(hI, sI)  # product tree of hI

            if not SCALED_REMAINDER_TREE:
                # Using remainder trees
                ptree_hI = poly_redc.reciprocal_tree(
                    {'rpoly': [1], 'rdeg': 0, 'fpoly': [1], 'fdeg': 0, 'a': 1},
                    2 * sJ + 1,
                    ptree_hI,
                    sI,
                )  # reciprocal of the root is used to compute sons' reciprocals

            else:
                # Using scaled remainder trees
                assert (2 * sJ - sI + 1) > sI
                ptree_hI['reciprocal'], ptree_hI['a'] = poly_redc.reciprocal(
                    ptree_hI['poly'][::-1], sI + 1, 2 * sJ - sI + 1
                )
                ptree_hI['scaled'], ptree_hI['as'] = (
                    list(ptree_hI['reciprocal'][:sI]),
                    ptree_hI['a'],
                )

            # Notice, 0 <= sK <= 1
            assert sK <= 1
            if sK == 1:
                K = [list(P2)]
            else:
                K = []

            return None

        # At this step, sI > 1
        assert sI > 1
        if sJ == 1:
            # This branch corresponds with global_L[i] = 11 and global_L[i] = 13
            Q = curve.xDBL(P, A)  # x([2]P)
            Q2 = curve.xDBL(Q, A)  # x([2]Q)

            J = [list(P)]

            I = [[0, 0]] * sI
            I[0] = list(Q)  # x(   Q)
            I[1] = curve.xADD(Q2, I[0], I[0])  # x([3]Q)
            for ii in range(2, sI, 1):
                I[ii] = curve.xADD(I[ii - 1], Q2, I[ii - 2])  # x([2**i + 1]Q)

            hI = [
                [fp.fp_sub(0, iP[0]), iP[1]] for iP in I
            ]  # we only need to negate x-coordinate of each point
            ptree_hI = poly_mul.product_tree(hI, sI)  # product tree of hI

            if not SCALED_REMAINDER_TREE:
                # Using remainder trees
                ptree_hI = poly_redc.reciprocal_tree(
                    {'rpoly': [1], 'rdeg': 0, 'fpoly': [1], 'fdeg': 0, 'a': 1},
                    2 * sJ + 1,
                    ptree_hI,
                    sI,
                )  # reciprocal of the root is used to compute sons' reciprocals

            else:
                # Using scaled remainder trees
                assert (2 * sJ - sI + 1) <= sI
                ptree_hI['scaled'], ptree_hI['as'] = poly_redc.reciprocal(
                    ptree_hI['poly'][::-1], sI + 1, sI
                )
                ptree_hI['reciprocal'], ptree_hI['a'] = (
                    list(ptree_hI['scaled'][: (2 * sJ - sI + 1)]),
                    ptree_hI['as'],
                )

            # Notice, 0 <= sK <= 1
            assert sK <= 1
            if sK == 1:
                K = [list(Q)]
            else:
                K = []

            return None

        # Now, we ensure sI >= sJ > 1
        assert sJ > 1

        # In other words, we can proceed by the general case

        # ------------------------------------------------
        # Computing [j]P for each j in {1, 3, ..., 2*sJ - 1}
        J = [[0, 0]] * sJ
        J[0] = list(P)  # x(   P)
        P2 = curve.xDBL(P, A)  # x([2]P)
        J[1] = curve.xADD(P2, J[0], J[0])  # x([3]P)
        for jj in range(2, sJ, 1):
            J[jj] = curve.xADD(J[jj - 1], P2, J[jj - 2])  # x([2*jj + 1]P)

        # -------------------------------------------------------
        # Computing [i]P for i in { (2*sJ) * (2i + 1) : 0 <= i < sI}
        bhalf_floor = sJ // 2
        bhalf_ceil = sJ - bhalf_floor
        P4 = curve.xDBL(P2, A)  # x([4]P)
        P2[0], P4[0] = fp.fp_cswap(P2[0], P4[0], sJ % 2)  # Constant-time swap
        P2[1], P4[1] = fp.fp_cswap(
            P2[1], P4[1], sJ % 2
        )  # x([4]P) <--- coditional swap ---> x([2]P)
        Q = curve.xADD(J[bhalf_ceil], J[bhalf_floor - 1], P2)  # Q := [2b]P
        P2[0], P4[0] = fp.fp_cswap(P2[0], P4[0], sJ % 2)  # Constant-time swap
        P2[1], P4[1] = fp.fp_cswap(
            P2[1], P4[1], sJ % 2
        )  # x([4]P) <--- coditional swap ---> x([2]P)

        I = [[0, 0]] * sI
        I[0] = list(Q)  # x(   Q)
        Q2 = curve.xDBL(Q, A)  # x([2]Q)
        I[1] = curve.xADD(Q2, I[0], I[0])  # x([3]Q)
        for ii in range(2, sI, 1):
            I[ii] = curve.xADD(I[ii - 1], Q2, I[ii - 2])  # x([2**i + 1]Q)

        # --------------------------------------------------------------
        # Computing [k]P for k in { 4*sJ*sI + 1, ..., l - 6, l - 4, l - 2}
        K = [[0, 0]] * sK

        if sK >= 1:
            K[0] = list(P2)  # x([l - 2]P) = x([-2]P) = x([2]P)

        if sK >= 2:
            K[1] = list(P4)  # x([l - 4]P) = x([-4]P) = x([4]P)

        for k in range(2, sK, 1):
            K[k] = curve.xADD(K[k - 1], P2, K[k - 2])

        # ------------------------------------------------------------------------------------------------------
        #                   ~~~~~~~~               ~~~~~~~~
        #                    |    |                 |    |
        # Computing h_I(W) = |    | (W - x([i]P)) = |    | (Zi * W - Xi) / Zi where x([i]P) = Xi/Zi
        #                    i in I                 i in I
        # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

        hI = [
            [fp.fp_sub(0, iP[0]), iP[1]] for iP in I
        ]  # we only need to negate x-coordinate of each point
        ptree_hI = poly_mul.product_tree(hI, sI)  # product tree of hI

        if not SCALED_REMAINDER_TREE:
            # Using scaled remainder trees
            ptree_hI = poly_redc.reciprocal_tree(
                {'rpoly': [1], 'rdeg': 0, 'fpoly': [1], 'fdeg': 0, 'a': 1},
                2 * sJ + 1,
                ptree_hI,
                sI,
            )  # reciprocal of the root is used to compute sons' reciprocals

        else:
            # Using scaled remainder trees
            if sI < (2 * sJ - sI + 1):
                ptree_hI['reciprocal'], ptree_hI['a'] = poly_redc.reciprocal(
                    ptree_hI['poly'][::-1], sI + 1, 2 * sJ - sI + 1
                )
                ptree_hI['scaled'], ptree_hI['as'] = (
                    list(ptree_hI['reciprocal'][:sI]),
                    ptree_hI['a'],
                )

            else:
                ptree_hI['scaled'], ptree_hI['as'] = poly_redc.reciprocal(
                    ptree_hI['poly'][::-1], sI + 1, sI
                )
                ptree_hI['reciprocal'], ptree_hI['a'] = (
                    list(ptree_hI['scaled'][: (2 * sJ - sI + 1)]),
                    ptree_hI['as'],
                )

        # Polynomial D_j of algorithm 2 from https://eprint.iacr.org/2020/341 is not required,
        # but we need some some squares and products determined by list J
        # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

        # Ensuring the cardinality of each ser coincide with the expected one
        assert len(I) == sI
        assert len(J) == sJ
        assert len(K) == sK

        return None


    # Next function perform algorithm 2 of https://eprint.iacr.org/2020/341 with input \alpha = 1 and \alpha = -1, and
    # then it computes the isogenous Montgomery curve coefficient
    def xISOG(A, i):

        nonlocal J
        nonlocal sJ
        nonlocal ptree_hI
        nonlocal sI
        nonlocal K
        nonlocal sK
        nonlocal XZJ4

        AA = fp.fp_add(A[0], A[0])  # 2A' + 4C
        AA = fp.fp_sub(AA, A[1])  # 2A'
        AA = fp.fp_add(AA, AA)  # 4A'

        # Polynomial D_j of algorithm 2 from https://eprint.iacr.org/2020/341 is not required,
        # but we need some some squares and products determined by list J
        # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

        SUB_SQUARED = [0 for j in range(0, sJ, 1)]  #
        ADD_SQUARED = [0 for j in range(0, sJ, 1)]  #

        # List XZJ4 is required for degree-l isogeny evaluations...
        XZJ4 = [0 for j in range(0, sJ, 1)]  # 2*Xj*Zj
        for j in range(0, sJ, 1):

            SUB_SQUARED[j] = fp.fp_sub(J[j][0], J[j][1])  # (Xj - Zj)
            SUB_SQUARED[j] = fp.fp_sqr(SUB_SQUARED[j])  # (Xj - Zj)^2

            ADD_SQUARED[j] = fp.fp_add(J[j][0], J[j][1])  # (Xj + Zj)
            ADD_SQUARED[j] = fp.fp_sqr(ADD_SQUARED[j])  # (Xj + Zj)^2

            XZJ4[j] = fp.fp_sub(SUB_SQUARED[j], ADD_SQUARED[j])  # -4*Xj*Zj

        # --------------------------------------------------------------------------------------------------
        #                   ~~~~~~~~
        #                    |    |
        # Computing E_J(W) = |    | [ F0(W, x([j]P)) * alpha^2 + F1(W, x([j]P)) * alpha + F2(W, x([j]P)) ]
        #                    j in J
        # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates
        # In particular, for a degree-l isogeny construction, we need alpha = 1 and alpha = -1

        # EJ_0 is the one determined by alpha = 1
        EJ_0 = [[0, 0, 0] for j in range(0, sJ, 1)]
        # EJ_1 is the one determined by alpha = -1
        EJ_1 = [[0, 0, 0] for j in range(0, sJ, 1)]

        for j in range(0, sJ, 1):

            # However, each SUB_SQUARED[j] and ADD_SQUARED[j] should be multiplied by C
            tadd = fp.fp_mul(ADD_SQUARED[j], A[1])
            tsub = fp.fp_mul(SUB_SQUARED[j], A[1])

            # We require the double of tadd and tsub
            tadd2 = fp.fp_add(tadd, tadd)
            tsub2 = fp.fp_add(tsub, tsub)

            t1 = fp.fp_mul(XZJ4[j], AA)  #       A *(-4*Xj*Zj)

            # Case alpha = 1
            linear = fp.fp_sub(
                t1, tadd2
            )  #       A *(-4*Xj*Zj)  - C * (2 * (Xj + Zj)^2)
            EJ_0[j] = [tsub, linear, tsub]

            # Case alpha = -1
            linear = fp.fp_sub(
                tsub2, t1
            )  #       C * (2 * (Xj - Zj)^2) - A *(-4*Xj*Zj)
            EJ_1[j] = [tadd, linear, tadd]

        # The faster way for multiplying is using a divide-and-conquer approach
        poly_EJ_0 = poly_mul.product_selfreciprocal_tree(EJ_0, sJ)[
            'poly'
        ]  # product tree of EJ_0 (we only require the root)
        poly_EJ_1 = poly_mul.product_selfreciprocal_tree(EJ_1, sJ)[
            'poly'
        ]  # product tree of EJ_1 (we only require the root)

        if not SCALED_REMAINDER_TREE:
            # Remainder tree computation
            remainders_EJ_0 = poly_redc.multieval_unscaled(
                poly_EJ_0, 2 * sJ + 1, ptree_hI, sI
            )
            remainders_EJ_1 = poly_redc.multieval_unscaled(
                poly_EJ_1, 2 * sJ + 1, ptree_hI, sI
            )

        else:
            # Approach using scaled remainder trees
            if ptree_hI != None:
                poly_EJ_0 = poly_redc.poly_redc(poly_EJ_0, 2 * sJ + 1, ptree_hI)
                fg_0 = poly_mul.poly_mul_middle(ptree_hI['scaled'], sI, poly_EJ_0[::-1], sI)
                remainders_EJ_0 = poly_redc.multieval_scaled(
                    fg_0[::-1], sI, [1] + [0] * (sI - 1), sI, ptree_hI, sI
                )

                poly_EJ_1 = poly_redc.poly_redc(poly_EJ_1, 2 * sJ + 1, ptree_hI)
                fg_1 = poly_mul.poly_mul_middle(ptree_hI['scaled'], sI, poly_EJ_1[::-1], sI)
                remainders_EJ_1 = poly_redc.multieval_scaled(
                    fg_1[::-1], sI, [1] + [0] * (sI - 1), sI, ptree_hI, sI
                )
            else:
                remainders_EJ_0 = []
                remainders_EJ_1 = []

        # Multipying all the remainders
        r0 = poly_mul.product(remainders_EJ_0, sI)
        r1 = poly_mul.product(remainders_EJ_1, sI)

        # ---------------------------------------------------------------------------------
        # Now, we proceed by computing the missing part which is determined by K
        # Notice, the denominators are the same and then they annulled between them
        # In other words, it is not required to compute the product of all Zk's with k In K

        # Case alpha = 1
        hK_0 = [[fp.fp_sub(K[k][1], K[k][0])] for k in range(0, sK, 1)]
        hK_0 = poly_mul.product(hK_0, sK)  # product of (Zk - Xk) for each k in K
        # Case alpha = -1
        hK_1 = [[fp.fp_add(K[k][1], K[k][0])] for k in range(0, sK, 1)]
        hK_1 = poly_mul.product(hK_1, sK)  # product of (Zk + Xk) for each k in K

        # --------------------------------------------------------------
        # Now, we have all the ingredients for computing the image curve
        A24m = fp.fp_sub(A[0], A[1])  # A' - 2C

        A24 = fp.fp_exp(A[0], global_L[i])  # (A' + 2C)^l
        A24m = fp.fp_exp(A24m, global_L[i])  # (A' - 2C)^l

        t24m = fp.fp_mul(
            hK_1, r1
        )  # output of algorithm 2 with alpha =-1 and without the demoninator
        t24m = fp.fp_sqr(t24m)  # raised at 2
        t24m = fp.fp_sqr(t24m)  # raised at 4
        t24m = fp.fp_sqr(t24m)  # raised at 8

        t24 = fp.fp_mul(
            hK_0, r0
        )  # output of algorithm 2 with alpha = 1 and without the demoninator
        t24 = fp.fp_sqr(t24)  # raised at 2
        t24 = fp.fp_sqr(t24)  # raised at 4
        t24 = fp.fp_sqr(t24)  # raised at 8

        A24 = fp.fp_mul(A24, t24m)
        A24m = fp.fp_mul(A24m, t24)

        # Now, we have d = (A24m / A24) where the image Montgomery cuve coefficient is
        #      B'   2*(1 + d)   2*(A24 + A24m)
        # B = ---- = --------- = --------------
        #      C      (1 - d)     (A24 - A24m)
        # However, we required B' + 2C = 4*A24 and 4C = 4 * (A24 - A24m)

        t24m = fp.fp_sub(A24, A24m)  #   (A24 - A24m)
        t24m = fp.fp_add(t24m, t24m)  # 2*(A24 - A24m)
        t24m = fp.fp_add(t24m, t24m)  # 4*(A24 - A24m)

        t24 = fp.fp_add(A24, A24)  # 2 * A24
        t24 = fp.fp_add(t24, t24)  # 4 * A24

        # return [t24, t24m], ptree_hI, XZJ4
        return [t24, t24m]


    def xEVAL(P, A):

        AA = fp.fp_add(A[0], A[0])  # 2A' + 4C
        AA = fp.fp_sub(AA, A[1])  # 2A'
        AA = fp.fp_add(AA, AA)  # 4A'

        # --------------------------------------------------------------------------------------------------
        #                   ~~~~~~~~
        #                    |    |
        # Computing E_J(W) = |    | [ F0(W, x([j]P)) * alpha^2 + F1(W, x([j]P)) * alpha + F2(W, x([j]P)) ]
        #                    j in J
        # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates
        # In particular, for a degree-l isogeny construction, we need alpha = 1 and alpha = -1

        # EJ_0 is the one determined by alpha = x
        EJ_0 = [[0, 0, 0] for j in range(0, sJ, 1)]
        # Notice, the corresponding EJ_1 that is determined by alpha = 1/x can be computed by using EJ_0

        XZ_add = fp.fp_add(P[0], P[1])  # X + Z
        XZ_sub = fp.fp_sub(P[0], P[1])  # X - Z

        AXZ2 = fp.fp_mul(P[0], P[1])  # X * Z
        t1 = fp.fp_sqr(P[0])  # X^2
        t2 = fp.fp_sqr(P[1])  # Z^2

        CX2Z2 = fp.fp_add(t1, t2)  #      X^2 + Z^2
        CX2Z2 = fp.fp_mul(CX2Z2, A[1])  # C * (X^2 + Z^2)

        AXZ2 = fp.fp_add(AXZ2, AXZ2)  #       2 * (X * Z)
        CXZ2 = fp.fp_mul(AXZ2, A[1])  # C  * [2 * (X * Z)]
        AXZ2 = fp.fp_mul(AXZ2, AA)  # A' * [2 * (X * Z)]

        for j in range(0, sJ, 1):

            XZj_add = fp.fp_add(J[j][0], J[j][1])  # Xj + Zj
            XZj_sub = fp.fp_sub(J[j][0], J[j][1])  # Xj - Zj

            t1 = fp.fp_mul(XZ_sub, XZj_add)  # (X - Z) * (Xj + Zj)
            t2 = fp.fp_mul(XZ_add, XZj_sub)  # (X + Z) * (Xj - Zj)

            # Computing the quadratic coefficient
            quadratic = fp.fp_sub(t1, t2)  #   2 * [(X*Zj) - (Z*Xj)]
            quadratic = fp.fp_sqr(quadratic)  # ( 2 * [(X*Zj) - (Z*Xj)] )^2
            quadratic = fp.fp_mul(A[1], quadratic)  # C * ( 2 * [(X*Zj) - (Z*Xj)] )^2

            # Computing the constant coefficient
            constant = fp.fp_add(t1, t2)  #   2 * [(X*Xj) - (Z*Zj)]
            constant = fp.fp_sqr(constant)  # ( 2 * [(X*Xj) - (Z*Zj)] )^2
            constant = fp.fp_mul(A[1], constant)  # C * ( 2 * [(X*Xj) - (Z*Zj)] )^2

            # Computing the linear coefficient
            # ----------------------------------------------------------------------------------------------------------
            # C * [ (-2*Xj*Zj)*(alpha^2 + 1) + (-2*alpha)*(Xj^2 + Zj^2)] + [A' * (-2*Xj*Zj) * (2*X*Z)] where alpha = X/Z
            t1 = fp.fp_add(J[j][0], J[j][1])  #     (Xj + Zj)
            t1 = fp.fp_sqr(t1)  #     (Xj + Zj)^2
            t1 = fp.fp_add(t1, t1)  # 2 * (Xj + Zj)^2
            t1 = fp.fp_add(
                t1, XZJ4[j]
            )  # 2 * (Xj + Zj)^2 - (4*Xj*Zj) := 2 * (Xj^2 + Zj^2)
            t1 = fp.fp_mul(t1, CXZ2)  # [2 * (Xj^2 + Zj^2)] * (2 * [ C * (X * Z)])

            t2 = fp.fp_mul(CX2Z2, XZJ4[j])  # [C * (X^2 + Z^2)] * (-4 * Xj * Zj)
            t1 = fp.fp_sub(
                t2, t1
            )  # [C * (X^2 + Z^2)] * (-4 * Xj * Zj) - [2 * (Xj^2 + Zj^2)] * (2 * [ C * (X * Z)])

            t2 = fp.fp_mul(AXZ2, XZJ4[j])  # (2 * [A' * (X * Z)]) * (-4 * Xj * Zj)
            linear = fp.fp_add(
                t1, t2
            )  # This is our desired equation but multiplied by 2
            linear = fp.fp_add(
                linear, linear
            )  # This is our desired equation but multiplied by 4
            # ----------------------------------------------------------------------------------------------------------

            # Case alpha = X / Z
            EJ_0[j] = [constant, linear, quadratic]

        # The faster way for multiplying is using a divide-and-conquer approach
        poly_EJ_0 = poly_mul.product_tree(EJ_0, sJ)[
            'poly'
        ]  # product tree of EJ_0 (we only require the root)
        poly_EJ_1 = list(
            poly_EJ_0[::-1]
        )  # product tree of EJ_1(x) = x^{2b + 1} EJ_0(1/X)

        if not SCALED_REMAINDER_TREE:
            # Remainder tree computation
            remainders_EJ_0 = poly_redc.multieval_unscaled(
                poly_EJ_0, 2 * sJ + 1, ptree_hI, sI
            )
            remainders_EJ_1 = poly_redc.multieval_unscaled(
                poly_EJ_1, 2 * sJ + 1, ptree_hI, sI
            )

        else:
            # Approach using scaled remainder trees
            if ptree_hI != None:
                poly_EJ_0 = poly_redc.poly_redc(poly_EJ_0, 2 * sJ + 1, ptree_hI)
                fg_0 = poly_mul.poly_mul_middle(ptree_hI['scaled'], sI, poly_EJ_0[::-1], sI)
                remainders_EJ_0 = poly_redc.multieval_scaled(
                    fg_0[::-1], sI, [1] + [0] * (sI - 1), sI, ptree_hI, sI
                )

                poly_EJ_1 = poly_redc.poly_redc(poly_EJ_1, 2 * sJ + 1, ptree_hI)
                fg_1 = poly_mul.poly_mul_middle(ptree_hI['scaled'], sI, poly_EJ_1[::-1], sI)
                remainders_EJ_1 = poly_redc.multieval_scaled(
                    fg_1[::-1], sI, [1] + [0] * (sI - 1), sI, ptree_hI, sI
                )
            else:
                remainders_EJ_0 = []
                remainders_EJ_1 = []

        # Multipying all the remainders
        r0 = poly_mul.product(remainders_EJ_0, sI)
        r1 = poly_mul.product(remainders_EJ_1, sI)

        # ---------------------------------------------------------------------------------
        # Now, we proceed by computing the missing part which is determined by K
        # Notice, the denominators are the same and then they annulled between them
        # In other words, it is not required to compute the product of all Zk's with k In K

        hK_0 = [[0]] * sK
        hK_1 = [[0]] * sK
        for k in range(0, sK, 1):

            XZk_add = fp.fp_add(K[k][0], K[k][1])  # Xk + Zk
            XZk_sub = fp.fp_sub(K[k][0], K[k][1])  # Xk - Zk
            t1 = fp.fp_mul(XZ_sub, XZk_add)  # (X - Z) * (Xk + Zk)
            t2 = fp.fp_mul(XZ_add, XZk_sub)  # (X + Z) * (Xk - Zk)

            # Case alpha = X/Z
            hK_0[k] = [fp.fp_sub(t1, t2)]  # 2 * [(X*Zk) - (Z*Xk)]

            # Case 1/alpha = Z/X
            hK_1[k] = [fp.fp_add(t1, t2)]  # 2 * [(X*Xk) - (Z*Zk)]

        hK_0 = poly_mul.product(hK_0, sK)  # product of (XZk - ZXk) for each k in K
        hK_1 = poly_mul.product(hK_1, sK)  # product of (XXk - ZZk) for each k in K

        # ---------------------------------------------------------------------------------
        # Now, unifying all the computations
        XX = fp.fp_mul(
            hK_1, r1
        )  # output of algorithm 2 with 1/alpha = Z/X and without the demoninator
        XX = fp.fp_sqr(XX)
        XX = fp.fp_mul(XX, P[0])

        ZZ = fp.fp_mul(
            hK_0, r0
        )  # output of algorithm 2 with alpha = X/Z and without the demoninator
        ZZ = fp.fp_sqr(ZZ)
        ZZ = fp.fp_mul(ZZ, P[1])

        return [XX, ZZ]


    # Get cost of the isogeny constructions and evaluations
    sI_list = None
    sJ_list = None


    def cISOG_and_cEVAL():

        nonlocal C_xISOG
        nonlocal C_xEVAL

        nonlocal sI_list
        nonlocal sJ_list

        if tuned:

            sI_list = []
            sJ_list = []
            f = open(ijk_data + prime)

            for i in range(0, n, 1):

                bc = f.readline()
                bc = [int(bci) for bci in bc.split()]
                sJ_list.append(bc[0])
                sI_list.append(bc[1])

            f.close()
        # First, we look for a full torsion point
        A = [2, 4]
        T_p, T_m = curve.full_torsion_points(A)

        for i in range(0, n, 1):

            if tuned:
                set_parameters_velu(sJ_list[i], sI_list[i], i)
            else:
                # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                # These paramters are required in KPs, xISOG, and xEVAL
                if global_L[i] == 3:
                    b = 0
                    c = 0
                else:
                    b = int(floor(sqrt(global_L[i] - 1) / 2.0))
                    c = int(floor((global_L[i] - 1.0) / (4.0 * b)))

                set_parameters_velu(b, c, i)

            # Getting an orderl-l point
            Tp = list(T_p)
            for j in range(i + 1, n, 1):
                Tp = curve.xMUL(Tp, A, j)

            # Cost of xISOG() and KPs()
            fp.set_zero_ops()
            KPs(Tp, A, i)
            t = fp.get_ops()
            C_xISOG[i] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            fp.set_zero_ops()
            B = xISOG(A, i)
            t = fp.get_ops()
            C_xISOG[i] += numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            # xEVAL: kernel point determined by the next isogeny evaluation
            fp.set_zero_ops()
            T_p = xEVAL(T_p, A)

            # Cost of xEVAL
            fp.set_zero_ops()
            T_m = xEVAL(T_m, A)
            t = fp.get_ops()
            C_xEVAL[i] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            # Updating the new next curve
            A = list(B)

        fp.set_zero_ops()
        assert curve.coeff(A) == 0x6
        return None


    # Now, we proceed to store all the correct costs
    cISOG_and_cEVAL()

    return attrdict(name='svelu', **locals())
