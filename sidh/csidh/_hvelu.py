from sidh.csidh._poly_mul import poly_mul_init, product_tree, product_selfreciprocal_tree, poly_mul_middle, product
from sidh.csidh._poly_redc import poly_redc_init, reciprocal, poly_redc, reciprocal_tree, multieval_scaled
from sidh.math import isequal, bitlength, hamming_weight
from sidh.constants import ijk_data

import numpy
from sympy import floor, sqrt, sign

class Hvelu(object):
    name = 'hvelu'
    def __init__(self, curve, _verbose):

        # Get cost of the isogeny constructions and evaluations
        self.sI_list = None
        self.sJ_list = None

        # Tradicional velu formulas will be used for l_i <= 101
        self.HYBRID_BOUND = 83

        # Global variables to be used in KPs, xISOG, and xEVAL

        # Here, J is a set of cardinality sJ
        self.J = None
        self.sJ = None

        # Here, ptree_I corresponds with the product tree determined by I, and I is a set of cardinality sJ
        self.ptree_hI = None
        self.sI = (None,)

        # Here, K is a set of cardinality sK
        self.K = None
        self.sK = None

        # An extra global variable which is used in xISOG and xEVAL
        self.XZJ4 = None

        self.SCALED_REMAINDER_TREE = True
        self.verbose = verbose = _verbose

        self.prime = curve.prime
        self.curve = curve
        self.fp = self.curve.fp
        poly_redc_init(self.fp)
        poly_mul_init(curve)
        self.global_L = self.curve.L

        self.C_xEVAL = list(
            map(self.cEVAL, self.global_L)
        )  # list of the costs of each degree-l isogeny evaluation
        self.C_xISOG = C_xISOG = list(
            map(self.cISOG, self.global_L)
        )  # list of the costs of each degree-l isogeny construction

        # Now, we proceed to store all the correct costs
        self.cISOG_and_cEVAL() # FIXME: not for general use, right?

    def cEVAL(self, l):
        numpy.array([2.0 * (l - 1.0), 2.0, (l + 1.0)])

    def cISOG(self, l):
        numpy.array(
            [
                (3.0 * l + 2.0 * hamming_weight(l) - 9.0 + isequal[l == 3] * 4.0),
                (l + 2.0 * bitlength(l) + 1.0 + isequal[l == 3] * 2.0),
                (3.0 * l - 7.0 + isequal[l == 3] * 6.0),
            ]
        )


    def yDBL(self, P, A):
        '''
        -------------------------------------------------------------------------
        yDBL()
        input : a projective Twisted Edwards y-coordinate point y(P) := YP/WP,
                and the projective Montgomery constants A24:= A + 2C and C24:=4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Twisted Edwards y-coordinate point y([2]P)
        -------------------------------------------------------------------------
        '''
        t_0 = self.fp.fp_sqr(P[0])
        t_1 = self.fp.fp_sqr(P[1])
        Z = self.fp.fp_mul(A[1], t_0)
        X = self.fp.fp_mul(Z, t_1)
        t_1 = self.fp.fp_sub(t_1, t_0)
        t_0 = self.fp.fp_mul(A[0], t_1)
        Z = self.fp.fp_add(Z, t_0)
        Z = self.fp.fp_mul(Z, t_1)
        return [ self.fp.fp_sub(X, Z), self.fp.fp_add(X, Z)]



    def yADD(self, P, Q, PQ):
        '''
        -------------------------------------------------------------------------
        yADD()
        input : the projective Twisted Edwards y-coordinate points y(P) := YP/WP,
                y(Q) := YQ/WQ, and y(P-Q) := YPQ/QPQ
        output: the projective Twisted Edwards y-coordinate point y(P+Q)
        -------------------------------------------------------------------------
        '''
        a = self.fp.fp_mul(P[1], Q[0])
        b = self.fp.fp_mul(P[0], Q[1])
        c = self.fp.fp_add(a, b)
        d = self.fp.fp_sub(a, b)
        c = self.fp.fp_sqr(c)
        d = self.fp.fp_sqr(d)

        xD = self.fp.fp_add(PQ[1], PQ[0])
        zD = self.fp.fp_sub(PQ[1], PQ[0])
        X = self.fp.fp_mul(zD, c)
        Z = self.fp.fp_mul(xD, d)
        return [ self.fp.fp_sub(X, Z), self.fp.fp_add(X, Z)]




    def KPs_t(self, P, A, i):
        '''
        -------------------------------------------------------------------------
        KPs()
        input : the projective Twisted Edwards y-coordinate points y(P) := YP/WP,
                the projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and a positive integer 0 <= i < n
        output: the list of projective Twisted Edwards y-coordinate points y(P),
                y([2]P), y([3]P), ..., and y([d_i]P) where l_i = 2 * d_i + 1
        -------------------------------------------------------------------------
        '''
        d = (self.fp.L[i] - 1) // 2

        self.K = [[0, 0] for j in range(d + 1)]
        self.K[0] = list([ self.fp.fp_sub(P[0], P[1]), self.fp.fp_add(P[0], P[1])])  # 2a
        self.K[1] = self.yDBL(self.K[0], A)  # 4M + 2S + 4a
        for j in range(2, d, 1):
            self.K[j] = self.yADD(self.K[j - 1], self.K[0], self.K[j - 2])  # 4M + 2S + 6a

        return self.K  # 2(l - 3)M + (l - 3)S + 3(l - 3)a




    def xISOG_t(self, A, i):
        '''
        ------------------------------------------------------------------------------
        xISOG()
        input : the projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, the list of projective Twisted 
                Edwards y-coordinate points y(P), y([2]P), y([3]P), ..., and y([d_i]P)
                where l_i = 2 * d_i + 1, and a positive integer 0 <= i < n
        output: the projective Montgomery constants a24:= a + 2c and c24:=4c where
                E': y^2 = x^3 + (a/c)*x^2 + x is a degree-(l_i) isogenous curve to E
        ------------------------------------------------------------------------------
        '''
        l = self.global_L[i]  # l
        bits_of_l = bitlength(l)  # Number of bits of L[i]
        d = (l - 1) // 2  # Here, l = 2d + 1

        By = self.K[0][0]
        Bz = self.K[0][1]
        for j in range(1, d, 1):

            By = self.fp.fp_mul(By, self.K[j][0])
            Bz = self.fp.fp_mul(Bz, self.K[j][1])

        bits_of_l -= 1
        constant_d_edwards = self.fp.fp_sub(A[0], A[1])  # 1a

        tmp_a = A[0]
        tmp_d = constant_d_edwards
        # left-to-right method for computing a^l and d^l
        for j in range(1, bits_of_l + 1):

            tmp_a = self.fp.fp_sqr(tmp_a)
            tmp_d = self.fp.fp_sqr(tmp_d)
            if ((l >> (bits_of_l - j)) & 1) != 0:

                tmp_a = self.fp.fp_mul(tmp_a, A[0])
                tmp_d = self.fp.fp_mul(tmp_d, constant_d_edwards)

        for j in range(3):

            By = self.fp.fp_sqr(By)
            Bz = self.fp.fp_sqr(Bz)

        C0 = self.fp.fp_mul(tmp_a, Bz)
        C1 = self.fp.fp_mul(tmp_d, By)
        C1 = self.fp.fp_sub(C0, C1)

        return [C0, C1]  # (l - 1 + 2*HW(l) - 2)M + 2(|l|_2 + 1)S + 2a



    def xEVAL_t(self, P, i):
        '''
        ------------------------------------------------------------------------------
        xEVAL()
        input : the projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, the list of projective Twisted
                Edwards y-coordinate points y(P), y([2]P), y([3]P), ..., and y([d_i]P)
                where l_i = 2 * d_i + 1, and a positive integer 0 <= i < n
        output: the projective Montgomery constants a24:= a + 2c and c24:=4c where
                E': y^2 = x^3 + (a/c)*x^2 + x is a degree-(l_i) isogenous curve to E
        ------------------------------------------------------------------------------
        '''
        d = (self.global_L[i] - 1) // 2  # Here, l = 2d + 1

        Q0 = self.fp.fp_add(P[0], P[1])
        Q1 = self.fp.fp_sub(P[0], P[1])
        R0, R1 = self.curve.CrissCross(self.K[0][1], self.K[0][0], Q0, Q1)
        for j in range(1, d, 1):

            T0, T1 = self.curve.CrissCross(self.K[j][1], self.K[j][0], Q0, Q1)
            R0 = self.fp.fp_mul(T0, R0)
            R1 = self.fp.fp_mul(T1, R1)

        R0 = self.fp.fp_sqr(R0)
        R1 = self.fp.fp_sqr(R1)
        X = self.fp.fp_mul(P[0], R0)
        Z = self.fp.fp_mul(P[1], R1)

        return [X, Z]  # 2(l - 1)M + 2S + (l + 1)a


    # Next functions is used for setting the cardinalities sI, sJ, and sK
    def set_parameters_velu(self, b, c, i):

        assert b <= c

        # At this step, everythin is correct
        self.sJ = b
        self.sI = c
        d = ((self.global_L[i] - 2 - 4 * b * c - 1) // 2) + 1
        assert d >= 0
        self.sK = d
        return None


    def print_parameters_velu(self):

        print("| sI: %3d, sJ: %3d, sK: %3d |" % (sI, sJ, sK), end="")
        return None


    # KPs computes x([i]P), x([j]P), and x([k]P) for each i in I, j in J, and j in J.
    # I, J, and K are defined according to Examples 4.7 and 4.12 of https://eprint.iacr.org/2020/341
    def KPs_s(self, P, A, i):

        # This functions computes all the independent data from the input of algorithm 2 of https://eprint.iacr.org/2020/341
        if self.sI == 0:

            # Case global_L[i] = 3 is super very special case (nothing to do)
            self.J = []
            self.ptree_hI = None
            self.K = [list(P)]
            # J, b, ptree_hI, c, K, d
            assert self.sJ == 0 and self.sI == 0 and self.sK == 1
            return None

        # We need to ensure sI is greater than or equal sJ
        assert self.sI >= self.sJ

        # Additionally, sK should be greater than or equal to zero. If that is not the case, then sJ and sI were badly chosen
        assert self.sK >= 0

        if self.sI == 1:
            # Branch corresponds with global_L[i] = 5 and global_L[i] = 7
            # Recall sJ > 0, then sJ = 1
            assert self.sJ == 1
            P2 = self.xDBL(P, A)

            self.J = [list(P)]

            I = [list(P2)]
            hI = [
                list([ self.fp.fp_sub(0, P2[0]), P2[1]])
            ]  # we only need to negate x-coordinate of each point
            self.ptree_hI = product_tree(hI, self.sI)  # product tree of hI

            if not self.SCALED_REMAINDER_TREE:
                # Using remainder trees
                self.ptree_hI = reciprocal_tree(
                    {'rpoly': [1], 'rdeg': 0, 'fpoly': [1], 'fdeg': 0, 'a': 1},
                    2 * self.sJ + 1,
                    self.ptree_hI,
                    self.sI,
                )  # reciprocal of the root is used to compute sons' reciprocals

            else:
                # Using scaled remainder trees
                assert (2 * self.sJ - self.sI + 1) > self.sI
                self.ptree_hI['reciprocal'], self.ptree_hI['a'] = reciprocal(
                    self.ptree_hI['poly'][::-1], self.sI + 1, 2 * self.sJ - self.sI + 1
                )
                self.ptree_hI['scaled'], self.ptree_hI['as'] = (
                    list(self.ptree_hI['reciprocal'][:self.sI]),
                    self.ptree_hI['a'],
                )

            # Notice, 0 <= sK <= 1
            assert self.sK <= 1
            if self.sK == 1:
                self.K = [list(P2)]
            else:
                self.K = []

            return None

        # At this step, sI > 1
        assert self.sI > 1
        if self.sJ == 1:
            # This branch corresponds with global_L[i] = 11 and global_L[i] = 13
            Q = self.xDBL(P, A)  # x([2]P)
            Q2 = self.xDBL(Q, A)  # x([2]Q)

            self.J = [list(P)]

            I = [[0, 0]] * sI
            I[0] = list(Q)  # x(   Q)
            I[1] = self.xADD(Q2, I[0], I[0])  # x([3]Q)
            for ii in range(2, self.sI, 1):
                I[ii] = self.xADD(I[ii - 1], Q2, I[ii - 2])  # x([2**i + 1]Q)

            hI = [
                [ self.fp.fp_sub(0, iP[0]), iP[1]] for iP in I
            ]  # we only need to negate x-coordinate of each point
            self.ptree_hI = product_tree(hI, self.sI)  # product tree of hI

            if not self.SCALED_REMAINDER_TREE:
                # Using remainder trees
                self.ptree_hI = reciprocal_tree(
                    {'rpoly': [1], 'rdeg': 0, 'fpoly': [1], 'fdeg': 0, 'a': 1},
                    2 * self.sJ + 1,
                    self.ptree_hI,
                    self.sI,
                )  # reciprocal of the root is used to compute sons' reciprocals

            else:
                # Using scaled remainder trees
                assert (2 * self.sJ - self.sI + 1) <= self.sI
                self.ptree_hI['scaled'], ptree_hI['as'] = reciprocal(
                    ptree_hI['poly'][::-1], self.sI + 1, self.sI
                )
                self.ptree_hI['reciprocal'], self.ptree_hI['a'] = (
                    list(self.ptree_hI['scaled'][: (2 * self.sJ - self.sI + 1)]),
                    self.ptree_hI['as'],
                )

            # Notice, 0 <= sK <= 1
            assert self.sK <= 1
            if self.sK == 1:
                self.K = [list(Q)]
            else:
                self.K = []

            return None

        # Now, we ensure sI >= sJ > 1
        assert self.sJ > 1

        # In other words, we can proceed by the general case

        # ------------------------------------------------
        # Computing [j]P for each j in {1, 3, ..., 2*sJ - 1}
        self.J = [[0, 0]] * self.sJ
        self.J[0] = list(P)  # x(   P)
        P2 = self.curve.xDBL(P, A)  # x([2]P)
        self.J[1] = self.curve.xADD(P2, self.J[0], self.J[0])  # x([3]P)
        for jj in range(2, self.sJ, 1):
            self.J[jj] = self.curve.xADD(self.J[jj - 1], P2, self.J[jj - 2])  # x([2*jj + 1]P)

        # -------------------------------------------------------
        # Computing [i]P for i in { (2*sJ) * (2i + 1) : 0 <= i < sI}
        bhalf_floor = self.sJ // 2
        bhalf_ceil = self.sJ - bhalf_floor
        P4 = self.curve.xDBL(P2, A)  # x([4]P)
        P2[0], P4[0] = self.fp.fp_cswap(P2[0], P4[0], self.sJ % 2)  # Constant-time swap
        P2[1], P4[1] = self.fp.fp_cswap(
            P2[1], P4[1], self.sJ % 2
        )  # x([4]P) <--- coditional swap ---> x([2]P)
        Q = self.curve.xADD(self.J[bhalf_ceil], self.J[bhalf_floor - 1], P2)  # Q := [2b]P
        P2[0], P4[0] = self.fp.fp_cswap(P2[0], P4[0], self.sJ % 2)  # Constant-time swap
        P2[1], P4[1] = self.fp.fp_cswap(
            P2[1], P4[1], self.sJ % 2
        )  # x([4]P) <--- coditional swap ---> x([2]P)

        I = [[0, 0]] * self.sI
        I[0] = list(Q)  # x(   Q)
        Q2 = self.curve.xDBL(Q, A)  # x([2]Q)
        I[1] = self.curve.xADD(Q2, I[0], I[0])  # x([3]Q)
        for ii in range(2, self.sI, 1):
            I[ii] = self.curve.xADD(I[ii - 1], Q2, I[ii - 2])  # x([2**i + 1]Q)

        # --------------------------------------------------------------
        # Computing [k]P for k in { 4*sJ*sI + 1, ..., l - 6, l - 4, l - 2}
        self.K = [[0, 0]] * self.sK

        if self.sK >= 1:
            self.K[0] = list(P2)  # x([l - 2]P) = x([-2]P) = x([2]P)

        if self.sK >= 2:
            self.K[1] = list(P4)  # x([l - 4]P) = x([-4]P) = x([4]P)

        for k in range(2, self.sK, 1):
            self.K[k] = self.curve.xADD(self.K[k - 1], P2, self.K[k - 2])

        # ------------------------------------------------------------------------------------------------------
        #                   ~~~~~~~~               ~~~~~~~~
        #                    |    |                 |    |
        # Computing h_I(W) = |    | (W - x([i]P)) = |    | (Zi * W - Xi) / Zi where x([i]P) = Xi/Zi
        #                    i in I                 i in I
        # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

        hI = [
            [ self.fp.fp_sub(0, iP[0]), iP[1]] for iP in I
        ]  # we only need to negate x-coordinate of each point
        self.ptree_hI = product_tree(hI, self.sI)  # product tree of hI

        if not self.SCALED_REMAINDER_TREE:
            # Using scaled remainder trees
            self.ptree_hI = reciprocal_tree(
                {'rpoly': [1], 'rdeg': 0, 'fpoly': [1], 'fdeg': 0, 'a': 1},
                2 * sJ + 1,
                self.ptree_hI,
                self.sI,
            )  # reciprocal of the root is used to compute sons' reciprocals

        else:
            # Using scaled remainder trees
            if self.sI < (2 * self.sJ - self.sI + 1):
                self.ptree_hI['reciprocal'], self.ptree_hI['a'] = reciprocal(
                    self.ptree_hI['poly'][::-1], self.sI + 1, 2 * self.sJ - self.sI + 1
                )
                self.ptree_hI['scaled'], self.ptree_hI['as'] = (
                    list(self.ptree_hI['reciprocal'][:self.sI]),
                    self.ptree_hI['a'],
                )

            else:
                self.ptree_hI['scaled'], self.ptree_hI['as'] = reciprocal(
                    self.ptree_hI['poly'][::-1], self.sI + 1, self.sI
                )
                self.ptree_hI['reciprocal'], self.ptree_hI['a'] = (
                    list(self.ptree_hI['scaled'][: (2 * self.sJ - self.sI + 1)]),
                    self.ptree_hI['as'],
                )

        # Polynomial D_j of algorithm 2 from https://eprint.iacr.org/2020/341 is not required,
        # but we need some some squares and products determined by list J
        # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

        # Ensuring the cardinality of each ser coincide with the expected one
        assert len(I) == self.sI
        assert len(self.J) == self.sJ
        assert len(self.K) == self.sK

        return None


    # Next function perform algorithm 2 of https://eprint.iacr.org/2020/341 with input \alpha = 1 and \alpha = -1, and
    # then it computes the isogenous Montgomery curve coefficient
    def xISOG_s(self, A, i):

        AA = self.fp.fp_add(A[0], A[0])  # 2A' + 4C
        AA = self.fp.fp_sub(AA, A[1])  # 2A'
        AA = self.fp.fp_add(AA, AA)  # 4A'

        # Polynomial D_j of algorithm 2 from https://eprint.iacr.org/2020/341 is not required,
        # but we need some some squares and products determined by list J
        # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

        SUB_SQUARED = [0 for j in range(0, self.sJ, 1)]  #
        ADD_SQUARED = [0 for j in range(0, self.sJ, 1)]  #

        # List XZJ4 is required for degree-l isogeny evaluations...
        self.XZJ4 = [0 for j in range(0, self.sJ, 1)]  # 2*Xj*Zj
        for j in range(0, self.sJ, 1):

            SUB_SQUARED[j] = self.fp.fp_sub(self.J[j][0], self.J[j][1])  # (Xj - Zj)
            SUB_SQUARED[j] = self.fp.fp_sqr(SUB_SQUARED[j])  # (Xj - Zj)^2

            ADD_SQUARED[j] = self.fp.fp_add(self.J[j][0], self.J[j][1])  # (Xj + Zj)
            ADD_SQUARED[j] = self.fp.fp_sqr(ADD_SQUARED[j])  # (Xj + Zj)^2

            self.XZJ4[j] = self.fp.fp_sub(SUB_SQUARED[j], ADD_SQUARED[j])  # -4*Xj*Zj

        # --------------------------------------------------------------------------------------------------
        #                   ~~~~~~~~
        #                    |    |
        # Computing E_J(W) = |    | [ F0(W, x([j]P)) * alpha^2 + F1(W, x([j]P)) * alpha + F2(W, x([j]P)) ]
        #                    j in J
        # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates
        # In particular, for a degree-l isogeny construction, we need alpha = 1 and alpha = -1

        # EJ_0 is the one determined by alpha = 1
        EJ_0 = [[0, 0, 0] for j in range(0, self.sJ, 1)]
        # EJ_1 is the one determined by alpha = -1
        EJ_1 = [[0, 0, 0] for j in range(0, self.sJ, 1)]

        for j in range(0, self.sJ, 1):

            # However, each SUB_SQUARED[j] and ADD_SQUARED[j] should be multiplied by C
            tadd = self.fp.fp_mul(ADD_SQUARED[j], A[1])
            tsub = self.fp.fp_mul(SUB_SQUARED[j], A[1])

            # We require the double of tadd and tsub
            tadd2 = self.fp.fp_add(tadd, tadd)
            tsub2 = self.fp.fp_add(tsub, tsub)

            t1 = self.fp.fp_mul(self.XZJ4[j], AA)  #       A *(-4*Xj*Zj)

            # Case alpha = 1
            linear = self.fp.fp_sub(
                t1, tadd2
            )  #       A *(-4*Xj*Zj)  - C * (2 * (Xj + Zj)^2)
            EJ_0[j] = [tsub, linear, tsub]

            # Case alpha = -1
            linear = self.fp.fp_sub(
                tsub2, t1
            )  #       C * (2 * (Xj - Zj)^2) - A *(-4*Xj*Zj)
            EJ_1[j] = [tadd, linear, tadd]

        # The faster way for multiplying is using a divide-and-conquer approach
        poly_EJ_0 = product_selfreciprocal_tree(EJ_0, self.sJ)[
            'poly'
        ]  # product tree of EJ_0 (we only require the root)
        poly_EJ_1 = product_selfreciprocal_tree(EJ_1, self.sJ)[
            'poly'
        ]  # product tree of EJ_1 (we only require the root)

        if not self.SCALED_REMAINDER_TREE:
            # Remainder tree computation
            remainders_EJ_0 = multieval_unscaled(
                poly_EJ_0, 2 * self.sJ + 1, self.ptree_hI, self.sI
            )
            remainders_EJ_1 = multieval_unscaled(
                poly_EJ_1, 2 * self.sJ + 1, self.ptree_hI, self.sI
            )

        else:
            # Approach using scaled remainder trees
            if self.ptree_hI != None:
                poly_EJ_0 = poly_redc(poly_EJ_0, 2 * self.sJ + 1, self.ptree_hI)
                fg_0 = poly_mul_middle(self.ptree_hI['scaled'], self.sI, poly_EJ_0[::-1], self.sI)
                remainders_EJ_0 = multieval_scaled(
                    fg_0[::-1], self.sI, [1] + [0] * (self.sI - 1), self.sI, self.ptree_hI, self.sI
                )

                poly_EJ_1 = poly_redc(poly_EJ_1, 2 * self.sJ + 1, self.ptree_hI)
                fg_1 = poly_mul_middle(self.ptree_hI['scaled'], self.sI, poly_EJ_1[::-1], self.sI)
                remainders_EJ_1 = multieval_scaled(
                    fg_1[::-1], self.sI, [1] + [0] * (self.sI - 1), self.sI, self.ptree_hI, self.sI
                )
            else:
                remainders_EJ_0 = []
                remainders_EJ_1 = []

        # Multipying all the remainders
        r0 = product(remainders_EJ_0, self.sI)
        r1 = product(remainders_EJ_1, self.sI)

        # ---------------------------------------------------------------------------------
        # Now, we proceed by computing the missing part which is determined by K
        # Notice, the denominators are the same and then they annulled between them
        # In other words, it is not required to compute the product of all Zk's with k In K

        # Case alpha = 1
        hK_0 = [[ self.fp.fp_sub(self.K[k][1], self.K[k][0])] for k in range(0, self.sK, 1)]
        hK_0 = product(hK_0, self.sK)  # product of (Zk - Xk) for each k in K
        # Case alpha = -1
        hK_1 = [[ self.fp.fp_add(self.K[k][1], self.K[k][0])] for k in range(0, self.sK, 1)]
        hK_1 = product(hK_1, self.sK)  # product of (Zk + Xk) for each k in K

        # --------------------------------------------------------------
        # Now, we have all the ingredients for computing the image curve
        A24m = self.fp.fp_sub(A[0], A[1])  # A' - 2C

        A24 = self.fp.fp_exp(A[0], self.global_L[i])  # (A' + 2C)^l
        A24m = self.fp.fp_exp(A24m, self.global_L[i])  # (A' - 2C)^l

        t24m = self.fp.fp_mul(
            hK_1, r1
        )  # output of algorithm 2 with alpha =-1 and without the demoninator
        t24m = self.fp.fp_sqr(t24m)  # raised at 2
        t24m = self.fp.fp_sqr(t24m)  # raised at 4
        t24m = self.fp.fp_sqr(t24m)  # raised at 8

        t24 = self.fp.fp_mul(
            hK_0, r0
        )  # output of algorithm 2 with alpha = 1 and without the demoninator
        t24 = self.fp.fp_sqr(t24)  # raised at 2
        t24 = self.fp.fp_sqr(t24)  # raised at 4
        t24 = self.fp.fp_sqr(t24)  # raised at 8

        A24 = self.fp.fp_mul(A24, t24m)
        A24m = self.fp.fp_mul(A24m, t24)

        # Now, we have d = (A24m / A24) where the image Montgomery cuve coefficient is
        #      B'   2*(1 + d)   2*(A24 + A24m)
        # B = ---- = --------- = --------------
        #      C      (1 - d)     (A24 - A24m)
        # However, we required B' + 2C = 4*A24 and 4C = 4 * (A24 - A24m)

        t24m = self.fp.fp_sub(A24, A24m)  #   (A24 - A24m)
        t24m = self.fp.fp_add(t24m, t24m)  # 2*(A24 - A24m)
        t24m = self.fp.fp_add(t24m, t24m)  # 4*(A24 - A24m)

        t24 = self.fp.fp_add(A24, A24)  # 2 * A24
        t24 = self.fp.fp_add(t24, t24)  # 4 * A24

        # return [t24, t24m], ptree_hI, XZJ4
        return [t24, t24m]


    def xEVAL_s(self, P, A):

        AA = self.fp.fp_add(A[0], A[0])  # 2A' + 4C
        AA = self.fp.fp_sub(AA, A[1])  # 2A'
        AA = self.fp.fp_add(AA, AA)  # 4A'

        # --------------------------------------------------------------------------------------------------
        #                   ~~~~~~~~
        #                    |    |
        # Computing E_J(W) = |    | [ F0(W, x([j]P)) * alpha^2 + F1(W, x([j]P)) * alpha + F2(W, x([j]P)) ]
        #                    j in J
        # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates
        # In particular, for a degree-l isogeny construction, we need alpha = 1 and alpha = -1

        # EJ_0 is the one determined by alpha = x
        EJ_0 = [[0, 0, 0] for j in range(0, self.sJ, 1)]
        # Notice, the corresponding EJ_1 that is determined by alpha = 1/x can be computed by using EJ_0

        XZ_add = self.fp.fp_add(P[0], P[1])  # X + Z
        XZ_sub = self.fp.fp_sub(P[0], P[1])  # X - Z

        AXZ2 = self.fp.fp_mul(P[0], P[1])  # X * Z
        t1 = self.fp.fp_sqr(P[0])  # X^2
        t2 = self.fp.fp_sqr(P[1])  # Z^2

        CX2Z2 = self.fp.fp_add(t1, t2)  #      X^2 + Z^2
        CX2Z2 = self.fp.fp_mul(CX2Z2, A[1])  # C * (X^2 + Z^2)

        AXZ2 = self.fp.fp_add(AXZ2, AXZ2)  #       2 * (X * Z)
        CXZ2 = self.fp.fp_mul(AXZ2, A[1])  # C  * [2 * (X * Z)]
        AXZ2 = self.fp.fp_mul(AXZ2, AA)  # A' * [2 * (X * Z)]

        for j in range(0, self.sJ, 1):

            XZj_add = self.fp.fp_add(self.J[j][0], self.J[j][1])  # Xj + Zj
            XZj_sub = self.fp.fp_sub(self.J[j][0], self.J[j][1])  # Xj - Zj

            t1 = self.fp.fp_mul(XZ_sub, XZj_add)  # (X - Z) * (Xj + Zj)
            t2 = self.fp.fp_mul(XZ_add, XZj_sub)  # (X + Z) * (Xj - Zj)

            # Computing the quadratic coefficient
            quadratic = self.fp.fp_sub(t1, t2)  #   2 * [(X*Zj) - (Z*Xj)]
            quadratic = self.fp.fp_sqr(quadratic)  # ( 2 * [(X*Zj) - (Z*Xj)] )^2
            quadratic = self.fp.fp_mul(A[1], quadratic)  # C * ( 2 * [(X*Zj) - (Z*Xj)] )^2

            # Computing the constant coefficient
            constant = self.fp.fp_add(t1, t2)  #   2 * [(X*Xj) - (Z*Zj)]
            constant = self.fp.fp_sqr(constant)  # ( 2 * [(X*Xj) - (Z*Zj)] )^2
            constant = self.fp.fp_mul(A[1], constant)  # C * ( 2 * [(X*Xj) - (Z*Zj)] )^2

            # Computing the linear coefficient
            # ----------------------------------------------------------------------------------------------------------
            # C * [ (-2*Xj*Zj)*(alpha^2 + 1) + (-2*alpha)*(Xj^2 + Zj^2)] + [A' * (-2*Xj*Zj) * (2*X*Z)] where alpha = X/Z
            t1 = self.fp.fp_add(self.J[j][0], self.J[j][1])  #     (Xj + Zj)
            t1 = self.fp.fp_sqr(t1)  #     (Xj + Zj)^2
            t1 = self.fp.fp_add(t1, t1)  # 2 * (Xj + Zj)^2
            t1 = self.fp.fp_add(
                t1, self.XZJ4[j]
            )  # 2 * (Xj + Zj)^2 - (4*Xj*Zj) := 2 * (Xj^2 + Zj^2)
            t1 = self.fp.fp_mul(t1, CXZ2)  # [2 * (Xj^2 + Zj^2)] * (2 * [ C * (X * Z)])

            t2 = self.fp.fp_mul(CX2Z2, self.XZJ4[j])  # [C * (X^2 + Z^2)] * (-4 * Xj * Zj)
            t1 = self.fp.fp_sub(
                t2, t1
            )  # [C * (X^2 + Z^2)] * (-4 * Xj * Zj) - [2 * (Xj^2 + Zj^2)] * (2 * [ C * (X * Z)])

            t2 = self.fp.fp_mul(AXZ2, self.XZJ4[j])  # (2 * [A' * (X * Z)]) * (-4 * Xj * Zj)
            linear = self.fp.fp_add(
                t1, t2
            )  # This is our desired equation but multiplied by 2
            linear = self.fp.fp_add(
                linear, linear
            )  # This is our desired equation but multiplied by 4
            # ----------------------------------------------------------------------------------------------------------

            # Case alpha = X / Z
            EJ_0[j] = [constant, linear, quadratic]

        # The faster way for multiplying is using a divide-and-conquer approach
        poly_EJ_0 = product_tree(EJ_0, self.sJ)[
            'poly'
        ]  # product tree of EJ_0 (we only require the root)
        poly_EJ_1 = list(
            poly_EJ_0[::-1]
        )  # product tree of EJ_1(x) = x^{2b + 1} EJ_0(1/X)

        if not self.SCALED_REMAINDER_TREE:
            # Remainder tree computation
            remainders_EJ_0 = multieval_unscaled(
                poly_EJ_0, 2 * self.sJ + 1, self.ptree_hI, self.sI
            )
            remainders_EJ_1 = multieval_unscaled(
                poly_EJ_1, 2 * self.sJ + 1, self.ptree_hI, self.sI
            )

        else:
            # Approach using scaled remainder trees
            if self.ptree_hI != None:
                poly_EJ_0 = poly_redc(poly_EJ_0, 2 * self.sJ + 1, self.ptree_hI)
                fg_0 = poly_mul_middle(self.ptree_hI['scaled'], self.sI, poly_EJ_0[::-1], self.sI)
                remainders_EJ_0 = multieval_scaled(
                    fg_0[::-1], self.sI, [1] + [0] * (self.sI - 1), self.sI, self.ptree_hI, self.sI
                )

                poly_EJ_1 = poly_redc(poly_EJ_1, 2 * self.sJ + 1, self.ptree_hI)
                fg_1 = poly_mul_middle(self.ptree_hI['scaled'], self.sI, poly_EJ_1[::-1], self.sI)
                remainders_EJ_1 = multieval_scaled(
                    fg_1[::-1], self.sI, [1] + [0] * (self.sI - 1), self.sI, self.ptree_hI, self.sI
                )
            else:
                remainders_EJ_0 = []
                remainders_EJ_1 = []

        # Multipying all the remainders
        r0 = product(remainders_EJ_0, self.sI)
        r1 = product(remainders_EJ_1, self.sI)

        # ---------------------------------------------------------------------------------
        # Now, we proceed by computing the missing part which is determined by K
        # Notice, the denominators are the same and then they annulled between them
        # In other words, it is not required to compute the product of all Zk's with k In K

        hK_0 = [[0]] * self.sK
        hK_1 = [[0]] * self.sK
        for k in range(0, self.sK, 1):

            XZk_add = self.fp.fp_add(self.K[k][0], self.K[k][1])  # Xk + Zk
            XZk_sub = self.fp.fp_sub(self.K[k][0], self.K[k][1])  # Xk - Zk
            t1 = self.fp.fp_mul(XZ_sub, XZk_add)  # (X - Z) * (Xk + Zk)
            t2 = self.fp.fp_mul(XZ_add, XZk_sub)  # (X + Z) * (Xk - Zk)

            # Case alpha = X/Z
            hK_0[k] = [ self.fp.fp_sub(t1, t2)]  # 2 * [(X*Zk) - (Z*Xk)]

            # Case 1/alpha = Z/X
            hK_1[k] = [ self.fp.fp_add(t1, t2)]  # 2 * [(X*Xk) - (Z*Zk)]

        hK_0 = product(hK_0, self.sK)  # product of (XZk - ZXk) for each k in K
        hK_1 = product(hK_1, self.sK)  # product of (XXk - ZZk) for each k in K

        # ---------------------------------------------------------------------------------
        # Now, unifying all the computations
        XX = self.fp.fp_mul(
            hK_1, r1
        )  # output of algorithm 2 with 1/alpha = Z/X and without the demoninator
        XX = self.fp.fp_sqr(XX)
        XX = self.fp.fp_mul(XX, P[0])

        ZZ = self.fp.fp_mul(
            hK_0, r0
        )  # output of algorithm 2 with alpha = X/Z and without the demoninator
        ZZ = self.fp.fp_sqr(ZZ)
        ZZ = self.fp.fp_mul(ZZ, P[1])

        return [XX, ZZ]


    def KPs(self, P, A, i):

        if self.global_L[i] <= self.HYBRID_BOUND:
            return self.KPs_t(P, A, i)
        else:
            return self.KPs_s(P, A, i)


    def xISOG(self, A, i):

        if self.global_L[i] <= self.HYBRID_BOUND:
            return self.xISOG_t(A, i)
        else:
            return self.xISOG_s(A, i)


    def xEVAL(self, P, v):

        if type(v) == int:
            return self.xEVAL_t(P, v)
        else:
            return self.xEVAL_s(P, v)


    def cISOG_and_cEVAL(self):

        if self.verbose:

            self.sI_list = []
            self.sJ_list = []
            f = open(ijk_data + self.prime)

            for i in range(0, self.fp.n, 1):

                bc = f.readline()
                bc = [int(bci) for bci in bc.split()]
                sJ_list.append(bc[0])
                sI_list.append(bc[1])

            f.close()

        # First, we look for a full torsion point
        A = [2, 4]
        T_p, T_m = self.curve.full_torsion_points(A)

        for i in range(0, self.fp.n, 1):

            if self.verbose:
                self.set_parameters_velu(sJ_list[i], sI_list[i], i)
            else:
                # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                # These paramters are required in KPs, xISOG, and xEVAL
                if self.global_L[i] == 3:
                    b = 0
                    c = 0
                else:
                    b = int(floor(sqrt(self.global_L[i] - 1) / 2.0))
                    c = int(floor((self.global_L[i] - 1.0) / (4.0 * b)))

                self.set_parameters_velu(b, c, i)

            # Getting an orderl-l point
            Tp = list(T_p)
            for j in range(i + 1, self.fp.n, 1):
                Tp = self.curve.xMUL(Tp, A, j)

            # Cost of xISOG() and KPs()
            self.fp.set_zero_ops()
            self.KPs(Tp, A, i)
            t = self.fp.get_ops()
            self.C_xISOG[i] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            self.fp.set_zero_ops()
            B = self.xISOG(A, i)
            t = self.fp.get_ops()
            self.C_xISOG[i] += numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            # xEVAL: kernel point determined by the next isogeny evaluation
            self.fp.set_zero_ops()
            if self.global_L[i] <= self.HYBRID_BOUND:
                T_p = self.xEVAL(T_p, i)
            else:
                T_p = self.xEVAL(T_p, A)

            # Cost of xEVAL
            self.fp.set_zero_ops()
            if self.global_L[i] <= self.HYBRID_BOUND:
                T_m = self.xEVAL(T_m, i)
            else:
                T_m = self.xEVAL(T_m, A)

            t = self.fp.get_ops()
            self.C_xEVAL[i] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            # Updating the new next curve
            A = list(B)

        self.fp.set_zero_ops()
        assert self.curve.coeff(A) == 0x6
        return None
