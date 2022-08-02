from random import SystemRandom
from pkg_resources import resource_filename

# Related polynomial multiplication procedures
from sibc.polymul import PolyMul
# Related polynomial reduction procedures
from sibc.polyredc import PolyRedc

from sibc.math import isequal, bitlength, hamming_weight, cswap
from sibc.constants import parameters

import numpy
from math import floor, sqrt

def doc(s):
    class __doc(object):
        def __init__(self,f):
            self.func = f
            self.desc = s
        def __call__(self,*args,**kwargs):
            return self.func(*args,**kwargs)
        def __repr__(self):
            return self.desc
    return __doc

def MontgomeryIsogeny(name : str, uninitialized = False):

    cutoff = 83
    cutoff_string = f' with cutoff ell <= {cutoff}' if name == 'hvelu' else ''
    NAME = 'Isogeny class using the %s Velu\'s formulae%s' % ({'tvelu':'tradicional', 'svelu':'square-root', 'hvelu':'hybrid'}[name], cutoff_string)
    @doc(NAME)
    class Formulae():

        def __init__(self, curve, tuned, multievaluation):
            self.name = name
            self.multievaluation_name = {True:'scaled', False:'unscaled'}[multievaluation]
            self.tuned_name = {True:'-tuned', False:''}[tuned]
            self.uninitialized = uninitialized

            # Get cost of the isogeny constructions and evaluations
            if tuned:
                # Reading tuned velusqrt parameters from the stored data
                self.sI_list = []
                self.sJ_list = []

                if curve.name in parameters['csidh'].keys():
                    algorithm = '/csidh'
                elif curve.name in parameters['bsidh'].keys():
                    algorithm = '/bsidh'
                else:
                    assert False, "only CSIDH and B-SIDH are currently implemented"

                path = resource_filename('sibc', "data/ijk/" + curve.model + algorithm + '/' + curve.name + '-' + self.multievaluation_name)
                f = open(path)

                for i in range(0, curve.n, 1):

                    bc = f.readline()
                    bc = [int(bci) for bci in bc.split()]
                    self.sJ_list.append(bc[0])
                    self.sI_list.append(bc[1])

                f.close()

            else:
                self.sI_list = None
                self.sJ_list = None

            # Tradicional velu formulas will be used for l_i <= self.HYBRID_BOUND
            self.HYBRID_BOUND = {'tvelu':max(curve.L), 'svelu':1, 'hvelu':cutoff}[name]

            # Global variables to be used in kps, xisog, and xeval

            # Here, J is a set of cardinality sJ
            self.J = None
            self.sJ = None

            # Here, ptree_I corresponds with the product tree determined by I, and I is a set of cardinality sJ
            self.ptree_hI = None
            self.sI = (None,)

            # Here, K is a set of cardinality sK
            self.K = None
            self.sK = None

            # An extra global variable which is used in xisog and xeval
            self.XZJ4 = None
            self.ADD_SQUARED = None
            self.SUB_SQUARED = None

            self.SCALED_MULTIEVALUATION = multievaluation
            self.tuned = tuned

            self.prime = curve.prime
            self.curve = curve
            self.field = self.curve.field
            self.L = self.curve.L

            self.poly_mul = PolyMul(self.field)
            self.poly_redc = PolyRedc(self.poly_mul)

            self.c_xeval = list(
                map(self.ceval, self.L)
            )  # list of the costs of each degree-l isogeny evaluation
            self.c_xisog = list(
                map(self.cisog, self.L)
            )  # list of the costs of each degree-l isogeny construction

            # Now, we proceed to store all the correct costs
            if name != 'tvelu' and uninitialized:
                print("// Precomputation cost regarding kps, xisog, and xeval of the velusqrt formulae")
                self.velusqrt_cost()

        def ceval(self, l):
            return numpy.array([2.0 * (l - 1.0), 2.0, (l + 1.0)])

        def cisog(self, l):
            return numpy.array(
                [
                    (
                        3.0 * l
                        + 2.0 * hamming_weight(l)
                        - 9.0
                        + isequal[l == 3] * 4.0
                    ),
                    (l + 2.0 * bitlength(l) + 1.0 + isequal[l == 3] * 2.0),
                    (3.0 * l - 7.0 + isequal[l == 3] * 6.0),
                ]
            )

        def ydbl(self, P, A):
            '''
            -------------------------------------------------------------------------
            ydbl()
            input : a projective Twisted Edwards y-coordinate point y(P) := YP/WP,
                    and the projective Montgomery constants A24:= A + 2C and C24:=4C
                    where E : y^2 = x^3 + (A/C)*x^2 + x
            output: the projective Twisted Edwards y-coordinate point y([2]P)
            -------------------------------------------------------------------------
            '''
            t_0 = (P[0] ** 2)
            t_1 = (P[1] ** 2)
            Z = (A[1] * t_0)
            X = (Z * t_1)
            t_1 -= t_0
            t_0 = (A[0] * t_1)
            Z += t_0
            Z *= t_1
            t_3 = X -Z # X - Z
            X += Z # X +Z
            return [t_3, X]

        def yadd(self, P, Q, PQ):
            '''
            -------------------------------------------------------------------------
            yadd()
            input : the projective Twisted Edwards y-coordinate points y(P) := YP/WP,
                    y(Q) := YQ/WQ, and y(P-Q) := YPQ/QPQ
            output: the projective Twisted Edwards y-coordinate point y(P+Q)
            -------------------------------------------------------------------------
            '''
            a = (P[1] * Q[0])
            b = (P[0] * Q[1])
            c = (a + b)
            a -= b # d = (a - b)
            c **= 2
            a **= 2 # d = d ** 2

            xD = (PQ[1] + PQ[0])
            zD = (PQ[1] - PQ[0])
            c *= zD # X
            Z = (xD * a) # xD * (a-b)**2
            t_3 = c - Z # X -Z
            c += Z # X + Z
            return [t_3, c]

        def kps_t(self, P, A, i):
            '''
            -------------------------------------------------------------------------
            kps_t()
            input : the projective Twisted Edwards y-coordinate points y(P) := YP/WP,
                    the projective Montgomery constants A24:= A + 2C and C24:=4C where
                    E : y^2 = x^3 + (A/C)*x^2 + x, and a positive integer 0 <= i < n
            output: the list of projective Twisted Edwards y-coordinate points y(P),
                    y([2]P), y([3]P), ..., and y([d_i]P) where l_i = 2 * d_i + 1
            -------------------------------------------------------------------------
            '''
            d = (self.curve.L[i] - 1) // 2

            self.K = [[0, 0] for j in range(d + 1)]
            self.K[0] = list(
                [(P[0] - P[1]), (P[0] + P[1])]
            )  # 2a
            self.K[1] = self.ydbl(self.K[0], A)  # 4M + 2S + 4a
            for j in range(2, d, 1):
                self.K[j] = self.yadd(
                    self.K[j - 1], self.K[0], self.K[j - 2]
                )  # 4M + 2S + 6a

            return self.K  # 2(l - 3)M + (l - 3)S + 3(l - 3)a

        def xisog_t(self, A, i):
            '''
            ------------------------------------------------------------------------------
            xisog_t()
            input : the projective Montgomery constants A24:= A + 2C and C24:=4C where
                    E : y^2 = x^3 + (A/C)*x^2 + x, the list of projective Twisted 
                    Edwards y-coordinate points y(P), y([2]P), y([3]P), ..., and y([d_i]P)
                    where l_i = 2 * d_i + 1, and a positive integer 0 <= i < n
            output: the projective Montgomery constants a24:= a + 2c and c24:=4c where
                    E': y^2 = x^3 + (a/c)*x^2 + x is a degree-(l_i) isogenous curve to E
            ------------------------------------------------------------------------------
            '''
            l = self.L[i]  # l
            bits_of_l = bitlength(l)  # Number of bits of L[i]
            d = (l - 1) // 2  # Here, l = 2d + 1

            # +0 to get a new FF allocation we can mutate inplace in
            # the for-loop below:
            By = self.K[0][0].copy()
            Bz = self.K[0][1].copy()
            for j in range(1, d, 1):

                By *= self.K[j][0] # By *= self.K[j][0]
                Bz *= self.K[j][1] # Bz *= self.K[j][1]

            bits_of_l -= 1
            constant_d_edwards = (A[0] - A[1])  # 1a

            tmp_a = A[0].copy()
            tmp_d = constant_d_edwards.copy()
            # left-to-right method for computing a^l and d^l
            for j in range(1, bits_of_l + 1):

                tmp_a **= 2
                tmp_d **= 2
                if ((l >> (bits_of_l - j)) & 1) != 0:

                    tmp_a *= A[0]
                    tmp_d *= constant_d_edwards

            for j in range(3):

                By **= 2
                Bz **= 2

            tmp_a *= Bz # C0 = (tmp_a * Bz)
            tmp_d *= By # C1
            C1 = (tmp_a - tmp_d) # C1 = C0 - C1

            return [tmp_a, C1]  # (l - 1 + 2*HW(l) - 2)M + 2(|l|_2 + 1)S + 2a

        def xeval_t(self, P, i):
            '''
            ------------------------------------------------------------------------------
            xeval_t()
            input : the projective Montgomery constants A24:= A + 2C and C24:=4C where
                    E : y^2 = x^3 + (A/C)*x^2 + x, the list of projective Twisted
                    Edwards y-coordinate points y(P), y([2]P), y([3]P), ..., and y([d_i]P)
                    where l_i = 2 * d_i + 1, and a positive integer 0 <= i < n
            output: the projective Montgomery constants a24:= a + 2c and c24:=4c where
                    E': y^2 = x^3 + (a/c)*x^2 + x is a degree-(l_i) isogenous curve to E
            ------------------------------------------------------------------------------
            '''
            d = (self.L[i] - 1) // 2  # Here, l = 2d + 1

            Q0 = (P[0] + P[1])
            Q1 = (P[0] - P[1])
            R0, R1 = self.curve.crisscross(self.K[0][1], self.K[0][0], Q0, Q1)
            for j in range(1, d, 1):

                T0, T1 = self.curve.crisscross(self.K[j][1], self.K[j][0], Q0, Q1)
                R0 *= T0
                R1 *= T1

            R0 **= 2
            R1 **= 2
            R0 *= P[0] # X
            R1 *= P[1] # Z

            return [R0, R1]  # 2(l - 1)M + 2S + (l + 1)a

        # Next functions is used for setting the cardinalities sI, sJ, and sK
        def set_parameters_velu(self, b, c, i):

            assert b <= c

            # At this step, everythin is correct
            self.sJ = b
            self.sI = c
            d = ((self.L[i] - 2 - 4 * b * c - 1) // 2) + 1
            assert d >= 0
            self.sK = d
            return None

        def print_parameters_velu(self):

            print(
                "| sI: %3d, sJ: %3d, sK: %3d |" % (self.sI, self.sJ, self.sK),
                end="",
            )
            return None

        # kps computes x([i]P), x([j]P), and x([k]P) for each i in I, j in J, and j in J.
        # I, J, and K are defined according to Examples 4.7 and 4.12 of https://eprint.iacr.org/2020/341
        def kps_s(self, P, A, i):

            # This functions computes all the independent data from the input of algorithm 2 of https://eprint.iacr.org/2020/341
            if self.sI == 0:

                # Case L[i] = 3 is super very special case (nothing to do)
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
                # Branch corresponds with L[i] = 5 and L[i] = 7
                # Recall sJ > 0, then sJ = 1
                assert self.sJ == 1
                P2 = self.curve.xdbl(P, A)

                self.J = [list(P)]

                I = [list(P2)]
                hI = [
                    list([(0 - P2[0]), P2[1]])
                ]  # we only need to negate x-coordinate of each point
                self.ptree_hI = self.poly_mul.product_tree(
                    hI, self.sI
                )  # product tree of hI

                if not self.SCALED_MULTIEVALUATION:
                    # Using remainder trees
                    self.ptree_hI = self.poly_redc.reciprocal_tree(
                        {'rpoly': [1], 'rdeg': 0, 'fpoly': [1], 'fdeg': 0, 'a': 1},
                        2 * self.sJ + 1,
                        self.ptree_hI,
                        self.sI,
                    )  # reciprocal of the root is used to compute sons' reciprocals

                else:
                    # Using scaled remainder trees
                    assert (2 * self.sJ - self.sI + 1) > self.sI
                    (
                        self.ptree_hI['reciprocal'],
                        self.ptree_hI['a'],
                    ) = self.poly_redc.reciprocal(
                        self.ptree_hI['poly'][::-1],
                        self.sI + 1,
                        2 * self.sJ - self.sI + 1,
                    )
                    self.ptree_hI['scaled'], self.ptree_hI['as'] = (
                        list(self.ptree_hI['reciprocal'][: self.sI]),
                        self.ptree_hI['a'],
                    )

                # Notice, 0 <= sK <= 1
                assert self.sK <= 1
                if self.sK == 1:
                    self.K = [list(P2)]
                else:
                    self.K = []

                self.SUB_SQUARED = [
                    self.J[0][0] - self.J[0][1]
                ]  # (Xj - Zj)
                self.SUB_SQUARED[0] **= 2  # (Xj - Zj)^2

                self.ADD_SQUARED = [
                    self.J[0][0] + self.J[0][1]
                ]  # (Xj + Zj)
                self.ADD_SQUARED[0] **= 2  # (Xj + Zj)^2

                self.XZJ4 = [
                    self.SUB_SQUARED[0] - self.ADD_SQUARED[0]
                ]  # -4*Xj*Zj
                return None

            # At this step, sI > 1
            assert self.sI > 1
            if self.sJ == 1:
                # This branch corresponds with L[i] = 11 and L[i] = 13
                Q = self.curve.xdbl(P, A)  # x([2]P)
                Q2 = self.curve.xdbl(Q, A)  # x([2]Q)

                self.J = [list(P)]

                I = [[0, 0]] * self.sI
                I[0] = list(Q)  # x(   Q)
                I[1] = self.curve.xadd(Q2, I[0], I[0])  # x([3]Q)
                for ii in range(2, self.sI, 1):
                    I[ii] = self.curve.xadd(I[ii - 1], Q2, I[ii - 2])  # x([2**i + 1]Q)

                hI = [
                    [(0 - iP[0]), iP[1]] for iP in I
                ]  # we only need to negate x-coordinate of each point
                self.ptree_hI = self.poly_mul.product_tree(
                    hI, self.sI
                )  # product tree of hI

                if not self.SCALED_MULTIEVALUATION:
                    # Using remainder trees
                    self.ptree_hI = self.poly_redc.reciprocal_tree(
                        {'rpoly': [1], 'rdeg': 0, 'fpoly': [1], 'fdeg': 0, 'a': 1},
                        2 * self.sJ + 1,
                        self.ptree_hI,
                        self.sI,
                    )  # reciprocal of the root is used to compute sons' reciprocals

                else:
                    # Using scaled remainder trees
                    assert (2 * self.sJ - self.sI + 1) <= self.sI
                    (
                        self.ptree_hI['scaled'],
                        self.ptree_hI['as'],
                    ) = self.poly_redc.reciprocal(
                        self.ptree_hI['poly'][::-1], self.sI + 1, self.sI
                    )
                    self.ptree_hI['reciprocal'], self.ptree_hI['a'] = (
                        list(
                            self.ptree_hI['scaled'][: (2 * self.sJ - self.sI + 1)]
                        ),
                        self.ptree_hI['as'],
                    )

                # Notice, 0 <= sK <= 1
                assert self.sK <= 1
                if self.sK == 1:
                    self.K = [list(Q)]
                else:
                    self.K = []

                self.SUB_SQUARED = [
                    self.J[0][0] - self.J[0][1]
                ]  # (Xj - Zj)
                self.SUB_SQUARED[0] **= 2  # (Xj - Zj)^2

                self.ADD_SQUARED = [
                    self.J[0][0] + self.J[0][1]
                ]  # (Xj + Zj)
                self.ADD_SQUARED[0] **= 2  # (Xj + Zj)^2

                self.XZJ4 = [
                    self.SUB_SQUARED[0] - self.ADD_SQUARED[0]
                ]  # -4*Xj*Zj
                return None

            # Now, we ensure sI >= sJ > 1
            assert self.sJ > 1

            # In other words, we can proceed by the general case

            # ------------------------------------------------
            # Computing [j]P for each j in {1, 3, ..., 2*sJ - 1}
            self.J = [[0, 0]] * self.sJ
            self.J[0] = list(P)  # x(   P)
            P2 = self.curve.xdbl(P, A)  # x([2]P)
            self.J[1] = self.curve.xadd(P2, self.J[0], self.J[0])  # x([3]P)
            for jj in range(2, self.sJ, 1):
                self.J[jj] = self.curve.xadd(
                    self.J[jj - 1], P2, self.J[jj - 2]
                )  # x([2*jj + 1]P)

            # -------------------------------------------------------
            # Computing [i]P for i in { (2*sJ) * (2i + 1) : 0 <= i < sI}
            bhalf_floor = self.sJ // 2
            bhalf_ceil = self.sJ - bhalf_floor
            P4 = self.curve.xdbl(P2, A)  # x([4]P)
            P2[0], P4[0] = cswap(
                P2[0], P4[0], self.sJ % 2
            )  # Constant-time swap
            P2[1], P4[1] = cswap(
                P2[1], P4[1], self.sJ % 2
            )  # x([4]P) <--- coditional swap ---> x([2]P)
            Q = self.curve.xadd(
                self.J[bhalf_ceil], self.J[bhalf_floor - 1], P2
            )  # Q := [2b]P
            P2[0], P4[0] = cswap(
                P2[0], P4[0], self.sJ % 2
            )  # Constant-time swap
            P2[1], P4[1] = cswap(
                P2[1], P4[1], self.sJ % 2
            )  # x([4]P) <--- coditional swap ---> x([2]P)

            I = [[0, 0]] * self.sI
            I[0] = list(Q)  # x(   Q)
            Q2 = self.curve.xdbl(Q, A)  # x([2]Q)
            I[1] = self.curve.xadd(Q2, I[0], I[0])  # x([3]Q)
            for ii in range(2, self.sI, 1):
                I[ii] = self.curve.xadd(I[ii - 1], Q2, I[ii - 2])  # x([2**i + 1]Q)

            # --------------------------------------------------------------
            # Computing [k]P for k in { 4*sJ*sI + 1, ..., l - 6, l - 4, l - 2}
            self.K = [[0, 0]] * self.sK

            if self.sK >= 1:
                self.K[0] = list(P2)  # x([l - 2]P) = x([-2]P) = x([2]P)

            if self.sK >= 2:
                self.K[1] = list(P4)  # x([l - 4]P) = x([-4]P) = x([4]P)

            for k in range(2, self.sK, 1):
                self.K[k] = self.curve.xadd(self.K[k - 1], P2, self.K[k - 2])

            # ------------------------------------------------------------------------------------------------------
            #                   ~~~~~~~~               ~~~~~~~~
            #                    |    |                 |    |
            # Computing h_I(W) = |    | (W - x([i]P)) = |    | (Zi * W - Xi) / Zi where x([i]P) = Xi/Zi
            #                    i in I                 i in I
            # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

            hI = [
                [(0 - iP[0]), iP[1]] for iP in I
            ]  # we only need to negate x-coordinate of each point
            self.ptree_hI = self.poly_mul.product_tree(
                hI, self.sI
            )  # product tree of hI

            if not self.SCALED_MULTIEVALUATION:
                # Using scaled remainder trees
                self.ptree_hI = self.poly_redc.reciprocal_tree(
                    {'rpoly': [1], 'rdeg': 0, 'fpoly': [1], 'fdeg': 0, 'a': 1},
                    2 * self.sJ + 1,
                    self.ptree_hI,
                    self.sI,
                )  # reciprocal of the root is used to compute sons' reciprocals

            else:
                # Using scaled remainder trees
                if self.sI < (2 * self.sJ - self.sI + 1):
                    (
                        self.ptree_hI['reciprocal'],
                        self.ptree_hI['a'],
                    ) = self.poly_redc.reciprocal(
                        self.ptree_hI['poly'][::-1],
                        self.sI + 1,
                        2 * self.sJ - self.sI + 1,
                    )
                    self.ptree_hI['scaled'], self.ptree_hI['as'] = (
                        list(self.ptree_hI['reciprocal'][: self.sI]),
                        self.ptree_hI['a'],
                    )

                else:
                    (
                        self.ptree_hI['scaled'],
                        self.ptree_hI['as'],
                    ) = self.poly_redc.reciprocal(
                        self.ptree_hI['poly'][::-1], self.sI + 1, self.sI
                    )
                    self.ptree_hI['reciprocal'], self.ptree_hI['a'] = (
                        list(
                            self.ptree_hI['scaled'][: (2 * self.sJ - self.sI + 1)]
                        ),
                        self.ptree_hI['as'],
                    )

            # Polynomial D_j of algorithm 2 from https://eprint.iacr.org/2020/341 is not required,
            # but we need some some squares and products determined by list J
            # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

            self.SUB_SQUARED = [0 for j in range(0, self.sJ, 1)]  #
            self.ADD_SQUARED = [0 for j in range(0, self.sJ, 1)]  #

            # List XZJ4 is required for degree-l isogeny evaluations...
            self.XZJ4 = [0 for j in range(0, self.sJ, 1)]  # 2*Xj*Zj
            for j in range(0, self.sJ, 1):

                self.SUB_SQUARED[j] = (
                    self.J[j][0] - self.J[j][1]
                )  # (Xj - Zj)
                self.SUB_SQUARED[j] **= 2  # (Xj - Zj)^2

                self.ADD_SQUARED[j] = (
                    self.J[j][0] + self.J[j][1]
                )  # (Xj + Zj)
                self.ADD_SQUARED[j] **= 2  # (Xj + Zj)^2

                self.XZJ4[j] = (
                    self.SUB_SQUARED[j] - self.ADD_SQUARED[j]
                )  # -4*Xj*Zj

            # Ensuring the cardinality of each ser coincide with the expected one
            assert len(I) == self.sI
            assert len(self.J) == self.sJ
            assert len(self.K) == self.sK

            return None

        # Next function perform algorithm 2 of https://eprint.iacr.org/2020/341 with input \alpha = 1 and \alpha = -1, and
        # then it computes the isogenous Montgomery curve coefficient
        def xisog_s(self, A, i):

            AA = (A[0] + A[0])  # 2A' + 4C
            AA -= A[1]  # 2A'
            AA += AA  # 4A'

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

                # However, each self.SUB_SQUARED[j] and self.ADD_SQUARED[j] should be multiplied by C
                tadd = (self.ADD_SQUARED[j] * A[1])
                tsub = (self.SUB_SQUARED[j] * A[1])

                # We require the double of tadd and tsub
                tadd2 = (tadd + tadd)
                tsub2 = (tsub + tsub)

                t1 = (self.XZJ4[j] * AA)  #       A *(-4*Xj*Zj)

                # Case alpha = 1
                linear = (
                    t1 - tadd2
                )  #       A *(-4*Xj*Zj)  - C * (2 * (Xj + Zj)^2)
                EJ_0[j] = [tsub, linear, tsub]

                # Case alpha = -1
                linear = (
                    tsub2 - t1
                )  #       C * (2 * (Xj - Zj)^2) - A *(-4*Xj*Zj)
                EJ_1[j] = [tadd, linear, tadd]

            # The faster way for multiplying is using a divide-and-conquer approach
            poly_EJ_0 = self.poly_mul.product_selfreciprocal_tree(EJ_0, self.sJ)[
                'poly'
            ]  # product tree of EJ_0 (we only require the root)
            poly_EJ_1 = self.poly_mul.product_selfreciprocal_tree(EJ_1, self.sJ)[
                'poly'
            ]  # product tree of EJ_1 (we only require the root)

            if not self.SCALED_MULTIEVALUATION:
                # Remainder tree computation
                remainders_EJ_0 = self.poly_redc.multieval_unscaled(
                    poly_EJ_0, 2 * self.sJ + 1, self.ptree_hI, self.sI
                )
                remainders_EJ_1 = self.poly_redc.multieval_unscaled(
                    poly_EJ_1, 2 * self.sJ + 1, self.ptree_hI, self.sI
                )

            else:
                # Approach using scaled remainder trees
                if self.ptree_hI != None:
                    poly_EJ_0 = self.poly_redc.poly_redc(
                        poly_EJ_0, 2 * self.sJ + 1, self.ptree_hI
                    )
                    fg_0 = self.poly_mul.poly_mul_middle(
                        self.ptree_hI['scaled'], self.sI, poly_EJ_0[::-1], self.sI
                    )
                    remainders_EJ_0 = self.poly_redc.multieval_scaled(
                        fg_0[::-1],
                        self.sI,
                        [1] + [0] * (self.sI - 1),
                        self.sI,
                        self.ptree_hI,
                        self.sI,
                    )

                    poly_EJ_1 = self.poly_redc.poly_redc(
                        poly_EJ_1, 2 * self.sJ + 1, self.ptree_hI
                    )
                    fg_1 = self.poly_mul.poly_mul_middle(
                        self.ptree_hI['scaled'], self.sI, poly_EJ_1[::-1], self.sI
                    )
                    remainders_EJ_1 = self.poly_redc.multieval_scaled(
                        fg_1[::-1],
                        self.sI,
                        [1] + [0] * (self.sI - 1),
                        self.sI,
                        self.ptree_hI,
                        self.sI,
                    )
                else:
                    remainders_EJ_0 = []
                    remainders_EJ_1 = []

            # Multipying all the remainders
            r0 = self.poly_mul.product(remainders_EJ_0, self.sI)
            r1 = self.poly_mul.product(remainders_EJ_1, self.sI)

            # ---------------------------------------------------------------------------------
            # Now, we proceed by computing the missing part which is determined by K
            # Notice, the denominators are the same and then they annulled between them
            # In other words, it is not required to compute the product of all Zk's with k In K

            # Case alpha = 1
            hK_0 = [
                [(self.K[k][1] - self.K[k][0])]
                for k in range(0, self.sK, 1)
            ]
            hK_0 = self.poly_mul.product(
                hK_0, self.sK
            )  # product of (Zk - Xk) for each k in K
            # Case alpha = -1
            hK_1 = [
                [(self.K[k][1] + self.K[k][0])]
                for k in range(0, self.sK, 1)
            ]
            hK_1 = self.poly_mul.product(
                hK_1, self.sK
            )  # product of (Zk + Xk) for each k in K

            # --------------------------------------------------------------
            # Now, we have all the ingredients for computing the image curve
            A24m = (A[0] - A[1])  # A' - 2C

            A24 = (A[0] ** self.L[i])  # (A' + 2C)^l
            A24m **= self.L[i]  # (A' - 2C)^l

            t24m = (
                r1 * hK_1
            )  # output of algorithm 2 with alpha =-1 and without the demoninator
            t24m **= 2  # raised at 2
            t24m **= 2  # raised at 4
            t24m **= 2  # raised at 8

            t24 = (
                r0 * hK_0
            )  # output of algorithm 2 with alpha = 1 and without the demoninator
            t24 **= 2  # raised at 2
            t24 **= 2  # raised at 4
            t24 **= 2  # raised at 8

            A24 *= t24m
            A24m *= t24

            # Now, we have d = (A24m / A24) where the image Montgomery cuve coefficient is
            #      B'   2*(1 + d)   2*(A24 + A24m)
            # B = ---- = --------- = --------------
            #      C      (1 - d)     (A24 - A24m)
            # However, we required B' + 2C = 4*A24 and 4C = 4 * (A24 - A24m)

            t24m = (A24 - A24m)  #   (A24 - A24m)
            t24m += t24m  # 2*(A24 - A24m)
            t24m += t24m  # 4*(A24 - A24m)

            t24 = (A24 + A24)  # 2 * A24
            t24 += t24  # 4 * A24

            # return [t24, t24m], ptree_hI, XZJ4
            return [t24, t24m]

        def xeval_s(self, P, A):

            AA = (A[0] + A[0])  # 2A' + 4C
            AA -= A[1]  # 2A'
            AA += AA  # 4A'

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

            XZ_add = (P[0] + P[1])  # X + Z
            XZ_sub = (P[0] - P[1])  # X - Z

            AXZ2 = (P[0] * P[1])  # X * Z
            t1 = (P[0] ** 2)  # X^2
            t2 = (P[1] ** 2)  # Z^2

            CX2Z2 = (t1 + t2)  #      X^2 + Z^2
            CX2Z2 *= A[1]  # C * (X^2 + Z^2)

            AXZ2 += AXZ2  #       2 * (X * Z)
            CXZ2 = (AXZ2 * A[1])  # C  * [2 * (X * Z)]
            AXZ2 *= AA  # A' * [2 * (X * Z)]

            for j in range(0, self.sJ, 1):

                XZj_add = (self.J[j][0] + self.J[j][1])  # Xj + Zj
                XZj_sub = (self.J[j][0] - self.J[j][1])  # Xj - Zj

                t1 = (XZ_sub * XZj_add)  # (X - Z) * (Xj + Zj)
                t2 = (XZ_add * XZj_sub)  # (X + Z) * (Xj - Zj)

                # Computing the quadratic coefficient
                quadratic = (t1 - t2)  #   2 * [(X*Zj) - (Z*Xj)]
                quadratic **= 2 # ( 2 * [(X*Zj) - (Z*Xj)] )^2
                quadratic *= A[1] # C * ( 2 * [(X*Zj) - (Z*Xj)] )^2

                # Computing the constant coefficient
                constant = t1 # safe because we reassign t1 in the next section
                constant += t2  #   2 * [(X*Xj) - (Z*Zj)]
                constant **= 2  # ( 2 * [(X*Xj) - (Z*Zj)] )^2
                constant *= A[1] # C * ( 2 * [(X*Xj) - (Z*Zj)] )^2

                # Computing the linear coefficient
                # ----------------------------------------------------------------------------------------------------------
                # C * [ (-2*Xj*Zj)*(alpha^2 + 1) + (-2*alpha)*(Xj^2 + Zj^2)] + [A' * (-2*Xj*Zj) * (2*X*Z)] where alpha = X/Z
                t1 = (self.J[j][0] + self.J[j][1])  #     (Xj + Zj)
                t1 **= 2  #     (Xj + Zj)^2
                t1 += t1  # 2 * (Xj + Zj)^2
                t1 += (self.XZJ4[j]
                       ) # 2 * (Xj + Zj)^2 - (4*Xj*Zj) := 2 * (Xj^2 + Zj^2)
                t1 *= (CXZ2
                       ) # [2 * (Xj^2 + Zj^2)] * (2 * [ C * (X * Z)])

                t2 = (
                    CX2Z2 * self.XZJ4[j]
                )  # [C * (X^2 + Z^2)] * (-4 * Xj * Zj)
                t1 = (
                    t2 - t1
                )  # [C * (X^2 + Z^2)] * (-4 * Xj * Zj) - [2 * (Xj^2 + Zj^2)] * (2 * [ C * (X * Z)])

                t2 = (
                    AXZ2 * self.XZJ4[j]
                )  # (2 * [A' * (X * Z)]) * (-4 * Xj * Zj)
                linear = (
                    t1 + t2
                )  # This is our desired equation but multiplied by 2
                linear += linear
                # This is our desired equation but multiplied by 4
                # ----------------------------------------------------------------------------------------------------------

                # Case alpha = X / Z
                EJ_0[j] = [constant, linear, quadratic]

            # The faster way for multiplying is using a divide-and-conquer approach
            poly_EJ_0 = self.poly_mul.product_tree(EJ_0, self.sJ)[
                'poly'
            ]  # product tree of EJ_0 (we only require the root)
            poly_EJ_1 = list(
                poly_EJ_0[::-1]
            )  # product tree of EJ_1(x) = x^{2b + 1} EJ_0(1/X)

            if not self.SCALED_MULTIEVALUATION:
                # Remainder tree computation
                remainders_EJ_0 = self.poly_redc.multieval_unscaled(
                    poly_EJ_0, 2 * self.sJ + 1, self.ptree_hI, self.sI
                )
                remainders_EJ_1 = self.poly_redc.multieval_unscaled(
                    poly_EJ_1, 2 * self.sJ + 1, self.ptree_hI, self.sI
                )

            else:
                # Approach using scaled remainder trees
                if self.ptree_hI != None:
                    poly_EJ_0 = self.poly_redc.poly_redc(
                        poly_EJ_0, 2 * self.sJ + 1, self.ptree_hI
                    )
                    fg_0 = self.poly_mul.poly_mul_middle(
                        self.ptree_hI['scaled'], self.sI, poly_EJ_0[::-1], self.sI
                    )
                    remainders_EJ_0 = self.poly_redc.multieval_scaled(
                        fg_0[::-1],
                        self.sI,
                        [1] + [0] * (self.sI - 1),
                        self.sI,
                        self.ptree_hI,
                        self.sI,
                    )

                    poly_EJ_1 = self.poly_redc.poly_redc(
                        poly_EJ_1, 2 * self.sJ + 1, self.ptree_hI
                    )
                    fg_1 = self.poly_mul.poly_mul_middle(
                        self.ptree_hI['scaled'], self.sI, poly_EJ_1[::-1], self.sI
                    )
                    remainders_EJ_1 = self.poly_redc.multieval_scaled(
                        fg_1[::-1],
                        self.sI,
                        [1] + [0] * (self.sI - 1),
                        self.sI,
                        self.ptree_hI,
                        self.sI,
                    )
                else:
                    remainders_EJ_0 = []
                    remainders_EJ_1 = []

            # Multipying all the remainders
            r0 = self.poly_mul.product(remainders_EJ_0, self.sI)
            r1 = self.poly_mul.product(remainders_EJ_1, self.sI)

            # ---------------------------------------------------------------------------------
            # Now, we proceed by computing the missing part which is determined by K
            # Notice, the denominators are the same and then they annulled between them
            # In other words, it is not required to compute the product of all Zk's with k In K

            hK_0 = [[0]] * self.sK
            hK_1 = [[0]] * self.sK
            for k in range(0, self.sK, 1):

                t1 = (self.K[k][0] + self.K[k][1])  # Xk + Zk
                t2 = (self.K[k][0] - self.K[k][1])  # Xk - Zk
                t1 *= XZ_sub  # (X - Z) * (Xk + Zk)
                t2 *= XZ_add  # (X + Z) * (Xk - Zk)

                # Case alpha = X/Z
                hK_0[k] = [(t1 - t2)]  # 2 * [(X*Zk) - (Z*Xk)]

                # Case 1/alpha = Z/X
                t1 += t2 # hK_1[k] = [t1 + t2]
                hK_1[k] = [t1]  # 2 * [(X*Xk) - (Z*Zk)]

            hK_0 = self.poly_mul.product(
                hK_0, self.sK
            )  # product of (XZk - ZXk) for each k in K
            hK_1 = self.poly_mul.product(
                hK_1, self.sK
            )  # product of (XXk - ZZk) for each k in K

            # ---------------------------------------------------------------------------------
            # Now, unifying all the computations
            XX = (
                r1 * hK_1
            )  # output of algorithm 2 with 1/alpha = Z/X and without the demoninator
            XX **= 2
            XX *= P[0]

            ZZ = (
                r0 * hK_0
            )  # output of algorithm 2 with alpha = X/Z and without the demoninator
            ZZ **= 2
            ZZ *= P[1]

            return [XX, ZZ]

        # Kernel computation for xeval_2 and xisog_2
        def kps_2(self, P):

            self.K = [None, None]
            self.K[0] = (P[0] + P[1])
            self.K[1] = (P[0] - P[1])
            return None

        # Degree-2 isogeny construction
        def xisog_2(self, P):

            A24 = (P[0] ** 2)
            C24 = (P[1] ** 2)
            A24 = (C24 - A24)
            return [A24, C24]

        # Degree-2 isogeny evluation
        def xeval_2(self, Q):

            t2 = (Q[0] + Q[1])
            t3 = (Q[0] - Q[1])
            t0 = (self.K[0] * t3)
            t1 = (self.K[1] * t2)
            t2 = (t0 + t1)
            t3 = (t0 - t1)
            XQ = (Q[0] * t2)
            ZQ = (Q[1] * t3)
            return [XQ, ZQ]

        # Kernel computation for xeval_4 and xisog_4
        def kps_4(self, P):
            self.K = [ None, None, None, None]

            self.K[1] = (P[0] - P[1])
            self.K[2] = (P[0] + P[1])
            self.K[0] = (P[1] ** 2)
            self.K[3] = (self.K[0] + self.K[0])
            self.K[0] = (self.K[3] + self.K[3])
            return None

        # Degree-4 isogeny construction
        def xisog_4(self, P):

            C24 = (self.K[3] ** 2)
            A24 = (P[0] ** 2)
            A24 += A24
            A24 **= 2
            return [A24, C24]


        # Degree-4 isogeny evaluation
        def xeval_4(self, Q):

            t0 = (Q[0] + Q[1])
            t1 = (Q[0] - Q[1])
            XQ = (t0 * self.K[1])
            ZQ = (t1 * self.K[2])
            t0 *= t1
            t0 *= self.K[0]
            t1 = (XQ + ZQ)
            ZQ = (XQ - ZQ)
            t1 **= 2
            ZQ **= 2
            XQ = (t0 + t1)
            t0 = (ZQ - t0)
            XQ *= t1
            ZQ *= t0
            return [XQ, ZQ]

        # Kernel computation for xeval_3 and xisog_3
        def kps_3(self, P):

            self.K = [None, None]
            self.K[0] = (P[0] - P[1])
            self.K[1] = (P[0] + P[1])
            return None

        # Degree-3 isogeny construction
        def xisog_3(self, P):

            # This function returns (A' + 2C, A' - 2C), which is required for tripling a point
            t0 = (self.K[0] ** 2)
            t1 = (self.K[1] ** 2)
            t2 = (t0 + t1)
            #assert((self.K[0] + self.K[1]) == (P[0] + P[0]))
            t3 = (self.K[0] + self.K[1])
            #t3 = (P[0] + P[0])
            t3 **= 2
            t3 -= t2
            t2 = (t1 + t3)
            t3 += t0
            t4 = (t3 + t0)
            t4 += t4
            t4 = (t1 + t4)
            A24m = (t2 * t4)
            t4 = (t1 + t2)
            t4 += t4
            t4 = (t0 + t4)
            A24p = (t3 * t4)
            return [A24p, A24m]

        # Degree-3 isogeny evaluation
        def xeval_3(self, Q):

            t0 = (Q[0] + Q[1])
            t1 = (Q[0] - Q[1])
            t0 *= self.K[0]
            t1 *= self.K[1]
            t2 = (t1 + t0)
            t0 = (t1 - t0)
            t2 **= 2
            t0 **= 2
            XQ = (Q[0] * t2)
            ZQ = (Q[1] * t0)
            return [XQ, ZQ]

        def kps(self, P, A, i):

            if self.L[i] == 4:
                return self.kps_4(P)

            if self.L[i] <= self.HYBRID_BOUND:
                return self.kps_t(P, A, i)
            else:
                return self.kps_s(P, A, i)

        def xisog(self, A, i):

            if self.L[i] == 4:
                return self.xisog_4(A)

            if self.L[i] <= self.HYBRID_BOUND:
                return self.xisog_t(A, i)
            else:
                return self.xisog_s(A, i)

        def xeval(self, P, v):

            if type(v) == int:
                if self.L[v] != 4:
                    return self.xeval_t(P, v)
                else:
                    return self.xeval_4(P)
            else:
                return self.xeval_s(P, v)

        def velusqrt_cost(self):

            if self.curve.name in parameters['csidh'].keys():
                # First, we look for a full torsion point
                A = [self.field(2), self.field(4)]
                T_p, T_m = self.curve.generators(A)

                for i in range(0, self.curve.n, 1):

                    if self.tuned:
                        self.set_parameters_velu(self.sJ_list[i], self.sI_list[i], i)
                    else:
                        # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                        # These paramters are required in kps, xisog, and xeval
                        if self.L[i] == 3:
                            b = 0
                            c = 0
                        else:
                            b = int(floor(sqrt(self.L[i] - 1) / 2.0))
                            c = int(floor((self.L[i] - 1.0) / (4.0 * b)))

                        self.set_parameters_velu(b, c, i)

                    # Getting an orderl-l point
                    Tp = list(T_p)
                    for j in range(i + 1, self.curve.n, 1):
                        Tp = self.curve.xmul(Tp, A, j)

                    # Cost of xisog() and kps()
                    self.field.init_runtime()
                    self.kps(Tp, A, i)
                    t = [self.field.fpmul, self.field.fpsqr, self.field.fpadd]
                    self.c_xisog[i] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

                    self.field.init_runtime()
                    B = self.xisog(A, i)
                    t = [self.field.fpmul, self.field.fpsqr, self.field.fpadd]
                    self.c_xisog[i] += numpy.array(
                        [t[0] * 1.0, t[1] * 1.0, t[2] * 1.0]
                    )

                    # xeval: kernel point determined by the next isogeny evaluation
                    self.field.init_runtime()
                    if self.L[i] <= self.HYBRID_BOUND:
                        T_p = self.xeval(T_p, i)
                    else:
                        T_p = self.xeval(T_p, A)

                    # Cost of xeval
                    self.field.init_runtime()
                    if self.L[i] <= self.HYBRID_BOUND:
                        T_m = self.xeval(T_m, i)
                    else:
                        T_m = self.xeval(T_m, A)

                    t = [self.field.fpmul, self.field.fpsqr, self.field.fpadd]
                    self.c_xeval[i] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

                    # Updating the new next curve
                    A = list(B)

                self.field.init_runtime()
                return None

            elif self.curve.name in parameters['bsidh'].keys():

                random = SystemRandom()

                # Reading public generators points
                f = open(resource_filename('sibc', 'data/gen/' + self.curve.model + '/bsidh/' + self.curve.name))

                # x(PA), x(QA) and x(PA - QA)
                PQA = f.readline()
                PQA = [int(x, 16) for x in PQA.split()]
                PA = [self.field(PQA[0:2]), self.field(1)]
                QA = [self.field(PQA[2:4]), self.field(1)]
                PQA= [self.field(PQA[4:6]), self.field(1)]

                # x(PB), x(QB) and x(PB - QB)
                PQB = f.readline()
                PQB = [int(x, 16) for x in PQB.split()]
                PB = [self.field(PQB[0:2]), self.field(1)]
                QB = [self.field(PQB[2:4]), self.field(1)]
                PQB= [self.field(PQB[4:6]), self.field(1)]

                f.close()

                # Three point ladder: case (p + 1)
                A = [ self.field(8), self.field(4)]
                S = list(PA)
                T = list(QA)
                ST = list(PQA)

                assert self.curve.isinfinity(S) == False
                assert self.curve.isinfinity(T) == False

                for i in range(0, self.curve.np, 1):
                    for idx in range(0, self.curve.Ep[i] - 1, 1):
                        S = self.curve.xmul(S, A, i)
                        T = self.curve.xmul(T, A, i)
                        ST = self.curve.xmul(ST, A, i)

                k = random.randint(0, self.field.basefield.p + 1)
                R = self.curve.Ladder3pt(k, S, T, ST, self.field(6))
                T_p = list(R)
                T_m = list(S)
                for idx in range(0, self.curve.np, 1):

                    # -------------------------------------------------------------
                    # Random kernel point
                    Tp = list(T_p)
                    for i in range(idx + 1, self.curve.np, 1):
                        Tp = self.curve.xmul(Tp, A, i)

                    if self.tuned:
                        self.set_parameters_velu(self.sJ_list[idx], self.sI_list[idx], idx)
                    else:
                        # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                        # These paramters are required in kps, xisog, and xeval
                        if self.L[idx] <= 4:
                            b = 0
                            c = 0
                        else:
                            b = int(floor(sqrt(self.L[idx] - 1) / 2.0))
                            c = int(floor((self.L[idx] - 1.0) / (4.0 * b)))

                        self.set_parameters_velu(b, c, idx)

                    # -------------------------------------------------------------
                    # kps procedure
                    self.field.basefield.init_runtime()
                    self.field.init_runtime()
                    self.kps(Tp, A, idx)
                    t = [self.field.basefield.fpmul, self.field.basefield.fpsqr, self.field.basefield.fpadd]
                    self.c_xisog[idx] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

                    # -------------------------------------------------------------
                    # xisog
                    self.field.basefield.init_runtime()
                    self.field.init_runtime()
                    Tp[0], A[0] = cswap(Tp[0], A[0], self.L[idx] == 4)
                    Tp[1], A[1] = cswap(Tp[1], A[1], self.L[idx] == 4)
                    B = self.xisog(A, idx)
                    Tp[0], A[0] = cswap(Tp[0], A[0], self.L[idx] == 4)
                    Tp[1], A[1] = cswap(Tp[1], A[1], self.L[idx] == 4)
                    t = [self.field.basefield.fpmul, self.field.basefield.fpsqr, self.field.basefield.fpadd]
                    self.c_xisog[idx] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

                    # -------------------------------------------------------------
                    # xeval: kernel point determined by the next isogeny evaluation
                    self.field.basefield.init_runtime()
                    self.field.init_runtime()
                    if self.L[idx] <= self.HYBRID_BOUND or self.L[idx] == 4:
                        T_p = self.xeval(T_p, idx)
                    else:
                        T_p = self.xeval(T_p, A)

                    # xeval bench
                    self.field.basefield.init_runtime()
                    self.field.init_runtime()
                    if self.L[idx] <= self.HYBRID_BOUND or self.L[idx] == 4:
                        T_m = self.xeval(T_m, idx)
                    else:
                        T_m = self.xeval(T_m, A)

                    t = [self.field.basefield.fpmul, self.field.basefield.fpsqr, self.field.basefield.fpadd]
                    self.c_xeval[idx] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

                    # assert(validate(B))
                    A = list(B)

                # Three point ladder: case (p - 1)
                A = [ self.field(8), self.field(4)]
                S = list(PB)
                T = list(QB)
                ST = list(PQB)

                assert self.curve.isinfinity(S) == False
                assert self.curve.isinfinity(T) == False

                for i in range(self.curve.np, self.curve.np + self.curve.nm, 1):
                    for idx in range(0, self.curve.Em[i - self.curve.np] - 1, 1):
                        S = self.curve.xmul(S, A, i)
                        T = self.curve.xmul(T, A, i)
                        ST= self.curve.xmul(ST, A, i)

                k = random.randint(0, self.field.basefield.p - 1)
                R = self.curve.Ladder3pt(k, S, T, ST, self.field(6))
                T_p = list(R)
                T_m = list(S)

                for idx in range(self.curve.np, self.curve.np + self.curve.nm, 1):

                    # -------------------------------------------------------------
                    # Random kernel point
                    Tp = list(T_p)
                    for i in range(idx + 1, self.curve.np + self.curve.nm, 1):
                        Tp = self.curve.xmul(Tp, A, i)

                    if self.tuned:
                        self.set_parameters_velu(self.sJ_list[idx], self.sI_list[idx], idx)
                    else:
                        # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                        # These paramters are required in kps, xisog, and xeval
                        if self.L[idx] == 3:
                            b = 0
                            c = 0
                        else:
                            b = int(floor(sqrt(self.L[idx] - 1) / 2.0))
                            c = int(floor((self.L[idx] - 1.0) / (4.0 * b)))

                        self.set_parameters_velu(b, c, idx)

                    # -------------------------------------------------------------
                    # kps procedure
                    self.field.basefield.init_runtime()
                    self.field.init_runtime()
                    self.kps(Tp, A, idx)
                    t = [self.field.basefield.fpmul, self.field.basefield.fpsqr, self.field.basefield.fpadd]
                    self.c_xisog[idx] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

                    # -------------------------------------------------------------
                    # xisog
                    self.field.basefield.init_runtime()
                    self.field.init_runtime()
                    B = self.xisog(A, idx)
                    t = [self.field.basefield.fpmul, self.field.basefield.fpsqr, self.field.basefield.fpadd]
                    self.c_xisog[idx] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

                    # -------------------------------------------------------------
                    # xeval: kernel point determined by the next isogeny evaluation
                    self.field.basefield.init_runtime()
                    self.field.init_runtime()
                    if self.L[idx] <= self.HYBRID_BOUND:
                        T_p = self.xeval(T_p, idx)
                    else:
                        T_p = self.xeval(T_p, A)

                    # xeval bench
                    self.field.basefield.init_runtime()
                    self.field.init_runtime()
                    if self.L[idx] <= self.HYBRID_BOUND:
                        T_m = self.xeval(T_m, idx)
                    else:
                        T_m = self.xeval(T_m, A)

                    t = [self.field.basefield.fpmul, self.field.basefield.fpsqr, self.field.basefield.fpadd]
                    self.c_xeval[idx] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

                    # assert(validate(B))
                    A = list(B)
            else:
                assert False, "only CSIDH and B-SIDH are currently implemented"


    Formulae.__name__ = NAME
    Formulae.name = name
    return Formulae
