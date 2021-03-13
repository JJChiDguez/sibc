import numpy
from sympy import floor, sqrt, sign
from random import SystemRandom

from sidh.bsidh.poly_redc import Poly_redc
from sidh.bsidh.poly_mul import Poly_mul
from sidh.math import hamming_weight, isequal, bitlength
from sidh.common import attrdict
from sidh.bsidh.montgomery import MontgomeryCurve
from pkg_resources import resource_filename


def Tvelu(curve):

    global_L = curve.global_L
    prime = curve.prime
    fp = curve.fp
    p = curve.p
    np = curve.np
    nm = curve.nm
    Em = curve.Em
    Ep = curve.Ep
    random = SystemRandom()
    poly_mul = Poly_mul(curve)
    poly_redc = Poly_redc(poly_mul)


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

    K = None



    def yDBL(P, A):
        '''
        -------------------------------------------------------------------------
        yDBL()
        input : a projective Twisted Edwards y-coordinate point y(P) := YP/WP,
                and the projective Montgomery constants A24:= A + 2C and C24:=4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Twisted Edwards y-coordinate point y([2]P)
        -------------------------------------------------------------------------
        '''
        t_0 = fp.fp2_sqr(P[0])
        t_1 = fp.fp2_sqr(P[1])
        Z = fp.fp2_mul(A[1], t_0)
        X = fp.fp2_mul(Z, t_1)
        t_1 = fp.fp2_sub(t_1, t_0)
        t_0 = fp.fp2_mul(A[0], t_1)
        Z = fp.fp2_add(Z, t_0)
        Z = fp.fp2_mul(Z, t_1)
        return [fp.fp2_sub(X, Z), fp.fp2_add(X, Z)]

    def yADD(P, Q, PQ):
        '''
        -------------------------------------------------------------------------
        yADD()
        input : the projective Twisted Edwards y-coordinate points y(P) := YP/WP,
                y(Q) := YQ/WQ, and y(P-Q) := YPQ/QPQ
        output: the projective Twisted Edwards y-coordinate point y(P+Q)
        -------------------------------------------------------------------------
        '''
        a = fp.fp2_mul(P[1], Q[0])
        b = fp.fp2_mul(P[0], Q[1])
        c = fp.fp2_add(a, b)
        d = fp.fp2_sub(a, b)
        c = fp.fp2_sqr(c)
        d = fp.fp2_sqr(d)

        xD = fp.fp2_add(PQ[1], PQ[0])
        zD = fp.fp2_sub(PQ[1], PQ[0])
        X = fp.fp2_mul(zD, c)
        Z = fp.fp2_mul(xD, d)
        return [fp.fp2_sub(X, Z), fp.fp2_add(X, Z)]

    def KPs_t(P, A, i):
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
        assert curve.isinfinity(P) == False
        nonlocal K
        d = (global_L[i] - 1) // 2

        K = [[[0, 0], [0, 0]] for j in range(d + 1)]
        K[0] = list([fp.fp2_sub(P[0], P[1]), fp.fp2_add(P[0], P[1])])  # 2a
        if global_L[i] == 3:
            K[1] = list([list(K[0]), list(K[1])])
        else:
            K[1] = yDBL(K[0], A)  # 4M + 2S + 4a

        for j in range(2, d, 1):
            K[j] = yADD(K[j - 1], K[0], K[j - 2])  # 4M + 2S + 6a

        return K  # 2(l - 3)M + (l - 3)S + 3(l - 3)a

    def xISOG_t(A, i):
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
        nonlocal K
        l = global_L[i]  # l
        bits_of_l = bitlength(l)  # Number of bits of L[i]
        d = (l - 1) // 2  # Here, l = 2d + 1

        By = K[0][0]
        Bz = K[0][1]
        for j in range(1, d, 1):

            By = fp.fp2_mul(By, K[j][0])
            Bz = fp.fp2_mul(Bz, K[j][1])

        bits_of_l -= 1
        constant_d_edwards = fp.fp2_sub(A[0], A[1])  # 1a

        tmp_a = A[0]
        tmp_d = constant_d_edwards
        # left-to-right method for computing a^l and d^l
        for j in range(1, bits_of_l + 1):

            tmp_a = fp.fp2_sqr(tmp_a)
            tmp_d = fp.fp2_sqr(tmp_d)
            if ((l >> (bits_of_l - j)) & 1) != 0:

                tmp_a = fp.fp2_mul(tmp_a, A[0])
                tmp_d = fp.fp2_mul(tmp_d, constant_d_edwards)

        for j in range(3):

            By = fp.fp2_sqr(By)
            Bz = fp.fp2_sqr(Bz)

        C0 = fp.fp2_mul(tmp_a, Bz)
        C1 = fp.fp2_mul(tmp_d, By)
        C1 = fp.fp2_sub(C0, C1)

        return [C0, C1]  # (l - 1 + 2*HW(l) - 2)M + 2(|l|_2 + 1)S + 2a

    def xEVAL_t(P, i):
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
        nonlocal K
        d = (global_L[i] - 1) // 2  # Here, l = 2d + 1

        Q0 = fp.fp2_add(P[0], P[1])
        Q1 = fp.fp2_sub(P[0], P[1])
        R0, R1 = curve.CrissCross(K[0][1], K[0][0], Q0, Q1)
        for j in range(1, d, 1):

            T0, T1 = curve.CrissCross(K[j][1], K[j][0], Q0, Q1)
            R0 = fp.fp2_mul(T0, R0)
            R1 = fp.fp2_mul(T1, R1)

        R0 = fp.fp2_sqr(R0)
        R1 = fp.fp2_sqr(R1)
        X = fp.fp2_mul(P[0], R0)
        Z = fp.fp2_mul(P[1], R1)

        return [X, Z]  # 2(l - 1)M + 2S + (l + 1)a

    def xISOG_4(P):
        """
        Degree-4 isogeny construction
        """
        nonlocal K
        K = [[0, 0], [0, 0], [0, 0]]

        K[1] = fp.fp2_sub(P[0], P[1])
        K[2] = fp.fp2_add(P[0], P[1])
        K[0] = fp.fp2_sqr(P[1])
        K[0] = fp.fp2_add(K[0], K[0])

        C24 = fp.fp2_sqr(K[0])
        K[0] = fp.fp2_add(K[0], K[0])
        A24 = fp.fp2_sqr(P[0])
        A24 = fp.fp2_add(A24, A24)
        A24 = fp.fp2_sqr(A24)
        return [A24, C24]

    def xEVAL_4(Q):
        """
        Degree-4 isogeny evaluation
        """
        t0 = fp.fp2_add(Q[0], Q[1])
        t1 = fp.fp2_sub(Q[0], Q[1])
        XQ = fp.fp2_mul(t0, K[1])
        ZQ = fp.fp2_mul(t1, K[2])
        t0 = fp.fp2_mul(t0, t1)
        t0 = fp.fp2_mul(t0, K[0])
        t1 = fp.fp2_add(XQ, ZQ)
        ZQ = fp.fp2_sub(XQ, ZQ)
        t1 = fp.fp2_sqr(t1)
        ZQ = fp.fp2_sqr(ZQ)
        XQ = fp.fp2_add(t0, t1)
        t0 = fp.fp2_sub(ZQ, t0)
        XQ = fp.fp2_mul(XQ, t1)
        ZQ = fp.fp2_mul(ZQ, t0)

        return [XQ, ZQ]

    def KPs(P, A, i):
        if global_L[i] != 4:
            return KPs_t(P, A, i)

    def xISOG(A, i):
        if global_L[i] != 4:
            return xISOG_t(A, i)
        else:
            # A should corresponds with an order-4 point
            return xISOG_4(A)

    def xEVAL(P, i):
        if global_L[i] != 4:
            return xEVAL_t(P, i)
        else:
            return xEVAL_4(P)

    def cISOG_and_cEVAL():
        """
        Get cost of the isogeny constructions and evaluations
        """
        nonlocal C_xISOG
        nonlocal C_xEVAL

        # E[p + 1]
        # First, we look for a full torsion point
        A = [[0x8, 0x0], [0x4, 0x0]]

        # Reading public generators points
        f = open(resource_filename(__name__, '../data/gen/' + prime))
        # x(PA), x(QA) and x(PA - QA)
        PQA = f.readline()
        PQA = [int(x, 16) for x in PQA.split()]
        PA = [list(PQA[0:2]), [0x1, 0x0]]
        QA = [list(PQA[2:4]), [0x1, 0x0]]
        PQA = [list(PQA[4:6]), [0x1, 0x0]]

        # x(PB), x(QB) and x(PB - QB)
        PQB = f.readline()
        PQB = [int(x, 16) for x in PQB.split()]
        PB = [list(PQB[0:2]), [0x1, 0x0]]
        QB = [list(PQB[2:4]), [0x1, 0x0]]
        PQB = [list(PQB[4:6]), [0x1, 0x0]]

        f.close()

        for i in range(0, Ep[0] - 1, 1):

            PA = curve.xMUL(PA, A, 0)
            QA = curve.xMUL(QA, A, 0)
            PQA = curve.xMUL(PQA, A, 0)

        # Random kernels for counting the
        T_p = curve.Ladder3pt(random.randint(0, p - 1), PA, QA, PQA, A)
        T_m = curve.Ladder3pt(random.randint(0, p - 1), PB, QB, PQB, A)

        for i in range(0, np, 1):
            for j in range(0, Ep[i] - 1, 1):
                T_p =curve. xMUL(T_p, A, i)

        for i in range(np, np + nm, 1):
            for j in range(0, Em[i - np] - 1, 1):
                T_m =curve. xMUL(T_m, A, i)

        for i in range(0, np, 1):

            # Getting an orderl-l point
            Tp = list(T_p)
            for j in range(i + 1, np, 1):
                Tp =curve. xMUL(Tp, A, j)

            # Cost of xISOG() and KPs()
            fp.fp.set_zero_ops()
            KPs(Tp, A, i)
            t = fp.fp.get_ops()
            C_xISOG[i] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            fp.fp.set_zero_ops()
            Tp[0], A[0] = fp.fp2_cswap(Tp[0], A[0], global_L[i] == 4)
            Tp[1], A[1] = fp.fp2_cswap(Tp[1], A[1], global_L[i] == 4)
            B = xISOG(A, i)
            Tp[0], A[0] = fp.fp2_cswap(Tp[0], A[0], global_L[i] == 4)
            Tp[1], A[1] = fp.fp2_cswap(Tp[1], A[1], global_L[i] == 4)
            t = fp.fp.get_ops()
            C_xISOG[i] += numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            # xEVAL: kernel point determined by the next isogeny evaluation
            fp.fp.set_zero_ops()
            T_p = xEVAL(T_p, i)

            # Cost of xEVAL
            fp.fp.set_zero_ops()
            T_m = xEVAL(T_m, i)
            t = fp.fp.get_ops()
            C_xEVAL[i] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            # Updating the new next curve
            A = list(B)

        # E[p - 1]
        # First, we look for a full torsion point
        A = [[0x8, 0x0], [0x4, 0x0]]
        T_m = curve.Ladder3pt(random.randint(0, p - 1), PA, QA, PQA, A)
        T_p = curve.Ladder3pt(random.randint(0, p - 1), PB, QB, PQB, A)

        for i in range(np, np + nm, 1):

            # Getting an orderl-l point
            Tp = list(T_p)
            for j in range(i + 1, np + nm, 1):
                Tp =curve. xMUL(Tp, A, j)

            # Cost of xISOG() and KPs()
            fp.fp.set_zero_ops()
            KPs(Tp, A, i)
            t = fp.fp.get_ops()
            C_xISOG[i] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            fp.fp.set_zero_ops()
            B = xISOG(A, i)
            t = fp.fp.get_ops()
            C_xISOG[i] += numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            # xEVAL: kernel point determined by the next isogeny evaluation
            fp.fp.set_zero_ops()
            T_p = xEVAL(T_p, i)

            # Cost of xEVAL
            fp.fp.set_zero_ops()
            T_m = xEVAL(T_m, i)
            t = fp.fp.get_ops()
            C_xEVAL[i] = numpy.array([t[0] * 1.0, t[1] * 1.0, t[2] * 1.0])

            # Updating the new next curve
            A = list(B)
        fp.fp.set_zero_ops()
        return None


    # Now, we proceed to store all the correct costs
    cISOG_and_cEVAL()
    return attrdict(name='tvelu', **locals())
