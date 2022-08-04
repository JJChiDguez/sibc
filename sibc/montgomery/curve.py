from random import SystemRandom
from pkg_resources import resource_filename

import numpy

# Loading the library corresponding to the arithmetic in GF(p)
from sibc.primefield import PrimeField
# Loading the library corresponding to the arithmetic in GF(p²)
from sibc.quadraticfield import QuadraticField

from functools import reduce
from copy import copy

from sibc.constants import bitlength, parameters
from sibc.common import attrdict
from sibc.math import cswap
from math import sqrt

def filename_to_list_of_lists_of_ints(path):
    res = []
    try:
        with open(path) as fh:
            for line in fh:
                res.append(list(map(int, line.split())))
    except:
        res = []
    return res

# MontgomeryCurve class determines the family of supersingular elliptic curves over GF(p)
def MontgomeryCurve(prime):

    name = prime
    model= 'montgomery'
    if prime in parameters['csidh'].keys():
        # CSIDH only requires the factorization of p + 1
        L = parameters['csidh'][prime]['L']
        n = parameters['csidh'][prime]['n']
        montgomery_n = n
        cofactor = parameters['csidh'][prime]['cofactor']
        p = parameters['csidh'][prime]['p']
        p_minus_one_halves = parameters['csidh'][prime]['p_minus_one_halves']
        validation_stop = sum([bitlength(l_i) for l_i in L]) / 2.0 + 2
        field = PrimeField(p)
        path = resource_filename('sibc', "data/sdacs/csidh/" + prime)

    elif prime in parameters['bsidh'].keys():
        # B-SIDH only requires the factorization of p + 1 and p - 1
        nm = parameters['bsidh'][prime]['nm']
        np = parameters['bsidh'][prime]['np']
        Lp = parameters['bsidh'][prime]['Lp']
        Lm = parameters['bsidh'][prime]['Lm']
        Ep = parameters['bsidh'][prime]['Ep']
        Em = parameters['bsidh'][prime]['Em']
        L = list(Lp + Lm)
        n = len(L)
        montgomery_n = n
        p = parameters['bsidh'][prime]['p']
        cofactor_p = (p + 1) // reduce(lambda x, y: (x * y), Lp)
        cofactor_m = (p - 1) // reduce(lambda x, y: (x * y), Lm)
        validation_stop = sum([bitlength(l_i) for l_i in Lp]) / 2.0 + 2
        p_minus_one_halves = parameters['bsidh'][prime]['p_minus_one_halves']
        p_minus_3_quarters = parameters['bsidh'][prime]['p_minus_3_quarters']
        field = QuadraticField(p)
        path = resource_filename('sibc', "data/sdacs/bsidh/" + prime)

    elif prime in parameters['sidh'].keys():
        # SIDH only requires the exponents two and three such that p + 1 = 2^two x 3^three
        p = parameters['sidh'][prime]['p']
        two =  parameters['sidh'][prime]['two']
        three =  parameters['sidh'][prime]['three']
        p_minus_one_halves = parameters['sidh'][prime]['p_minus_one_halves']
        p_minus_3_quarters = parameters['sidh'][prime]['p_minus_3_quarters']
        field = QuadraticField(p)
        # Next path give a fixed list of sdacs (only for three and four),
        # but it is not required in sidh configuration
        L = [3,4]
        n = 2
        montgomery_n = n
        path = resource_filename('sibc', "data/sdacs/sidh")

    else:
        assert False, "only csidh, bsidh, and sidh are currently implemented"

    # Shortest Differential Addition Chains (SDACs) for each l_i
    SDACS = filename_to_list_of_lists_of_ints(path)
    assert len(SDACS) > 0, f'Not precomputed sdacs for {prime} prime'
    SDACS_LENGTH = list(map(len, SDACS))
    SDACS_REVERSED = list(map(lambda x:x[::-1], SDACS))

    cmul = lambda l: numpy.array(
        [
            4.0 * (SDACS_LENGTH[L.index(l)] + 2),
            2.0 * (SDACS_LENGTH[L.index(l)] + 2),
            6.0 * (SDACS_LENGTH[L.index(l)] + 2) - 2.0,
        ]
    )
    c_xmul = list(map(cmul, L))  # list of the costs of each [l]P

    SQR = 1.00
    ADD = 0.00
    def measure(x):
        """
        Field additions, multiplications, and squarings
        SQR = 1.00
        # In F_p, we have SQR_{F_p} = SQR x MUL_{F_p}
        ADD = 0.00
        # In F_p, we have ADD_{F_p} = ADD x MUL_{F_p}
        """
        return x[0] + SQR * x[1] + ADD * x[2]

    random = SystemRandom()

    def jinvariant(A):
        """
        -------------------------------------------------------------------------
        jinvariant()
        input : projective Montgomery constants A24 := A + 2C and C24 := 4C where
                E : y^2 = x^3 + (A/C)*x^2 + x
        output: the j-invariant of E
        -------------------------------------------------------------------------
        """
        A4_squared = (A[0] + A[0])              # (2 * A24)
        A4_squared -= A[1]        # (2 * A24) - C24
        A4_squared += A4_squared  # 4*A = 2[(2 * A24) - C24]

        # Now, we have A = A' / C' := A4 / A[1] = (4*A) / (4*C)
        A4_squared = (A4_squared ** 2)  # (A')^2
        C4_squared = (A[1] ** 2)        # (C')^2
        t = (C4_squared + C4_squared)   # 2 * [(C')^2]

        num = (C4_squared + t)      # 3 * [(C')^2]
        num = (A4_squared - num)    # (A')^2 - 3 * [(C')^2]
        s = (num ** 2)              # { (A')^2 - 3 * [(C')^2] }^2
        num *= s             # { (A')^2 - 3 * [(C')^2] }^3

        C4_squared = (C4_squared ** 2)  # (C')^4
        den = (t + t)                   # 4 * [(C')^2]
        den = (A4_squared - den)        # (A')^2 - 4 * [(C')^2]
        den *= C4_squared        # {(A')^2 - 4 * [(C')^2] } * [(C')^4]
        den **= -1               # 1 / {(A')^2 - 4 * [(C')^2] } * [(C')^4]

        num *= den  # j := { (A')^2 - 3 * [(C')^2] }^3 / {(A')^2 - 4 * [(C')^2] } * [(C')^4]
        num += num  #   2*j
        num += num  #   4*j
        num += num  #   8*j
        num += num  #  16*j
        num += num  #  32*j
        num += num  #  64*j
        num += num  # 128*j
        num += num  # 256*j
        return num

    def elligator(A):
        """ elligator() samples two points on E[pi + 1] or E[pi - 1] """
        Ap = (A[0] + A[0])
        Ap -= A[1]
        Ap += Ap
        Cp = A[1]

        # TODO random.randint is no good entropy source here
        u = field(random.randint(2, p_minus_one_halves))
        u_squared = (u ** 2)

        u_squared_plus_one = (u_squared + 1)
        u_squared_minus_one = (u_squared - 1)

        C_times_u_squared_minus_one = (Cp * u_squared_minus_one)
        AC_times_u_squared_minus_one = (Ap * C_times_u_squared_minus_one)

        tmp = (Ap ** 2)
        tmp *= u_squared
        aux = (C_times_u_squared_minus_one ** 2)
        tmp += aux
        tmp *= AC_times_u_squared_minus_one

        alpha, beta = field(0), u
        alpha, beta = cswap(alpha, beta, tmp == 0)
        u_squared_plus_one *= alpha
        alpha *= C_times_u_squared_minus_one

        Tp_X = (Ap + alpha)
        Tm_X = (Ap * u_squared)
        Tm_X += alpha
        Tm_X = (-Tm_X)

        tmp += u_squared_plus_one
        Tp_X, Tm_X = cswap(Tp_X, Tm_X, not tmp.issquare())

        return (
            [Tp_X, C_times_u_squared_minus_one],
            [Tm_X, C_times_u_squared_minus_one],
        )

    def affine_to_projective(affine):
        """
        affine_to_projective()
        input : the affine Montgomery coefficient A=A'/C with C=1
        output: projective Montgomery constants A24 := A' + 2C and C24 := 4C
                where E : y^2 = x^3 + (A'/C)*x^2 + x
        """
        return [affine + field(2), field(4)]

    def coeff(A):
        """
        ----------------------------------------------------------------------
        coeff()
        input : projective Montgomery constants A24 := A + 2C and C24 := 4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the affine Montgomery coefficient A/C
        ----------------------------------------------------------------------
        """
        output = (A[0] + A[0])         # (2 * A24)
        output -= A[1]                 # (2 * A24) - C24
        C24_inv = (A[1] ** -1)         # 1 / (C24)
        output += output               # 4*A = 2[(2 * A24) - C24]
        output *= C24_inv              # A/C = 2[(2 * A24) - C24] / C24

        return output

    def get_A(P, Q, PQ):
        """
        ----------------------------------------------------------------------
        coeff()
        input : the affine Montgomery x-coordinate points x(P) := (XP : YP), 
                x(Q) := (XQ : ZQ), and x(P - Q) := (XPQ : ZPQ) on the curve 
                E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery coefficient (A + 2C : 4C)
        ----------------------------------------------------------------------
        """
        t0 = (P[0] + P[1])   # XP + ZP
        t1 = (Q[0] + Q[1])   # XQ + ZQ

        t = (t0 * t1)        # (XP + ZP) * (XQ + ZQ)
        XPXQ = (P[0] * Q[0]) # XP * XQ
        ZPZQ = (P[1] * Q[1]) # ZP * ZQ

        t -= XPXQ
        t -= ZPZQ           # XPZQ + ZPXQ
        s = (XPXQ - ZPZQ)   # XPXQ - ZPZQ

        t0 = (t *  PQ[0])   # (XPZQ + ZPXQ) * XPQ
        t1 = (s *  PQ[1])   # (XPXQ - ZPZQ) * ZPQ
        t0 += t1            # (XPZQ + ZPXQ) * XPQ + (XPXQ - ZPZQ) * ZPQ
        t0 **= 2            # [(XPZQ + ZPXQ) * XPQ + (XPXQ - ZPZQ) * ZPQ] ^ 2

        t1 = (t * PQ[1])    # (XPZQ + ZPXQ) * ZPQ
        s = (ZPZQ * PQ[0])  # ZPZQ * XPQ
        t1 += s             # (XPZQ + ZPXQ) * ZPQ + ZPZQ * XPQ
        s = (XPXQ * PQ[0])  # (XPXQ) * XPQ
        s += s              # 2 * [(XPXQ) * XPQ]
        s += s              # 4 * [(XPXQ) * XPQ]
        t1 *= s             # [(XPZQ + ZPXQ) * ZPQ + ZPZQ * XPQ] * (4 * [(XPXQ) * XPQ])

        t = (ZPZQ * PQ[1])  # ZPZQ * ZPQ

        XPXQ = (t0 - t1)    # [(XPZQ + ZPXQ) * XPQ + (XPXQ - ZPZQ) * ZPQ] ^ 2 - [(XPZQ + ZPXQ) * ZPQ + ZPZQ * XPQ] * (4 * [(XPXQ) * XPQ])
        ZPZQ = (s * t)      # (4 * [(XPXQ) * XPQ]) * (ZPZQ * ZPQ)

        # ---
        B1 = (ZPZQ + ZPZQ)  # 2C
        B0 = (XPXQ + B1)    # A + 2C
        B1 += B1      # 4C
        return [B0, B1]

    def isinfinity(P):
        """ isinfinity(P) determines if x(P) := (XP : ZP) = (1 : 0) """
        return P[1] == 0

    def isequal(P, Q):
        """ isequal(P, Q) determines if x(P) = x(Q) """
        return (P[0] * Q[1]) == (P[1] * Q[0])

    def xdbl(P, A):
        """
        ----------------------------------------------------------------------
        xdbl()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, and
                the  projective Montgomery constants A24:= A + 2C and C24:=4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery x-coordinate point x([2]P)
        ----------------------------------------------------------------------
        """
        t_0 = (P[0] - P[1])
        t_0 **= 2
        Z = (A[1] * t_0)
        t_1 = (P[0] + P[1])
        t_1 **= 2
        X = (Z * t_1)
        t_1 -= t_0
        Z += (A[0] * t_1)
        Z *= t_1

        return [X, Z]
    try: # use specialized xdbl from field if available
        xdbl = field.xdbl
    except: pass

    def xadd(P, Q, PQ):
        """
        ----------------------------------------------------------------------
        xadd()
        input : the projective Montgomery x-coordinate points x(P) := XP/ZP,
                x(Q) := XQ/ZQ, and x(P-Q) := XPQ/ZPQ
        output: the projective Montgomery x-coordinate point x(P+Q)
        ----------------------------------------------------------------------
        """
        a = (P[0] + P[1])
        b = (P[0] - P[1])
        b *= (Q[0] + Q[1]) # c = (Q[0] + Q[1])
        a *= (Q[0] - Q[1]) # d = (Q[0] - Q[1])
        c = (a + b)
        a -= b # a = d
        c **= 2
        c *= PQ[1] # X = c
        a **= 2
        a *= PQ[0] # Z = d
        return [c, a]
    try: # use specialized xadd from the field if available
        xadd = field.xadd
    except: pass

    def xdbladd(P, Q, PQ, A):
        """
        ----------------------------------------------------------------------
        xdbladd()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, x(Q) := XQ/ZQ, 
                and x(P-Q) := XPQ/ZPQ, and the projective Montgomery constant (a24 : 1) 
                = (A + 2C: 4C)  where E : y^2 = x^3 + (A/C)x^2 + x
        output: the projective Montgomery x-coordinate point x([2]P), x([P+Q])
        ----------------------------------------------------------------------
        """
        t0 = (P[0] + P[1])
        t1 = (P[0] - P[1])
        X2 = (t0 ** 2)
        t2 = (Q[0] - Q[1])
        X3 = (Q[0] + Q[1])
        t0 *= t2
        Z2 = (t1 ** 2)
        # ---
        t1 *= X3
        t2 = (X2 - Z2)
        X2 *= Z2
        X3 = (A * t2)
        Z3 = (t0 - t1)
        Z2 += X3
        X3 = (t0 + t1)
        # ---
        Z2 *= t2
        Z3 **= 2
        X3 **= 2
        Z3 *= PQ[0]
        X3 *= PQ[1]
        return [X2, Z2], [X3, Z3]

    def xtpl(P, A):
        """
        ----------------------------------------------------------------------
        xtpl()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, and
                the  projective Montgomery constants A24p:=A+2C and A24m:=A-2C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery x-coordinate point x([3]P)
        ----------------------------------------------------------------------
        """
        # ---
        t0 = (P[0] - P[1])
        t2 = (t0 ** 2)
        t1 = (P[0] + P[1])
        t3 = (t1 ** 2) # (P0*P1)^2
        t4 = (t1 + t0) # P0*P1 + (P0-P1)^2
        t7 = (t1 - t0) # P0*P1 - (P0-P1)^2
        # ---
        t11 = (t4 ** 2) # (P0*P1 + (P0-P1)^2)^
        t11 -= t3
        t11 -= t2
        t5 = (t3 * A[0])
        t3 *= t5
        t6 = (t2 * A[1])
        # ---
        t2 *= t6
        t8 = (t2 - t3)
        t5 -= t6
        t11 *= t5
        t10 = (t8 + t11)
        t10 **= 2
        # ---
        t10 *= t4 # X = t10 * (P0*P1 + (P0-P1)^2)
        t8 -= t11 # Z
        t8 **= 2
        t8 *= t7 # Z
        return [t10,t8] # X, Z

    def xdbl_twice(P, A):
        """
        xdbl_twice()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, and
                the  projective Montgomery constants A24:= A + 2C and C24:=4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery x-coordinate point x([4]P)
        ----------------------------------------------------------------------
        """
        return xdbl(xdbl(P,A), A)

    def xmul(P, A, j):
        """
        ----------------------------------------------------------------------
        xmul()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and an positive integer j
        output: the projective Montgomery x-coordinate point x([L[j]]P)
        ----------------------------------------------------------------------
        """
        P2 = xdbl(P, A)
        R = [P, P2, xadd(P2, P, P)]
        try: # fast-track for PrimeField
            P[0].x
            def isinf(POINT): return POINT[1].x == 0
        except:
            isinf = isinfinity

        for sdac in SDACS_REVERSED[j]:

            #XXXX isinfinity inlined for speed here:
            if isinf(R[sdac]): # if isinfinity(R[sdac]):
                R[:] = R[sdac ^ 1], R[2], xdbl(R[2], A)
            else:
                R[:] = R[sdac ^ 1], R[2], xadd(R[2], R[sdac ^ 1], R[sdac])

        return R[2]

    def Ladder3pt(m, P, Q, PQ, A):
        """
        ----------------------------------------------------------------------
        Ladder3pt()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP,
                x(Q) := XQ/ZQ, and x(P-Q) := XPQ/ZPQ, the affine Montgomery
                constant A where E : y^2 = x^3 + Ax^2 + x, and a positive
                integer m
        output: the projective Montgomery x-coordinate point x(P + [m]Q)
        ----------------------------------------------------------------------
        """
        X0 = list([field(Q[0]), field(Q[1])])
        X1 = list([field(P[0]), field(P[1])])
        X2 = list([field(PQ[0]), field(PQ[1])])
        t = 0x1
        for i in range(0, bitlength(p), 1):
            X1, X2 = cswap(X1, X2, t & m == 0)
            X0, X1 = xdbladd(X0, X1, X2, A)
            X1, X2 = cswap(X1, X2, t & m == 0)
            t <<= 1
        return X1

    # Golden ration is used in prac algorithm
    phi = (1.0 + sqrt(0.5)) / 2.0
    phi_nu, phi_de = phi.as_integer_ratio()
    def euclid2d(m, n, P, Q, PQ, A):
        """ The 2-dimensional scalar pseudomultiplication: x([r]P + [s - r]P) with r = s / {Golden Ratio}') """

        s0, s1 = m+0, n+0
        x0  = list(P)
        x1  = list(Q)
        diff= list(PQ)

        while s0 != 0:
            if s1 < s0:
                x0, x1 = cswap(x0, x1, 1)
                s0, s1 = s1, s0
            if s1 <= 4*s0:
                # Fibonacci step
                x   = list(x0)
                x0  = xadd(x1, x0, diff)
                diff= list(x)
                s1 -= s0
            elif (s0 % 2) == (s1 % 2):
                x0 = xadd(x1, x0, diff)
                x1 = xdbl(x1, A)
                s1 -= s0
                s1 //= 2
            elif (s1 % 2) == 0:
                diff=xadd(x1, diff, x0)
                x1  =xdbl(x1, A)
                s1 //= 2
            else:
                diff= xadd(x0, diff, x1)
                x0  = xdbl(x0, A)
                s0 //= 2

        while s1 % 2 == 0:
            x1 = xdbl(x1, A)
            s1 //= 2

        if s1 > 1:
            # Ladder step on the missing part: x0 will correspond with Ladder(x1)
            diff= list(x1)
            x0  = xdbl(x1, A)
            s1_binary = bin(s1)[2:][::-1]
            s1_length = len(s1_binary)
            for i in range(s1_length - 2, -1, -1):
                x0, x1 = cswap(x0, x1, int(s1_binary[i + 1]) ^ int(s1_binary[i]))
                x1 = xadd(x0, x1, diff)
                x0 = xdbl(x0, A)

            x0, x1 = cswap(x0, x1, int(s1_binary[0]))
        else:
            # In this case, the output should correspond with x1, thus we swap to x0
            x0, x1 = cswap(x0, x1, 1)

        return x0

    def prac(k, P, A):
        """ PRAC algorithm: (simplified) 1-D Euclidean pseudomultiplication """

        s = k
        infty = [field(1), field(0)]  # Point at infinity

        # Reducing the problem from k = 2^i x s to s
        x = list(P)
        while s & 1 == 0: # % 2
            x = xdbl(x, A)
            s //= 2

        r = (s * phi_nu) // phi_de
        x = euclid2d(r, s - r, x, x, infty, A)
        return x

    def cofactor_multiples(P, A, points):
        """
        ----------------------------------------------------------------------
        cofactor_multiples()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and subset of |[0, n]|
        output: the projective Montgomery x-coordinate points x([(p+1) / l_0]P),
                x([(p+1) / l_1]P), ..., x([(p+1) / l_{n-1}]P).
        ----------------------------------------------------------------------
        """
        n = len(points)
        if n == 1:
            # In this recursion level we have an order-l point
            return [P]
        elif n > 0:
            # We proceed by applying a divide-and-conquer procedure
            h = n // 2
            if h > 0:

                # 1st half
                first_half = []
                second_P = P
                for j in range(h):
                    second_P = xmul(second_P, A, points[j])
                    first_half.append(points[j])

                # 2nd half
                second_half = []
                first_P = P
                for j in range(h, n):
                    first_P = xmul(first_P, A, points[j])
                    second_half.append(points[j])

                return cofactor_multiples(first_P, A, first_half) + cofactor_multiples(
                    second_P, A, second_half
                )

            return []

    def isfullorder(seq):
        """ isfullorder() checks if the point at infinity belongs or not to a given list """
        tmp = [not isinfinity(seq_i) for seq_i in seq]
        return reduce(lambda x, y: (x and y), tmp)

    def generators(A):
        """ generators() looks for two full-order poinst on E[pi + 1] or E[pi - 1] """
        output = [[0, 0], [0, 0]]
        while [0, 0] in output:

            T_p, T_m = elligator(A)
            T_p = prac(cofactor, T_p, A)
            if isfullorder(cofactor_multiples(T_p, A, range(0, n, 1))) and output[
                0
            ] == [0, 0]:
                output[0] = list(T_p)

            T_m = prac(cofactor, T_m, A)
            if isfullorder(cofactor_multiples(T_m, A, range(0, n, 1))) and output[
                1
            ] == [0, 0]:
                output[1] = list(T_m)

        return output[0], output[1]

    def crisscross(alpha, beta, gamma, delta):
        """ crisscross() computes a*c + b*d, and a*c - b*d """
        t_1 = (alpha * delta)
        t_2 = (beta * gamma)
        #return (t_1 + t_2), (t_1 - t_2)
        # shave off a FF allocation: ##
        t_3 = t_1.copy() # object.__new__(t_1.__class__); t_3.x = t_1.x #      ## copy(t_1)
        t_1 += t_2                    ##
        t_3 -= t_2                   ##
        return t_1, t_3 # (t_1 + t_2), (t_1 - t_2)

    if prime in parameters['csidh'].keys():
        def issupersingular(A):
            """ issupersingular() verifies supersingularity """
            T_p, _ = elligator(A) # T_p is always in GF(p), and thus has torsion (p+1)
            T_p = prac(cofactor, T_p, A)
            assert n == montgomery_n
            P = cofactor_multiples(T_p, A, range(0, montgomery_n, 1))

            assert n == montgomery_n
            bits_of_the_order = 0
            for i in range(0, montgomery_n, 1):

                if isinfinity(P[i]) == False:

                    Q = xmul(P[i], A, i)

                    if isinfinity(Q) == False:
                        return False

                    bits_of_the_order += bitlength(L[i])
                    if bits_of_the_order > validation_stop:
                        return True

    elif prime in parameters['bsidh'].keys():
        def issupersingular(A):
            """ issupersingular() verifies supersingularity """
            P, _ = elligator(A) # T_p is always in GF(p²)

            # Is P a torsion-(p + 1)?
            T_p = prac(cofactor_p, P, A)
            Tp = cofactor_multiples(T_p, A, range(0, np, 1))
            Tp = [xmul(Tp[i], A, i) for i in range(0, np, 1)]
            Tp = [isinfinity(Tp_i) for Tp_i in Tp]

            assert n == montgomery_n
            # Is P a torsion-(p - 1)?
            T_m = prac(cofactor_m, P, A)
            Tm = cofactor_multiples(T_m, A, range(np, n, 1))
            Tm = [xmul(Tm[i - np], A, i) for i in range(np, n, 1)]
            Tm = [isinfinity(Tm_i) for Tm_i in Tm]

            return reduce(lambda x, y: (x and y), Tp) or reduce(lambda x, y: (x and y), Tm)
    elif prime in parameters['sidh'].keys():
        def issupersingular(A):
            """ issupersingular() verifies supersingularity """
            x = field(random.randint(1, p))
            T = [x, field(1)]
            for i in range(0, two):
                T = xdbl(T, A)
            for i in range(0, three):
                T = xtpl(T, A)
            
            return isinfinity(T)
    else:
        assert False, "only CSIDH, CSIDH, and B-SIDH are currently implemented"

    return attrdict(**locals())
