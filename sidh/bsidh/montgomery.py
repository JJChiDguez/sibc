import math
from random import SystemRandom
from pkg_resources import resource_filename

from functools import reduce
import numpy

# Loading the library corresponding to the arithmetic in F_{p^2}
from sidh.fp2 import F_p2
from sidh.constants import sdacs_data, bitlength, parameters
from sidh.math import hamming_weight
from sidh.common import attrdict

bitlength = lambda x: len(bin(x)[2:])  # number of bits
hamming_weight = lambda x: bin(x).count("1")
# hamming weight: number of bits equal 1

def filename_to_list_of_lists_of_ints(path):
    res = []
    try:
        with open(path) as fh:
            for line in fh:
                res.append(list(map(int, line.split())))
    except:
        res = []
    return res


def write_list_of_lists_of_ints_to_file(path, data):
    with open(path, 'w') as fh:
        for line in data:
            fh.writelines(' '.join(str(v) for v in line))
        fh.writelines('')


def MontgomeryCurve(prime):
    A = parameters['bsidh']['A']
    nm = parameters['bsidh'][prime]['nm']
    np = parameters['bsidh'][prime]['np']
    Lp = parameters['bsidh'][prime]['Lp']
    Lm = parameters['bsidh'][prime]['Lm']
    Ep = parameters['bsidh'][prime]['Ep']
    Em = parameters['bsidh'][prime]['Em']
    p = parameters['bsidh'][prime]['p']


    fp = F_p2(p)

    SQR = 1.00
    # In F_p, we have SQR_{F_p} = SQR x MUL_{F_p}
    ADD = 0.00
    # In F_p, we have ADD_{F_p} = ADD x MUL_{F_p}
    measure = lambda x: (x[0] + SQR * x[1] + ADD * x[2])

    def dacs(l, r0, r1, r2, chain):
        '''
        dacs()
        inputs: a small odd prime number l, three integer numbers, and a list
        output: all the differential additions chains corresponding with the input l

        NOTE: this is a recursive approach
        '''
        if r2 == l:

            return [(chain, r2)]
        elif r2 < l and len(chain) <= 1.5 * log(l, 2):

            return dacs(l, r0, r2, r2 + r0, chain + [1]) + dacs(
                l, r1, r2, r2 + r1, chain + [0]
            )
        else:
            return []

    def sdac(l):
        '''
        sdac()
        input: a small odd prime number l
        output: the shortest differential additions chains corresponding with the input l

        NOTE: this function uses a recursive function
        '''
        all_dacs = dacs(l, 1, 2, 3, [])
        return min(all_dacs, key=lambda t: len(t[0]))[0]

    path = resource_filename(__name__, "../data/sdacs/" + prime)
    SDACS = filename_to_list_of_lists_of_ints(path)
    if len(SDACS) == 0:
        print("// SDAC's to be computed")
        SDACS = list(
            map(sdac, Lp)
        )  # Shortest Differential Addition Chains for each small odd prime l in Lp
        # Next nm SDACS
        SDACS = SDACS + list(
            map(sdac, Lm)
        )  # Shortest Differential Addition Chains for each small odd prime l in Lm
        #print("// Storing SDAC's in a file")
        #write_list_of_lists_of_ints_to_file(path, SDACS)
    SDACS_LENGTH = list(map(len, SDACS))

       ##   if sys.argv[0] != 'header.py':
       ##       print("// Shortest Differential Addition Chains (SDAC) for each l_i;")
       #try:

       #    # List of Small odd primes, L := [l_0, ..., l_{n-1}]
       #    f = open(sdacs_data + setting.prime)
       #    if sys.argv[0] != 'header.py':
       #        print("// SDAC's to be read from a file")

       #    SDACS = []
       #    SDACS_LENGTH = []
       #    for i in range(0, np + nm, 1):

       #        tmp = f.readline()
       #        tmp = [int(b) for b in tmp.split()]
       #        SDACS.append(tmp)
       #        SDACS_LENGTH.append(len(tmp))

       #    f.close()

       #except IOError:

       #    if sys.argv[0] != 'header.py':
       #        print("// SDAC's to be computed")

       #    # First np SDACS
       #    SDACS = list(
       #        map(sdac, Lp)
       #    )  # Shortest Differential Addition Chains for each small odd prime l in Lp
       #    SDACS_LENGTH = [len(sdacs) for sdacs in SDACS]
       #    # Next nm SDACS
       #    SDACS = SDACS + list(
       #        map(sdac, Lm)
       #    )  # Shortest Differential Addition Chains for each small odd prime l in Lm
       #    SDACS_LENGTH = [len(sdacs) for sdacs in SDACS]

       #    # Storing the SDACS to be used
       #    print("// Storing SDAC's in a file")
       #    f = open(sdacs_data + setting.prime, 'w')
       #    for i in range(0, np + nm, 1):
       #        f.writelines(' '.join([str(tmp) for tmp in SDACS[i]]) + '\n')

       #    f.close()

    # List of Small Isogeny Degree, Lp := [l_0, ..., l_{n-1}] [We need to include case l=2 and l=4]
    SIDp = reduce(lambda x, y: x + y, [[Lp[i]] * Ep[i] for i in range(0, np, 1)])
    # List of Small Isogeny Degree, Lm := [l_0, ..., l_{n-1}]
    SIDm = reduce(lambda x, y: x + y, [[Lm[i]] * Em[i] for i in range(0, nm, 1)])

    SID = list(SIDp + SIDm)
    n = len(SID)

    L = global_L = list(Lp + Lm)

    cMUL = lambda l: numpy.array(
        [
            4.0 * (SDACS_LENGTH[global_L.index(l)] + 2),
            2.0 * (SDACS_LENGTH[global_L.index(l)] + 2),
            6.0 * (SDACS_LENGTH[global_L.index(l)] + 2) - 2.0,
        ]
    )
    C_xMUL = list(map(cMUL, global_L))  # list of the costs of each [l]P


    def coeff(A):
        '''
        ----------------------------------------------------------------------
        coeff()
        input : projective Montgomery constants A24 := A + 2C and C24 := 4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the affine Montgomery coefficient A/C
        ----------------------------------------------------------------------
        '''
        output = fp.fp2_add(A[0], A[0])  # (2 * A24)
        output = fp.fp2_sub(output, A[1])  # (2 * A24) - C24
        C24_inv = fp.fp2_inv(A[1])  # 1 / (C24)
        output = fp.fp2_add(output, output)  # 4*A = 2[(2 * A24) - C24]
        output = fp.fp2_mul(output, C24_inv)  # A/C = 2[(2 * A24) - C24] / C24
        return output

    def jinvariant(A):
        '''
        ----------------------------------------------------------------------
        jinvariant() input : projective Montgomery constants A24 := A + 2C and C24
        := 4C where E : y^2 = x^3 + (A/C)*x^2 + x output: the j-invariant of E
        ----------------------------------------------------------------------
        '''
        A4_squared = fp.fp2_add(A[0], A[0])  # (2 * A24)
        A4_squared = fp.fp2_sub(A4_squared, A[1])  # (2 * A24) - C24
        A4_squared = fp.fp2_add(A4_squared, A4_squared)  # 4*A = 2[(2 * A24) - C24]

        # Now, we have A = A' / C' := A4 / A[1] = (4*A) / (4*C)
        A4_squared = fp.fp2_sqr(A4_squared)  # (A')^2
        C4_squared = fp.fp2_sqr(A[1])  # (C')^2
        t = fp.fp2_add(C4_squared, C4_squared)  # 2 * [(C')^2]

        num = fp.fp2_add(C4_squared, t)  # 3 * [(C')^2]
        num = fp.fp2_sub(A4_squared, num)  # (A')^2 - 3 * [(C')^2]
        s = fp.fp2_sqr(num)  # { (A')^2 - 3 * [(C')^2] }^2
        num = fp.fp2_mul(num, s)  # { (A')^2 - 3 * [(C')^2] }^3

        C4_squared = fp.fp2_sqr(C4_squared)  # (C')^4
        den = fp.fp2_add(t, t)  # 4 * [(C')^2]
        den = fp.fp2_sub(A4_squared, den)  # (A')^2 - 4 * [(C')^2]
        den = fp.fp2_mul(den, C4_squared)  # {(A')^2 - 4 * [(C')^2] } * [(C')^4]
        den = fp.fp2_inv(den)  # 1 / {(A')^2 - 4 * [(C')^2] } * [(C')^4]

        num = fp.fp2_mul(
            num, den
        )  # j := { (A')^2 - 3 * [(C')^2] }^3 / {(A')^2 - 4 * [(C')^2] } * [(C')^4]
        num = fp.fp2_add(num, num)  #   2*j
        num = fp.fp2_add(num, num)  #   4*j
        num = fp.fp2_add(num, num)  #   8*j
        num = fp.fp2_add(num, num)  #  16*j
        num = fp.fp2_add(num, num)  #  32*j
        num = fp.fp2_add(num, num)  #  64*j
        num = fp.fp2_add(num, num)  # 128*j
        num = fp.fp2_add(num, num)  # 256*j
        return num

    def isinfinity(P):
        """
        isinfinity(P) determines if x(P) := (XP : ZP) = (1 : 0)
        """
        return (P[1][0] == 0) and (P[1][1] == 0)

    def areequal(P, Q):
        """
        areequal(P, Q) determines if x(P) = x(Q)
        """
        XPZQ = fp.fp2_mul(P[0], Q[1])
        ZPXQ = fp.fp2_mul(P[1], Q[0])
        return (XPZQ[0] == ZPXQ[0]) and (XPZQ[0] == ZPXQ[0])

    def xDBL(P, A):
        '''
        ----------------------------------------------------------------------
        xDBL()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, and
                the  projective Montgomery constants A24:= A + 2C and C24:=4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery x-coordinate point x([2]P)
        ----------------------------------------------------------------------
        '''
        t_0 = fp.fp2_sub(P[0], P[1])
        t_1 = fp.fp2_add(P[0], P[1])
        t_0 = fp.fp2_sqr(t_0)
        t_1 = fp.fp2_sqr(t_1)
        Z = fp.fp2_mul(A[1], t_0)
        X = fp.fp2_mul(Z, t_1)
        t_1 = fp.fp2_sub(t_1, t_0)
        t_0 = fp.fp2_mul(A[0], t_1)
        Z = fp.fp2_add(Z, t_0)
        Z = fp.fp2_mul(Z, t_1)
        return [X, Z]

    def xADD(P, Q, PQ):
        '''
        ----------------------------------------------------------------------
        xADD()
        input : the projective Montgomery x-coordinate points x(P) := XP/ZP,
                x(Q) := XQ/ZQ, and x(P-Q) := XPQ/ZPQ
        output: the projective Montgomery x-coordinate point x(P+Q)
        ----------------------------------------------------------------------
        '''
        a = fp.fp2_add(P[0], P[1])
        b = fp.fp2_sub(P[0], P[1])
        c = fp.fp2_add(Q[0], Q[1])
        d = fp.fp2_sub(Q[0], Q[1])
        a = fp.fp2_mul(a, d)
        b = fp.fp2_mul(b, c)
        c = fp.fp2_add(a, b)
        d = fp.fp2_sub(a, b)
        c = fp.fp2_sqr(c)
        d = fp.fp2_sqr(d)
        X = fp.fp2_mul(PQ[1], c)
        Z = fp.fp2_mul(PQ[0], d)
        return [X, Z]

    def xDBLADD(P, Q, PQ, A):
        """
        Next function computes x([2]P) and x(P + Q)
        """

        # In C-code implementations this can be optimized
        S = xADD(P, Q, PQ)
        T = xDBL(P, A)
        return T, S

    # Modificar esta parte para usar cadenas de addicion
    def xMUL(P, A, j):
        '''
        ----------------------------------------------------------------------
        xMUL()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and an positive integer j
        output: the projective Montgomery x-coordinate point x([L[j]]P)
        ----------------------------------------------------------------------
        '''

        P2 = xDBL(P, A)
        R = [P, P2, xADD(P2, P, P)]

        for i in range(SDACS_LENGTH[j] - 1, -1, -1):

            T = xADD(R[2], R[SDACS[j][i] ^ 1], R[SDACS[j][i]])

            R[0] = list(R[SDACS[j][i] ^ 1])
            R[1] = list(R[2])
            R[2] = list(T)

        return R[2]

    def Ladder3pt(m, P, Q, PQ, A):
        """
        Next function computes x(P + [m]Q)
        """
        X0 = list([list(Q[0]), list(Q[1])])
        X1 = list([list(P[0]), list(P[1])])
        X2 = list([list(PQ[0]), list(PQ[1])])
        t = 0x1
        for i in range(0, bitlength(p), 1):
            # In C-code implementations this branch should be implemented with cswap's
            if t & m != 0:
                X0, X1 = xDBLADD(X0, X1, X2, A)
            else:
                X0, X2 = xDBLADD(X0, X2, X1, A)
            t <<= 1
        return X1

    def prime_factors(P, A, points):
        """
        ----------------------------------------------------------------------
        prime_factors()
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
                    second_P = xMUL(second_P, A, points[j])
                    first_half.append(points[j])
                # 2nd half
                second_half = []
                first_P = P
                for j in range(h, n):
                    first_P = xMUL(first_P, A, points[j])
                    second_half.append(points[j])

                return prime_factors(first_P, A, first_half) + prime_factors(
                    second_P, A, second_half
                )
            return []

    def random_affine_point(A):
        """
        Random affine point (x,y) in E : y^2 = x^3 + Ax^2 + x
        """
        return "Not implemented yet!"

    def difference_point(P, Q, A):
        """
        Next function computes x(P - Q) giving two affine points P and Q
        """
        return "Not implemented yet!"

    def isfull_order(seq):
        tmp = [not isinfinity(seq_i) for seq_i in seq]
        return reduce(lambda x, y: (x and y), tmp)

    def full_torsion_points(A):
        return "Not implemented yet!"

    def CrissCross(alpha, beta, gamma, delta):
        t_1 = fp.fp2_mul(alpha, delta)
        t_2 = fp.fp2_mul(beta, gamma)
        return fp.fp2_add(t_1, t_2), fp.fp2_sub(t_1, t_2)

    return attrdict(**locals())
