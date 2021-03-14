import math
from random import SystemRandom

import numpy

# Loading the library corresponding to the arithmetic in GF(p)
from sidh.primefield import PrimeField, cswap   # <--- cswap should be later imported from common
#from sidh.quadraticfield import QuadraticField # <--- fp2 classes not integrated for now

from functools import reduce

from sidh.constants import sdacs_data, bitlength, parameters
from sidh.common import attrdict


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
            fh.writelines(' '.join(str(v) for v in v in line))
        fh.writelines()


def MontgomeryCurve(prime, style):

    if True:  # algorithm == 'csidh':
        # this is montgomery.py currently
        A = parameters['csidh']['A']
        L = parameters['csidh'][prime]['L']
        n = parameters['csidh'][prime]['n']
        exponent_of_two = parameters['csidh'][prime]['exponent_of_two']

        p = (2 ** (exponent_of_two)) * reduce(
            lambda x, y: (x * y), L
        ) - 1  # p := 4 * l_0 * ... * l_n - 1
        # p_minus_one_halves = (p - 1) // 2  # (p - 1) / 2
        p_minus_one_halves = parameters['csidh'][prime]['p_minus_one_halves']
        validation_stop = sum([bitlength(l_i) for l_i in L]) / 2.0 + 2
        FiniteFieldClass = PrimeField
    else:
        assert False, "bsidh not refactored yet"
        # this is for a possible future where we have a unified montgomery.py
        # for csidh and bsidh. the code in this else block was moved from fp.py

        f = open(sop_data + prime)

        # The prime to be used
        p = f.readline()
        self.p = int(p, 16)

        # List corresponding (p + 1)
        Lp = f.readline()
        Lp = [int(lp) for lp in Lp.split()]
        # exponent_of_twop = Lp[0]
        # Lp = Lp[1:]
        Ep = f.readline()
        Ep = [int(ep) for ep in Ep.split()]
        assert len(Ep) == len(Lp)
        np = len(Lp)

        # List corresponding (p - 1)
        Lm = f.readline()
        Lm = [int(lm) for lm in Lm.split()]
        Em = f.readline()
        Em = [int(em) for em in Em.split()]
        assert len(Em) == len(Lm)
        nm = len(Lm)

        f.close()

        # pp = (2**exponent_of_twop) * reduce(lambda x,y : (x*y), [ Lp[i]**Ep[i] for i in range(0, np, 1)  ])
        pp = reduce(
            lambda x, y: (x * y), [Lp[i] ** Ep[i] for i in range(0, np, 1)]
        )
        pm = reduce(
            lambda x, y: (x * y), [Lm[i] ** Em[i] for i in range(0, nm, 1)]
        )
        assert (p + 1) % pp == 0
        assert (p - 1) % pm == 0

        if p % 4 == 1:
            print("// Case p = 1 mod 4 is not implemented yet!")
            exit(-1)

        p_minus_one_halves = (p - 1) // 2
        p_minus_3_quarters = (p - 3) // 4
        #FiniteFieldClass = QuadraticField

    ff = FiniteFieldClass(p)
    # print("// Shortest Differential Addition Chains (SDAC) for each l_i;")
    # List of Small odd primes, L := [l_0, ..., l_{n-1}]
    # print("// SDAC's to be read from a file")
    path = sdacs_data + prime
    SDACS = filename_to_list_of_lists_of_ints(path)
    if len(SDACS) == 0:
        print("// SDAC's to be computed")
        SDACS = generate_sdacs(L)
        print("// Storing SDAC's in a file")
        write_list_of_lists_of_ints_to_file(path, SDACS)
    SDACS_LENGTH = list(map(len, SDACS))

    cMUL = lambda l: numpy.array(
        [
            4.0 * (SDACS_LENGTH[L.index(l)] + 2),
            2.0 * (SDACS_LENGTH[L.index(l)] + 2),
            6.0 * (SDACS_LENGTH[L.index(l)] + 2) - 2.0,
        ]
    )
    C_xmul = list(map(cMUL, L))  # list of the costs of each [l]P

    SQR = 1.00
    ADD = 0.00

    random = SystemRandom()

    def elligator(A):

        Ap = (A[0] + A[0])
        Ap = (Ap - A[1])
        Ap = (Ap + Ap)
        Cp = A[1]

        u = ff(random.randint(2, p_minus_one_halves))
        u_squared = (u ** 2)

        u_squared_plus_one = (u_squared + 1)
        u_squared_minus_one = (u_squared - 1)

        C_times_u_squared_minus_one = (Cp * u_squared_minus_one)
        AC_times_u_squared_minus_one = (Ap * C_times_u_squared_minus_one)

        tmp = (Ap ** 2)
        tmp = (tmp * u_squared)
        aux = (C_times_u_squared_minus_one ** 2)
        tmp = (tmp + aux)
        tmp = (AC_times_u_squared_minus_one * tmp)

        alpha, beta = 0, u
        alpha, beta = cswap(alpha, beta, tmp == 0)
        u_squared_plus_one = (alpha * u_squared_plus_one)
        alpha = (alpha * C_times_u_squared_minus_one)

        Tp_X = (Ap + alpha)
        Tm_X = (Ap * u_squared)
        Tm_X = (Tm_X + alpha)
        Tm_X = (-Tm_X)

        tmp = (tmp + u_squared_plus_one)
        Tp_X, Tm_X = cswap(Tp_X, Tm_X, not tmp.issquare())

        return (
            [Tp_X, C_times_u_squared_minus_one],
            [Tm_X, C_times_u_squared_minus_one],
        )

    def generate_sdacs(L):
        return list(
            map(sdac, L)
        )  # Shortest Differential Addition Chains for each small odd prime l in L

    def measure(x):
        """
        SQR = 1.00
        # In F_p, we have SQR_{F_p} = SQR x MUL_{F_p}
        ADD = 0.00
        # In F_p, we have ADD_{F_p} = ADD x MUL_{F_p}
        """
        return x[0] + SQR * x[1] + ADD * x[2]

    def dacs(l, r0, r1, r2, chain):
        '''
        dacs()
        inputs: a small odd prime number l, three integer numbers, and a list
        output: all the differential additions chains corresponding with the input l

        NOTE: this is a recursive approach
        '''
        if r2 == l:

            return [(chain, r2)]
        elif r2 < l and len(chain) <= 1.5 * math.log(l, 2):

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

    def affine_to_projective(affine):
        """
        input : the affine Montgomery coefficient A=A'/C with C=1
        output: projective Montgomery constants A24 := A' + 2C and C24 := 4C
                where E : y^2 = x^3 + (A'/C)*x^2 + x
        """
        return [affine + ff(2), ff(4)]

    def coeff(A):
        '''
        ----------------------------------------------------------------------
        coeff()
        input : projective Montgomery constants A24 := A + 2C and C24 := 4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the affine Montgomery coefficient A/C
        ----------------------------------------------------------------------
        '''
        output = (A[0] + A[0])         # (2 * A24)
        output = (output - A[1])       # (2 * A24) - C24
        C24_inv = (A[1] ** -1)         # 1 / (C24)
        output = (output + output)     # 4*A = 2[(2 * A24) - C24]
        output = (output * C24_inv)    # A/C = 2[(2 * A24) - C24] / C24

        return output

    def isinfinity(P):
        """
        isinfinity(P) determines if x(P) := (XP : ZP) = (1 : 0)
        """
        return P[1] == 0

    def isequal(P, Q):
        """ isequal(P, Q) determines if x(P) = x(Q) """
        return (P[0] * Q[1]) == (P[1] * Q[0])

    def xdbl(P, A):
        '''
        ----------------------------------------------------------------------
        xdbl()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, and
                the  projective Montgomery constants A24:= A + 2C and C24:=4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery x-coordinate point x([2]P)
        ----------------------------------------------------------------------
        '''
        t_0 = (P[0] - P[1])
        t_1 = (P[0] + P[1])
        t_0 = (t_0 ** 2)
        t_1 = (t_1 ** 2)
        Z = (A[1] * t_0)
        X = (Z * t_1)
        t_1 = (t_1 - t_0)
        t_0 = (A[0] * t_1)
        Z = (Z + t_0)
        Z = (Z * t_1)

        return [X, Z]

    def xadd(P, Q, PQ):
        '''
        ----------------------------------------------------------------------
        xadd()
        input : the projective Montgomery x-coordinate points x(P) := XP/ZP,
                x(Q) := XQ/ZQ, and x(P-Q) := XPQ/ZPQ
        output: the projective Montgomery x-coordinate point x(P+Q)
        ----------------------------------------------------------------------
        '''
        a = (P[0] + P[1])
        b = (P[0] - P[1])
        c = (Q[0] + Q[1])
        d = (Q[0] - Q[1])
        a = (a * d)
        b = (b * c)
        c = (a + b)
        d = (a - b)
        c = (c ** 2)
        d = (d ** 2)
        X = (PQ[1] * c)
        Z = (PQ[0] * d)
        return [X, Z]

    # Modificar esta parte para usar cadenas de addicion
    def xmul(P, A, j):
        '''
        ----------------------------------------------------------------------
        xmul()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and an positive integer j
        output: the projective Montgomery x-coordinate point x([L[j]]P)
        ----------------------------------------------------------------------
        '''
        P2 = xdbl(P, A)
        R = [P, P2, xadd(P2, P, P)]

        for i in range(SDACS_LENGTH[j] - 1, -1, -1):

            if isinfinity(R[SDACS[j][i]]):
                T = xdbl(R[2], A)
            else:
                T = xadd(R[2], R[SDACS[j][i] ^ 1], R[SDACS[j][i]])

            R[0] = list(R[SDACS[j][i] ^ 1])
            R[1] = list(R[2])
            R[2] = list(T)

        return R[2]

    def cofactor_multiples(P, A, points):
        '''
        ----------------------------------------------------------------------
        cofactor_multiples()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and subset of |[0, n]|
        output: the projective Montgomery x-coordinate points x([(p+1) / l_0]P),
                x([(p+1) / l_1]P), ..., x([(p+1) / l_{n-1}]P).
        ----------------------------------------------------------------------
        '''
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
        tmp = [not isinfinity(seq_i) for seq_i in seq]
        return reduce(lambda x, y: (x and y), tmp)

    def generators(A):

        output = [[0, 0], [0, 0]] if style != 'wd1' else [[0, 0], [1, 1]]
        while [0, 0] in output:

            T_p, T_m = elligator(A)
            for i in range(0, exponent_of_two, 1):
                T_p = xdbl(T_p, A)

            if isfull_order(cofactor_multiples(T_p, A, range(0, n, 1))) and output[
                0
            ] == [0, 0]:
                output[0] = list(T_p)

            if style != 'wd1':
                for i in range(0, exponent_of_two, 1):
                    T_m = xdbl(T_m, A)
                if isfull_order(
                    cofactor_multiples(T_m, A, range(0, n, 1))
                ) and output[1] == [0, 0]:
                    output[1] = list(T_m)

        return output[0], output[1]

    def crisscross(alpha, beta, gamma, delta):

        t_1 = (alpha * delta)
        t_2 = (beta * gamma)
        return (t_1 + t_2), (t_1 - t_2)

    def issupersingular(A):

        while True:

            T_p, _ = elligator(A)
            for i in range(0, exponent_of_two, 1):
                T_p = xdbl(T_p, A)

            P = cofactor_multiples(T_p, A, range(0, n, 1))

            bits_of_the_order = 0
            for i in range(0, n, 1):

                if isinfinity(P[i]) == False:

                    Q = xmul(P[i], A, i)

                    if isinfinity(Q) == False:
                        return False

                    bits_of_the_order += bitlength(L[i])
                    if bits_of_the_order > validation_stop:
                        return True

    return attrdict(**locals())
