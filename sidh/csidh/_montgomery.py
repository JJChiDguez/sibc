import math
import random

#
import numpy

# Loading the library corresponding to the arithmetic in F_p
from sidh._fp import *
from sidh.constants import sdacs_data, bitlength, parameters
from sidh._math import hamming_weight

class MontgomeryLadder(object):
    def __init__(self, prime, style):
        self.style = style
        self.prime = prime
        self.fp = F_p('csidh', prime)
        self.L = L = parameters['csidh'][prime]['L']
        self.A = parameters['csidh']['A']
        print("// Shortest Differential Addition Chains (SDAC) for each l_i;")
        # List of Small odd primes, L := [l_0, ..., l_{n-1}]
        print("// SDAC's to be read from a file")
        path = sdacs_data + prime
        self.SDACS = self.filename_to_list_of_lists_of_ints(path)
        if len(self.SDACS) == 0:
            print("// SDAC's to be computed")
            self.SDACS = generate_sdacs(L)
            print("// Storing SDAC's in a file")
            write_list_of_lists_of_ints_to_file(path, self.SDACS)
        self.SDACS_LENGTH = list(map(len, self.SDACS))

        cMUL = lambda l: numpy.array(
            [
                4.0 * (self.SDACS_LENGTH[L.index(l)] + 2),
                2.0 * (self.SDACS_LENGTH[L.index(l)] + 2),
                6.0 * (self.SDACS_LENGTH[L.index(l)] + 2) - 2.0,
            ]
        )
        self.C_xMUL = list(map(cMUL, L))  # list of the costs of each [l]P

    def elligator(self, A):

        Ap = self.fp.fp_add(A[0], A[0])
        Ap = self.fp.fp_sub(Ap, A[1])
        Ap = self.fp.fp_add(Ap, Ap)
        Cp = A[1]

        u = random.randint(2, self.fp.p_minus_one_halves)
        u_squared = self.fp.fp_sqr(u)

        u_squared_plus_one = self.fp.fp_add(u_squared, 1)
        u_squared_minus_one = self.fp.fp_sub(u_squared, 1)

        C_times_u_squared_minus_one = self.fp.fp_mul(Cp, u_squared_minus_one)
        AC_times_u_squared_minus_one = self.fp.fp_mul(Ap, C_times_u_squared_minus_one)

        tmp = self.fp.fp_sqr(Ap)
        tmp = self.fp.fp_mul(tmp, u_squared)
        aux = self.fp.fp_sqr(C_times_u_squared_minus_one)
        tmp = self.fp.fp_add(tmp, aux)
        tmp = self.fp.fp_mul(AC_times_u_squared_minus_one, tmp)

        alpha, beta = 0, u
        alpha, beta = self.fp.fp_cswap(alpha, beta, tmp == 0)
        u_squared_plus_one = self.fp.fp_mul(alpha, u_squared_plus_one)
        alpha = self.fp.fp_mul(alpha, C_times_u_squared_minus_one)

        Tp_X = self.fp.fp_add(Ap, alpha)
        Tm_X = self.fp.fp_mul(Ap, u_squared)
        Tm_X = self.fp.fp_add(Tm_X, alpha)
        Tm_X = self.fp.fp_sub(0, Tm_X)

        tmp = self.fp.fp_add(tmp, u_squared_plus_one)
        Tp_X, Tm_X = self.fp.fp_cswap(Tp_X, Tm_X, (1 - jacobi(tmp, self.fp.p)) // 2)

        return (
            [Tp_X, C_times_u_squared_minus_one],
            [Tm_X, C_times_u_squared_minus_one],
        )



    def filename_to_list_of_lists_of_ints(self, path):
        res = []
        try:
            with open(path) as fh:
                for line in fh:
                    res.append(list(map(int, line.split())))
        except:
            res = []
        return res

    def write_list_of_lists_of_ints_to_file(self, path, data):
        with open(path, 'w') as fh:
            for line in data:
                fh.writelines(' '.join(str(v) for v in v in line))
            fh.writelines()

    def generate_sdacs(self, L):
        return list(
            map(sdac, L)
        )  # Shortest Differential Addition Chains for each small odd prime l in L


    def measure(self, x):
        """
        SQR = 1.00
        # In F_p, we have SQR_{F_p} = SQR x MUL_{F_p}
        ADD = 0.00
        # In F_p, we have ADD_{F_p} = ADD x MUL_{F_p}
        """
        SQR = 1.00
        ADD = 0.00
        return (x[0] + SQR * x[1] + ADD * x[2])

    def dacs(self, l, r0, r1, r2, chain):
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

    def sdac(self, l):
        '''
        sdac()
        input: a small odd prime number l
        output: the shortest differential additions chains corresponding with the input l

        NOTE: this function uses a recursive function
        '''
        all_dacs = dacs(l, 1, 2, 3, [])
        return min(all_dacs, key=lambda t: len(t[0]))[0]

    def affine_to_projective(self, affine):
        """
        input : the affine Montgomery coefficient A/C
        output: projective Montgomery constants A24 := A + 2C and C24 := 4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        """
        res = [ 0, 0 ]
        res[0] = self.fp.fp_add(affine, self.A[0])
        res[1] = self.fp.fp_add(0, self.A[1])
        return res

    def coeff(self, A):
        '''
        ----------------------------------------------------------------------
        coeff()
        input : projective Montgomery constants A24 := A + 2C and C24 := 4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the affine Montgomery coefficient A/C
        ----------------------------------------------------------------------
        '''
        output = self.fp.fp_add(A[0], A[0])  # (2 * A24)
        output = self.fp.fp_sub(output, A[1])  # (2 * A24) - C24
        C24_inv = self.fp.fp_inv(A[1])  # 1 / (C24)
        output = self.fp.fp_add(output, output)  # 4*A = 2[(2 * A24) - C24]
        output = self.fp.fp_mul(output, C24_inv)  # A/C = 2[(2 * A24) - C24] / C24

        return output

    def isinfinity(self, P):
        """
        isinfinity(P) determines if x(P) := (XP : ZP) = (1 : 0)
        """
        return P[1] == 0


    def areequal(self, P, Q):
        """ areequal(P, Q) determines if x(P) = x(Q) """
        return self.fp.fp_mul(P[0], Q[1]) == self.fp.fp_mul(P[1], Q[0])

    def xDBL(self, P, A):
        '''
        ----------------------------------------------------------------------
        xDBL()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, and
                the  projective Montgomery constants A24:= A + 2C and C24:=4C
                where E : y^2 = x^3 + (A/C)*x^2 + x
        output: the projective Montgomery x-coordinate point x([2]P)
        ----------------------------------------------------------------------
        '''
        t_0 = self.fp.fp_sub(P[0], P[1])
        t_1 = self.fp.fp_add(P[0], P[1])
        t_0 = self.fp.fp_sqr(t_0)
        t_1 = self.fp.fp_sqr(t_1)
        Z = self.fp.fp_mul(A[1], t_0)
        X = self.fp.fp_mul(Z, t_1)
        t_1 = self.fp.fp_sub(t_1, t_0)
        t_0 = self.fp.fp_mul(A[0], t_1)
        Z = self.fp.fp_add(Z, t_0)
        Z = self.fp.fp_mul(Z, t_1)

        return [X, Z]




    def xADD(self, P, Q, PQ):
        '''
        ----------------------------------------------------------------------
        xADD()
        input : the projective Montgomery x-coordinate points x(P) := XP/ZP,
                x(Q) := XQ/ZQ, and x(P-Q) := XPQ/ZPQ
        output: the projective Montgomery x-coordinate point x(P+Q)
        ----------------------------------------------------------------------
        '''
        a = self.fp.fp_add(P[0], P[1])
        b = self.fp.fp_sub(P[0], P[1])
        c = self.fp.fp_add(Q[0], Q[1])
        d = self.fp.fp_sub(Q[0], Q[1])
        a = self.fp.fp_mul(a, d)
        b = self.fp.fp_mul(b, c)
        c = self.fp.fp_add(a, b)
        d = self.fp.fp_sub(a, b)
        c = self.fp.fp_sqr(c)
        d = self.fp.fp_sqr(d)
        X = self.fp.fp_mul(PQ[1], c)
        Z = self.fp.fp_mul(PQ[0], d)
        return [X, Z]


    # Modificar esta parte para usar cadenas de addicion
    def xMUL(self, P, A, j):
        '''
        ----------------------------------------------------------------------
        xMUL()
        input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
                projective Montgomery constants A24:= A + 2C and C24:=4C where
                E : y^2 = x^3 + (A/C)*x^2 + x, and an positive integer j
        output: the projective Montgomery x-coordinate point x([L[j]]P)
        ----------------------------------------------------------------------
        '''
        P2 = self.xDBL(P, A)
        R = [P, P2, self.xADD(P2, P, P)]

        for i in range(self.SDACS_LENGTH[j] - 1, -1, -1):

            if self.isinfinity(R[self.SDACS[j][i]]):
                T = self.xDBL(R[2], A)
            else:
                T = self.xADD(R[2], R[self.SDACS[j][i] ^ 1], R[self.SDACS[j][i]])

            R[0] = list(R[self.SDACS[j][i] ^ 1])
            R[1] = list(R[2])
            R[2] = list(T)

        return R[2]




    def prime_factors(self, P, A, points):
        '''
        ----------------------------------------------------------------------
        prime_factors()
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
                    second_P = self.xMUL(second_P, A, points[j])
                    first_half.append(points[j])

                # 2nd half
                second_half = []
                first_P = P
                for j in range(h, n):
                    first_P = self.xMUL(first_P, A, points[j])
                    second_half.append(points[j])

                return self.prime_factors(first_P, A, first_half) + self.prime_factors(
                    second_P, A, second_half
                )

            return []


    def isfull_order(self, seq):
        tmp = [not self.isinfinity(seq_i) for seq_i in seq]
        return reduce(lambda x, y: (x and y), tmp)


    def full_torsion_points(self, A):

        if self.style != 'wd1':
            output = [[0, 0], [0, 0]]
        else:
            output = [[0, 0], [1, 1]]

        while [0, 0] in output:

            T_p, T_m = self.elligator(A)
            for i in range(0, self.fp.exponent_of_two, 1):
                T_p = self.xDBL(T_p, A)

            if self.isfull_order(self.prime_factors(T_p, A, range(0, self.fp.n, 1))) and output[
                0
            ] == [0, 0]:
                output[0] = list(T_p)

            if self.style != 'wd1':
                for i in range(0, self.fp.exponent_of_two, 1):
                    T_m = self.xDBL(T_m, A)
                if self.isfull_order(self.prime_factors(T_m, A, range(0, self.fp.n, 1))) and output[
                    1
                ] == [0, 0]:
                    output[1] = list(T_m)

        return output[0], output[1]


    def CrissCross(self, alpha, beta, gamma, delta):

        t_1 = self.fp.fp_mul(alpha, delta)
        t_2 = self.fp.fp_mul(beta, gamma)
        return self.fp.fp_add(t_1, t_2), self.fp.fp_sub(t_1, t_2)


    def validate(self, A):

        while True:

            T_p, _ = self.elligator(A)
            for i in range(0, self.fp.exponent_of_two, 1):
                T_p = self.xDBL(T_p, A)

            P = self.prime_factors(T_p, A, range(0, self.fp.n, 1))

            bits_of_the_order = 0
            for i in range(0, self.fp.n, 1):

                if self.isinfinity(P[i]) == False:

                    Q = self.xMUL(P[i], A, i)

                    if self.isinfinity(Q) == False:
                        return False

                    bits_of_the_order += bitlength(self.L[i])
                    if bits_of_the_order > self.fp.validation_stop:
                        return True


