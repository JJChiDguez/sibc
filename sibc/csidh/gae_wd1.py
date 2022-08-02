from random import SystemRandom
import numpy
from math import floor, ceil, sqrt, log

from sibc.math import isequal, bitlength, hamming_weight, cswap, sign
from sibc.constants import parameters
from sibc.common import geometric_serie, rounds

class CachedIndexTuple(tuple):
    def __init__(self, container):
        self.cached_index = dict((b,a) for a,b in enumerate(container))
        self.index = self.cached_index.__getitem__

class Gae_wd1(object):
    def __init__(self, prime, tuned, curve, formula):
        self.random = SystemRandom()

        # In order to achieve efficiency, the optimal strategies and their cost are saved in two global dictionaries (hash tables)
        self.S = {1: {}}  # Initialization of each strategy
        self.C = {1: {}}  # Initialization of the costs: 0.

        self.curve = curve
        self.prime = prime
        self.formula = formula
        self.formula.L = CachedIndexTuple(self.formula.L)
        self.formula_name = formula.name
        self.field = self.curve.field
        self.tuned = tuned
        self.L = CachedIndexTuple(parameters['csidh'][prime]['L'])
        self.m = parameters['csidh'][prime]['wd1']['m']
        self.c_xmul = self.curve.c_xmul
        n = parameters['csidh'][prime]['n']
        for i in range(n):
            self.S[1][tuple([self.L[i]])] = []
            # Strategy with a list with only one element (a small odd prime number l_i)
            self.C[1][tuple([self.L[i]])] = (
                self.formula.c_xisog[i] - self.curve.c_xmul[i] + 2.0 * numpy.array([4.0, 2.0, 6.0])
                )
            # For catching the weigth of horizontal edges of the form [(0,j),(0,j+1)]
        for i in range(2, n + 1):
            self.C[i] = {}
            self.S[i] = {}
        ######################################################################################################################
        # Next functions are used for computing optimal bounds
        self.basis = numpy.eye(n, dtype=int)

        self.C_out, self.L_out, self.R_out, self.S_out, self.r_out = self.strategy_block_cost(
            self.L[::-1], self.m[::-1]
        )
        self.temporal_m = list(set(self.m))

    def GAE_at_0(self, exp):
        if (len(self.temporal_m) == 1) or (
            (len(self.temporal_m) == 2) and (0 in self.temporal_m)
        ):
            return self.GAE(
                [self.curve.field(2), self.curve.field(4)],
                exp,
                [self.L_out[0]],
                [self.R_out[0]],
                [self.S_out[0]],
                [self.temporal_m[-1]],
                self.m,
            )
        else:
            return self.GAE(
                [self.curve.field(2), self.curve.field(4)],
                exp, self.L_out, self.R_out, self.S_out, self.r_out, self.m
            )

    def GAE_at_A(self, exp, A):
        assert self.curve.issupersingular(A), "non-supersingular input curve"
        if (len(self.temporal_m) == 1) or (
            (len(self.temporal_m) == 2) and (0 in self.temporal_m)
        ):
            ss = self.GAE(
                A,
                exp,
                [self.L_out[0]],
                [self.R_out[0]],
                [self.S_out[0]],
                [self.temporal_m[-1]],
                self.m,
            )
        else:
            ss = self.GAE(A, exp, self.L_out, self.R_out, self.S_out, self.r_out, self.m)
        return ss

    def random_exponents(self, m=None):
        if m is None:
            m = self.m
        """
        random_exponents(m) implements an uniform random sample from S(m_1) x S(m_2) x ... x S(m_n)
        """
        return [
            self.random.randint(0, m_i)
            for m_i in self.m
        ]

    def print_exponents(self, label : str, exp):
        print("%s := ( %s );" % (label, ' '.join(map(str, exp))) )
        return None

    def security(self, M, n):
        """
        security()
        inputs : the list M of maximum number of degree-(l_i) isogeny constructions
                 to be performed, and the length of M
        output : bits of security that M brings
        """
        return sum(list([log(M[i] + 1, 2) for i in range(n)]))

    def dynamic_programming_algorithm(self, L, n):
        """
        dynamic_programming_algorithm():
        inputs: the list of small odd primes to be processed and its length
        output: the optimal strategy and its cost of the input list of small odd primes
        """
        # If the approach uses dummy operations, to set DUMMY = 2.0;
        # otherwise, to set DUMMY = 1.0 (dummy free approach);

        if len(L) != n:

            # If the list of prime numbers doesn't have size n, then we return [],-1
            print(
                "error:\tthe list of prime numbers has different size from %d."
                % n
            )
            return [], -1
        else:

            # Assuming #L = n, we proceed.
            get_neighboring_sets = lambda L, k: [
                tuple(L[i : i + k]) for i in range(n - k + 1)
            ]  # This function computes all the k-tuple: (l_1, l_2, ..., l_{k)),
            # (l_2, l_3, ..., l_{k+1)), ..., (l_{n-k}, l_{n-k+1, ..., l_{n)).
            for i in range(2, n + 1):

                # TODOXXX see the loop unrolling in gae_df.py
                for Tuple in get_neighboring_sets(L, i):

                    if self.C[i].get(Tuple) is None:

                        alpha = [
                            (
                                b,
                                self.C[len(Tuple[:b])][Tuple[:b]]
                                + self.C[  # Subtriangle on the right side with b leaves
                                    len(Tuple[b:])
                                ][
                                    Tuple[b:]
                                ]
                                + 1.0  # Subtriangle on the left side with (i - b) leaves
                                * sum(
                                    [
                                        self.c_xmul[
                                            self.formula.L.index(t)
                                        ]
                                        for t in Tuple[:b]
                                    ]
                                )
                                + 1.0  # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with b leaves
                                * sum(
                                    [
                                        self.formula.c_xeval[
                                            self.formula.L.index(t)
                                        ]
                                        for t in Tuple[b:]
                                    ]
                                )
                                + 1.0  # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                                * sum(
                                    [
                                        self.c_xmul[
                                            self.formula.L.index(t)
                                        ]
                                        for t in Tuple[b:]
                                    ]
                                ),  # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                            )
                            for b in range(1, i)
                        ]
                        b, self.C[i][Tuple] = min(
                            alpha, key=lambda t: self.curve.measure(t[1])
                        )  # We save the minimal cost corresponding to the triangle with leaves Tuple
                        self.S[i][Tuple] = (
                            [b]
                            + self.S[i - b][Tuple[b:]]
                            + self.S[b][Tuple[:b]]
                        )  # We save the optimal strategy corresponding to the triangle with leaves Tuple

            return (
                self.S[n][tuple(L)],
                self.C[n][tuple(L)]
                + self.curve.c_xmul[self.formula.L.index(L[0])]
                - 2.0 * numpy.array([4.0, 2.0, 6.0]),
            )  # The weight of the horizontal edges [(0,n-1),(0,n)] must be equal to c_xisog[self.formula.L.index(L[0])].

    def evaluate_strategy(self, E, P, L, strategy, n, m, e):
        """
        evaluate_strategy():
        inputs : a projective Montgomery constants A24:= A + 2C and C24:=4C where E : y^2 = x^3 + (A/C)*x^2 + x,
                 a list of two projective Montgomery x-coordinate torsion-(l_1 x ... l_n) points x(T_+) := XP/ZP
                 and x(T_-) := XM/ZM, the list of small odd primes [l_1, l_2, ..., l_n], an strategy and length
                 of the given list of small odd primes, maximum number of degree-l_i isogeny constructions, and
                 the secret integer vector to be evaluated
        output : the projective Montgomery constants a24:= a + 2c and c24:=4c where E': y^2 = x^3 + (a/c)*x^2 + x
                 is E / <P>, the new maximum and current number of degree-l_i isogeny constructions to be performed
                 after the strategy evaluation.

        NOTE: T_+ belongs to E[pi - 1] and T_- belongs to E[pi + 1]. In particular, P = [T_-, T_+]
        """
        v = list(m)
        u = list(e)
        ramifications = []
        moves = [
            0
        ]  # moves: this list determines whether an isogeny construction must be performed
        k = 0  # k: current element of the strategy

        ramifications.append(P)
        E_i = list(E)
        for i in range(len(strategy)):

            pos = self.formula.L.index(
                L[n - 1 - i]
            )  # Current element of self.formula.L to be required

            # Reaching the vertex (n - 1 - i, i)

            # Vertical edges
            prev = sum(moves)
            while prev < (n - 1 - i):

                moves.append(
                    strategy[k]
                )  # Number of vertical edges to be performed
                T = list(ramifications[-1])  # New ramification
                for j in range(prev, prev + strategy[k], 1):
                    T = self.curve.xmul(T, E_i, self.formula.L.index(L[j]))

                ramifications.append(T)
                prev += strategy[k]
                k += 1

            # At this point, vertex (n - 1 - i, i) has been reached
            if (
                v[pos] > 0
            ):  # Maximum number of degree-l_{pos} isogeny constructions?

                # At this step, ramifications[-1] is the i-th leaf
                if self.curve.isinfinity(ramifications[-1]) == False:

                    # Dummy or NOT Dummy degree-(l_{n-1-i}) isogeny construction, that's the question?
                    b_i = isequal[u[pos] == 0]

                    ramifications[-1][0], ramifications[0][0] = cswap(
                        ramifications[-1][0], ramifications[0][0], b_i
                    )
                    ramifications[-1][1], ramifications[0][1] = cswap(
                        ramifications[-1][1], ramifications[0][1], b_i
                    )

                    if self.formula_name != 'tvelu':
                        # This branchs corresponds with the use of the new velu's formulaes

                        if self.tuned:
                            self.formula.set_parameters_velu(
                                self.formula.sJ_list[pos],
                                self.formula.sI_list[pos],
                                pos
                            )

                        else:
                            # -------------------------------------------------------------
                            # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                            # These paramters are required in KPs, xISOG, and xEVAL
                            if self.formula.L[pos] == 3:
                                b = 0
                                c = 0
                            else:
                                b = int(floor(sqrt(self.formula.L[pos] - 1) / 2.0))
                                c = int(
                                    floor((self.formula.L[pos] - 1.0) / (4.0 * b))
                                )

                            self.formula.set_parameters_velu(b, c, pos)

                        if (
                            self.formula_name == 'hvelu'
                            and self.formula.L[pos] <= self.formula.HYBRID_BOUND
                        ):
                            K = self.formula.kps(ramifications[-1], E_i, pos)
                        else:
                            self.formula.kps(ramifications[-1], E_i, pos)
                            T = self.curve.xmul(ramifications[-1], E_i, pos)

                    else:
                        K = self.formula.kps(ramifications[-1], E_i, pos)

                    ramifications[-1][0], ramifications[0][0] = cswap(
                        ramifications[-1][0], ramifications[0][0], b_i
                    )
                    ramifications[-1][1], ramifications[0][1] = cswap(
                        ramifications[-1][1], ramifications[0][1], b_i
                    )

                    # New isogeny construction
                    C_i = self.formula.xisog(E_i, pos)

                    # Next, the horizontal edge [(0,i),(0,i+1)] is performed
                    if self.formula_name == 'tvelu' or (
                        self.formula_name == 'hvelu'
                        and self.formula.L[pos] <= self.formula.HYBRID_BOUND
                    ):
                        d_i = (self.formula.L[pos] - 1) // 2
                        mask = isequal[
                            (self.formula.L[pos] == 3)
                        ]  # catching special case when l = 3

                        Z = self.formula.yadd(
                            K[(d_i + mask) - 1], K[0], K[(d_i + mask) - 2]
                        )  # y([d_i + 1]K[0])
                        Z[0], K[d_i][0] = cswap(
                            Z[0], K[d_i][0], mask ^ 1
                        )
                        Z[1], K[d_i][1] = cswap(
                            Z[1], K[d_i][1], mask ^ 1
                        )

                        T = self.formula.yadd(
                            K[d_i], K[d_i - 1], K[0]
                        )  # y([2*d_i + 1]K[0]) := y([l_i]K[0])
                        T = [
                            (T[1] + T[0]),
                            (T[1] - T[0]),
                        ]  # x([l_i]K[0])

                    if self.formula_name == 'tvelu' or (
                        self.formula_name == 'hvelu'
                        and self.formula.L[pos] <= self.formula.HYBRID_BOUND
                    ):
                        ramifications[0] = self.formula.xeval(ramifications[0], pos)
                    else:
                        ramifications[0] = self.formula.xeval(ramifications[0], E_i)

                    T[0], ramifications[0][0] = cswap(
                        T[0], ramifications[0][0], b_i
                    )
                    T[1], ramifications[0][1] = cswap(
                        T[1], ramifications[0][1], b_i
                    )

                    # The remainder horizontal edges are performed
                    for j in range(1, len(moves) - 1, 1):

                        T = self.curve.xmul(ramifications[j], E_i, pos)

                        if self.formula_name == 'tvelu' or (
                            self.formula_name == 'hvelu'
                            and self.formula.L[pos] <= self.formula.HYBRID_BOUND
                        ):
                            ramifications[j] = self.formula.xeval(
                                ramifications[j], pos
                            )
                        else:
                            ramifications[j] = self.formula.xeval(
                                ramifications[j], E_i
                            )

                        T[0], ramifications[j][0] = cswap(
                            T[0], ramifications[j][0], b_i
                        )
                        T[1], ramifications[j][1] = cswap(
                            T[1], ramifications[j][1], b_i
                        )

                    C_i[0], E_i[0] = cswap(C_i[0], E_i[0], b_i ^ 1)
                    C_i[1], E_i[1] = cswap(C_i[1], E_i[1], b_i ^ 1)

                    v[pos] -= 1
                    u[pos] -= b_i ^ 1
            else:
                for j in range(0, len(moves) - 1, 1):

                    ramifications[j] = self.curve.xmul(ramifications[j], E_i, pos)

            moves.pop()
            ramifications.pop()

        pos = self.formula.L.index(
            L[0]
        )  # Current element of self.formula.L to be required
        if self.curve.isinfinity(ramifications[0]) == False:

            if m[pos] > 0:

                b_i = isequal[e[pos] == 0]

                if self.formula_name != 'tvelu':
                    # This branchs corresponds with the use of the new velu's formulaes

                    if self.tuned:
                        self.formula.set_parameters_velu(
                            self.formula.sJ_list[pos],
                            self.formula.sI_list[pos],
                            pos
                        )

                    else:
                        # -------------------------------------------------------------
                        # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                        # These paramters are required in KPs, xISOG, and xEVAL
                        if self.formula.L[pos] == 3:
                            b = 0
                            c = 0
                        else:
                            b = int(floor(sqrt(self.formula.L[pos] - 1) / 2.0))
                            c = int(floor((self.formula.L[pos] - 1.0) / (4.0 * b)))

                        self.formula.set_parameters_velu(b, c, pos)

                    self.formula.kps(ramifications[0], E_i, pos)

                else:
                    self.formula.kps(ramifications[0], E_i, pos)

                C_i = self.formula.xisog(E_i, pos)
                C_i[0], E_i[0] = cswap(C_i[0], E_i[0], b_i ^ 1)
                C_i[1], E_i[1] = cswap(C_i[1], E_i[1], b_i ^ 1)

                v[pos] -= 1
                u[pos] -= b_i ^ 1

        return E_i, v, u

    def GAE(self, A, e, L, R, St, r, m):
        """
        GAE():
        inputs : a projective Montgomery constants A24:= A + 2C and C24:=4C where E : y^2 = x^3 + (A/C)*x^2 + x,
                 the secret integer vector to be evaluated, all the sublist of small odd primes to be required and
                 their complement of each sublist, the (optimal) strategies corresponding to each sublist and their
                 number of times to be evaluated, and the maximum number of degree-l_i isogeny constructions
        output : the projective Montgomery constants a24:= a + 2c and c24:=4c where E': y^2 = x^3 + (a/c)*x^2 + x
                 corresponds with the image of the group action evaluation.

        NOTE: GAE comes from Group Action Evaluation, and the input sublists are determined by the output of function
              rounds().
              THIS IS THE IMPLEMENTATION OF OUR PROPOSED STRATEGY METHOD.
        """
        E_k = list(A)
        n = len(L)

        for j in range(0, n, 1):

            for k in range(0, r[j], 1):

                T_p, T_m = self.curve.elligator(E_k)
                T_p = self.curve.prac(self.curve.cofactor, T_p, E_k)

                for l in R[j]:
                    T_p = self.curve.xmul(
                        T_p, E_k, self.formula.L.index(l)
                    )

                E_k, m, e = self.evaluate_strategy(
                    E_k,
                    T_p,
                    L[j],
                    St[j],
                    len(L[j]),
                    m,
                    e,
                )

        # Multiplicative strategy on the set of unreached small odd prime numbers
        unreached_sop = [
            self.formula.L[i]
            for i in range(len(self.formula.L))
            if m[i] > 0
        ]
        remainder_sop = [
            l for l in self.formula.L if l not in unreached_sop
        ]

        while len(unreached_sop) > 0:

            T_p, T_m = self.curve.elligator(E_k)
            T_p = self.curve.prac(self.curve.cofactor, T_p, E_k)

            for l in remainder_sop:
                T_p = self.curve.xmul(T_p, E_k, self.formula.L.index(l))

            current_n = len(unreached_sop)
            E_k, m, e = self.evaluate_strategy(
                E_k,
                T_p,
                unreached_sop,
                list(range(current_n - 1, 0, -1)),
                current_n,
                m,
                e,
            )

            # If the maximum of degree-(l_k) has been reached then the current batch (and its complement) must be updated
            tmp_unreached = [
                unreached_sop[k]
                for k in range(current_n)
                if m[self.formula.L.index(unreached_sop[k])] > 0
            ]
            tmp_remainder = [
                unreached_sop[k]
                for k in range(current_n)
                if m[self.formula.L.index(unreached_sop[k])] == 0
            ]

            unreached_sop = list(
                tmp_unreached
            )  # Removing elements from the batch
            remainder_sop = (
                remainder_sop + tmp_remainder
            )  # Adding elements to the complement of the batch

        return E_k

    def strategy_block_cost(self, L, e):
        """
        Next function computes the expected cost of our approach by assuming we have full torsion points
        """
        elligator_cost = numpy.array([7.0, 3.0, 10.0])  # Elligator cost
        mul_fp_by_cofactor = (
            numpy.array([4.0, 2.0, 4.0]) * int(ceil(1.64 * bitlength(self.curve.cofactor)))
        )  # Cost of computing x([self.curve.cofactor]P)

        n = len(L)
        e_prime = [geometric_serie(e[k], L[k]) for k in range(n)]

        tmp_r, tmp_Ls, tmp_Cs = rounds(e_prime, n)

        C_e = numpy.array([0.0, 0.0, 0.0])
        S_out = []
        L_out = []
        R_out = []
        for j in range(len(tmp_r)):

            R_out.append([L[k] for k in tmp_Cs[j]])
            L_out.append([L[k] for k in tmp_Ls[j]])

            bo_C = 1.0 * sum(
                [
                    self.c_xmul[self.formula.L.index(L[k])]
                    for k in tmp_Cs[j]
                ]
            )
            S_tmp, go_C = self.dynamic_programming_algorithm(
                [L[k] for k in tmp_Ls[j]], len(tmp_Ls[j])
            )

            S_out.append(S_tmp)
            C_e = (
                C_e
                + (go_C + bo_C + elligator_cost + 1.0 * mul_fp_by_cofactor)
                * tmp_r[j]
            )

        return C_e, L_out, R_out, S_out, tmp_r
