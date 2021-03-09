import random
import numpy
from sympy import symbols, floor, sqrt, sign

from sidh._math import isequal, bitlength, hamming_weight
from sidh.constants import parameters

# In order to achieve efficiency, the optimal strategies and their cost are saved in two global dictionaries (hash tables)
S = {1: {}}  # Initialization of each strategy
C = {1: {}}  # Initialization of the costs: 0.
L = None
basis = None
C_xMUL = None

class Gae_df(object):
    def __init__(self, prime, verbose, curve, formula):
        global S
        global C
        global L
        global basis
        global C_xMUL
        self.curve = curve
        self.prime = prime
        self.formula = formula
        self.formula_name = formula.name
        self.fp = self.curve.fp
        self.verbose = verbose
        self.L = L = parameters['csidh'][prime]['L']
        self.m = parameters['csidh'][prime]['df']['m']
        C_xMUL = self.curve.C_xMUL
        n = parameters['csidh'][prime]['n']
        for i in range(n):
            S[1][tuple([L[i]])] = []
            # Strategy with a list with only one element (a small odd prime number l_i)
            C[1][tuple([L[i]])] = self.formula.C_xISOG[i]
            # For catching the weigth of horizontal edges of the form [(0,j),(0,j+1)]
        for i in range(2, n + 1):
            C[i] = {}
            S[i] = {}
        ######################################################################################################################
        # Next functions are used for computing optimal bounds
        #self.basis = basis = numpy.eye(n, dtype=int)

    def pubkey(self, sk):
        pass

    def dh(self, sk, pk):
        assert self.curve.validate(pk), "public key does not validate"
        temporal_m = list(set(self.m))
        C_out, L_out, R_out, S_out, r_out = self.strategy_block_cost(L[::-1], self.m[::-1])
        ss = self.GAE(
            pk,
            sk,
            [L_out[0]],
            [R_out[0]],
            [S_out[0]],
            [temporal_m[-1]],
            self.m,
        )
        return self.curve.coeff(ss)

    def random_key(self):
        """
        random_key(m) implements an uniform random sample from S(m_1) x S(m_2) x ... x S(m_n)
        """
        return [2 * (random.randint(0, m_i) - (m_i // 2)) - (m_i % 2) for m_i in self.m ]

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
        global S, C
        # If the approach uses dummy operations, to set DUMMY = 2.0;
        # otherwise, to set DUMMY = 1.0 (dummy free approach);

        if len(L) != n:

            # If the list of prime numbers doesn't have size n, then we return [],-1
            print(
                "error:\tthe list of prime numbers has different size from %d." % n
            )
            return [], -1
        else:

            # Assuming #L = n, we proceed.
            get_neighboring_sets = lambda L, k: [
                tuple(L[i : i + k]) for i in range(n - k + 1)
            ]  # This function computes all the k-tuple: (l_1, l_2, ..., l_{k)),
            # (l_2, l_3, ..., l_{k+1)), ..., (l_{n-k}, l_{n-k+1, ..., l_{n)).
            for i in range(2, n + 1):

                for Tuple in get_neighboring_sets(L, i):

                    if C[i].get(Tuple) is None:

                        alpha = [
                            (
                                b,
                                C[len(Tuple[:b])][Tuple[:b]]
                                + C[  # Subtriangle on the right side with b leaves
                                    len(Tuple[b:])
                                ][Tuple[b:]]
                                + 2.0  # Subtriangle on the left side with (i - b) leaves
                                * sum(
                                    [C_xMUL[self.formula.global_L.index(t)] for t in Tuple[:b]]
                                )
                                + 2.0  # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with b leaves
                                * sum(
                                    [self.formula.C_xEVAL[self.formula.global_L.index(t)] for t in Tuple[b:]]
                                )
                                + 1.0  # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                                * sum(
                                    [C_xMUL[self.formula.global_L.index(t)] for t in Tuple[b:]]
                                ),  # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                            )
                            for b in range(1, i - 1)
                        ] + [
                            (
                                i - 1,
                                C[i - 1][Tuple[: (i - 1)]]
                                + C[  # Subtriangle on the right side with (i - 1) leaves
                                    1
                                ][
                                    Tuple[(i - 1) :]
                                ]
                                + 1.0  # Subtriangle on the left side with 1 leaf (only one vertex)
                                * sum(
                                    [
                                        C_xMUL[self.formula.global_L.index(t)]
                                        for t in Tuple[: (i - 1)]
                                    ]
                                )
                                + 2.0  # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with 1 leaf
                                * self.formula.C_xEVAL[self.formula.global_L.index(Tuple[i - 1])]
                                + 1.0  # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - 1) leaves
                                * C_xMUL[
                                    self.formula.global_L.index(Tuple[i - 1])
                                ],  # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - 1) leaves
                            )
                        ]
                        b, C[i][Tuple] = min(
                            alpha, key=lambda t: self.curve.measure(t[1])
                        )  # We save the minimal cost corresponding to the triangle with leaves Tuple
                        S[i][Tuple] = (
                            [b] + S[i - b][Tuple[b:]] + S[b][Tuple[:b]]
                        )  # We save the optimal strategy corresponding to the triangle with leaves Tuple

            return (
                S[n][tuple(L)],
                C[n][tuple(L)],
            )  # The weight of the horizontal edges [(0,n-1),(0,n)] must be equal to C_xISOG[self.formula.global_L.index(L[0])].



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
        t = 1  # sign to be flip

        ramifications.append(
            list(P)
        )  # list of pair of points (T_-, T_+) such that T_- in E[\pi + 1] and T_+ in E[\pi - 1]
        E_i = list(E)
        for i in range(len(strategy)):

            pos = self.formula.global_L.index(
                L[n - 1 - i]
            )  # Current element of self.formula.global_L to be required

            # Reaching the vertex (n - 1 - i, i)

            # Vertical edges
            prev = sum(moves)
            while (prev + strategy[k]) < (n - 1 - i):

                moves.append(
                    strategy[k]
                )  # Number of vertical edges to be performed
                T = list([list(ramifications[-1][0]), list(ramifications[-1][1])])
                for j in range(prev, prev + strategy[k], 1):
                    T = list(
                        [
                             self.curve.xMUL(T[0], E_i, self.formula.global_L.index(L[j])),
                             self.curve.xMUL(T[1], E_i, self.formula.global_L.index(L[j])),
                        ]
                    )

                ramifications.append(list([list(T[0]), list(T[1])]))
                prev += strategy[k]
                k += 1

            # Vertical edges without ramifications
            s_i = sign(e[pos])  # Sign of e[pos]
            c_i = (s_i + 1) // 2  # Constant-swap of T_+ and T_-

            if prev < (n - 1 - i):

                moves.append(
                    strategy[k]
                )  # Number of vertical edges (without ramifications) to be performed
                T = list([list(ramifications[-1][0]), list(ramifications[-1][1])])

                T[0][0], T[1][0] = self.fp.fp_cswap(T[0][0], T[1][0], c_i)
                T[0][1], T[1][1] = self.fp.fp_cswap(T[0][1], T[1][1], c_i)
                for j in range(prev, prev + strategy[k], 1):
                    T[0] = list(
                         self.curve.xMUL(T[0], E_i, self.formula.global_L.index(L[j]))
                    )  # A single scalar multiplication is required

                T[0][0], T[1][0] = self.fp.fp_cswap(T[0][0], T[1][0], c_i)
                T[0][1], T[1][1] = self.fp.fp_cswap(T[0][1], T[1][1], c_i)

                ramifications.append(list([list(T[0]), list(T[1])]))
                prev += strategy[k]
                k += 1

            # At this point, vertex (n - 1 - i, i) has been reached
            if (
                v[pos] > 0
            ):  # Maximum number of degree-l_{pos} isogeny constructions?

                # At this step, ramifications[-1] is the i-th leaf
                # Swap the torsion-(l_{n-j}) elliptic curve points on the current leaf
                #                           T_- <- ramifications[-1][0]
                #                           T_+ <- ramifications[-1][1]
                ramifications[-1][0][0], ramifications[-1][1][0] = self.fp.fp_cswap(
                    ramifications[-1][0][0], ramifications[-1][1][0], c_i
                )  # XT_+ <-> XT_-
                ramifications[-1][0][1], ramifications[-1][1][1] = self.fp.fp_cswap(
                    ramifications[-1][0][1], ramifications[-1][1][1], c_i
                )  # ZT_+ <-> ZT_-

                if  self.curve.isinfinity(ramifications[-1][0]) == False:

                    if self.formula_name != 'tvelu':
                        # This branchs corresponds with the use of the new velu's formulaes

                        E_prev = list(E_i)
                        if self.verbose:
                            self.formula.set_parameters_velu(sJ_list[pos], sI_list[pos], pos)

                        else:
                            # -------------------------------------------------------------
                            # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                            # These paramters are required in KPs, xISOG, and xEVAL
                            if self.formula.global_L[pos] == 3:
                                b = 0
                                c = 0
                            else:
                                b = int(floor(sqrt(self.formula.global_L[pos] - 1) / 2.0))
                                c = int(floor((self.formula.global_L[pos] - 1.0) / (4.0 * b)))

                            self.formula.set_parameters_velu(b, c, pos)

                    self.formula.KPs(ramifications[-1][0], E_i, pos)

                    # New isogeny construction
                    E_i = self.formula.xISOG(E_i, pos)

                    # Horizontal edges are performed
                    for j in range(0, len(moves) - 1, 1):

                        # Swap points in E[\pi - 1] and E[\pi + 1]
                        ramifications[j][0][0], ramifications[j][1][0] = self.fp.fp_cswap(
                            ramifications[j][0][0], ramifications[j][1][0], c_i
                        )
                        ramifications[j][0][1], ramifications[j][1][1] = self.fp.fp_cswap(
                            ramifications[j][0][1], ramifications[j][1][1], c_i
                        )

                        if self.formula_name == 'tvelu' or (
                            self.formula_name == 'hvelu'
                            and self.formula.global_L[pos] <= self.formula.HYBRID_BOUND
                        ):
                            # This branchs corresponds with the use of the tradicional velu's formulaes
                            # T_- (or T_)
                            ramifications[j][0] = self.formula.xEVAL(ramifications[j][0], pos)
                            # T_+ or (T_+)
                            ramifications[j][1] = self.formula.xEVAL(ramifications[j][1], pos)
                        else:
                            # This branchs corresponds with the use of the new velu's formulaes
                            # T_- (or T_)
                            ramifications[j][0] = self.formula.xEVAL(
                                ramifications[j][0], E_prev
                            )
                            # T_+ or (T_+)
                            ramifications[j][1] = self.formula.xEVAL(
                                ramifications[j][1], E_prev
                            )

                        # Multiplying by the current small odd prime number in order to decrease the order of both points T_+ and T_-
                        ramifications[j][1] =  self.curve.xMUL(ramifications[j][1], E_i, pos)

                        # Undo swap of points in E[\pi - 1] and E[\pi + 1]
                        ramifications[j][0][0], ramifications[j][1][0] = self.fp.fp_cswap(
                            ramifications[j][0][0], ramifications[j][1][0], c_i
                        )
                        ramifications[j][0][1], ramifications[j][1][1] = self.fp.fp_cswap(
                            ramifications[j][0][1], ramifications[j][1][1], c_i
                        )

                    b_i = (
                        isequal[u[pos] == 1] ^ isequal[u[pos] == -1]
                    )  # 0 if u[pos] != +1,-1; otherwise, 1 [This ask must be performed in constant-time]
                    v[pos] -= 1
                    u[pos] -= s_i * (
                        1 + b_i
                    )  # reduced by 1 unit if it was different from +1 and -1; otherwise, reduced by 2 units

                else:

                    for j in range(0, len(moves) - 1, 1):

                        ramifications[j][0][0], ramifications[j][1][0] = self.fp.fp_cswap(
                            ramifications[j][0][0], ramifications[j][1][0], c_i
                        )
                        ramifications[j][0][1], ramifications[j][1][1] = self.fp.fp_cswap(
                            ramifications[j][0][1], ramifications[j][1][1], c_i
                        )
                        ramifications[j][1] =  self.curve.xMUL(ramifications[j][1], E_i, pos)
                        ramifications[j][0][0], ramifications[j][1][0] = self.fp.fp_cswap(
                            ramifications[j][0][0], ramifications[j][1][0], c_i
                        )
                        ramifications[j][0][1], ramifications[j][1][1] = self.fp.fp_cswap(
                            ramifications[j][0][1], ramifications[j][1][1], c_i
                        )
            else:
                # This branch only depends on randomness

                # At this step, ramifications[-1] is the i-th leaf
                ramifications[-1][0][0], ramifications[-1][1][0] = self.fp.fp_cswap(
                    ramifications[-1][0][0], ramifications[-1][1][0], c_i
                )
                ramifications[-1][0][1], ramifications[-1][1][1] = self.fp.fp_cswap(
                    ramifications[-1][0][1], ramifications[-1][1][1], c_i
                )

                if  self.curve.isinfinity(ramifications[-1][0]) == False:

                    for j in range(0, len(moves) - 1, 1):

                        ramifications[j][0] =  self.curve.xMUL(ramifications[j][0], E_i, pos)
                        ramifications[j][1] =  self.curve.xMUL(ramifications[j][1], E_i, pos)
                else:

                    for j in range(0, len(moves) - 1, 1):

                        ramifications[j][0][0], ramifications[j][1][0] = self.fp.fp_cswap(
                            ramifications[j][0][0], ramifications[j][1][0], c_i
                        )
                        ramifications[j][0][1], ramifications[j][1][1] = self.fp.fp_cswap(
                            ramifications[j][0][1], ramifications[j][1][1], c_i
                        )
                        ramifications[j][1] =  self.curve.xMUL(ramifications[j][1], E_i, pos)
                        ramifications[j][0][0], ramifications[j][1][0] = self.fp.fp_cswap(
                            ramifications[j][0][0], ramifications[j][1][0], c_i
                        )
                        ramifications[j][0][1], ramifications[j][1][1] = self.fp.fp_cswap(
                            ramifications[j][0][1], ramifications[j][1][1], c_i
                        )

            moves.pop()
            ramifications.pop()

        pos = self.formula.global_L.index(L[0])  # Current element of self.formula.global_L to be required
        s_i = sign(e[pos])  # Sign of e[pos]
        c_i = (s_i + 1) // 2  # Constant-swap of T_+ and T_-

        # T_+ or T_- ?
        # Root
        ramifications[0][0][0], ramifications[0][1][0] = self.fp.fp_cswap(
            ramifications[0][0][0], ramifications[0][1][0], c_i
        )
        ramifications[0][0][1], ramifications[0][1][1] = self.fp.fp_cswap(
            ramifications[0][0][1], ramifications[0][1][1], c_i
        )

        if  self.curve.isinfinity(ramifications[0][0]) == False:

            if m[pos] > 0:
                if self.formula_name != 'tvelu':
                    # This branchs corresponds with the use of the new velu's formulaes

                    if self.verbose:
                        self.formula.set_parameters_velu(sJ_list[pos], sI_list[pos], pos)

                    else:
                        # -------------------------------------------------------------
                        # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                        # These paramters are required in KPs, xISOG, and xEVAL
                        if self.formula.global_L[pos] == 3:
                            b = 0
                            c = 0
                        else:
                            b = int(floor(sqrt(self.formula.global_L[pos] - 1) / 2.0))
                            c = int(floor((self.formula.global_L[pos] - 1.0) / (4.0 * b)))

                        self.formula.set_parameters_velu(b, c, pos)

                self.formula.KPs(ramifications[0][0], E_i, pos)
                E_i = self.formula.xISOG(E_i, pos)

                b_i = (
                    isequal[u[pos] == 1] ^ isequal[u[pos] == -1]
                )  # 0 if u[pos] != +1,-1; otherwise, 1 [This ask must be performed in constant-time]
                v[pos] -= 1
                u[pos] -= s_i * (
                    1 + b_i
                )  # reduced by 1 unit if it was different from +1 and -1; otherwise, reduced by 2 units

        return E_i, v, u




    def geometric_serie(self, m, l):
        """
        geometric_serie()
        inputs: and integer m, and a prime number l
        output: the nearest integer to
                      l
                m x -----
                    l - 1
        """
        l_float = float(l)
        m_float = float(m)
        return floor((m_float * l_float) / (l_float - 1.0) + 0.5)


    def filtered(self, List, sublist):
        """
        filtered()
        inputs : a list L and a sublist SL of L
        output : L \ SL
        """
        return [e for e in List if e not in sublist]

    def rounds(self, e, n):
        """
        rounds()
        inputs : an integer vector (maximum number of isogeny constructions to be performed),
                 and the length of the vector
        output : the subset of (indexes of) the small odd primes that determines the optimal
                 strategy to be used, the number of times that each strategy will be used, and
                 the complement of each subset (with respect to the set of all the small odd primes)
        """
        tmp_N = range(n)
        tmp_e = list(e)
        rounds_out = []
        sublists_L = []
        sublists_C = []
        while [e_i for e_i in tmp_e if e_i > 0] != []:
            e_min = min([e_i for e_i in tmp_e if e_i > 0])
            rounds_out.append(e_min)
            sublists_L.append([i for i in tmp_N if tmp_e[i] >= e_min])
            sublists_C.append(self.filtered(tmp_N, sublists_L[len(sublists_L) - 1]))
            tmp_e = [(tmp_e[i] - e_min) for i in tmp_N]

        return rounds_out, sublists_L, sublists_C

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
                for ii in range(0, self.curve.fp.exponent_of_two, 1):
                    T_p =  self.curve.xDBL(T_p, E_k)
                    T_m =  self.curve.xDBL(T_m, E_k)

                for l in R[j]:
                    T_p =  self.curve.xMUL(T_p, E_k, self.formula.global_L.index(l))
                    T_m =  self.curve.xMUL(T_m, E_k, self.formula.global_L.index(l))

                E_k, m, e = self.evaluate_strategy(
                    E_k, list([list(T_m), list(T_p)]), L[j], St[j], len(L[j]), m, e
                )

        # Multiplicative strategy on the set of unreached small odd prime numbers
        unreached_sop = [self.formula.global_L[i] for i in range(len(self.formula.global_L)) if m[i] > 0]
        remainder_sop = [l for l in self.formula.global_L if l not in unreached_sop]

        while len(unreached_sop) > 0:

            T_p, T_m = self.curve.elligator(E_k)
            for ii in range(0, self.curve.fp.exponent_of_two, 1):
                T_p =  self.curve.xDBL(T_p, E_k)
                T_m =  self.curve.xDBL(T_m, E_k)

            for l in remainder_sop:
                T_p =  self.curve.xMUL(T_p, E_k, self.formula.global_L.index(l))
                T_m =  self.curve.xMUL(T_m, E_k, self.formula.global_L.index(l))

            current_n = len(unreached_sop)
            E_k, m, e = self.evaluate_strategy(
                E_k,
                list([list(T_m), list(T_p)]),
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
                if m[self.formula.global_L.index(unreached_sop[k])] > 0
            ]
            tmp_remainder = [
                unreached_sop[k]
                for k in range(current_n)
                if m[self.formula.global_L.index(unreached_sop[k])] == 0
            ]

            unreached_sop = list(tmp_unreached)  # Removing elements from the batch
            remainder_sop = (
                remainder_sop + tmp_remainder
            )  # Adding elements to the complement of the batch

        return E_k



    def strategy_block_cost(self, L, e):
        """
        Next function computes the expected cost of our approach by assuming we have full torsion points
        """
        elligator_cost = numpy.array([7.0, 3.0, 10.0])  # Elligator cost
        mul_fp_by_four = (
            numpy.array([4.0, 2.0, 4.0]) * self.curve.fp.exponent_of_two
        )  # Cost of computing x([2^exponent_of_two]P)

        n = len(L)
        e_prime = [self.geometric_serie(e[k], L[k]) for k in range(n)]

        tmp_r, tmp_Ls, tmp_Cs = self.rounds(e_prime, n)

        C_e = numpy.array([0.0, 0.0, 0.0])
        S_out = []
        L_out = []
        R_out = []
        for j in range(len(tmp_r)):

            R_out.append([L[k] for k in tmp_Cs[j]])
            L_out.append([L[k] for k in tmp_Ls[j]])

            bo_C = 2.0 * sum([C_xMUL[self.formula.global_L.index(L[k])] for k in tmp_Cs[j]])
            S_tmp, go_C = self.dynamic_programming_algorithm(
                [L[k] for k in tmp_Ls[j]], len(tmp_Ls[j])
            )

            S_out.append(S_tmp)
            C_e = (
                C_e
                + (go_C + bo_C + elligator_cost + 2.0 * mul_fp_by_four) * tmp_r[j]
            )

        return C_e, L_out, R_out, S_out, tmp_r
