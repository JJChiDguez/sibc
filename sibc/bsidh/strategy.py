from random import SystemRandom
from pkg_resources import resource_filename

from functools import reduce
from math import floor, sqrt

from sibc.math import isequal, bitlength, hamming_weight, cswap, sign
from sibc.constants import parameters

class Strategy(object):
    def __init__(self, prime, tuned, curve, formula):
        self.random = SystemRandom()

        # In order to achieve efficiency, the optimal strategies and their cost are saved in two global dictionaries (hash tables)
        self.S = {1: {}}  # Initialization of each strategy
        self.C = {1: {}}  # Initialization of the costs: 0.

        self.curve = curve
        self.prime = prime
        self.formula = formula
        self.formula_name = formula.name
        self.field = self.curve.field
        self.tuned = tuned
        self.L = curve.L
        self.SIDp = reduce(lambda x, y: x + y, [[self.curve.Lp[i]] * self.curve.Ep[i] for i in range(0, self.curve.np, 1)])
        self.SIDm = reduce(lambda x, y: x + y, [[self.curve.Lm[i]] * self.curve.Em[i] for i in range(0, self.curve.nm, 1)])

        self.c_xmul = self.curve.c_xmul
        n = curve.n
        for i in range(n):
            self.S[1][tuple([self.L[i]])] = []
            # Strategy with a list with only one element (a small odd prime number l_i)
            self.C[1][tuple([self.L[i]])] = self.formula.c_xisog[i]
            # For catching the weigth of horizontal edges of the form [(0,j),(0,j+1)]
        for i in range(2, len(self.SIDp) + len(self.SIDm) + 1):
            self.C[i] = {}
            self.S[i] = {}

        # Reading public generators points
        f = open(resource_filename('sibc', 'data/gen/' + curve.model + '/bsidh/' + prime))

        # x(PA), x(QA) and x(PA - QA)
        PQA = f.readline()
        PQA = [int(x, 16) for x in PQA.split()]
        self.PA = self.field(PQA[0:2])
        self.QA = self.field(PQA[2:4])
        self.PQA= self.field(PQA[4:6])

        # x(PB), x(QB) and x(PB - QB)
        PQB = f.readline()
        PQB = [int(x, 16) for x in PQB.split()]
        self.PB = self.field(PQB[0:2])
        self.QB = self.field(PQB[2:4])
        self.PQB= self.field(PQB[4:6])

        f.close()

        if not formula.uninitialized:
            file_path = (
                "data/strategies/"
                + curve.model
                + '/'
                + 'bsidh/bsidh'
                + '-'
                + prime
                + '-'
                + formula.name
                + '-'
                + formula.multievaluation_name
                + formula.tuned_name
            )
            file_path = resource_filename('sibc', file_path)
            f = open(file_path)
            # Corresponding to the list of Small Isogeny Degree, Lp := [l_0, ...,
            # l_{n-1}] [We need to include case l=2 and l=4]
            tmp = f.readline()
            tmp = [int(b) for b in tmp.split()]
            self.Sp = list(tmp)
            # Corresponding to the list of Small Isogeny Degree, Lm := [l_0, ...,
            # l_{n-1}]
            tmp = f.readline()
            tmp = [int(b) for b in tmp.split()]
            self.Sm = list(tmp)
            f.close()
        ######################################################################################################################

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
                                ),
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
                self.C[n][tuple(L)],
            )  # The weight of the horizontal edges [(0,n-1),(0,n)] must be equal to c_xisog[self.formula.L.index(L[0])].

    def random_scalar_A(self): return self.random.randint(0, self.curve.p + 1)
    def random_scalar_B(self): return self.random.randint(0, self.curve.p - 1)

    def strategy_at_6_A(self, sk_a):
        # a24 = (A + 2) / 4 = (6 + 2) / 4 = 2
        Ra = self.curve.Ladder3pt(
            sk_a,
            [self.PA, self.field(1)],
            [self.QA, self.field(1)],
            [self.PQA, self.field(1)],
            self.field(2)
        )
        _, P, Q, PQ = self.evaluate_strategy(
            True,
            [self.PB, self.field(1)],
            [self.QB, self.field(1)],
            [self.PQB, self.field(1)],
            [self.field(8), self.field(4)],
            Ra,
            self.SIDp[::-1],
            self.Sp,
            len(self.SIDp)
        )
        # --- from proj to affine
        inv = (P[1] * Q[1])
        inv *= PQ[1]
        inv **= -1
        # --- P
        PB_a = (Q[1] * PQ[1])
        PB_a *= inv
        PB_a *= P[0]
        # --- Q
        QB_a = (P[1] * PQ[1])
        QB_a *= inv
        QB_a *= Q[0]
        # --- PQ
        PQB_a = (P[1] * Q[1])
        PQB_a *= inv
        PQB_a *= PQ[0]
        return (PB_a, QB_a, PQB_a)

    def strategy_at_6_B(self, sk_b):
        # a24 = (A + 2) / 4 = (6 + 2) / 4 = 2
        Rb = self.curve.Ladder3pt(
            sk_b,
            [self.PB, self.field(1)],
            [self.QB, self.field(1)],
            [self.PQB, self.field(1)],
            self.field(2)
        )
        _, P, Q, PQ = self.evaluate_strategy(
            True,
            [self.PA, self.field(1)],
            [self.QA, self.field(1)],
            [self.PQA, self.field(1)],
            [self.field(8), self.field(4)],
            Rb,
            self.SIDm[::-1],
            self.Sm,
            len(self.SIDm)
        )
        # --- from proj to affine
        inv = (P[1] * Q[1])
        inv *= PQ[1]
        inv **= -1
        # --- P
        PA_b = (Q[1] * PQ[1])
        PA_b *= inv
        PA_b *= P[0]
        # --- Q
        QA_b = (P[1] * PQ[1])
        QA_b *= inv
        QA_b *= Q[0]
        # --- PQ
        PQA_b = (P[1] * Q[1])
        PQA_b *= inv
        PQA_b *= PQ[0]
        return (PA_b, QA_b, PQA_b)

    def strategy_A(self, sk_a, pk_b):
        (PA_b, QA_b, PQA_b) = pk_b
        A = self.curve.get_A(
            [PA_b, self.field(1)],
            [QA_b, self.field(1)],
            [PQA_b, self.field(1)]
        )
        #assert self.curve.issupersingular(A), "non-supersingular input curve"
        a24 = A[1] ** -1
        a24 *= A[0]
        RB_a = self.curve.Ladder3pt(
            sk_a,
            [PA_b, self.field(1)],
            [QA_b, self.field(1)],
            [PQA_b, self.field(1)],
            a24
        )
        ss_a, _, _, _ = self.evaluate_strategy(
            False,
            None,
            None,
            None,
            A,
            RB_a,
            self.SIDp[::-1],
            self.Sp,
            len(self.SIDp)
        )
        return ss_a

    def strategy_B(self, sk_b, pk_a):
        (PB_a, QB_a, PQB_a) = pk_a
        A = self.curve.get_A(
            [PB_a, self.field(1)],
            [QB_a, self.field(1)],
            [PQB_a, self.field(1)]
        )
        #assert self.curve.issupersingular(A), "non-supersingular input curve"
        a24 = A[1] ** -1
        a24 *= A[0]
        RA_b = self.curve.Ladder3pt(
            sk_b,
            [PB_a, self.field(1)],
            [QB_a, self.field(1)],
            [PQB_a, self.field(1)],
            a24
        )
        ss_b, _, _, _ = self.evaluate_strategy(
            False,
            None,
            None,
            None,
            A,
            RA_b,
            self.SIDm[::-1],
            self.Sm,
            len(self.SIDm)
        )
        return ss_b

    def evaluate_strategy(self, EVAL, S_in, T_in, ST_in, E, P, L, strategy, n):
        '''
        evaluate_strategy():
                 primes;
        output : the projective Montgomery constants a24:= a + 2c and c24:=4c where E': y^2 = x^3 + (a/c)*x^2 + x
                     is E / <P>
        '''

        ramifications = []
        moves = [
            0
        ]  # moves: this list determines whether an isogeny construction must be performed
        k = 0  # k: current element of the strategy

        ramifications.append(list(P))
        E_i = list(E)

        if EVAL:
            # Public points to be evaluated
            S_out = list(
                S_in
            )  # x(S) should be a torsion point with not common factors in L
            T_out = list(
                T_in
            )  # x(T) should be a torsion point with not common factors in L
            ST_out = list(
                ST_in
            )  # x(S - T) should be a torsion point with not common factors in L
        else:
            S_out = None
            T_out = None
            ST_out = None

        assert len(strategy) == (n - 1)
        for i in range(len(strategy)):

            pos = self.L.index(
                L[n - 1 - i]
            )  # Current element of self.L to be required

            # Reaching the vertex (n - 1 - i, i)
            # Vertical edges (scalar multiplications)
            prev = sum(moves)
            while prev < (n - 1 - i):

                moves.append(
                    strategy[k]
                )  # Number of vertical edges to be performed
                T = list(ramifications[-1])  # New ramification
                for j in range(prev, prev + strategy[k], 1):
                    T = self.curve.xmul(T, E_i, self.L.index(L[j]))

                ramifications.append(list(T))
                prev += strategy[k]
                k += 1

            # Deciding which velu variant will be used
            if self.formula_name != 'tvelu':
                # This branchs corresponds with the use of the new velu's formulaes

                if self.tuned:
                    self.formula.set_parameters_velu(self.formula.sJ_list[pos], self.formula.sI_list[pos], pos)

                else:
                    # -------------------------------------------------------------
                    # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                    # These paramters are required in self.formula.kps, self.formula.xisog, and self.formula.xeval
                    if self.L[pos] <= 4:
                        b = 0
                        c = 0
                    else:
                        b = int(floor(sqrt(self.L[pos] - 1) / 2.0))
                        c = int(floor((self.L[pos] - 1.0) / (4.0 * b)))

                    if self.formula_name != 'tvelu':
                        self.formula.set_parameters_velu(b, c, pos)

            # Kernel Points computation
            self.formula.kps(ramifications[-1], E_i, pos)

            if not EVAL:
                # Isogeny construction
                ramifications[-1][0], E_i[0] = cswap(
                    ramifications[-1][0], E_i[0], self.L[pos] == 4
                )
                ramifications[-1][1], E_i[1] = cswap(
                    ramifications[-1][1], E_i[1], self.L[pos] == 4
                )
                C_i = self.formula.xisog(E_i, pos)
                ramifications[-1][0], E_i[0] = cswap(
                    ramifications[-1][0], E_i[0], self.L[pos] == 4
                )
                ramifications[-1][1], E_i[1] = cswap(
                    ramifications[-1][1], E_i[1], self.L[pos] == 4
                )               

            # Now, we proceed by perform horizontal edges (isogeny evaluations)
            for j in range(0, len(moves) - 1, 1):

                if (
                    self.formula_name == 'tvelu'
                    or (
                        self.formula_name == 'hvelu'
                        and self.L[pos] <= self.formula.HYBRID_BOUND
                    )
                    or (self.L[pos] == 4)
                ):
                    ramifications[j] = self.formula.xeval(ramifications[j], pos)
                else:
                    ramifications[j] = self.formula.xeval(ramifications[j], E_i)

            if EVAL:
                # Evaluating public points
                if (
                    self.formula_name == 'tvelu'
                    or (
                        self.formula_name == 'hvelu'
                        and self.L[pos] <= self.formula.HYBRID_BOUND
                    )
                    or (self.L[pos] == 4)
                ):

                    S_out = self.formula.xeval(S_out, pos)
                    T_out = self.formula.xeval(T_out, pos)
                    ST_out = self.formula.xeval(ST_out, pos)
                else:

                    S_out = self.formula.xeval(S_out, E_i)
                    T_out = self.formula.xeval(T_out, E_i)
                    ST_out = self.formula.xeval(ST_out, E_i)

                C_i = self.curve.get_A(S_out, T_out, ST_out)

            # Updating the Montogmery curve coefficients
            E_i = [self.field(C_i[0]), self.field(C_i[1])]

            moves.pop()
            ramifications.pop()

        pos = self.L.index(L[0])  # Current element of self.L to be required

        if self.formula_name != 'tvelu':
            # This branchs corresponds with the use of the new velu's formulaes

            if self.tuned:
                self.formula.set_parameters_velu(self.formula.sJ_list[pos], self.formula.sI_list[pos], pos)

            else:
                # -------------------------------------------------------------
                # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                # These paramters are required in self.formula.kps, self.formula.xisog, and self.formula.xeval
                if self.L[pos] <= 4:
                    b = 0
                    c = 0
                else:
                    b = int(floor(sqrt(self.L[pos] - 1) / 2.0))
                    c = int(floor((self.L[pos] - 1.0) / (4.0 * b)))

                self.formula.set_parameters_velu(b, c, pos)

        # Kernel Points computations
        self.formula.kps(ramifications[0], E_i, pos)

        if not EVAL:
            # Isogeny construction
            ramifications[0][0], E_i[0] = cswap(
                ramifications[0][0], E_i[0], self.L[pos] == 4
            )
            ramifications[0][1], E_i[1] = cswap(
                ramifications[0][1], E_i[1], self.L[pos] == 4
            )
            C_i = self.formula.xisog(E_i, pos)
            ramifications[0][0], E_i[0] = cswap(
                ramifications[0][0], E_i[0], self.L[pos] == 4
            )
            ramifications[0][1], E_i[1] = cswap(
                ramifications[0][1], E_i[1], self.L[pos] == 4
            )

        else:
            
            # Evaluating public points
            if (
                self.formula_name == 'tvelu'
                or (self.formula_name == 'hvelu' and self.L[pos] <= self.formula.HYBRID_BOUND)
                or (self.L[pos] == 4)
            ):

                S_out = self.formula.xeval(S_out, pos)
                T_out = self.formula.xeval(T_out, pos)
                ST_out = self.formula.xeval(ST_out, pos)

            else:

                S_out = self.formula.xeval(S_out, E_i)
                T_out = self.formula.xeval(T_out, E_i)
                ST_out = self.formula.xeval(ST_out, E_i)

            C_i = self.curve.get_A(S_out, T_out, ST_out)

        # Updating the Montogmery curve coefficients
        E_i = [self.field(C_i[0]), self.field(C_i[1])]

        return E_i, S_out, T_out, ST_out

