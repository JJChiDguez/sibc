from random import SystemRandom
from pkg_resources import resource_filename

from functools import reduce
from math import floor, sqrt

from sibc.math import isequal, bitlength, hamming_weight, cswap, sign
from sibc.constants import parameters

class Strategy(object):
    def __init__(self, prime, curve, formula):
        self.random = SystemRandom()

        self.curve = curve
        self.prime = prime
        self.formula = formula
        self.formula_name = formula.name
        self.field = self.curve.field
        self.L = curve.L
        self.two = self.curve.two
        self.three = self.curve.three

        # Reading public generators points
        f = open(resource_filename('sibc', 'data/gen/' + curve.model + '/sidh/' + prime))

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
        # To check the stored strategies format/name-file
        if not formula.uninitialized:
            file_path = (
                "data/strategies/"
                + curve.model
                + '/sidh/'
                + 'sidh'
                + '-'
                + prime
            )
            file_path = resource_filename('sibc', file_path)
            f = open(file_path)
            # Corresponding to the list of Small Isogeny Degree, Lp := [l_0, ...,
            # l_{n-1}] [We need to include case l=2 and l=4]
            tmp = f.readline()
            tmp = [int(b) for b in tmp.split()]
            self.S2 = list(tmp)
            # Corresponding to the list of Small Isogeny Degree, Lm := [l_0, ...,
            # l_{n-1}]
            tmp = f.readline()
            tmp = [int(b) for b in tmp.split()]
            self.S3 = list(tmp)
            f.close()
        ######################################################################################################################

    def dynamic_programming_algorithm(self, ell, e):
        """
        dynamic_programming_algorithm():
        inputs: the small prime (either 2 or 3) and its exponent
        output: the optimal strategy and its cost
        """
        # My intuition says, e = 1 mod 2 can be improved
        n = {2:(e//2), 3:e}[ell]
        p = 12          # Both of x([4]P) and x([3]P) cost 4M + 8S, assuming M=S we have a cost of 12 multiplications
        q = {2:8, 3:6}[ell]  # cost of a degree-4 and degree-3 isogeny evaluation (6M + 2S + 6a and 4M + 2S + 4a, respectively)

        S = {1:[]}
        C = {1:0 }
        for i in range(2, n+1):
            b, cost = min(((b, C[i-b] + C[b] + b*p + (i-b)*q) for b in range(1,i)), key=lambda t: t[1])
            S[i] = [b] + S[i-b] + S[b]
            C[i] = cost

        return S[n], C[n]

    def random_scalar_A(self): return self.random.randint(0, 2 ** self.curve.two)
    def random_scalar_B(self): return self.random.randint(0, 3 ** self.curve.three)

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
            2,
            self.S2,
            self.two
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
        # A + 2C = 8, and A - 2C = 4
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
            3,
            self.S3,
            self.three
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
            2,
            self.S2,
            self.two
        )
        return ss_a

    def strategy_B(self, sk_b, pk_a):
        (PB_a, QB_a, PQB_a) = pk_a
        # Notice, get_A() returns (A'+2C:4C)
        A = self.curve.get_A(
            [PB_a, self.field(1)],
            [QB_a, self.field(1)],
            [PQB_a, self.field(1)]
        )
        # Recall, 3-isogenies work with (A'+2C: A'-2C) = (A[0] : A[0] - A[1])
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
            [A[0], A[0] - A[1]],
            RA_b,
            3,
            self.S3,
            self.three
        )
        # Now, we map from (A'+2C:A'-2C) into (A'+2C:4C) = (ss_b[0] : ss_b[0] - ss_b[1])
        return [ss_b[0], ss_b[0] - ss_b[1]]

    def evaluate_strategy(self, EVAL, S_in, T_in, ST_in, E, P, ell, strategy, n):
        '''
        evaluate_strategy():
                 primes;
        output : the projective Montgomery constants a24:= a + 2c and c24:=4c where E': y^2 = x^3 + (a/c)*x^2 + x
                 is E / <P>
        '''
        exp = {2:n//2, 3:n}[ell]
        xmul = {2:self.curve.xdbl_twice, 3:self.curve.xtpl}[ell]        # quadrupling or tripling of points
        kps = {2:self.formula.kps_4, 3:self.formula.kps_3}[ell]         # kernel computation (cost of two additions)
        xisog = {2:self.formula.xisog_4, 3:self.formula.xisog_3}[ell]   # degree-2 or degree-3 isogeny construction
        xeval = {2:self.formula.xeval_4, 3:self.formula.xeval_3}[ell]   # degree-2 or degree-3 isogeny evaluation

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

        assert( (exp - 1) == len(strategy))
        if (ell == 2) and (n % 2 == 1):
            T = list(ramifications[0])  # New ramification
            for j in range(0, n - 1,1):
                T = self.curve.xdbl(T, E_i)
            
            self.formula.kps_2(T)
            E_i = self.formula.xisog_2(T)
            ramifications[0] = self.formula.xeval_2(ramifications[0])

            if EVAL:
                # Evaluating public points
                S_out = self.formula.xeval_2(S_out)
                T_out = self.formula.xeval_2(T_out)
                ST_out = self.formula.xeval_2(ST_out)

        for i in range(0, exp - 1, 1):

            # Reaching the vertex (n - 1 - i, i)
            # Vertical edges (scalar multiplications)
            prev = sum(moves)
            while prev < (exp - 1 - i):

                moves.append(
                    strategy[k]
                )  # Number of vertical edges to be performed
                T = list(ramifications[-1])  # New ramification
                for j in range(prev, prev + strategy[k], 1):
                    T = xmul(T, E_i)

                ramifications.append(list(T))
                prev += strategy[k]
                k += 1

            # Kernel Points computation
            kps(ramifications[-1])
            # Isogeny construction
            E_i = xisog(ramifications[-1])
            # Now, we proceed by perform horizontal edges (isogeny evaluations)
            for j in range(0, len(moves) - 1, 1):
                ramifications[j] = xeval(ramifications[j])

            if EVAL:
                # Evaluating public points
                S_out = xeval(S_out)
                T_out = xeval(T_out)
                ST_out = xeval(ST_out)

            moves.pop()
            ramifications.pop()

        # Kernel Points computations
        kps(ramifications[0])
        # Isogeny construction
        E_i = xisog(ramifications[0])
        if EVAL:
            # Evaluating public points
            S_out = xeval(S_out)
            T_out = xeval(T_out)
            ST_out = xeval(ST_out)

        return E_i, S_out, T_out, ST_out
