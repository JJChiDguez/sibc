from random import SystemRandom
import numpy
from sympy import symbols, floor, sqrt, sign, log
from pkg_resources import resource_filename

from sidh.math import isequal, bitlength, hamming_weight
from sidh.constants import parameters
from sidh.common import attrdict


def Gae(prime, tuned, curve, formula):
    '''
    >>> from sidh.bsidh.strategy import Gae
    >>> from sidh.bsidh.montgomery import MontgomeryCurve
    >>> from sidh.bsidh.hvelu import Hvelu
    >>> import pdb
    >>> prime = 'b2'
    >>> tuned = False
    >>> multievaluation = False
    >>> curve = MontgomeryCurve(prime)
    >>> formula = Hvelu(curve, tuned, multievaluation)
    >>> gae = Gae(prime, tuned, curve, formula)
    >>> sk_a = 0x642DCCC20D71FAFDFBA18D94E19777E8601494E718CB04E330C4BFE0181C209
    >>> sk_b = 0x135EB57C05FD58E80531C7CDDE36F2A7BBA88C55E8A70A0A97917D554AFA1EB6
    >>> pk_a = gae.pubkey_A(sk_a)
    >>> a_curve = curve.coeff(pk_a)
    >>> pk_b = gae.pubkey_B(sk_b)
    >>> b_curve = curve.coeff(pk_b)
    >>> ss_a = gae.dh_A(sk_a, pk_b)
    >>> curve_ss_a = curve.coeff(ss_a)
    >>> ss_b = gae.dh_B(sk_b, pk_a)
    >>> curve_ss_b = curve.coeff(ss_b)
    >>> assert curve_ss_a == curve_ss_b
    '''
    fp = curve.fp
    global_L = curve.L
    n = curve.n
    m = curve.p
    curve = curve
    random = SystemRandom()
    tuned_name = ('-classical','-suitable')[tuned]
    SQR, ADD = curve.SQR, curve.ADD

    SIDp = curve.SIDp
    SIDm = curve.SIDm

    f_name = 'data/strategies/bsidh-'+prime+'-'+formula.name+tuned_name
    f = open(resource_filename('sidh', f_name))
    # Corresponding to the list of Small Isogeny Degree, Lp := [l_0, ...,
    # l_{n-1}] [We need to include case l=2 and l=4]
    tmp = f.readline()
    tmp = [int(b) for b in tmp.split()]
    Sp = list(tmp)
    # Corresponding to the list of Small Isogeny Degree, Lm := [l_0, ...,
    # l_{n-1}]
    tmp = f.readline()
    tmp = [int(b) for b in tmp.split()]
    Sm = list(tmp)
    f.close()

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

    # These are for nonlocal in the pk/dh functions
    PA_b, QA_b, PQA_b = None, None, None
    PB_a, QB_a, PQB_a = None, None, None


    A = [[0x8, 0x0], [0x4, 0x0]]
    a = curve.coeff(A)

    # random_key() implements an uniform random integer sample [this functions should be modified]
    random = SystemRandom()
    random_key_A = lambda m=m: random.randint(0, m+1)
    random_key_B = lambda m=m: random.randint(0, m-1)
    random_key = lambda m=m: random.randint(0, m)

    # In order to achieve efficiency, the optimal strategies and their cost are saved in two global dictionaries (hash tables)
    S = {1: {}}  # Initialization of each strategy
    C = {1: {}}  # Initialization of the costs: 0.

    for i in range(n):
        j = global_L.index(curve.SID[i])
        S[1][tuple([global_L[j]])] = []
        # Strategy with a list with only one element (a small odd prime number l_i)
        C[1][tuple([global_L[j]])] = formula.C_xISOG[
            j
        ]  # Degree-l_i isogeny construction cost

    for i in range(2, n + 1):

        C[i] = {}
        S[i] = {}

    def dynamic_programming_algorithm(L, n):
        '''
        dynamic_programming_algorithm():
        inputs: the list of small odd primes to be processed and its length
        output: the optimal strategy and its cost of the input list of small odd primes
        '''
        nonlocal S, C
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
                                + C[  # Subtriangle on the left side with b leaves
                                    len(Tuple[b:])
                                ][Tuple[b:]]
                                + 1.0  # Subtriangle on the right side with (i - b) leaves
                                * sum(
                                    [curve.C_xMUL[global_L.index(t)] for t in Tuple[:b]]
                                )
                                + 1.0  # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with b leaves
                                * sum(
                                    [formula.C_xEVAL[global_L.index(t)] for t in Tuple[b:]]
                                ),  # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                            )
                            for b in range(1, i)
                        ]
                        b, C[i][Tuple] = min(
                            alpha, key=lambda t: curve.measure(t[1])
                        )  # We save the minimal cost corresponding to the triangle with leaves Tuple
                        S[i][Tuple] = (
                            [b] + S[i - b][Tuple[b:]] + S[b][Tuple[:b]]
                        )  # We save the optimal strategy corresponding to the triangle with leaves Tuple

            return S[n][tuple(L)], C[n][tuple(L)]  #

    def evaluate_strategy(EVAL, S_in, T_in, ST_in, E, P, L, strategy, n):
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

            pos = global_L.index(
                L[n - 1 - i]
            )  # Current element of global_L to be required

            # Reaching the vertex (n - 1 - i, i)
            # Vertical edges (scalar multiplications)
            prev = sum(moves)
            while prev < (n - 1 - i):

                moves.append(
                    strategy[k]
                )  # Number of vertical edges to be performed
                T = list(ramifications[-1])  # New ramification
                for j in range(prev, prev + strategy[k], 1):
                    T = curve.xMUL(T, E_i, global_L.index(L[j]))

                ramifications.append(list(T))
                prev += strategy[k]
                k += 1

            # Deciding which velu variant will be used
            if formula != 'tvelu':
                # This branchs corresponds with the use of the new velu's formulaes

                if tuned:
                    formula.set_parameters_velu(formula.sJ_list[pos], formula.sI_list[pos], pos)

                else:
                    # -------------------------------------------------------------
                    # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                    # These paramters are required in formula.KPs, formula.xISOG, and formula.xEVAL
                    if global_L[pos] <= 4:
                        b = 0
                        c = 0
                    else:
                        b = int(floor(sqrt(global_L[pos] - 1) / 2.0))
                        c = int(floor((global_L[pos] - 1.0) / (4.0 * b)))

                    formula.set_parameters_velu(b, c, pos)

            # Kernel Points computation
            formula.KPs(ramifications[-1], E_i, pos)

            # Isogeny construction
            ramifications[-1][0], E_i[0] = fp.fp2_cswap(
                ramifications[-1][0], E_i[0], global_L[pos] == 4
            )
            ramifications[-1][1], E_i[1] = fp.fp2_cswap(
                ramifications[-1][1], E_i[1], global_L[pos] == 4
            )
            C_i = formula.xISOG(E_i, pos)
            ramifications[-1][0], E_i[0] = fp.fp2_cswap(
                ramifications[-1][0], E_i[0], global_L[pos] == 4
            )
            ramifications[-1][1], E_i[1] = fp.fp2_cswap(
                ramifications[-1][1], E_i[1], global_L[pos] == 4
            )

            # Now, we proceed by perform horizontal edges (isogeny evaluations)
            for j in range(0, len(moves) - 1, 1):

                if (
                    formula.name == 'tvelu'
                    or (
                        formula.name == 'hvelu'
                        and global_L[pos] <= formula.HYBRID_BOUND
                    )
                    or (global_L[pos] == 4)
                ):
                    ramifications[j] = formula.xEVAL(ramifications[j], pos)
                else:
                    ramifications[j] = formula.xEVAL(ramifications[j], E_i)

            if EVAL:
                # Evaluating public points
                if (
                    formula.name == 'tvelu'
                    or (
                        formula.name == 'hvelu'
                        and global_L[pos] <= formula.HYBRID_BOUND
                    )
                    or (global_L[pos] == 4)
                ):

                    S_out = formula.xEVAL(S_out, pos)
                    T_out = formula.xEVAL(T_out, pos)
                    ST_out = formula.xEVAL(ST_out, pos)
                else:

                    S_out = formula.xEVAL(S_out, E_i)
                    T_out = formula.xEVAL(T_out, E_i)
                    ST_out = formula.xEVAL(ST_out, E_i)

            # Updating the Montogmery curve coefficients
            E_i = [list(C_i[0]), list(C_i[1])]

            moves.pop()
            ramifications.pop()

        pos = global_L.index(L[0])  # Current element of global_L to be required

        if formula.name != 'tvelu':
            # This branchs corresponds with the use of the new velu's formulaes

            if tuned:
                formula.set_parameters_velu(formula.sJ_list[pos], formula.sI_list[pos], pos)

            else:
                # -------------------------------------------------------------
                # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                # These paramters are required in formula.KPs, formula.xISOG, and formula.xEVAL
                if global_L[pos] <= 4:
                    b = 0
                    c = 0
                else:
                    b = int(floor(sqrt(global_L[pos] - 1) / 2.0))
                    c = int(floor((global_L[pos] - 1.0) / (4.0 * b)))

                formula.set_parameters_velu(b, c, pos)

        # Kernel Points computations
        formula.KPs(ramifications[0], E_i, pos)

        # Isogeny construction
        ramifications[0][0], E_i[0] = fp.fp2_cswap(
            ramifications[0][0], E_i[0], global_L[pos] == 4
        )
        ramifications[0][1], E_i[1] = fp.fp2_cswap(
            ramifications[0][1], E_i[1], global_L[pos] == 4
        )
        C_i = formula.xISOG(E_i, pos)
        ramifications[0][0], E_i[0] = fp.fp2_cswap(
            ramifications[0][0], E_i[0], global_L[pos] == 4
        )
        ramifications[0][1], E_i[1] = fp.fp2_cswap(
            ramifications[0][1], E_i[1], global_L[pos] == 4
        )

        if EVAL:
            # Evaluating public points
            if (
                formula.name == 'tvelu'
                or (formula.name == 'hvelu' and global_L[pos] <= formula.HYBRID_BOUND)
                or (global_L[pos] == 4)
            ):

                S_out = formula.xEVAL(S_out, pos)
                T_out = formula.xEVAL(T_out, pos)
                ST_out = formula.xEVAL(ST_out, pos)

            else:

                S_out = formula.xEVAL(S_out, E_i)
                T_out = formula.xEVAL(T_out, E_i)
                ST_out = formula.xEVAL(ST_out, E_i)

        # Updating the Montogmery curve coefficients
        E_i = [list(C_i[0]), list(C_i[1])]

        return E_i, S_out, T_out, ST_out

    def pubkey_A(sk_a):
        nonlocal PB_a, QB_a, PQB_a
        Ra = curve.Ladder3pt(sk_a, PA, QA, PQA, A)
        pk_a, PB_a, QB_a, PQB_a = evaluate_strategy( True, PB, QB, PQB, A,
                Ra, curve.SIDp[::-1], Sp, len(curve.SIDp))
        #a_curve = curve.coeff(pk_a)
        return pk_a

    def pubkey_B(sk_b):
        nonlocal PA_b, QA_b, PQA_b
        Rb = curve.Ladder3pt(sk_b, PB, QB, PQB, A)
        pk_b, PA_b, QA_b, PQA_b = evaluate_strategy( True, PA, QA, PQA, A,
                Rb, curve.SIDm[::-1], Sm, len(curve.SIDm))
        #b_curve = curve.coeff(pk_b)
        return pk_b

    def dh_A(sk_a, pk_b):
        # sk here is alice's secret key
        # pk_b here is from bob (not processed by coeff)
        nonlocal PA_b, QA_b, PQA_b
        RB_a = curve.Ladder3pt(sk_a, PA_b, QA_b, PQA_b, pk_b)
        ss_a, _, _, _ = evaluate_strategy( False, PB, QB, PQB, pk_b,
                RB_a, curve.SIDp[::-1], Sp, len(curve.SIDp))
        #ss_a_curve = curve.coeff(ss_a)
        #return ss_a_curve
        return ss_a

    def dh_B(sk_b, pk_a):
        # sk_b here is bob's secret key
        # pk_a here is from alice (not processed by coeff)
        nonlocal PB_a, QB_a, PQB_a
        RA_b = curve.Ladder3pt(sk_b, PB_a, QB_a, PQB_a, pk_a)
        ss_b, _, _, _ = evaluate_strategy( False, PA, QA, PQA, pk_a,
                RA_b, curve.SIDm[::-1], Sm, len(curve.SIDm))
        #ss_b_curve = curve.coeff(ss_b)
        #return ss_b_curve
        return ss_b

    return attrdict(locals())

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
