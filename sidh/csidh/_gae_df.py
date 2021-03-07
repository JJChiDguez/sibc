#if __formulaes == 'tvelu':
#    from sidh.csidh.tvelu import *
#elif __formulaes == 'svelu':
#    from sidh.csidh.svelu import *
#else:
#    from sidh.csidh.hvelu import *

import numpy
from sympy import symbols, floor, sqrt, sign

from sidh._math import isequal, bitlength, hamming_weight
from sidh.csidh._hvelu import Hvelu, affine_to_projective, measure, elligator, xDBL, xMUL, set_parameters_velu, KPs, xISOG, HYBRID_BOUND, xEVAL
from sidh.csidh._montgomery import isinfinity, coeff, validate
from sidh._fp import fp_cswap
from sidh.constants import parameters

# In order to achieve efficiency, the optimal strategies and their cost are saved in two global dictionaries (hash tables)
S = {1: {}}  # Initialization of each strategy
C = {1: {}}  # Initialization of the costs: 0.
L = None
n = None
basis = None
hvelu = None
C_xMUL = None

__formulaes ='hvelu'
__verbose = False

class Gae_df(object):
    """

    >>> sk_a = [ 9,  2,-14,  3, 11,-12,  4,  4,  0, 20,  8,  7,-12,-20, 23, 15,  5,  3, 15,-19,  7,-17,-19,  1, -2, 14,  0, -9,  2,  4, 11,  2,  7,  9,  9, -1,  5, -7,  5,  4, -4,  6,  6,  7, -8, -2,  0,  2, -6,  5, -2,  0, -2,  4, -5, -1, -5,  3,  3,  5,  3, -5,  5,  3, -2,  2, -4, -2,  0, -2,  2,  0, 2, -3 ]
    >>> pk_a = [0x2b84079dfe34daca1000b23de50ed8e581c2f8f174fdbcb03d3ca6c96361e731152f45bdd4959832de406e8f3c4f0b4949c4826af82b3a7362677a458196fbcf, 0x76599305d04fb32f8907b2279ee35be99786da7c055e4a922712a6b546d457fa8db9529049bbe500be47e23dae04ecd34e043264a02bb1917dfdf9fa540f233]
    >>> sk_b = [ 3,-16,-10,  1, 15, 20,-20,-22,-16,-22,  0,-19,  6, -4, -9, 13,-11, 13,-13, -1, 23, 21, -5, 13, -4, -2, 12, 15, -4,-10, -5,  0, 11,  1, -1, -1,  7,  1, -3,  6,  0,  2, -4, -5,  0,  2, -4, -2, -4, -5,  6,  2, -6, -4,  5, -5,  5, -3,  1,  3, -1, -5,  3, -5, -4,  2,  4,  2,  2,  4,  0, -2, 0, -3 ]
    >>> pk_b = [ 0x5e2fd48334e49b47bb754f88b3345c77d604eb1fadc29b35b931724459143abde2f22346b595e3b161d80f3659870f16d4983bfa58f5f2f9718d3b375c21d65c, 0x314b346a927c21051052b790809d895627ed8fbe4f008408d361223a97556ec0e6d3b544b0898daffcdbff5c5b409ccb5cc9e2edc95504fca54318071e28e054 ]
    >>> ss = 0x1ADB783878BA330BB2A842E7F8B3392329A2CD3B407900E4CF6A8F13B744BFFEFF617BDE2CEBBB9CE97D32BC6FC1BCE2D88381B03B3E13CFF0651EEA82D02937
    >>> g = Gae_df('p512', 'df', False)
    // Shortest Differential Addition Chains (SDAC) for each l_i;
    // SDAC's to be read from a file
    >>> g.dh(sk_a, pk_b) == g.dh(sk_b, pk_a) == ss
    True
    >>> hex(coeff(pk_a)).upper()[2:]
    '27AD85DDC08BF510F08A8562BA4909803675536A0BCE6250E3BED4A9401AC123FE75C18866625E9FCFCAF03D0927ED46665E153E786244DAAAC9A83075060C82'
    >>> hex(coeff(pk_b)).upper()[2:]
    '181C39753CCB4D3358E32B4471EE73EDC568846CA3B0B17571A09BD7373B4658251ADF466FF1FFB29D89B382184703C708F71497611A4B643BD984D847F3A430'
    >>> coeff(affine_to_projective(coeff(pk_a))) == coeff(pk_a)
    True
    >>> g.dh(sk_b, affine_to_projective(coeff(pk_a))) == ss
    True
    >>> list(map(hex,affine_to_projective(coeff(pk_a))))
    ['0x27ad85ddc08bf510f08a8562ba4909803675536a0bce6250e3bed4a9401ac123fe75c18866625e9fcfcaf03d0927ed46665e153e786244daaac9a83075060c84', '0x4']
    """
    def __init__(self, prime, style, verbose):
        global S
        global C
        global L
        global n
        global basis
        global hvelu
        global C_xMUL
        self.prime = prime
        self.hvelu = hvelu = Hvelu(prime, style, verbose)
        self.L = L = parameters['csidh'][prime]['L']
        self.n = n = parameters['csidh'][prime]['n']
        C_xMUL = hvelu.ML.C_xMUL
        for i in range(n):
            S[1][tuple([L[i]])] = []
            # Strategy with a list with only one element (a small odd prime number l_i)
            C[1][tuple([L[i]])] = hvelu.C_xISOG[i]
            # For catching the weigth of horizontal edges of the form [(0,j),(0,j+1)]
        for i in range(2, n + 1):
            C[i] = {}
            S[i] = {}
        ######################################################################################################################
        # Next functions are used for computing optimal bounds
        #self.basis = basis = numpy.eye(n, dtype=int)
        self.m = parameters['csidh'][prime]['df']['m']
        self.temporal_m = list(set(self.m))

    def dh(self, sk, pk):
        assert validate(pk), "public key does not validate"
        C_out, L_out, R_out, S_out, r_out = strategy_block_cost(L[::-1], self.m[::-1])
        ss = GAE(
            pk,
            sk,
            [L_out[0]],
            [R_out[0]],
            [S_out[0]],
            [self.temporal_m[-1]],
            parameters['csidh'][self.prime]['df']['m'],
        )
        return coeff(ss)

def random_key(m):
    """
    random_key(m) implements an uniform random sample from S(m_1) x S(m_2) x ... x S(m_n)
    """
    return [2 * (random.randint(0, m_i) - (m_i // 2)) - (m_i % 2) for m_i in m ]

def security(M, n):
    """
    security()
    inputs : the list M of maximum number of degree-(l_i) isogeny constructions
             to be performed, and the length of M
    output : bits of security that M brings
    """
    return sum(list([log(M[i] + 1, 2) for i in range(n)]))


def dynamic_programming_algorithm(L, n):
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
                                [C_xMUL[hvelu.global_L.index(t)] for t in Tuple[:b]]
                            )
                            + 2.0  # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with b leaves
                            * sum(
                                [hvelu.C_xEVAL[hvelu.global_L.index(t)] for t in Tuple[b:]]
                            )
                            + 1.0  # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                            * sum(
                                [C_xMUL[hvelu.global_L.index(t)] for t in Tuple[b:]]
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
                                    C_xMUL[hvelu.global_L.index(t)]
                                    for t in Tuple[: (i - 1)]
                                ]
                            )
                            + 2.0  # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with 1 leaf
                            * hvelu.C_xEVAL[hvelu.global_L.index(Tuple[i - 1])]
                            + 1.0  # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - 1) leaves
                            * C_xMUL[
                                hvelu.global_L.index(Tuple[i - 1])
                            ],  # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - 1) leaves
                        )
                    ]
                    b, C[i][Tuple] = min(
                        alpha, key=lambda t: measure(t[1])
                    )  # We save the minimal cost corresponding to the triangle with leaves Tuple
                    S[i][Tuple] = (
                        [b] + S[i - b][Tuple[b:]] + S[b][Tuple[:b]]
                    )  # We save the optimal strategy corresponding to the triangle with leaves Tuple

        return (
            S[n][tuple(L)],
            C[n][tuple(L)],
        )  # The weight of the horizontal edges [(0,n-1),(0,n)] must be equal to C_xISOG[hvelu.global_L.index(L[0])].



def evaluate_strategy(E, P, L, strategy, n, m, e):
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

        pos = hvelu.global_L.index(
            L[n - 1 - i]
        )  # Current element of hvelu.global_L to be required

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
                        xMUL(T[0], E_i, hvelu.global_L.index(L[j])),
                        xMUL(T[1], E_i, hvelu.global_L.index(L[j])),
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

            T[0][0], T[1][0] = fp_cswap(T[0][0], T[1][0], c_i)
            T[0][1], T[1][1] = fp_cswap(T[0][1], T[1][1], c_i)
            for j in range(prev, prev + strategy[k], 1):
                T[0] = list(
                    xMUL(T[0], E_i, hvelu.global_L.index(L[j]))
                )  # A single scalar multiplication is required

            T[0][0], T[1][0] = fp_cswap(T[0][0], T[1][0], c_i)
            T[0][1], T[1][1] = fp_cswap(T[0][1], T[1][1], c_i)

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
            ramifications[-1][0][0], ramifications[-1][1][0] = fp_cswap(
                ramifications[-1][0][0], ramifications[-1][1][0], c_i
            )  # XT_+ <-> XT_-
            ramifications[-1][0][1], ramifications[-1][1][1] = fp_cswap(
                ramifications[-1][0][1], ramifications[-1][1][1], c_i
            )  # ZT_+ <-> ZT_-

            if isinfinity(ramifications[-1][0]) == False:

                if __formulaes != 'tvelu':
                    # This branchs corresponds with the use of the new velu's formulaes

                    E_prev = list(E_i)
                    if __verbose:
                        set_parameters_velu(sJ_list[pos], sI_list[pos], pos)

                    else:
                        # -------------------------------------------------------------
                        # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                        # These paramters are required in KPs, xISOG, and xEVAL
                        if hvelu.global_L[pos] == 3:
                            b = 0
                            c = 0
                        else:
                            b = int(floor(sqrt(hvelu.global_L[pos] - 1) / 2.0))
                            c = int(floor((hvelu.global_L[pos] - 1.0) / (4.0 * b)))

                        set_parameters_velu(b, c, pos)

                KPs(ramifications[-1][0], E_i, pos)

                # New isogeny construction
                E_i = xISOG(E_i, pos)

                # Horizontal edges are performed
                for j in range(0, len(moves) - 1, 1):

                    # Swap points in E[\pi - 1] and E[\pi + 1]
                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(
                        ramifications[j][0][0], ramifications[j][1][0], c_i
                    )
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(
                        ramifications[j][0][1], ramifications[j][1][1], c_i
                    )

                    if __formulaes == 'tvelu' or (
                        __formulaes == 'hvelu'
                        and hvelu.global_L[pos] <= HYBRID_BOUND
                    ):
                        # This branchs corresponds with the use of the tradicional velu's formulaes
                        # T_- (or T_)
                        ramifications[j][0] = xEVAL(ramifications[j][0], pos)
                        # T_+ or (T_+)
                        ramifications[j][1] = xEVAL(ramifications[j][1], pos)
                    else:
                        # This branchs corresponds with the use of the new velu's formulaes
                        # T_- (or T_)
                        ramifications[j][0] = xEVAL(
                            ramifications[j][0], E_prev
                        )
                        # T_+ or (T_+)
                        ramifications[j][1] = xEVAL(
                            ramifications[j][1], E_prev
                        )

                    # Multiplying by the current small odd prime number in order to decrease the order of both points T_+ and T_-
                    ramifications[j][1] = xMUL(ramifications[j][1], E_i, pos)

                    # Undo swap of points in E[\pi - 1] and E[\pi + 1]
                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(
                        ramifications[j][0][0], ramifications[j][1][0], c_i
                    )
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(
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

                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(
                        ramifications[j][0][0], ramifications[j][1][0], c_i
                    )
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(
                        ramifications[j][0][1], ramifications[j][1][1], c_i
                    )
                    ramifications[j][1] = xMUL(ramifications[j][1], E_i, pos)
                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(
                        ramifications[j][0][0], ramifications[j][1][0], c_i
                    )
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(
                        ramifications[j][0][1], ramifications[j][1][1], c_i
                    )
        else:
            # This branch only depends on randomness

            # At this step, ramifications[-1] is the i-th leaf
            ramifications[-1][0][0], ramifications[-1][1][0] = fp_cswap(
                ramifications[-1][0][0], ramifications[-1][1][0], c_i
            )
            ramifications[-1][0][1], ramifications[-1][1][1] = fp_cswap(
                ramifications[-1][0][1], ramifications[-1][1][1], c_i
            )

            if isinfinity(ramifications[-1][0]) == False:

                for j in range(0, len(moves) - 1, 1):

                    ramifications[j][0] = xMUL(ramifications[j][0], E_i, pos)
                    ramifications[j][1] = xMUL(ramifications[j][1], E_i, pos)
            else:

                for j in range(0, len(moves) - 1, 1):

                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(
                        ramifications[j][0][0], ramifications[j][1][0], c_i
                    )
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(
                        ramifications[j][0][1], ramifications[j][1][1], c_i
                    )
                    ramifications[j][1] = xMUL(ramifications[j][1], E_i, pos)
                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(
                        ramifications[j][0][0], ramifications[j][1][0], c_i
                    )
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(
                        ramifications[j][0][1], ramifications[j][1][1], c_i
                    )

        moves.pop()
        ramifications.pop()

    pos = hvelu.global_L.index(L[0])  # Current element of hvelu.global_L to be required
    s_i = sign(e[pos])  # Sign of e[pos]
    c_i = (s_i + 1) // 2  # Constant-swap of T_+ and T_-

    # T_+ or T_- ?
    # Root
    ramifications[0][0][0], ramifications[0][1][0] = fp_cswap(
        ramifications[0][0][0], ramifications[0][1][0], c_i
    )
    ramifications[0][0][1], ramifications[0][1][1] = fp_cswap(
        ramifications[0][0][1], ramifications[0][1][1], c_i
    )

    if isinfinity(ramifications[0][0]) == False:

        if m[pos] > 0:
            if __formulaes != 'tvelu':
                # This branchs corresponds with the use of the new velu's formulaes

                if __verbose:
                    set_parameters_velu(sJ_list[pos], sI_list[pos], pos)

                else:
                    # -------------------------------------------------------------
                    # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                    # These paramters are required in KPs, xISOG, and xEVAL
                    if hvelu.global_L[pos] == 3:
                        b = 0
                        c = 0
                    else:
                        b = int(floor(sqrt(hvelu.global_L[pos] - 1) / 2.0))
                        c = int(floor((hvelu.global_L[pos] - 1.0) / (4.0 * b)))

                    set_parameters_velu(b, c, pos)

            KPs(ramifications[0][0], E_i, pos)
            E_i = xISOG(E_i, pos)

            b_i = (
                isequal[u[pos] == 1] ^ isequal[u[pos] == -1]
            )  # 0 if u[pos] != +1,-1; otherwise, 1 [This ask must be performed in constant-time]
            v[pos] -= 1
            u[pos] -= s_i * (
                1 + b_i
            )  # reduced by 1 unit if it was different from +1 and -1; otherwise, reduced by 2 units

    return E_i, v, u




def geometric_serie(m, l):
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


def filtered(List, sublist):
    """
    filtered()
    inputs : a list L and a sublist SL of L
    output : L \ SL
    """
    return [e for e in List if e not in sublist]

def rounds(e, n):
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
        sublists_C.append(filtered(tmp_N, sublists_L[len(sublists_L) - 1]))
        tmp_e = [(tmp_e[i] - e_min) for i in tmp_N]

    return rounds_out, sublists_L, sublists_C

def GAE(A, e, L, R, St, r, m):
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

            T_p, T_m = elligator(E_k)
            for ii in range(0, hvelu.ML.fp.exponent_of_two, 1):
                T_p = xDBL(T_p, E_k)
                T_m = xDBL(T_m, E_k)

            for l in R[j]:
                T_p = xMUL(T_p, E_k, hvelu.global_L.index(l))
                T_m = xMUL(T_m, E_k, hvelu.global_L.index(l))

            E_k, m, e = evaluate_strategy(
                E_k, list([list(T_m), list(T_p)]), L[j], St[j], len(L[j]), m, e
            )

    # Multiplicative strategy on the set of unreached small odd prime numbers
    unreached_sop = [hvelu.global_L[i] for i in range(len(hvelu.global_L)) if m[i] > 0]
    remainder_sop = [l for l in hvelu.global_L if l not in unreached_sop]

    while len(unreached_sop) > 0:

        T_p, T_m = elligator(E_k)
        for ii in range(0, hvelu.ML.fp.exponent_of_two, 1):
            T_p = xDBL(T_p, E_k)
            T_m = xDBL(T_m, E_k)

        for l in remainder_sop:
            T_p = xMUL(T_p, E_k, hvelu.global_L.index(l))
            T_m = xMUL(T_m, E_k, hvelu.global_L.index(l))

        current_n = len(unreached_sop)
        E_k, m, e = evaluate_strategy(
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
            if m[hvelu.global_L.index(unreached_sop[k])] > 0
        ]
        tmp_remainder = [
            unreached_sop[k]
            for k in range(current_n)
            if m[hvelu.global_L.index(unreached_sop[k])] == 0
        ]

        unreached_sop = list(tmp_unreached)  # Removing elements from the batch
        remainder_sop = (
            remainder_sop + tmp_remainder
        )  # Adding elements to the complement of the batch

    return E_k



def strategy_block_cost(L, e):
    """
    Next function computes the expected cost of our approach by assuming we have full torsion points
    """
    elligator_cost = numpy.array([7.0, 3.0, 10.0])  # Elligator cost
    mul_fp_by_four = (
        numpy.array([4.0, 2.0, 4.0]) * hvelu.ML.fp.exponent_of_two
    )  # Cost of computing x([2^exponent_of_two]P)

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

        bo_C = 2.0 * sum([C_xMUL[hvelu.global_L.index(L[k])] for k in tmp_Cs[j]])
        S_tmp, go_C = dynamic_programming_algorithm(
            [L[k] for k in tmp_Ls[j]], len(tmp_Ls[j])
        )

        S_out.append(S_tmp)
        C_e = (
            C_e
            + (go_C + bo_C + elligator_cost + 2.0 * mul_fp_by_four) * tmp_r[j]
        )

    return C_e, L_out, R_out, S_out, tmp_r

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
