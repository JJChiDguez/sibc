from sidh.framework import *

if setting.formulaes == 'tvelu':
    from sidh.bsidh.tvelu import *

elif setting.formulaes == 'svelu':
    from sidh.bsidh.svelu import *

else:
    from sidh.bsidh.hvelu import *


# random_key() implements an uniform random integer sample [this functions should be modified]
random_key = lambda m: random.randint(0, m)

# In order to achieve efficiency, the optimal strategies and their cost are saved in two global dictionaries (hash tables)
S = { 1: {} }    # Initialization of each strategy
C = { 1: {} }    # Initialization of the costs: 0.

for i in range(n):

    j = global_L.index(SID[i])
    S[1][tuple([global_L[j]])] = [];        # Strategy with a list with only one element (a small odd prime number l_i)
    C[1][tuple([global_L[j]])] = C_xISOG[j] # Degree-l_i isogeny construction cost

for i in range(2, n + 1):

    C[i] = { }
    S[i] = { }

'''
    dynamic_programming_algorithm():
    inputs: the list of small odd primes to be processed and its length
    output: the optimal strategy and its cost of the input list of small odd primes
'''
def dynamic_programming_algorithm(L, n):
    global S, C
    # If the approach uses dummy operations, to set DUMMY = 2.0;
    # otherwise, to set DUMMY = 1.0 (dummy free approach);

    if( len(L) != n ):

        # If the list of prime numbers doesn't have size n, then we return [],-1
        print("error:\tthe list of prime numbers has different size from %d." % n);
        return [], -1;
    else:

        # Assuming #L = n, we proceed.
        get_neighboring_sets = lambda L, k: [ tuple(L[i:i+k]) for i in range(n-k+1)] # This function computes all the k-tuple: (l_1, l_2, ..., l_{k)),
                                                                                     # (l_2, l_3, ..., l_{k+1)), ..., (l_{n-k}, l_{n-k+1, ..., l_{n)).
        for i in range(2, n+1):
 
            for Tuple in get_neighboring_sets(L, i):

                if C[i].get(Tuple) is None:

                    alpha = [ (b, 
                        C[len(Tuple[:b])][Tuple[:b]] +       # Subtriangle on the left side with b leaves
                        C[len(Tuple[b:])][Tuple[b:]] +       # Subtriangle on the right side with (i - b) leaves
                        1.0*sum([ C_xMUL[global_L.index(t)] for t in Tuple[:b] ])   +   # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with b leaves
                        1.0*sum([ C_xEVAL[global_L.index(t)] for t in Tuple[b:] ])      # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                        ) for b in range(1, i)
                        ]
                    b, C[i][Tuple] = min(alpha, key=lambda t: measure(t[1]))     # We save the minimal cost corresponding to the triangle with leaves Tuple
                    S[i][Tuple] = [b] + S[i - b][Tuple[b:]] + S[b][Tuple[:b]]    # We save the optimal strategy corresponding to the triangle with leaves Tuple

        return S[n][tuple(L)], C[n][tuple(L)]   #


'''
    evaluate_strategy():
    inputs : a projective Montgomery constants A24:= A + 2C and C24:=4C where E : y^2 = x^3 + (A/C)*x^2 + x,
             a projective Montgomery x-coordinate torsion-(l_1 x ... l_n) point x(P) := XP/ZP, the list of
             small odd primes [l_1, l_2, ..., l_n], an strategy and length of the given list of small odd 
             primes;
    output : the projective Montgomery constants a24:= a + 2c and c24:=4c where E': y^2 = x^3 + (a/c)*x^2 + x
             is E / <P>
'''
def evaluate_strategy(EVAL, S_in, T_in, ST_in, E, P, L, strategy, n):

    ramifications = []
    moves = [0]         # moves: this list determines whether an isogeny construction must be performed
    k = 0               # k: current element of the strategy

    ramifications.append(list(P))
    E_i = list(E)

    if EVAL:
        # Public points to be evaluated
        S_out = list(S_in)      # x(S) should be a torsion point with not common factors in L
        T_out = list(T_in)      # x(T) should be a torsion point with not common factors in L
        ST_out= list(ST_in)     # x(S - T) should be a torsion point with not common factors in L
    else:
        S_out = None
        T_out = None
        ST_out= None

    assert(len(strategy) == (n - 1))
    for i in range(len(strategy)):

        pos = global_L.index(L[n - 1 - i])      # Current element of global_L to be required

        # Reaching the vertex (n - 1 - i, i)
        # Vertical edges (scalar multiplications)
        prev = sum(moves)
        while prev < (n - 1 - i):

            moves.append(strategy[k])                   # Number of vertical edges to be performed
            T = list( ramifications[-1] )               # New ramification
            for j in range(prev, prev + strategy[k], 1):
                T = xMUL(T, E_i, global_L.index(L[j]))

            ramifications.append(list(T))
            prev += strategy[k]
            k += 1

        # Deciding which velu variant will be used
        if setting.formulaes != 'tvelu':
            # This branchs corresponds with the use of the new velu's formulaes

            if setting.verbose:
                set_parameters_velu(sJ_list[pos], sI_list[pos], pos)

            else:
                # -------------------------------------------------------------
                # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
                # These paramters are required in KPs, xISOG, and xEVAL
                if global_L[pos] <= 4:
                    b = 0
                    c = 0
                else:
                    b = int(floor( sqrt(global_L[pos] - 1) / 2.0) )
                    c = int(floor( (global_L[pos] - 1.0) / (4.0*b) ))

                set_parameters_velu(b, c, pos)

        # Kernel Points computation
        KPs(ramifications[-1], E_i, pos)

        # Isogeny construction
        ramifications[-1][0], E_i[0] = fp2_cswap(ramifications[-1][0], E_i[0], global_L[pos] == 4)
        ramifications[-1][1], E_i[1] = fp2_cswap(ramifications[-1][1], E_i[1], global_L[pos] == 4)
        C_i = xISOG(E_i, pos)
        ramifications[-1][0], E_i[0] = fp2_cswap(ramifications[-1][0], E_i[0], global_L[pos] == 4)
        ramifications[-1][1], E_i[1] = fp2_cswap(ramifications[-1][1], E_i[1], global_L[pos] == 4)

        # Now, we proceed by perform horizontal edges (isogeny evaluations)
        for j in range(0, len(moves) - 1, 1):

            if setting.formulaes == 'tvelu' or (setting.formulaes == 'hvelu' and global_L[pos] <= HYBRID_BOUND) or (global_L[pos] == 4):
                ramifications[j] = xEVAL(ramifications[j], pos)
            else:
                ramifications[j] = xEVAL(ramifications[j], E_i)

        if EVAL:
            # Evaluating public points
            if setting.formulaes == 'tvelu' or (setting.formulaes == 'hvelu' and global_L[pos] <= HYBRID_BOUND) or (global_L[pos] == 4):

                S_out = xEVAL(S_out, pos)
                T_out = xEVAL(T_out, pos)
                ST_out= xEVAL(ST_out,pos)
            else:
            
                S_out = xEVAL(S_out, E_i)
                T_out = xEVAL(T_out, E_i)
                ST_out= xEVAL(ST_out,E_i)

        # Updating the Montogmery curve coefficients
        E_i = [list(C_i[0]), list(C_i[1])]

        moves.pop()
        ramifications.pop()

    pos = global_L.index(L[0])                                      # Current element of global_L to be required

    if setting.formulaes != 'tvelu':
        # This branchs corresponds with the use of the new velu's formulaes

        if setting.verbose:
            set_parameters_velu(sJ_list[pos], sI_list[pos], pos)

        else:
            # -------------------------------------------------------------
            # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
            # These paramters are required in KPs, xISOG, and xEVAL
            if global_L[pos] <= 4:
                b = 0
                c = 0
            else:
                b = int(floor( sqrt(global_L[pos] - 1) / 2.0) )
                c = int(floor( (global_L[pos] - 1.0) / (4.0*b) ))

            set_parameters_velu(b, c, pos)

    # Kernel Points computations
    KPs(ramifications[0], E_i, pos)

    # Isogeny construction
    ramifications[0][0], E_i[0] = fp2_cswap(ramifications[0][0], E_i[0], global_L[pos] == 4)
    ramifications[0][1], E_i[1] = fp2_cswap(ramifications[0][1], E_i[1], global_L[pos] == 4)
    C_i = xISOG(E_i, pos)
    ramifications[0][0], E_i[0] = fp2_cswap(ramifications[0][0], E_i[0], global_L[pos] == 4)
    ramifications[0][1], E_i[1] = fp2_cswap(ramifications[0][1], E_i[1], global_L[pos] == 4)

    if EVAL:
        # Evaluating public points
        if setting.formulaes == 'tvelu' or (setting.formulaes == 'hvelu' and global_L[pos] <= HYBRID_BOUND) or (global_L[pos] == 4):

            S_out = xEVAL(S_out, pos)
            T_out = xEVAL(T_out, pos)
            ST_out= xEVAL(ST_out,pos)

        else:
            
            S_out = xEVAL(S_out, E_i)
            T_out = xEVAL(T_out, E_i)
            ST_out= xEVAL(ST_out,E_i)

    # Updating the Montogmery curve coefficients
    E_i = [list(C_i[0]), list(C_i[1])]

    return E_i, S_out, T_out, ST_out
