from sidh.bsidh.poly_redc import *

cEVAL = lambda l: numpy.array( [2.0*(l - 1.0), 2.0, (l + 1.0)] )
cISOG = lambda l: numpy.array( [(3.0*l + 2.0*hamming_weight(l) - 9.0 + isequal[l == 3]*4.0), (l + 2.0*bitlength(l) + 1.0 + isequal[l == 3]*2.0), (3.0*l - 7.0 + isequal[l == 3]*6.0)] )

C_xEVAL = list(map(cEVAL, global_L))   # list of the costs of each degree-l isogeny evaluation
C_xISOG = list(map(cISOG, global_L))   # list of the costs of each degree-l isogeny construction

# Global variables to be used in KPs, xISOG, and xEVAL

# Here, J is a set of cardinality sJ
J = None
sJ = None

# Here, ptree_I corresponds with the product tree determined by I, and I is a set of cardinality sJ
ptree_hI = None
sI = None, 

# Here, K is a set of cardinality sK
K = None
sK = None

# An extra global variable which is used in xISOG and xEVAL
XZJ4 = None

SCALED_REMAINDER_TREE = True

''' -------------------------------------------------------------------------
    yDBL()
    input : a projective Twisted Edwards y-coordinate point y(P) := YP/WP,
            and the projective Montgomery constants A24:= A + 2C and C24:=4C
            where E : y^2 = x^3 + (A/C)*x^2 + x
    output: the projective Twisted Edwards y-coordinate point y([2]P)
    ------------------------------------------------------------------------- '''
def yDBL(P, A):

    t_0 = fp2_sqr(P[0])
    t_1 = fp2_sqr(P[1])
    Z = fp2_mul(A[1], t_0);
    X = fp2_mul(Z, t_1);
    t_1 = fp2_sub(t_1, t_0);
    t_0 = fp2_mul(A[0], t_1);
    Z = fp2_add(Z, t_0);
    Z = fp2_mul(Z, t_1);
    return [fp2_sub(X,Z), fp2_add(X,Z)]

''' -------------------------------------------------------------------------
    yADD()
    input : the projective Twisted Edwards y-coordinate points y(P) := YP/WP,
            y(Q) := YQ/WQ, and y(P-Q) := YPQ/QPQ
    output: the projective Twisted Edwards y-coordinate point y(P+Q)
    ------------------------------------------------------------------------- '''
def yADD(P, Q, PQ):

    a = fp2_mul(P[1], Q[0])
    b = fp2_mul(P[0], Q[1])
    c = fp2_add(a, b)
    d = fp2_sub(a, b)
    c = fp2_sqr(c)
    d = fp2_sqr(d)

    xD = fp2_add(PQ[1], PQ[0])
    zD = fp2_sub(PQ[1], PQ[0])
    X = fp2_mul(zD, c)
    Z = fp2_mul(xD, d)
    return [fp2_sub(X,Z), fp2_add(X,Z)]

''' -------------------------------------------------------------------------
    KPs()
    input : the projective Twisted Edwards y-coordinate points y(P) := YP/WP,
            the projective Montgomery constants A24:= A + 2C and C24:=4C where 
            E : y^2 = x^3 + (A/C)*x^2 + x, and a positive integer 0 <= i < n
    output: the list of projective Twisted Edwards y-coordinate points y(P),
            y([2]P), y([3]P), ..., and y([d_i]P) where l_i = 2 * d_i + 1
    ------------------------------------------------------------------------- '''
def KPs_t(P, A, i):

    assert(isinfinity(P) == False)
    global K
    d = (global_L[i] - 1) // 2

    K = [ [[0,0], [0,0]] for j in range(d + 1) ]
    K[0] = list([fp2_sub(P[0], P[1]), fp2_add(P[0],P[1])])  # 2a
    if global_L[i] == 3:
        K[1] = list([ list(K[0]), list(K[1]) ])
    else:
        K[1] = yDBL(K[0], A)                                # 4M + 2S + 4a

    for j in range(2, d, 1):
        K[j] = yADD(K[j - 1], K[0], K[j - 2])               # 4M + 2S + 6a

    return K                                                # 2(l - 3)M + (l - 3)S + 3(l - 3)a

''' ------------------------------------------------------------------------------
    xISOG()
    input : the projective Montgomery constants A24:= A + 2C and C24:=4C where
            E : y^2 = x^3 + (A/C)*x^2 + x, the list of projective Twisted 
            Edwards y-coordinate points y(P), y([2]P), y([3]P), ..., and y([d_i]P)
            where l_i = 2 * d_i + 1, and a positive integer 0 <= i < n
    output: the projective Montgomery constants a24:= a + 2c and c24:=4c where
            E': y^2 = x^3 + (a/c)*x^2 + x is a degree-(l_i) isogenous curve to E
    ------------------------------------------------------------------------------ '''
def xISOG_t(A, i):

    global K
    l = global_L[i]             # l
    bits_of_l = bitlength(l)    # Number of bits of L[i]
    d = (l - 1) // 2            # Here, l = 2d + 1

    By = K[0][0]
    Bz = K[0][1]
    for j in range(1, d, 1):

        By = fp2_mul(By, K[j][0])
        Bz = fp2_mul(Bz, K[j][1])

    bits_of_l -= 1
    constant_d_edwards = fp2_sub(A[0], A[1])     # 1a

    tmp_a = A[0]
    tmp_d = constant_d_edwards
    # left-to-right method for computing a^l and d^l
    for j in range(1, bits_of_l + 1):

        tmp_a = fp2_sqr(tmp_a)
        tmp_d = fp2_sqr(tmp_d)
        if( ( (l >> (bits_of_l - j)) & 1 ) != 0 ):

            tmp_a = fp2_mul(tmp_a, A[0])
            tmp_d = fp2_mul(tmp_d, constant_d_edwards)

    for j in range(3):

        By = fp2_sqr(By)
        Bz = fp2_sqr(Bz)

    C0 = fp2_mul(tmp_a, Bz)
    C1 = fp2_mul(tmp_d, By)
    C1 = fp2_sub(C0, C1)

    return [C0, C1]                             # (l - 1 + 2*HW(l) - 2)M + 2(|l|_2 + 1)S + 2a

''' ------------------------------------------------------------------------------
    xEVAL()
    input : the projective Montgomery constants A24:= A + 2C and C24:=4C where
            E : y^2 = x^3 + (A/C)*x^2 + x, the list of projective Twisted 
            Edwards y-coordinate points y(P), y([2]P), y([3]P), ..., and y([d_i]P)
            where l_i = 2 * d_i + 1, and a positive integer 0 <= i < n
    output: the projective Montgomery constants a24:= a + 2c and c24:=4c where
            E': y^2 = x^3 + (a/c)*x^2 + x is a degree-(l_i) isogenous curve to E
    ------------------------------------------------------------------------------ '''
def xEVAL_t(P, i):

    global K
    d = (global_L[i] - 1) // 2                  # Here, l = 2d + 1

    Q0 = fp2_add(P[0], P[1])
    Q1 = fp2_sub(P[0], P[1])
    R0, R1 = CrissCross(K[0][1], K[0][0], Q0, Q1)
    for j in range(1, d, 1):

        T0, T1 = CrissCross(K[j][1], K[j][0], Q0, Q1)
        R0 = fp2_mul(T0, R0)
        R1 = fp2_mul(T1, R1)

    R0 = fp2_sqr(R0)
    R1 = fp2_sqr(R1)
    X = fp2_mul(P[0], R0)
    Z = fp2_mul(P[1], R1)

    return [X, Z]                               # 2(l - 1)M + 2S + (l + 1)a

# Next functions is used for setting the cardinalities sI, sJ, and sK
def set_parameters_velu(b, c, i):


    global sJ
    global sI
    global sK

    assert(b <= c)
    
    # At this step, everythin is correct
    sJ = b
    sI = c
    d = ((global_L[i] - 2 - 4*b*c - 1) // 2) + 1
    assert(d >= 0)
    sK = d
    return None

def print_parameters_velu():

    print("| sI: %3d, sJ: %3d, sK: %3d |" % (sI, sJ, sK), end="")
    return None

# KPs computes x([i]P), x([j]P), and x([k]P) for each i in I, j in J, and j in J.
# I, J, and K are defined according to Examples 4.7 and 4.12 of https://eprint.iacr.org/2020/341
def KPs_s(P, A, i):

    # Global variables to be used
    global J
    global sJ
    global ptree_hI
    global sI
    global K
    global sK

    # This functions computes all the independent data from the input of algorithm 2 of https://eprint.iacr.org/2020/341
    if sI == 0:

        # Case global_L[i] = 3 is super very special case (nothing to do)
        J = []
        ptree_hI = None
        K = [list(P)]
        #J, b, ptree_hI, c, K, d
        assert(sJ == 0 and sI == 0 and sK == 1)
        return None

    # We need to ensure sI is greater than or equal sJ
    assert(sI >= sJ)
    
    # Additionally, sK should be greater than or equal to zero. If that is not the case, then sJ and sI were badly chosen
    assert(sK >= 0)


    if sI == 1:
        # Branch corresponds with global_L[i] = 5 and global_L[i] = 7
        # Recall sJ > 0, then sJ = 1
        assert(sJ == 1)
        P2 = xDBL(P, A)

        J = [ list(P) ]

        I = [ list(P2) ]
        hI = [ list([fp2_sub([0,0], P2[0]), P2[1]]) ]            # we only need to negate x-coordinate of each point
        ptree_hI = product_tree(hI, sI)                      # product tree of hI

        if not SCALED_REMAINDER_TREE:
            # Using remainder trees
            ptree_hI = reciprocal_tree({'rpoly':[[1,0]], 'rdeg':0, 'fpoly':[[1,0]], 'fdeg':0, 'a':[1,0]}, 2*sJ + 1, ptree_hI, sI)   # reciprocal of the root is used to compute sons' reciprocals

        else:
            # Using scaled remainder trees
            assert( (2*sJ - sI + 1) > sI )
            ptree_hI['reciprocal'], ptree_hI['a'] = reciprocal(ptree_hI['poly'][::-1], sI + 1, 2*sJ - sI + 1)
            ptree_hI['scaled'], ptree_hI['as'] = list(ptree_hI['reciprocal'][:sI]), ptree_hI['a']

        # Notice, 0 <= sK <= 1
        assert(sK <= 1)
        if sK == 1:
            K = [ list(P2) ]
        else:
            K = []

        return None

    # At this step, sI > 1
    assert(sI > 1)
    if sJ == 1:
        # This branch corresponds with global_L[i] = 11 and global_L[i] = 13
        Q = xDBL(P, A)                              # x([2]P)
        Q2= xDBL(Q, A)                              # x([2]Q)

        J = [ list(P) ]

        I = [ [[0,0],[0,0]] ] * sI
        I[0] = list(Q)                              # x(   Q)
        I[1] = xADD(Q2, I[0], I[0])                 # x([3]Q)
        for ii in range(2, sI, 1):
            I[ii] = xADD(I[ii - 1], Q2, I[ii -2])   # x([2**i + 1]Q)

        hI = [ [fp2_sub([0,0], iP[0]), iP[1]] for iP in I ]      # we only need to negate x-coordinate of each point
        ptree_hI = product_tree(hI, sI)                      # product tree of hI

        if not SCALED_REMAINDER_TREE:
            # Using remainder trees
            ptree_hI = reciprocal_tree({'rpoly':[[1,0]], 'rdeg':0, 'fpoly':[[1,0]], 'fdeg':0, 'a':[1,0]}, 2*sJ + 1, ptree_hI, sI)   # reciprocal of the root is used to compute sons' reciprocals

        else:
            # Using scaled remainder trees
            assert( (2*sJ - sI + 1) <= sI )
            ptree_hI['scaled'], ptree_hI['as'] = reciprocal(ptree_hI['poly'][::-1], sI + 1, sI)
            ptree_hI['reciprocal'], ptree_hI['a'] = list(ptree_hI['scaled'][:(2*sJ - sI + 1)]), ptree_hI['as']

        # Notice, 0 <= sK <= 1
        assert(sK <= 1)
        if sK == 1:
            K = [ list(Q) ]
        else:
            K = []

        return None

    # Now, we ensure sI >= sJ > 1
    assert(sJ > 1)

    # In other words, we can proceed by the general case

    # ------------------------------------------------
    # Computing [j]P for each j in {1, 3, ..., 2*sJ - 1}
    J = [ [[0,0],[0,0]] ] * sJ 
    J[0]= list(P)                               # x(   P)
    P2  = xDBL(P, A)                            # x([2]P)
    J[1] = xADD(P2, J[0], J[0])                 # x([3]P)
    for jj in range(2, sJ, 1):
        J[jj] = xADD(J[jj - 1], P2, J[jj - 2])  # x([2*jj + 1]P)

    # -------------------------------------------------------
    # Computing [i]P for i in { (2*sJ) * (2i + 1) : 0 <= i < sI}
    bhalf_floor= sJ // 2
    bhalf_ceil = sJ - bhalf_floor
    P4 = xDBL(P2, A)                                    # x([4]P)
    P2[0], P4[0] = fp2_cswap(P2[0], P4[0], sJ % 2)        # Constant-time swap
    P2[1], P4[1] = fp2_cswap(P2[1], P4[1], sJ % 2)        # x([4]P) <--- coditional swap ---> x([2]P)
    Q = xADD(J[bhalf_ceil], J[bhalf_floor - 1], P2)     # Q := [2b]P
    P2[0], P4[0] = fp2_cswap(P2[0], P4[0], sJ % 2)        # Constant-time swap
    P2[1], P4[1] = fp2_cswap(P2[1], P4[1], sJ % 2)        # x([4]P) <--- coditional swap ---> x([2]P)

    I = [ [0,0] ] * sI
    I[0] = list(Q)                              # x(   Q)
    Q2 = xDBL(Q,A)                              # x([2]Q)
    I[1] = xADD(Q2, I[0], I[0])                 # x([3]Q)
    for ii in range(2, sI, 1):
        I[ii] = xADD(I[ii - 1], Q2, I[ii -2])   # x([2**i + 1]Q)

    # --------------------------------------------------------------
    # Computing [k]P for k in { 4*sJ*sI + 1, ..., l - 6, l - 4, l - 2}
    K = [ [[0,0], [0,0]] ] * sK

    if sK >= 1:
        K[0] = list(P2)                             # x([l - 2]P) = x([-2]P) = x([2]P)
   
    if sK >= 2:
        K[1] = list(P4)                             # x([l - 4]P) = x([-4]P) = x([4]P)

    for k in range(2, sK, 1):
        K[k] = xADD(K[k - 1], P2, K[k - 2])

    # ------------------------------------------------------------------------------------------------------
    #                   ~~~~~~~~               ~~~~~~~~
    #                    |    |                 |    |
    # Computing h_I(W) = |    | (W - x([i]P)) = |    | (Zi * W - Xi) / Zi where x([i]P) = Xi/Zi
    #                    i in I                 i in I
    # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

    hI = [ [fp2_sub([0,0], iP[0]), iP[1]] for iP in I ]      # we only need to negate x-coordinate of each point
    ptree_hI = product_tree(hI, sI)                      # product tree of hI

    if not SCALED_REMAINDER_TREE: 
        # Using scaled remainder trees
        ptree_hI = reciprocal_tree({'rpoly':[[1,0]], 'rdeg':0, 'fpoly':[[1,0]], 'fdeg':0, 'a':[1,0]}, 2*sJ + 1, ptree_hI, sI)   # reciprocal of the root is used to compute sons' reciprocals

    else:
        # Using scaled remainder trees
        if sI < (2*sJ - sI + 1):
            ptree_hI['reciprocal'], ptree_hI['a'] = reciprocal(ptree_hI['poly'][::-1], sI + 1, 2*sJ - sI + 1)
            ptree_hI['scaled'], ptree_hI['as'] = list(ptree_hI['reciprocal'][:sI]), ptree_hI['a']

        else:
            ptree_hI['scaled'], ptree_hI['as'] = reciprocal(ptree_hI['poly'][::-1], sI + 1, sI)
            ptree_hI['reciprocal'], ptree_hI['a'] = list(ptree_hI['scaled'][:(2*sJ - sI + 1)]), ptree_hI['as']

    # Polynomial D_j of algorithm 2 from https://eprint.iacr.org/2020/341 is not required, 
    # but we need some some squares and products determined by list J
    # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

    # Ensuring the cardinality of each ser coincide with the expected one
    assert(len(I) == sI)
    assert(len(J) == sJ)
    assert(len(K) == sK)

    return None

# Next function perform algorithm 2 of https://eprint.iacr.org/2020/341 with input \alpha = 1 and \alpha = -1, and
# then it computes the isogenous Montgomery curve coefficient
def xISOG_s(A, i):

    global J
    global sJ
    global ptree_hI
    global sI
    global K
    global sK
    global XZJ4

    AA = fp2_add(A[0], A[0])     # 2A' + 4C
    AA = fp2_sub(AA, A[1])       # 2A'
    AA = fp2_add(AA, AA)         # 4A'

    # Polynomial D_j of algorithm 2 from https://eprint.iacr.org/2020/341 is not required, 
    # but we need some some squares and products determined by list J
    # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

    SUB_SQUARED = [ 0 for j in range(0, sJ, 1) ]             # 
    ADD_SQUARED = [ 0 for j in range(0, sJ, 1) ]             # 

    # List XZJ4 is required for degree-l isogeny evaluations...
    XZJ4 = [ 0 for j in range(0, sJ, 1) ]                         # 2*Xj*Zj
    for j in range(0, sJ, 1):

        SUB_SQUARED[j] = fp2_sub(J[j][0], J[j][1])               # (Xj - Zj)
        SUB_SQUARED[j] = fp2_sqr(SUB_SQUARED[j])                 # (Xj - Zj)^2
        
        ADD_SQUARED[j] = fp2_add(J[j][0], J[j][1])               # (Xj + Zj)
        ADD_SQUARED[j] = fp2_sqr(ADD_SQUARED[j])                 # (Xj + Zj)^2

        XZJ4[j] = fp2_sub(SUB_SQUARED[j], ADD_SQUARED[j])        # -4*Xj*Zj

    # --------------------------------------------------------------------------------------------------
    #                   ~~~~~~~~
    #                    |    | 
    # Computing E_J(W) = |    | [ F0(W, x([j]P)) * alpha^2 + F1(W, x([j]P)) * alpha + F2(W, x([j]P)) ]
    #                    j in J 
    # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates
    # In particular, for a degree-l isogeny construction, we need alpha = 1 and alpha = -1

    # EJ_0 is the one determined by alpha = 1
    EJ_0 = [ [[0,0],[0,0],[0,0]] for j in range(0, sJ, 1) ]
    # EJ_1 is the one determined by alpha = -1
    EJ_1 = [ [[0,0],[0,0],[0,0]] for j in range(0, sJ, 1) ]

    for j in range(0, sJ, 1):

        # However, each SUB_SQUARED[j] and ADD_SQUARED[j] should be multiplied by C
        tadd = fp2_mul(ADD_SQUARED[j], A[1])
        tsub = fp2_mul(SUB_SQUARED[j], A[1])

        # We require the double of tadd and tsub
        tadd2= fp2_add(tadd, tadd)
        tsub2= fp2_add(tsub, tsub)

        t1 = fp2_mul(XZJ4[j], AA)                    #       A *(-4*Xj*Zj)

        # Case alpha = 1
        linear = fp2_sub(t1, tadd2)                  #       A *(-4*Xj*Zj)  - C * (2 * (Xj + Zj)^2)
        EJ_0[j] = [tsub, linear, tsub]

        # Case alpha = -1
        linear = fp2_sub(tsub2, t1)                  #       C * (2 * (Xj - Zj)^2) - A *(-4*Xj*Zj)
        EJ_1[j] = [tadd, linear, tadd] 

    # The faster way for multiplying is using a divide-and-conquer approach
    poly_EJ_0 = product_selfreciprocal_tree(EJ_0, sJ)['poly']       # product tree of EJ_0 (we only require the root)
    poly_EJ_1 = product_selfreciprocal_tree(EJ_1, sJ)['poly']       # product tree of EJ_1 (we only require the root)
   
    if not SCALED_REMAINDER_TREE:
        # Remainder tree computation
        remainders_EJ_0 = multieval_unscaled(poly_EJ_0, 2*sJ + 1,  ptree_hI, sI)
        remainders_EJ_1 = multieval_unscaled(poly_EJ_1, 2*sJ + 1,  ptree_hI, sI)

    else:
        # Approach using scaled remainder trees
        if ptree_hI != None:
            poly_EJ_0 = poly_redc(poly_EJ_0, 2*sJ + 1, ptree_hI)
            fg_0 = poly_mul_middle(ptree_hI['scaled'], sI, poly_EJ_0[::-1], sI)
            remainders_EJ_0 = multieval_scaled(fg_0[::-1], sI, [[1,0]] + [[0,0]]*(sI - 1), sI, ptree_hI, sI)
    
            poly_EJ_1 = poly_redc(poly_EJ_1, 2*sJ + 1, ptree_hI)
            fg_1 = poly_mul_middle(ptree_hI['scaled'], sI, poly_EJ_1[::-1], sI)
            remainders_EJ_1 = multieval_scaled(fg_1[::-1], sI, [[1,0]] + [[0,0]]*(sI - 1), sI, ptree_hI, sI)
        else:
            remainders_EJ_0 = []
            remainders_EJ_1 = []

    # Multipying all the remainders
    r0 = product(remainders_EJ_0, sI)
    r1 = product(remainders_EJ_1, sI)
    
    # ---------------------------------------------------------------------------------
    # Now, we proceed by computing the missing part which is determined by K
    # Notice, the denominators are the same and then they annulled between them
    # In other words, it is not required to compute the product of all Zk's with k In K

    # Case alpha = 1
    hK_0 = [ [fp2_sub(K[k][1], K[k][0])] for k in range(0, sK, 1 )]
    hK_0 = product(hK_0, sK)                     # product of (Zk - Xk) for each k in K
    # Case alpha = -1
    hK_1 = [ [fp2_add(K[k][1], K[k][0])] for k in range(0, sK, 1 )]
    hK_1 = product(hK_1, sK)                     # product of (Zk + Xk) for each k in K 
    
    # --------------------------------------------------------------
    # Now, we have all the ingredients for computing the image curve
    A24m = fp2_sub(A[0], A[1])                    # A' - 2C

    A24 = fp2_exp(A[0], global_L[i])                    # (A' + 2C)^l
    A24m= fp2_exp(A24m, global_L[i])                    # (A' - 2C)^l

    t24m= fp2_mul(hK_1, r1)                      # output of algorithm 2 with alpha =-1 and without the demoninator
    t24m= fp2_sqr(t24m)                          # raised at 2
    t24m= fp2_sqr(t24m)                          # raised at 4
    t24m= fp2_sqr(t24m)                          # raised at 8

    t24 = fp2_mul(hK_0, r0)                      # output of algorithm 2 with alpha = 1 and without the demoninator 
    t24 = fp2_sqr(t24)                           # raised at 2
    t24 = fp2_sqr(t24)                           # raised at 4
    t24 = fp2_sqr(t24)                           # raised at 8

    A24 = fp2_mul(A24, t24m)
    A24m= fp2_mul(A24m, t24)

    # Now, we have d = (A24m / A24) where the image Montgomery cuve coefficient is
    #      B'   2*(1 + d)   2*(A24 + A24m)
    # B = ---- = --------- = --------------
    #      C      (1 - d)     (A24 - A24m)
    # However, we required B' + 2C = 4*A24 and 4C = 4 * (A24 - A24m)

    t24m = fp2_sub(A24, A24m)                    #   (A24 - A24m)
    t24m = fp2_add(t24m, t24m)                   # 2*(A24 - A24m)
    t24m = fp2_add(t24m, t24m)                   # 4*(A24 - A24m)

    t24 = fp2_add(A24, A24)                      # 2 * A24
    t24 = fp2_add(t24, t24)                      # 4 * A24
    
    #return [t24, t24m], ptree_hI, XZJ4
    return [t24, t24m]

def xEVAL_s(P, A):

    AA = fp2_add(A[0], A[0])     # 2A' + 4C
    AA = fp2_sub(AA, A[1])       # 2A'
    AA = fp2_add(AA, AA)         # 4A'

    # --------------------------------------------------------------------------------------------------
    #                   ~~~~~~~~
    #                    |    | 
    # Computing E_J(W) = |    | [ F0(W, x([j]P)) * alpha^2 + F1(W, x([j]P)) * alpha + F2(W, x([j]P)) ]
    #                    j in J 
    # In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates
    # In particular, for a degree-l isogeny construction, we need alpha = 1 and alpha = -1

    # EJ_0 is the one determined by alpha = x
    EJ_0 = [ [[0,0],[0,0],[0,0]] for j in range(0, sJ, 1) ]
    # Notice, the corresponding EJ_1 that is determined by alpha = 1/x can be computed by using EJ_0

    XZ_add = fp2_add(P[0], P[1])                 # X + Z
    XZ_sub = fp2_sub(P[0], P[1])                 # X - Z
    
    AXZ2 = fp2_mul(P[0], P[1])                   # X * Z
    t1 = fp2_sqr(P[0])                           # X^2
    t2 = fp2_sqr(P[1])                           # Z^2

    CX2Z2 = fp2_add(t1, t2)                      #      X^2 + Z^2
    CX2Z2 = fp2_mul(CX2Z2, A[1])                 # C * (X^2 + Z^2)

    AXZ2 = fp2_add(AXZ2, AXZ2)                   #       2 * (X * Z)
    CXZ2 = fp2_mul(AXZ2, A[1])                   # C  * [2 * (X * Z)]
    AXZ2 = fp2_mul(AXZ2, AA)                     # A' * [2 * (X * Z)]

    for j in range(0, sJ, 1):

        XZj_add = fp2_add(J[j][0], J[j][1])      # Xj + Zj
        XZj_sub = fp2_sub(J[j][0], J[j][1])      # Xj - Zj

        t1 = fp2_mul(XZ_sub, XZj_add)            # (X - Z) * (Xj + Zj)
        t2 = fp2_mul(XZ_add, XZj_sub)            # (X + Z) * (Xj - Zj)

        # Computing the quadratic coefficient
        quadratic = fp2_sub(t1, t2)              #   2 * [(X*Zj) - (Z*Xj)]
        quadratic = fp2_sqr(quadratic)           # ( 2 * [(X*Zj) - (Z*Xj)] )^2
        quadratic = fp2_mul(A[1], quadratic)     # C * ( 2 * [(X*Zj) - (Z*Xj)] )^2

        # Computing the constant coefficient
        constant = fp2_add(t1, t2)               #   2 * [(X*Xj) - (Z*Zj)]
        constant = fp2_sqr(constant)             # ( 2 * [(X*Xj) - (Z*Zj)] )^2
        constant = fp2_mul(A[1], constant)       # C * ( 2 * [(X*Xj) - (Z*Zj)] )^2

        # Computing the linear coefficient
        # ----------------------------------------------------------------------------------------------------------
        # C * [ (-2*Xj*Zj)*(alpha^2 + 1) + (-2*alpha)*(Xj^2 + Zj^2)] + [A' * (-2*Xj*Zj) * (2*X*Z)] where alpha = X/Z
        t1 = fp2_add(J[j][0], J[j][1])           #     (Xj + Zj)
        t1 = fp2_sqr(t1)                         #     (Xj + Zj)^2
        t1 = fp2_add(t1, t1)                     # 2 * (Xj + Zj)^2
        t1 = fp2_add(t1, XZJ4[j])                # 2 * (Xj + Zj)^2 - (4*Xj*Zj) := 2 * (Xj^2 + Zj^2)
        t1 = fp2_mul(t1, CXZ2)                   # [2 * (Xj^2 + Zj^2)] * (2 * [ C * (X * Z)])

        t2 = fp2_mul(CX2Z2, XZJ4[j])             # [C * (X^2 + Z^2)] * (-4 * Xj * Zj)
        t1 = fp2_sub(t2, t1)                     # [C * (X^2 + Z^2)] * (-4 * Xj * Zj) - [2 * (Xj^2 + Zj^2)] * (2 * [ C * (X * Z)])

        t2 = fp2_mul(AXZ2, XZJ4[j])              # (2 * [A' * (X * Z)]) * (-4 * Xj * Zj)
        linear = fp2_add(t1, t2)                 # This is our desired equation but multiplied by 2
        linear = fp2_add(linear, linear)         # This is our desired equation but multiplied by 4
        # ----------------------------------------------------------------------------------------------------------

        # Case alpha = X / Z
        EJ_0[j] = [constant, linear, quadratic] 

    # The faster way for multiplying is using a divide-and-conquer approach
    poly_EJ_0 = product_tree(EJ_0, sJ)['poly']       # product tree of EJ_0 (we only require the root)
    poly_EJ_1 = list(poly_EJ_0[::-1])               # product tree of EJ_1(x) = x^{2b + 1} EJ_0(1/X)

    if not SCALED_REMAINDER_TREE:
        # Remainder tree computation
        remainders_EJ_0 = multieval_unscaled(poly_EJ_0, 2*sJ + 1,  ptree_hI, sI)
        remainders_EJ_1 = multieval_unscaled(poly_EJ_1, 2*sJ + 1,  ptree_hI, sI)

    else:
        # Approach using scaled remainder trees
        if ptree_hI != None:
            poly_EJ_0 = poly_redc(poly_EJ_0, 2*sJ + 1, ptree_hI)
            fg_0 = poly_mul_middle(ptree_hI['scaled'], sI, poly_EJ_0[::-1], sI)
            remainders_EJ_0 = multieval_scaled(fg_0[::-1], sI, [[1,0]] + [[0,0]]*(sI - 1), sI, ptree_hI, sI)
    
            poly_EJ_1 = poly_redc(poly_EJ_1, 2*sJ + 1, ptree_hI)
            fg_1 = poly_mul_middle(ptree_hI['scaled'], sI, poly_EJ_1[::-1], sI)
            remainders_EJ_1 = multieval_scaled(fg_1[::-1], sI, [[1,0]] + [[0,0]]*(sI - 1), sI, ptree_hI, sI)
        else:
            remainders_EJ_0 = []
            remainders_EJ_1 = []

    # Multipying all the remainders
    r0 = product(remainders_EJ_0, sI)
    r1 = product(remainders_EJ_1, sI)
    
    # ---------------------------------------------------------------------------------
    # Now, we proceed by computing the missing part which is determined by K
    # Notice, the denominators are the same and then they annulled between them
    # In other words, it is not required to compute the product of all Zk's with k In K

    hK_0 = [ [[0,0]] ] * sK
    hK_1 = [ [[0,0]] ] * sK
    for k in range(0, sK, 1):

        XZk_add = fp2_add(K[k][0], K[k][1])      # Xk + Zk
        XZk_sub = fp2_sub(K[k][0], K[k][1])      # Xk - Zk
        t1 = fp2_mul(XZ_sub, XZk_add)            # (X - Z) * (Xk + Zk)
        t2 = fp2_mul(XZ_add, XZk_sub)            # (X + Z) * (Xk - Zk)

        # Case alpha = X/Z
        hK_0[k] = [fp2_sub(t1, t2)]              # 2 * [(X*Zk) - (Z*Xk)]

        # Case 1/alpha = Z/X
        hK_1[k] = [fp2_add(t1, t2)]              # 2 * [(X*Xk) - (Z*Zk)]

    hK_0 = product(hK_0, sK)                     # product of (XZk - ZXk) for each k in K
    hK_1 = product(hK_1, sK)                     # product of (XXk - ZZk) for each k in K

    # ---------------------------------------------------------------------------------
    # Now, unifying all the computations
    XX = fp2_mul(hK_1, r1)                       # output of algorithm 2 with 1/alpha = Z/X and without the demoninator
    XX = fp2_sqr(XX)
    XX = fp2_mul(XX, P[0])

    ZZ = fp2_mul(hK_0, r0)                       # output of algorithm 2 with alpha = X/Z and without the demoninator
    ZZ = fp2_sqr(ZZ)
    ZZ = fp2_mul(ZZ, P[1])

    return [XX, ZZ]

# Degree-4 isogeny construction
def xISOG_4(P):

    global K
    K = [ [0,0], [0,0], [0,0] ]

    K[1] = fp2_sub(P[0], P[1])
    K[2] = fp2_add(P[0], P[1])
    K[0] = fp2_sqr(P[1])
    K[0] = fp2_add(K[0], K[0])

    C24 = fp2_sqr(K[0])
    K[0]= fp2_add(K[0], K[0])
    A24 = fp2_sqr(P[0])
    A24 = fp2_add(A24, A24)
    A24 = fp2_sqr(A24)
    return [A24, C24]

# Degree-4 isogeny evaluation
def xEVAL_4(Q):

    t0 = fp2_add(Q[0], Q[1])
    t1 = fp2_sub(Q[0], Q[1])
    XQ = fp2_mul(t0, K[1])
    ZQ = fp2_mul(t1, K[2])
    t0 = fp2_mul(t0, t1)
    t0 = fp2_mul(t0, K[0])
    t1 = fp2_add(XQ, ZQ)
    ZQ = fp2_sub(XQ, ZQ)
    t1 = fp2_sqr(t1)
    ZQ = fp2_sqr(ZQ)
    XQ = fp2_add(t0, t1)
    t0 = fp2_sub(ZQ, t0)
    XQ = fp2_mul(XQ, t1)
    ZQ = fp2_mul(ZQ, t0)

    return [XQ, ZQ]

# Tradicional velu formulas will be used for l_i <= 101
HYBRID_BOUND = 83

def KPs(P, A, i):

    if global_L[i] <= HYBRID_BOUND and (global_L[i] != 4):
        return KPs_t(P, A, i)
    else:
        return KPs_s(P, A, i)

def xISOG(A, i):

    if global_L[i] <= HYBRID_BOUND:

        if global_L[i] != 4:
            return xISOG_t(A, i)
        else:
            # A should corresponds with an order-4 point
            return xISOG_4(A)

    else:
        return xISOG_s(A, i)

def xEVAL(P, v):

    if type(v) == int:

        if global_L[v] != 4:
            return xEVAL_t(P, v)
        else:
            return xEVAL_4(P)

    else:
        return xEVAL_s(P, v)

# Get cost of the isogeny constructions and evaluations
sI_list = None
sJ_list = None
def cISOG_and_cEVAL():

    global C_xISOG
    global C_xEVAL

    global sI_list
    global sJ_list

    if setting.verbose:

        sI_list = []
        sJ_list = []
        f = open(ijk_data + setting.prime)

        for i in range(0, np + nm, 1):

            bc = f.readline()
            bc = [ int(bci) for bci in bc.split() ]
            sJ_list.append(bc[0])
            sI_list.append(bc[1])

        f.close()

    # E[p + 1]
    # First, we look for a full torsion point
    A = [ [0x8, 0x0], [0x4, 0x0] ]

    # Reading public generators points
    f = open(gen_data + setting.prime)
    # x(PA), x(QA) and x(PA - QA)
    PQA = f.readline()
    PQA = [ int(x, 16) for x in PQA.split() ]
    PA  = [list(PQA[0:2]), [0x1,0x0]]
    QA  = [list(PQA[2:4]), [0x1,0x0]]
    PQA = [list(PQA[4:6]), [0x1,0x0]]

    # x(PB), x(QB) and x(PB - QB)
    PQB = f.readline()
    PQB = [ int(x, 16) for x in PQB.split() ]
    PB  = [list(PQB[0:2]), [0x1,0x0]]
    QB  = [list(PQB[2:4]), [0x1,0x0]]
    PQB = [list(PQB[4:6]), [0x1,0x0]]

    f.close()

    for i in range(0, Ep[0] - 1, 1):

        PA = xMUL(PA, A, 0)
        QA = xMUL(QA, A, 0)
        PQA= xMUL(PQA, A, 0)

    # Random kernels for counting the 
    T_p = Ladder3pt(random.randint(0, p - 1), PA, QA, PQA, A)
    T_m = Ladder3pt(random.randint(0, p - 1), PB, QB, PQB, A)

    for i in range(0, np, 1):
        for j in range(0, Ep[i] - 1, 1):
            T_p = xMUL(T_p, A, i)

    for i in range(np, np + nm, 1):
        for j in range(0, Em[i - np] - 1, 1):
            T_m = xMUL(T_m, A, i)

    for i in range(0, np, 1):

        if setting.verbose:
            set_parameters_velu(sJ_list[i], sI_list[i], i)

        else:
            # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
            # These paramters are required in KPs, xISOG, and xEVAL
            if global_L[i] <= 4:
                b = 0
                c = 0
            else:
                b = int(floor( sqrt(global_L[i] - 1) / 2.0) )
                c = int(floor( (global_L[i] - 1.0) / (4.0*b) ))

            set_parameters_velu(b, c, i)

        # Getting an orderl-l point
        Tp = list(T_p)
        for j in range(i + 1, np, 1):
            Tp = xMUL(Tp, A, j)

        # Cost of xISOG() and KPs()
        set_zero_ops()
        KPs(Tp, A, i)
        t = get_ops()
        C_xISOG[i] = numpy.array( [ t[0] * 1.0, t[1] * 1.0, t[2] * 1.0] )

        set_zero_ops()
        Tp[0], A[0] = fp2_cswap(Tp[0], A[0], global_L[i] == 4)
        Tp[1], A[1] = fp2_cswap(Tp[1], A[1], global_L[i] == 4)
        B = xISOG(A, i)
        Tp[0], A[0] = fp2_cswap(Tp[0], A[0], global_L[i] == 4)
        Tp[1], A[1] = fp2_cswap(Tp[1], A[1], global_L[i] == 4)
        t = get_ops()
        C_xISOG[i] += numpy.array( [ t[0] * 1.0, t[1] * 1.0, t[2] * 1.0] )

        # xEVAL: kernel point determined by the next isogeny evaluation
        set_zero_ops()
        if global_L[i] <= HYBRID_BOUND:
            T_p = xEVAL(T_p, i)
        else:
            T_p = xEVAL(T_p, A)

        # Cost of xEVAL
        set_zero_ops()
        if global_L[i] <= HYBRID_BOUND:
            T_m = xEVAL(T_m, i)
        else:
            T_m = xEVAL(T_m, A)

        t = get_ops()
        C_xEVAL[i] = numpy.array( [ t[0] * 1.0, t[1] * 1.0, t[2] * 1.0] )

        # Updating the new next curve
        A = list(B)

    # E[p - 1]
    # First, we look for a full torsion point
    A = [ [0x8, 0x0], [0x4, 0x0] ]
    T_m = Ladder3pt(random.randint(0, p - 1), PA, QA, PQA, A)
    T_p = Ladder3pt(random.randint(0, p - 1), PB, QB, PQB, A)

    for i in range(np, np + nm, 1):

        if setting.verbose:
            set_parameters_velu(sJ_list[i], sI_list[i], i)

        else:
            # Parameters sJ and sI correspond with the parameters b and b' from example 4.12 of https://eprint.iacr.org/2020/341
            # These paramters are required in KPs, xISOG, and xEVAL
            if global_L[i] == 3:
                b = 0
                c = 0
            else:
                b = int(floor( sqrt(global_L[i] - 1) / 2.0) )
                c = int(floor( (global_L[i] - 1.0) / (4.0*b) ))

            set_parameters_velu(b, c, i)

        # Getting an orderl-l point
        Tp = list(T_p)
        for j in range(i + 1, np + nm, 1):
            Tp = xMUL(Tp, A, j)

        # Cost of xISOG() and KPs()
        set_zero_ops()
        KPs(Tp, A, i)
        t = get_ops()
        C_xISOG[i] = numpy.array( [ t[0] * 1.0, t[1] * 1.0, t[2] * 1.0] )

        set_zero_ops()
        B = xISOG(A, i)
        t = get_ops()
        C_xISOG[i] += numpy.array( [ t[0] * 1.0, t[1] * 1.0, t[2] * 1.0] )

        # xEVAL: kernel point determined by the next isogeny evaluation
        set_zero_ops()
        if global_L[i] <= HYBRID_BOUND:
            T_p = xEVAL(T_p, i)
        else:
            T_p = xEVAL(T_p, A)

        # Cost of xEVAL
        set_zero_ops()
        if global_L[i] <= HYBRID_BOUND:
            T_m = xEVAL(T_m, i)
        else:
            T_m = xEVAL(T_m, A)

        t = get_ops()
        C_xEVAL[i] = numpy.array( [ t[0] * 1.0, t[1] * 1.0, t[2] * 1.0] )

        # Updating the new next curve
        A = list(B)
    set_zero_ops()
    return None

# Now, we proceed to store all the correct costs
cISOG_and_cEVAL()
