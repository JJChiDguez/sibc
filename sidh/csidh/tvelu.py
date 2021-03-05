from sidh.csidh.montgomery import *

cEVAL = lambda l: numpy.array([2.0 * (l - 1.0), 2.0, (l + 1.0)])
cISOG = lambda l: numpy.array(
    [
        (3.0 * l + 2.0 * hamming_weight(l) - 9.0 + isequal[l == 3] * 4.0),
        (l + 2.0 * bitlength(l) + 1.0 + isequal[l == 3] * 2.0),
        (3.0 * l - 7.0 + isequal[l == 3] * 6.0),
    ]
)

C_xEVAL = list(
    map(cEVAL, global_L)
)  # list of the costs of each degree-l isogeny evaluation
C_xISOG = list(
    map(cISOG, global_L)
)  # list of the costs of each degree-l isogeny construction

K = None

''' -------------------------------------------------------------------------
    yDBL()
    input : a projective Twisted Edwards y-coordinate point y(P) := YP/WP,
            and the projective Montgomery constants A24:= A + 2C and C24:=4C
            where E : y^2 = x^3 + (A/C)*x^2 + x
    output: the projective Twisted Edwards y-coordinate point y([2]P)
    ------------------------------------------------------------------------- '''


def yDBL(P, A):

    t_0 = fp_sqr(P[0])
    t_1 = fp_sqr(P[1])
    Z = fp_mul(A[1], t_0)
    X = fp_mul(Z, t_1)
    t_1 = fp_sub(t_1, t_0)
    t_0 = fp_mul(A[0], t_1)
    Z = fp_add(Z, t_0)
    Z = fp_mul(Z, t_1)
    return [fp_sub(X, Z), fp_add(X, Z)]


''' -------------------------------------------------------------------------
    yADD()
    input : the projective Twisted Edwards y-coordinate points y(P) := YP/WP,
            y(Q) := YQ/WQ, and y(P-Q) := YPQ/QPQ
    output: the projective Twisted Edwards y-coordinate point y(P+Q)
    ------------------------------------------------------------------------- '''


def yADD(P, Q, PQ):

    a = fp_mul(P[1], Q[0])
    b = fp_mul(P[0], Q[1])
    c = fp_add(a, b)
    d = fp_sub(a, b)
    c = fp_sqr(c)
    d = fp_sqr(d)

    xD = fp_add(PQ[1], PQ[0])
    zD = fp_sub(PQ[1], PQ[0])
    X = fp_mul(zD, c)
    Z = fp_mul(xD, d)
    return [fp_sub(X, Z), fp_add(X, Z)]


''' -------------------------------------------------------------------------
    KPs()
    input : the projective Twisted Edwards y-coordinate points y(P) := YP/WP,
            the projective Montgomery constants A24:= A + 2C and C24:=4C where 
            E : y^2 = x^3 + (A/C)*x^2 + x, and a positive integer 0 <= i < n
    output: the list of projective Twisted Edwards y-coordinate points y(P),
            y([2]P), y([3]P), ..., and y([d_i]P) where l_i = 2 * d_i + 1
    ------------------------------------------------------------------------- '''


def KPs(P, A, i):

    global K
    d = (global_L[i] - 1) // 2

    K = [[0, 0] for j in range(d + 1)]
    K[0] = list([fp_sub(P[0], P[1]), fp_add(P[0], P[1])])  # 2a
    K[1] = yDBL(K[0], A)  # 4M + 2S + 4a

    for j in range(2, d, 1):
        K[j] = yADD(K[j - 1], K[0], K[j - 2])  # 4M + 2S + 6a

    return K  # 2(l - 3)M + (l - 3)S + 3(l - 3)a


''' ------------------------------------------------------------------------------
    xISOG()
    input : the projective Montgomery constants A24:= A + 2C and C24:=4C where
            E : y^2 = x^3 + (A/C)*x^2 + x, the list of projective Twisted 
            Edwards y-coordinate points y(P), y([2]P), y([3]P), ..., and y([d_i]P)
            where l_i = 2 * d_i + 1, and a positive integer 0 <= i < n
    output: the projective Montgomery constants a24:= a + 2c and c24:=4c where
            E': y^2 = x^3 + (a/c)*x^2 + x is a degree-(l_i) isogenous curve to E
    ------------------------------------------------------------------------------ '''


def xISOG(A, i):

    global K
    l = global_L[i]  # l
    bits_of_l = bitlength(l)  # Number of bits of L[i]
    d = (l - 1) // 2  # Here, l = 2d + 1

    By = K[0][0]
    Bz = K[0][1]
    for j in range(1, d, 1):

        By = fp_mul(By, K[j][0])
        Bz = fp_mul(Bz, K[j][1])

    bits_of_l -= 1
    constant_d_edwards = fp_sub(A[0], A[1])  # 1a

    tmp_a = A[0]
    tmp_d = constant_d_edwards
    # left-to-right method for computing a^l and d^l
    for j in range(1, bits_of_l + 1):

        tmp_a = fp_sqr(tmp_a)
        tmp_d = fp_sqr(tmp_d)
        if ((l >> (bits_of_l - j)) & 1) != 0:

            tmp_a = fp_mul(tmp_a, A[0])
            tmp_d = fp_mul(tmp_d, constant_d_edwards)

    for j in range(3):

        By = fp_sqr(By)
        Bz = fp_sqr(Bz)

    C0 = fp_mul(tmp_a, Bz)
    C1 = fp_mul(tmp_d, By)
    C1 = fp_sub(C0, C1)

    return [C0, C1]  # (l - 1 + 2*HW(l) - 2)M + 2(|l|_2 + 1)S + 2a


''' ------------------------------------------------------------------------------
    xEVAL()
    input : the projective Montgomery constants A24:= A + 2C and C24:=4C where
            E : y^2 = x^3 + (A/C)*x^2 + x, the list of projective Twisted 
            Edwards y-coordinate points y(P), y([2]P), y([3]P), ..., and y([d_i]P)
            where l_i = 2 * d_i + 1, and a positive integer 0 <= i < n
    output: the projective Montgomery constants a24:= a + 2c and c24:=4c where
            E': y^2 = x^3 + (a/c)*x^2 + x is a degree-(l_i) isogenous curve to E
    ------------------------------------------------------------------------------ '''


def xEVAL(P, i):

    global K
    d = (global_L[i] - 1) // 2  # Here, l = 2d + 1

    Q0 = fp_add(P[0], P[1])
    Q1 = fp_sub(P[0], P[1])
    R0, R1 = CrissCross(K[0][1], K[0][0], Q0, Q1)
    for j in range(1, d, 1):

        T0, T1 = CrissCross(K[j][1], K[j][0], Q0, Q1)
        R0 = fp_mul(T0, R0)
        R1 = fp_mul(T1, R1)

    R0 = fp_sqr(R0)
    R1 = fp_sqr(R1)
    X = fp_mul(P[0], R0)
    Z = fp_mul(P[1], R1)

    return [X, Z]  # 2(l - 1)M + 2S + (l + 1)a
