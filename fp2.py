from fp import *

# Addition in the quadratic extension field of fp
def fp2_add(a, b):
    return [ fp_add(a[0], b[0]), fp_add(a[1], b[1]) ]

# Substraction in the quadratic extension field of fp
def fp2_sub(a, b):
    return [ fp_sub(a[0], b[0]), fp_sub(a[1], b[1]) ]

# Product in the quadratic extension field of fp
def fp2_mul(a, b):

    z0 = fp_add(a[0], a[1])
    z1 = fp_add(b[0], b[1])

    t = fp_mul(z0, z1)
    z2 = fp_mul(a[0], b[0])
    z3 = fp_mul(a[1], b[1])

    c = [0, 0]
    c[0] = fp_sub(z2, z3)
    c[1] = fp_sub(t, z2)
    c[1] = fp_sub(c[1], z3)
    return c

# Squaring in the quadratic extension field of fp
def fp2_sqr(a):

    z0 = fp_add(a[0], a[0])
    z1 = fp_add(a[0], a[1])
    z2 = fp_sub(a[0], a[1])

    b = [0,0]
    b[0] = fp_mul(z1, z2)
    b[1] = fp_mul(z0, a[1])
    return b

# Inverse in the quadratic extension field of fp
def fp2_inv(a):

    N0 = fp_sqr(a[0])
    N1 = fp_sqr(a[1])

    S1 = fp_add(N0, N1)
    S1 = fp_inv(S1)

    b = [0, 0]
    S2 = fp_sub(0, a[1])
    b[0] = fp_mul(S1, a[0])
    b[1] = fp_mul(S1, S2)
    return b

# Exponentiation in the quadratic extension field of fp
def fp2_exp(a, e):

    bits_of_e = bitlength(e)
    bits_of_e -= 1
    tmp_a = list(a)
    # left-to-right method for computing a^e
    for j in range(1, bits_of_e + 1):

        tmp_a = fp2_sqr(tmp_a)
        if( ( (e >> (bits_of_e - j)) & 1 ) != 0 ):
            tmp_a = fp2_mul(tmp_a, a)

    return tmp_a

# Sqrt (if exists) in the quadratic extension field of fp
def fp2_issquare(a):

    a1 = fp2_exp(a, p_minus_3_quarters)
    alpha = fp2_sqr(a1)
    alpha = fp2_mul(alpha, a)

    alpha_conjugated = list(alpha)
    alpha_conjugated[1] = fp_sub(0, alpha_conjugated[1])

    a0 = fp2_mul(alpha, alpha_conjugated)
    if a0[1] == 0 and a0[0] == (p - 1):
        # a doesn't have sqrt in fp2
        return False, None

    x0 = fp2_mul(a1, a)
    if alpha[1] == 0 and alpha[0] == (p - 1):
        return True, [fp_sub(0, x0[1]), x0[0]]

    else:

        alpha[0] = fp_add(alpha[0], 1)
        b = fp2_exp(alpha, p_minus_one_halves)
        b = fp2_mul(b, x0)
        return True, b

def fp2_cswap(x, y, b):
    z0, w0 = fp_cswap(x[0], y[0], b)
    z1, w1 = fp_cswap(x[1], y[1], b)
    return [z0, z1], [w0, w1]
"""
# Tests
print("p := 0x%X;" % p)
print("fp := GF(p);")
print("P<x> := PolynomialRing(fp);");
print("fp2<i> := ext<fp | x^2 + 1>;")
print("P<x> := PolynomialRing(fp2);");

a = [random.randint(0, p - 1), random.randint(0, p - 1)]
b = [random.randint(0, p - 1), random.randint(0, p - 1)]
print("a := 0x%X + 0x%X * i;" % (a[0], a[1]))
print("b := 0x%X + 0x%X * i;" % (b[0], b[1]))

e = random.randint(0, p - 1)
print("e := 0x%X;" % e)

c = fp2_add(a, b)
print("(0x%X + 0x%X * i) eq (a+b);" % (c[0], c[1]))
c = fp2_sub(a, b)
print("(0x%X + 0x%X * i) eq (a-b);" % (c[0], c[1]))
c = fp2_mul(a, b)
print("(0x%X + 0x%X * i) eq (a*b);" % (c[0], c[1]))
c = fp2_sqr(a)
print("(0x%X + 0x%X * i) eq (a^2);" % (c[0], c[1]))
c = fp2_inv(a)
print("(0x%X + 0x%X * i) eq (1/a);" % (c[0], c[1]))
c = fp2_exp(a, e)
print("(0x%X + 0x%X * i) eq (a^e);" % (c[0], c[1]))
t, c = fp2_issquare(a)
if t:
    print("(0x%X + 0x%X * i)^2 eq a;" % (c[0], c[1]))
"""
