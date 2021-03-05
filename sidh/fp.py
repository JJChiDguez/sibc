from sidh.framework import *

if (
    setting.formulaes != 'tvelu'
    and setting.formulaes != 'svelu'
    and setting.formulaes != 'hvelu'
) or (setting.algorithm != 'csidh' and setting.algorithm != 'bsidh'):

    print("  ,-~~-.___.          ")
    print(
        " / |  '     \\        SYNTAX ERROR ..., run python3 %s -h for help"
        % sys.argv[0]
    )
    print("(  )         0        ")
    print(" \_/-, ,----'         ")
    print("    ====           // ")
    print("   /  \-'~;    /~~~(O)")
    print("  /  __/~|   /       |")
    print("=(  _____| (_________|")
    exit(2)

bitlength = lambda x: len(bin(x)[2:])  # number of bits
hamming_weight = lambda x: bin(x).count("1")
# hamming weight: number of bits equal 1

sign = lambda x: (1, -1)[x < 0]  # Sign of an integer
isequal = {True: 1, False: 0}  # Simulating constant-time integer comparison

# --------------------------------------------------------------------------------------------------------------------------------
try:

    if setting.algorithm == 'csidh':

        # List of Small odd primes, L := [l_0, ..., l_{n-1}]
        f = open(sop_data + setting.prime)
        L = f.read()
        L = [int(l) for l in L.split()]
        exponent_of_two = L[0]  #   Exponent of cofactor-2
        L = list(L[1:])  #   Small Odd Primes l_i's
        n = len(L)  #   Number of l_i's to be used
        f.close()

    else:

        f = open(sop_data + setting.prime)

        # The prime to be used
        p = f.readline()
        p = int(p, 16)

        # List corresponding (p + 1)
        Lp = f.readline()
        Lp = [int(lp) for lp in Lp.split()]
        # exponent_of_twop = Lp[0]
        # Lp = Lp[1:]
        Ep = f.readline()
        Ep = [int(ep) for ep in Ep.split()]
        assert len(Ep) == len(Lp)
        np = len(Lp)

        # List corresponding (p - 1)
        Lm = f.readline()
        Lm = [int(lm) for lm in Lm.split()]
        Em = f.readline()
        Em = [int(em) for em in Em.split()]
        assert len(Em) == len(Lm)
        nm = len(Lm)

        f.close()

except IOError:

    print("  ,-~~-.___.          ")
    print(" / |  '     \\        SYNTAX ERROR ..., file not accessible")
    print("(  )         0        ")
    print(" \_/-, ,----'         ")
    print("    ====           // ")
    print("   /  \-'~;    /~~~(O)")
    print("  /  __/~|   /       |")
    print("=(  _____| (_________|")
    exit(2)

if setting.algorithm == 'csidh':

    p = (2 ** (exponent_of_two)) * reduce(
        lambda x, y: (x * y), L
    ) - 1  # p := 4 * l_0 * ... * l_n - 1
    p_minus_one_halves = (p - 1) // 2  # (p - 1) / 2

    validation_stop = sum([bitlength(l_i) for l_i in L]) / 2.0 + 2
else:

    # pp = (2**exponent_of_twop) * reduce(lambda x,y : (x*y), [ Lp[i]**Ep[i] for i in range(0, np, 1)  ])
    pp = reduce(
        lambda x, y: (x * y), [Lp[i] ** Ep[i] for i in range(0, np, 1)]
    )
    pm = reduce(
        lambda x, y: (x * y), [Lm[i] ** Em[i] for i in range(0, nm, 1)]
    )
    assert (p + 1) % pp == 0
    assert (p - 1) % pm == 0

    if p % 4 == 1:
        print("// Case p = 1 mod 4 is not implemented yet!")
        exit(-1)

    p_minus_one_halves = (p - 1) // 2
    p_minus_3_quarters = (p - 3) // 4

# --------------------------------------------------------------------------------------------------------------------------------
# Checking if p is composite


def _try_composite(a, d, n, s):
    if pow(a, d, n) == 1:
        return False
    for i in range(s):
        if pow(a, 2 ** i * d, n) == n - 1:
            return False
    return True  # n  is definitely composite


def is_prime(n):
    """
    Miller-Rabin primality test.
 
    A return value of False means n is certainly not prime. A return value of
    True means n is very likely a prime.
    """
    if n != int(n):
        return False
    n = int(n)
    # Miller-Rabin test for prime
    if n == 0 or n == 1 or n == 4 or n == 6 or n == 8 or n == 9:
        return False

    if n == 2 or n == 3 or n == 5 or n == 7:
        return True
    s = 0
    d = n - 1
    while d % 2 == 0:
        d >>= 1
        s += 1
    assert 2 ** s * d == n - 1

    def trial_composite(a):
        if pow(a, d, n) == 1:
            return False
        for i in range(s):
            if pow(a, 2 ** i * d, n) == n - 1:
                return False
        return True

    bar = Bar('// Primality test on p', max=128)
    for i in range(128):  # number of trials
        a = random.randrange(2, n)
        if trial_composite(a):
            return False
        bar.next()
    bar.finish()

    return True


if is_prime(p) == False:
    print(
        "[ERROR]\tThe integer number p := 4 * l_1, * ... * l_n - 1 is not prime where L := ",
        L,
    )
    exit(3)

"""
if( (sys.argv[0] != 'header.py') ):
    print("/*")
    print("The prime number to be used has the following form\n")
    print("           %3d" % n)
    print("        ~~~~~~~~~")
    print("         |     | ")
    print("p := 4 x |     | l_j   -   1, where each l_j is a small odd prime")
    print("           j=1  ")
    print("*/\n")
"""

# At this point, p is very likely a prime. Thus, we can continue
# --------------------------------------------------------------------------------------------------------------------------------

# Jacobi symbol used for checking if an integer has square-root in fp
def jacobi(a, n):

    assert n > a > 0 and n % 2 == 1
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            r = n % 8
            if r == 3 or r == 5:
                t = -t
        a, n = n, a
        if a % 4 == n % 4 == 3:
            t = -t
        a %= n
    if n == 1:
        return t
    else:
        return 0


# Extended GCD
def xgcd(aa, bb):
    lastremainder, remainder = abs(aa), abs(bb)
    x, lastx, y, lasty = 0, 1, 1, 0
    while remainder:
        lastremainder, (quotient, remainder) = (
            remainder,
            divmod(lastremainder, remainder),
        )
        x, lastx = lastx - quotient * x, x
        y, lasty = lasty - quotient * y, y

    return (
        lastremainder,
        lastx * (-1 if aa < 0 else 1),
        lasty * (-1 if bb < 0 else 1),
    )


# counters for field operations performed
fpadd = 0
fpsqr = 0
fpmul = 0


def set_zero_ops():

    global fpadd
    global fpsqr
    global fpmul
    fpadd -= fpadd
    fpsqr -= fpsqr
    fpmul -= fpmul


def show_ops(label, a, b, flag):

    global fpadd
    global fpsqr
    global fpmul

    print("| %s: %7dM + %7dS + %7da" % (label, fpmul, fpsqr, fpadd), end="\t")

    return None


def get_ops():

    global fpadd
    global fpsqr
    global fpmul

    return [fpmul, fpsqr, fpadd]


# Modular inverse
def fp_inv(a):

    g, x, y = xgcd(a, p)
    # if g != 1:
    # 	raise ValueError
    return x % p


# Modular addition
def fp_add(a, b):

    global fpadd
    fpadd += 1
    return (a + b) % p


# Modular substraction
def fp_sub(a, b):

    global fpadd
    fpadd += 1
    return (a - b) % p


# Modular multiplication
def fp_mul(a, b):

    global fpmul
    fpmul += 1
    return (a * b) % p


# Modular squaring
def fp_sqr(a):

    global fpsqr
    fpsqr += 1
    # print(fpsqr)
    return (a ** 2) % p


# constant-time swap
def fp_cswap(x, y, b):

    z = list([x, y])
    z = list(z[:: (1 - 2 * b)])
    return z[0], z[1]


# Modular exponentiation
def fp_exp(a, e):

    bits_of_e = bitlength(e)
    bits_of_e -= 1
    tmp_a = a
    # left-to-right method for computing a^e
    for j in range(1, bits_of_e + 1):

        tmp_a = fp_sqr(tmp_a)
        if ((e >> (bits_of_e - j)) & 1) != 0:
            tmp_a = fp_mul(tmp_a, a)

    return tmp_a


# --------------------------------------------------------------------------------------------------------------------------------
'''
    chunks()
    inputs: a string, a list, and the maximum  number of elements in each chunk
    -----
    NOTE: This function divide the input list into len(L) / k chunks.
'''
chunks = (
    lambda NAME, L, n: [NAME + ' =\t{']
    + [
        '\t' + ','.join(list(map(format, L[i * n : (i + 1) * n], ['3d'] * n)))
        for i in range((len(L) + n - 1) // n)
    ]
    + ['\t};']
)
'''
    printl()
    inputs: a string, a list, and the maximum number k of elements in each chunk
    -----
    NOTE: this function prints a given list by chunks of size k.
'''


def printl(NAME, L, k):

    to_print = chunks(NAME, L, k)
    print(to_print[0])
    for i in range(1, len(to_print) - 2):
        print(to_print[i] + ",")

    print(to_print[len(to_print) - 2])
    print(to_print[len(to_print) - 1])
