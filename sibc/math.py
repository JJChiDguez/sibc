from random import SystemRandom
#from progress.bar import Bar

# number of bits, use builtin int.bit_length if present:
bitlength = getattr(int, 'bit_length', lambda x: len(bin(x)[2:]))

# python3.10 has builtin popcount aka hamming weight aka int.bit_count:
hamming_weight = getattr(int, 'bit_count', lambda x: bin(x).count(r'1'))
# hamming weight: number of bits equal 1

sign = lambda x: (1, -1)[x < 0]  # Sign of an integer
isequal = {True: 1, False: 0}  # Simulating constant-time integer comparison
random = SystemRandom()

# constant-time swap (this function should be modified into bit opterations for ensrue constant-time)
def cswap(x, y, b):
    z = list([x, y])
    z = list(z[:: (1 - 2 * b)])
    return z[0], z[1]

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

    #bar = Bar('// Primality test on p', max=128)
    iters = 1
    for i in range(iters):  # number of trials
        a = random.randrange(2, n)
        if trial_composite(a):
            return False
        #bar.next()
    #bar.finish()

    return True


def Brent_Cycle_finding_method(N : int):
    random = SystemRandom()
    if N % 2 == 0:
        return 2
    y, c, m = random.randint(1, N - 1), random.randint(1, N - 1), random.randint(1, N - 1)
    g, r, q = 1, 1, 1
    while g == 1:
        x = y
        for i in range(r):
            y = ((y * y) % N + c) % N
        k = 0
        while k < r and g == 1:
            ys = y
            for i in range(min(m, r - k)):
                y = ((y * y) % N + c) % N
                q = q * (abs(x - y)) % N
            g = gcd(q, N)
            k += m
        r *= 2
    if g == N:
        while True:
            ys = ((ys * ys) % N + c) % N
            g = gcd(abs(x - ys), N)
            if g > 1:
                break
    return g

def Factorint(n : int):

    # Random factorization procedure by using Brent's Cycle finding method
    if n == 1:
        return {}
    if is_prime(n):
        return {str(n): 1}
    else:
        m = Brent_Cycle_finding_method(n)   # Factoring n
        factors_0 = Factorint(m)        # Factoring m
        factors_1 = Factorint(n // m)   # Factoring n // m

        # Unifying results
        tmp = {k: factors_0[k] for k in factors_0.keys() - factors_1}
        tmp = dict({k: factors_1[k] for k in factors_1.keys() - factors_0}, **tmp)
        tmp = dict({k: factors_0[k] + factors_1[k] for k in factors_0.keys() & factors_1}, **tmp)
        return tmp

def Factorization(n : int):
    """
    Integer factorization using the Brent Cycle finding method

    ...
    Parameters
    ----------
        - positive integer n
    Returns
    -------
        - Factorization of n
    Notes
    -----
        - The output is given as a dictionary where the keys are the prime factors 
        and the values their multiplicities
    """
    return dict(sorted(Factorint(n).items(), key = lambda kv: int(kv[0])))
