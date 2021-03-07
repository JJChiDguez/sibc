import random
from progress.bar import Bar

bitlength = lambda x: len(bin(x)[2:])  # number of bits
hamming_weight = lambda x: bin(x).count("1")
# hamming weight: number of bits equal 1

sign = lambda x: (1, -1)[x < 0]  # Sign of an integer
isequal = {True: 1, False: 0}  # Simulating constant-time integer comparison

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


