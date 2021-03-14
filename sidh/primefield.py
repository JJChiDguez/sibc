from random import SystemRandom
from sidh.math import bitlength, hamming_weight, is_prime
from functools import wraps

# ---> Next three functions should be moved to common.py
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
            k = k + m
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
		m = Brent_Cycle_finding_method(n)	# Factoring n
		factors_0 = Factorint(m)		# Factoring m
		factors_1 = Factorint(n // m)	# Factoring n // m

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

# <---

def check(func):
	@wraps(func)
	def method(self, other):
		if type(self) is not type(other):
			other = self.__class__(other)
		else:
			if self.field.p != other.field.p:
				raise ValueError
		return func(self, other)

	return method

def doc(s):
	class __doc(object):
		def __init__(self,f):
			self.func = f
			self.desc = s
		def __call__(self,*args,**kwargs):
			return self.func(*args,**kwargs)
		def __repr__(self):
			return self.desc
	return __doc

def PrimeField(p : int):
	"""
	Prime Field class constructor

	...
	Parameters
	----------
		- prime number p
	Returns
	-------
		- Prime Field class of characteristic p
	"""

	if not is_prime(p):
		# To remove the use of progress.bar in is_prime() ... and maybe to renamed as IsPrime()
		raise TypeError(f'The integer {p} is not a prime number, and thus it does not allow to construct a prime field')

	NAME = 'Prime Field GF(p) of characteristic p = 0x%X' % (p)
	@doc(NAME)
	class FiniteField():

		def __init__(self, x):
			self.x = x % p if isinstance(x, int) else x.x
			self.field = FiniteField

		def __invert__(self): return self.inverse()
		@check
		def __add__(self, other): self.field.fpadd += 1; return FiniteField(self.x + other.x)
		@check
		def __radd__(self, other): return self + other
		@check
		def __sub__(self, other): self.field.fpadd += 1; return FiniteField(self.x - other.x)
		@check
		def __rsub__(self, other): return -self + other
		@check
		def __mul__(self, other): self.field.fpmul += 1; return FiniteField(self.x * other.x)
		@check
		def __rmul__(self, other): return self * other
		@check
		def __truediv__(self, other): return self * ~other
		@check
		def __rtruediv__(self, other): return ~self * other
		@check
		def __floordiv__(self, other): return self * ~other
		@check
		def __rfloordiv__(self, other): return ~self * other
		@check
		def __div__(self, other): return self * ~other
		@check
		def __rdiv__(self, other): return ~self * other

		def __neg__(self): return FiniteField(-self.x)

		@check
		def __eq__(self, other): return isinstance(other, self.__class__) and self.x == other.x

		def __abs__(self):
			"""
			Signed representation of a prime field element

			...
			Parameters
			----------
				- self an element of a Prime Field
			Returns
			-------
				- an integer x belonging to |[ (-p+1)/2 .. (p+1)/2 ]| such that self = x modulo the characteristic of the Prime Field
			-----
			Usage: self.abs()

			"""
			neg = (-self)
			if self.x < neg.x:
				return self.x
			else:
				return -neg.x

		def __str__(self): return str(self.x)
		def __repr__(self): return str(self.x)
 
		def __divmod__(self, divisor):
			q,r = divmod(self.x, divisor.x)
			return (FiniteField(q), FiniteField(r))

		def __pow__(self, e):
			"""
			Exponentiation

			...
			Parameters
			----------
				- self, which is an element of a Prime Field
				- an integer e
			Returns
			-------
				- self raised to e
			-----
			Usage:
				- self.pow(e)
				- self ** e
			Notes
			-----
				- This is a constant-time implementation by using the left-to-right method
				- It allows negative exponents, but any exponent is expected to belong to |[ 0 .. p - 1 ]|
			"""
			if e == 0:
				return FiniteField(1)

			elif e < 0:
				return self.inverse() ** (-e)

			else:
				self.field.fpsqr += (bitlength(e) - 1)
				self.field.fpmul += (hamming_weight(e) - 1)
				return FiniteField(pow(self.x, e, self.field.p))

		def issquare(self):
			"""
			Checking if a given element is a quadratic residue

			...
			Parameters
			----------
				- self, which is an element of a Prime Field
			Returns
			-------
				- True if self is a quadratic residue; otherwise, False
			-----
			Usage:
				- self.issquare()
			Notes
			-----
				- This is a constant-time implementation by rasing to (p - 1) / 2
				- In other words, this function determines if the input has square-root in the Prime Field
			"""
			return self ** ((self.field.p - 1) // 2) == 1
 
		def inverse(self):
			"""
			Multiplicative inverse computation

			...
			Parameters
			----------
				- self, which is an element of a Prime Field
			Returns
			-------
				- the multiplivative inverse of self
			-----
			Usage:
				- self.inverse()
				- self ** -1
				- 1 / self, which performs an extra field multiplication
			Notes
			-----
				- This is a constant-time implementation by raising to (p - 2)
			"""
			return self ** (self.field.p - 2)					# constant-time

		def sqrt(self):
			"""
			Square-root computation by using the Tonelli-Shanks algorithm

			...
			Parameters
			----------
				- self, which is an element of a Prime Field
			Returns
			-------
				- a square-root of self
			-----
			Usage:
				- self.sqrt()
			Notes
			-----
				- This is a non-constant-time implementation but it is only used on public data
			"""
			
			if self == 0:
				return self

			if not self.issquare():
				raise TypeError(f'The element {self} does not have square-root in the prime field {FiniteField.__name__}')

			if self.field.p % 4 == 3:
				return self ** int((self.field.p + 1) // 4)

			q = self.field.p - 1
			s = 0
			while q % 2 == 0:
				q //= 2
				s += 1

			z = self.__class__(2)
			while z.issquare():
				z += 1

			m = s
			c = z ** int(q)
			t = self ** int(q)
			r_exp = (q + 1) // 2
			r = self ** int(r_exp)

			while t != 1:
				i = 1
				while not (t ** (2 ** i)) == 1:
					i += 1
				two_exp = m - (i + 1)
				b = c ** (self.__class__(2) ** two_exp).x
				m = i
				c = b ** 2
				t *= c
				r *= b
			return r
 
	FiniteField.fpadd = 0  # Number of field additions performed
	FiniteField.fpsqr = 0  # Number of field squarings performed
	FiniteField.fpmul = 0  # Number of field multiplications performed
	FiniteField.p = p
	FiniteField.__name__ = NAME

	return FiniteField

# ---> Next three functions should be moved to common.py

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
# <---

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
