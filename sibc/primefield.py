from random import SystemRandom
from sibc.math import bitlength, hamming_weight, is_prime
from sibc.common import check, doc

def init_runtime(field):
	field.fpadd = 0
	field.fpsqr = 0
	field.fpmul = 0

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

		def __str__(self): return hex(self.x)
		def __repr__(self): return hex(self.x)
 
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
				raise TypeError(f'The element {self} does not have square-root in the {FiniteField.__name__}')

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

	FiniteField.show_runtime = lambda label: print(
		"| %s: %7dM + %7dS + %7da"
		% (label, FiniteField.fpmul, FiniteField.fpsqr, FiniteField.fpadd),
		end="\t",
    )
	FiniteField.init_runtime = lambda: init_runtime(FiniteField)

	return FiniteField