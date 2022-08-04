from random import SystemRandom
from sibc.math import bitlength, hamming_weight, is_prime
from sibc.common import check, doc
from copy import copy

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
	#@doc(NAME) # XXX this results in a 25% speedup
	class FiniteField:
		__slots__ = (r'x',)

		def __init__(self, x):
			# we have 2 mil instantiations with int during a DH
			# and zero FiniteField(FiniteField), so we optimize
			# for fast construction time; you can use
			# FiniteField(int(ff)) if you need that.
			try:
				self.x = x.x % p
			except:
				self.x = int.__mod__(x,p)

		def __invert__(self): return self.inverse()
		#@check
		def __add__(self, other):
			FiniteField.fpadd += 1
			try:
				return FiniteField(self.x + other.x)
			except AttributeError:
				return FiniteField(int.__add__(other, self.x)) #ret.x = int.__add__(other, self.x)
			#ret.x %= p
			#return ret
		#@check
		def __radd__(self, other): return self + other

		def __iadd__(self, other):
			FiniteField.fpadd += 1
			try:
				self.x = (self.x + other.x) % p
				return self
			except:
				self.x = int.__add__(self.x, other) % p
				return self

		#@check
		def __sub__(self, other):
			FiniteField.fpadd += 1
			try:
				return FiniteField(self.x - other.x)
			except AttributeError:
				return FiniteField(int.__sub__(self.x,other))

		def __isub__(self, other):
			FiniteField.fpadd += 1
			try:
				self.x = (self.x - other.x) % p
				return self
			except AttributeError:
				self.x = int.__sub__(self.x, other) % p
				return self

		#@check
		def __rsub__(self, other):
			ret = -self
			ret += other
			return ret

		#@check
		def __mul__(self, other):
			FiniteField.fpmul += 1
			try:
				return FiniteField(self.x * other.x)
			except:
				return FiniteField(int.__mul__(self.x,other))

		def __imul__(self, other):
			try:
				self.x = (self.x*other.x) % p
			except AttributeError:
				self.x = other.__mul__(self.x) % p
			FiniteField.fpmul += 1
			return self

		#@check
		def __rmul__(self, other): return self * other
		#@check
		def __truediv__(self, other): return self * ~other

		def __itruediv__(self, other):
			return self.__truediv__(self, other)

		#@check
		def __rtruediv__(self, other): return ~self * other
		#@check
		def __floordiv__(self, other): return self * ~other
		def __ifloordiv__(self, other):
			self *= other.inverse()
			return self
		#@check
		def __rfloordiv__(self, other): return ~self * other
		#@check
		def __div__(self, other): return self * ~other
		def __idiv__(self, other):
			self.x *= ~other
			self.x %= p
			return self
		#@check
		def __rdiv__(self, other): return ~self * other

		def __neg__(self): return FiniteField(-self.x)

		#@check
		def __eq__(self, other):
			try:
				assert other.__class__ is self.__class__
				return self.x == other.x
			except:
				assert other.__class__ is int
				return int.__eq__(self.x, other)

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
				FiniteField.fpsqr += (bitlength(e) - 1)
				FiniteField.fpmul += (hamming_weight(e) - 1)
				return FiniteField(pow(self.x, e, FiniteField.p))
		def __ipow__(self, e):
			if e > 0:
				FiniteField.fpsqr += (bitlength(e) - 1)
				# hamming weight of 1 and 2 is 1, subtracting 1
				# means we're adding 0 to the counter.
				# by special-casing that here we can avoid
				# ~300k function calls per CSIDH:
				if 2 < e:
					FiniteField.fpmul += (hamming_weight(e) - 1)
				self.x = pow(self.x, e, FiniteField.p)
				return self
			elif e < 0:
				self.x = pow(self.x,
					e % (FiniteField.p_minus_one),
					FiniteField.p)
				return self
			else: # e == 0
				self.x = 1
				return self


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
			return self ** ((FiniteField.p_minus_one) // 2) == 1
 
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
			return self ** (FiniteField.p - 2)					# constant-time

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

			if FiniteField.p % 4 == 3:
				return self ** int((FiniteField.p + 1) // 4)

			q = FiniteField.p_minus_one
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

		def copy(self):
			ret = object.__new__(FiniteField)
			ret.x = self.x
			return ret

		@staticmethod
		def xadd(P, Q, PQ):
			"""
			----------------------------------------------------------------------
			xadd()
			input : the projective Montgomery x-coordinate points x(P) := XP/ZP,
					x(Q) := XQ/ZQ, and x(P-Q) := XPQ/ZPQ
			output: the projective Montgomery x-coordinate point x(P+Q)
			----------------------------------------------------------------------
			"""
			FiniteField.fpadd += 6
			a = (P[0].x + P[1].x)
			b = (P[0].x - P[1].x)
			b *= (Q[0].x + Q[1].x) # c = (Q[0] + Q[1])
			a *= (Q[0].x - Q[1].x) # d = (Q[0] - Q[1])
			c = (a + b) % p
			a -= b # a = d
			FiniteField.fpmul += 4
			c = pow(c, 2, p)
			a = pow(a, 2, p)
			return [PQ[1] * c, PQ[0] * a] # X,Z

		@staticmethod
		def xdbl(P, A):
			"""
			----------------------------------------------------------------------
			xdbl()
			input : a projective Montgomery x-coordinate point x(P) := XP/ZP, and
			the  projective Montgomery constants A24:= A + 2C and C24:=4C
			where E : y^2 = x^3 + (A/C)*x^2 + x
			output: the projective Montgomery x-coordinate point x([2]P)
			----------------------------------------------------------------------
			"""
			t_0 = pow(P[0].x - P[1].x, 2, p)
			#Z = (A[1] * t_0)
			Z = A[1].x * t_0
			t_1 = pow(P[0].x + P[1].x, 2, p)
			X = FiniteField(Z * t_1)
			t_1 -= t_0
			Z += A[0].x * t_1
			Z = FiniteField(Z * t_1)
			return [X, Z]

	FiniteField.fpadd = 0  # Number of field additions performed
	FiniteField.fpsqr = 0  # Number of field squarings performed
	FiniteField.fpmul = 0  # Number of field multiplications performed
	FiniteField.p = p
	FiniteField.p_minus_one = p -1
	FiniteField.__name__ = NAME

	FiniteField.show_runtime = lambda label: print(
		"| %s: %7dM + %7dS + %7da"
		% (label, FiniteField.fpmul, FiniteField.fpsqr, FiniteField.fpadd),
		end="\t",
    )
	FiniteField.init_runtime = lambda: init_runtime(FiniteField)

	return FiniteField
