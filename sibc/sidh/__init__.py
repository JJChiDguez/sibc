from hashlib import shake_256
from random import SystemRandom

from sibc.montgomery.curve import MontgomeryCurve
from sibc.montgomery.isogeny import MontgomeryIsogeny
from sibc.sidh.strategy import Strategy

from sibc.constants import parameters
from sibc.common import attrdict

default_parameters = dict(
    curvemodel='montgomery',
    prime='p434',
    uninitialized=False,
    verbose=False
)

class SIDH(object):
    """

    SIDH

    Here is one group action test with random keys:

    >>> from sibc.sidh import SIDH
    >>> sidh = SIDH('montgomery', 'p434', False, False)
    >>> sk_a, sk_b = sidh.secret_key_a(), sidh.secret_key_b()
    >>> pk_a, pk_b = sidh.public_key_a(sk_a), sidh.public_key_b(sk_b)
    >>> ss_a, ss_b = sidh.dh_a(sk_a, pk_b), sidh.dh_b(sk_b, pk_a)
    >>> ss_a == ss_b
    True

    >>> from sibc.sidh import SIDH, default_parameters
    >>> sidh = SIDH(**default_parameters)
    >>> sk_a, sk_b = sidh.secret_key_a(), sidh.secret_key_b()
    >>> pk_a, pk_b = sidh.public_key_a(sk_a), sidh.public_key_b(sk_b)
    >>> ss_a, ss_b = sidh.dh_a(sk_a, pk_b), sidh.dh_b(sk_b, pk_a)
    >>> ss_a == ss_b
    True

    >>> from sibc.sidh import SIDH, default_parameters
    >>> sidh = SIDH(**default_parameters)
    >>> sk_a, pk_a = sidh.keygen_a()
    >>> sk_b, pk_b = sidh.keygen_b()
    >>> ss_a, ss_b = sidh.dh_a(sk_a, pk_b), sidh.dh_b(sk_b, pk_a)
    >>> ss_a == ss_b
    True
    
    Other tests which were previously here are now in the test directory.

    """

    def __init__(
        self,
        curvemodel,
        prime,
        uninitialized,
        verbose
    ):

        self.params = attrdict(parameters['sidh'][prime])
        self.p_bytes = (self.params.p_bits + 8 - (self.params.p_bits % 8)) // 8
        self.prime = prime
        self.uninitialized = uninitialized
        self.verbose = verbose

        random = SystemRandom()

        if curvemodel == 'montgomery':
            self.isogeny = MontgomeryIsogeny('tvelu', uninitialized = self.uninitialized)
            self.curve = MontgomeryCurve(prime)
            self.field = self.curve.field
            self.basefield = self.curve.field.basefield
        else:
            self.curve = None
            raise NotImplemented

        self.formula = self.isogeny(self.curve, False , False)

        if self.formula is not None and self.curve is not None:
            self.strategy = Strategy(prime, self.curve, self.formula)
        else:
            self.strategy = None
            raise NotImplemented

    def secret_key_a(self):
        k = self.strategy.random_scalar_A()
        return k.to_bytes(length=self.p_bytes, byteorder='little')

    def secret_key_b(self):
        k = self.strategy.random_scalar_B()
        return k.to_bytes(length=self.p_bytes, byteorder='little')

    def public_key_a(self, sk):
        sk = int.from_bytes(sk, byteorder='little')
        x, y, z = self.strategy.strategy_at_6_A(sk)
        pk = x.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        x.im.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        y.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        y.im.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        z.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        z.im.x.to_bytes(length=self.p_bytes, byteorder='little')
        return pk

    def public_key_b(self, sk):
        sk = int.from_bytes(sk, byteorder='little')
        x, y, z = self.strategy.strategy_at_6_B(sk)
        pk = x.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        x.im.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        y.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        y.im.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        z.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        z.im.x.to_bytes(length=self.p_bytes, byteorder='little')
        return pk

    def dh_a(self, sk, pk):
        sk = int.from_bytes(sk, byteorder='little')
        # --- P
        P_re = int.from_bytes(pk[0:self.p_bytes], byteorder='little')
        P_im = int.from_bytes(pk[self.p_bytes:(2*self.p_bytes)], byteorder='little')
        # --- Q
        Q_re = int.from_bytes(pk[(2*self.p_bytes):(3*self.p_bytes)], byteorder='little')
        Q_im = int.from_bytes(pk[(3*self.p_bytes):(4*self.p_bytes)], byteorder='little')
        # --- Q
        PQ_re = int.from_bytes(pk[(4*self.p_bytes):(5*self.p_bytes)], byteorder='little')
        PQ_im = int.from_bytes(pk[(5*self.p_bytes):(6*self.p_bytes)], byteorder='little')

        pk = (self.field([P_re, P_im]), self.field([Q_re, Q_im]), self.field([PQ_re, PQ_im]))
        ss = self.strategy.strategy_A(sk, pk)
        ss = self.curve.jinvariant(ss)
        x = ss.re.x.to_bytes(length=self.p_bytes, byteorder='little')
        y = ss.im.x.to_bytes(length=self.p_bytes, byteorder='little')
        return x + y

    def dh_b(self, sk, pk):
        sk = int.from_bytes(sk, byteorder='little')
        # --- P
        P_re = int.from_bytes(pk[0:self.p_bytes], byteorder='little')
        P_im = int.from_bytes(pk[self.p_bytes:(2*self.p_bytes)], byteorder='little')
        # --- Q
        Q_re = int.from_bytes(pk[(2*self.p_bytes):(3*self.p_bytes)], byteorder='little')
        Q_im = int.from_bytes(pk[(3*self.p_bytes):(4*self.p_bytes)], byteorder='little')
        # --- Q
        PQ_re = int.from_bytes(pk[(4*self.p_bytes):(5*self.p_bytes)], byteorder='little')
        PQ_im = int.from_bytes(pk[(5*self.p_bytes):(6*self.p_bytes)], byteorder='little')

        pk = (self.field([P_re, P_im]), self.field([Q_re, Q_im]), self.field([PQ_re, PQ_im]))
        ss = self.strategy.strategy_B(sk, pk)
        ss = self.curve.jinvariant(ss)
        x = ss.re.x.to_bytes(length=self.p_bytes, byteorder='little')
        y = ss.im.x.to_bytes(length=self.p_bytes, byteorder='little')
        return x + y

    def derive_a(self, sk, pk):
        return self.dh_a(sk, pk)

    def derive_b(self, sk, pk):
        return self.dh_b(sk, pk)

    def keygen_a(self):
        sk = self.strategy.random_scalar_A()
        x, y, z = self.strategy.strategy_at_6_A(sk)
        sk = sk.to_bytes(length=self.p_bytes, byteorder='little')
        pk = x.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        x.im.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        y.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        y.im.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        z.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        z.im.x.to_bytes(length=self.p_bytes, byteorder='little')
        return sk, pk

    def keygen_b(self):
        sk = self.strategy.random_scalar_B()
        x, y, z = self.strategy.strategy_at_6_B(sk)
        sk = sk.to_bytes(length=self.p_bytes, byteorder='little')
        pk = x.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        x.im.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        y.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        y.im.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        z.re.x.to_bytes(length=self.p_bytes, byteorder='little') +\
        z.im.x.to_bytes(length=self.p_bytes, byteorder='little')
        return sk, pk

class SIKE(object):
    """

    SIKE

    Here is one group action test with random keys:

    >>> from sibc.sidh import SIKE, default_parameters
    >>> sike = SIKE(**default_parameters)
    >>> s, sk3, pk3 = sike.KeyGen()
    >>> c, K = sike.Encaps(pk3)
    >>> K_ = sike.Decaps((s, sk3, pk3), c)
    >>> K == K_
    True

    >>> sike503 = SIKE('montgomery', 'p503', False, False)
    >>> s, sk3, pk3 = sike503.KeyGen()
    >>> c, K = sike503.Encaps(pk3)
    >>> K_ = sike503.Decaps((s, sk3, pk3), c)
    >>> K == K_
    True

    >>> sike610 = SIKE('montgomery', 'p610', False, False)
    >>> s, sk3, pk3 = sike610.KeyGen()
    >>> c, K = sike610.Encaps(pk3)
    >>> K_ = sike610.Decaps((s, sk3, pk3), c)
    >>> K == K_
    True

    >>> sike751 = SIKE('montgomery', 'p751', False, False)
    >>> s, sk3, pk3 = sike751.KeyGen()
    >>> c, K = sike751.Encaps(pk3)
    >>> K_ = sike751.Decaps((s, sk3, pk3), c)
    >>> K == K_
    True

    Other tests which were previously here are now in the test directory.

    """

    def __init__(
            self,
            curvemodel,
            prime,
            uninitialized,
            verbose
    ):
        self.sidh = SIDH(curvemodel, prime, uninitialized, verbose)
        self.n = {'p434':128, 'p503':128, 'p610':192, 'p751':256}[prime]
        self.n_bytes = self.n // 8
        self.k_bytes = self.n // 8

    def bytes_xor(self, a, b):
        c = b''
        for i in range(0, self.n_bytes, 1):
            c += (a[i] ^ b[i]).to_bytes(length=1, byteorder='little')
        return c

    def Gen(self):
        return self.sidh.keygen_b()

    def Enc(self, pk3, m, r):
        c0 = self.sidh.public_key_a(r)
        j = self.sidh.derive_a(r, pk3)
        h = shake_256(j).digest(self.n_bytes)
        c1 = self.bytes_xor(h, m)
        return c0, c1

    def Dec(self, sk3, c):
        (c0, c1) = c
        j = self.sidh.derive_b(sk3, c0)
        h = shake_256(j).digest(self.n_bytes)
        m = self.bytes_xor(h, c1)
        return m

    def KeyGen(self):
        sk3, pk3 = self.Gen()
        s = self.sidh.strategy.random.randint(0, 2**self.n).to_bytes(length=self.n_bytes, byteorder='little')
        return s, sk3, pk3

    def Encaps(self, pk3):
        m = self.sidh.strategy.random.randint(0, 2**self.n).to_bytes(length=self.n_bytes, byteorder='little')
        r = shake_256(m + pk3).digest(self.sidh.strategy.two // 8)
        c0, c1 = self.Enc(pk3, m, r)
        K = shake_256(m + c0 + c1).digest(self.k_bytes)
        return (c0, c1), K

    def Decaps(self, parameters, c):
        (s, sk3, pk3) = parameters
        (c0, c1) = c
        m_ = self.Dec(sk3, c)
        r_ = shake_256(m_ + pk3).digest(self.sidh.strategy.two // 8)
        c0_ = self.sidh.public_key_a(r_)
        if c0_ == c0:
            K = shake_256(m_ + c0 + c1).digest(self.k_bytes)
        else:
            K = shake_256(s + c0 + c1).digest(self.k_bytes)

        return K

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
