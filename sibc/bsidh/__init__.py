from hashlib import shake_256
from random import SystemRandom

from sibc.montgomery.curve import MontgomeryCurve
from sibc.montgomery.isogeny import MontgomeryIsogeny
from sibc.bsidh.strategy import Strategy

from sibc.constants import parameters
from sibc.common import attrdict

default_parameters = dict(
    curvemodel='montgomery',
    prime='p253',
    formula='hvelu',
    tuned=True,
    uninitialized=False,
    multievaluation=False,
    verbose=False
)

class BSIDH(object):
    """

    BSIDH

    Here is one isogeny strategy evaluation test with random keys:

    >>> from sibc.bsidh import BSIDH
    >>> bsidh_hvelu = BSIDH('montgomery', 'p253', 'hvelu', True, False, False, False)
    >>> sk_a, sk_b = bsidh_hvelu.secret_key_a(), bsidh_hvelu.secret_key_b()
    >>> pk_a, pk_b = bsidh_hvelu.public_key_a(sk_a), bsidh_hvelu.public_key_b(sk_b)
    >>> ss_a, ss_b = bsidh_hvelu.dh_a(sk_a, pk_b), bsidh_hvelu.dh_b(sk_b, pk_a)
    >>> ss_a == ss_b
    True

    >>> from sibc.bsidh import BSIDH, default_parameters
    >>> bsidh = BSIDH(**default_parameters)
    >>> sk_a, sk_b = bsidh.secret_key_a(), bsidh.secret_key_b()
    >>> pk_a, pk_b = bsidh.public_key_a(sk_a), bsidh.public_key_b(sk_b)
    >>> ss_a, ss_b = bsidh.dh_a(sk_a, pk_b), bsidh.dh_b(sk_b, pk_a)
    >>> ss_a == ss_b
    True

    >>> from sibc.bsidh import BSIDH, default_parameters
    >>> bsidh = BSIDH(**default_parameters)
    >>> sk_a, pk_a = bsidh.keygen_a()
    >>> sk_b, pk_b = bsidh.keygen_b()
    >>> ss_a, ss_b = bsidh.derive_a(sk_a, pk_b), bsidh.derive_b(sk_b, pk_a)
    >>> ss_a == ss_b
    True
    
    Other tests which were previously here are now in the test directory.

    """

    def __init__(
        self,
        curvemodel,
        prime,
        formula,
        tuned,
        multievaluation,
        uninitialized,
        verbose
    ):

        self.params = attrdict(parameters['bsidh'][prime])
        self.p_bytes = (self.params.p_bits + 8 - (self.params.p_bits % 8)) // 8
        self.prime = prime
        self.tuned = tuned
        self.uninitialized = uninitialized
        self.multievaluation = multievaluation
        self.verbose = verbose

        random = SystemRandom()

        if curvemodel == 'montgomery':
            self.isogeny = MontgomeryIsogeny(formula, uninitialized = self.uninitialized)
            self.curve = MontgomeryCurve(prime)
            self.field = self.curve.field
            self.basefield = self.curve.field.basefield
        else:
            self.curve = None
            raise NotImplemented

        self.formula = self.isogeny(self.curve, self.tuned, self.multievaluation)

        if self.formula is not None and self.curve is not None:
            self.strategy = Strategy(prime, self.tuned, self.curve, self.formula)
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
        ss = self.curve.coeff(ss)
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
        ss = self.curve.coeff(ss)
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

class BSIKE(object):
    """

    SIKE

    Here is one group action test with random keys:

    >>> from sibc.bsidh import BSIKE, default_parameters
    >>> bsike = BSIKE(**default_parameters)
    >>> s, sk3, pk3 = bsike.KeyGen()
    >>> c, K = bsike.Encaps(pk3)
    >>> K_ = bsike.Decaps((s, sk3, pk3), c)
    >>> K == K_
    True

    >>> bsike255 = BSIKE('montgomery', 'p255', 'hvelu', True, False, False, False)
    >>> s, sk3, pk3 = bsike255.KeyGen()
    >>> c, K = bsike255.Encaps(pk3)
    >>> K_ = bsike255.Decaps((s, sk3, pk3), c)
    >>> K == K_
    True

    >>> bsike247 = BSIKE('montgomery', 'p247', 'hvelu', True, False, False, False)
    >>> s, sk3, pk3 = bsike247.KeyGen()
    >>> c, K = bsike247.Encaps(pk3)
    >>> K_ = bsike247.Decaps((s, sk3, pk3), c)
    >>> K == K_
    True

    >>> bsike237 = BSIKE('montgomery', 'p237', 'hvelu', True, False, False, False)
    >>> s, sk3, pk3 = bsike237.KeyGen()
    >>> c, K = bsike237.Encaps(pk3)
    >>> K_ = bsike237.Decaps((s, sk3, pk3), c)
    >>> K == K_
    True

    >>> bsike257 = BSIKE('montgomery', 'p257', 'hvelu', True, False, False, False)
    >>> s, sk3, pk3 = bsike257.KeyGen()
    >>> c, K = bsike257.Encaps(pk3)
    >>> K_ = bsike257.Decaps((s, sk3, pk3), c)
    >>> K == K_
    True

    Other tests which were previously here are now in the test directory.

    """

    def __init__(
            self,
            curvemodel,
            prime,
            formula,
            tuned,
            multievaluation,
            uninitialized,
            verbose
    ):
        self.bsidh = BSIDH(curvemodel, prime, formula, tuned, multievaluation, uninitialized, verbose)
        # Currently, BSIDH only has NIST LEVEL 1 of security
        self.n = {'p253':128, 'p255':128, 'p247':128, 'p237':128, 'p257':128}[prime]
        self.n_bytes = self.n // 8
        self.k_bytes = self.n // 8

    def bytes_xor(self, a, b):
        c = b''
        for i in range(0, self.n_bytes, 1):
            c += (a[i] ^ b[i]).to_bytes(length=1, byteorder='little')
        return c

    def Gen(self):
        return self.bsidh.keygen_b()

    def Enc(self, pk3, m, r):
        c0 = self.bsidh.public_key_a(r)
        j = self.bsidh.derive_a(r, pk3)
        h = shake_256(j).digest(self.n_bytes)
        c1 = self.bytes_xor(h, m)
        return c0, c1

    def Dec(self, sk3, c):
        (c0, c1) = c
        j = self.bsidh.derive_b(sk3, c0)
        h = shake_256(j).digest(self.n_bytes)
        m = self.bytes_xor(h, c1)
        return m

    def KeyGen(self):
        sk3, pk3 = self.Gen()
        s = self.bsidh.strategy.random.randint(0, 2**self.n).to_bytes(length=self.n_bytes, byteorder='little')
        return s, sk3, pk3

    def Encaps(self, pk3):
        m = self.bsidh.strategy.random.randint(0, 2**self.n).to_bytes(length=self.n_bytes, byteorder='little')
        r = shake_256(m + pk3).digest(self.bsidh.p_bytes)
        c0, c1 = self.Enc(pk3, m, r)
        K = shake_256(m + c0 + c1).digest(self.k_bytes)
        return (c0, c1), K

    def Decaps(self, parameters, c):
        (s, sk3, pk3) = parameters
        (c0, c1) = c
        m_ = self.Dec(sk3, c)
        r_ = shake_256(m_ + pk3).digest(self.bsidh.p_bytes)
        c0_ = self.bsidh.public_key_a(r_)
        if c0_ == c0:
            K = shake_256(m_ + c0 + c1).digest(self.k_bytes)
        else:
            K = shake_256(s + c0 + c1).digest(self.k_bytes)

        return K

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
