from struct import pack, unpack
from random import SystemRandom

from sibc.montgomery.curve import MontgomeryCurve
from sibc.montgomery.isogeny import MontgomeryIsogeny
from sibc.bsidh.strategy import Strategy

from sibc.constants import parameters
from sibc.common import attrdict

default_parameters = dict(
    curvemodel='montgomery',
    prime='b2',
    formula='hvelu',
    tuned=True,
    uninitialized=False,
    multievaluation=False,
    verbose=False
)

class BSIDH(object):
    """

    BSIDH

    Here is one group action test with random keys:

    >>> bsidh_hvelu = BSIDH('montgomery', 'b2', 'hvelu', True, False, False, False)
    >>> sk_a, sk_b = bsidh_hvelu.secret_key_a(), bsidh_hvelu.secret_key_b()
    >>> pk_a, pk_b = bsidh_hvelu.public_key_a(sk_a), bsidh_hvelu.public_key_b(sk_b)
    >>> ss_a, ss_b = bsidh_hvelu.dh_a(sk_a, pk_b), bsidh_hvelu.dh_b(sk_b, pk_a)
    >>> ss_a == ss_b
    True

    >>> from sibc.bsidh import BSIDH, default_parameters
    >>> b = BSIDH(**default_parameters)
    >>> sk_a, sk_b = b.secret_key_a(), b.secret_key_b()
    >>> pk_a, pk_b = b.public_key_a(sk_a), b.public_key_b(sk_b)
    >>> ss_a, ss_b = b.dh_a(sk_a, pk_b), b.dh_b(sk_b, pk_a)
    >>> ss_a == ss_b
    True

    >>> from sibc.bsidh import BSIDH, default_parameters
    >>> b = BSIDH(**default_parameters)
    >>> sk_a, pk_a = b.keygen_a()
    >>> sk_b, pk_b = b.keygen_b()
    >>> ss_a, ss_b = b.dh_a(sk_a, pk_b), b.dh_b(sk_b, pk_a)
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
        self.p_bytes = (self.params.p_bits + (self.params.p_bits % 8)) // 8
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

    def derive_b(self, sk, pk):
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

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
