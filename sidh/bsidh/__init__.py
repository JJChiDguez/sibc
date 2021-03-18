from struct import pack, unpack
from random import SystemRandom

from sidh.montgomery.curve import MontgomeryCurve
from sidh.montgomery.isogeny import MontgomeryIsogeny
from sidh.bsidh.strategy import Strategy

from sidh.constants import parameters
from sidh.common import attrdict

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

    >>> bsidh_tvelu = BSIDH('montgomery', 'b2', 'hvelu', True, False, False, False)
    >>> sk_a, sk_b = bsidh_tvelu.secret_key_a(), bsidh_tvelu.secret_key_b()
    >>> pk_a, pk_b = bsidh_tvelu.public_key_a(sk_a), bsidh_tvelu.public_key_b(sk_b)
    >>> curve_ss_a, curve_ss_b = bsidh_tvelu.dh_a(sk_a, pk_b), bsidh_tvelu.dh_b(sk_b, pk_a)
    >>> curve_ss_a == curve_ss_b
    True

    >>> from sidh.bsidh import BSIDH, default_parameters
    >>> b = BSIDH(**default_parameters)
    >>> sk_a, sk_b = b.secret_key_a(), b.secret_key_b()
    >>> pk_a, pk_b = b.public_key_a(sk_a), b.public_key_b(sk_b)
    >>> curve_ss_a, curve_ss_b = b.dh_a(sk_a, pk_b), b.dh_b(sk_b, pk_a)
    >>> curve_ss_a == curve_ss_b
    True

    Other tests which were previously here are now in the test directory.

    """

    def __init__(
        self,
        curvemodel,
        prime,
        formula,
        tuned,
        uninitialized,
        multievaluation,
        verbose
    ):

        self.params = attrdict(parameters['bsidh'][prime])
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
        return k.to_bytes(length=32, byteorder='little')

    def secret_key_b(self):
        k = self.strategy.random_scalar_B()
        return k.to_bytes(length=32, byteorder='little')

    def public_key_a(self, sk):
        # To be modified: public key will corresponds with the image point if PB, QB, and PQB
        sk = int.from_bytes(sk, byteorder='little')
        x, y = self.strategy.strategy_at_6_A(sk)
        a = x.re.x; b = x.im.x;
        c = y.re.x; d = y.im.x;
        pk = a.to_bytes(length=32, byteorder='little') + b.to_bytes(length=32, byteorder='little') +\
            c.to_bytes(length=32, byteorder='little') + d.to_bytes(length=32, byteorder='little')
        return pk

    def public_key_b(self, sk):
        # To be modified: public key will corresponds with the image point if PA, QA, and PQA
        sk = int.from_bytes(sk, byteorder='little')
        x, y = self.strategy.strategy_at_6_B(sk)
        a = x.re.x; b = x.im.x;
        c = y.re.x; d = y.im.x;
        e = a.to_bytes(length=32, byteorder='little') + b.to_bytes(length=32, byteorder='little') +\
            c.to_bytes(length=32, byteorder='little') + d.to_bytes(length=32, byteorder='little')
        return e

    def dh_a(self, sk, pk):
        sk = int.from_bytes(sk, byteorder='little')
        a, b = int.from_bytes(pk[0:32], byteorder='little'), int.from_bytes(pk[32:64], byteorder='little')
        c, d = int.from_bytes(pk[64:96], byteorder='little'), int.from_bytes(pk[96:128], byteorder='little')
        pk = [self.field([a, b]), self.field([c, d])]
        ss = self.strategy.strategy_A(sk, pk)
        curve_ss_a = self.curve.coeff(ss)
        x, y = curve_ss_a.re, curve_ss_a.im
        x = x.x.to_bytes(length=32, byteorder='little')
        y = y.x.to_bytes(length=32, byteorder='little')
        return x + y

    def dh_b(self, sk, pk):
        sk = int.from_bytes(sk, byteorder='little')
        a, b = int.from_bytes(pk[0:32], byteorder='little'), int.from_bytes(pk[32:64], byteorder='little')
        c, d = int.from_bytes(pk[64:96], byteorder='little'), int.from_bytes(pk[96:128], byteorder='little')
        pk = [self.field([a, b]), self.field([c, d])]
        ss = self.strategy.strategy_B(sk, pk)
        curve_ss_b = self.curve.coeff(ss)
        x, y = curve_ss_b.re, curve_ss_b.im
        x = x.x.to_bytes(length=32, byteorder='little')
        y = y.x.to_bytes(length=32, byteorder='little')
        return x + y

    def derive_a(self, sk, pk):
        sk = int.from_bytes(sk, byteorder='little')
        a, b = int.from_bytes(pk[0:32], byteorder='little'), int.from_bytes(pk[32:64], byteorder='little')
        c, d = int.from_bytes(pk[64:96], byteorder='little'), int.from_bytes(pk[96:128], byteorder='little')
        pk = [self.field([a, b]), self.field([c, d])]
        ss = self.strategy.strategy_A(sk, pk)
        curve_ss_a = self.curve.coeff(ss)
        x, y = curve_ss_a.re, curve_ss_a.im
        x = x.x.to_bytes(length=32, byteorder='little')
        y = y.x.to_bytes(length=32, byteorder='little')
        return x + y

    def derive_b(self, sk, pk):
        sk = int.from_bytes(sk, byteorder='little')
        a, b = int.from_bytes(pk[0:32], byteorder='little'), int.from_bytes(pk[32:64], byteorder='little')
        c, d = int.from_bytes(pk[64:96], byteorder='little'), int.from_bytes(pk[96:128], byteorder='little')
        pk = [self.field([a, b]), self.field([c, d])]
        ss = self.strategy.strategy_B(sk, pk)
        curve_ss_b = self.curve.coeff(ss)
        x, y = curve_ss_b.re, curve_ss_b.im
        x = x.x.to_bytes(length=32, byteorder='little')
        y = y.x.to_bytes(length=32, byteorder='little')
        return x + y

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
