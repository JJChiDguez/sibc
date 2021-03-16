from struct import pack, unpack

from sympy import symbols, floor, sqrt, sign
from pkg_resources import resource_filename
from random import SystemRandom

from sidh.bsidh.hvelu import Hvelu
from sidh.bsidh.tvelu import Tvelu
from sidh.bsidh.svelu import Svelu
from sidh.bsidh.montgomery import MontgomeryCurve
from sidh.bsidh.strategy import Gae
from sidh.constants import parameters
from sidh.common import attrdict

default_parameters = dict(curvemodel='montgomery', prime='b2',
                        formula='hvelu', tuned=False,
                        multievaluation=False, verbose=False)

class BSIDH(object):
    """

    BSIDH

    Here is one group action test with random keys:

    >>> bsidh_tvelu = BSIDH('montgomery', 'b2', 'hvelu', False, False, False)
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

    def __init__(self, curvemodel, prime, formula, tuned, multievaluation,
                 verbose):

        self.params = attrdict(parameters['bsidh'][prime])
        self.prime = prime
        self.tuned = tuned
        self.multievaluation = multievaluation
        self.verbose = verbose

        random = SystemRandom()

        if curvemodel == 'montgomery':
            self.curve = MontgomeryCurve(prime)
            self.A = self.curve.A
            self.fp = self.curve.fp
        else:
            self.curve = None
            raise NotImplemented

        if formula == 'hvelu':
            self.formula = Hvelu(self.curve, self.tuned, self.multievaluation)
        elif formula == 'tvelu':
            self.formula = Tvelu(self.curve)
        elif formula == 'svelu':
            self.formula = Svelu(self.curve, self.tuned, self.multievaluation)
        else:
            self.formula = None
            raise NotImplemented

        if self.formula is not None and self.curve is not None:
            self.gae = Gae(prime, self.tuned, self.curve, self.formula)
        else:
            self.gae = None
            raise NotImplemented

    def secret_key_a(self):
        k = self.gae.random_key_A()
        return k.to_bytes(length=32, byteorder='little')

    def secret_key_b(self):
        k = self.gae.random_key_B()
        return k.to_bytes(length=32, byteorder='little')

    def public_key_a(self, sk):
        sk = int.from_bytes(sk, byteorder='little')
        x, y = self.gae.pubkey_A(sk)
        a = x[0]; b = x[1];
        c = y[0]; d = y[1];
        pk = a.to_bytes(length=32, byteorder='little') + b.to_bytes(length=32, byteorder='little') +\
            c.to_bytes(length=32, byteorder='little') + d.to_bytes(length=32, byteorder='little')
        return pk

    def public_key_b(self, sk):
        sk = int.from_bytes(sk, byteorder='little')
        x, y = self.gae.pubkey_B(sk)
        a = x[0]; b = x[1];
        c = y[0]; d = y[1];
        e = a.to_bytes(length=32, byteorder='little') + b.to_bytes(length=32, byteorder='little') +\
            c.to_bytes(length=32, byteorder='little') + d.to_bytes(length=32, byteorder='little')
        return e

    def dh_a(self, sk, pk):
        sk = int.from_bytes(sk, byteorder='little')
        a, b = int.from_bytes(pk[0:32], byteorder='little'), int.from_bytes(pk[32:64], byteorder='little')
        c, d = int.from_bytes(pk[64:96], byteorder='little'), int.from_bytes(pk[96:128], byteorder='little')
        pk = [(a, b), (c, d)]
        ss = self.gae.dh_A(sk, pk)
        curve_ss_a = self.curve.coeff(ss)
        x, y = curve_ss_a
        x = x.to_bytes(length=32, byteorder='little')
        y = y.to_bytes(length=32, byteorder='little')
        return x + y

    def dh_b(self, sk, pk):
        sk = int.from_bytes(sk, byteorder='little')
        a, b = int.from_bytes(pk[0:32], byteorder='little'), int.from_bytes(pk[32:64], byteorder='little')
        c, d = int.from_bytes(pk[64:96], byteorder='little'), int.from_bytes(pk[96:128], byteorder='little')
        pk = [(a, b), (c, d)]
        ss = self.gae.dh_B(sk, pk)
        curve_ss_b = self.curve.coeff(ss)
        x, y = curve_ss_b
        x = x.to_bytes(length=32, byteorder='little')
        y = y.to_bytes(length=32, byteorder='little')
        return x + y

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
