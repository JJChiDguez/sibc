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
    >>> sk_a, sk_b = bsidh_tvelu.secret_key(), bsidh_tvelu.secret_key()
    >>> pk_a, pk_b = bsidh_tvelu.public_key_a(sk_a), bsidh_tvelu.public_key_b(sk_b)
    >>> bsidh_tvelu.dh_a(sk_a, pk_b) == bsidh_tvelu.dh_b(sk_b, pk_a)
    True

    >>> from sidh.bsidh import BSIDH, default_parameters
    >>> b = BSIDH(**default_parameters)
    >>> # alice generates a key
    >>> alice_secret_key = b.secret_key()
    >>> alice_public_key = b.public_key(alice_secret_key)
    >>> # bob generates a key
    >>> bob_secret_key = b.secret_key()
    >>> bob_public_key = b.public_key(bob_secret_key)
    >>> # if either alice or bob use their secret key with the other's respective
    >>> # public key, the resulting shared secrets are the same
    >>> shared_secret_alice = b.dh(alice_secret_key, bob_public_key)
    >>> shared_secret_bob = b.dh(bob_secret_key, alice_public_key)
    >>> # Alice and bob produce an identical shared secret
    >>> shared_secret_alice == shared_secret_bob
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

    def secret_key(self):
        k = self.gae.random_key()
        return k.to_bytes(length=32, byteorder='little')

    def public_key_a(self, sk):
        sk = int.from_bytes(sk, byteorder='little')
        x, y = self.gae.pubkey_A(sk)
        return x.to_bytes(length=32, byteorder='little')+y.to_bytes(length=32, byteorder='little')

    def public_key_b(self, sk):
        sk = int.from_bytes(sk, byteorder='little')
        x, y = self.gae.pubkey_B(sk)
        return x.to_bytes(length=32, byteorder='little')+y.to_bytes(length=32, byteorder='little')

    def dh_a(self, sk, pk):
        sk = unpack('<{}b'.format(len(sk)), sk)
        pk = x, y = int.from_bytes(pk[:32], 'little'), int.from_bytes(pk[32:], 'little')
        ss = self.curve.coeff(self.gae.dh_A(sk, pk)).to_bytes(length=32, byteorder='little')
        return ss

    def dh_b(self, sk, pk):
        sk = unpack('<{}b'.format(len(sk)), sk)
        pk = x, y = int.from_bytes(pk[:32], 'little'), int.from_bytes(pk[32:], 'little')
        ss = self.curve.coeff(self.gae.dh_B(sk, pk)).to_bytes(length=32, byteorder='little')
        return ss

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
