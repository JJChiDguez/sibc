from struct import pack, unpack

from sidh.bsidh.hvelu import Hvelu
from sidh.bsidh.tvelu import Tvelu
from sidh.bsidh.svelu import Svelu
from sidh.bsidh.montgomery import MontgomeryCurve
from sidh.bsidh.strategy import Gae
from sidh.constants import parameters
from sidh.common import attrdict

default_parameters = dict(curvemodel='montgomery', prime='b6',
                        formula='hvelu', exponent=2, tuned=False,
                        multievaluation=False, verbose=False)

class BSIDH(object):
    """

    BSIDH

    Here is one group action test with random keys:

    >>> bsidh_tvelu = BSIDH('montgomery', 'b6', 'tvelu', False, False, False)
    >>> sk_a, sk_b = bsidh_tvelu.secret_key(), bsidh_tvelu.secret_key()
    >>> pk_a, pk_b = bsidh_tvelu.public_key(sk_a), bsidh_tvelu.public_key(sk_b)
    >>> bsidh_tvelu.dh(sk_a, pk_b) == bsidh_tvelu.dh(sk_b, pk_a)
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
                 verbose, exponent):
        self.curvemodel = curvemodel
        self.prime = prime
        self.params = attrdict(parameters['bsidh'][prime])
        self.tuned = tuned
        self.multievaluation = multievaluation

        if self.curvemodel == 'montgomery':
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

        self.gae = Gae(prime, self.tuned, self.curve, self.formula)

    def dh(self, sk, pk):
        sk = unpack('<{}b'.format(len(sk)), sk)
        pk = int.from_bytes(pk, 'little')
        pk = self.curve.affine_to_projective(pk)
        ss = self.curve.coeff(self.gae.dh(sk, pk)).to_bytes(length=64, byteorder='little')
        return ss

    def secret_key(self):
        k = self.gae.random_key()
        return pack('<{}b'.format(len(k)), *k)

    def public_key(self, sk):
        sk = unpack('<{}b'.format(len(sk)), sk)
        xy = self.gae.pubkey(sk)
        x = self.curve.coeff(xy)
        # this implies a y of 4 on the receiver side
        return x.to_bytes(length=64, byteorder='little')


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
