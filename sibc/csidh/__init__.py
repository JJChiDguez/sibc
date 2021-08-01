from struct import pack, unpack

from sibc.montgomery.curve import MontgomeryCurve
from sibc.montgomery.isogeny import MontgomeryIsogeny
from sibc.csidh.gae_df import Gae_df
from sibc.csidh.gae_wd1 import Gae_wd1
from sibc.csidh.gae_wd2 import Gae_wd2

from sibc.constants import parameters
from sibc.common import attrdict

default_parameters = dict(
    curvemodel='montgomery',
    prime='p512',
    formula='hvelu',
    style='df',
    exponent=10,
    tuned=True,
    uninitialized=False,
    multievaluation=False,
    verbose=False,
)


class CSIDH(object):
    """

    CSIDH

    Here is one group action test with random keys:

    >>> from sibc.csidh import CSIDH
    >>> csidh_tvelu_wd1 = CSIDH('montgomery', 'p512', 'tvelu', 'wd1', 10, False, False, False, False)
    >>> sk_a, sk_b = csidh_tvelu_wd1.secret_key(), csidh_tvelu_wd1.secret_key()
    >>> pk_a, pk_b = csidh_tvelu_wd1.public_key(sk_a), csidh_tvelu_wd1.public_key(sk_b)
    >>> csidh_tvelu_wd1.dh(sk_a, pk_b) == csidh_tvelu_wd1.dh(sk_b, pk_a)
    True

    >>> from sibc.csidh import CSIDH, default_parameters
    >>> csidh = CSIDH(**default_parameters)
    >>> # alice generates a key
    >>> alice_secret_key = csidh.secret_key()
    >>> alice_public_key = csidh.public_key(alice_secret_key)
    >>> # bob generates a key
    >>> bob_secret_key = csidh.secret_key()
    >>> bob_public_key = csidh.public_key(bob_secret_key)
    >>> # if either alice or bob use their secret key with the other's respective
    >>> # public key, the resulting shared secrets are the same
    >>> shared_secret_alice = csidh.dh(alice_secret_key, bob_public_key)
    >>> shared_secret_bob = csidh.dh(bob_secret_key, alice_public_key)
    >>> # Alice and bob produce an identical shared secret
    >>> shared_secret_alice == shared_secret_bob
    True

    >>> from sibc.csidh import CSIDH, default_parameters
    >>> csidh = CSIDH(**default_parameters)
    >>> # alice generates a key
    >>> alice_secret_key, alice_public_key = csidh.keygen()
    >>> # bob generates a key
    >>> bob_secret_key, bob_public_key = csidh.keygen()
    >>> # if either alice or bob use their secret key with the other's respective
    >>> # public key, the resulting shared secrets are the same
    >>> shared_secret_alice = csidh.derive(alice_secret_key, bob_public_key)
    >>> shared_secret_bob = csidh.derive(bob_secret_key, alice_public_key)
    >>> # Alice and bob produce an identical shared secret
    >>> shared_secret_alice == shared_secret_bob
    True

    Other tests which were previously here are now in the test directory.

    """

    def __init__(
        self,
        curvemodel,
        prime,
        formula,
        style,
        exponent,
        tuned,
        multievaluation,
        uninitialized,
        verbose,
    ):
        self.curvemodel = curvemodel
        self.prime = prime
        self.style = style
        self._exponent = exponent
        self.tuned = tuned
        self.uninitialized = uninitialized
        self.multievaluation = multievaluation
        self.params = attrdict(parameters['csidh'][prime])
        self.params.update(self.params[style])
        self.p_bytes = (self.params.p_bits + 8 - (self.params.p_bits % 8)) // 8

        if self.curvemodel == 'montgomery':
            self.isogeny = MontgomeryIsogeny(formula, uninitialized = self.uninitialized)
            self.curve = MontgomeryCurve(prime)
            self.field = self.curve.field
        else:
            self.curve = None
            raise NotImplemented

        self.formula = self.isogeny(self.curve, self.tuned, self.multievaluation)

        if self.style == 'df':
            self.gae = Gae_df(prime, self.tuned, self.curve, self.formula)
        elif self.style == 'wd1':
            self.gae = Gae_wd1(prime, self.tuned, self.curve, self.formula)
        elif self.style == 'wd2':
            self.gae = Gae_wd2(prime, self.tuned, self.curve, self.formula)
        else:
            self.gae = NotImplemented

    def dh(self, sk, pk):
        sk = unpack('<{}b'.format(len(sk)), sk)
        pk = int.from_bytes(pk, 'little')
        pk = self.curve.affine_to_projective(pk)
        ss = self.curve.coeff(self.gae.GAE_at_A(sk, pk)).x.to_bytes(
            length=self.p_bytes, byteorder='little'
        )
        return ss

    def derive(self, sk, pk):
        return self.dh(sk, pk)

    def secret_key(self):
        k = self.gae.random_exponents()
        return pack('<{}b'.format(len(k)), *k)

    def public_key(self, sk):
        sk = unpack('<{}b'.format(len(sk)), sk)
        xy = self.gae.GAE_at_0(sk)
        pk = self.curve.coeff(xy)
        # this implies a y of 4 on the receiver side
        return pk.x.to_bytes(length=self.p_bytes, byteorder='little')

    def keygen(self):
        sk = self.gae.random_exponents()
        xy = self.gae.GAE_at_0(sk)
        pk = self.curve.coeff(xy)
        # this implies a y of 4 on the receiver side
        return pack('<{}b'.format(len(sk)), *sk), pk.x.to_bytes(length=self.p_bytes, byteorder='little')

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
