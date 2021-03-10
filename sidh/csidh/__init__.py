#x = dict(df=gae_df, wd1=gae_wd1, wd2=gae_wd2)
from sidh.csidh.gae_df import Gae_df
from sidh.csidh.gae_wd1 import Gae_wd1
from sidh.csidh.gae_wd2 import Gae_wd2
from sidh.csidh.hvelu import Hvelu
from sidh.csidh.tvelu import Tvelu
from sidh.csidh.svelu import Svelu
from sidh.csidh.montgomery import MontgomeryCurve
from sidh.constants import parameters
from sidh.common import attrdict

class CSIDH(object):
    """

    CSIDH

    Here is one group action test with random keys:

    >>> csidh_tvelu_wd1 = CSIDH('montgomery', 'p512', 'tvelu', 'wd1', False, 2)
    >>> sk_a, sk_b = csidh_tvelu_wd1.random_key(), csidh_tvelu_wd1.random_key()
    >>> pk_a, pk_b = csidh_tvelu_wd1.pubkey(sk_a), csidh_tvelu_wd1.pubkey(sk_b)
    >>> csidh_tvelu_wd1.dh(sk_a, pk_b) == csidh_tvelu_wd1.dh(sk_b, pk_a)
    True

    Other tests which were previously here are now in the test directory.

    """

    def __init__(self, curvemodel, prime, formula, style, verbose, exponent):
        self.curvemodel = curvemodel
        self.prime = prime
        self.style = style
        self._exponent = exponent
        self.fp = None
        self.params = attrdict(parameters['csidh'][prime])
        self.params.update(self.params[style])

        # Where do we do our math? On a curve!
        if self.curvemodel == 'montgomery':
            self.curve = MontgomeryCurve(prime, style)
            self.fp = self.curve.fp
        else:
            self.curve = None
            raise NotImplemented

        # Formulas for our algorithm
        if formula == 'hvelu':
            self.formula = Hvelu(self.curve, verbose)
        elif formula == 'tvelu':
            self.formula = Tvelu(self.curve)
        elif formula == 'svelu':
            self.formula = Svelu(self.curve, verbose)

        # Side channel protection styles
        if self.style == 'df':
            self.gae = Gae_df(prime, verbose, self.curve, self.formula)
        elif self.style == 'wd1':
            self.gae = Gae_wd1(prime, verbose, self.curve, self.formula)
        elif self.style == 'wd2':
            self.gae = Gae_wd2(prime, verbose, self.curve, self.formula)
        else:
            self.gae = NotImplemented

    def dh(self, sk, pk):
        return self.gae.dh(sk, pk)

    def random_key(self):
        return self.gae.random_key()

    def pubkey(self, sk):
        return self.gae.pubkey(sk)

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
