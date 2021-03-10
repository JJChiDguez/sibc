"""
This module creates test classes for all 54 permutations of prime, formula,
style, and verbose.

These obviously take a long time to run.

You can run the tests with pytest and limit which are run using the -k option.
For instance, to only run the test_genpubvalidate test on p512 with df and
tvelu, you could run this:
    pytest -k 'p512 and df and tvelu and validate'
"""

from unittest import TestCase
from sidh.csidh import CSIDH

PRIMES = ('p512', 'p1024', 'p1792')

FORMULAS = ('tvelu', 'svelu', 'hvelu')

STYLES = ('df', 'wd1', 'wd2')


class CSIDH_gae_test_base(object):
    def setUp(self):
        self.c = CSIDH(
            'montgomery', self.prime, self.formula, self.style, self.verbose, 2
        )

    def test_genpubvalidate(self):
        sk = self.c.gae.random_key()
        pk = self.c.gae.pubkey(sk)
        self.assertTrue(self.c.curve.validate(pk))

    def test_group_action_with_random_keys(self):
        sk_a, sk_b = self.c.gae.random_key(), self.c.gae.random_key()
        pk_a, pk_b = self.c.gae.pubkey(sk_a), self.c.gae.pubkey(sk_b)
        self.assertEqual(self.c.gae.dh(sk_a, pk_b), self.c.gae.dh(sk_b, pk_a))


for prime in PRIMES:
    for formula in FORMULAS:
        for style in STYLES:
            for verbose in (False, True):

                class cls(CSIDH_gae_test_base, TestCase):
                    formula = formula
                    style = style
                    prime = prime
                    verbose = verbose

                globals()[
                    'csidh_gae_'
                    + '_'.join(
                        [
                            prime,
                            formula,
                            style,
                            ('classical', 'suitable')[verbose],
                        ]
                    )
                ] = cls
del cls
