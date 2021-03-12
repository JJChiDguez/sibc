"""
This module creates test classes for all 108 permutations of prime, formula,
style, tuned, and multievaluation.

These obviously take a long time to run. pytest-xdist is useful for speeding up
the test run with parallelized tests.

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
            'montgomery',
            self.prime,
            self.formula,
            self.style,
            self.exponent,
            self.tuned,
            self.multievaluation,
            self.verbose,
        )

    def test_group_action_with_random_keys(self):
        sk_a, sk_b = self.c.gae.random_key(), self.c.gae.random_key()
        pk_a, pk_b = self.c.gae.pubkey(sk_a), self.c.gae.pubkey(sk_b)
        self.assertTrue(self.c.curve.validate(pk_a))
        self.assertTrue(self.c.curve.validate(pk_b))
        ss_a = self.c.gae.dh(sk_a, pk_b)
        ss_b = self.c.gae.dh(sk_b, pk_a)
        self.assertNotEqual(ss_a, ss_b)
        self.assertEqual(self.c.curve.coeff(ss_a), self.c.curve.coeff(ss_b))


for prime in PRIMES:
    for formula in FORMULAS:
        for style in STYLES:
            for tuned in (False, True):
                for multievaluation in (
                    (False,) if formula == "tvelu" else (False, True)
                ):

                    class cls(CSIDH_gae_test_base, TestCase):
                        formula = formula
                        style = style
                        prime = prime
                        tuned = tuned
                        exponent = 2
                        multievaluation = multievaluation
                        verbose = False

                    globals()[
                        'csidh_gae_'
                        + '_'.join(
                            [
                                prime,
                                formula,
                                style,
                                ('classical', 'suitable')[tuned],
                                ('no-multievaluation', 'multievaluation')[
                                    multievaluation
                                ],
                            ]
                        )
                    ] = cls
del cls
