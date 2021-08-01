"""
This module creates test classes for all 81 permutations of prime, formula,
style, tuned, and multievaluation.

These obviously take a long time to run. pytest-xdist is useful for speeding up
the test run with parallelized tests.

You can run the tests with pytest and limit which are run using the -k option.
For instance, to only run the test_genpubvalidate test on p512 with df and
tvelu, you could run this:
    pytest -k 'p512 and df and tvelu and validate'
"""

from unittest import TestCase
from sibc.csidh import CSIDH

configuration=(
    ('p512', 'df', '10'),
    ('p512', 'wd1', '10'),
    ('p512', 'wd2', '5'),
    ('p1024', 'df', '3'),
    ('p1024', 'wd1', '3'),
    ('p1024', 'wd2', '2'),
    ('p1792', 'df', '2'),
    ('p1792', 'wd1', '2'),
    ('p1792', 'wd2', '1'),
    ('p2048', 'df', '1'),
    ('p2048', 'wd1', '1'),
    ('p2048', 'wd2', '1'),
    ('p4096', 'df', '1'),
    ('p4096', 'wd1', '1'),
    ('p4096', 'wd2', '1'),
    ('p5120', 'df', '1'),
    ('p5120', 'wd1', '1'),
    ('p5120', 'wd2', '1'),
    ('p6144', 'df', '1'),
    ('p6144', 'wd1', '1'),
    ('p6144', 'wd2', '1'),
    ('p8192', 'df', '1'),
    ('p8192', 'wd1', '1'),
    ('p8192', 'wd2', '1'),
    ('p9216', 'df', '1'),
    ('p9216', 'wd1', '1'),
    ('p9216', 'wd2', '1'),
)

FORMULAS = ('tvelu', 'svelu', 'hvelu')

class CSIDH_gae_test_base(object):
    def setUp(self):
        self.c = CSIDH(
            'montgomery',
            prime=self.prime,
            formula=self.formula,
            style=self.style,
            exponent=self.exponent,
            tuned=self.tuned,
            multievaluation=self.multievaluation,
            uninitialized=self.uninitialized,
            verbose=self.verbose,
        )

    def test_group_action_with_random_exponents(self):
        sk_a, sk_b = self.c.gae.random_exponents(), self.c.gae.random_exponents()
        pk_a, pk_b = self.c.gae.GAE_at_0(sk_a), self.c.gae.GAE_at_0(sk_b)
        self.assertTrue(self.c.curve.issupersingular(pk_a))
        self.assertTrue(self.c.curve.issupersingular(pk_b))
        ss_a = self.c.gae.GAE_at_A(sk_a, pk_b)
        ss_b = self.c.gae.GAE_at_A(sk_b, pk_a)
        self.assertNotEqual(ss_a, ss_b)
        self.assertEqual(self.c.curve.coeff(ss_a), self.c.curve.coeff(ss_b))


for prime, style, exponent in configuration:
    for formula in FORMULAS:
        for tuned in (False, True):
            for multievaluation in (False, True):
                if formula == "tvelu" and (tuned or multievaluation):
                    # tvelu doesn't have tuned or multievaluation modes
                    continue

                class cls(CSIDH_gae_test_base, TestCase):
                    prime = prime
                    formula = formula
                    style = style
                    exponent = exponent
                    tuned = tuned
                    multievaluation = multievaluation
                    uninitialized = False
                    verbose = False

                globals()[
                    'csidh_gae_'
                    + '_'.join(
                        [
                            prime,
                            formula,
                            style,
                            ('', 'tuned')[tuned],
                            ('unscaled', 'scaled')[
                                multievaluation
                            ],
                        ]
                    )
                ] = cls
del cls
