"""
This module creates test classes for all permutations of prime, formula, tuned,
and multievaluation.
"""

from unittest import TestCase
from sidh.bsidh import BSIDH

PRIMES = ('b2',)

FORMULAS = ('tvelu', 'svelu', 'hvelu')


class BSIDH_gae_test_base(object):
    def setUp(self):
        self.algo = BSIDH(
            curvemodel='montgomery',
            prime=self.prime,
            formula=self.formula,
            uninitialized=self.uninitialized,
            tuned=self.tuned,
            multievaluation=self.multievaluation,
            verbose=self.verbose,
        )

    def test_group_action_with_random_keys(self):
        sk_a, sk_b = self.algo.gae.random_scalary_A(), self.algo.gae.random_key_B()
        pk_a, pk_b = self.algo.gae.strategy_at_6_A(sk_a), self.algo.gae.strategy_at_6_B(sk_b)
        # bsidh does not have a validate function on the montgomery curve
        self.assertTrue(self.algo.curve.issupersingular(pk_a))
        self.assertTrue(self.algo.curve.issupersingular(pk_b))
        a_curve = self.algo.curve.coeff(pk_a)
        b_curve = self.algo.curve.coeff(pk_b)
        ss_a = self.algo.gae.strategy_A(sk_a, pk_b)
        curve_ss_a = self.algo.curve.coeff(ss_a)
        ss_b = self.algo.gae.strategy_B(sk_b, pk_a)
        curve_ss_b = self.algo.curve.coeff(ss_b)
        self.assertNotEqual(ss_a, ss_b)
        self.assertEqual(curve_ss_a, curve_ss_b)


for prime in PRIMES:
    for formula in FORMULAS:
        for tuned in (False, True):
            for multievaluation in (False, True):
                if formula == "tvelu" and (tuned or multievaluation):
                    # tvelu doesn't have tuned or multievaluation modes
                    continue

                class cls(BSIDH_gae_test_base, TestCase):
                    formula = formula
                    prime = prime
                    tuned = tuned
                    uninitialized = False
                    multievaluation = multievaluation
                    verbose = False

                globals()[
                    'bsidh_gae_'
                    + '_'.join(
                        [
                            prime,
                            formula,
                            ('classical', 'suitable')[tuned],
                            ('no-multievaluation', 'multievaluation')[
                                multievaluation
                            ],
                        ]
                    )
                ] = cls
del cls
#assert (
#    len([algo for algo in dir() if algo.startswith('bsidh_gae')]) == 81
#), "unexpected number of permutations"
