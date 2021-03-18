"""
This module creates test classes for all permutations of prime, formula, tuned,
and multievaluation.
"""

from unittest import TestCase
from sidh.bsidh import BSIDH

PRIMES = ('b2',)

FORMULAS = ('tvelu', 'svelu', 'hvelu')


class BSIDH_strategy_test_base(object):
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

    def test_strategy_evaluation_with_random_keys(self):
        sk_a, sk_b = self.algo.strategy.random_scalar_A(), self.algo.strategy.random_scalar_B()
        pk_a, pk_b = self.algo.strategy.strategy_at_6_A(sk_a), self.algo.strategy.strategy_at_6_B(sk_b)
        # bsidh does not have a validate function on the montgomery curve
        self.assertTrue(self.algo.curve.issupersingular(pk_a))
        self.assertTrue(self.algo.curve.issupersingular(pk_b))
        a_curve = self.algo.curve.coeff(pk_a)
        b_curve = self.algo.curve.coeff(pk_b)
        ss_a = self.algo.strategy.strategy_A(sk_a, pk_b)
        curve_ss_a = self.algo.curve.coeff(ss_a)
        ss_b = self.algo.strategy.strategy_B(sk_b, pk_a)
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

                class cls(BSIDH_strategy_test_base, TestCase):
                    formula = formula
                    prime = prime
                    tuned = tuned
                    uninitialized = False
                    multievaluation = multievaluation
                    verbose = False

                globals()[
                    'bsidh_strategy_'
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
#    len([algo for algo in dir() if algo.startswith('bsidh_strategy')]) == 81
#), "unexpected number of permutations"
