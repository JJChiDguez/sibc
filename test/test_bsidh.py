"""
This module creates test classes for all permutations of prime, formula, tuned,
and multievaluation.
"""

from unittest import TestCase
from sibc.bsidh import BSIDH

PRIMES = ('p253', 'p255', 'p247', 'p237', 'p257')

FORMULAS = ('tvelu', 'svelu', 'hvelu')


class BSIDH_strategy_test_base(object):
    def setUp(self):
        self.algo = BSIDH(
            curvemodel='montgomery',
            prime=self.prime,
            formula=self.formula,
            tuned=self.tuned,
            multievaluation=self.multievaluation,
            uninitialized=self.uninitialized,
            verbose=self.verbose,
        )

    def test_strategy_evaluation_with_random_keys(self):
        sk_a, sk_b = self.algo.strategy.random_scalar_A(), self.algo.strategy.random_scalar_B()
        pk_a, pk_b = self.algo.strategy.strategy_at_6_A(sk_a), self.algo.strategy.strategy_at_6_B(sk_b)
        # bsidh does not have a validate function on the montgomery curve
        (PA_b, QA_b, PQA_b) = pk_b
        self.assertTrue(self.algo.curve.issupersingular(self.algo.curve.get_A(
            [PA_b, self.algo.field(1)],
            [QA_b, self.algo.field(1)],
            [PQA_b, self.algo.field(1)]
        )))
        (PB_a, QB_a, PQB_a) = pk_a
        self.assertTrue(self.algo.curve.issupersingular(self.algo.curve.get_A(
            [PB_a, self.algo.field(1)],
            [QB_a, self.algo.field(1)],
            [PQB_a, self.algo.field(1)]
        )))
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
                    multievaluation = multievaluation
                    uninitialized = False
                    verbose = False

                globals()[
                    'bsidh_strategy_evaluation_'
                    + '_'.join(
                        [
                            prime,
                            formula,
                            ('', 'tuned')[tuned],
                            ('unscaled', 'scaled')[
                                multievaluation
                            ],
                        ]
                    )
                ] = cls
del cls
