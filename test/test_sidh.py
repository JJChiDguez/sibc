"""
This module creates test classes for all permutations of prime.
"""

from unittest import TestCase
from sibc.sidh import SIDH

PRIMES = ('p434', 'p503', 'p610', 'p751')

class SIDH_strategy_test_base(object):
    def setUp(self):
        self.algo = SIDH(
            curvemodel='montgomery',
            prime=self.prime,
            uninitialized=self.uninitialized,
            verbose=self.verbose,
        )

    def test_strategy_evaluation_with_random_keys(self):
        sk_a, sk_b = self.algo.strategy.random_scalar_A(), self.algo.strategy.random_scalar_B()
        pk_a, pk_b = self.algo.strategy.strategy_at_6_A(sk_a), self.algo.strategy.strategy_at_6_B(sk_b)
        # sidh does not have a validate function on the montgomery curve
        ss_a = self.algo.strategy.strategy_A(sk_a, pk_b)
        curve_ss_a = self.algo.curve.coeff(ss_a)
        ss_b = self.algo.strategy.strategy_B(sk_b, pk_a)
        curve_ss_b = self.algo.curve.coeff(ss_b)
        self.assertNotEqual(ss_a, ss_b)
        self.assertEqual(curve_ss_a, curve_ss_b)


for prime in PRIMES:
    class cls(SIDH_strategy_test_base, TestCase):
        prime = prime
        uninitialized = False
        verbose = False

    globals()[
        'sidh_strategy_evaluation_'
        + '_'.join(
            [
                prime,
            ]
        )
    ] = cls
del cls
