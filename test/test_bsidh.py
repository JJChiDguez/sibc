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
            exponent=self.exponent,
            tuned=self.tuned,
            multievaluation=self.multievaluation,
            verbose=self.verbose,
        )

    def test_group_action_with_random_keys(self):
        sk_a, sk_b = self.algo.gae.random_key(), self.algo.gae.random_key()
        pk_a, pk_b = self.algo.gae.pubkey(sk_a), self.algo.gae.pubkey(sk_b)
        # bsidh does not have a validate function on the montgomery curve
        # self.assertTrue(self.algo.curve.validate(pk_a))
        # self.assertTrue(self.algo.curve.validate(pk_b))
        ss_a = self.algo.gae.dh(sk_a, pk_b)
        ss_b = self.algo.gae.dh(sk_b, pk_a)
        self.assertNotEqual(ss_a, ss_b)
        self.assertEqual(self.algo.curve.coeff(ss_a), self.algo.curve.coeff(ss_b))


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
                    exponent = 2
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
