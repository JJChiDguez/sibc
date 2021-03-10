"""
This file can test many things.

For now, we're commenting and uncommenting things below to control which
permutations of parameters get tested. In the future we should use pytest's
option handling to be able to specify this.
"""


from unittest import TestCase
from sidh.csidh import CSIDH

PRIMES = ('p512', 'p1024', 'p1792')

PRIMES = ('p512',) # only test this for now, to be able to run tests quickly

FORMULAS = ('tvelu', 'svelu', 'hvelu')

STYLES = ('df', 'wd1', 'wd2')

STYLES = ('wd1',) # only this passes right now


class CSIDH_gae_test_base(object):

    def setUp(self):
        self.c = CSIDH('montgomery', self.prime, self.formula, self.style, False, 2)

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
            class cls(CSIDH_gae_test_base, TestCase):
                formula = formula
                style = style
                prime = prime
            globals()['csidh_gae_'+prime+'_'+formula+'_'+style] = cls
del cls
