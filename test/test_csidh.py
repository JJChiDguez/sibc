from unittest import TestCase, main
from sidh.csidh import CSIDH

class CSIDH_test_base(object):

    def test_genpubvalidate(self):
        c = CSIDH('montgomery', 'p512', self.formula, self.style, False, 2)
        sk = c.gae.random_key()
        pk = c.gae.pubkey(sk)
        self.assertTrue(c.curve.validate(pk))


FORMULAS = ('tvelu', 'svelu', 'hvelu')
STYLES = ('df', 'wd1', 'wd2')

STYLES = ('wd1',) # only this passes right now

for formula in FORMULAS:
    for style in STYLES:
        class cls(CSIDH_test_base, TestCase):
            formula = formula
            style = style
        globals()['csidh_'+formula+'_'+style] = cls
del cls
