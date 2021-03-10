from unittest import TestCase, main
from sidh.csidh import CSIDH

class TestCSIDH(TestCase):

#    def setUp(self):

    def test_genpubvalidate(self):
        c = CSIDH('montgomery', 'p512', 'hvelu', 'wd1', False, 2)
        sk = c.gae.random_key()
        pk = c.gae.pubkey(sk)
        self.assertTrue(c.curve.validate(pk))

