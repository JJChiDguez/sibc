from unittest import TestCase
from sidh.csidh import CSIDH
from sidh.common import attrdict


class known_df_p512(object):

    """
    Known-answer tests

    We currently only have these for p512 df keys.

    These keys and answers came from running csidh/main.py (which generated
    these keys) prior to the big refactoring.

    To run this test for just one formula, you can run commands like this:
        pytest -k 'known_df and hvelu'
    """

    keys = attrdict(
        A=attrdict(
            sk=[
                # fmt: off
                9, 2, -14, 3, 11, -12, 4, 4, 0, 20, 8, 7, -12, -20, 23, 15, 5,
                3, 15, -19, 7, -17, -19, 1, -2, 14, 0, -9, 2, 4, 11, 2, 7, 9,
                9, -1, 5, -7, 5, 4, -4, 6, 6, 7, -8, -2, 0, 2, -6, 5, -2, 0,
                -2, 4, -5, -1, -5, 3, 3, 5, 3, -5, 5, 3, -2, 2, -4, -2, 0, -2,
                2, 0, 2, -3,
                # fmt: on
            ],
            pk=[
                0x2B84079DFE34DACA1000B23DE50ED8E581C2F8F174FDBCB03D3CA6C96361E731152F45BDD4959832DE406E8F3C4F0B4949C4826AF82B3A7362677A458196FBCF,
                0x76599305D04FB32F8907B2279EE35BE99786DA7C055E4A922712A6B546D457FA8DB9529049BBE500BE47E23DAE04ECD34E043264A02BB1917DFDF9FA540F233,
            ],
            compressed=0x27AD85DDC08BF510F08A8562BA4909803675536A0BCE6250E3BED4A9401AC123FE75C18866625E9FCFCAF03D0927ED46665E153E786244DAAAC9A83075060C82,
        ),
        B=attrdict(
            sk=[
                # fmt: off
                3, -16, -10, 1, 15, 20, -20, -22, -16, -22, 0, -19, 6, -4, -9,
                13, -11, 13, -13, -1, 23, 21, -5, 13, -4, -2, 12, 15, -4, -10,
                -5, 0, 11, 1, -1, -1, 7, 1, -3, 6, 0, 2, -4, -5, 0, 2, -4, -2,
                -4, -5, 6, 2, -6, -4, 5, -5, 5, -3, 1, 3, -1, -5, 3, -5, -4, 2,
                4, 2, 2, 4, 0, -2, 0, -3,
                # fmt: on
            ],
            pk=[
                0x5E2FD48334E49B47BB754F88B3345C77D604EB1FADC29B35B931724459143ABDE2F22346B595E3B161D80F3659870F16D4983BFA58F5F2F9718D3B375C21D65C,
                0x314B346A927C21051052B790809D895627ED8FBE4F008408D361223A97556EC0E6D3B544B0898DAFFCDBFF5C5B409CCB5CC9E2EDC95504FCA54318071E28E054,
            ],
            compressed=0x181C39753CCB4D3358E32B4471EE73EDC568846CA3B0B17571A09BD7373B4658251ADF466FF1FFB29D89B382184703C708F71497611A4B643BD984D847F3A430,
        ),
    )

    ss = 0x1ADB783878BA330BB2A842E7F8B3392329A2CD3B407900E4CF6A8F13B744BFFEFF617BDE2CEBBB9CE97D32BC6FC1BCE2D88381B03B3E13CFF0651EEA82D02937

    prime = 'p512'
    formula = NotImplemented
    style = 'df'
    tuned = False
    exponent = 2
    multievaluation = False
    verbose = False

    def setUp(self):
        self.csidh = CSIDH(
            'montgomery',
            self.prime,
            self.formula,
            self.style,
            self.tuned,
            self.exponent,
            self.multievaluation,
            self.verbose,
        )

        self.coeff = self.csidh.curve.coeff

    def test_dh_AB(self):
        self.assertEqual(
            self.coeff(self.csidh.gae.dh(self.keys.A.sk, self.keys.B.pk)),
            self.ss,
        )

    def test_dh_BA(self):
        self.assertEqual(
            self.coeff(self.csidh.gae.dh(self.keys.B.sk, self.keys.A.pk)),
            self.ss,
        )

    def test_compress(self):
        for keys in self.keys.values():
            self.assertEqual(keys.compressed, self.csidh.curve.coeff(keys.pk))

    def test_compress_roundtrip(self):
        compress, uncompress = (
            self.csidh.curve.coeff,
            self.csidh.curve.affine_to_projective,
        )
        for keys in self.keys.values():
            self.assertEqual(
                keys.compressed, compress(uncompress(compress(keys.pk)))
            )


class Test_known_df_p512_hvelu(known_df_p512, TestCase):
    formula = 'hvelu'


class Test_known_df_p512_tvelu(known_df_p512, TestCase):
    formula = 'tvelu'


class Test_known_df_p512_svelu(known_df_p512, TestCase):
    formula = 'svelu'
