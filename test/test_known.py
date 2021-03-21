from unittest import TestCase
from sibc.csidh import CSIDH
from sibc.common import attrdict


class known_df_p512(object):

    """
    Known-answer tests

    We currently only have these for p512 df keys.

    These keys and answers came from running csidh/main.py from 
    https://github.com/JJChiDguez/csidh_withstrategies (which
    generated these keys).

    To run this test for just one formula, you can run commands like this:
        pytest -k 'known_df and hvelu'
    """

    keys = attrdict(
        A=attrdict(
            sk=[7, -14, -16, -5, 7, -4, 6, -16, 14, 22, -10, 5, 16, 8, -17, -1, -1, -11, 5, 11, 1, -1, 7, -5, 10, -2, 14, -13, 12, -4, -11, 10, -9, -9, 5, 9, 1, 9, 9, 0, 6, -4, 6, -5, 4, -6, 2, 4, 6, 3, 0, -2, 0, -6, -5, -3, 5, 5, -3, 5, 1, -1, 3, -1, 4, -2, 0, 0, 2, 4, 0, 4, 4, 1],
            pk=[4728146032486912801094583991353286261973546366143541397245686753762901316062339836438041892681193229023769378118684082801187071058443826977747221505877393, 2983658557129372264342771993184706542050709103347982710108512572730239907028070100828017140583214399358736276490324811325106128107134002761476490169757932],
            compressed=0x3BD2E2F7B2EA49932420DAF81189C350519FF0E29589587C0344111D4822C80AD88100DAC094898DE9FF615A61929B2A6D23C1C149CB031E47AE3BC4E4DC6A19,
        ),
        B=attrdict(
            sk=[3, 14, 18, -15, 17, 22, -12, -14, -20, 10, 10, 3, 10, -8, -17, 13, -15, 19, -19, -21, -15, 23, 11, 1, -2, 10, 6, -3, 0, 6, 9, 12, -9, -9, 3, 5, -1, -9, 7, -2, -8, 8, 0, -3, 4, -2, 0, -2, -6, -7, 0, 2, 6, -4, -5, 3, 5, -5, 1, -3, -5, -1, -1, -3, 2, 4, 0, 0, 4, 2, 2, 4, -2, 3],
            pk=[3960926043577541377267819439096804391618603251899293401108425870931130716835203327490747655445345747395550952611253698265199934334080097314192061511711009, 487723576446903129146235324971287599085008588633854300660046476997824105838685867154322957184759627700657558727262867917195404714619211343006256912881095],
            compressed=0x3E533EEDF73446A8FCDAE46A97E4855C62E7A67B6F5CBE2710087C8FFC3C63A4BFA6B365553FB23219955B13A99207721F936E790795CB0AE09D206655197063,
        ),
    )

    ss = 0x13F022ADE3CC59641C069BF9FE798A6AACE097DF451EBDF95DC016978E494748BD9935F311F4D403AFBC96021FD19EBB7381D077D2A1A8D59FC81F97508B4F5A

    prime = 'p512'
    formula = NotImplemented
    style = 'df'
    tuned = False
    exponent = 10
    multievaluation = False
    uninitialized = False
    verbose = False

    def setUp(self):
        self.csidh = CSIDH(
            'montgomery',
            prime=self.prime,
            formula=self.formula,
            style=self.style,
            exponent=self.exponent,
            tuned=self.tuned,
            uninitialized=self.uninitialized,
            multievaluation=self.multievaluation,
            verbose=self.verbose,
        )

        self.coeff = self.csidh.curve.coeff
        self.field = self.csidh.field

    def test_dh_AB(self):
        self.assertEqual(
            self.coeff(self.csidh.gae.GAE_at_A(self.keys.A.sk, [self.field(pk) for pk in self.keys.B.pk])).x,
            self.ss,
        )

    def test_dh_BA(self):
        self.assertEqual(
            self.coeff(self.csidh.gae.GAE_at_A(self.keys.B.sk, [self.field(pk) for pk in self.keys.A.pk])).x,
            self.ss,
        )

    def test_compress(self):
        for keys in self.keys.values():
            self.assertEqual(keys.compressed, self.csidh.curve.coeff([self.field(pk) for pk in keys.pk]).x)

    def test_compress_roundtrip(self):
        compress, uncompress = (
            self.csidh.curve.coeff,
            self.csidh.curve.affine_to_projective,
        )
        for keys in self.keys.values():
            self.assertEqual(
                keys.compressed, compress(uncompress(compress([self.field(pk) for pk in keys.pk]))).x
            )


class Test_known_df_p512_hvelu(known_df_p512, TestCase):
    formula = 'hvelu'


class Test_known_df_p512_tvelu(known_df_p512, TestCase):
    formula = 'tvelu'


class Test_known_df_p512_svelu(known_df_p512, TestCase):
    formula = 'svelu'
