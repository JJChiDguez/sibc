#x = dict(df=gae_df, wd1=gae_wd1, wd2=gae_wd2)
from sidh.csidh._gae_df import Gae_df
from sidh.csidh._hvelu import Hvelu
from sidh.csidh._tvelu import Tvelu
from sidh.csidh._montgomery import MontgomeryLadder

class CSIDH(object):
    """

    CSIDH

    >>> sk_a = [ 9,  2,-14,  3, 11,-12,  4,  4,  0, 20,  8,  7,-12,-20, 23, 15,  5,  3, 15,-19,  7,-17,-19,  1, -2, 14,  0, -9,  2,  4, 11,  2,  7,  9,  9, -1,  5, -7,  5,  4, -4,  6,  6,  7, -8, -2,  0,  2, -6,  5, -2,  0, -2,  4, -5, -1, -5,  3,  3,  5,  3, -5,  5,  3, -2,  2, -4, -2,  0, -2,  2,  0, 2, -3 ]
    >>> pk_a = [0x2b84079dfe34daca1000b23de50ed8e581c2f8f174fdbcb03d3ca6c96361e731152f45bdd4959832de406e8f3c4f0b4949c4826af82b3a7362677a458196fbcf, 0x76599305d04fb32f8907b2279ee35be99786da7c055e4a922712a6b546d457fa8db9529049bbe500be47e23dae04ecd34e043264a02bb1917dfdf9fa540f233]
    >>> sk_b = [ 3,-16,-10,  1, 15, 20,-20,-22,-16,-22,  0,-19,  6, -4, -9, 13,-11, 13,-13, -1, 23, 21, -5, 13, -4, -2, 12, 15, -4,-10, -5,  0, 11,  1, -1, -1,  7,  1, -3,  6,  0,  2, -4, -5,  0,  2, -4, -2, -4, -5,  6,  2, -6, -4,  5, -5,  5, -3,  1,  3, -1, -5,  3, -5, -4,  2,  4,  2,  2,  4,  0, -2, 0, -3 ]
    >>> pk_b = [ 0x5e2fd48334e49b47bb754f88b3345c77d604eb1fadc29b35b931724459143abde2f22346b595e3b161d80f3659870f16d4983bfa58f5f2f9718d3b375c21d65c, 0x314b346a927c21051052b790809d895627ed8fbe4f008408d361223a97556ec0e6d3b544b0898daffcdbff5c5b409ccb5cc9e2edc95504fca54318071e28e054 ]
    >>> ss = 0x1ADB783878BA330BB2A842E7F8B3392329A2CD3B407900E4CF6A8F13B744BFFEFF617BDE2CEBBB9CE97D32BC6FC1BCE2D88381B03B3E13CFF0651EEA82D02937
    >>> csidh = CSIDH('montgomery', 'p512', 'tvelu', 'df', False, 2)
    // Shortest Differential Addition Chains (SDAC) for each l_i;
    // SDAC's to be read from a file
    >>> csidh.dh(sk_a, pk_b) == csidh.dh(sk_b, pk_a) == ss
    True
    >>> coeff, affine_to_projective = csidh.curve.coeff, csidh.curve.affine_to_projective
    >>> hex(coeff(pk_a)).upper()[2:]
    '27AD85DDC08BF510F08A8562BA4909803675536A0BCE6250E3BED4A9401AC123FE75C18866625E9FCFCAF03D0927ED46665E153E786244DAAAC9A83075060C82'
    >>> hex(coeff(pk_b)).upper()[2:]
    '181C39753CCB4D3358E32B4471EE73EDC568846CA3B0B17571A09BD7373B4658251ADF466FF1FFB29D89B382184703C708F71497611A4B643BD984D847F3A430'
    >>> coeff(affine_to_projective(coeff(pk_a))) == coeff(pk_a)
    True
    >>> csidh.dh(sk_b, affine_to_projective(coeff(pk_a))) == ss
    True
    >>> list(map(hex,affine_to_projective(coeff(pk_a))))
    ['0x27ad85ddc08bf510f08a8562ba4909803675536a0bce6250e3bed4a9401ac123fe75c18866625e9fcfcaf03d0927ed46665e153e786244daaac9a83075060c84', '0x4']
    """

    def __init__(self, curvemodel, prime, formula, style, verbose, exponent):
        self.curvemodel = curvemodel
        self.prime = prime
        self.style = style
        self._exponent = exponent
        self.fp = None

        # Where do we do our math? On a curve!
        if self.curvemodel == 'montgomery':
            self.curve = MontgomeryLadder(prime, style)
            self.fp = self.curve.fp
        else:
            self.curve = None
            raise NotImplemented

        # Formulas for ???
        if formula == 'hvelu':
            self.formula = Hvelu(self.curve, verbose)
        elif formula == 'tvelu':
            self.formula = Tvelu(self.curve)
        else:
            self.formula = NotImplemented

        # Side channel protection styles
        if self.style == 'df':
            self.gae = Gae_df(prime, style, verbose, self.curve, self.formula)
        else:
            self.gae = None

    def dh(self, sk, pk):
        return self.gae.dh(sk, pk)

    def pubkey(self, sk):
        return self.gae.pubkey(sk)

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
