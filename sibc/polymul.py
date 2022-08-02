from sibc.constants import parameters
# Next library is used for printing polynomials
from math import floor, sqrt, log

from sibc.common import attrdict
from functools import reduce
from copy import copy as pythoncopy

def PolyMul(field, maxdeg = None, mindeg = 64):

    # ----------------------------------------------------------------------------------

    try:
        copy = field.copy
    except:
        copy = pythoncopy

    # Table of 2 raised to a negative power
    if maxdeg != None:
        inverse_of_two = (field(2) ** -1)
        max_exp = 2 * int(floor(sqrt(maxdeg)))
        negative_powers_of_two = dict()
        negative_powers_of_two[1] = 1

        j = 1
        j_next = 1
        for i in range(0, max_exp, 1):
            j_next *= 2
            negative_powers_of_two[j_next] = (
                negative_powers_of_two[j] * inverse_of_two
            )
            j = j_next

    # ----------------------------------------------------------------------------------
    # Here, a polynomial is represented as a list of integers.
    def print_poly(h, lenh):

        # Just a pretty way for printing polynomials
        print(reduce( lambda x,y: x + ' + ' + y, [ f'({h[i]}) * (x ^ {i})' for i in range(0, lenh, 1)]))
        return None

    def karatsuba_mul(f, flen, g, glen):
        """
        Karatsuba style multiplication of two polynomials
        """

        if flen < glen:
            return poly_mul(g, glen, f, flen)

        # At this step, we ensure flen >= glen
        if [] == g or [] == f:

            # Multiplication by 0, here the zero polynomial is represented as the empty list
            return []

        if glen == 1:
            # XXX 54k
            # Multipication by a constant
            return [(f[i] * g[0]) for i in range(0, flen, 1)]

        # At this step, we ensure flen >= glen >= 2
        if flen == 2:
            # XXX 57k
            # Multiplication of linear polynomials over fp
            # Thus c(x) = (a * b)(x) is a quadratic polynomial
            c = [0, 0, 0]

            c[0] = (f[0] * g[0])  # coeff of x^0
            c[2] = (f[1] * g[1])  # coeff of x^2
            f01 = (f[0] + f[1])
            g01 = (g[0] + g[1])
            c[1] = (f01 * g01)
            c[1] -= c[0]
            c[1] -= c[2]  # coeff of x^1
            return c

        if flen == 3:

            # f(x) is a quadratic polynomial
            if glen == 2:

                # g(x) is a linear polynomial
                c = list(
                    [0, 0, 0, 0]
                )  # Thus c(x) = (f * g)(x) is a cubic polynomial

                c[0] = (f[0] * g[0])  # coeff of x^0
                c[2] = (f[1] * g[1])
                f01 = (f[0] + f[1])
                g01 = (g[0] + g[1])
                c[1] = (f01 * g01)
                c[1] -= c[0]
                c[1] -= c[2]  # coeff of x^1
                c[3] = (f[2] * g[1])  # coeff of x^3
                f2g0 = (f[2] * g[0])
                c[2] += f2g0  # coeff of x^2
                return c

            if glen == 3:
                # XXX 24k but recursion
                # g(x) is a a quadratic polynomial

                karatsuba_0 = karatsuba_mul(
                    f[:1], 1, g[:1], 1
                )  # field multiplication
                karatsuba_1 = karatsuba_mul(
                    f[1:flen], 2, g[1:glen], 2
                )  # multiplication of two linear polynomials

                # Middle part
                f_01 = (f[0] + f[1])
                f_02 = (f[0] + f[2])
                g_01 = (g[0] + g[1])
                g_02 = (g[0] + g[2])

                f_01 *= g_01 # t_01
                f_02 *= g_02 # t_02

                f_01 -= karatsuba_0[0] # t_01 - karat # l_coeff
                f_01 -= karatsuba_1[0] # l_coeff
                f_02 -= karatsuba_0[0] # t_02 - # q_coeff
                f_02 -= karatsuba_1[2]
                f_02 += karatsuba_1[0] # we don't use karatsuba_1[0] XXX

                return list(karatsuba_0 + [f_01, f_02] + karatsuba_1[1:])

        else:
            # XXX 20k but worth looking at
            # Multiplication of polynomials of degree >= 3 over fp
            nf = flen // 2
            mf = flen - nf
            f_low = f[:nf]
            f_high = f[nf:flen]

            if glen <= nf:
                # XXX 700 neglible
                c0 = karatsuba_mul(f_low, nf, g[:glen], glen)
                c1 = karatsuba_mul(f_high, mf, g[:glen], glen)
                return (
                    c0[:nf]
                    + [
                        (c0[nf + i] + c1[i])
                        for i in range(0, glen - 1, 1)
                    ]
                    + c1[(glen - 1) :]
                )

            mg = glen - nf
            g_low = g[:nf]
            g_high = g[nf:glen]

            f_mid = [
                (f_low[i] + f_high[i]) for i in range(0, nf, 1)
            ] + f_high[nf:]

            fg_low = karatsuba_mul(f_low, nf, g_low, nf)
            fg_high = karatsuba_mul(f_high, mf, g_high, mg)

            if mg < mf:

                g_mid = [
                    (g_low[i] + g_high[i]) for i in range(0, mg, 1)
                ] + g_low[mg:]
                fg_mid = poly_mul(f_mid, mf, g_mid, nf)

            else:

                g_mid = [
                    (g_low[i] + g_high[i]) for i in range(0, nf, 1)
                ] + g_high[nf:]
                fg_mid = poly_mul(f_mid, mf, g_mid, mg)

            for i in range(0, mf + mg - 1, 1):
                fg_mid[i] -= fg_high[i]
            for i in range(0, 2 * nf - 1, 1):
                fg_mid[i] -= fg_low[i]

            ret = []
            ret += fg_low[:nf]
            for i in range(0, nf - 1, 1):
                fg_low[nf + i] += fg_mid[i]
                ret.append(fg_low[nf + i])
            ret.append(fg_mid[nf - 1])
            for i in range(0, mf - 1):
                fg_mid[nf + i] += fg_high[i]
                ret.append(fg_mid[nf + i])
            ret += fg_high[(mf - 1) :]

            return ret

    def qring_mul(f, g, e):

        n = len(f)
        m = len(f[0])

        if (n == 1) and (m <= 8):

            # The product f[1] * g[1] is directly computed in B = fp[x] / (x^m + 1)
            h = karatsuba_mul(f[0], m, g[0], m)  # Polynomial multiplication

            # Next, we reduce h modulo (x^m + 1)
            h_0 = h[:m]
            h_1 = h[m:]
            h_1 += [0] * (len(h_0) - len(h_1))

            return [[(h_0[i] - h_1[i]) for i in range(0, m, 1)]]

        elif (n == 1) and (m > 8):

            # We need to  embed f[1] and g[1] in B[y]/(y^n2+1) being
            # B = Fp[x]/(x^m2+1) where n2, andm2 are computed as follows
            s = int(floor(log(m, 2))) // 2
            m2 = 2 ** (s + 1)
            n2 = (2 * m) // m2
            hm2 = m2 // 2  # notice that n2*hm2 = m

            # Now, F and G will corresponds with the  images of f[1] and g[1],
            # respectively, in B2[y]/(y^n2+1); here, B2 = Fp[x]/(x^hm2+1).
            # Notice B2[y]/(y^n2+1) is 'contained' in B[y]/(y^n2+1)
            F = []
            G = []
            ind = -1

            for j in range(0, n2, 1):

                Fj = []
                Gj = []

                for i in range(0, hm2, 1):
                    ind += 1
                    Fj.append(f[0][ind])
                    Gj.append(g[0][ind])

                for i in range(0, hm2, 1):
                    Fj.append(0)
                    Gj.append(0)

                F.append(Fj)
                G.append(Gj)

            # Next, we recursively multiply F and G but now in B[y]/(y^n2+1) = B[y]/(y^n2-x^m2)
            H = qring_mul(F, G, m2)
            # At this point, H is equal to  n2*F*G, so we need to divide by n2
            # in2 =  fp.fp_inv(n2)   # This inverse can be omited (we only care about polynomials with the same roots)

            # H lives in B[y]/(y^n2+1), 'contained' in Fp[x,y]/(y^n2+1), so to recover
            # the corresponding element h in Fp[x][x^m+1], in H we substitute y by x^hm2.
            h = list([0] * m)
            degy = 0
            func = {0:(lambda x,y: x+y), 1:lambda x,y: x-y}
            for j in range(0, n2, 1):

                deg = degy
                for i in range(0, m2, 1):

                    # for a polynomial position (l,k) (that is position [j-1,i-1] when viewed
                    # as sequences), deg is the degree of x^k*y^l after substitution, namely,
                    # deg = l*hm2 + k. Sometimes we have deg>m-1, this is why case we have to
                    # take residues (deg mod m) and also flip the sign of H(l,k) when deg>m-1
                    q, r = divmod(deg, m)
                    # (-1)^q determines if  or  will be needed
                    hr = H[j][i]
                    h[r] = func[q % 2](h[r], hr)
                    deg += 1

                degy = degy + hm2

            return [
                [
                    (h[i] * negative_powers_of_two[n2])
                    for i in range(0, m, 1)
                ]
            ]

        else:

            # The strategy proceeds by spliting the degree-n problem over B=Fp[x]/(x^m+1) into
            # two degree-n/2 over B, and then each of the two problems is handled recursively
            # until reaching degree-0 products over B

            # However, f*g is to be computed in B[y]/(y^n-x^e) = B[y]/[(y^n2)^2-(x^e2)^2]
            n2 = n // 2
            e2 = e // 2

            F1, F2 = [], []
            G1, G2 = [], []
            # (-1)^q determines if an addition or substraction will be needed
            func = {0:(lambda x,y: x+y), 1:lambda x,y: x-y}
            for j in range(0, n2, 1):

                Fj1, Fj2 = [], []
                Gj1, Gj2 = [], []
                for i in range(0, m, 1):

                    q, r = divmod(i - e2, m)
                    sgn = q % 2
                    # ---
                    Fj1.append(func[sgn](f[j][i], f[n2 + j][r]))
                    Gj1.append(func[sgn](g[j][i], g[n2 + j][r]))
                    # ---
                    Fj2.append(func[1 - sgn](f[j][i], f[n2 + j][r]))
                    Gj2.append(func[1 - sgn](g[j][i], g[n2 + j][r]))

                F1.append(Fj1)
                F2.append(Fj2)
                G1.append(Gj1)
                G2.append(Gj2)

            # Next, we recursively multiply F1 and G1, and also F2 and G2
            H1 = qring_mul(F1, G1, e2)
            H2 = qring_mul(F2, G2, m + e2)

            h1, h2 = [], []
            for j in range(0, n2, 1):

                h1j, h2j = [], []
                for i in range(0, m, 1):

                    q, r = divmod(i - m + e2, m)
                    sgn = q % 2
                    h1j.append((H1[j][i] + H2[j][i]))
                    h2j.append(
                        func[1 - sgn](0, (H1[j][r] - H2[j][r]))
                    )

                h1.append(h1j)
                h2.append(h2j)

            return h1 + h2

    if maxdeg != None:
        def poly_mul(f, flen, g, glen):

            ff = list(f) + [0] * (glen - flen)
            gg = list(g) + [0] * (flen - glen)
            # assert(len(ff) == len(gg))
            hlen = len(ff)

            # hlen > 0 implies `fast` multiplication in fp[x]/(x^n + 1) will be always used
            if hlen > mindeg:
                n = 2 ** len(bin(2 * hlen - 1)[2:])
                fg_qring = qring_mul(
                    [ff + [0] * (n - hlen)], [gg + [0] * (n - hlen)], 0
                )

                # return [ (fg_qring[0][i], in2) for i in range(0, flen + glen, 1) ]
                return fg_qring[0][: (flen + glen - 1)]

            else:
                # This branch is only for checking if karatsuba is better(?)
                return karatsuba_mul(f, flen, g, glen)
    else:
        poly_mul = karatsuba_mul # def poly_mul(f, flen, g, glen): return karatsuba_mul(f, flen, g, glen)

    def poly_mul_modxn(n, f, flen, g, glen):
        """
        Next function computes f*g mod x^n but using the same idea from the C-code
        implementation of https://eprint.iacr.org/2020/341
        """

        if flen < glen:
            return poly_mul_modxn(n, g, glen, f, flen)

        # ------------------------------------
        # At this step, we ensure flen >= glen

        # We require f mod x^n and g mod x^n
        if glen > n:
            return poly_mul_modxn(n, f[:n], n, g[:n], n)

        if flen > n:
            return poly_mul_modxn(n, f[:n], n, g, glen)

        if n == 0:
            # This case return zero as the empty list, that is []
            return []

        if [] == g or [] == f:

            # Multiplication by cero (we need to ensure the output has length n, the zero polynomial corresponds a list with n zeroes
            return [0] * n

        if n == 1:
            # Only the constant coefficients are required
            return [(f[0] * g[0])]

        if n >= (flen + glen - 1):

            # Karatsuba multiplication with no possible savings
            return poly_mul(f, flen, g, glen) + [0] * (n - flen - glen + 1)

        # At this step we ensure n <  (flen + glen - 1)
        if glen == 1:

            return [(f[i] * g[0]) for i in range(0, n, 1)]

        # Here, flen >= glen >= 2
        if n == 2:

            # Multiplication of two linear polynomials modulo x^2
            # And thus, the cuadratic coefficient is not required
            f0g0 = (f[0] * g[0])
            f0g1 = (f[0] * g[1])
            #f1g0 = (f[1] * g[0])
            #return [f0g0, (f0g1 + f1g0)]
            f0g1 += (f[1] * g[0]) # f0g0+f1g0
            return [f0g0, f0g1]

        if n == 3:

            if glen == 2:

                # Multiplication modulo x^3 of a linear polynomial with another of degree at least 2
                f0g0 = (f[0] * g[0])
                f1g1 = (f[1] * g[1])

                t01 = (f[0] + f[1])
                t01 *= (g[0] + g[1])
                t01 -= f0g0
                t01 -= f1g1
                #f2g0 = (f[2] * g[0])
                #return [f0g0, t01, (f1g1 + f2g0)]
                f1g1 += (f[2] * g[0]) # f1g1 + f2g0
                return [f0g0, t01, f1g1]

            if glen == 3:
                c00 = (f[0] * g[0])  # coeff of x^0

                c11 = (f[1] * g[1])
                c01 = (f[0] + f[1]) # f01
                c01 *= (g[0] + g[1]) # f01*g01
                c01 -= c00
                c01 -= c11

                c02 = (f[0] + f[2]) # f02
                g02 = (g[0] + g[2])
                c02 *= g02 # f02 * g02
                c02 -= c00
                c02 -= (f[2] * g[2])
                c02 += c11
                return [c00, c01, c02]

        if n == 4 and glen == 4:

            # This special case can be used as the general case it is expensive. Maybe with anoth bound for n should fine
            #S = n
            #for i in range(0, n // 2, 1):
            #    S = S + n - 1 - 2 * i
            # XXX: since n == 4 here, we can precompute this, but with
            # variable n we should use Gauss summation for closed-form.
            # n + sum([n -1 -2*i for i in range(0, n // 2, 1)])
            # 4 + sum([4 -1 -2*i for i in range(0, 4 // 2, 1)])
            # 4 + sum([3,1]) === 8
            S = 8

            c = list([0] * S)
            nf = n // 2  # floor of n/2
            nc = -(-n // 2)  # Ceiling of n/2

            for i in range(0, n, 1):
                c[i] = (f[i] * g[i])
            k = n
            for i in range(0, nf, 1):
                for j in range(0, n - 2 * i - 1, 1):
                    c[k] = (f[i] + f[i + j + 1])  # these overlap and
                    c[k] *= (g[i] + g[i + j + 1]) # could be cached XXX
                    c[k] -= (c[i] + c[i + j + 1])
                    k += 1

            c[n - 1] = copy(c[0]) # XXX
            for i in range(1, nf, 1):
                for j in range(1, n - 2 * i, 1):
                    c[n + 2 * i - 1 + j] += c[(1 + i) * n - i ** 2 - 1 + j]

            for i in range(1, nc, 1):
                c[n + 2 * i - 1] += c[i]

            delta = n - 1
            return c[delta : (n + delta)]

        # The goal is to divide the multiplication using karatsuba style multiplication
        # but saving multiplications (some higher monomials can be omited)
        l1 = n // 2  # floor(n / 2)
        l0 = n - l1  # ceil( n / 2)

        # We proceed by spliting f(x) as f0(x^2) + x * f1(x^2)
        f1len = flen // 2  # floor(flen / 2)
        f0len = flen - f1len  # ceil( flen / 2)
        f0 = [f[2 * i + 0] for i in range(0, f0len, 1)]
        f1 = [f[2 * i + 1] for i in range(0, f1len, 1)]

        # We proceed by spliting g(x) as g0(x^2) + x * g1(x^2)
        g1len = glen // 2  # floor(glen / 2)
        g0len = glen - g1len  # ceil( geln / 2)
        g0 = g[ : 2*g0len :2 ] # [g[2 * i + 0] for i in range(0, g0len, 1)]
        g1 = [g[2 * i + 1] for i in range(0, g1len, 1)]
        # TODO check which of these is faster --^

        # Middle part like karatsuba
        f01 = [(f0[i] + f1[i]) for i in range(0, f1len, 1)] + f0[
            f1len:
        ]
        g01 = [(g0[i] + g1[i]) for i in range(0, g1len, 1)] + g0[
            g1len:
        ]

        # product f0 * g0 must has degree at most x^l0
        n0 = f0len + g0len - 1
        if n0 > l0:
            n0 = l0
        fg_0 = poly_mul_modxn(n0, f0, f0len, g0, g0len)

        # product f01 * g01 must has degree at most x^l1
        n01 = f0len + g0len - 1
        if n01 > l1:
            n01 = l1
        fg_01 = poly_mul_modxn(n01, f01, f0len, g01, g0len)

        # product f1 * g1 must has degree at most x^l1
        n1 = f1len + g1len - 1
        if n1 > l1:
            n1 = l1
        fg_1 = poly_mul_modxn(n1, f1, f1len, g1, g1len)

        # Computing the middle part
        for i in range(n01):
            fg_01[i] -= fg_0[i]
        for i in range(n1):
            fg_01[i] -= fg_1[i]

        # Unifying the computations
        fg = [0]*n
        # fgs are interleaved: [fg0.0, fg01.0, fg0.1, fg01.1]
        fg[:2*n0:2] = fg_0[:n0]
        fg[1:2*n01:2] = fg_01[:n01]

        for i in range(0, n1 - 1, 1):
            fg[2 * i + 2] += fg_1[i]

        if 2 * n1 < n:
            fg[2 * n1] += fg_1[n1 - 1]

        return fg

    def quasi_poly_mul_middle(g, glen, f, flen):
        """
        Next function computes the central part polynomial of f * g mod (x^m - 1)
        Next two functions are based on the work title "The Middle Product
        Algorithm, I." by G. Hanrot. M. Quercia, and P. Zimmermann
        """

        #    if not glen == len(g):
        #        import pdb; pdb.set_trace();
        assert glen == len(g)
        assert flen == len(f)

        if glen == 0:
            return []

        if glen == 1:
            return [(g[0] * f[0])]

        glen0 = glen // 2  # floor(glen / 2)
        glen1 = glen - glen0  # ceil(glen / 2)

        A = poly_mul_middle(
            g[glen0:],
            glen1,
            [
                (f[i] + f[i + glen1])
                for i in range(0, 2 * glen1 - 1, 1)
            ],
            2 * glen1 - 1,
        )

        if glen % 2 == 0:
            B = poly_mul_middle(
                [(g[glen1 + i] - g[i]) for i in range(0, glen0, 1)],
                glen0,
                f[glen1 : (3 * glen1 - 1)],
                2 * glen1 - 1,
            )
        else:
            B = poly_mul_middle(
                [g[glen0]]
                + [(g[glen1 + i] - g[i]) for i in range(0, glen0, 1)],
                glen1,
                f[glen1 : (3 * glen1 - 1)],
                2 * glen1 - 1,
            )

        C = poly_mul_middle(
            g[:glen0],
            glen0,
            [
                (f[glen1 + i] + f[2 * glen1 + i])
                for i in range(0, 2 * glen0 - 1, 1)
            ],
            2 * glen0 - 1,
        )

        assert len(A) == glen1
        assert len(B) == glen1
        assert len(C) == glen0
        return [(A[i] - B[i]) for i in range(0, glen1, 1)] + [
            (C[i] + B[i]) for i in range(0, glen0, 1)
        ]

    def poly_mul_middle(g, glen, f, flen):
        """
        Next function computes the central part polynomial of f * g mod (x^m - 1)
        This functions is an extension of quasi_poly_mul_middle()
        """

        if (glen == (flen // 2)) and (flen % 2 == 0):
            # Here deg(f) = 2d - 1 and deg(g) = d - 1
            return quasi_poly_mul_middle(g, glen, f[1:] + [0], flen)

        elif (glen == ((flen // 2) + 1)) and (flen % 2 == 0):
            # Here deg(f) = 2d - 1 and deg(g) = d
            return quasi_poly_mul_middle(g, glen, [0] + f, flen + 1)

        elif (glen == ((flen // 2) + 1)) and (flen % 2 == 1):
            # Here deg(f) = 2d and deg(g) = d
            return quasi_poly_mul_middle(g, glen, f, flen)

        else:

            # Unbalanced case when glen < floor(flen/2) and ceil(flen/2) < glen
            # This case only happens when using scaled remainder trees
            if glen == 0:
                return []

            # low part
            F0len = flen - glen
            F0 = f[:F0len]
            F0 = F0[::-1]
            G = g[::-1]
            fg_low = poly_mul_modxn(glen - 1, F0, F0len, G, glen)
            fg_low = fg_low[::-1]

            # high part
            F1 = f[F0len:]
            fg_high = poly_mul_modxn(glen, F1, glen, g, glen)

            return [
                (fg_low[i] + fg_high[i]) for i in range(0, glen - 1, 1)
            ] + [fg_high[glen - 1]]

    def poly_mul_selfreciprocal(g, glen, f, flen):
        """
        Next function computes the product of two polynomials g(x) annd f(x) such that
        g(x) = x^dg * g(1/x) and f(x) = x^df * f(1/x) where deg g = dg, and deg f = df
        """
        if glen == 0 and flen == 0:
            return []

        if glen == 1 and flen == 1:
            return [(g[0] * f[0])]

        if glen == 2 and flen == 2:

            h = [0] * 3
            h[0] = (g[0] * f[0])
            h[1] = (h[0] + h[0])
            h[2] = h[0]
            return h

        if glen == 3 and flen == 3:

            h = [0] * 5
            h[0] = (g[0] * f[0])
            h[2] = (g[1] * f[1])
            h[1] = (g[0] + g[1]) # g01
            h[1] *=(f[0] + f[1]) # g01f01
            h[2] += h[0]
            h[1] -= h[2]
            h[2] += h[0]
            h[3] = h[1]
            h[4] = h[0]
            return h

        if glen == 4 and flen == 4:

            h = [0] * 7
            h[0] = (g[0] * f[0])
            h[3] = (g[1] * f[1])
            h[2] = g[0] + g[1]
            h[2] *= f[0] + f[1]
            h[2] -= h[0]
            h[1] = (h[2] - h[3])
            h[3] += h[0]
            h[3] += h[3]
            h[4] = h[2]
            h[5] = h[1]
            h[6] = h[0]
            return h

        if glen == 5 and flen == 5:
            # XXX 4k, not an issue.
            h = [0] * 9
            g10 = (g[1] - g[0])
            f01 = (f[0] - f[1])
            h[1] = (g10 * f01)
            g20 = (g[2] - g[0])
            f02 = (f[0] - f[2])
            h[2] = (g20 * f02)
            g21 = (g[2] - g[1])
            f12 = (f[1] - f[2])
            h[3] = (g21 * f12)

            h[0] = (g[0] * f[0])
            g1f1 = (g[1] * f[1])
            g2f2 = (g[2] * f[2])

            t = (g1f1 + h[0])
            h[1] += t
            h[3] += h[1]
            h[4] = (t + g2f2)
            h[4] += t
            h[2] += t
            h[2] += g2f2
            h[3] += g1f1
            h[3] += g2f2

            h[5] = h[3]
            h[6] = h[2]
            h[7] = h[1]
            h[8] = h[0]
            return h

        copy=g[0].__class__ # caching the class constructor seems to be the most efficient way to copy these.

        if glen == flen and (glen % 2) == 1:

            # Case two degree polynomials of same even degree
            len0 = (glen + 1) // 2
            len1 = glen // 2

            g0 = [0] * len0
            f0 = [0] * len0
            g1 = [0] * len1
            f1 = [0] * len1

            # We proceed by applying the same idea as in poly_mul_modxn
            # ---
            for i in range(0, len0, 1):
                g0[i] = copy(g[2 * i])
                f0[i] = copy(f[2 * i])

            h0 = poly_mul_selfreciprocal(g0, len0, f0, len0)

            # ---
            for i in range(0, len1, 1):
                g1[i] = copy(g[2 * i + 1])
                f1[i] = copy(f[2 * i + 1])

            h1 = poly_mul_selfreciprocal(g1, len1, f1, len1)

            # ---
            for i in range(0, len1, 1):
                g0[i] += g1[i]
                f0[i] += f1[i]

            for i in range(0, len1, 1):
                g0[i + 1] += g1[i]
                f0[i + 1] += f1[i]

            h01 = poly_mul_selfreciprocal(g0, len0, f0, len0)

            # Mixing the computations
            for i in range(0, 2 * len0 - 1, 1):
                h01[i] = h01[i] - h0[i] # h01[i] -= h0[i]

            for i in range(0, 2 * len1 - 1, 1):
                h01[i] = h01[i] - h1[i] # h01[i] -= h1[i]

            for i in range(0, 2 * len1 - 1, 1):
                h01[i + 1] -= h1[i]

            for i in range(0, 2 * len1 - 1, 1):
                h01[i + 1] -= h1[i]

            for i in range(0, 2 * len1 - 1, 1):
                h01[i + 2] -= h1[i]

            for i in range(1, 2 * len0 - 1, 1):
                h01[i] = (h01[i] - h01[i - 1])

            # Unifying the computations
            hlen = 2 * glen - 1
            h = [0] * hlen

            for i in range(0, 2 * len0 - 1, 1):
                h[2 * i] = copy(h0[i])

            for i in range(0, 2 * len0 - 2, 1):
                h[2 * i + 1] = copy(h01[i])

            for i in range(0, 2 * len1 - 1, 1):
                h[2 * i + 2] += h1[i]

            return h

        if glen == flen and (glen % 2) == 0:

            # Case two degree polynomials of same odd degree

            half = glen // 2
            h0 = poly_mul(g[:half], half, f[:half], half)
            h1 = poly_mul(g[:half], half, f[half:], half)

            h = [0] * (2 * glen - 1)
            for i in range(0, glen - 1, 1):
                h[i] += h0[i]

            for i in range(0, glen - 1, 1):
                h[2 * glen - 2 - i] += h0[i]

            for i in range(0, glen - 1, 1):
                h[half + i] += h1[i]

            for i in range(0, glen - 1, 1):
                h[glen + half - 2 - i] += h1[i]

            return h

        # General case
        hlen = glen + flen - 1
        m = (hlen + 1) // 2
        h = poly_mul_modxn(m, g, glen, f, flen)
        h += [0] * (hlen - m)
        for i in range(m, hlen, 1):
            h[i] = h[hlen - 1 - i]

        return h

    def product_tree(f, n):
        """
        Product tree of given list of polynomials
        """

        if n == 0:
            return {'left': None, 'right': None, 'poly': [1], 'deg': 0}

        if n == 1:

            # No multiplication is required
            return {
                'left': None,
                'right': None,
                'poly': f[0],
                'deg': len(f[0]) - 1,
            }

        else:

            m = n - (n // 2)
            left = product_tree(f[:m], m)
            right = product_tree(f[m:], n - m)
            return {
                'left': left,
                'right': right,
                'poly': poly_mul(
                    left['poly'],
                    left['deg'] + 1,
                    right['poly'],
                    right['deg'] + 1,
                ),
                'deg': left['deg'] + right['deg'],
            }

    def product_selfreciprocal_tree(f, n):
        """
        Product tree of given list of self reciprocal polynomials
        """

        if n == 0:
            return {'left': None, 'right': None, 'poly': [1], 'deg': 0}

        if n == 1:

            # No multiplication is required
            return {
                'left': None,
                'right': None,
                'poly': f[0],
                'deg': len(f[0]) - 1,
            }

        else:

            m = n - (n // 2)
            left = product_selfreciprocal_tree(f[:m], m)
            right = product_selfreciprocal_tree(f[m:], n - m)
            return {
                'left': left,
                'right': right,
                'poly': poly_mul_selfreciprocal(
                    left['poly'],
                    left['deg'] + 1,
                    right['poly'],
                    right['deg'] + 1,
                ),
                'deg': left['deg'] + right['deg'],
            }

            # Now, the resultant of two polinomials corresponds with product of all g(x) mod f_i(x) when each f_i is a linear polynomial

    def product(list_g_mod_f, n):

        # When finishing tests, remove the below asserts
        if n == 0:
            return 1

        assert len(list_g_mod_f) == n

        out = copy(list_g_mod_f[0][0])
        for j in range(1, n, 1):

            out *= list_g_mod_f[j][0]

        return out

    return attrdict(locals())
