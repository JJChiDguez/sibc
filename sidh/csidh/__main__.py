from logging import getLogger
import click
from click.exceptions import Exit
from sidh.framework import *
from .constants import parameters


#        from sidh.csidh.gae_wd1 import *
#            from sidh.csidh.tvelu import *
#                from sidh.csidh.poly_redc import *
#                    from sidh.csidh.poly_mul import *
#                        from sidh.csidh.montgomery import *
#                            from sidh.fp import *
#                        from sympy import *


def csidh_dh(GAE, L, m, pk, sk):
    """
    >>> sk_a = [  9,  2,-14,  3, 11,-12,  4,  4,  0, 20,  8,  7,-12,-20, 23, 15, 5,  3, 15,-19,  7,-17,-19,  1, -2, 14,  0, -9,  2,  4, 11,  2,  7,  9,  9, -1,  5, -7,  5,  4, -4,  6,  6,  7, -8, -2,  0,  2, -6,  5, -2,  0, -2,  4, -5, -1, -5,  3,  3,  5,  3, -5,  5,  3, -2,  2, -4, -2,  0, -2,  2,  0, 2, -3 ]
    >>> pk_a = 0x27AD85DDC08BF510F08A8562BA4909803675536A0BCE6250E3BED4A9401AC123FE75C18866625E9FCFCAF03D0927ED46665E153E786244DAAAC9A83075060C82
    >>> sk_b = [  3,-16,-10,  1, 15, 20,-20,-22,-16,-22,  0,-19,  6, -4, -9, 13,-11, 13,-13, -1, 23, 21, -5, 13, -4, -2, 12, 15, -4,-10, -5,  0, 11,  1, -1, -1,  7,  1, -3,  6,  0,  2, -4, -5,  0,  2, -4, -2, -4, -5,  6,  2, -6, -4,  5, -5,  5, -3,  1,  3, -1, -5, 3, -5, -4,  2,  4,  2,  2,  4,  0, -2, 0, -3 ]
    >>> pk_b = 0x181C39753CCB4D3358E32B4471EE73EDC568846CA3B0B17571A09BD7373B4658251ADF466FF1FFB29D89B382184703C708F71497611A4B643BD984D847F3A430
    >>> ss = 0x1ADB783878BA330BB2A842E7F8B3392329A2CD3B407900E4CF6A8F13B744BFFEFF617BDE2CEBBB9CE97D32BC6FC1BCE2D88381B03B3E13CFF0651EEA82D02937
    >>> csidh_dh(sk_a, pk_b) == csidh(sk_b, pk_a) == ss
    True
    """
    C_out, L_out, R_out, S_out, r_out = strategy_block_cost(L[::-1], m[::-1])

def get_strategy_params(prime, style):

    else:
        delta = None
    return delta, m, sigma, kappa

def load_strategy(m, n):
    temporal_m = list(set(m))
    if len(temporal_m) > 1:
        # Maximum number of degree-(l_i) isogeny constructions is m_i (different for each l_i)
        LABEL_m = 'different_bounds'
    else:
        # Maximum number of degree-(l_i) isogeny constructions is m (the same for each l_i)
        LABEL_m = 'with_same_bounds'

    if setting.verbose:
        verb = '-suitable'
    else:
        verb = '-classical'

    try:

        # List of Small Odd Primes, L := [l_0, ..., l_{n-1}]
        m_prime = [geometric_serie(m[k], L[k]) for k in range(n)]
        r_out, L_out, R_out = rounds(m_prime[::-1], n)
        for j in range(0, len(r_out), 1):

            R_out[j] = list([L[::-1][k] for k in R_out[j]])
            L_out[j] = list([L[::-1][k] for k in L_out[j]])

        f = open(
            strategy_data
            + setting.algorithm
            + '-'
            + setting.prime
            + '-'
            + setting.style
            + '-'
            + setting.formulaes
            + '-'
            + LABEL_m
            + verb
        )
        print("// Strategies to be read from a file")
        S_out = []
        for i in range(0, len(r_out), 1):

            tmp = f.readline()
            tmp = [int(b) for b in tmp.split()]
            S_out.append(tmp)

        f.close()

    except IOError:
        compute_strategy()

def compute_strategy():
    print("// Strategies to be computed")
    C_out, L_out, R_out, S_out, r_out = strategy_block_cost(L[::-1], m[::-1])
    f = open(
        strategy_data
        + setting.algorithm
        + '-'
        + setting.prime
        + '-'
        + setting.style
        + '-'
        + setting.formulaes
        + '-'
        + LABEL_m
        + verb,
        'w',
    )
    for i in range(0, len(r_out)):

        f.writelines(' '.join([str(tmp) for tmp in S_out[i]]) + '\n')

    f.close()

    print(
        "// All the experiments are assuming S = %1.6f x M and a = %1.6f x M. The measures are given in millions of field operations.\n"
        % (SQR, ADD)
    )


    ''' -------------------------------------------------------------------------------------
        Framework
        ------------------------------------------------------------------------------------- '''

    print("p := 0x%X;" % p)
    print("fp := GF(p);")
    print("P<x> := PolynomialRing(fp);")
    print("fp2<i> := ext<fp | x^2 + 1>;")
    print("P<x> := PolynomialRing(fp2);")

    A = [2, 4]
    print("public_coeff := 0x%X;\n" % coeff(A))

    ''' -------------------------------------------------------------------------------------
        Main
        ------------------------------------------------------------------------------------- '''

    print("// Maximum number of degree-(\ell_i) isogeny constructions: m_i")
    print("/*")
    printl("m", m, n // 3)
    print("*/")
    print(
        "// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )
    print("// Public Key Generation")
    print(
        "// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
    )

    # ------------------------------------------------------------------------- Alice
    set_zero_ops()
    public_validation = validate(A)
    assert public_validation

    #a_private = random_key(m)
    a_private = [ 9,  2,-14,  3, 11,-12,  4,  4,  0, 20,  8,  7,-12,-20, 23, 15,
            5,  3, 15,-19,  7,-17,-19,  1, -2, 14,  0, -9,  2,  4, 11,  2,  7,  9,
            9, -1,  5, -7,  5,  4, -4,  6,  6,  7, -8, -2,  0,  2, -6,  5, -2,  0,
            -2,  4, -5, -1, -5,  3,  3,  5,  3, -5,  5,  3, -2,  2, -4, -2,  0, -2,
            2,  0, 2, -3 ]

    print("// Private key corresponding to Alice")
    print("/*")
    printl("sk_a", a_private, n // 3)
    print("*/")

    RUNNING_TIME = get_ops()
    print("// Public key corresponding to Alice")
    print(
        "// public key validation cost  :\t%2.3fM + %2.3fS + %2.3fa = %2.3fM,"
        % (
            RUNNING_TIME[0] / (10.0 ** 6),
            RUNNING_TIME[1] / (10.0 ** 6),
            RUNNING_TIME[2] / (10.0 ** 6),
            measure(RUNNING_TIME) / (10.0 ** 6),
        )
    )

    set_zero_ops()
    if (len(temporal_m) == 1) or ((len(temporal_m) == 2) and (0 in temporal_m)):
        a_public = GAE(
            A, a_private, [L_out[0]], [R_out[0]], [S_out[0]], [temporal_m[-1]], m
        )
    else:
        a_public = GAE(A, a_private, L_out, R_out, S_out, r_out, m)

    RUNNING_TIME = get_ops()
    print(
        "// group action evaluation cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
        % (
            RUNNING_TIME[0] / (10.0 ** 6),
            RUNNING_TIME[1] / (10.0 ** 6),
            RUNNING_TIME[2] / (10.0 ** 6),
            measure(RUNNING_TIME) / (10.0 ** 6),
        )
    )
    print("pk_a := 0x%X;\n" % coeff(a_public))

    # ------------------------------------------------------------------------- Bob
    set_zero_ops()

    #b_private = random_key(m)
    b_private = [ 3,-16,-10,  1, 15, 20,-20,-22,-16,-22,  0,-19,  6, -4, -9,
            13,-11, 13,-13, -1, 23, 21, -5, 13, -4, -2, 12, 15, -4,-10, -5,  0, 11,
            1, -1, -1,  7,  1, -3,  6,  0,  2, -4, -5,  0,  2, -4, -2, -4, -5,  6,
            2, -6, -4,  5, -5,  5, -3,  1,  3, -1, -5,  3, -5, -4,  2,  4,  2,  2,
            4,  0, -2, 0, -3 ]

    print("// Private key corresponding to Bob")
    print("/*")
    printl("sk_b", b_private, n // 3)
    print("*/")

    set_zero_ops()
    public_validation = validate(A)
    assert public_validation
    RUNNING_TIME = get_ops()
    print("// Public key corresponding to Bob")
    print(
        "// public key validation cost  :\t%2.3fM + %2.3fS + %2.3fa = %2.3fM,"
        % (
            RUNNING_TIME[0] / (10.0 ** 6),
            RUNNING_TIME[1] / (10.0 ** 6),
            RUNNING_TIME[2] / (10.0 ** 6),
            measure(RUNNING_TIME) / (10.0 ** 6),
        )
    )

    set_zero_ops()
    if (len(temporal_m) == 1) or ((len(temporal_m) == 2) and (0 in temporal_m)):
        b_public = GAE(
            A, b_private, [L_out[0]], [R_out[0]], [S_out[0]], [temporal_m[-1]], m
        )
    else:
        b_public = GAE(A, b_private, L_out, R_out, S_out, r_out, m)

    RUNNING_TIME = get_ops()
    print(
        "// group action evaluation cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
        % (
            RUNNING_TIME[0] / (10.0 ** 6),
            RUNNING_TIME[1] / (10.0 ** 6),
            RUNNING_TIME[2] / (10.0 ** 6),
            measure(RUNNING_TIME) / (10.0 ** 6),
        )
    )
    print("pk_b := 0x%X;" % coeff(b_public))

    print(
        "\n// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )
    print("// Secret Sharing Computation")
    print(
        "// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
    )

    # ------------------------------------------------------------------------- Alice
    set_zero_ops()
    public_validation = validate(b_public)
    assert public_validation
    RUNNING_TIME = get_ops()
    print("// Secret sharing corresponding to Alice")
    print(
        "// public key validation cost  :\t%2.3fM + %2.3fS + %2.3fa = %2.3fM,"
        % (
            RUNNING_TIME[0] / (10.0 ** 6),
            RUNNING_TIME[1] / (10.0 ** 6),
            RUNNING_TIME[2] / (10.0 ** 6),
            measure(RUNNING_TIME) / (10.0 ** 6),
        )
    )

    set_zero_ops()
    if (len(temporal_m) == 1) or ((len(temporal_m) == 2) and (0 in temporal_m)):
        ss_a = GAE(
            b_public,
            a_private,
            [L_out[0]],
            [R_out[0]],
            [S_out[0]],
            [temporal_m[-1]],
            m,
        )
    else:
        ss_a = GAE(b_public, a_private, L_out, R_out, S_out, r_out, m)

    RUNNING_TIME = get_ops()

    print(
        "// group action evaluation cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
        % (
            RUNNING_TIME[0] / (10.0 ** 6),
            RUNNING_TIME[1] / (10.0 ** 6),
            RUNNING_TIME[2] / (10.0 ** 6),
            measure(RUNNING_TIME) / (10.0 ** 6),
        )
    )
    print("ss_a := 0x%X;\n" % coeff(ss_a))
    print("expected: 0x1ADB783878BA330BB2A842E7F8B3392329A2CD3B407900E4CF6A8F13B744BFFEFF617BDE2CEBBB9CE97D32BC6FC1BCE2D88381B03B3E13CFF0651EEA82D02937")

    # ------------------------------------------------------------------------- Bob
    set_zero_ops()
    public_validation = validate(a_public)
    assert public_validation
    RUNNING_TIME = get_ops()
    print("// Secret sharing corresponding to Bob")
    print(
        "// public key validation cost  :\t%2.3fM + %2.3fS + %2.3fa = %2.3fM,"
        % (
            RUNNING_TIME[0] / (10.0 ** 6),
            RUNNING_TIME[1] / (10.0 ** 6),
            RUNNING_TIME[2] / (10.0 ** 6),
            measure(RUNNING_TIME) / (10.0 ** 6),
        )
    )

    set_zero_ops()
    if (len(temporal_m) == 1) or ((len(temporal_m) == 2) and (0 in temporal_m)):
        ss_b = GAE(
            a_public,
            b_private,
            [L_out[0]],
            [R_out[0]],
            [S_out[0]],
            [temporal_m[-1]],
            m,
        )
    else:
        ss_b = GAE(a_public, b_private, L_out, R_out, S_out, r_out, m)

    RUNNING_TIME = get_ops()
    print(
        "// group action evaluation cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;"
        % (
            RUNNING_TIME[0] / (10.0 ** 6),
            RUNNING_TIME[1] / (10.0 ** 6),
            RUNNING_TIME[2] / (10.0 ** 6),
            measure(RUNNING_TIME) / (10.0 ** 6),
        )
    )
    print("ss_b := 0x%X;\n" % coeff(ss_b))


    print(
        "\n// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    )
    if coeff(ss_a) == coeff(ss_b):
        print('\x1b[0;30;43m' + '\"Successfully passed!\";' + '\x1b[0m')
    else:
        print(
            '\x1b[0;30;41m'
            + '\"Great Scott!... The sky is falling. NOT PASSED!!!\"'
            + '\x1b[0m'
        )

@click.version_option()
@click.command(short_help="csidh utility for key generation and key derivation")
@click.option(
    "-p",
    "--prime",
    type=click.Choice(['b2', 'b3', 'b5', 'b6', 'p1024', 'p1792', 'p512', 'sv']),
    default="p512",
)
@click.option(
    "-f",
    "--formula",
    type=click.Choice(['tvelu', 'svelu', 'hvelu']),
    default='hvelu',
)
@click.option(
    "-a",
    "--algorithm",
    type=click.Choice(["csidh", "bsidh"]),
    default='csidh',
)
@click.option(
    "-s",
    "--style",
    type=click.Choice(['wd1', 'wd2', 'df']),
    default='df',
)
#   @click.option(
#       "-e",
#       "--exponent",
#       type=click.Choice(['2', '3']),
#       default='2',
#   )
#   @click.option(
#       "-m",
#       "--multievaluation",
#       type=click.Choice(['???', '????']),
#       default='???',
#   )
@click.option(
    "-c",
    "--curvemodel",
    type=click.Choice(['edwards', 'montgomery']),
    default='montgomery',
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Not the kind of verbosity you might expect",
)
def main(prime, formula, algorithm, style, exponent, verbose, multievaluation, curvemodel):
    """
    sidh-csidh-util with reasonable defaults

    eg: sidh -p p512 -f svelu -a csidh -b 1 -s df -v

    """
    if formula == 'hvelu':
        from sidh.csidh._gae_df import gae_df
        gae = Gae_df(prime, style, verbose)
#    csidh_dh(pk, sk)
    raise Exit(0)

if __name__ == '__main__':
    main()
