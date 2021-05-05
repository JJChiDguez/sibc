import click

from sibc.common import attrdict, printl

@click.command()
@click.pass_context
def csidh_sdacs(ctx):
    """ SDACs as C-code headers files """
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
    L = algo.params.L
    n = algo.params.n
    m = algo.params.m
    exponent_of_two = (len(bin(algo.curve.cofactor)[2:]) - 1)
    validation_stop = algo.curve.validation_stop
    SDACS = algo.curve.SDACS

    # ---
    k = 3  # Number of rows (format of the list)
    # ---

    print("#ifndef _SDACS_%s_H_" % setting.prime)
    print("#define _SDACS_%s_H_" % setting.prime)
    print("")
    assert n == len(L)

    if (2**exponent_of_two ) == algo.curve.cofactor:
        print(
            "#define cofactor %d\t// Exponent of 2 in the factorization of (p + 1) "
            % exponent_of_two
        )
        print("")
    print("#ifdef _MONT_C_CODE_")

    print(
        "#define UPPER_BOUND %d\t// Bits of 4 * sqrt( [p + 1] / [2^e] )"
        % validation_stop
    )
    print(
        "\n// Recall, the list of Small Odd Primes (SOPs) is stored such that l_0 < l_1 < ... < l_{n-1}"
    )

    assert n == len(SDACS)
    print(
        "\n// The list of Shortest Differential Addition Chains (SDACs) corresponding with each l_i"
    )
    printl(
        "static int LENGTHS[]", [len(sdac_i) for sdac_i in SDACS], n // k + 1
    )
    print("")
    for i in range(0, n, 1):
        if len(SDACS[i]) == 0:
            print("static char SDAC%d[1];" % i)
        else:
            print(
                'static char SDAC'
                + str(i)
                + '[] = \"'
                + ''.join([str(sdac_ij) for sdac_ij in SDACS[i]])
                + '\";'
            )

    SDACS_string = "static char *SDACs[N] = {\n\t"

    for i in range(0, n - 1, 1):
        if (i + 1) % (n // k + 1) == 0:
            SDACS_string = SDACS_string + "SDAC%d, \n\t" % (i)
        else:
            SDACS_string = SDACS_string + "SDAC%d, " % (i)

    SDACS_string = SDACS_string + "SDAC%d\n\t};" % (n - 1)
    print("")
    print(SDACS_string)
    print("#endif")
    print(
        "\n#endif /* required framework for the SDACs, which is used in CSIDH-%s */"
        % setting.prime[1:]
    )
    return attrdict(name='csidh-sdacs', **locals())
