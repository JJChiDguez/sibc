import click
import numpy

from sibc.common import attrdict, geometric_serie, rounds, printl
from pkg_resources import resource_filename

@click.command()
@click.pass_context
def csidh_header(ctx):
    """ Optimal strategies as C-code headers files """
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
    L = algo.params.L
    n = algo.params.n
    m = algo.params.m

    strategy_block_cost = algo.gae.strategy_block_cost
    basis = numpy.eye(n, dtype=int)
    measure = algo.curve.measure

    # ==========================================================================

    if algo.formula.name != 'tvelu':
        set_parameters_velu = algo.formula.set_parameters_velu

    temporal_m = list(set(m))
    if len(temporal_m) > 1:
        # Maximum number of degree-(l_i) isogeny constructions is m_i (different for each l_i)
        bounds = '-diffbounds'
    else:
        # Maximum number of degree-(l_i) isogeny constructions is m (the same for each l_i)
        bounds = '-samebounds'

    # List of Small Odd Primes, L := [l_0, ..., l_{n-1}]
    m_prime = [geometric_serie(m[k], L[k]) for k in range(n)]
    r_out, L_out, R_out = rounds(m_prime[::-1], n)
    for j in range(0, len(r_out), 1):

        R_out[j] = list([L[::-1][k] for k in R_out[j]])
        L_out[j] = list([L[::-1][k] for k in L_out[j]])

    file_path = (
        "data/strategies/"
        + algo.curve.model
        + '/'
        + 'csidh/csidh'
        + '-'
        + setting.prime
        + '-'
        + setting.style
        + '-e'
        + setting.exponent
        + bounds
        + '-'
        + setting.formula
        + '-'
        + algo.formula.multievaluation_name
        + algo.formula.tuned_name
    )
    file_path = resource_filename('sibc', file_path)
    f = open(file_path)
    S_out = []
    for i in range(0, len(r_out), 1):

        tmp = f.readline()
        tmp = [int(b) for b in tmp.split()]
        S_out.append(tmp)

    f.close()

    if (len(set(m)) == 1) or ((len(set(m)) == 2) and (0 in set(m))):
        L_out = list([L_out[0]])
        R_out = list([R_out[0]])
        S_out = list([S_out[0]])
        r_out = list([r_out[0]])
    # ---
    k = 3  # Number of rows (format of the list)
    # ---

    print("#ifndef _STRATEGIES_H_")
    print("#define _STRATEGIES_H_\n")
    print("#include <inttypes.h>\n\n")
    print(
        "// This script assumes the C-code implementation has the list of Small Odd Primes (SOPs) stored such that l_0 < l_1 < ... < l_{n-1}"
    )
    print("// Recall, the strategies process from small SOPs to large SOPs.\n")

    for i in range(len(r_out)):
        print(
            "// -----------------------------------------------------------------------------------------------------------------------------------"
        )
        print("// Strategy number %d\n" % (i))
        L_string = "static uint32_t L%d[] " % (i)
        R_string = "static uint32_t W%d[] " % (i)
        S_string = "static uint32_t S%d[] " % (i)

        printl(
            L_string,
            [L.index(l_i) for l_i in L_out[i]],
            len(L_out[i]) // k + 1,
        )
        if R_out[i] != []:
            printl(
                R_string,
                [L.index(r_i) for r_i in R_out[i]],
                len(R_out[i]) // k + 1,
            )
        else:
            print("static uint32_t W%d[1];" % i)
        if S_out[i] != []:
            printl(S_string, S_out[i], len(S_out[i]) // k + 1)
        else:
            print("static uint32_t S%d[1];" % i)

    print("\n")
    print(
        "// -----------------------------------------------------------------------------------------------------------------------------------"
    )
    print(
        "// -----------------------------------------------------------------------------------------------------------------------------------"
    )
    print("#define NUMBER_OF_DIFFERENT_STRATEGIES  %d" % len(L_out))
    print("")
    L_string = (
        "static uint32_t *L_STRATEGY[NUMBER_OF_DIFFERENT_STRATEGIES] = {\n\t"
    )
    R_string = (
        "static uint32_t *W_STRATEGY[NUMBER_OF_DIFFERENT_STRATEGIES] = {\n\t"
    )
    S_string = "static uint32_t *S[NUMBER_OF_DIFFERENT_STRATEGIES] = {\n\t"

    tmp_sizes = "static uint32_t NUMBER_OF_PRIMES[] = {\n\t"
    tmp_round = "static uint8_t ROUNDS[] = {\n\t"

    for i in range(len(L_out) - 1):

        L_string = L_string + "L%d, " % (i)
        R_string = R_string + "W%d, " % (i)
        S_string = S_string + "S%d, " % (i)

        tmp_sizes = tmp_sizes + "%3d," % (len(L_out[i]))
        tmp_round = tmp_round + "%3d," % (r_out[i])

    L_string = L_string + "L%d\n\t};" % (len(L_out) - 1)
    R_string = R_string + "W%d\n\t};" % (len(L_out) - 1)
    S_string = S_string + "S%d\n\t};" % (len(L_out) - 1)

    tmp_sizes = tmp_sizes + "%3d\n\t};" % (len(L_out[len(L_out) - 1]))
    tmp_round = tmp_round + "%3d\n\t};" % (r_out[len(L_out) - 1])

    print(
        "// L_STRATEGY[i] determines the small odd primes l_i per each strategy"
    )
    print(L_string)
    print("\n// W_STRATEGY[i] determines L - L_STRATEGY[i]")
    print(R_string)
    print(
        "\n// S_STRATEGY[i] determines the optimal strategy for L_STRATEGY[i]"
    )
    print(S_string)

    print("\n// Number of primes for each strategy")
    print(tmp_sizes)
    print("")
    print("// Number of rounds per each different strategy")
    print(tmp_round)

    print("")
    print("// Maximum number of degree-(l_i) isogeny constructions")
    printl("static uint8_t M[]", m, n // k + 1)

    STYLE_NAME = {
        'wd2': 'OAYT-style',
        'wd1': 'MCR-style',
        'df': 'dummy-free-style',
    }[setting.style]
    print(
        "\n#endif /* required framework for the strategies to be used in CSIDH-%s using %s */"
        % (setting.prime[1:], STYLE_NAME)
    )
    # Need empty line at EOF for -Wnewline-eof
    print("") 
    return attrdict(name='header', **locals())
