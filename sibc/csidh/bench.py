from pkg_resources import resource_filename
import click

from progress.bar import Bar
import statistics

from sibc.common import attrdict, geometric_serie, rounds
from sibc.constants import strategy_data

@click.command()
@click.pass_context
def csidh_bench(ctx):
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
    coeff = algo.curve.coeff
    p = algo.params.p
    L = algo.params.L
    m = algo.params.m
    n = algo.params.n
    SQR, ADD = algo.curve.SQR, algo.curve.ADD
    init_runtime = algo.field.init_runtime
    validate = algo.curve.issupersingular
    measure = algo.curve.measure
    GAE_at_0 = algo.gae.GAE_at_0
    GAE_at_A = algo.gae.GAE_at_A
    strategy_block_cost = algo.gae.strategy_block_cost
    random_exponents = algo.gae.random_exponents
    print_exponents = algo.gae.print_exponents

    if algo.formula.name != 'tvelu':
        set_parameters_velu = algo.formula.set_parameters_velu

    temporal_m = list(set(m))
    if len(temporal_m) > 1:
        # Maximum number of degree-(l_i) isogeny constructions is m_i (different for each l_i)
        LABEL_m = 'different_bounds'
    else:
        # Maximum number of degree-(l_i) isogeny constructions is m (the same for each l_i)
        LABEL_m = 'with_same_bounds'

    if setting.tuned:
        verb = '-suitable'
    else:
        verb = '-classical'

    # List of Small Odd Primes, L := [l_0, ..., l_{n-1}]
    m_prime = [geometric_serie(m[k], L[k]) for k in range(n)]
    r_out, L_out, R_out = rounds(m_prime[::-1], n)
    for j in range(0, len(r_out), 1):

        R_out[j] = list([L[::-1][k] for k in R_out[j]])
        L_out[j] = list([L[::-1][k] for k in L_out[j]])

    file_path = (
        "data/strategies/"
        + setting.algorithm
        + '-'
        + setting.prime
        + '-'
        + setting.style
        + '-'
        + setting.formula
        + '-'
        + LABEL_m
        + verb
    )
    file_path = resource_filename('sibc', file_path)
    try:
        f = open(file_path)
        print("// Strategies to be read from a file")
        S_out = []
        for i in range(0, len(r_out), 1):

            tmp = f.readline()
            tmp = [int(b) for b in tmp.split()]
            S_out.append(tmp)

        f.close()

    except IOError:

        print("// Strategies to be computed")
        C_out, L_out, R_out, S_out, r_out = strategy_block_cost(
            L[::-1], m[::-1]
        )
        f = open(file_path, 'w')
        for i in range(0, len(r_out)):

            f.writelines(' '.join([str(tmp) for tmp in S_out[i]]) + '\n')

        f.close()

    print(
        "// The running time is assuming S = %1.2f x M and a = %1.2f x M, and giving in millions of field operations.\n"
        % (SQR, ADD)
    )

    ''' -------------------------------------------------------------------------------------
        Main
        ------------------------------------------------------------------------------------- '''

    e = random_exponents(m)
    A = [algo.curve.field(2), algo.curve.field(4)]

    tmp = list(set(m))
    SAMPLE = [[0.0, 0.0, 0.0]] * setting.benchmark
    SAMPLE_VALIDATE = [[0.0, 0.0, 0.0]] * setting.benchmark
    bar = Bar('// Running experiments', max=setting.benchmark)

    for main_i in range(setting.benchmark):

        e = random_exponents(m)
        init_runtime()
        B = GAE_at_A(e, A)
        SAMPLE[main_i] = [algo.field.fpmul, algo.field.fpsqr, algo.field.fpadd]

        init_runtime()
        V = validate(B)
        assert V
        SAMPLE_VALIDATE[main_i] = [algo.field.fpmul, algo.field.fpsqr, algo.field.fpadd]

        bar.next()
    bar.finish()

    AVERAGE = [
        statistics.mean([ops[0] for ops in SAMPLE]),
        statistics.mean([ops[1] for ops in SAMPLE]),
        statistics.mean([ops[2] for ops in SAMPLE]),
    ]
    print(
        "// Average running time (GAE + key validation):\t%2.3fM + %2.3fS + %2.3fa := %2.3fM"
        % (
            AVERAGE[0] / (10.0 ** 6),
            AVERAGE[1] / (10.0 ** 6),
            AVERAGE[2] / (10.0 ** 6),
            measure(AVERAGE) / (10.0 ** 6),
        )
    )

    AVERAGE = [
        statistics.mean([ops[0] for ops in SAMPLE_VALIDATE]),
        statistics.mean([ops[1] for ops in SAMPLE_VALIDATE]),
        statistics.mean([ops[2] for ops in SAMPLE_VALIDATE]),
    ]
    print(
        "// Average running time (validate):\t\t%2.3fM + %2.3fS + %2.3fa := %2.3fM\n"
        % (
            AVERAGE[0] / (10.0 ** 6),
            AVERAGE[1] / (10.0 ** 6),
            AVERAGE[2] / (10.0 ** 6),
            measure(AVERAGE) / (10.0 ** 6),
        )
    )
