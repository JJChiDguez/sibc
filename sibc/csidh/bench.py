from pkg_resources import resource_filename
import click

from progress.bar import Bar
import statistics

from sibc.common import attrdict, geometric_serie, rounds

@click.command()
@click.pass_context
def csidh_bench(ctx):
    """ Average GF(p)-operation cost of a GAE """
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
    random_exponents = algo.gae.random_exponents
    print_exponents = algo.gae.print_exponents

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
