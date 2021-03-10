import click
from progress.bar import Bar
from sidh.constants import strategy_data
from sidh.csidh import CSIDH
import statistics

@click.command()
@click.pass_context
def csidh_bench(ctx):
    setting = ctx.meta['sidh.kwargs']
    self = setting['algo']
    n = self.params.n
    m = self.params.m
    L = self.params.L

    if len(set(m)) > 1:
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
        m_prime = [self.gae.geometric_serie(m[k], L[k]) for k in range(n)]
        r_out, L_out, R_out = self.gae.rounds(m_prime[::-1], n)
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
            + setting.formula
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

        print("// Strategies to be computed")
        C_out, L_out, R_out, S_out, r_out = self.gae.strategy_block_cost(L[::-1], m[::-1])
        f = open(
            strategy_data
            + setting.algorithm
            + '-'
            + setting.prime
            + '-'
            + setting.style
            + '-'
            + setting.formula
            + '-'
            + LABEL_m
            + verb,
            'w',
        )
        for i in range(0, len(r_out)):

            f.writelines(' '.join([str(tmp) for tmp in S_out[i]]) + '\n')

        f.close()

    print(
        "// All the experiments are assuming S = %2.3fM and a = %2.3fM. The measures are given in millions of field operations.\n"
        % (self.curve.SQR, self.curve.ADD)
    )

    ''' -------------------------------------------------------------------------------------
        Framework
        ------------------------------------------------------------------------------------- '''

    # print("p := 0x%X;" % p)
    # print("fp := GF(p);")
    # print("P<x> := PolynomialRing(fp);");
    # print("fp2<i> := ext<fp | x^2 + 1>;")
    # print("P<x> := PolynomialRing(fp2);");

    A = self.curve.A
    assert A == [2, 4]
    # print("public_coeff := 0x%X;\n" % coeff(A))

    ''' -------------------------------------------------------------------------------------
        Main
        ------------------------------------------------------------------------------------- '''

    e = self.gae.random_key(m)
    B = list(A)

    tmp = list(set(m))
    SAMPLE = [[0.0, 0.0, 0.0]] * setting.benchmark
    SAMPLE_VALIDATE = [[0.0, 0.0, 0.0]] * setting.benchmark
    bar = Bar('// Running experiments', max=setting.benchmark)

    for main_i in range(setting.benchmark):

        e = self.gae.random_key(m)

        if (len(tmp) == 1) or ((len(tmp) == 2) and (0 in tmp)):

            self.fp.set_zero_ops()
            B = self.gae.GAE(B, e, [L_out[0]], [R_out[0]], [S_out[0]], [tmp[-1]], m)
            SAMPLE[main_i] = self.fp.get_ops()

        else:

            self.fp.set_zero_ops()
            B = self.gae.GAE(B, e, L_out, R_out, S_out, r_out, m)
            SAMPLE[main_i] = self.fp.get_ops()

        self.fp.set_zero_ops()
        V = self.curve.validate(B)
        assert V
        SAMPLE_VALIDATE[main_i] = self.fp.get_ops()

        # print("Random(EllipticCurve(x^3 + 0x%X * x^2 + x)) * (p+1);" % coeff(B))

        bar.next()
    bar.finish()

    AVERAGE = [
        statistics.mean([ops[0] for ops in SAMPLE]),
        statistics.mean([ops[1] for ops in SAMPLE]),
        statistics.mean([ops[2] for ops in SAMPLE]),
    ]
    print(
        "// Average number of field operations (GAE):\t\t%2.3fM + %2.3fS + %2.3fa := %2.3fM"
        % (
            AVERAGE[0] / (10.0 ** 6),
            AVERAGE[1] / (10.0 ** 6),
            AVERAGE[2] / (10.0 ** 6),
            self.curve.measure(AVERAGE) / (10.0 ** 6),
        )
    )

    AVERAGE = [
        statistics.mean([ops[0] for ops in SAMPLE_VALIDATE]),
        statistics.mean([ops[1] for ops in SAMPLE_VALIDATE]),
        statistics.mean([ops[2] for ops in SAMPLE_VALIDATE]),
    ]
    print(
        "// Average number of field operations (validate):\t%2.3fM + %2.3fS + %2.3fa := %2.3fM\n"
        % (
            AVERAGE[0] / (10.0 ** 6),
            AVERAGE[1] / (10.0 ** 6),
            AVERAGE[2] / (10.0 ** 6),
            self.curve.measure(AVERAGE) / (10.0 ** 6),
        )
    )

