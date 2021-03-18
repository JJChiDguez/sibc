import click

from sidh.common import attrdict
from sidh.constants import strategy_data


@click.command()
@click.pass_context
def csidh_strategy(ctx):
    algo = ctx.meta['sidh.kwargs']['algo']
    setting = ctx.meta['sidh.kwargs']
    L = algo.params.L
    n = algo.params.n
    m = algo.params.m
    delta = algo.params.delta
    geometric_serie = algo.gae.geometric_serie
    strategy_block_cost = algo.gae.strategy_block_cost
    rounds = algo.gae.rounds
    sigma, kappa = algo.params.sigma, algo.params.kappa

    # ==========================================================================

    if len(set(m)) > 1:
        # Maximum number of degree-(l_i) isogeny constructions is m_i (different for each l_i)
        LABEL_m = 'different_bounds'
    else:
        # Maximum number of degree-(l_i) isogeny constructions is m (the same for each l_i)
        LABEL_m = 'with_same_bounds'

    if setting.tuned:
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
        C_out, L_out, R_out, S_out, r_out = strategy_block_cost(
            L[::-1], m[::-1]
        )
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

    filename = (
        setting.algorithm
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
    return attrdict(name='csidh-strategy', **locals())