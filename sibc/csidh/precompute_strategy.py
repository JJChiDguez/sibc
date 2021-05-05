from pkg_resources import resource_filename
import click

from sibc.common import attrdict, geometric_serie, rounds

@click.command()
@click.pass_context
def csidh_precompute_strategy(ctx):
    """ Precomputation of optimal strategies """
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

    assert setting.uninitialized, 'option -u (--uninitialized) is required!'

    if algo.formula.name != 'tvelu':
        set_parameters_velu = algo.formula.set_parameters_velu

    temporal_m = list(set(m))
    if len(temporal_m) > 1:
        # Maximum number of degree-(l_i) isogeny constructions is m_i (different for each l_i)
        bounds = '-diffbounds'
    else:
        # Maximum number of degree-(l_i) isogeny constructions is m (the same for each l_i)
        bounds = '-samebounds'

    tuned = {True:'-tuned', False:''}[setting.tuned]
    multievaluation = {True:'scaled', False:'unscaled'}[setting.multievaluation]
    # List of Small Odd Primes, L := [l_0, ..., l_{n-1}]
    m_prime = [geometric_serie(m[k], L[k]) for k in range(n)]
    r_out, L_out, R_out = rounds(m_prime[::-1], n)
    for j in range(0, len(r_out), 1):

        R_out[j] = list([L[::-1][k] for k in R_out[j]])
        L_out[j] = list([L[::-1][k] for k in L_out[j]])

    file_path = (
        "data/strategies/"
        + algo.curve.model
        + '/csidh/'
        + 'csidh'
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
        + multievaluation
        + tuned
    )
    file_path = resource_filename('sibc', file_path)
    print("// Strategies to be computed")
    C_out, L_out, R_out, S_out, r_out = strategy_block_cost(
        L[::-1], m[::-1]
    )
    f = open(file_path, 'w')
    for i in range(0, len(r_out)):

        f.writelines(' '.join([str(tmp) for tmp in S_out[i]]) + '\n')

    f.close()

    return attrdict(name='csidh-precompute-strategy', **locals())