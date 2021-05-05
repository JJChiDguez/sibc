import click
from pkg_resources import resource_filename

from sibc.common import attrdict


@click.command()
@click.pass_context
def bsidh_precompute_strategy(ctx):
    """ Precomputation of optimal strategies """
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
    SIDp = algo.strategy.SIDp
    SIDm = algo.strategy.SIDm

    assert setting.uninitialized, 'option -u (--uninitialized) is required!'

    tuned = {True:'-tuned', False:''}[setting.tuned]
    multievaluation = {True:'scaled', False:'unscaled'}[setting.multievaluation]
    file_path = (
        "data/strategies/"
        + algo.curve.model
        + '/bsidh/'
        + 'bsidh'
        + '-'
        + setting.prime
        + '-'
        + setting.formula
        + '-'
        + multievaluation
        + tuned
    )
    file_path = resource_filename('sibc', file_path)
    print("// Strategies to be computed")
    # List of Small Isogeny Degree, Lp := [l_0, ..., l_{n-1}] [We need to
    # include case l=2 and l=4]
    Sp, Cp = algo.strategy.dynamic_programming_algorithm(SIDp[::-1], len(SIDp))
    # List of Small Isogeny Degree, Lm := [l_0, ..., l_{n-1}]
    Sm, Cm = algo.strategy.dynamic_programming_algorithm(SIDm[::-1], len(SIDm))
    f = open(file_path, 'w')
    f.writelines(' '.join([str(tmp) for tmp in Sp]) + '\n')
    f.writelines(' '.join([str(tmp) for tmp in Sm]) + '\n')
    f.close()

    return attrdict(name='bsidh-precompute-strategy', **locals())