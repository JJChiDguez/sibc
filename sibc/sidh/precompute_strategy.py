import click
from pkg_resources import resource_filename

from sibc.common import attrdict


@click.command()
@click.pass_context
def sidh_precompute_strategy(ctx):
    """ Precomputation of optimal strategies """
    algo = ctx.meta['sibc.kwargs']['algo']
    setting = ctx.meta['sibc.kwargs']
    e2 = algo.strategy.two
    e3 = algo.strategy.three

    assert setting.uninitialized, 'option -u (--uninitialized) is required!'

    file_path = (
        "data/strategies/"
        + algo.curve.model
        + '/sidh/'
        + 'sidh'
        + '-'
        + setting.prime
    )
    file_path = resource_filename('sibc', file_path)
    print("// Strategies to be computed")
    # List of Small Isogeny Degree, Lp := [l_0, ..., l_{n-1}] [We need to
    # include case l=2 and l=4]
    S2, C2 = algo.strategy.dynamic_programming_algorithm(2, e2)
    # List of Small Isogeny Degree, Lm := [l_0, ..., l_{n-1}]
    S3, C3 = algo.strategy.dynamic_programming_algorithm(3, e3)
    f = open(file_path, 'w')
    f.writelines(' '.join([str(tmp) for tmp in S2]) + '\n')
    f.writelines(' '.join([str(tmp) for tmp in S3]) + '\n')
    f.close()

    return attrdict(name='sidh-precompute-strategy', **locals())