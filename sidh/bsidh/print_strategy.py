import click

from sidh.common import attrdict
from sidh.constants import strategy_data
from sidh.bsidh.strategy import Gae

from pkg_resources import resource_filename


@click.command()
@click.pass_context
def bsidh_strategy(ctx):
    algo = ctx.meta['sidh.kwargs']['algo']
    setting = ctx.meta['sidh.kwargs']
    tuned_name = ('-classical','-suitable')[setting.tuned]
    SIDp = algo.curve.SIDp
    SIDm = algo.curve.SIDm

    f_name = 'data/strategies/'+setting.algorithm+'-'+setting.prime+'-'+algo.formula.name+tuned_name
    try:
        f = open(resource_filename('sidh', f_name))
        print("// Strategies to be read from a file")
        # Corresponding to the list of Small Isogeny Degree, Lp := [l_0, ...,
        # l_{n-1}] [We need to include case l=2 and l=4]
        tmp = f.readline()
        tmp = [int(b) for b in tmp.split()]
        Sp = list(tmp)
        # Corresponding to the list of Small Isogeny Degree, Lm := [l_0, ..., l_{n-1}]
        tmp = f.readline()
        tmp = [int(b) for b in tmp.split()]
        Sm = list(tmp)
        f.close()
    except IOError:
        print("// Strategies to be computed")
        # List of Small Isogeny Degree, Lp := [l_0, ..., l_{n-1}] [We need to
        # include case l=2 and l=4]
        Sp, Cp = algo.gae.dynamic_programming_algorithm(SIDp[::-1], len(SIDp))
        # List of Small Isogeny Degree, Lm := [l_0, ..., l_{n-1}]
        Sm, Cm = algo.gae.dynamic_programming_algorithm(SIDm[::-1], len(SIDm))

        f_name = 'data/strategies/x'+setting.algorithm+'-'+setting.prime+'-'+algo.formula.name+tuned_name
        f = open(f_name, 'w',)

        f.writelines(' '.join([str(tmp) for tmp in Sp]) + '\n')
        f.writelines(' '.join([str(tmp) for tmp in Sm]) + '\n')

        f.close()

    # List of strategies
    S_out = [list(Sp), list(Sm)]
    filename = (
        setting.algorithm + '-' + setting.prime + '-' + setting.formula + tuned_name
    )
