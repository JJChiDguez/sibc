from logging import getLogger
import click
from base64 import b64encode, b64decode
from click.exceptions import Exit
# csidh
from sibc.csidh.bench import csidh_bench
from sibc.csidh.bounds import csidh_bounds
from sibc.csidh.header import csidh_header
from sibc.csidh.ijk import csidh_ijk
from sibc.csidh.precompute_parameters import csidh_precompute_parameters
from sibc.csidh.precompute_strategy import csidh_precompute_strategy
from sibc.csidh.sdacs import csidh_sdacs
from sibc.csidh.test import csidh_test
from sibc.csidh.main import csidh_main
# bsidh
from sibc.bsidh.test import bsidh_test
from sibc.bsidh.precompute_strategy import bsidh_precompute_strategy
from sibc.bsidh.precompute_parameters import bsidh_precompute_parameters
from sibc.bsidh.main import bsidh_main

from sibc.plot_strategy import plot_strategy
from sibc.timing import print_timing
from sibc.common import attrdict
from sibc.constants import parameters


@click.version_option()
@click.group(short_help="sibc utility")
@click.option(
    "-p",
    "--prime",
    type=click.Choice(
        ['b2', 'b3', 'b5', 'b6', 'p1024', 'p1792', 'p512', 's1']
    ),
    default="p512",
    show_default=True,
)
@click.option(
    "-f",
    "--formula",
    type=click.Choice(['tvelu', 'svelu', 'hvelu']),
    default='hvelu',
    show_default=True,
)
@click.option(
    "-a",
    "--algorithm",
    type=click.Choice(["csidh", "bsidh"]),
    default='csidh',
    show_default=True,
)
@click.option(
    "-s",
    "--style",
    type=click.Choice(['wd1', 'wd2', 'df']),
    default='df',
    show_default=True,
)
@click.option(
    "-e",
    "--exponent",
    type=click.Choice([f'{i}' for i in range(1, 11, 1)]),
    default='10',
    show_default=True,
)
@click.option(
    "-m", "--multievaluation", is_flag=True, help="", show_default=True,
)
@click.option(
    "-c",
    "--curvemodel",
    type=click.Choice(['edwards', 'montgomery']),
    default='montgomery',
    show_default=True,
)
@click.option(
    "-b", "--benchmark", default=128, show_default=True,
)
@click.option(
    "-t", "--tuned", is_flag=True, help="", show_default=True,
)
@click.option(
    "-u", "--uninitialized", is_flag=True, help="", show_default=True,
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Not the kind of verbosity you might expect",
    show_default=True,
)
@click.pass_context
def main(ctx, **kwargs):
    """

    \b
      ,-~~-.___.          
     / |  '     \\        
    (  )         0        
     \_/-, ,----'         
        ====           // 
       /  \-'~;    /~~~(O)
      /  __/~|   /       |
    =(  _____| (_________|

    """
    algo_args = kwargs.copy()
    algorithm = algo_args.pop('algorithm')
    algo_args.pop('benchmark')
    if algorithm == 'csidh':
        from sibc.csidh import CSIDH

        algo = CSIDH(**algo_args)
    elif algorithm == 'bsidh':
        from sibc.bsidh import BSIDH
        algo_args.pop('style')
        algo_args.pop('exponent')
        algo = BSIDH(**algo_args)
    else:
        click.echo('algorithm not implemented')
        raise Exit(1)
    kwargs['algo'] = algo
    ctx.meta['sibc.kwargs'] = attrdict(kwargs)


@main.command()
@click.pass_context
def csidh_genkey(ctx):
    "Generate random CSIDH secret key"
    algo = ctx.meta['sibc.kwargs']['algo']
    click.echo(b64encode(algo.secret_key()))


@main.command()
@click.argument('secret_key', type=click.File())
@click.pass_context
def csidh_pubkey(ctx, secret_key):
    "Derive CSIDH public key from CSIDH secret key"
    algo = ctx.meta['sibc.kwargs']['algo']
    click.echo(b64encode(algo.public_key(b64decode(secret_key.read()))))


@main.command()
@click.argument('secret_key', type=click.File())
@click.argument('public_key', type=str)
@click.pass_context
def csidh_dh(ctx, secret_key, public_key):
    "Derive shared secret key from CSIDH sk, CSIDH pk"
    algo = ctx.meta['sibc.kwargs']['algo']
    click.echo(
        b64encode(algo.dh(b64decode(secret_key.read()), b64decode(public_key)))
    )


main.add_command(print_timing)
main.add_command(plot_strategy)
# csidh
main.add_command(csidh_bench)
main.add_command(csidh_bounds)
main.add_command(csidh_header)
main.add_command(csidh_ijk)
main.add_command(csidh_precompute_parameters)
main.add_command(csidh_precompute_strategy)
main.add_command(csidh_sdacs)
main.add_command(csidh_test)
main.add_command(csidh_main)
# bsidh
main.add_command(bsidh_test)
main.add_command(bsidh_precompute_strategy)
main.add_command(bsidh_precompute_parameters)
main.add_command(bsidh_main)

if __name__ == '__main__':
    main()
