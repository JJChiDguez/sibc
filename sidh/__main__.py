from logging import getLogger
import click
from base64 import b64encode, b64decode
from click.exceptions import Exit
from sidh.csidh.bench import csidh_bench
from sidh.csidh.bounds import csidh_bounds
from sidh.csidh.header import csidh_header
from sidh.csidh.ijk import csidh_ijk
from sidh.csidh.print_parameters import csidh_parameters
from sidh.csidh.print_strategy import csidh_strategy
from sidh.csidh.sdacs import csidh_sdacs
from sidh.csidh.suitable_bounds import csidh_suitable_bounds
from sidh.csidh.test import csidh_test
from sidh.csidh.main import csidh_main
from sidh.bsidh.test import bsidh_test
from sidh.bsidh.print_strategy import bsidh_strategy
from sidh.bsidh.print_parameters import bsidh_parameters
from sidh.bsidh.main import bsidh_main
from sidh.printstrategy import print_strategy
from sidh.timing import print_timing
from sidh.common import attrdict
from sidh.constants import parameters


@click.version_option()
@click.group(short_help="sidh utility")
@click.option(
    "-p",
    "--prime",
    type=click.Choice(
        ['b2', 'b3', 'b5', 'b6', 'p1024', 'p1792', 'p512', 'sv']
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
    type=click.Choice(['5','10']),
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
        from sidh.csidh import CSIDH

        algo = CSIDH(**algo_args)
    elif algorithm == 'bsidh':
        from sidh.bsidh import BSIDH
        algo_args.pop('style')
        algo_args.pop('exponent')
        algo = BSIDH(**algo_args)
    else:
        click.echo('algorithm not implemented')
        raise Exit(1)
    kwargs['algo'] = algo
    ctx.meta['sidh.kwargs'] = attrdict(kwargs)


@main.command()
@click.pass_context
def genkey(ctx):
    "Generate a secret key"
    algo = ctx.meta['sidh.kwargs']['algo']
    click.echo(b64encode(algo.secret_key()))


@main.command()
@click.argument('secret_key', type=click.File())
@click.pass_context
def pubkey(ctx, secret_key):
    "Generate a public key from a secret key"
    algo = ctx.meta['sidh.kwargs']['algo']
    click.echo(b64encode(algo.public_key(b64decode(secret_key.read()))))


@main.command()
@click.argument('secret_key', type=click.File())
@click.argument('public_key', type=str)
@click.pass_context
def dh(ctx, secret_key, public_key):
    "Generate a shared secret between a secret key and a public key"
    algo = ctx.meta['sidh.kwargs']['algo']
    click.echo(
        b64encode(algo.dh(b64decode(secret_key.read()), b64decode(public_key)))
    )


main.add_command(print_timing)
main.add_command(print_strategy)
main.add_command(csidh_bench)
main.add_command(csidh_bounds)
main.add_command(csidh_header)
main.add_command(csidh_ijk)
main.add_command(csidh_parameters)
main.add_command(csidh_strategy)
main.add_command(csidh_sdacs)
main.add_command(csidh_suitable_bounds)
main.add_command(csidh_test)
main.add_command(csidh_main)
main.add_command(bsidh_test)
main.add_command(bsidh_strategy)
main.add_command(bsidh_parameters)
main.add_command(bsidh_main)

if __name__ == '__main__':
    main()
