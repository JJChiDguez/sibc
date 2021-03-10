from logging import getLogger
import click
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
from sidh.printstrategy import print_strategy
from sidh.common import attrdict
from sidh.constants import parameters

@click.version_option()
@click.group(short_help="sidh utility")
@click.option(
    "-p",
    "--prime",
    type=click.Choice(['b2', 'b3', 'b5', 'b6', 'p1024', 'p1792', 'p512', 'sv']),
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
   type=click.Choice(['2']),
   default='2',
    show_default=True,
)
#   @click.option(
#       "-m",
#       "--multievaluation",
#       type=click.Choice(['???', '????']),
#       default='???',
#   )
@click.option(
    "-c",
    "--curvemodel",
    type=click.Choice(['edwards', 'montgomery']),
    default='montgomery',
    show_default=True,
)
@click.option(
    "-b",
    "--benchmark",
    default=128,
    show_default=True,
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
        click.echo('BSIDH not yet implemented; try again later')
        raise Exit(1)
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
    click.echo(algo.random_key())

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

if __name__ == '__main__':
    main()
