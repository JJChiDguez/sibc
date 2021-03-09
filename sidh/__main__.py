from logging import getLogger
import click
from click.exceptions import Exit
from .constants import parameters
from sidh.csidh.bounds import bounds
from sidh.common import attrdict

@click.version_option()
@click.group(short_help="sidh utility")
@click.option(
    "-p",
    "--prime",
    type=click.Choice(['b2', 'b3', 'b5', 'b6', 'p1024', 'p1792', 'p512', 'sv']),
    default="p512",
)
@click.option(
    "-f",
    "--formula",
    type=click.Choice(['tvelu', 'svelu', 'hvelu']),
    default='hvelu',
)
@click.option(
    "-a",
    "--algorithm",
    type=click.Choice(["csidh", "bsidh"]),
    default='csidh',
)
@click.option(
    "-s",
    "--style",
    type=click.Choice(['wd1', 'wd2', 'df']),
    default='df',
)
@click.option(
   "-e",
   "--exponent",
   type=click.Choice(['2', '3']),
   default='2',
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
)
@click.option(
    "-b",
    "--benchmark",
    default=128,
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Not the kind of verbosity you might expect",
)
@click.pass_context
def main(ctx, **kwargs):
    """
    sidh main stub
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
    algo = ctx.meta['sidh.kwargs']['algo']
    click.echo(algo.random_key())

main.add_command(bounds)

if __name__ == '__main__':
    main()
