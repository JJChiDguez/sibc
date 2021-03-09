from logging import getLogger
import click
from click.exceptions import Exit
from .constants import parameters

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
#   @click.option(
#       "-e",
#       "--exponent",
#       type=click.Choice(['2', '3']),
#       default='2',
#   )
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
    ctx.meta['sidh.kwargs'] = kwargs

@main.command()
@click.pass_context
def genkey(ctx):
    kwargs = dict(ctx.meta['sidh.kwargs'])
    algorithm = kwargs.pop('algorithm')
    _ = kwargs.pop('multievaluation')
    if algorithm == 'bsidh':
        click.echo('BSIDH not yet implemented; try again later')
        raise Exit(1)
    elif algorithm == 'csidh':
        from sidh.csidh import CSIDH as algo
    else:
        click.echo('algorithm not implemented')
        raise Exit(1)
    click.echo(algo(**kwargs).random_key())


if __name__ == '__main__':
    main()
