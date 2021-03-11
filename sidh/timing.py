import click
import timeit
from sidh.csidh import CSIDH

@click.command()
@click.pass_context
def print_timing(ctx):
    setting = ctx.meta['sidh.kwargs']
    setting.pop('algo')
    setting.pop('algorithm')
    rounds = setting.pop('benchmark')
    click.echo("Running ({} rounds):".format(rounds))
    from sidh.csidh import CSIDH
    c = CSIDH(**setting)
    s = """sk = c.random_key()"""
    click.echo(s)
    click.echo("with garbage collection:")
    print(timeit.timeit(s, number=rounds, globals=locals(), setup='import gc;gc.enable()'))
    click.echo("without garbage collection:")
    print(timeit.timeit(s, number=rounds, globals=locals()))
    s += """\npk = c.pubkey(sk)"""
    click.echo("with garbage collection:")
    click.echo(s)
    print(timeit.timeit(s, number=rounds, globals=locals(), setup='import gc;gc.enable()'))
    click.echo("without garbage collection:")
    print(timeit.timeit(s, number=rounds, globals=locals()))
