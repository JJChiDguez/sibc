import click
import timeit
from sibc.csidh import CSIDH

@click.command()
@click.pass_context
def print_timing(ctx):
    setting = ctx.meta['sibc.kwargs']
    setting.pop('algo')
    setting.pop('algorithm')
    rounds = setting.pop('benchmark')
    try:
        from cpuinfo import get_cpu_info
        info = get_cpu_info()
        click.echo("python: {}".format(info['python_version']))
        click.echo("arch: {}".format(info['arch']))
        click.echo("speed: {}".format(str(info['hz_actual'])))
        click.echo("cpu: {}".format(info['brand']))
    except:
        pass
    click.echo("Running ({} rounds):".format(rounds))
    from sibc.csidh import CSIDH
    c = CSIDH(**setting)
    s = """sk = c.secret_key()"""
    click.echo(s)
    click.echo("with garbage collection:")
    print(timeit.timeit(s, number=rounds, globals=locals(), setup='import gc;gc.enable()'))
    click.echo("without garbage collection:")
    print(timeit.timeit(s, number=rounds, globals=locals()))
    s += """\npk = c.public_key(sk)"""
    click.echo("with garbage collection:")
    click.echo(s)
    print(timeit.timeit(s, number=rounds, globals=locals(), setup='import gc;gc.enable()'))
    click.echo("without garbage collection:")
    print(timeit.timeit(s, number=rounds, globals=locals()))
