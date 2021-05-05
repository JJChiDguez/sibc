import math
import click
from click.exceptions import Exit
from pkg_resources import resource_string, resource_filename

def write_list_of_lists_of_ints_to_file(path, data):
    with open(path, 'w') as fh:
        for line in data:
            fh.writelines(' '.join(str(v) for v in line))
            fh.writelines('\n')

@click.command()
@click.option(
    "-p",
    "--prime",
    type=str,
    default=None,
    show_default=True,
    help='Prime number description given as either csidh or bsidh format (see README.md file for more details)',
)
@click.option(
    "-a",
    "--algorithm",
    type=click.Choice(["csidh", "bsidh"]),
    default='csidh',
    show_default=True,
)
def main(prime, algorithm):
    assert prime != None, 'argument --prime is required'
    """ Computing and storing sdacs """
    f = open(resource_filename('sibc', 'data/sop/' + algorithm + '/' + prime))
    if algorithm == 'csidh':
        # CSIDH only requires the factorization of p + 1
        L = f.readline()
        # The first value in L corresponds with the cofactor h of (p+1), which is not required here
        L = [int(l) for l in L.split()][1:]
        n = len(L)

    elif algorithm == 'bsidh':
        # B-SIDH only requires the factorization of p + 1 and p - 1
        # The prime to be used
        p = f.readline()
        p = int(p, 16)

        # List corresponding (p + 1)
        Lp = f.readline()
        Lp = [int(lp) for lp in Lp.split()]
        # exponent_of_twop = Lp[0]
        # Lp = Lp[1:]
        Ep = f.readline()
        Ep = [int(ep) for ep in Ep.split()]
        assert len(Ep) == len(Lp)
        np = len(Lp)

        # List corresponding (p - 1)
        Lm = f.readline()
        Lm = [int(lm) for lm in Lm.split()]
        Em = f.readline()
        Em = [int(em) for em in Em.split()]
        assert len(Em) == len(Lm)
        nm = len(Lm)

        L = list(Lp + Lm)
        n = len(L)

    else:
        click.echo("only csidh and bsidh are currently implemented")
        raise Exit(1)

    f.close()

    def dacs(l, r0, r1, r2, chain):
        """
        dacs()
        inputs: a small odd prime number l, three integer numbers, and a list
        output: all the differential additions chains corresponding with the input l

        NOTE: this is a recursive approach
        """
        if r2 == l:

            return [(chain, r2)]
        elif r2 < l and len(chain) <= 1.5 * math.log(l, 2):

            return dacs(l, r0, r2, r2 + r0, chain + [1]) + dacs(
                l, r1, r2, r2 + r1, chain + [0]
            )
        else:
            return []

    def sdac(l):
        """
        sdac()
        input: a small odd prime number l
        output: the shortest differential additions chains corresponding with the input l

        NOTE: this function uses a recursive function
        """
        all_dacs = dacs(l, 1, 2, 3, [])
        return min(all_dacs, key=lambda t: len(t[0]))[0]

    def generate_sdacs(L):
        """ Shortest Differential Addition Chains for each small odd prime l in L """
        return list(
            map(sdac, L)
        )

    # Shortest Differential Addition Chains (SDACs) for each l_i
    path = resource_filename('sibc', "data/sdacs/" + algorithm + '/' + prime)
    print("// Computing sdacs")
    SDACS = generate_sdacs(L)
    print("// Storing sdacs into a file")
    write_list_of_lists_of_ints_to_file(path, SDACS)

if __name__ == "__main__":
    main()