import os
from functools import reduce
from .math import bitlength, is_prime

base_path = "/usr/share/python3-sidh/data/"
if not os.path.exists(base_path) and os.path.exists('./data'):
    # this allows running it from the repo without installing on the system, for now
    # FIXME: use pkg_resources to locate data
    base_path = "./data/"
strategy_data = base_path + "/strategies/"
sdacs_data = base_path + "/sdacs/"
sop_data = base_path + "/sop/"
ijk_data = base_path + "/ijk/"
gen_data = base_path + "/gen/"
tmp_dir = "./"  # Drop the files in the current working directory


def csidh_get_sop_from_disk(prime):
    assert prime in ('p512', 'p1024', 'p1792'), "unsupported prime for csidh"
    # List of Small odd primes, L := [l_0, ..., l_{n-1}]
    f = open(sop_data + prime)
    L = f.read()
    L = [int(l) for l in L.split()]
    exponent_of_two = L[0]  #   Exponent of cofactor-2
    L = list(L[1:])  #   Small Odd Primes l_i's
    n = len(L)  #   Number of l_i's to be used
    f.close()
    p = (2 ** (exponent_of_two)) * reduce(
        lambda x, y: (x * y), L
    ) - 1  # p := 4 * l_0 * ... * l_n - 1
    p_minus_one_halves = (p - 1) // 2  # (p - 1) / 2

    validation_stop = sum([bitlength(l_i) for l_i in L]) / 2.0 + 2
    #    assert is_prime(p), "[ERROR]\tThe integer number p := 4 * l_1, * ... * l_n - 1 is not prime where L := %s" % (L,) #FIXME
    return dict(
        L=L,
        exponent_of_two=exponent_of_two,
        n=n,
        p=p,
        p_minus_one_halves=p_minus_one_halves,
        validation_stop=validation_stop,
        p_bits=int(prime[1:]),
    )


#   """
#   if( (sys.argv[0] != 'header.py') ):
#       print("/*")
#       print("The prime number to be used has the following form\n")
#       print("           %3d" % n)
#       print("        ~~~~~~~~~")
#       print("         |     | ")
#       print("p := 4 x |     | l_j   -   1, where each l_j is a small odd prime")
#       print("           j=1  ")
#       print("*/\n")
#   """


def bsidh_get_sop_from_disk(prime):
    f = open(sop_data + prime)

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
    f.close()
    # pp = (2**exponent_of_twop) * reduce(lambda x,y : (x*y), [ Lp[i]**Ep[i] for i in range(0, np, 1)  ])
    pp = reduce(
        lambda x, y: (x * y), [Lp[i] ** Ep[i] for i in range(0, np, 1)]
    )
    pm = reduce(
        lambda x, y: (x * y), [Lm[i] ** Em[i] for i in range(0, nm, 1)]
    )
    assert (p + 1) % pp == 0
    assert (p - 1) % pm == 0

    if p % 4 == 1:
        print("// Case p = 1 mod 4 is not implemented yet!")
        exit(-1)

    p_minus_one_halves = (p - 1) // 2
    p_minus_3_quarters = (p - 3) // 4
    return  # FIXME do the bsidh part


parameters = dict(
    csidh=dict(
        A=[2, 4],
        p512=dict(
            wd1=dict(
                m=
                # fmt: off
                    [ 15, 18, 20, 21, 21, 22, 22, 22, 22, 22, 22, 19, 20, 22, 23, 23,
                    23, 23, 23, 23, 23, 21, 23, 20, 16, 16, 16, 15, 14, 12, 13, 12,
                    11, 11, 10, 10, 9, 9, 9, 8, 8, 8, 8, 7, 7, 7, 6, 6, 6, 6, 6, 6,
                    6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    3, ]
                ,
                # fmt: on
                sigma=5,
                kappa=11,
                delta=1,
            ),
            wd2=dict(
                m=
                # fmt: off
                    [ 7, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11,
                    11, 11, 11, 11, 9, 11, 9, 8, 8, 8, 7, 7, 7, 7, 6, 6, 5, 5, 5,
                    5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                    3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, ]
                ,
                # fmt: on
                sigma=3,
                kappa=8,
                delta=2,
            ),
            df=dict(
                m=
                # fmt: off
                    [ 15, 18, 20, 21, 21, 22, 22, 22, 22, 22, 22, 19, 20, 22, 23, 23,
                    23, 23, 23, 23, 23, 23, 23, 19, 16, 16, 16, 15, 14, 12, 13, 12,
                    11, 11, 11, 9, 9, 9, 9, 8, 8, 8, 8, 7, 8, 6, 6, 6, 6, 7, 6, 6,
                    6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                    3, ]
                ,
                # fmt: on
                sigma=5,
                kappa=11,
                delta=1,
            ),
            **csidh_get_sop_from_disk('p512')
        ),
        p1024=dict(
            wd1=dict(
                m=
                # fmt: off
                    [ 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                    5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3,
                    3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, ]
                ,
                # fmt: on
                sigma=5,
                kappa=4,
                delta=1,
            ),
            wd2=dict(
                m=
                # fmt: off
                    [ 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, ]
                ,
                # fmt: on
                sigma=3,
                kappa=5,
                delta=2,
            ),
            df=dict(
                m=
                # fmt: off
                    [ 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                    5, 6, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3,
                    3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, ]
                ,
                # fmt: on
                sigma=5,
                kappa=4,
                delta=1,
            ),
            **csidh_get_sop_from_disk('p1024')
        ),
        p1792=dict(
            wd1=dict(
                m=
                # fmt: off
                    [ 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                    5, 5, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ]
                ,
                # fmt: on
                sigma=5,
                kappa=4,
                delta=1,
            ),
            wd2=dict(
                m=[1]
                * 207,  # 207 is the n value returned by csidh_get_sop_from_disk('p1792')
                sigma=3,
                kappa=5,
                delta=2,
            ),
            df=dict(
                m=
                # fmt: off
                    [ 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                    5, 5, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ]
                ,
                # fmt: on
                sigma=5,
                kappa=4,
                delta=1,
            ),
            **csidh_get_sop_from_disk('p1792')
        ),
    ),
    bsidh=dict(),  # FIXME
)
