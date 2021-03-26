import os
from functools import reduce
from .math import bitlength, is_prime
from pkg_resources import resource_string, resource_filename

def csidh_get_sop_from_disk(prime):
    assert prime in ('p512', 'p1024', 'p1792'), "unsupported prime for csidh"
    # List of Small odd primes, L := [l_0, ..., l_{n-1}]
    L = resource_string(__name__, "data/sop/" + prime)
    L = [int(l) for l in L.split()]
    cofactor = L[0]
    L = list(L[1:])  #   Small Odd Primes l_i's
    n = len(L)  #   Number of l_i's to be used
    p = cofactor * reduce(
        lambda x, y: (x * y), L
    ) - 1  # p := cofactor * l_0 * ... * l_n - 1
    p_minus_one_halves = (p - 1) // 2  # (p - 1) / 2

    validation_stop = sum([bitlength(l_i) for l_i in L]) / 2.0 + 2
    return dict(
        L=L,
        cofactor=cofactor,
        n=n,
        p=p,
        p_minus_one_halves=p_minus_one_halves,
        validation_stop=validation_stop,
        p_bits=bitlength(p),
    )

def csidh_get_exp_from_disk(exponent, style, security, attack, n):

    m = resource_string(__name__, "data/exponents/" + attack + '/' + security + '/' + style + '/e' + exponent)
    m  = [int(mi) for mi in m.split()]
    assert len(m) <= n, 'Not enough number of small odd prime numbers in the factorization of p + 1'
    return m + [0]*(n - len(m))

def bsidh_get_sop_from_disk(prime):
    f = open(resource_filename('sibc', 'data/sop/'+ prime))

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
    return dict(
        nm=nm,
        np=np,
        Lm=Lm,
        Lp=Lp,
        Ep=Ep,
        Em=Em,
        pp=pp,
        pm=pm,
        p=p,
        p_minus_one_halves=p_minus_one_halves,
        p_bits=bitlength(p),
        p_minus_3_quarters=p_minus_3_quarters,
    )


parameters = dict(
    csidh=dict(
        p512=dict(
            wd2=dict(
                m=csidh_get_exp_from_disk('5', 'wd2', '128', 'mitm', 74),
            ),
            wd1=dict(
                m=csidh_get_exp_from_disk('10', 'wd1', '128', 'mitm', 74),
            ),
            df=dict(
                m=csidh_get_exp_from_disk('10', 'df', '128', 'mitm', 74),
            ),
            **csidh_get_sop_from_disk('p512')
        ),
        p1024=dict(
            wd2=dict(
                m=csidh_get_exp_from_disk('2', 'wd2', '128', 'mitm', 130),
            ),
            wd1=dict(
                m=csidh_get_exp_from_disk('3', 'wd1', '128', 'mitm', 130),
            ),
            df=dict(
                m=csidh_get_exp_from_disk('3', 'df', '128', 'mitm', 130),
            ),
            **csidh_get_sop_from_disk('p1024')
        ),
        p1792=dict(
            wd2=dict(
                m=csidh_get_exp_from_disk('1', 'wd2', '128', 'mitm', 207),
            ),
            wd1=dict(
                m=csidh_get_exp_from_disk('2', 'wd1', '128', 'mitm', 207),
            ),
            df=dict(
                m=csidh_get_exp_from_disk('2', 'df', '128', 'mitm', 207),
            ),
            **csidh_get_sop_from_disk('p1792')
        ),
    ),
    bsidh=dict(
        b2=dict(
            **bsidh_get_sop_from_disk('b2')
        ),
        b3=dict(
            **bsidh_get_sop_from_disk('b3')
        ),
        b5=dict(
            **bsidh_get_sop_from_disk('b5')
        ),
        b6=dict(
            **bsidh_get_sop_from_disk('b6')
        ),
        s1=dict(
            **bsidh_get_sop_from_disk('s1')
        ),
    ),
)
