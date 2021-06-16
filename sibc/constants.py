import os
from functools import reduce
from .math import bitlength, is_prime
from pkg_resources import resource_string, resource_filename

# When including a new csidh-prime, you need to update the next list
# The csidh format in the following dictionary is
# PRIME:[
#   EXPONENT,
#   STYLE,
#   SECURITY (KEYSPACE SIZE IN BITS),
#   ATTACK (either mitm or vow-gcs),
#   NUMBER n OF SMALL ODD PRIMES IN THE FACTORIZATION OF p + 1
# ]
csidh_primes = {
    'p512':[
        ('5', 'wd2', '128', 'mitm', 74),
        ('10', 'wd1', '128', 'mitm', 74),
        ('10', 'df', '128', 'mitm', 74)
    ],
    'p1024':[
        ('2', 'wd2', '128', 'mitm', 130),
        ('3', 'wd1', '128', 'mitm', 130),
        ('3', 'df', '128', 'mitm', 130)
    ],
    'p1792':[
        ('1', 'wd2', '128', 'mitm', 207),
        ('2', 'wd1', '128', 'mitm', 207),
        ('2', 'df', '128', 'mitm', 207)
    ],
    'p2048':[
        ('1', 'wd2', '128', 'vow-gcs', 221),
        ('1', 'wd1', '128', 'vow-gcs', 221),
        ('1', 'df', '128', 'vow-gcs', 221)
    ],
    'p4096':[
        ('1', 'wd2', '128', 'vow-gcs', 221),
        ('1', 'wd1', '128', 'vow-gcs', 221),
        ('1', 'df', '128', 'vow-gcs', 221)
    ],
    'p5120':[
        ('1', 'wd2', '128', 'mitm', 256),
        ('1', 'wd1', '128', 'mitm', 256),
        ('1', 'df', '128', 'mitm', 256)
    ],
    'p6144':[
        ('1', 'wd2', '192', 'vow-gcs', 306),
        ('1', 'wd1', '192', 'vow-gcs', 306),
        ('1', 'df', '192', 'vow-gcs', 306)
    ],
    'p8192':[
        ('1', 'wd2', '192', 'vow-gcs', 306),
        ('1', 'wd1', '192', 'vow-gcs', 306),
        ('1', 'df', '192', 'vow-gcs', 306)
    ],
    'p9216':[
        ('1', 'wd2', '192', 'mitm', 384),
        ('1', 'wd1', '192', 'mitm', 384),
        ('1', 'df', '192', 'mitm', 384)
    ]
}

# When including a new bsidh-prime, you need to update the next list
bsidh_primes = (
    'p253',
    'p255',
    'p247',
    'p237',
    'p257'
)

# When including a new [sidh/sike]-prime, you need to update the next list
sidh_primes = (
    'p434',
    'p503',
    'p610',
    'p751'
)

def csidh_get_sop_from_disk(prime):
    assert prime in csidh_primes.keys(), "unsupported prime for csidh"
    # List of Small odd primes, L := [l_0, ..., l_{n-1}]
    L = resource_string(__name__, "data/sop/csidh/" + prime)
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
    assert prime in bsidh_primes, "unsupported prime for bsidh"
    f = open(resource_filename('sibc', 'data/sop/bsidh/'+ prime))

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

def sidh_get_sop_from_disk(prime):
    assert prime in sidh_primes, "unsupported prime for sidh/sike"
    f = open(resource_filename('sibc', 'data/sop/sidh/'+ prime))

    # p + 1 = 2^two x 3^three
    Exponents = f.readline()
    two, three = [int(e) for e in Exponents.split()]
    p = (2 ** two) * (3 ** three) - 1

    if p % 4 == 1:
        print("// Case p = 1 mod 4 is not implemented yet!")
        exit(-1)

    p_minus_one_halves = (p - 1) // 2
    p_minus_3_quarters = (p - 3) // 4
    return dict(
        two=two,
        three=three,
        p=p,
        p_minus_one_halves=p_minus_one_halves,
        p_bits=bitlength(p),
        p_minus_3_quarters=p_minus_3_quarters,
    )

parameters = dict(
    csidh=dict({
        PRIME:{
            **{style:{
                'm':csidh_get_exp_from_disk(exp, style, sec, attack, n)
                } for (exp, style, sec, attack, n) in csidh_primes[PRIME]
            },
            **csidh_get_sop_from_disk(PRIME)
        } for PRIME in csidh_primes.keys()
    }),
    bsidh=dict({
        PRIME:dict(**bsidh_get_sop_from_disk(PRIME)) for PRIME in bsidh_primes
    }),
    sidh=dict({
        PRIME:dict(**sidh_get_sop_from_disk(PRIME)) for PRIME in sidh_primes
    })
)