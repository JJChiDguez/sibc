# **S**upersingular **I**sogeny-**B**ased **C**ryptographic constructions

<img align="center" src="https://github.com/JJChiDguez/sibc/blob/master/sibc-logo.png">

This repository includes a python-code library named **sibc**, which allows a user-friendly
interface to deal with some isogeny-based cryptographic primitive.

The current version of the **sibc** library has integrated SIDH, CSIDH, and B-SIDH schemes using traditional
and velusqrt formulae on Montgomery curve x-only projective coordinates; in particular, **sibc** allows
working with prime and quadratic field classes that permit operating field elements as integers.
Additionally, the cryptographic primitives are implemented in constant-time concerning the number
of field operations, where a constant-time procedure refers to its running time does not depend on
the input or it possibly does from randomness as CSIDH does.

It is worthing to mention, the library is constantly extended, and some signature schemes will be integrated into the **sibc** library.

The **sibc** library aims to allow isogeny-contributors for building new primitives with a constant-time nature.


> :warning: There is a new devastating attack against *SIDH* and *SIKE* by Castryck & Decru. Currently, there are two public implementations of the Castryck-Decru attack:
> 
> 1. [**Magma** code](https://homes.esat.kuleuven.be/~wcastryc/) from [Castryck-Decru preprint](https://eprint.iacr.org/2022/975), and
> 2. [**Sagemath** code](https://github.com/jack4818/Castryck-Decru-SageMath) translation from the Magma code, by Giacomo Pope.
>
> :warning: The attack does extend to B-SIDH and B-SIKE.
> 
> :exclamation: The attack does not apply to CSIDH.

## Installation 

Install the `sibc` module which provides the `sibc` program:

```
sudo python3 setup.py install
```

## For development:
For this installation method, any further modification in sibc directory will be reflect when running sibc library.

```
# Installing required package
# Before running the following commands, ensure you have the lastest version of pip
pip3 install dh click numpy progress matplotlib networkx stdeb setuptools-scm setuptools

# only pip3 install cpuinfo is missing for macOS (to be fixed in coming versions)
pip3 install pytest pytest-xdist

# Installing the library
sudo pip3 install -e .
```

### Debian package build

To build a package for Debian or Ubuntu, we suggest the use of stdeb:

```
sudo apt install -y dh-python python3-click python3-progress\
  python3-numpy python3-matplotlib python3-networkx \
  python3-stdeb python3-setuptools-scm python3-setuptools python3-cpuinfo
python3 setup.py bdist_deb
sudo dpkg -i deb_dist/python3-sibc_0.0.1-1_all.deb
```

## Usage 
The syntax compilation can be viewed by running one of the following three commands:
```bash
# Corresponding with the key-exchange protocol
sibc --help

# Corresponding with benchmarking (only for CSIDH, which has a variable
# running-time cost independent from the key)
sibc csidh-bench

# Corresponding with the costs of KPs (Kernel Point computation), xISOG
# (isogeny construction), and xEVAL (isogeny evaluation)
sibc csidh-test
```

Usage for the `sibc` tool:
```bash

Usage: sibc [OPTIONS] COMMAND [ARGS]...

    ,-~~-.___.          
   / |  '     \        
  (  )         0        
   \_/-, ,----'         
      ====           // 
     /  \-'~;    /~~~(O)
    /  __/~|   /       |
  =(  _____| (_________|

Options:
  -p, --prime [p434|p503|p610|p751|p253|p255|p247|p237|p257|p512|p1024|p1792|p2048|p4096|p5120|p6144|p8192|p9216]
                                  [default: p512]
  -f, --formula [tvelu|svelu|hvelu]
                                  [default: hvelu]
  -a, --algorithm [sidh|sike|csidh|bsidh|bsike]
                                  [default: csidh]
  -s, --style [wd1|wd2|df]        [default: df]
  -e, --exponent [1|2|3|4|5|6|7|8|9|10]
                                  [default: 10]
  -m, --multievaluation           [default: False]
  -c, --curvemodel [edwards|montgomery]
                                  [default: montgomery]
  -b, --benchmark INTEGER         [default: 128]
  -t, --tuned                     [default: False]
  -u, --uninitialized             [default: False]
  -v, --verbose                   Not the kind of verbosity you might expect
                                  [default: False]
  --version                       Show the version and exit.
  --help                          Show this message and exit.

Commands:
  bsidh-main                   Random instance example of a key-exchange
  bsidh-precompute-parameters  Precomputation of tuned velusqrt parameters
  bsidh-precompute-strategy    Precomputation of optimal strategies
  bsidh-test                   GF(p²)-operation cost of kps, xisog, and...
  csidh-bench                  Average GF(p)-operation cost of a GAE
  csidh-bounds                 Greedy-based search of optimal exponents
  csidh-dh                     Derive shared secret key from CSIDH sk,...
  csidh-genkey                 Generate random CSIDH secret key
  csidh-header                 Optimal strategies as C-code headers files
  csidh-ijk                    Velusqrt parameters as C-code headers files
  csidh-main                   Random instance example of a key-exchange
  csidh-precompute-parameters  Precomputation of tuned velusqrt parameters
  csidh-precompute-strategy    Precomputation of optimal strategies
  csidh-pubkey                 Derive CSIDH public key from CSIDH secret key
  csidh-sdacs                  SDACs as C-code headers files
  csidh-test                   GF(p)-operation cost of kps, xisog, and xeval
  decaps                       (B)SIKE decapsulation
  encaps                       (B)SIKE encapsulation
  keygen                       Generate random (B)SIKE secret and public...
  plot-strategy                draw strategy graphs as a subgraph...
  print-timing
  sidh-precompute-strategy     Precomputation of optimal strategies
```

## SIDH cryptographic API

`CSIDH`, `BSIDH`, `SIDH`, `SIKE`, and `BSIKE` objects are available from the `sibc` package and module.

Automatically generated documentation is available with pydoc after `sibc` is
installed:
```
pydoc3 sibc.csidh
pydoc3 sibc.bsidh
pydoc3 sibc.sidh
```

### Command Line Interface: examples

```bash
# CSIDH
sk_a="$(sibc csidh-genkey)"
pk_a="$(echo "$sk_a"|sibc csidh-pubkey -)"
sk_b="$(sibc csidh-genkey)"
pk_b="$(echo "$sk_b"|sibc csidh-pubkey -)"
ss_a="$(echo "$sk_a"|sibc csidh-dh - "$pk_b")"
ss_b="$(echo "$sk_b"|sibc csidh-dh - "$pk_a")"
echo $ss_a
echo $ss_b

# SIKE
sk="$(sibc -a sidh -p p434 keygen)"
pk3=`echo "${sk}" | tail -n1`
ck="$(echo "$pk3"|sibc -a sidh -p p434 encaps -)"
c0=`echo "${ck}" | head -1`
c1=`echo "${ck}" | tail -2 | head -1`
K=`echo "${ck}" | tail -n1`
K_="$(echo "$sk"|sibc -a sidh -p p434 decaps - "$c0 $c1")"
echo $K
echo $K_

# BSIKE
sk="$(sibc -a bsidh -p p253 keygen)"
pk3=`echo "${sk}" | tail -n1`
ck="$(echo "$pk3"|sibc -a bsidh -p p253 encaps -)"
c0=`echo "${ck}" | head -1`
c1=`echo "${ck}" | tail -2 | head -1`
K=`echo "${ck}" | tail -n1`
K_="$(echo "$sk"|sibc -a bsidh -p p253 decaps - "$c0 $c1")"
echo $K
echo $K_
```

### Basic shared secret generation example with CSIDH
```python3
from sibc.csidh import CSIDH, default_parameters
csidh = CSIDH(**default_parameters)

# alice generates a key
alice_secret_key = csidh.secret_key()
alice_public_key = csidh.public_key(alice_secret_key)

# bob generates a key
bob_secret_key = csidh.secret_key()
bob_public_key = csidh.public_key(bob_secret_key)

# if either alice or bob use their secret key with the other's respective
# public key, the resulting shared secrets are the same
shared_secret_alice = csidh.dh(alice_secret_key, bob_public_key)
shared_secret_bob = csidh.dh(bob_secret_key, alice_public_key)

# Alice and bob produce an identical shared secret
assert shared_secret_alice == shared_secret_bob
```

### Basic shared secret generation example with BSIDH
```python3
from sibc.bsidh import BSIDH, default_parameters
bsidh = BSIDH(**default_parameters)
sk_a, pk_a = bsidh.keygen_a()
sk_b, pk_b = bsidh.keygen_b()
ss_a, ss_b = bsidh.derive_a(sk_a, pk_b), bsidh.derive_b(sk_b, pk_a)
ss_a == ss_b
```

### Basic example with BSIKE (BSIDH + key encapsulation)
```python3
from sibc.bsidh import BSIKE, default_parameters
bsike = BSIKE(**default_parameters)
s, sk3, pk3 = bsike.KeyGen()
c, K = bsike.Encaps(pk3)
K_ = bsike.Decaps((s, sk3, pk3), c)
K == K_

bsike255 = BSIKE('montgomery', 'p255', 'hvelu', True, False, False, False)
s, sk3, pk3 = bsike255.KeyGen()
c, K = bsike255.Encaps(pk3)
K_ = bsike255.Decaps((s, sk3, pk3), c)
K == K_
```

### Basic shared secret generation example with SIDH
```python3
from sibc.sidh import SIDH, default_parameters
sidh = SIDH(**default_parameters)
sk_a, pk_a = sidh.keygen_a()
sk_b, pk_b = sidh.keygen_b()
ss_a, ss_b = sidh.dh_a(sk_a, pk_b), sidh.dh_b(sk_b, pk_a)
ss_a == ss_b
```

### Basic example with SIKE (SIDH + key encapsulation)
```python3
from sibc.sidh import SIKE, default_parameters
sike = SIKE(**default_parameters)
s, sk3, pk3 = sike.KeyGen()
c, K = sike.Encaps(pk3)
K_ = sike.Decaps((s, sk3, pk3), c)
K == K_

sike503 = SIKE('montgomery', 'p503', False, False)
s, sk3, pk3 = sike503.KeyGen()
c, K = sike503.Encaps(pk3)
K_ = sike503.Decaps((s, sk3, pk3), c)
K == K_
```

## Adding new primes

The field characteristic `p` should be stored in directory `data/sop/`, and
CSIDH and BSIDH have different structures (see below):

```bash
# CSIDH format (here p = cofactor * l_1 * .... l_n - 1)
cofactor l_1 l_2 ... l_n

# BSIDH format
Hexadecimal representation of the prime p
4 l_1 l_2 ... l_n
c e_1 e_2 ... e_n
l'_1 l'_2 ... l'_m
e'_1 e'_2 ... e'_m

# SIDH format: p = 2^{e_2} * 3^{e_3} - 1
e_2 e_3
```

For the case of BSIDH, `M := (4^c * l_1^{e_1} * l_2^{e_2} * ... * l_n^{e_n})`
must divide `(p + 1)`, and `N := (l'_1^{e'_1} * l'_2^{e'_2} * ... *
l'_n^{e'_n})` must divide `(p-1)`. Additionally, the order-`M` generators `PA`,
`QA` and `PQA := PA - QA` should be stored in directory `gen/` as projective
x-coordinate points. Similarly, the order-`N` generators `PB`, `QB` and `PQB :=
PB - QB` also should be stored it the same directory. Both 3-tuples of points
must be stored in a single file with the following syntax:

```bash
Re(x(PA)) Im(x(PA)) Re(x(QA)) Im(x(QA)) Re(x(PQA)) Im(x(PQA))
Re(x(PB)) Im(x(PB)) Re(x(QB)) Im(x(QB)) Re(x(PQB)) Im(x(PQB))
```

where `Re(X)` and `Im(X)` denote the real and imaginary parts of X with respect
to `F_p[i]/(i^2 + 1)`, respectively. Moreover, all the above twelve integers
should be stored in hexadecimal.

For SIDH, generators have order either `M=2^{e_2}` or `N=3^{e_3}`.

## Examples

We summarize some examples of runs of the `sibc` tool as follows:

```bash
# CSIDH
# A single random intances of a key exchange
sibc -p p512 -f hvelu -a csidh -s df -e 10 csidh-main
sibc -p p512 -f hvelu -a csidh -s df -e 10 -m csidh-main
sibc -p p512 -f hvelu -a csidh -s df -e 10 -t csidh-main
sibc -p p512 -f hvelu -a csidh -s df -e 10 -m -t csidh-main
# Average GF(p)-operation cost of 64 random instances
sibc -p p512 -f hvelu -a csidh -s df -m -e 10 -b 64 csidh-bench
sibc -p p512 -f hvelu -a csidh -s df -m -e 10 -b 64 -m csidh-bench
sibc -p p512 -f hvelu -a csidh -s df -m -e 10 -b 64 -t csidh-bench
sibc -p p512 -f hvelu -a csidh -s df -m -e 10 -b 64 -t -m csidh-bench
# GF(p)-operation cost of kps, xisog, and xeval blocks
sibc -p p512 -f hvelu -a csidh csidh-test
sibc -p p512 -f hvelu -a csidh -m csidh-test
sibc -p p512 -f hvelu -a csidh -t csidh-test
sibc -p p512 -f hvelu -a csidh -t -m csidh-test

# BSIDH
# A single random intances of a key exchange
sibc -p p253 -f hvelu -a bsidh bsidh-main
sibc -p p253 -f hvelu -a bsidh -m bsidh-main
sibc -p p253 -f hvelu -a bsidh -t bsidh-main
sibc -p p253 -f hvelu -a bsidh -t -m bsidh-main
# GF(p²)-operation cost of kps, xisog, and xeval blocks
sibc -p p253 -f tvelu -a bsidh bsidh-test
sibc -p p253 -f svelu -a bsidh bsidh-test
sibc -p p253 -f hvelu -a bsidh -t bsidh-test
```

Remark, our implementation allows us to plot each optimal strategy required (only tested in Linux machines):

```bash
# CSIDH
sibc -p p512 -f tvelu -a csidh -s df -e 10 plot-strategy
sibc -p p512 -f svelu -a csidh -s wd1 -e 10 plot-strategy
sibc -p p512 -f hvelu -a csidh -s wd2 -e 5 plot-strategy

# BSIDH
sibc -p p253 -f hvelu -a bsidh plot-strategy
sibc -p p253 -f hvelu -a bsidh -m plot-strategy
sibc -p p253 -f hvelu -a bsidh -t -m plot-strategy

# SIDH
sibc -p p434 -a sidh plot-strategy
sibc -p p503 -a sidh plot-strategy
sibc -p p610 -a sidh plot-strategy
sibc -p p751 -a sidh plot-strategy
```

Additionally, one can created files with extension `.h` that includes all the
required variables in a the sdacs, strategies, and velusqrt (at least for CSIDH
implementations).

```bash
# Suitable bounds search with e = 10.
sibc -a csidh -p p512 -s df -f hvelu -e 10 -u csidh-bounds # The greedy-based algorithm on a large searching space, it could take hours or even days!: option -u is required
# SDACs (options -s and -e do not affect the output)
sibc -p p512 -f hvelu -a csidh -s df -e 10 csidh-sdacs
# Optimal sizes of I, J, and K required in velusqrt (options -s and -e do not affect the output)
sibc -p p512 -f hvelu -a csidh -s df -e 10 -t csidh-ijk  # option -t is required
# Optimal strategies
sibc -p p512 -f hvelu -a csidh -s df -e 10 -t csidh-header
```

## BSIDH primes
Currently only `p253`, `p255`, `p247`, `p237`, and `p257` are implemented and tested in the current API.
Extending this to other primes is straight-forward.

## SIDH primes
Currently only `p434`, `p503`, `p610`, and `p751` are implemented and tested in the current API.
Extending this to other primes is straight-forward.

## Precomputing data for a new prime instances

Generating new data can be easily done by adding and running to either `misc/create-csidh-data.sh` or
`misc/create-bsidh-data.sh`. The new prime number description should b stored as previously mentioned.

```bash
bash misc/create-csidh-data.sh
bash misc/create-bsidh-data.sh
```

Also, you can do it manually by doing something as follows:

```bash
# CSIDH
sibc-precompute-sdacs -p p512 -a csidh # SDACs
sibc -p p512 -f hvelu -a csidh -m -t -u csidh-precompute-parameters # Tuned velusqrt parameters
sibc -p p512 -f hvelu -a csidh -s df -m -u csidh-precompute-strategy # Strategies
# BSIDH
sibc-precompute-sdacs -p p253 -a bsidh # SDACs
sudo sibc -p p253 -f svelu -a bsidh -u bsidh-precompute-parameters # Tuned velusqrt parameters: the option -u is required
sudo sibc -p p253 -f svelu -a bsidh -u bsidh-precompute-strategy # Strategies
# SIDH
sibc -p p434 -a sidh -u sidh-precompute-strategy # Strategies
sibc -p p503 -a sidh -u sidh-precompute-strategy # Strategies
sibc -p p610 -a sidh -u sidh-precompute-strategy # Strategies
sibc -p p751 -a sidh -u sidh-precompute-strategy # Strategies
```

Furthermore, you can create tests by running `bash misc/create-tests.sh` and `bash misc/test-cli.sh`
(only `CSIDH`, `BSIKE`, and `SIKE` are handled by now).

If either a new prime instace or primitive is included, then you should add it to misc directory.
New primitives require new bash scripts. SIDH instances are simple, only one configuration (for now),
and thus you can omit adding them in `misc` directory.
 
## Remarks

The primes labeled as `p253`, `p255`, `p247`, and `p237` correspond with the examples 2, 3, 5, and 6 from
[B-SIDH paper](https://eprint.iacr.org/2019/1145), respectively. In particular,  `p257` denotes the
prime number given in [velusqrt paper](https://eprint.iacr.org/2020/341). The field airthmetic is
centered on primes `p = 3 mod 4`. Multiplying and squaring in `GF(p²) = GF(p)[u]/(u^2 + 1)` have a
cost of 3 and 2 multiplications in `GF(p)`, respectively.

The current implementation does not have implemented the B-SIDH key validation, it will be included
in the next library version.

### Adding new prime instances
When adding a new isogeny-based instances (the prime number, and public parameters) should be included
in parameter list of `sibc/__main__.py` (click option `-p`, `--prime`) by modifying `sibc/constants.py`
(updating the lists/dictionary `[csidh/bsidh]_primes`). If a new primitive is included, then you need
to update `sibc/__main__.py` by extending the click options `-p` (`--prime`) and `-a` (`--algorithm`),
and also to include its branch in `sibc/montgomery/curve.py` and `sibc/montgomery/isogeny.py` files.

## Changes

Significant changes are listed in the [CHANGELOG](CHANGELOG.md) file. Future integrations/modifications
are listed in the [TODOLIST](TODOLIST.md) file.

## Authors

1. **_Gora Adj_** <gora.adj@udl.cat>, <gora.adj@gmail.com>;
2. **_Jesús-Javier Chi-Domínguez_** <jesus.dominguez@tii.ae>, <chidoys@gmail.com>; and
3. **_Francisco Rodríguez-Henríquez_** <francisco@cs.cinvestav.mx>.

### Main cryptographer contributors (collaborators):

1. **_Jorge Ch&aacute;vez-Saab_** <jchavez@computacion.cs.cinvestav.mx>, <jorgechavezsaab@gmail.com>

### Logo creator contributor

1. **_Fabiola-Argentina Hern&aacute;ndez-Torres_** <farg.cls@outlook.com>

### User-friendly interface contributors

1. Jacob Appelbaum <j.appelbaum@tue.nl>; and
2. Leif Ryge <leif@synthesize.us>.



All contributors are listed in the [CONTRIBUTORS](CONTRIBUTORS) file.

### How to cite this library

```latex
@misc{sibc,
  author = {{Gora Adj} and {Jes\'us-Javier Chi-Dom\'inguez} and {Francisco Rodr\'iguez-Henr\'iquez}},
  title = {{SIBC} Python library},
  year = {2021},
  howpublished = {\url{https://github.com/JJChiDguez/sibc/}}
}
```

## License

This project is licensed under the GNU general public license - see the
[LICENSE](LICENSE) file for details.

## Funding

This project has initially received funding from the European Research Council (ERC) under the
European Union's Horizon 2020 research and innovation programme (grant agreement No 804476),
while the second author was doing a postdoctoral research stay in Tampere University.

The third author received partial funds from the Mexican Science council CONACyT project 313572.
