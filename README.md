# **S**upersingular **I**sogeny-**B**ased **C**ryptography constructions

This repository includes a python-code library named **sibc**, which allows an user-
friendly usage to deal with some isogeny-based cryptographic primitive.

The current version of **sibc** library has integrated CSIDH and B-SIDH schemes using
traditional and velusqrt formulae on Montgomery curve x-only projective coordinates.

The current version allows working with prime and quadratic field classes that permit
operating field elements as integers. Moreover, the current cryptographic primitives
are implemented in constant-time concerning the number of field operations. Here, a
constant-time algorithm means its running time does not depend on the input or it possibly
does from randomness as CSIDH does.

## Installation 

Install the `sibc` module which provides the `sibc` program:
```
sudo python3 setup.py install
```

For development:
```
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
  -p, --prime [b2|b3|b5|b6|p1024|p1792|p512|s1]
                                  [default: p512]
  -f, --formula [tvelu|svelu|hvelu]
                                  [default: hvelu]
  -a, --algorithm [csidh|bsidh]   [default: csidh]
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
  csidh-dh                     Derive shared secret key from CSIDH sk,
                               CSIDH...
  csidh-genkey                 Generate random CSIDH secret key
  csidh-header                 Optimal strategies as C-code headers files
  csidh-ijk                    Velusqrt parameters as C-code headers files
  csidh-main                   Random instance example of a key-exchange
  csidh-precompute-parameters  Precomputation of tuned velusqrt parameters
  csidh-precompute-strategy    Precomputation of optimal strategies
  csidh-pubkey                 Derive CSIDH public key from CSIDH secret key
  csidh-sdacs                  SDACs as C-code headers files
  csidh-test                   GF(p)-operation cost of kps, xisog, and xeval
  plot-strategy                draw strategy graphs as a subgraph Discrete...
  print-timing
```

## SIDH cryptographic API

CSIDH and BSIDH objects are available from the `sibc` package and module.

Automatically generated documentation is available with pydoc after `sibc` is
installed:
```
pydoc3 sibc.csidh
pydoc3 sibc.bsidh
```

### Basic shared secret generation example with CSIDH
```python3
from sibc.csidh import CSIDH, default_parameters
c = CSIDH(**default_parameters)

# alice generates a key
alice_secret_key = c.secret_key()
alice_public_key = c.public_key(alice_secret_key)

# bob generates a key
bob_secret_key = c.secret_key()
bob_public_key = c.public_key(bob_secret_key)

# if either alice or bob use their secret key with the other's respective
# public key, the resulting shared secrets are the same
shared_secret_alice = c.dh(alice_secret_key, bob_public_key)
shared_secret_bob = c.dh(bob_secret_key, alice_public_key)

# Alice and bob produce an identical shared secret
assert shared_secret_alice == shared_secret_bob
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
should be stored in hexadecimal."

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
sibc -p b2 -f hvelu -a bsidh bsidh-main
sibc -p b2 -f hvelu -a bsidh -m bsidh-main
sibc -p b2 -f hvelu -a bsidh -t bsidh-main
sibc -p b2 -f hvelu -a bsidh -t -m bsidh-main
# GF(p²)-operation cost of kps, xisog, and xeval blocks
sibc -p b2 -f tvelu -a bsidh bsidh-test
sibc -p b2 -f svelu -a bsidh bsidh-test
sibc -p b2 -f hvelu -a bsidh -t bsidh-test
```

Remark, our implementation allows us to plot each optimal strategy required:

```bash
# CSIDH
sibc -p p512 -f tvelu -a csidh -s df -e 10 plot-strategy
sibc -p p512 -f svelu -a csidh -s wd1 -e 10 plot-strategy
sibc -p p512 -f hvelu -a csidh -s wd2 -e 5 plot-strategy

# BSIDH
sibc -p b2 -f hvelu -a bsidh plot-strategy
sibc -p b2 -f hvelu -a bsidh -m plot-strategy
sibc -p b2 -f hvelu -a bsidh -t -m plot-strategy
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
Currently only `b2`, `b3`, `b5`, `b6`, and `s1` are implemented and tested in the current API.
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
sibc -p p512 -f hvelu -a csidh -m -t csidh-precompute-parameters # Tuned velusqrt parameters
sibc -p p512 -f hvelu -a csidh -s df -m csidh-precompute-strategy # Strategies
# BSIDH
sibc-precompute-sdacs -p b2 -a bsidh # SDACs
sudo sibc -p b2 -f svelu -a bsidh -u bsidh-precompute-parameters # Tuned velusqrt parameters: the option -u is required
sudo sibc -p b2 -f svelu -a bsidh -u bsidh-precompute-strategy # Strategies
```

Furthermore, you can create tests by running `bash misc/create-tests.sh` and `bash misc/test-cli.sh`
(only csidh is handled by now).

## Remarks

The primes labeled as `b2`, `b3`, `b5`, and `b6` correspond with the examples 2, 3, 5, and 6 from
[B-SIDH paper](https://eprint.iacr.org/2019/1145), respectively. In particular,  `s1` denotes the
prime number given in [velusqrt paper](https://eprint.iacr.org/2020/341). The field airthmetic is
centered on primes `p = 3 mod 4`. Multiplying and squaring in `GF(p²) = GF(p)[u]/(u^2 + 1)` have a
cost of 3 and 2 multiplications in `GF(p)`, respectively.

The current implementation does not have implemented the B-SIDH key validation, it will be included
in the next library version.

## Changes

Significant changes are listed in the [CHANGELOG](CHANGELOG.md) file. Future integrations/modifications
are listed in the [TODOLIST](TODOLIST.md) file.

## Authors

1. **Gora Adj** <gora.adj@gmail.com>, <gora.adj@udl.cat>;
2. **_Jesús-Javier Chi-Domínguez_** <jesus.chidominguez@tuni.fi>, <chidoys@gmail.com>, <jjchi@computacion.cs.cinvestav.mx>; and
3. **_Francisco Rodríguez-Henríquez_** <francisco@cs.cinvestav.mx>.

### Main contributors

1. Jacob Appelbaum <j.appelbaum@tue.nl>; and
2. Leif Ryge <leif@synthesize.us>.

All contributors are listed in the [CONTRIBUTORS](CONTRIBUTORS) file.

## License

This project is licensed under the GNU general public license - see the
[LICENSE](LICENSE) file for details.

## Funding

This project has initially received funding from the European Research Council (ERC) under the
European Union's Horizon 2020 research and innovation programme (grant agreement No 804476).

The third author received partial funds from the Mexican Science council CONACyT project 313572.
