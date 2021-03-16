# On the Velu's formulae and its applications to CSIDH and BSIDH constant-time implementations

At a combined computational expense of about *6l* field operations, Velu's
formulae are used to construct and evaluate degree-*l* isogenies in the vast
majority of isogeny-based primitive implementations. Recently, Bernstein, de
Feo, Leroux and Smith introduced an new approach for solving this same problem
at a reduced cost of just *O(sqrt(l))* field operations. In this work, we
present a concrete computational analysis of these novel formulae, along with
several algorithmic tricks that helped us to slightly, but noticeably, reduce
their practical cost.

## Installation 

Install the `sibc` module which provides the `sibc` program:
```
sudo python3 setup.py install
```

For development:
```
sudo pip install -e . 
```

### Debian package build

To build a package for Debian or Ubuntu, we suggest the use of stdeb:

```
sudo apt install -y dh-python python3-click  python3-sympy  python3-progress\
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


$ sibc --help
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
  -p, --prime [b2|b3|b5|b6|p1024|p1792|p512|sv]
                                  [default: p512]
  -f, --formula [tvelu|svelu|hvelu]
                                  [default: hvelu]
  -a, --algorithm [csidh|bsidh]   [default: csidh]
  -s, --style [wd1|wd2|df]        [default: df]
  -e, --exponent [5|10]           [default: 10]
  -m, --multievaluation           [default: False]
  -c, --curvemodel [edwards|montgomery]
                                  [default: montgomery]
  -b, --benchmark INTEGER         [default: 128]
  -t, --tuned                     [default: False]
  -v, --verbose                   Not the kind of verbosity you might expect
                                  [default: False]

  --version                       Show the version and exit.
  --help                          Show this message and exit.

Commands:
  bsidh-main
  bsidh-parameters
  bsidh-strategy
  bsidh-test
  csidh-bench
  csidh-bounds
  csidh-dh               Derive shared secret key from CSIDH sk, CSIDH pk
  csidh-genkey           Generate random CSIDH secret key
  csidh-header
  csidh-ijk
  csidh-main
  csidh-parameters
  csidh-pubkey           Derive CSIDH public key from CSIDH secret key
  csidh-sdacs
  csidh-strategy
  csidh-suitable-bounds
  csidh-test
  print-strategy         draw graphs
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
# CSIDH format (here p = 2^c * l_1 * .... l_n - 1)
c l_1 l_2 ... l_n

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
sibc -p p1024 -f tvelu -a csidh -s df csidh-main
sibc -p p512 -f svelu -a csidh -s wd2 csidh-main
sibc -p p1792 -f hvelu -a csidh -s wd1 -v csidh-main

sibc -p p512 -f hvelu -a csidh -s wd2 -b 1024 -v csidh-bench 
sibc -p p512 -f hvelu -a csidh -s wd1 -b 1024 -v csidh-bench 
sibc -p p512 -f hvelu -a csidh -s df  -b 1024 -v csidh-bench 

sibc -p p1792 -f tvelu -a csidh csidh-test
sibc -p p1792 -f svelu -a csidh csidh-test
sibc -p p1792 -f hvelu -a csidh csidh-test

# BSIDH
sibc -p b2 -f tvelu -a bsidh bsidh-test
sibc -p b2 -f svelu -a bsidh bsidh-test
sibc -p b2 -f hvelu -a bsidh -t bsidh-test
```

Remark, our implementation allows us to plot each optimal strategy required:

```bash
# CSIDH
sibc -p p1024 -f tvelu -a csidh -s df     print-strategy
sibc -p p512 -f svelu -a csidh -s wd2     print-strategy
sibc -p p1792 -f hvelu -a csidh -s wd1 -v print-strategy

# BSIDH
sibc -a bsidh -p b2 -f tvelu bsidh-strategy
sibc -a bsidh -t -p b2 -f svelu bsidh-strategy
sibc -a bsidh -t -p b2 -f hvelu bsidh-strategy
```

Additionally, one can created files with extension `.h` that includes all the
required variables in a the sdacs, strategies, and velusqrt (at least for CSIDH
implementations).

```bash
# Suitable bounds with m = 5
sibc -a csidh -p p512 -s wd2 -f hvelu -b 5 csidh-bounds
# Strategies also with m = 5
sibc -a csidh -p p512 -s wd2 -f hvelu -v -b 5 csidh-header
# SDACs (the flag -s doesn't affects the output)
sibc -a csidh -p p512 -s wd2 -f hvelu -v csidh-sdacs
# Optimal sizes of I, J, and K required in velusqrt (the flag -s doesn't affects the output)
sibc -a csidh -p p512 -s wd2 -f hvelu -v csidh-ijk
```

## BSIDH primes
Currently only `b2` is implemented and tested in the current API. Extending
this to other primes is straight-forward.

## Remarks

The primes labeled as `b2`, `b3`, `b5`, and `b6` correspond with the examples
2, 3, 5, and 6 of , respectively. In particular, we focused on primes such that
`p = 3 mod 4`. Additionally, the product and squaring in `F_p[i]/(i^2 + 1)`
were implemented using 3 and 2 products in `F_p`, respectively.

## Changes

Significant changes are listed in the [CHANGELOG](CHANGELOG.md) file.

## Authors

1. **Gora Adj** <gora.adj@gmail.com,gora.adj@udl.cat>,
2. **_Jesús-Javier Chi-Domínguez_** <jesus.chidominguez@tuni.fi>, <chidoys@gmail.com>, <jjchi@computacion.cs.cinvestav.mx>, and
3. **_Francisco Rodríguez-Henríquez_** <francisco@cs.cinvestav.mx>.

Additional contributors are listed in the [CONTRIBUTORS](CONTRIBUTORS) file.

## License

This project is licensed under the GNU general public license - see the
[LICENSE](LICENSE) file for details.

## Funding

This project has received funding from the European Research Council (ERC)
under the European Union's Horizon 2020 research and innovation programme
(grant agreement No 804476). 

The third author received partial funds from the Mexican Science council
CONACyT project 313572.
