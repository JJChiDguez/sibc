# On the Velu's formulae and its applications to CSIDH and B-SIDH constant-time implementations


At a combined computational expense of about *6l* field operations, Velu's formulae are used to construct and evaluate degree-*l* isogenies in the vast majority of isogeny-based primitive implementations. Recently, Bernstein, de Feo, Leroux and Smith introduced an new approach for solving this same problem at a reduced cost of just *O(sqrt(l))* field operations. In this work, we present a concrete computational analysis of these novel formulae, along with several algorithmic tricks that helped us to slightly, but noticeably, reduce their practical cost.


## Installation 

Install the `sidh` module which provides the `sidh` program:
```
sudo python3 setup.py install
```

For development:
```
sudo pip install -e . 
```

## Usage 
The syntax compilation can be viewed by running one of the following three commands:
```bash
# Corresponding with the key-exchange protocol
sidh --help
# Corresponding with benchmarking (only for CSIDH, which has a variable running-time cost independent from the key)
sidh csidh-bench
# Corresponding with the costs of KPs (Kernel Point computation), xISOG (isogeny construction), and xEVAL (isogeny evaluation)
sidh csidh-test
```

Usage for the `sidh` tool:
```bash
$ sidh --help
Usage: sidh [OPTIONS] COMMAND [ARGS]...

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
  -e, --exponent [2]              [default: 2]
  -c, --curvemodel [edwards|montgomery]
                                  [default: montgomery]
  -b, --benchmark INTEGER         [default: 128]
  -v, --verbose                   Not the kind of verbosity you might expect
                                  [default: False]

  --version                       Show the version and exit.
  --help                          Show this message and exit.

Commands:
  csidh-bench
  csidh-bounds
  csidh-header
  csidh-ijk
  csidh-main
  csidh-parameters
  csidh-sdacs
  csidh-strategy
  csidh-suitable-bounds
  csidh-test
  genkey                 Generate a secret key
  print-strategy         draw graphs
```

## sidh api

Currently only CSIDH is available as a library function.

### Basic shared secret generation
```python3
from sidh.csidh import CSIDH, default_parameters
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

The field characteristic `p` should be stored in directory `data/sop/`, and CSIDH and BSIDH have different structures (see below):

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

For the case of BSIDH, `M := (4^c * l_1^{e_1} * l_2^{e_2} * ... * l_n^{e_n})` must divide `(p + 1)`, and `N := (l'_1^{e'_1} * l'_2^{e'_2} * ... * l'_n^{e'_n})` must divide `(p-1)`. Additionally, the order-`M` generators `PA`, `QA` and `PQA := PA - QA` should be stored in directory `gen/` as projective x-coordinate points. Similarly, the order-`N` generators `PB`, `QB` and `PQB := PB - QB` also should be stored it the same directory. Both 3-tuples of points must be stored in a single file with the following syntax:

```bash
Re(x(PA)) Im(x(PA)) Re(x(QA)) Im(x(QA)) Re(x(PQA)) Im(x(PQA))
Re(x(PB)) Im(x(PB)) Re(x(QB)) Im(x(QB)) Re(x(PQB)) Im(x(PQB))
```

where `Re(X)` and `Im(X)` denote the real and imaginary parts of X with respect to `F_p[i]/(i^2 + 1)`, respectively. Moreover, all the above twelve integers should be stored in hexadecimal."

## Examples

We summarize some examples of runs as follows:

```bash
# CSIDH
sidh -p p1024 -f tvelu -a csidh -s df csidh-main
sidh -p p512 -f svelu -a csidh -s wd2 csidh-main
sidh -p p1792 -f hvelu -a csidh -s wd1 -v csidh-main

sidh -p p512 -f hvelu -a csidh -s wd2 -b 1024 -v csidh-bench 
sidh -p p512 -f hvelu -a csidh -s wd1 -b 1024 -v csidh-bench 
sidh -p p512 -f hvelu -a csidh -s df  -b 1024 -v csidh-bench 

sidh -p p1792 -f tvelu -a csidh csidh-test
sidh -p p1792 -f svelu -a csidh csidh-test
sidh -p p1792 -f hvelu -a csidh csidh-test

# BSIDH (all of these are currently not ported over to the click tui)
sidh -p b6 -f tvelu -a bsidh
sidh -p b5 -f svelu -a bsidh  <!- data missing for this ( /usr/share/python3-sidh/data/strategies/bsidh-b5-svelu-classical )
sidh -p b2 -f hvelu -a bsidh -v

sidh-test -p b6 -f tvelu -a bsidh
sidh-test -p b6 -f svelu -a bsidh
sidh-test -p b6 -f hvelu -a bsidh

```

Remark, our implementation allows us to plot each optimal strategy required:

```bash
# CSIDH
sidh -p p1024 -f tvelu -a csidh -s df     print-strategy
sidh -p p512 -f svelu -a csidh -s wd2     print-strategy
sidh -p p1792 -f hvelu -a csidh -s wd1 -v print-strategy

# BSIDH
sidh-print-strategy -p b6 -f tvelu -a bsidh
sidh-print-strategy -p b5 -f svelu -a bsidh <! data missing for this ( /usr/share/python3-sidh/data/strategies/bsidh-b5-svelu-classical )
sidh-print-strategy -p b2 -f hvelu -a bsidh -v
```

Additionally, one can created files with extension `.h` that includes all the required variables in a the sdacs, strategies, and velusqrt (at least for CSIDH implementations).

```bash
# Suitable bounds with m = 5
sidh -a csidh -p p512 -s wd2 -f hvelu -b 5 csidh-bounds
# Strategies also with m = 5
sidh -a csidh -p p512 -s wd2 -f hvelu -v -b 5 csidh-header
# SDACs (the flag -s doesn't affects the output)
sidh -a csidh -p p512 -s wd2 -f hvelu -v csidh-sdacs
# Optimal sizes of I, J, and K required in velusqrt (the flag -s doesn't affects the output)
sidh -a csidh -p p512 -s wd2 -f hvelu -v csidh-ijk
```

## Remarks

The primes labeled as `b2`, `b3`, `b5`, and `b6` correspond with the examples 2, 3, 5, and 6 of , respectively. In particular, we focused on primes such that `p = 3 mod 4`. Additionally, the product and squaring in `F_p[i]/(i^2 + 1)` were implemented using 3 and 2 products in `F_p`, respectively.

## Changes

Significant changes are listed in the [CHANGELOG](CHANGELOG.md) file.

## Authors

1. **Gora Adj** <gora.adj@gmail.com,gora.adj@udl.cat>,
2. **_Jesús-Javier Chi-Domínguez_** <jesus.chidominguez@tuni.fi>, <chidoys@gmail.com>, <jjchi@computacion.cs.cinvestav.mx>, and
3. **_Francisco Rodríguez-Henríquez_** <francisco@cs.cinvestav.mx>.

Additional contributors are listed in the [CONTRIBUTORS](CONTRIBUTORS) file.

## License

This project is licensed under the GNU general public license - see the [LICENSE](LICENSE) file for details

## Funding

This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 804476). 

The third author received partial funds from the Mexican Science council CONACyT project 313572.
