
<a name="v1.0.4"></a>
## [v1.0.4](https://github.com/JJChiDguez/sibc/compare/v1.0.3...v1.0.4)

> 2022-08-05

### Docs

* update README.md

### Pull request

* optimize csidh 

### Fix

* correct bounds concerning CSIDH-512 for wd1, wd2, and df; files modified               sibc/data/exponents/mitm/128/df/e10               sibc/data/exponents/mitm/128/wd1/e10               sibc/data/exponents/mitm/128/wd2/e5


<a name="v1.0.3"></a>
## [v1.0.3](https://github.com/JJChiDguez/sibc/compare/v1.0.1...v1.0.3)

> 2021-11-09

### Docs

* small tweaks in the README.md and setup.py files
* fix citation (now using and instead of commas)
* updated setup.py description
* updated CHANGELOG.md and TODOLIST.md

### Fix

* update bsidh-test to last public parameter representation
* strategy plotting is working on macOS


<a name="v1.0.1"></a>
## [v1.0.1](https://github.com/JJChiDguez/sibc/compare/v0.0.1...v1.0.1)

> 2021-08-02

### Docs

* updating README.md file
* adding missing bracket in how to cite it
* updating the README.md file Including how to cite the library (bibtex format), and fixing/updating authors email
* updating the readme.md
* extending and updating the README.md file (included logo, and a better sibc description)
* updating url in setup.py
* updating TODOLIST.md, and including commented example using keygen() procedure for csidh

### Feat

* added cli for sike and bsike, updated README
* added sidh and sike procedures updated README.md file, and extended tests; SIDH and SIKE tests are simple. removed duplicated code in csidh and bsidh, local tests correctly passed =D
* including sidh primes, sidh is not currently working (missing basis generation)      Wait for the coming updates (sidh will be soonly avalaible)
* renamed bsidh primes according to https://eprint.iacr.org/2020/1109, improving the integration of new [csidh/bsidh]-prime instances (by moving into sibc/parameters.py). In fact, a new instances should be included in [csidh/bsidh]_primes list/dictionary.
* including degree-3 isogeny procedures Now, we can play with kps_3, xisog_3, and xeval_3; Small modification in xtpl, which is now taking as input (A + 2C : A - 2C)
* including degree-2 isogeny procedures Now, we can play with kps_2, xisog_2, and xeval_2
* small modification in the file description of csidh-primes (assuming first integer is the cofactor instead of the exponent of 2 such that p = cofactor * prod(l_i's) - 1)
* BSIDH class is almost complete (still missing complete validation). However, bsidh public key is now corresponding with the affine x-coordinates of the images point.
* including recovery Montgomery curve coefficient BREAKING CHANGE: now xdbladd and thus Ladder3pt assumes the input curve coefficient is in affine coords

### Fix

* now, [bsidh,csidh]/precompute_parameters.py is correctly working
* correct license name in classifiers
* bsidh/test.py is now working with new xdbladd() and Ladder3pt() version (small issue solved)

### Recfactor

* integrating Ladder3pt() and xdbladd() taking the affine parameter a24 = (A + 2C)/(4C)

### Refact

* bsidh and csidh have their own subdirectories for the data

### Refactor

* now csidh allows any kind of prime name. FOR EXAMPLE, including a new prime named as p1024k74 is possible
* bsidh handles FIRST strategy evaluation by using get_A() instead of xisog(). Best perfomance by using projective get_A()


<a name="v0.0.1"></a>
## v0.0.1

> 2021-03-23

### Constants

* use pkg_resources

### Feat

* updating CHANGELOG.md; at this point we have the first complete working version of the code
* moving data from {ijk,gen,strategy}/montgomery to {ijk,gen,strategy}/montgomery, for an easy future integration of different curve models
* updating README.md and TODOLIST files
* updating path file required in csidh/header.py
* adding small help descriptions for the missing terminal commands
* removing outdated console commands, and including sibc-precompute-sdacs as a new command
* new isolated file for precomputing sdacs

### Fix

* data generation according curvemodel=montgomery
* small modification in the bsidh/__init__.py description
* giving the right ordering for the comments in __init__.py and test_xxxx.py files
* plotting of strategies working for both of csidh and bsidh
* integrating the new library name (sibc)
* README.md and CONTRIBUTORS according to the new library name
* .gitignore according to the new library name
* create-tests.sh and test-cli.sh; sidh-autosearch has been deleted

### Reactor

* cleaning the code (removing not used variables from sibc/constants.py file)

### Refactor

* moving strategy_block_cost() in the __init__() of Gae_xx().
* exponent bounds for csidh are now read them from disk (cleaning the code)
* updating/adding precomputed data for both of csidh and csidh primitives
* adating the tests according current file description, and adding scripts for precomputing data of csidh and bsidh
* insolating strategy precomputation
* IsogenyCurve class reads tuned velusqrt parameters from files with suffix either -scaled or -unscaled
* MontgomeryCurve class assumes sdacs were precomputed before

### Update

* adding sibc-0.0.1.tar.gz deb_dist/ to .gitignore
* adding sibc-installed.txt to .gitignore
* TODOLIST (anyother improvements [?])

