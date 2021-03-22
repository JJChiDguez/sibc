#  (2021-03-22)


### Bug Fixes

* .gitignore according to the new library name ([57d74f9](https://github.com/JJChiDguez/velusqrt/commit/57d74f9e573ac7693e2e5a58d97aca26224ce10c))
* create-tests.sh and test-cli.sh; sidh-autosearch has been deleted ([431c5bd](https://github.com/JJChiDguez/velusqrt/commit/431c5bdf49744a930141b6daf984d68b2199eaf8))
* data generation according curvemodel=montgomery ([c53e24d](https://github.com/JJChiDguez/velusqrt/commit/c53e24daa4b0fcbf281b783695741ead13e13a19))
* giving the right ordering for the comments in __init__.py and test_xxxx.py files ([ff9060a](https://github.com/JJChiDguez/velusqrt/commit/ff9060a146dddad05176bea6d8f0051dec8d853f))
* integrating the new library name (sibc) ([504369c](https://github.com/JJChiDguez/velusqrt/commit/504369c0b7969adb222f0866f752f69548accff1))
* plotting of strategies working for both of csidh and bsidh ([bb77b6e](https://github.com/JJChiDguez/velusqrt/commit/bb77b6ed85906a2ae095056b510d76e5a3a4652c))
* README.md and CONTRIBUTORS according to the new library name ([c638d73](https://github.com/JJChiDguez/velusqrt/commit/c638d730d9978845944a55c44d94c7b3f3e725c3))
* small modification in the bsidh/__init__.py description ([751d1b0](https://github.com/JJChiDguez/velusqrt/commit/751d1b0f4e3abb3a1efba96f466c2abf22388665))


### Features

* adding small help descriptions for the missing terminal commands ([fabd8f9](https://github.com/JJChiDguez/velusqrt/commit/fabd8f9b9ca5ddb3cdd7dcc1589114c67e660554))
* moving data from {ijk,gen,strategy}/montgomery to {ijk,gen,strategy}/montgomery, for an easy future integration of different curve models ([aba25d4](https://github.com/JJChiDguez/velusqrt/commit/aba25d4e31c60026283418d4ae33f83df962d637))
* new isolated file for precomputing sdacs ([c5065b6](https://github.com/JJChiDguez/velusqrt/commit/c5065b6287839d371e9620446f2d368cebe8c9c1))
* removing outdated console commands, and including sibc-precompute-sdacs as a new command ([1df88cb](https://github.com/JJChiDguez/velusqrt/commit/1df88cbe706e8923da8950394ca72ca38e3a8ec5))
* updating path file required in csidh/header.py ([a5f060c](https://github.com/JJChiDguez/velusqrt/commit/a5f060c8866a54729e671a96133e70fdb8c2aae5))
* updating README.md and TODOLIST files ([968aefa](https://github.com/JJChiDguez/velusqrt/commit/968aefab37d743e979db5f754a59c8900dedd4bb))

### Refactors

* package as ~~sidh~~ (now named as `sibc`) with submodules csidh and bsidh
  * setup.py includes the precomputed data files for most supported parameters
* refactor project to use click module for argument handling
* general refactor to object oriented interface
* programs such as ~~sidh~~ (now named as `sibc`) are automatically generated at install time from
  setup.py entrypoints.
* ~~sidh~~ (now named as `sibc`) now has the following subcommands:
  * bench, csidh-bounds, csidh-header, csidh-ijk, csidh-main,
    ~~csidh-parameters~~, csidh-sdacs, ~~csidh-strategy~~, ~~csidh-suitable-bounds~~,
    csidh-test, genkey
* pytest now supports dynamic running of tests:
  * the full test suite takes ~175 minutes on an i7-9750H system
  * csidh supports p512, p1024, p1792 - with each supported style and formula
* gae library does not use coeff on gae dh method output anymore
  * dh now returns projective coordinates
* CSIDH library api exposes CSIDH.secret_key(), CSIDH.public_key(), CSIDH.dh()
  * these methods that consume and return bytes objects
* ~~sidh~~ (now named as `sibc`) cli tools: pubkey, genkey, dh
  * these tools produce and consume base64 encoded byte values
* CSIDH object and cli tool now have `--tuned` and `--multievaluation` options
