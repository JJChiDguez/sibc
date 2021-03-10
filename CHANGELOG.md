- package as sidh with submodules csidh and bsidh
-- setup.py includes the precomputed data files for most supported parameters
- refactor project to use click module for argument handling
- general refactor to object oriented interface
- programs such as `sidh` are automatically generated at install time from
  setup.py entrypoints.
- `sidh` now has the following subcommands:
  - bench, csidh-bounds, csidh-header, csidh-ijk, csidh-main,
    csidh-parameters, csidh-sdacs, csidh-strategy, csidh-suitable-bounds,
    csidh-test, genkey
- pytest now supports dynamic running of tests:
  - the full test suite takes ~175 minutes on an i7-9750H system
  - csidh supports p512, p1024, p1792 - with each supported style and formula
