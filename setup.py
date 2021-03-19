#!/usr/bin/env python3
import setuptools
import os

try:
    from stdeb.command.sdist_dsc import sdist_dsc
    from stdeb.command.bdist_deb import bdist_deb
except ImportError:
    sdist_dsc = None
    bdist_deb = None


__version__ = '0.0.1'

if os.path.exists('requirements.txt'):
    with open("requirements.txt", "r") as obj:
        requirements = obj.read().splitlines()
else:
    requirements = []

with open("README.md", "r") as obj:
    long_description = obj.read()

setuptools.setup(
    name="sibc",
    version=__version__,
    author="JJChiDguez",
    author_email="chidoys@gmail.com",
    description=("Supersingular Isogeny-Based Cryptography constructions: currently, csidh and bsidh are implemented by using traditional and velusqrt formulae on Montgomery curve x-only projective coordinates"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="GPLv3",
    url="https://github.com/JJChiDguez/velusqrt",
    packages=setuptools.find_packages(),
    keywords="csidh, bsidh, sibc, encryption",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    entry_points={
        "console_scripts": [
            "sibc=sibc.__main__:main",
            "sibc-bench=sibc.csidh.bench:main",
            "sibc-test=sibc.csidh.test:main",
            "sibc-plot-strategy=sibc.plot_strategy:main",
            "sibc-precompute-parameters=sibc.csidh.parameters:main",
            "sibc-bounds=sibc.csidh.bounds:main",
            "sibc-header=sibc.csidh.header:main",
            "sibc-sdacs=sibc.csidh.sdacs:main",
            "sibc-ijk=sibc.csidh.ijk:main",
            "sibc-csidh-util=sibc.csidh.util:main",
        ]
    },
    zip_safe=False,
    install_requires=requirements,
    include_package_data=True,
    cmdclass=dict(
        bdist_deb=bdist_deb,
        sdist_dsc=sdist_dsc,
    ),
    package_data = {
            "" : [ "data/*/*" ]
        }
)
