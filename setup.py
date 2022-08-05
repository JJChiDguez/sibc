#!/usr/bin/env python3
import setuptools
import os

try:
    from stdeb.command.sdist_dsc import sdist_dsc
    from stdeb.command.bdist_deb import bdist_deb
except ImportError:
    sdist_dsc = None
    bdist_deb = None


__version__ = '1.0.4'

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
    description=("Supersingular Isogeny-Based Cryptography constructions: currently, csidh, bsidh, sidh, and sike are implemented by using traditional and velusqrt formulae on Montgomery curve x-only projective coordinates"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="GPLv3",
    url="https://github.com/JJChiDguez/sibc",
    packages=setuptools.find_packages(),
    keywords="csidh, bsidh, sidh, sike, sibc, encryption",
    classifiers=[
        "Programming Language :: Python :: 3",
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    entry_points={
        "console_scripts": [
            "sibc=sibc.__main__:main",
            "sibc-precompute-sdacs=sibc.precompute_sdacs:main", 
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
            "" : [ "sibc/data/*" ]
        }
)
