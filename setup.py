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
    name="sidh",
    version=__version__,
    author="JJChiDguez",
    author_email="chidoys@gmail.com",
    description=("Supersingular Isogeny Diffie Hellman (SIDH) constructions: csidh and bsidh using velusqrt formulae"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="GPLv3",
    url="https://github.com/JJChiDguez/velusqrt",
    packages=setuptools.find_packages(),
    keywords="csidh, bsidh, sidh, encryption",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    entry_points={
        "console_scripts": [
            "sidh=sidh.__main__:main",
            "sidh-bench=sidh.csidh.bench:main",
            "sidh-test=sidh.csidh.test:main",
            "sidh-print-strategy=sidh.printstrategy:main",
            "sidh-parameters=sidh.csidh.parameters:main",
            "sidh-bounds=sidh.csidh.bounds:main",
            "sidh-header=sidh.csidh.header:main",
            "sidh-sdacs=sidh.csidh.sdacs:main",
            "sidh-ijk=sidh.csidh.ijk:main",
            "sidh-csidh-util=sidh.csidh.util:main",
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
