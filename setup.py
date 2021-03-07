#!/usr/bin/env python3
import setuptools
import os

__version__ = '0.0.1'

if os.path.exists('requirements.txt'):
    with open("requirements.txt", "r") as obj:
        requirements = obj.read().splitlines()
else:
    requirements = []

with open("README.md", "r") as obj:
    long_description = obj.read()

strategy_files = ["data/strategies/" + strategy for strategy in os.listdir('data/strategies/')]
gen_files = ["data/gen/" + gen for gen in os.listdir('data/gen/')]
ijk_files = ["data/ijk/" + ijk for ijk in os.listdir('data/ijk/')]
sdacs_files = ["data/sdacs/" + sdacs for sdacs in os.listdir('data/sdacs/')]
sop_files = ["data/sop/" + sop for sop in os.listdir('data/sop/')]
figure_files = ["data/figures/" + figure for figure in os.listdir('data/figures/')]
data_files = [
        ( "/usr/share/python3-sidh/data/gen/", gen_files,),
        ( "/usr/share/python3-sidh/data/ijk/", ijk_files,),
        ( "/usr/share/python3-sidh/data/sdacs/", sdacs_files),
        ( "/usr/share/python3-sidh/data/sop/", sop_files),
        ( "/usr/share/python3-sidh/data/strategies/", strategy_files),
        ( "/usr/share/python3-sidh/data/figure/", figure_files,),
        ( "/usr/bin/", ["misc/sidh-autosearch"]),
    ]

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
            "sidh=sidh.main:main",
            "sidh-bench=sidh.bench:main",
            "sidh-test=sidh.test:main",
            "sidh-print-strategy=sidh.printstrategy:main",
            "sidh-parameters=sidh.parameters:main",
            "sidh-bounds=sidh.bounds:main",
            "sidh-header=sidh.header:main",
            "sidh-sdacs=sidh.sdacs:main",
            "sidh-ijk=sidh.ijk:main",
            "sidh-csidh-util=sidh.csidh.util:main",
        ]
    },
    zip_safe=False,
    install_requires=requirements,
    include_package_data=True,
    data_files=data_files,
)
