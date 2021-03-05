from progress.bar import Bar
from functools import reduce
import sys
import random
from math import ceil, floor, log, sqrt, pi
import numpy
import statistics
import argparse

base_path = "/usr/share/python3-sidh/data/"
strategy_data = base_path + "/strategies/"
sdacs_data = base_path + "/sdacs/"
sop_data = base_path + "/sop/"
ijk_data = base_path + "/ijk/"
gen_data = base_path + "/gen/"
tmp_dir = "./"  # Drop the files in the current working directory
# --------------------------------------------------------------------------------------------------------------------------------
def getinputs(args=sys.argv[1:]):

    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument(
        "-p",
        "--prime",
        help="prime number configuration should be stored in pSUFFIX (sop folder is taken as default).",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--formulaes",
        help="traditional (tvelu), sqrt (svelu), or hybrid (hvelu) velu formulaes to be used.",
        required=True,
    )
    parser.add_argument(
        "-a", "--algorithm", help="bsidh or csidh algorithm", required=True
    )
    parser.add_argument(
        "-s",
        "--style",
        help="style to be used: wd1 (with dummy operations and a single torsion point), wd2 (with dummy operations and a two torsion point), or df (dummy-free approach).",
    )
    parser.add_argument(
        "-b",
        "--benchmark",
        type=int,
        help="number of experiments to be used in the benchmark.",
        default=128,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest='verbose',
        action='store_true',
        help="Verbose mode.",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    options = parser.parse_args(args)
    return options


# Inputs
setting = getinputs(sys.argv[1:])
