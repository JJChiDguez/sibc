from sidh.framework import *

if setting.algorithm == 'csidh':
    from sidh.csidh.ijk import *

elif setting.algorithm == 'bsidh':
    print("Not implemented yet!")
    exit(0)
    from sidh.bsidh.ijk import *

else:
    # This case shouldn't happen (only bsidh or csidh are allowed as inputs)
    exit(-1)


def main():
    pass
