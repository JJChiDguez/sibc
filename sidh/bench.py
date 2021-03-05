from sidh.framework import *

if setting.algorithm == 'csidh':
    from sidh.csidh.bench import *

elif setting.algorithm == 'bsidh':

    print("// Not required, BSIDH has a fixed constant-time implementation")
    exit(0)

else:
    # This case shouldn't happen (only bsidh or csidh are allowed as inputs)
    exit(-1)


def main():
    pass
