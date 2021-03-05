from sidh.framework import *

if setting.algorithm == 'csidh':
    from sidh.csidh.parameters import *

elif setting.algorithm == 'bsidh':
    from sidh.bsidh.parameters import *

else:
    # This case shouldn't happen (only bsidh or csidh are allowed as inputs)
    exit(-1)

def main():
    pass
