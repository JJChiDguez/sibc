from framework import *

if setting.algorithm == 'csidh':
    from csidh.parameters import *

elif setting.algorithm == 'bsidh':
    from bsidh.parameters import *

else:
    # This case shouldn't happen (only bsidh or csidh are allowed as inputs)
    exit(-1)

