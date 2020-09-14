from framework import *

if setting.algorithm == 'csidh':
    from csidh.test import *

elif setting.algorithm == 'bsidh':
    from bsidh.test import *

else:
    # This case shouldn't happen (only bsidh or csidh are allowed as inputs)
    exit(-1)

