from framework import *

if setting.algorithm == 'csidh':
    from csidh.main import *

elif setting.algorithm == 'bsidh':
    from bsidh.main import *

else:
    # This case shouldn't happen (only bsidh or csidh are allowed as inputs)
    exit(-1)
