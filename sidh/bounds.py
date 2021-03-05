from sidh.framework import *

if setting.algorithm == 'csidh':
    from sidh.csidh.bounds import *

else:
    print("Only csidh requires this script")
    exit(-1)
