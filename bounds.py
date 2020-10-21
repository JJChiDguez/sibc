from framework import *

if setting.algorithm == 'csidh':
    from csidh.bounds import *

else:
    print("Only csidh requires this script")
    exit(-1)
