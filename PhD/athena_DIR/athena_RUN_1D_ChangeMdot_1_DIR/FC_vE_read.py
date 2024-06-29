import numpy as np
import math
import sys

filename = 'FC_vE_FLAG.txt'
FC_vE    = float(sys.argv[1]) # cm/s
FC_vE_FLAG = int(-1)
if(FC_vE > 0.0):
	FC_vE_FLAG = int(1)
	print('(FC_ve_read.py) Vr(Final active zone) > 0.0')
	print('(FC_ve_read.py) Vr = ' +str(FC_vE)+' cm/s')

c=[FC_vE_FLAG]
with open(filename,"w") as file:
	for x in zip(c):
		file.write(str(*x))
