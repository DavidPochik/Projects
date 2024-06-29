import numpy as np
import math
import sys
sys.path.insert(0, '/home/pochik.1/athena-accretion/mhd-winds-new-EOS/vis/python/')
import athena_read

FILE_name_f = "accretion.prim."+sys.argv[1]+".athdf"
FILE_name_g = "accretion.uov."+sys.argv[1]+".athdf"
R_PNS = sys.argv[2]
f_data = FILE_name_f
g_data = FILE_name_g
filename = 'vr'+str(R_PNS)+'km.txt'

data_prim = athena_read.athdf(f_data)
data_uov  = athena_read.athdf(g_data)

r1   = data_prim['x1f']
vr1  = (data_prim['vel1'])[0]

r    = r1[:]
vr   = vr1[0,:]
vrIC = vr[0]
rIC  = r[0]

print('vr at R_PNS = '+str(R_PNS)+' km is ' + str(vrIC) +' cm/s')
print('r at R_PNS = '+str(R_PNS)+' km is ' + str(rIC) + ' cm')

c=[vrIC]
with open(filename,"w") as file:
	for x in zip(c):
		file.write(str(*x))
