#import h5py
#import matplotlib.pyplot as plt
import numpy as np
import math
import sys
pythonpath=str(sys.argv[5])
sys.path.insert(0, pythonpath)
import athena_read

# inputs
rPNS         = float(sys.argv[1]) # in km
outputnumber = str(sys.argv[2])
Lnu          = str(sys.argv[3])
Mdot         = float(sys.argv[4])

FILE_name_f = "accretion.prim."+outputnumber+".athdf"
FILE_name_g = "accretion.uov."+outputnumber+".athdf"
IC_name     = 'Lnu_'+str(Lnu)+'e51_Mdot_'+str(Mdot)+'_RPNS_'+str(int(rPNS))+'km_extrapolate.txt'
f_data = FILE_name_f
g_data = FILE_name_g

data_prim = athena_read.athdf(f_data)
data_uov  = athena_read.athdf(g_data)

r1        = data_prim['x1v']
vr1       = (data_prim['vel1'])[0]
P1        = (data_prim['press'])[0]
rho1      = (data_prim['rho'])[0]
T         = (data_uov['dt1'])[0]
elecfrac1 = (data_prim['r0'])[0]

rho  = rho1[0,:]
vr1  = vr1[0,:]
T    = T[0,:]
r    = r1[:]
Ye   = elecfrac1[0,:]
nR   = np.size(r)
for i in range(0,nR):
	if(r[i]/1.0e5>float(rPNS)):
		iNew = i
		print('New inner boundary r = ' + str(r[i]/1.0e5) + ' km')
		print('iNew                 = ' + str(iNew))
		break

r_interp   = np.linspace(r[int(iNew)],r[nR-1],int(nR)).T
rho_interp = np.interp(r_interp,r,rho)
vr_interp  = np.interp(r_interp,r,vr1)
T_interp   = np.interp(r_interp,r,T)
Ye_interp  = np.interp(r_interp,r,Ye)

c=[r_interp, rho_interp, vr_interp, T_interp, Ye_interp]
with open(IC_name, "w") as file:
    for x in zip(*c):
        file.write("{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\n".format(*x))

