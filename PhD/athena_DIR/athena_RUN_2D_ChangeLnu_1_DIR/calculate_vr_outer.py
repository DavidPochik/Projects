#import h5py
#import matplotlib.pyplot as plt
import numpy as np
import math
import sys
pythonpath=str(sys.argv[3])
sys.path.insert(0, pythonpath)
import athena_read

fileDIR     = str(sys.argv[1])
file_name   = 'vr_FLAG.txt'
FILE_name_f = fileDIR+"/accretion.prim."+sys.argv[2]+".athdf"
f_data      = FILE_name_f
data_prim   = athena_read.athdf(f_data)
r           = data_prim['x1v']
theta       = data_prim['x2v']
vr          = (data_prim['vel1'])[0]
nR          = np.size(r)
nTheta      = np.size(theta)

vr_average = 0.0
for i in range(0,nTheta-1):
	delta_Theta = theta[i+1] - theta[i]
	vr_average += 0.5 * 0.25 * (np.sin(theta[i]) * vr[i,nR-1] + np.sin(theta[i+1]) * vr[i+1,nR-1])

vr_FLAG = int(-1)
if(vr_average>0.0):
	vr_FLAG = 1
	print('(calculate_vr_outer.py) <Vr(Final active zone)> > 0.0')
	print('(calculate_vr_outer.py) <Vr>                    = ' +str(vr_average)+' cm/s')

c=[vr_FLAG]
with open(file_name, "w") as file:
    for x in zip(c):
        file.write(str(*x)) 


