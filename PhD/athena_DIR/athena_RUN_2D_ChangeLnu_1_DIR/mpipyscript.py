from mpi4py import MPI
import numpy as np
import math
import h5py
import sys
import os
#sys.path.insert(0, '/users/PAS2055/dpochik/test/athena_CURRENT_WORKING_DIRECTORY_DIR/athena_WORKING_DIR/athena_RUN_2D_DIR/vis/python')
pythonpath = str(sys.argv[10])
sys.path.insert(0,pythonpath)
import athena_read
#import athena_read
import cmath
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
from matplotlib import ticker, cm
import matplotlib as mpl
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from scipy.interpolate import interp2d
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
try:
    import helmeos
except ImportError:
    try:
        from .. import helmeos
    except ImportError:
        import sys
        sys.path.append(os.path.dirname(os.path.abspath(__file__)))
        import helmeos

fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "helm_table.dat")
by_filename = helmeos.HelmTable(fn=fn, temp_n=201, dens_n=541)

DIR_name    = str(sys.argv[1])
start_index = int(sys.argv[2])
max_index   = int(sys.argv[3])
nTheta      = int(sys.argv[4])
nR          = int(sys.argv[5])
ts          = float(sys.argv[6])
rmax        = float(sys.argv[7])
name        = str(sys.argv[8])
Mdot_set    = float(sys.argv[9])
Mach_set    = float(sys.argv[11])
SrMax       = float(sys.argv[12])
nT          = max_index - start_index + 1
t0          = 0.0
time        = np.linspace(t0, ts, nT).T

comm         = MPI.COMM_WORLD
rank         = comm.Get_rank()
size         = comm.Get_size()
istart       = start_index
nmax         = max_index 
num_per_rank = int(round((nmax-istart)/size))
lower_ind    = istart+rank*num_per_rank
upper_ind    = istart+(rank + 1)*num_per_rank
if (rank==size-1):
	upper_ind = nmax+1

rho    = np.zeros(shape=(upper_ind-lower_ind,nTheta,nR))
vr     = np.zeros(shape=(upper_ind-lower_ind,nTheta,nR))
vtheta = np.zeros(shape=(upper_ind-lower_ind,nTheta,nR))
S      = np.zeros(shape=(upper_ind-lower_ind,nTheta,nR))
Ye     = np.zeros(shape=(upper_ind-lower_ind,nTheta,nR))
T      = np.zeros(shape=(upper_ind-lower_ind,nTheta,nR))
P      = np.zeros(shape=(upper_ind-lower_ind,nTheta,nR))
cs     = np.zeros(shape=(upper_ind-lower_ind,nTheta,nR))
qdot   = np.zeros(shape=(upper_ind-lower_ind,nTheta,nR))
gamma  = np.zeros(shape=(upper_ind-lower_ind,nTheta,nR))

Na     = 6.02214e23      # Avogadro's constant
mb     = 1.0 / Na        # Baryon mass in g
kb     = 1.380649e-16    # Boltzmann constant in erg/K
abar   = 1.0
c      = 3e10            # cm/s
k      = 1.380649e-16    # erg/K
hbar   = 1.0546e-27      # erg s
kb_MeV = 8.61733326e-11 # MeV/K

comm.Barrier()
# Read in data
for i in range(lower_ind,upper_ind):
	if(i<10):
	        n=DIR_name+"/accretion.prim.0000"+str(i)+".athdf"
	        n1=DIR_name+"/accretion.uov.0000"+str(i)+".athdf"
	if(i>=10 and i<100):
	        n=DIR_name+"/accretion.prim.000"+str(i)+".athdf"
	        n1=DIR_name+"/accretion.uov.000"+str(i)+".athdf"
	if(i>=100 and i<1000):
	        n=DIR_name+"/accretion.prim.00"+str(i)+".athdf"
	        n1=DIR_name+"/accretion.uov.00"+str(i)+".athdf"
	if(i>=1000):
	        n=DIR_name+"/accretion.prim.0"+str(i)+".athdf"
	        n1=DIR_name+"/accretion.uov.0"+str(i)+".athdf"

	data1                  = athena_read.athdf(n)
	data2                  = athena_read.athdf(n1)
	rhod                   = (data1['rho'])[0]	
	vrd                    = (data1['vel1'])[0]
	vthetad                = (data1['vel2'])[0]
	yed                    = (data1['r0'])[0]
	tempd                  = (data2['dt1'])[0]
	pd                     = (data1['press'])[0]
	csd                    = (data2['dt3'])[0]
	qdotd                  = (data2['dt2'])[0]

	abar    = 1.0
	EOSData = by_filename.eos_DT(rhod, tempd, abar, yed)
	
	rho[i-lower_ind][:][:]    = rhod[0:nTheta][0:nR]                          # g/cm^3
	vr[i-lower_ind][:][:]     = vrd[0:nTheta][0:nR]                           # cm/s
	vtheta[i-lower_ind][:][:] = vthetad[0:nTheta][0:nR]
	S[i-lower_ind][:][:]      = (EOSData['stot'] * (mb / kb))[0:nTheta][0:nR] # kb/baryon
	Ye[i-lower_ind][:][:]     = yed[0:nTheta][0:nR]
	T[i-lower_ind][:][:]      = tempd[0:nTheta][0:nR]                         # K
	P[i-lower_ind][:][:]      = pd[0:nTheta][0:nR]                            # erg/cm^3
	cs[i-lower_ind][:][:]     = csd[0:nTheta][0:nR]                           # cm/s
	qdot[i-lower_ind][:][:]   = qdotd[0:nTheta][0:nR]                         # erg/s/g
	gamma[i-lower_ind][:][:]  = (EOSData['gam1'])[0:nTheta][0:nR] 
	
#collect from all ranks
iSize = upper_ind - lower_ind
rho_gather_T    = []
vr_gather_T     = []
vtheta_gather_T = []
S_gather_T      = []
Ye_gather_T     = []
T_gather_T      = []
P_gather_T      = []
cs_gather_T     = []
qdot_gather_T   = []
gamma_gather_T  = []
for i in range(0,iSize):
	comm.Barrier()
	rho_gather    = comm.gather(rho[i,0:nTheta,0:nR],root=0)
	comm.Barrier()
	vr_gather     = comm.gather(vr[i,0:nTheta,0:nR],root=0)
	comm.Barrier()
	vtheta_gather = comm.gather(vtheta[i,0:nTheta,0:nR],root=0)
	comm.Barrier()
	S_gather      = comm.gather(S[i,0:nTheta,0:nR],root=0)
	comm.Barrier()
	Ye_gather     = comm.gather(Ye[i,0:nTheta,0:nR],root=0)
	comm.Barrier()
	T_gather      = comm.gather(T[i,0:nTheta,0:nR],root=0)
	comm.Barrier()
	P_gather      = comm.gather(P[i,0:nTheta,0:nR],root=0)
	comm.Barrier()
	cs_gather     = comm.gather(cs[i,0:nTheta,0:nR],root=0)
	comm.Barrier()
	qdot_gather   = comm.gather(qdot[i,0:nTheta,0:nR],root=0)
	comm.Barrier()
	gamma_gather  = comm.gather(gamma[i,0:nTheta,0:nR],root=0)
	comm.Barrier()

	rho_gather_T.append(rho_gather)
	vr_gather_T.append(vr_gather)
	vtheta_gather_T.append(vtheta_gather)
	S_gather_T.append(S_gather)
	Ye_gather_T.append(Ye_gather)
	T_gather_T.append(T_gather)
	P_gather_T.append(P_gather)
	cs_gather_T.append(cs_gather)
	qdot_gather_T.append(qdot_gather)
	gamma_gather_T.append(gamma_gather)
	
comm.Barrier()
if(rank==0):
	# reshape data into more 'intuitive' arrays
	rho_gather_T    = np.array(rho_gather_T)
	vr_gather_T     = np.array(vr_gather_T)
	vtheta_gather_T = np.array(vtheta_gather_T)
	S_gather_T      = np.array(S_gather_T)
	Ye_gather_T     = np.array(Ye_gather_T)
	T_gather_T      = np.array(T_gather_T)
	P_gather_T      = np.array(P_gather_T)
	cs_gather_T     = np.array(cs_gather_T)
	qdot_gather_T   = np.array(qdot_gather_T)
	gamma_gather_T  = np.array(gamma_gather_T)

	rho_array    = np.zeros(shape=(nT,nTheta,nR))
	vr_array     = np.zeros(shape=(nT,nTheta,nR))
	vtheta_array = np.zeros(shape=(nT,nTheta,nR))
	S_array      = np.zeros(shape=(nT,nTheta,nR))
	Ye_array     = np.zeros(shape=(nT,nTheta,nR))
	T_array      = np.zeros(shape=(nT,nTheta,nR))
	P_array      = np.zeros(shape=(nT,nTheta,nR))
	cs_array     = np.zeros(shape=(nT,nTheta,nR))
	qdot_array   = np.zeros(shape=(nT,nTheta,nR))
	gamma_array  = np.zeros(shape=(nT,nTheta,nR))
	for i in range(0,iSize):
		for k in range(0,nTheta):
			for j in range(0,nR):
				for l in range(0,size):
					rho_array[l+size*i][k][j]    = rho_gather_T[i][l][k][j]	
					vr_array[l+size*i][k][j]     = vr_gather_T[i][l][k][j]
					vtheta_array[l+size*i][k][j] = vtheta_gather_T[i][l][k][j]
					S_array[l+size*i][k][j]      = S_gather_T[i][l][k][j]
					Ye_array[l+size*i][k][j]     = Ye_gather_T[i][l][k][j]
					T_array[l+size*i][k][j]      = T_gather_T[i][l][k][j]
					P_array[l+size*i][k][j]      = P_gather_T[i][l][k][j]
					cs_array[l+size*i][k][j]     = cs_gather_T[i][l][k][j]
					qdot_array[l+size*i][k][j]   = qdot_gather_T[i][l][k][j]
					gamma_array[l+size*i][k][j]  = gamma_gather_T[i][l][k][j]
	
	def averaging():
		r       = data1['x1v']
		theta   = data1['x2v']
		nR      = np.size(r)
		nTheta  = np.size(theta)

		# time average
		rho_avg_t   = np.zeros(shape=(nTheta,nR))
		vr_avg_t    = np.zeros(shape=(nTheta,nR)) 
		S_avg_t     = np.zeros(shape=(nTheta,nR)) 
		Ye_avg_t    = np.zeros(shape=(nTheta,nR)) 
		T_avg_t     = np.zeros(shape=(nTheta,nR)) 
		P_avg_t     = np.zeros(shape=(nTheta,nR)) 
		cs_avg_t    = np.zeros(shape=(nTheta,nR)) 
		qdot_avg_t  = np.zeros(shape=(nTheta,nR)) 
		gamma_avg_t = np.zeros(shape=(nTheta,nR)) 
		DELTA_T     = time[nT-1] - time[0]
		for i in range(0,nR):
			for j in range(0,nTheta):
				for k in range(0,nT-1):
					delta_T          = time[k+1] - time[k]
					rho_avg_t[j,i]   += (1.0 / DELTA_T) * (0.5) * delta_T * (rho_array[k,j,i]   + rho_array[k+1,j,i])
					vr_avg_t[j,i]    += (1.0 / DELTA_T) * (0.5) * delta_T * (vr_array[k,j,i]    + vr_array[k+1,j,i])
					S_avg_t[j,i]     += (1.0 / DELTA_T) * (0.5) * delta_T * (S_array[k,j,i]     + S_array[k+1,j,i])
					Ye_avg_t[j,i]    += (1.0 / DELTA_T) * (0.5) * delta_T * (Ye_array[k,j,i]    + Ye_array[k+1,j,i])
					T_avg_t[j,i]     += (1.0 / DELTA_T) * (0.5) * delta_T * (T_array[k,j,i]     + T_array[k+1,j,i])
					P_avg_t[j,i]     += (1.0 / DELTA_T) * (0.5) * delta_T * (P_array[k,j,i]     + P_array[k+1,j,i])
					cs_avg_t[j,i]    += (1.0 / DELTA_T) * (0.5) * delta_T * (cs_array[k,j,i]    + cs_array[k+1,j,i])
					qdot_avg_t[j,i]  += (1.0 / DELTA_T) * (0.5) * delta_T * (qdot_array[k,j,i]  + qdot_array[k+1,j,i])
					gamma_avg_t[j,i] += (1.0 / DELTA_T) * (0.5) * delta_T * (gamma_array[k,j,i] + gamma_array[k+1,j,i])
	
		# theta average
		rho_avg_theta    = np.zeros(shape=(nT,nR))
		vr_avg_theta     = np.zeros(shape=(nT,nR)) 
		vtheta_avg_theta = np.zeros(shape=(nT,nR))
		S_avg_theta      = np.zeros(shape=(nT,nR)) 
		Ye_avg_theta     = np.zeros(shape=(nT,nR)) 
		T_avg_theta      = np.zeros(shape=(nT,nR)) 
		P_avg_theta      = np.zeros(shape=(nT,nR)) 
		cs_avg_theta     = np.zeros(shape=(nT,nR)) 
		qdot_avg_theta   = np.zeros(shape=(nT,nR)) 
		gamma_avg_theta  = np.zeros(shape=(nT,nR)) 
		Mdot_avg_theta   = np.zeros(shape=(nT,nR))
		for i in range(0,nR):
			for k in range(0,nT):
				for j in range(0,nTheta-1):
					delta_Theta           = theta[j+1] - theta[j]
					rho_avg_theta[k,i]   += (1.0 / 4.0) * delta_Theta * (rho_array[k,j,i]       * np.sin(theta[j]) + \
						                               		     rho_array[k,j+1,i]     * np.sin(theta[j+1]))
					vr_avg_theta[k,i]    += (1.0 / 4.0) * delta_Theta * (vr_array[k,j,i]        * np.sin(theta[j]) + \
						                 		             vr_array[k,j+1,i]      * np.sin(theta[j+1]))
					vtheta_avg_theta[k,i] += (1.0 / 4.0) * delta_Theta * (vtheta_array[k,j,i]   * np.sin(theta[j]) + \
											      vtheta_array[k,j+1,i] * np.sin(theta[j+1]))
					S_avg_theta[k,i]     += (1.0 / 4.0) * delta_Theta * (S_array[k,j,i]         * np.sin(theta[j]) + \
						                 		             S_array[k,j+1,i]       * np.sin(theta[j+1]))
					Ye_avg_theta[k,i]    += (1.0 / 4.0) * delta_Theta * (Ye_array[k,j,i]        * np.sin(theta[j]) + \
						                 		             Ye_array[k,j+1,i]      * np.sin(theta[j+1]))
					T_avg_theta[k,i]     += (1.0 / 4.0) * delta_Theta * (T_array[k,j,i]         * np.sin(theta[j]) + \
						                 		             T_array[k,j+1,i]       * np.sin(theta[j+1]))
					P_avg_theta[k,i]     += (1.0 / 4.0) * delta_Theta * (P_array[k,j,i]         * np.sin(theta[j]) + \
						                 		             P_array[k,j+1,i]       * np.sin(theta[j+1]))
					cs_avg_theta[k,i]    += (1.0 / 4.0) * delta_Theta * (cs_array[k,j,i]        * np.sin(theta[j]) + \
							         		             cs_array[k,j+1,i]      * np.sin(theta[j+1]))
					qdot_avg_theta[k,i]  += (1.0 / 4.0) * delta_Theta * (qdot_array[k,j,i]      * np.sin(theta[j]) + \
										             qdot_array[k,j+1,i]    * np.sin(theta[j+1]))
					gamma_avg_theta[k,i] += (1.0 / 4.0) * delta_Theta * (gamma_array[k,j,i]     * np.sin(theta[j]) + \
										             gamma_array[k,j+1,i]   * np.sin(theta[j+1]))	

				Mdot_avg_theta[k,i] = 4.0 * np.pi * r[i]**2 * rho_avg_theta[k,i] * np.abs(vr_avg_theta[k,i])

		# time and theta average
		rho_avg_t_theta   = np.zeros(nR)
		vr_avg_t_theta    = np.zeros(nR)
		S_avg_t_theta     = np.zeros(nR)
		Ye_avg_t_theta    = np.zeros(nR)
		T_avg_t_theta     = np.zeros(nR)
		P_avg_t_theta     = np.zeros(nR)
		cs_avg_t_theta    = np.zeros(nR)
		qdot_avg_t_theta  = np.zeros(nR)
		gamma_avg_t_theta = np.zeros(nR)
	
		for i in range(0,nR):
			for j in range(0,nTheta-1):
				delta_Theta = theta[j+1] - theta[j]
				rho_avg_t_theta[i]   += (1.0 / 4.0) * delta_Theta * (rho_avg_t[j,i]     * np.sin(theta[j]) + \
					                 		             rho_avg_t[j+1,i]   * np.sin(theta[j+1]))
				vr_avg_t_theta[i]    += (1.0 / 4.0) * delta_Theta * (vr_avg_t[j,i]      * np.sin(theta[j]) + \
					                 		             vr_avg_t[j+1,i]    * np.sin(theta[j+1]))
				S_avg_t_theta[i]     += (1.0 / 4.0) * delta_Theta * (S_avg_t[j,i]       * np.sin(theta[j]) + \
					                 		             S_avg_t[j+1,i]     * np.sin(theta[j+1]))
				Ye_avg_t_theta[i]    += (1.0 / 4.0) * delta_Theta * (Ye_avg_t[j,i]      * np.sin(theta[j]) + \
					                 		             Ye_avg_t[j+1,i]    * np.sin(theta[j+1]))
				T_avg_t_theta[i]     += (1.0 / 4.0) * delta_Theta * (T_avg_t[j,i]       * np.sin(theta[j]) + \
					                 		             T_avg_t[j+1,i]     * np.sin(theta[j+1]))
				P_avg_t_theta[i]     += (1.0 / 4.0) * delta_Theta * (P_avg_t[j,i]       * np.sin(theta[j]) + \
					                 		             P_avg_t[j+1,i]     * np.sin(theta[j+1]))
				cs_avg_t_theta[i]    += (1.0 / 4.0) * delta_Theta * (cs_avg_t[j,i]      * np.sin(theta[j]) + \
						         		             cs_avg_t[j+1,i]    * np.sin(theta[j+1]))
				qdot_avg_t_theta[i]  += (1.0 / 4.0) * delta_Theta * (qdot_avg_t[j,i]    * np.sin(theta[j]) + \
									             qdot_avg_t[j+1,i]  * np.sin(theta[j+1]))
				gamma_avg_t_theta[i] += (1.0 / 4.0) * delta_Theta * (gamma_avg_t[j,i]   * np.sin(theta[j]) + \
									             gamma_avg_t[j+1,i] * np.sin(theta[j+1]))

	
		return[rho_avg_t, vr_avg_t, S_avg_t, Ye_avg_t, T_avg_t, P_avg_t, cs_avg_t, qdot_avg_t, gamma_avg_t, \
		       rho_avg_theta, vr_avg_theta, vtheta_avg_theta, S_avg_theta, Ye_avg_theta, T_avg_theta, P_avg_theta, cs_avg_theta, qdot_avg_theta, gamma_avg_theta, \
		       Mdot_avg_theta, \
 		       rho_avg_t_theta, vr_avg_t_theta, S_avg_t_theta, Ye_avg_t_theta, T_avg_t_theta, P_avg_t_theta, cs_avg_t_theta, qdot_avg_t_theta, gamma_avg_t_theta]


	def BruntVaisaila_Chi_of_t(rho_avg, vr_avg, T_avg, P_avg, S_avg, qdot_avg, ye_avg, time):
		# Quantities are angle-averaged
		kb    = 8.617333262145*pow(10,-11) # MeV K^-1
		mb    = 1.66053872801*pow(10,-24)  # g
		G     = 6.6743*pow(10,-8)          # cm^3 g^-1 s^-1
		M     = 1.4*2*pow(10,33)           # g
		mu    = G*M
		r     = data1['x1v']
		nR    = np.size(r)
		imod  = 2
	
		g       = np.zeros(nR)
		drhodP  = np.zeros(shape=(nT,nR))
		dPdYe   = np.zeros(shape=(nT,nR))
		dPdS    = np.zeros(shape=(nT,nR))
		dSdr    = np.zeros(shape=(nT,nR))
		dYedr   = np.zeros(shape=(nT,nR))
		OmegaSq = np.zeros(shape=(nT,nR))
		Omega_I = np.zeros(shape=(nT,nR))
	
		print('(BruntVaisaila_Chi_of_t) Calculating derivatives...')
		for i in range(0,nT):	
			for j in range(0,nR):
				g[j] = mu / (r[i]**2)
				if(j==0):
					drhodP[i,j]  = (rho_avg[i,j+1] - rho_avg[i,j]) / (P_avg[i,j+1]  - P_avg[i,j])
					dPdYe[i,j]   = (P_avg[i,j+1]   - P_avg[i,j])   / (ye_avg[i,j+1] - ye_avg[i,j])
					dPdS[i,j]    = (P_avg[i,j+1]   - P_avg[i,j])   / (S_avg[i,j+1]  - S_avg[i,j])
					dSdr[i,j]    = (S_avg[i,j+1]   - S_avg[i,j])   / (r[j+1]        - r[j])
					dYedr[i,j]   = (ye_avg[i,j+1]  - ye_avg[i,j])  / (r[j+1]        - r[j])
				if(j==nR-1):
					drhodP[i,j]  = (rho_avg[i,j] - rho_avg[i,j-1]) / (P_avg[i,j]  - P_avg[i,j-1])
					dPdYe[i,j]   = (P_avg[i,j]   - P_avg[i,j-1])   / (ye_avg[i,j] - ye_avg[i,j-1])
					dPdS[i,j]    = (P_avg[i,j]   - P_avg[i,j-1])   / (S_avg[i,j]  - S_avg[i,j-1])
					dSdr[i,j]    = (S_avg[i,j]   - S_avg[i,j-1])   / (r[j]        - r[j-1])
					dYedr[i,j]   = (ye_avg[i,j]  - ye_avg[i,j-1])  / (r[j]        - r[j-1])
				else:
					drhodP[i,j]  = (rho_avg[i,j+1] - rho_avg[i,j-1]) / (P_avg[i,j+1]  - P_avg[i,j-1])
					dPdYe[i,j]   = (P_avg[i,j+1]   - P_avg[i,j-1])   / (ye_avg[i,j+1] - ye_avg[i,j-1])
					dPdS[i,j]    = (P_avg[i,j+1]   - P_avg[i,j-1])   / (S_avg[i,j+1]  - S_avg[i,j-1])
					dSdr[i,j]    = (S_avg[i,j+1]   - S_avg[i,j-1])   / (r[j+1]        - r[j-1])
					dYedr[i,j]   = (ye_avg[i,j+1]  - ye_avg[i,j-1])  / (r[j+1]        - r[j-1])
	
	
				OmegaSq[i,j] = g[i] / rho_avg[i,j] * drhodP[i,j] * (dPdS[i,j]*dSdr[i,j] + dPdYe[i,j]*dYedr[i,j])
				Omega_I[i,j] = np.imag(cmath.sqrt(OmegaSq[i,j]))
	
		rgain   = np.zeros(nT)
		Mdot    = np.zeros(shape=(nT,nR))
		i_rgain = np.zeros(nT)
		for j in range(0,nT):
			counter = 0
			for i in range(0,nR):
				Mdot[j,i] = 4.0 * math.pi * r[i]**2 * rho_avg[j,i] * np.abs(vr_avg[j,i])
				if(qdot_avg[j,i]>0.0 and counter==0):
					rgain[j] = r[i]
					i_rgain[j] = int(i)
					counter += 1
	
		MdotMax  = np.zeros(nT)
		rshock   = np.zeros(nT)
		i_rshock = np.zeros(nT)
		for j in range(0,nT):
			MdotMax[j] = np.max(Mdot[j,imod:int(nR)-imod])
			for i in range(imod,int(nR)-imod):
				if(Mdot[j,i]==MdotMax[j]):
					rshock[j]   = r[i]
					i_rshock[j] = int(i)
	
		Chi = np.zeros(nT)
		for j in range(0,nT):
			for i in range(int(i_rgain[j]),int(i_rshock[j])-1):
				deltaR  = r[i+1] - r[i]
				Chi[j] += 0.5 * deltaR * (Omega_I[j,i] / np.abs(vr_avg[j,i]) + Omega_I[j,i+1] / np.abs(vr_avg[j,i+1]))	

		Chi_t_avg    = 0.0
		Rshock_avg_t = 0.0
		Rgain_avg_t  = 0.0
		DELTA_T      = time[nT-1] - time[0]
		for j in range(0,nT-1):
			delta_T       = time[j+1] - time[j]
			Chi_t_avg    += (1.0 / DELTA_T) * (0.5) * delta_T * (Chi[j]    + Chi[j+1])
			Rshock_avg_t += (1.0 / DELTA_T) * (0.5) * delta_T * (rshock[j] + rshock[j+1])
			Rgain_avg_t  += (1.0 / DELTA_T) * (0.5) * delta_T * (rgain[j]  + rgain[j+1])

		Chi_t_avg_array = np.zeros(nT)
		Rshock_avg_array = np.zeros(nT)
		Rgain_avg_array  = np.zeros(nT)
		for j in range(0,nT):
			Chi_t_avg_array[j]  = Chi_t_avg 
			Rshock_avg_array[j] = Rshock_avg_t
			Rgain_avg_array[j]  = Rgain_avg_t
	
		return[Chi, Chi_t_avg_array, Chi_t_avg, rgain, rshock, Rshock_avg_array, Rgain_avg_array]

	def timescales(rho_avg, vr_avg, P_avg, qdot_avg):
		# Quantities are angle-averaged
		r      = data1['x1v']
		theta  = data1['x2v']
		nR     = np.size(r)
		nTheta = np.size(theta)
		nT     = np.size(time)
		H      = np.zeros(shape=(nT,nR))
		for j in range(0,nR):
			for k in range(0,nT):
				if(j==0):
					H[k,j] = pow(np.abs((np.log(rho_avg[k,j+1]) - np.log(rho_avg[k,j]))   / (r[j+1] - r[j])),-1.0)
				elif(j==nR-1):
					H[k,j] = pow(np.abs((np.log(rho_avg[k,j])   - np.log(rho_avg[k,j-1])) / (r[j]   - r[j-1])),-1.0)
				else:
					H[k,j] = pow(np.abs((np.log(rho_avg[k,j+1]) - np.log(rho_avg[k,j-1])) / (r[j+1] - r[j-1])),-1.0)

		Tadv       = np.abs(H / vr_avg)
		Theat      = np.abs((P_avg / rho_avg) * (1.0 / qdot_avg))
		Ratio      = Tadv / Theat
		return[Tadv, Theat, Ratio]

	def ReynoldsStressTensor(rhoR, vrR, vthetaR, csR, T):
		# rho, vr, and vtheta are angle-averaged
		# Domains
		R         = data1['x1v']
		r         = data1['x1v']
		DELTA_t   = time[nT-1]-time[0]
		# Reynolds stress terms
		R00       = np.zeros(nR) # r-r
		R01       = np.zeros(nR) # r-theta
		R10       = np.zeros(nR) # theta-r
		R11       = np.zeros(nR) # theta-theta
		R00_g     = np.zeros(nR)
		R01_g     = np.zeros(nR) # r-theta
		R10_g     = np.zeros(nR) # theta-r
		R11_g     = np.zeros(nR) # theta-theta
		# Velocity terms
		v0v0      = np.zeros(nR)
		v0v1      = np.zeros(nR) # v0v1 = v1v0
		v1v1      = np.zeros(nR)
		v0        = np.zeros(nR)
		v1        = np.zeros(nR)	
		# Density term
		rho_bar   = np.zeros(nR)
		cssq_tavg = np.zeros(nR)
		v_esc     = np.zeros(nR)
		# Turbulence diagnostic quantities (Raives et al. 2021)
		xi_turb   = np.zeros(nR)
		xi_th     = np.zeros(nR)
		xi_eff    = np.zeros(nR)
		alpha     = np.zeros(nR)
		K         = np.zeros(nR)
		# Physical constants
		G         = 6.6743*pow(10,-8) # cm^3 g^-1 s^-1
		M         = 1.4*2*pow(10,33)  # g
		mu        = G*M
		
		print('(ReynoldsStressTensor) Calculating Rij elements...')
		# Calculate Rij(r) quantites (equation 24 in Raives et al. 2021)
		for k in range(0,nR):
			v0v0_x_rho = 0.0
			v0v1_x_rho = 0.0
			v1v1_x_rho = 0.0
			v0_x_rho   = 0.0
			v1_x_rho   = 0.0
			v_esc[k]   = np.sqrt(2.0 * mu / r[k])
			# Time integration
			for i in range(0,nT-1):	
				delta_t        = T[i+1] - T[i]
				rho_bar[k]    += (1.0 / DELTA_t) * 0.5 * delta_t * (rhoR[i,k] + rhoR[i+1,k])
				v0v0_x_rho    += (1.0 / DELTA_t) * 0.5 * delta_t * (rhoR[i,k] * vrR[i,k]     * vrR[i,k]     + rhoR[i+1,k] * vrR[i+1,k]     * vrR[i+1,k])
				v0v1_x_rho    += (1.0 / DELTA_t) * 0.5 * delta_t * (rhoR[i,k] * vrR[i,k]     * vthetaR[i,k] + rhoR[i+1,k] * vrR[i+1,k]     * vthetaR[i+1,k])
				v1v1_x_rho    += (1.0 / DELTA_t) * 0.5 * delta_t * (rhoR[i,k] * vthetaR[i,k] * vthetaR[i,k] + rhoR[i+1,k] * vthetaR[i+1,k] * vthetaR[i+1,k])
				v0_x_rho      += (1.0 / DELTA_t) * 0.5 * delta_t * (rhoR[i,k] * vrR[i,k]                    + rhoR[i+1,k] * vrR[i+1,k])
				v1_x_rho      += (1.0 / DELTA_t) * 0.5 * delta_t * (rhoR[i,k] * vthetaR[i,k]                + rhoR[i+1,k] * vthetaR[i+1,k])
				cssq_tavg[k]  += (1.0 / DELTA_t) * 0.5 * delta_t * (csR[i,k]**2                             + csR[i+1,k]**2)
	
			# Divide by rho_bar
			v0v0[k] = v0v0_x_rho / rho_bar[k]
			v0v1[k] = v0v1_x_rho / rho_bar[k]
			v1v1[k] = v1v1_x_rho / rho_bar[k]
			v0[k]   = v0_x_rho   / rho_bar[k]
			v1[k]   = v1_x_rho   / rho_bar[k]
			# Calculate Reynolds stress tensor (units of erg/cm^3)
			R00[k]  = rho_bar[k] * (v0v0[k] - v0[k] * v0[k])
			R01[k]  = rho_bar[k] * (v0v1[k] - v0[k] * v1[k])
			R10[k]  = R01[k]
			R11[k]  = rho_bar[k] * (v1v1[k] - v1[k] * v1[k])

			# (units of erg/g)
			R00_g[k] = R00[k] / rho_bar[k]
			R01_g[k] = R01[k] / rho_bar[k]
			R10_g[k] = R01_g[k]
			R11_g[k] = R11[k] / rho_bar[k]
			# Diagnostics
			alpha[k]   = R11[k] / R00[k]
			if(R00[k]==0.0):
				print('')
				print('R00['+str(k)+']   = ' + str(R00[k]))
				print('v0v0['+str(k)+']  = ' + str(v0v0[k]))
				print('v0*v0['+str(k)+'] = ' + str(v0[k] * v0[k]))
				print('R10['+str(k)+']   = ' + str(R10[k]))
				print('R11['+str(k)+']   = ' + str(R11[k]))
				print('alpha['+str(k)+'] = ' + str(alpha[k]))
				print(' ')
			K[k]       = (1.0 + alpha[k]) * R00[k] / rho_bar[k] # Kinetic energy, cm^2/s^2
			xi_th[k]   = cssq_tavg[k] / v_esc[k]**2
			xi_turb[k] = 2.0 * K[k] / v_esc[k]**2
			xi_eff[k]  = xi_th[k] + 0.5 * xi_turb[k]
	
	
		return[R00, R01, R10, R11, R00_g, R01_g, R10_g, R11_g, v0v0, v0v1, v1v1, v0, v1, rho_bar, alpha, K, xi_th, xi_turb, xi_eff]
	

	def MachNumber(VrM, CsM):
		nT   = np.size(time)
		nR   = np.size(R)
		mach = np.zeros(nT)
		vr_o = np.zeros(nT)
		cs_o = np.zeros(nT)
		for i in range(0,nT):
			mach[i] = np.abs(VrM[i,nR-1] / CsM[i,nR-1])
			vr_o[i] = np.abs(VrM[i,nR-1])
			cs_o[i] = np.abs(CsM[i,nR-1])
	
		return [mach, vr_o, cs_o]

	print('Calculating averages...')	
	R = data1['x1v']
	[rho_avg_t, vr_avg_t, S_avg_t, Ye_avg_t, T_avg_t, P_avg_t, cs_avg_t, qdot_avg_t, gamma_avg_t, \
	rho_avg_theta, vr_avg_theta, vtheta_avg_theta, S_avg_theta, Ye_avg_theta, T_avg_theta, P_avg_theta, cs_avg_theta, qdot_avg_theta, gamma_avg_theta, Mdot_avg_theta, \
	rho_avg_t_theta, vr_avg_t_theta, S_avg_t_theta, Ye_avg_t_theta, T_avg_t_theta, P_avg_t_theta, cs_avg_t_theta, qdot_avg_t_theta, gamma_avg_t_theta] = \
	averaging()	
	[Chi, Chi_t_avg, Chi_t_avg_value, Rgain, Rshock, Rshock_avg_t, Rgain_avg_t] = \
	BruntVaisaila_Chi_of_t(rho_avg_theta, vr_avg_theta, T_avg_theta, P_avg_theta, S_avg_theta, qdot_avg_theta, Ye_avg_theta, time)
	[Tadv, Theat, Ratio] = timescales(rho_avg_theta, vr_avg_theta, P_avg_theta, qdot_avg_theta)
	[R00, R01, R10, R11, R00_g, R01_g, R10_g, R11_g, v0v0, v0v1, v1v1, v0, v1, rho_bar, alpha, K, xi_th, xi_turb, xi_eff] = \
	ReynoldsStressTensor(rho_avg_theta, vr_avg_theta, vtheta_avg_theta, cs_avg_theta, time)
	[MachNum, Vr_Outside, Cs_Outside] = MachNumber(vr_avg_theta, cs_avg_theta)
	

	######### Entropy-Vr plot parameters ###########
	# LHS, entropy
	min_S = np.amin(S_array)
	max_S = np.amax(S_array)
	# RHS, vr
	min_vr = -2.75*pow(10,9)
	max_vr = 2.75*pow(10,9)
	# colorbar?	
	zc   = 8e9
	
	arr = []
	for i in range(istart,nmax+1):
		if(i%1==0):
			arr.append(int(i))	

	fig,axs = plt.subplots(1,1,figsize=(10,10),subplot_kw=dict(projection='polar'))
	bx      = fig.add_axes([0.005,0.1,0.015,0.72])
	cx      = fig.add_axes([0.93,0.1,0.015,0.72])
	bx.yaxis.set_ticks_position("right")
	cx.yaxis.set_ticks_position("left")
	norm1   = mpl.colors.LogNorm(vmin=min_S, vmax=max_S)
	hot     = mpl.colormaps['hot']
	seismic = mpl.colormaps['seismic']
	cb1     = mpl.colorbar.ColorbarBase(bx,cmap=hot,norm=norm1,orientation='vertical',label=r"$S$ [ $k_{b}$ baryon$^{-1}$] ]")
	norm2   = mpl.colors.Normalize(vmin=min_vr, vmax=max_vr,clip=True)
	cb2     = mpl.colorbar.ColorbarBase(cx,cmap=seismic,norm=norm2,orientation='vertical',label=r"$v_{r}$ [ cm/s ]")	

	######### Gamma-T plot parameters ###########
	# LHS, gamma
	min_gamma = np.amin(gamma_array)
	max_gamma = np.amax(gamma_array)
	# RHS, T
	min_T = np.amin(T_array) * kb_MeV
	max_T = np.amax(T_array) * kb_MeV
	# colorbar?	
	zc2   = 8e9
	
	arr2 = []
	for i in range(istart,nmax+1):
		if(i%1==0):
			arr2.append(int(i))	

	fig2,axs2 = plt.subplots(1,1,figsize=(10,10),subplot_kw=dict(projection='polar'))
	bx2       = fig2.add_axes([0.005,0.1,0.015,0.72])
	cx2       = fig2.add_axes([0.93,0.1,0.015,0.72])
	bx2.yaxis.set_ticks_position("right")
	cx2.yaxis.set_ticks_position("left")
	norm1_2   = mpl.colors.Normalize(vmin=min_gamma, vmax=max_gamma)	
	inferno   = mpl.colormaps['inferno']
	cb1_2     = mpl.colorbar.ColorbarBase(bx2,cmap=hot,norm=norm1_2,orientation='vertical',label=r"$\gamma$")
	norm2_2   = mpl.colors.LogNorm(vmin=min_T, vmax=max_T,clip=True)
	cb2_2     = mpl.colorbar.ColorbarBase(cx2,cmap=inferno,norm=norm2_2,orientation='vertical',label=r"$T$ [ MeV ]")

	######### Entropy-Vr reduced rmax plot parameters ###########
	# LHS, entropy
	rtest = data1['x1v']
	n_rtest = np.size(rtest)
	i_SrMax = int(n_rtest)
	for l in range(0,int(np.size(rtest))):
		if(rtest[l]>SrMax):
			i_SrMax = l
			break

	print('entropy shape [][][]  = ' + str(np.shape(S_array[:][:][0:i_SrMax])))
	print('entropy shape [:,:,:] = ' + str(np.shape(S_array[:,:,0:i_SrMax])))
	min_S2 = np.amin(S_array[:,:,0:i_SrMax])
	max_S2 = np.amax(S_array[:,:,0:i_SrMax])
	# Using same vr and colorbar options as before
	
	arr3 = []
	for i in range(istart,nmax+1):
		if(i%1==0):
			arr3.append(int(i))	

	fig3,axs3 = plt.subplots(1,1,figsize=(10,10),subplot_kw=dict(projection='polar'))
	bx3       = fig3.add_axes([0.005,0.1,0.015,0.72])
	cx3       = fig3.add_axes([0.93,0.1,0.015,0.72])
	bx3.yaxis.set_ticks_position("right")
	cx3.yaxis.set_ticks_position("left")
	norm1b   = mpl.colors.LogNorm(vmin=min_S2, vmax=max_S2)	
	cb1b     = mpl.colorbar.ColorbarBase(bx3,cmap=hot,norm=norm1b,orientation='vertical',label=r"$S$ [ $k_{b}$ baryon$^{-1}$] ]")
	norm2b   = mpl.colors.Normalize(vmin=min_vr, vmax=max_vr,clip=True)
	cb2b     = mpl.colorbar.ColorbarBase(cx3,cmap=seismic,norm=norm2b,orientation='vertical',label=r"$v_{r}$ [ cm/s ]")

	def animate_Multi(i):
		print('(animate multipanel) i = ' + str(i))
	
		tformat  = '{:.2f}'.format
		plt.suptitle('t = '+tformat(i*float(ts)/float(max_index) * 1.0e3)+' ms')
		plt.tight_layout(pad=0.1)
		fig_multi.subplots_adjust(top=0.88)
		r = data1['x1v']
	
		for j in range(0,len(lines)):
			if(j==0):
				lines[j].set_data(r/1.0e5,rho_avg_theta[i,:])
			if(j==1):
				lines[j].set_data(r/1.0e5,T_avg_theta[i,:]*kb_MeV)
			if(j==2):
				lines[j].set_data(r/1.0e5,P_avg_theta[i,:])
			if(j==3):
				lines[j].set_data(r/1.0e5,np.abs(vr_avg_theta[i,:]))
			if(j==4):
				lines[j].set_data(r/1.0e5,cs_avg_theta[i,:])
			if(j==5):
				lines[j].set_data(r/1.0e5,Mdot_avg_theta[i,:]/1.989e33)
			if(j==6):
				lines[j].set_data(r/1.0e5,Mdot_EXP)
			if(j==7):
				lines[j].set_data(r/1.0e5,qdot_avg_theta[i,:]/1.0e21)
		return lines

	def animate_Tau(i):
		print('(animate tau) i = ' + str(i))
	
		tformat  = '{:.2f}'.format	
		plt.suptitle('t = '+tformat(i*float(ts)/float(max_index) * 1.0e3)+' ms')
		plt.tight_layout(pad=0.1)
		fig_tau.subplots_adjust(top=0.88)
		r = data1['x1v']
	
		for j in range(0,len(lines2)):
			if(j==0):
				lines2[j].set_data(r/1.0e5,Tadv[i,:])
			if(j==1):
				lines2[j].set_data(r/1.0e5,Theat[i,:])
		return lines2

	def animate_Entropy_avg_theta(i):
		print('(animate S avg theta) i = ' + str(i))
	
		tformat  = '{:.2f}'.format	
		plt.suptitle('t = '+tformat(i*float(ts)/float(max_index) * 1.0e3)+' ms')
		plt.tight_layout(pad=0.1)
		fig_entropy.subplots_adjust(top=0.88)
		r = data1['x1v']
	
		for j in range(0,len(linesS)):
			if(j==0):
				linesS[j].set_data(r/1.0e5,S_avg_theta[i,:])	
		return linesS

	########## Multipanel plot parameters ###########################
	fig_multi, [(ax1, ax2, ax3), (ax4, ax5, ax6)] = plt.subplots(2,3)
	lcolor  = 'b'
	kcolor  = 'black'
	line1,  = ax1.loglog([],   [], lw=1.75, color=lcolor)
	line2,  = ax2.loglog([],   [], lw=1.75, color=lcolor)
	line3,  = ax3.loglog([],   [], lw=1.75, color=lcolor)
	line4,  = ax4.loglog([],   [], lw=1.75, color=lcolor)
	line4a, = ax4.loglog([],   [], lw=1.75, color=kcolor,linestyle='dashed')
	line5,  = ax5.loglog([],   [], lw=1.75, color=lcolor)
	line5a, = ax5.loglog([],   [], lw=1.75, color=kcolor, linestyle='dashed')
	line6,  = ax6.semilogx([], [], lw=1.75, color=lcolor)


	rmin = min(data1['x1v']) # min radius in cm	
	# Density
	ax1.set_xlim(rmin/1.0e5,rmax/1.0e5)
	ax1.set_ylim(np.amin(rho_avg_theta),np.amax(rho_avg_theta))
	ax1.set_xlabel(r'$r$ [km]')
	ax1.set_ylabel(r'$\rho$ [g/cm$^{3}$]')
	
	# Temperature
	ax2.set_xlim(rmin/1.0e5,rmax/1.0e5)
	ax2.set_ylim(np.amin(T_avg_theta*kb_MeV),np.amax(T_avg_theta*kb_MeV))
	ax2.set_xlabel(r'$r$ [km]')
	ax2.set_ylabel(r'$T$ [MeV]')
	
	# Pressure
	ax3.set_xlim(rmin/1.0e5,rmax/1.0e5)
	ax3.set_ylim(np.amin(P_avg_theta),np.amax(P_avg_theta))
	ax3.set_xlabel(r'$r$ [km]')
	ax3.set_ylabel(r'$P$ [erg/cm$^{3}$]')

	velmin = np.amin(np.abs(vr_avg_theta))
	velmax = np.amax(np.abs(vr_avg_theta))
	csmin  = np.amin(cs_avg_theta)
	csmax  = np.amax(cs_avg_theta)
	minValues = [velmin, csmin]
	maxValues = [velmax, csmax]
	vplotmin = np.min(minValues)
	vplotmax = np.max(maxValues)
	
	# Velocity
	ax4.set_xlim(rmin/1.0e5,rmax/1.0e5)
	ax4.set_ylim(vplotmin,vplotmax)
	ax4.set_xlabel(r'$r$ [km]')
	ax4.set_ylabel(r'$\left|v_{r}\right|$ (blue), $C_{s}$ (black) [cm/s$^{1}$]')
	
	# Mdot
	ax5.set_xlim(rmin/1.0e5,rmax/1.0e5)
	ax5.set_ylim(np.amin(Mdot_avg_theta)/(1.989e33),np.amax(Mdot_avg_theta)/(1.989e33))
	ax5.set_xlabel(r'$r$ [km]')
	ax5.set_ylabel(r'$\dot{M}$ [M$_{\odot}$ s$^{-1}$]')
	
	# Total Heating
	ax6.set_xlim(rmin/1.0e5,rmax/1.0e5)
	ax6.set_ylim(np.amin(qdot_avg_theta)/(1.0e21),np.amax(qdot_avg_theta)/(1.0e21))
	ax6.set_xlabel(r'$r$ [km]')
	ax6.set_ylabel(r'$\dot{q}$ [$10^{21}$ erg/g/s]')

	Mdot_EXP = np.zeros(nR)
	for j in range(0,nR):
		Mdot_EXP[j] = Mdot_set  # Msun/s	

	Mach_EXP = np.zeros(nT)	
	frac     = np.zeros(nT)
	ones     = np.ones(nT)
	for j in range(0,nT):
		Mach_EXP[j] = Mach_set
		frac[j]     = np.abs(MachNum[j] - Mach_EXP[j]) / np.abs(Mach_EXP[j])

	lines = [line1, line2, line3, line4, line4a, line5, line5a, line6]
	arrM  = []
	k     = start_index
	while(int(k)<=int(max_index)):
		if(int(k)%1==0):
			arrM.append(int(k))
		k=int(k)+1

	def init():
		for i in range(0,len(lines)):
			lines[i].set_data([],[])
		return lines

	print('######### Generating multi-plot #########')
	animM = FuncAnimation(fig_multi, animate_Multi, init_func=init,
                               frames=arrM, interval=100, blit=True,repeat_delay=200)
	print(arrM)
	animM.save(name+'_Multipanel.gif', writer='imagemagick')

	############ Tadv and Theat plot parameters ###########
	fig_tau    = plt.figure()
	Tadvmin    = np.amin(Tadv)
	Tadvmax    = np.amax(Tadv)
	Theatmin   = np.amin(Theat)
	Theatmax   = np.amax(Theat)
	minTValues = [Tadvmin, Theatmin]
	maxTValues = [Tadvmax, Theatmax]
	Tplotmin   = np.min(minTValues)
	Tplotmax   = np.max(maxTValues)
	print('Tplotmin = ' + str(Tplotmin))
	print('Tplotmax = ' + str(Tplotmax))	
	axT        = \
		plt.axes(xlim=(rmin/1.0e5,rmax/1.0e5), ylim=(Tplotmin, Tplotmax),xlabel=r"$R$ [km]",ylabel=r"$\tau_{\mathrm{adv}}$ (red), $\tau_{\mathrm{heat}}$ (blue)")
	line1a     = axT.loglog([],   [], lw=1.75, color='red', label=r'$\tau_{\mathrm{adv}}$')[0]
	line1b     = axT.loglog([],   [], lw=1.75, color='blue', linestyle='dashed', label=r'$\tau_{\mathrm{heat}}$')[0]	
	lines2     = [line1a, line1b]
	arrT       = []

	j = start_index
	while(int(j)<=int(max_index)):
		if(int(j)%1==0):
			arrT.append(int(j))
		j=int(j)+1
	
	def init2():
		for i in range(0,len(lines2)):
			lines2[i].set_data([],[])
		return lines2

	print('######### Generating Tau #########')
	animT = FuncAnimation(fig_tau, animate_Tau, init_func=init2,
                               frames=arrT, interval=100, blit=True,repeat_delay=200)
	print(arrT)
	animT.save(name+'_Tadv_Theat.gif', writer='imagemagick')

	############ S angle average plot parameters ###########
	fig_entropy = plt.figure()
	entropy_min = np.amin(S_avg_theta)
	entropy_max = np.amax(S_avg_theta)
	axEntropy   = \
		plt.axes(xlim=(rmin/1.0e5,rmax/1.0e5), ylim=(entropy_min, entropy_max),xlabel=r"$R$ [km]",ylabel=r"$\left<S\right>_{\theta}$ [k$_{\mathrm{b}}$ baryon$^{-1}$]")
	line1s = axEntropy.loglog([],   [], lw=1.75, color='red')[0]	
	linesS = [line1s]
	arrEnt   = []

	j = start_index
	while(int(j)<=int(max_index)):
		if(int(j)%1==0):
			arrEnt.append(int(j))
		j=int(j)+1
	
	def initEnt():
		for i in range(0,len(linesS)):
			linesS[i].set_data([],[])
		return linesS

	print('######### Generating Entropy angle average plot #########')
	animEnt = FuncAnimation(fig_entropy, animate_Entropy_avg_theta, init_func=initEnt,
                               frames=arrEnt, interval=100, blit=True,repeat_delay=200)
	print(arrEnt)
	animEnt.save(name+'_Entropy_theta_avg.gif', writer='imagemagick')
	
	def animate_S_vr(i):
		print('(animate S-vr) i = ' + str(i))
		axs.clear()	

		r1      = data1['x1f']
		theta1  = data1['x2f']	
		Radius  = data1['x1v']
		Angle   = data1['x2v']
		ntheta1 = np.size(theta1)
		vr1     = (data1['vel1'])[0]
		
		index=0
		for j in range(0, len(r1)):
			if(r1[j]>=rmax):
				index=j
				break
		if(rmax==1e9):
			index=len(r1)-2
		
		Vr         = vr_array[i,:,:index+1]#vr_array[i][:][:index+1]
		Rho        = rho_array[i,:,:index+1]	#rho_array[i][:][:index+1]	
		Entropy    = S_array[i,:,:index+1]#S_array[i][:][:index+1]
		r1         = r1[:index+2]
		r2, theta2 = np.meshgrid(r1, theta1)		
		p1         = axs.pcolormesh(theta2,r2,Entropy[0:ntheta1,0:index+1],norm=norm1,cmap='hot') # works with index+2 
		p2         = axs.pcolormesh((-1.0)*theta2,r2,Vr[0:ntheta1,0:index+1],norm=norm2,cmap='seismic') # works with index+2 	
		st_nth     = 32
		st_nr      = 32
		st_pts     = np.linspace(data1['x2f'][1], data1['x2f'][-2], st_nth)
		r0         = data1['x1v'][0]
		st_pts     = np.array([st_pts, [r0] * st_nr]).T
		interp_nth = len(data1['x2v'])
		theta      = np.linspace(0, np.pi, interp_nth)
		interp_nr  = int(np.ceil(rmax / (data1['x1f'][1] - data1['x1f'][0])))	
		r          = np.linspace((data1['x1v'])[0], rmax, interp_nr)	
		axs.set_theta_offset(0.5*np.pi)
		axs.set_ylim([0.0,np.max(r)])	
		axs.set_title('Physical time = %f s'%(i*ts/nmax))

	def animate_S_vr_SrMax(i):
		print('(animate S-vr, reduced rmax) i = ' + str(i))
		axs3.clear()	

		r1      = data1['x1f']
		theta1  = data1['x2f']	
		Radius  = data1['x1v']
		Angle   = data1['x2v']
		ntheta1 = np.size(theta1)
		vr1     = (data1['vel1'])[0]
		
		index=0
		for j in range(0, len(r1)):
			if(r1[j]>=SrMax):
				index=j
				break
		if(SrMax==1e9):
			index=len(r1)-2
		
		Vr         = vr_array[i,:,:index+1]#vr_array[i][:][:index+1]
		Rho        = rho_array[i,:,:index+1]	#rho_array[i][:][:index+1]	
		Entropy    = S_array[i,:,:index+1]#S_array[i][:][:index+1]
		r1         = r1[:index+2]	
		r2, theta2 = np.meshgrid(r1, theta1)		
		p1         = axs3.pcolormesh(theta2,r2,Entropy[0:ntheta1,0:index+1],norm=norm1b,cmap='hot') # works with index+2 
		p2         = axs3.pcolormesh((-1.0)*theta2,r2,Vr[0:ntheta1,0:index+1],norm=norm2b,cmap='seismic') # works with index+2 	
		st_nth     = 32
		st_nr      = 32
		st_pts     = np.linspace(data1['x2f'][1], data1['x2f'][-2], st_nth)
		r0         = data1['x1v'][0]
		st_pts     = np.array([st_pts, [r0] * st_nr]).T
		interp_nth = len(data1['x2v'])
		theta      = np.linspace(0, np.pi, interp_nth)
		interp_nr  = int(np.ceil(SrMax / (data1['x1f'][1] - data1['x1f'][0])))	
		r          = np.linspace((data1['x1v'])[0], SrMax, interp_nr)	
		axs3.set_theta_offset(0.5*np.pi)
		axs3.set_ylim([0.0,np.max(r)])	
		axs3.set_title('Physical time = %f s'%(i*ts/nmax))

	def animate_Gamma_T(i):
		print('(animate gamma-T) i = ' + str(i))
		axs2.clear()	

		r1      = data1['x1f']
		theta1  = data1['x2f']
		Radius  = data1['x1v']
		Angle   = data1['x2v']
		vr1     = (data1['vel1'])[0]
		
		index=0
		for j in range(0, len(r1)):
			if(r1[j]>=rmax):
				index=j
				break
		if(rmax==1e9):
			index=len(r1)-2
		
		Gamma      = gamma_array[i][:][:index+1]
		Temp       = T_array[i][:][:index+1] * kb_MeV	

		r1         = r1[:index+2]
		r2, theta2 = np.meshgrid(r1, theta1)
		p1         = axs2.pcolormesh(theta2,r2,Gamma,norm=norm1_2,cmap='hot')
		p2         = axs2.pcolormesh((-1.0)*theta2,r2,Temp,norm=norm2_2,cmap='inferno')		
		st_nth     = 32
		st_nr      = 32
		st_pts     = np.linspace(data1['x2f'][1], data1['x2f'][-2], st_nth)
		r0         = data1['x1v'][0]
		st_pts     = np.array([st_pts, [r0] * st_nr]).T
		interp_nth = len(data1['x2v'])
		theta      = np.linspace(0, np.pi, interp_nth)
		interp_nr  = int(np.ceil(rmax / (data1['x1f'][1] - data1['x1f'][0])))	
		r          = np.linspace((data1['x1v'])[0], rmax, interp_nr)	
		axs2.set_theta_offset(0.5*np.pi)
		axs2.set_ylim([0.0,np.max(r)])	
		axs2.set_title('Physical time = %f s'%(i*ts/nmax))

	print('########## Generating S-vr plot #########')
	anim1 = animation.FuncAnimation(fig,animate_S_vr,frames=arr,repeat=False)
	print(arr)
	anim1.save(name+'_entropy_vr.gif',writer='imagemagick')	

	print('########## Generating Gamma-T plot #########')
	anim2 = animation.FuncAnimation(fig2,animate_Gamma_T,frames=arr2,repeat=False)
	print(arr2)
	anim2.save(name+'_Gamma_T.gif',writer='imagemagick')				

	print('########## Generating S-vr (reduced rmax) plot #########')
	anim3 = animation.FuncAnimation(fig3,animate_S_vr_SrMax,frames=arr3,repeat=False)
	print(arr3)
	anim3.save(name+'_entropy_vr_SrMax.gif',writer='imagemagick')	

	print('######## Generating Chi(t) plot #########')
	plt.figure(0)
	plt.clf()
	plt.plot(time,Chi,linewidth=2.0,linestyle='solid',color='black',label=r'$\left< \chi \right>_{\theta}$')
	plt.plot(time,Chi_t_avg,linewidth=1.5,linestyle='dashed',color='red',label=r'$\left< \chi \right>_{\theta,t}$')
	plt.xlabel(r'$t$ [s]')
	plt.ylabel(r'$\left< \chi \right>_{\theta}$')
	plt.title(r'$\left< \chi \right>_{\theta,t} = $' + str(Chi_t_avg_value))
	plt.legend()
	plt.tight_layout()
	plt.savefig(name+'_Chi_t.png')		

	print('######## Generating Rshock and Rgain plot #########')
	plt.figure(0)
	plt.clf()
	plt.plot(time,Rshock,linewidth=2.0,linestyle='solid',color='black',label=r'$\left< R_{\mathrm{shock}} \right>_{\theta}$')
	plt.plot(time,Rgain, linewidth=2.0,linestyle='solid',color='red',  label=r'$\left< R_{\mathrm{gain}}  \right>_{\theta}$')
	plt.plot(time,Rshock_avg_t,linewidth=1.5,linestyle='dashed',color='black',label=r'$\left< R_{\mathrm{shock}} \right>_{\theta,t} = $' + \
		str(Rshock_avg_t[0]/1.0e5) + ' km')
	plt.plot(time,Rgain_avg_t, linewidth=1.5,linestyle='dashed',color='red',  label=r'$\left< R_{\mathrm{gain}}  \right>_{\theta,t} = $' + \
		str(Rgain_avg_t[0]/1.0e5) + ' km')
	plt.xlabel(r'$t$ [s]')
	plt.ylabel(r'$[km]$')
	plt.legend()
	plt.tight_layout()
	plt.savefig(name+'_Rshock_Rgain.png')		


	print('######## Generating Rij plot #########')
	r = data1['x1v']
	plt.figure(0)
	plt.clf()
	plt.loglog(r,R00_g,linewidth=2.0,linestyle='solid',color='black',label=r'$\left<R_{0,0}\right>/\left< \rho \right>$')
	plt.loglog(r,R01_g,linewidth=2.0,linestyle='solid',color='red',  label=r'$\left<R_{0,1}\right>/\left< \rho \right>$')
	plt.loglog(r,R11_g,linewidth=2.0,linestyle='solid',color='blue', label=r'$\left<R_{1,1}\right>/\left< \rho \right>$')	
	plt.xlabel(r'$R$ [km]')
	plt.ylabel(r'$\left<R_{i,j}\right>$ [erg/g]')
	plt.legend()
	plt.tight_layout()
	plt.savefig(name+'_Rij_by_rho.png')		

	print('######## Generating cs^2/vesc^2 plot #########')
	r = data1['x1v']
	plt.figure(0)
	plt.clf()
	plt.loglog(r,xi_th,linewidth=2.0,linestyle='solid',color='black')
	plt.xlabel(r'$R$ [km]')
	plt.ylabel(r'$\left< \xi \right>_{t,\theta}$')
	plt.legend()
	plt.tight_layout()
	plt.savefig(name+'_cssq_vescsq.png')		

	print('######## Generating alpha plot #########')
	r = data1['x1v']
	plt.figure(0)
	plt.clf()
	plt.loglog(r,alpha,linewidth=2.0,linestyle='solid',color='black')
	plt.xlabel(r'$R$ [km]')
	plt.ylabel(r'$\alpha = R_{11} / R_{00}$')
	plt.legend()
	plt.tight_layout()
	plt.savefig(name+'_alpha.png')

	print('######## Generating Mach plot #########')
	r = data1['x1v']
	plt.figure(0)
	plt.clf()
	plt.loglog(time,frac,linewidth=2.0,linestyle='solid',color='black',label=r'$\mathcal{M}$ [code]')
#	plt.loglog(time,ones,linewidth=1.5,linestyle='dashed',color='black')
	plt.xlabel(r'$time$ [s]')
	plt.ylabel(r'$\left(\mathcal{M}_{\mathrm{code}} - \mathcal{M}_{\mathrm{expected}} \right)/ \mathcal{M}_{\mathrm{expected}}$')
#	plt.legend()
	plt.tight_layout()
	plt.savefig(name+'_Mach.png')
