import h5py
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
sys.path.insert(0, '/home/dpochik/mhd-winds-new-EOS/vis/python')
import athena_read
import cmath

DIR_name1 = str(sys.argv[1])
DIR_name2 = str(sys.argv[1])
nmin      = int(sys.argv[2])      # Initial output file number
nmax      = int(sys.argv[3])      # Final output file number
nT        = int(nmax-nmin+1)
minTime   = 0.0                   # Start time of simulation (seconds)
maxTime   = float(sys.argv[4])    # End time of simulation (seconds)
minR      = float(sys.argv[5])    # Minimum Radius in cm
maxR      = float(sys.argv[6])    # Maximum Radius in cm

#'/mnt/d/Research_Stuff/Automated_Analysis_DIR/2D_Data_DIR/Lnu_10_Mdot_0pt4_256x128_6_19_2023_DIR/data_DIR'

nprim   = DIR_name1 + "/accretion.prim.00000.athdf"
datap   = athena_read.athdf(nprim)
r       = datap['x1v']
theta   = datap['x2v']
radRes  = np.size(r)     # radial resolution.
thetRes = np.size(theta) # angular resolution.
Time    = np.linspace(minTime,maxTime,nT).T

iModifier = int(2)       # Removes endpoints in domain while finding max Mdot

def ReadData(mini,maxi,minT,maxT,DIR1,DIR2,nR,nTheta,i_mod):
	# Time domain
	nT   = int(maxi-mini+1)
	Time = np.linspace(minT,maxT,nT)

	# 3D arrays
	rho_3D    = np.zeros(shape=(nR,nTheta,nT))
	vr_3D     = np.zeros(shape=(nR,nTheta,nT))
	vtheta_3D = np.zeros(shape=(nR,nTheta,nT))
	cs_3D     = np.zeros(shape=(nR,nTheta,nT))
	T_3D      = np.zeros(shape=(nR,nTheta,nT))
	P_3D      = np.zeros(shape=(nR,nTheta,nT))
	S_3D      = np.zeros(shape=(nR,nTheta,nT))
	Ye_3D     = np.zeros(shape=(nR,nTheta,nT))
	qdot_3D   = np.zeros(shape=(nR,nTheta,nT))

	# Theta average quantities
	rho_theta_avg    = np.zeros(shape=(nR,nT))
	vr_theta_avg     = np.zeros(shape=(nR,nT))
	vtheta_theta_avg = np.zeros(shape=(nR,nT))
	cs_theta_avg     = np.zeros(shape=(nR,nT))
	T_theta_avg      = np.zeros(shape=(nR,nT))
	P_theta_avg      = np.zeros(shape=(nR,nT))
	qdot_theta_avg   = np.zeros(shape=(nR,nT))
	S_theta_avg      = np.zeros(shape=(nR,nT))
	Ye_theta_avg     = np.zeros(shape=(nR,nT))
	R_shock          = np.zeros(maxi-mini+1)
	R_gain           = np.zeros(maxi-mini+1)

	for i in range(mini,maxi+1):
		if(i<10):
			nprim = DIR1+"/accretion.prim.0000"+str(i)+".athdf"
			nuov  = DIR2+"/accretion.uov.0000"+str(i)+".athdf"
		if(i>=10 and i<100):
			nprim = DIR1+"/accretion.prim.000"+str(i)+".athdf"
			nuov  = DIR2+"/accretion.uov.000"+str(i)+".athdf"
		if(i>=100 and i<1000):
			nprim = DIR1+"/accretion.prim.00"+str(i)+".athdf"
			nuov  = DIR2+"/accretion.uov.00"+str(i)+".athdf"
		if(i>=1000):
			nprim = DIR1+"/accretion.prim.0"+str(i)+".athdf"
			nuov  = DIR2+"/accretion.uov.0"+str(i)+".athdf"

		print('(PlotEverything) i = ' + str(i) )
		datap   = athena_read.athdf(nprim)
		datau   = athena_read.athdf(nuov)
		r       = datap['x1v']
		theta   = datap['x2v']
		vr      = (datap['vel1'])[0]
		vthet   = (datap['vel2'])[0]
		vphi    = (datap['vel3'])[0]
		rho     = (datap['rho'])[0]
		qdot    = (datau['dt2'])[0]
		cs      = (datau['dt3'])[0]
		P       = (datap['press'])[0]
		Ye      = (datap['r0'])[0]
		T       = (datau['dt1'])[0]
		pi      = math.pi
#		nR      = np.size(rho[0,:])
#		nTheta  = np.size(rho[:,0])

		# Constants
		kb      = 8.617333262145e-11         # MeV K^-1
		kb_erg  = 1.380649e-16               # erg/K
		mb      = 1.66053872801*pow(10,-24)  # g
		G       = 6.6743*pow(10,-8)          # cm^3 g^-1 s^-1
		c       = 2.99792458e10              # cm/s
		M       = 1.4*2*pow(10,33)           # g
		hbar    = 1.0546e-27                 # erg s
		mu      = G*M
		gamma   = 4.0/3.0
		pi      = math.pi

		T_MeV = T * kb
		# Build 3D arrays
		for k in range(0,nR):
			for j in range(0,nTheta):
				S_3D[k,j,i] = (11.0 * pi**2 / 45.0) * (1.0 / pow(hbar * c , 3)) * pow(T[j,k] * kb_erg,3) / (rho[j,k] / mb)
				rho_3D[k,j,i]    = rho[j,k]
				vr_3D[k,j,i]     = vr[j,k]
				vtheta_3D[k,j,i] = vthet[j,k]
				cs_3D[k,j,i]     = cs[j,k]
				T_3D[k,j,i]      = T[j,k]
				P_3D[k,j,i]      = P[j,k]
				Ye_3D[k,j,i]     = Ye[j,k]
				qdot_3D[k,j,i]   = qdot[j,k]

		Mdot     = np.zeros(nR)
		# Calculate angle averages with trapezoid rule (2pi/4pi * 1/2 = 1/4 prefactor)
		for k in range(0,nR):
			for j in range(0,nTheta-1):
				deltaTheta              = theta[j+1]-theta[j]
				rho_theta_avg[k,i]     += 1.0 / (4.0) * deltaTheta * (rho_3D[k,j,i]    * np.sin(theta[j]) + rho_3D[k,j+1,i]    * np.sin(theta[j+1]))
				vr_theta_avg[k,i]      += 1.0 / (4.0) * deltaTheta * (vr_3D[k,j,i]     * np.sin(theta[j]) + vr_3D[k,j+1,i]     * np.sin(theta[j+1]))
				vtheta_theta_avg[k,i]  += 1.0 / (4.0) * deltaTheta * (vtheta_3D[k,j,i] * np.sin(theta[j]) + vtheta_3D[k,j+1,i] * np.sin(theta[j+1]))
				cs_theta_avg[k,i]      += 1.0 / (4.0) * deltaTheta * (cs_3D[k,j,i]     * np.sin(theta[j]) + cs_3D[k,j+1,i]     * np.sin(theta[j+1]))
				T_theta_avg[k,i]       += 1.0 / (4.0) * deltaTheta * (T_3D[k,j,i]      * np.sin(theta[j]) + T_3D[k,j+1,i]      * np.sin(theta[j+1]))
				P_theta_avg[k,i]       += 1.0 / (4.0) * deltaTheta * (P_3D[k,j,i]      * np.sin(theta[j]) + P_3D[k,j+1,i]      * np.sin(theta[j+1]))
				S_theta_avg[k,i]       += 1.0 / (4.0) * deltaTheta * (S_3D[k,j,i]      * np.sin(theta[j]) + S_3D[k,j+1,i]      * np.sin(theta[j+1]))
				qdot_theta_avg[k,i]    += 1.0 / (4.0) * deltaTheta * (qdot_3D[k,j,i]   * np.sin(theta[j]) + qdot_3D[k,j+1,i]   * np.sin(theta[j+1]))
				Ye_theta_avg[k,i]      += 1.0 / (4.0) * deltaTheta * (Ye_3D[k,j,i]     * np.sin(theta[j]) + Ye_3D[k,j+1,i]     * np.sin(theta[j+1]))

			Mdot[k] = 4.0 * pi * r[k]**2 * rho_theta_avg[k,i] * np.abs(vr_theta_avg[k,i])

		Mdotmax = np.max(Mdot[i_mod:nR-i_mod])
		count_s = 0
		count_g = 0
		for k in range(i_mod,nR-i_mod+1):
			if(Mdot[k]==Mdotmax and count_s==0):
				r_shock  = r[k]
				i_shock  = k
				count_s +=1
			if(qdot_theta_avg[k,i]>0.0 and count_g==0):
				r_gain   = r[k]
				i_gain   = k
				count_g += 1

		R_shock[i] = r_shock
		R_gain[i]  = r_gain

	rho_t_avg    = np.zeros(shape=(nR,nTheta))
	vr_t_avg     = np.zeros(shape=(nR,nTheta))
	T_t_avg      = np.zeros(shape=(nR,nTheta))
	P_t_avg      = np.zeros(shape=(nR,nTheta))
	S_t_avg      = np.zeros(shape=(nR,nTheta))
	Ye_t_avg     = np.zeros(shape=(nR,nTheta))
	qdot_t_avg   = np.zeros(shape=(nR,nTheta))
	DELTAT       = Time[nT-1] - Time[0]
	print('(PlotEverything) Time averaging quantities...')
	for k in range(0,nR):
		for j in range(0,nTheta):
			for i in range(0,nT-1):
				deltaT           = Time[i+1] - Time[i]
				rho_t_avg[k,j]  += (1.0 / DELTAT) * 0.5 * deltaT * (rho_3D[k,j,i]  + rho_3D[k,j,i+1])
				vr_t_avg[k,j]   += (1.0 / DELTAT) * 0.5 * deltaT * (vr_3D[k,j,i]   + vr_3D[k,j,i+1])
				T_t_avg[k,j]    += (1.0 / DELTAT) * 0.5 * deltaT * (T_3D[k,j,i]    + T_3D[k,j,i+1])
				P_t_avg[k,j]    += (1.0 / DELTAT) * 0.5 * deltaT * (P_3D[k,j,i]    + P_3D[k,j,i+1])
				S_t_avg[k,j]    += (1.0 / DELTAT) * 0.5 * deltaT * (S_3D[k,j,i]    + S_3D[k,j,i+1])
				Ye_t_avg[k,j]   += (1.0 / DELTAT) * 0.5 * deltaT * (Ye_3D[k,j,i]   + Ye_3D[k,j,i+1])
				qdot_t_avg[k,j] += (1.0 / DELTAT) * 0.5 * deltaT * (qdot_3D[k,j,i] + qdot_3D[k,j,i+1])

	rho_avg  = np.zeros(nR)
	vr_avg   = np.zeros(nR)
	T_avg    = np.zeros(nR)
	P_avg    = np.zeros(nR)
	S_avg    = np.zeros(nR)
	Ye_avg   = np.zeros(nR)
	qdot_avg = np.zeros(nR)
	print('(PlotEverything) Angle averaging quantities...')
	for k in range(0,nR):
		for j in range(0,nTheta-1):
			deltaTheta   = theta[j+1] - theta[j]
			rho_avg[k]  += 1.0 / (4.0) * deltaTheta * (rho_t_avg[k,j]  * np.sin(theta[j]) + rho_t_avg[k,j+1]  * np.sin(theta[j+1]))
			vr_avg[k]   += 1.0 / (4.0) * deltaTheta * (vr_t_avg[k,j]   * np.sin(theta[j]) + vr_t_avg[k,j+1]   * np.sin(theta[j+1]))
			T_avg[k]    += 1.0 / (4.0) * deltaTheta * (T_t_avg[k,j]    * np.sin(theta[j]) + T_t_avg[k,j+1]    * np.sin(theta[j+1]))
			P_avg[k]    += 1.0 / (4.0) * deltaTheta * (P_t_avg[k,j]    * np.sin(theta[j]) + P_t_avg[k,j+1]    * np.sin(theta[j+1]))
			S_avg[k]    += 1.0 / (4.0) * deltaTheta * (S_t_avg[k,j]    * np.sin(theta[j]) + S_t_avg[k,j+1]    * np.sin(theta[j+1]))
			Ye_avg[k]   += 1.0 / (4.0) * deltaTheta * (Ye_t_avg[k,j]   * np.sin(theta[j]) + Ye_t_avg[k,j+1]   * np.sin(theta[j+1]))
			qdot_avg[k] += 1.0 / (4.0) * deltaTheta * (qdot_t_avg[k,j] * np.sin(theta[j]) + qdot_t_avg[k,j+1] * np.sin(theta[j+1]))


	return[rho_theta_avg, vr_theta_avg, vtheta_theta_avg, cs_theta_avg, T_theta_avg, P_theta_avg, S_theta_avg, qdot_theta_avg, Ye_theta_avg, R_shock, R_gain, rho_avg, vr_avg, T_avg, P_avg, S_avg, Ye_avg, qdot_avg, Time, r]

def BruntVaisaila(rho, vr, T, P, S, qdot, ye, time, r, i_mod):
	# Quantities are time- and angle-averaged
	kb      = 8.617333262145*pow(10,-11) # MeV K^-1
	mb      = 1.66053872801*pow(10,-24)  # g
	G       = 6.6743*pow(10,-8)          # cm^3 g^-1 s^-1
	M       = 1.4*2*pow(10,33)           # g
	mu      = G*M
	nR      = np.size(r)
	nT      = np.size(time)
	g       = np.zeros(nR)
	drhodP  = np.zeros(nR)
	dPdYe   = np.zeros(nR)
	dPdS    = np.zeros(nR)
	dSdr    = np.zeros(nR)
	dYedr   = np.zeros(nR)
	OmegaSq = np.zeros(nR)
	Omega_I = np.zeros(nR)
	Mdot    = np.zeros(nR)
	print('(BruntVaisaila) Calculating derivatives...')
	for i in range(0,nR):
		g[i] = mu / (r[i]**2)
		Mdot[i] = 4.0 * math.pi * r[i]**2 * rho[i] * np.abs(vr[i])
		if(i==0):
			drhodP[i]  = (rho[i+1] - rho[i]) / (P[i+1]  - P[i])
			dPdYe[i]   = (P[i+1]   - P[i])   / (ye[i+1] - ye[i])
			dPdS[i]    = (P[i+1]   - P[i])   / (S[i+1]  - S[i])
			dSdr[i]    = (S[i+1]   - S[i])   / (r[i+1]  - r[i])
			dYedr[i]   = (ye[i+1]  - ye[i])  / (r[i+1]  - r[i])
		if(i==int(nR-1)):
			drhodP[i]  = (rho[i] - rho[i-1]) / (P[i]  - P[i-1])
			dPdYe[i]   = (P[i]   - P[i-1])   / (ye[i] - ye[i-1])
			dPdS[i]    = (P[i]   - P[i-1])   / (S[i]  - S[i-1])
			dSdr[i]    = (S[i]   - S[i-1])   / (r[i]  - r[i-1])
			dYedr[i]   = (ye[i]  - ye[i-1])  / (r[i]  - r[i-1])
		else:
			drhodP[i]  = (rho[i+1] - rho[i-1]) / (P[i+1]  - P[i-1])
			dPdYe[i]   = (P[i+1]   - P[i-1])   / (ye[i+1] - ye[i-1])
			dPdS[i]    = (P[i+1]   - P[i-1])   / (S[i+1]  - S[i-1])
			dSdr[i]    = (S[i+1]   - S[i-1])   / (r[i+1]  - r[i-1])
			dYedr[i]   = (ye[i+1]  - ye[i-1])  / (r[i+1]  - r[i-1])

			OmegaSq[i] = g[i] / rho[i] * drhodP[i] * (dPdS[i] * dSdr[i] + dPdYe[i] * dYedr[i])
			Omega_I[i] = np.imag(cmath.sqrt(OmegaSq[i]))

	counter = 0
	for i in range(0,nR):
		if(qdot[i]>0.0 and counter==0):
			rgain   = r[i]
			i_rgain = int(i)
			counter += 1

	MdotMax = np.max(Mdot[i_mod:np.size(r)-i_mod])
	for i in range(i_mod,int(np.size(r))-i_mod):
		if(Mdot[i]==MdotMax):
			rshock   = r[i]
			i_rshock = int(i)

	Chi = 0.0
	for i in range(int(i_rgain),int(i_rshock)):
		deltaR  = r[i+1] - r[i]
		Chi    += 0.5 * deltaR * (Omega_I[i] / np.abs(vr[i]) + Omega_I[i+1] / np.abs(vr[i+1]))

	return[Chi]

def BruntVaisaila_Chi_of_t(rho_avg, vr_avg, T_avg, P_avg, S_avg, qdot_avg, ye_avg, time, r, imod):
	# Quantities are angle-averaged
	kb    = 8.617333262145*pow(10,-11) # MeV K^-1
	mb    = 1.66053872801*pow(10,-24)  # g
	G     = 6.6743*pow(10,-8)          # cm^3 g^-1 s^-1
	M     = 1.4*2*pow(10,33)           # g
	mu    = G*M

	g       = np.zeros(np.size(r))
	drhodP  = np.zeros(shape=(np.size(r),np.size(time)))
	dPdYe   = np.zeros(shape=(np.size(r),np.size(time)))
	dPdS    = np.zeros(shape=(np.size(r),np.size(time)))
	dSdr    = np.zeros(shape=(np.size(r),np.size(time)))
	dYedr   = np.zeros(shape=(np.size(r),np.size(time)))
	OmegaSq = np.zeros(shape=(np.size(r),np.size(time)))
	Omega_I = np.zeros(shape=(np.size(r),np.size(time)))

	print('(BruntVaisaila_Chi_of_t) Calculating derivatives...')
	for i in range(0,int(np.size(r))):
		g[i] = mu / (r[i]**2)
		for j in range(0,int(np.size(time))):
			if(i==0):
				drhodP[i,j]  = (rho_avg[i+1,j] - rho_avg[i,j]) / (P_avg[i+1,j]  - P_avg[i,j])
				dPdYe[i,j]   = (P_avg[i+1,j]   - P_avg[i,j])   / (ye_avg[i+1,j] - ye_avg[i,j])
				dPdS[i,j]    = (P_avg[i+1,j]   - P_avg[i,j])   / (S_avg[i+1,j]  - S_avg[i,j])
				dSdr[i,j]    = (S_avg[i+1,j]   - S_avg[i,j])   / (r[i+1]        - r[i])
				dYedr[i,j]   = (ye_avg[i+1,j]  - ye_avg[i,j])  / (r[i+1]        - r[i])
			if(i==int(np.size(r))-1):
				drhodP[i,j]  = (rho_avg[i,j] - rho_avg[i-1,j]) / (P_avg[i,j]  - P_avg[i-1,j])
				dPdYe[i,j]   = (P_avg[i,j]   - P_avg[i-1,j])   / (ye_avg[i,j] - ye_avg[i-1,j])
				dPdS[i,j]    = (P_avg[i,j]   - P_avg[i-1,j])   / (S_avg[i,j]  - S_avg[i-1,j])
				dSdr[i,j]    = (S_avg[i,j]   - S_avg[i-1,j])   / (r[i]        - r[i-1])
				dYedr[i,j]   = (ye_avg[i,j]  - ye_avg[i-1,j])  / (r[i]        - r[i-1])
			else:
				drhodP[i,j]  = (rho_avg[i+1,j] - rho_avg[i-1,j]) / (P_avg[i+1,j]  - P_avg[i-1,j])
				dPdYe[i,j]   = (P_avg[i+1,j]   - P_avg[i-1,j])   / (ye_avg[i+1,j] - ye_avg[i-1,j])
				dPdS[i,j]    = (P_avg[i+1,j]   - P_avg[i-1,j])   / (S_avg[i+1,j]  - S_avg[i-1,j])
				dSdr[i,j]    = (S_avg[i+1,j]   - S_avg[i-1,j])   / (r[i+1]        - r[i-1])
				dYedr[i,j]   = (ye_avg[i+1,j]  - ye_avg[i-1,j])  / (r[i+1]        - r[i-1])


			OmegaSq[i,j] = g[i] / rho_avg[i,j] * drhodP[i,j] * (dPdS[i,j]*dSdr[i,j] + dPdYe[i,j]*dYedr[i,j])
			Omega_I[i,j] = np.imag(cmath.sqrt(OmegaSq[i,j]))

	rgain   = np.zeros(np.size(time))
	Mdot    = np.zeros(shape=(np.size(r),np.size(time)))
	i_rgain = np.zeros(np.size(time))
	for j in range(0,int(np.size(time))):
		counter = 0
		for i in range(0,int(np.size(r))):
			Mdot[i,j] = 4.0 * math.pi * r[i]**2 * rho_avg[i,j] * np.abs(vr_avg[i,j])
			if(qdot_avg[i,j]>0.0 and counter==0):
				rgain[j] = r[i]
				i_rgain[j] = int(i)
				counter += 1

	MdotMax  = np.zeros(np.size(time))
	rshock   = np.zeros(np.size(time))
	i_rshock = np.zeros(np.size(time))
	for j in range(0,int(np.size(time))):
		MdotMax[j] = np.max(Mdot[imod:np.size(r)-imod,j])
		for i in range(imod,int(np.size(r))-imod):
			if(Mdot[i,j]==MdotMax[j]):
				rshock[j]   = r[i]
				i_rshock[j] = int(i)

	Chi = np.zeros(np.size(time))
	for j in range(0,int(np.size(time))):
		for i in range(int(i_rgain[j]),int(i_rshock[j])):
			deltaR  = r[i+1] - r[i]
			Chi[j] += 0.5 * deltaR * (Omega_I[i,j] / np.abs(vr_avg[i,j]) + Omega_I[i+1,j] / np.abs(vr_avg[i+1,j]))


	return[Chi]

def ReynoldsStressTensor(rho, vr, vtheta, cs, T, R):
	# rho, vr, and vtheta are angle-averaged
	# Domains
	nT        = int(np.size(T))
	nR        = int(np.size(R))
	DELTA_t   = T[nT-1]-T[0]
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
		cssq_tavg  = 0.0
		v_esc      = np.sqrt(2.0 * mu / r[k])
		# Time integration
		for i in range(0,nT-1):
			delta_t     = T[i+1] - T[i]
			rho_bar[k] += (1.0 / DELTA_t) * 0.5 * delta_t * (rho[k,i] + rho[k,i+1])
			v0v0_x_rho += (1.0 / DELTA_t) * 0.5 * delta_t * (rho[k,i] * vr[k,i]     * vr[k,i]     + rho[k,i+1] * vr[k,i+1]     * vr[k,i+1])
			v0v1_x_rho += (1.0 / DELTA_t) * 0.5 * delta_t * (rho[k,i] * vr[k,i]     * vtheta[k,i] + rho[k,i+1] * vr[k,i+1]     * vtheta[k,i+1])
			v1v1_x_rho += (1.0 / DELTA_t) * 0.5 * delta_t * (rho[k,i] * vtheta[k,i] * vtheta[k,i] + rho[k,i+1] * vtheta[k,i+1] * vtheta[k,i+1])
			v0_x_rho   += (1.0 / DELTA_t) * 0.5 * delta_t * (rho[k,i] * vr[k,i]                   + rho[k,i+1] * vr[k,i+1])
			v1_x_rho   += (1.0 / DELTA_t) * 0.5 * delta_t * (rho[k,i] * vtheta[k,i]               + rho[k,i+1] * vtheta[k,i+1])
			cssq_tavg  += (1.0 / DELTA_t) * 0.5 * delta_t * (cs[k,i]**2                           + cs[k,i+1]**2)

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
		K[k]       = (1.0 + alpha[k]) * R00[k] / rho_bar[k] # Kinetic energy, cm^2/s^2
		xi_th[k]   = cssq_tavg / v_esc**2
		xi_turb[k] = 2.0 * K[k] / v_esc**2
		xi_eff[k]  = xi_th[k] + 0.5 * xi_turb[k]


	return[R00, R01, R10, R11, R00_g, R01_g, R10_g, R11_g, v0v0, v0v1, v1v1, v0, v1, rho_bar, alpha, K, xi_th, xi_turb, xi_eff]

def AdiabaticIndex(r, t, rho, P, qdot):
	# Calculate effective adiabatic index gamma(r) using time- and angle-averaged quantities
	# Quantities passed in are angle-averaged
	nT      = np.size(t)
	nR      = np.size(r)
	gamma   = np.zeros(nR)
	rho_t   = np.zeros(nR)
	P_t     = np.zeros(nR)
	qdot_t  = np.zeros(nR)
	DELTA_t = t[nT-1] - t[0]

	# take time-averages of rho and P, which are already angle-averaged
	print('(AdiabaticIndex) Calculating Gamma...')
	count = 0
	for k in range(0,nR):
		for i in range(0,nT-1):
			delta_t    = t[i+1] - t[i]
			rho_t[k]  += (1.0 / DELTA_t) * 0.5 * delta_t * (rho[k,i]  + rho[k,i+1])
			P_t[k]    += (1.0 / DELTA_t) * 0.5 * delta_t * (P[k,i]    + P[k,i+1])
			qdot_t[k] += (1.0 / DELTA_t) * 0.5 * delta_t * (qdot[k,i] + qdot[k,i+1])
		# Find time- and angle-averaged gain radius + index
		if(qdot_t[k]>0.0 and count==0):
			r_gain = r[k]
			i_gain = k
			count += 1


	# Calculate Gamma(r) = dln(P)/dln(rho)|_s
	for k in range(0,nR):
		if(k==0):
			gamma[k] = (np.log(P_t[k+1]) - np.log(P_t[k]))   / (np.log(rho_t[k+1]) - np.log(rho_t[k]))
		elif(k==nR-1):
			gamma[k] = (np.log(P_t[k])   - np.log(P_t[k-1])) / (np.log(rho_t[k])   - np.log(rho_t[k-1]))
		else:
			gamma[k] = (np.log(P_t[k+1]) - np.log(P_t[k-1])) / (np.log(rho_t[k+1]) - np.log(rho_t[k-1]))

	return[gamma, r_gain, i_gain]

def MachNumber(vr, Cs, R, time):
	nT   = np.size(time)
	nR   = np.size(R)
	mach = np.zeros(nT)
	vr_o = np.zeros(nT)
	cs_o = np.zeros(nT)
	for i in range(0,nT):
		mach[i] = np.abs(vr[nR-1,i] / Cs[nR-1,i])
		vr_o[i] = np.abs(vr[nR-1,i])
		cs_o[i] = np.abs(Cs[nR-1,i])

	return [mach, vr_o, cs_o]

[rho_theta_A, vr_theta_A, vtheta_theta_A, cs_theta_A, T_theta_A, P_theta_A, S_theta_A, qdot_theta_A, Ye_theta_A, Rshock, Rgain, rho_theta_t_A, vr_theta_t_A, T_theta_t_A, P_theta_t_A, S_theta_t_A, Ye_theta_t_A, qdot_theta_t_A, t, R] = ReadData(nmin,nmax,minTime,maxTime,DIR_name1,DIR_name2,radRes,thetRes,iModifier)

# --- FUNCTIONS ---

[chi]      = BruntVaisaila(rho_theta_t_A, vr_theta_t_A, T_theta_t_A, P_theta_t_A, S_theta_t_A, qdot_theta_t_A, Ye_theta_t_A, t, R, iModifier)

[chi_of_t] = BruntVaisaila_Chi_of_t(rho_theta_A, vr_theta_A, T_theta_A, P_theta_A, S_theta_A, qdot_theta_A, Ye_theta_A, t, R, iModifier)

[r00, r01, r10, r11, r00_g, r01_g, r10_g, r11_g, V0V0, V0V1, V1V1, V0, V1, Rhobar, Alpha, KE, Xi_th, Xi_turb, Xi_eff]= ReynoldsStressTensor(rho_theta_A, vr_theta_A, vtheta_theta_A, cs_theta_A, t, R)

[Gamma, RGain, iGain] = AdiabaticIndex(R, t, rho_theta_A, P_theta_A, qdot_theta_A)

[Mach, VrO, CsO] = MachNumber(vr_theta_A, cs_theta_A, R, Time)

# --- <Chi>_{t,theta} ---
print('---Chi derived from time- and angle-averaged quantities (Iwakami et al. 2014)---')
print('Chi = ' + str(chi))

# --- Chi(t) ---
name='Chi_of_t'
plt.figure(0)
plt.plot(Time,chi_of_t)
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$\chi$')
#plt.show()
plt.savefig(name+'.png')

# --- Mach number ---
name='MachNumber'
plt.figure(1)
plt.plot(Time,Mach)
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$\left|v_{r}/c_{s}\right|$')
#plt.show()
plt.savefig(name+'.png')

# --- Rshock and Rgain
DELTA_t = Time[nT-1]-Time[0]
R_shock_avg = 0.0
R_gain_avg  = 0.0
for i in range(nmin,nmax):
	delta_t      = Time[i+1] - Time[i]
	R_shock_avg += (1.0 / DELTA_t) * 0.5 * delta_t * (Rshock[i] + Rshock[i+1])
	R_gain_avg  += (1.0 / DELTA_t) * 0.5 * delta_t * (Rgain[i]  + Rgain[i+1])

R_shock_avg_ARRAY = np.zeros(nT)
R_gain_avg_ARRAY  = np.zeros(nT)
for i in range(nmin,nmax+1):
	R_shock_avg_ARRAY[i] = R_shock_avg
	R_gain_avg_ARRAY[i]  = R_gain_avg

name='Rshock_Rgain_2D'
plt.figure(2)
plt.plot(Time,Rshock/1.0e5,color='black',linewidth=2,label=r'$R_{\mathrm{shock}}$')
plt.plot(Time,R_shock_avg_ARRAY/1.0e5,color='black',linestyle='dashed',label=r'$\left<r_{\mathrm{shock}}\right>_{t}$ = '+'{:.1e}'.format(R_shock_avg/1.0e5) + ' km')
plt.plot(Time,Rgain/1.0e5,color='red',linewidth=2,label=r'$R_{\mathrm{gain}}$')
plt.plot(Time,R_gain_avg_ARRAY/1.0e5,color='red',linestyle='dashed',label=r'$\left<r_{\mathrm{gain}}\right>_{t}$ = '+'{:.1e}'.format(R_gain_avg/1.0e5) + ' km')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$[km]$')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend()
#plt.show()
plt.savefig(name+'.png')

# --- Reynolds Stress Stuff ---
critical = np.zeros(np.size(R))
for i in range(0,np.size(R)):
	critical[i] = 3.0 / 16.0 * Gamma[iGain]

name1='Rij_2D'
plt.figure(3)
plt.plot(r/1.0e5,np.abs(r00),label=r'$R00$ ($r$-$r$)')
plt.plot(r/1.0e5,np.abs(r01),label=r'$R01$ ($r$-$\theta$)')
plt.plot(r/1.0e5,np.abs(r10),label=r'$R10$ ($\theta$-$r$)')
plt.plot(r/1.0e5,np.abs(r11),label=r'$R11$ ($\theta$-$\theta$)')
plt.xlabel(r'$R$ [km]')
plt.ylabel(r'$R_{i,j}$ [erg/cm$^{3}$]')
plt.xlim(minR/1.0e5, maxR/1.0e5)
plt.legend()
plt.savefig(name1+'.png')
#plt.show()

name2='Rij_by_rho_2D'
plt.figure(4)
plt.plot(r/1.0e5,np.abs(r00_g),label=r'$R00$ ($r$-$r$)')
plt.plot(r/1.0e5,np.abs(r01_g),label=r'$R01$ ($r$-$\theta$)')
plt.plot(r/1.0e5,np.abs(r10_g),label=r'$R10$ ($\theta$-$r$)')
plt.plot(r/1.0e5,np.abs(r11_g),label=r'$R11$ ($\theta$-$\theta$)')
plt.xlabel(r'$R$ [km]')
plt.ylabel(r'${R_{i,j}}/{\bar{\rho}}$ [erg/g]')
plt.xlim(minR/1.0e5, maxR/1.0e5)
plt.legend()
plt.savefig(name2+'.png')
#plt.show()

name3='VV_terms_2D'
plt.figure(5)
plt.plot(r/1.0e5,V0V0,label=r'$\overline{v_{r}v_{r}}/\bar{\rho}$')
plt.plot(r/1.0e5,V0V1,label=r'$\overline{v_{r}v_{\theta}}/\bar{\rho}$')
plt.plot(r/1.0e5,V1V1,label=r'$\overline{v_{\theta}v_{\theta}}/\bar{\rho}$')
plt.plot(r/1.0e5,V0*V0,label=r'$\bar{v_{r}}\bar{v_{r}}/(\bar{\rho})^2$')
plt.plot(r/1.0e5,V0*V1,label=r'$\bar{v_{r}}\bar{v_{\theta}}/(\bar{\rho})^2$')
plt.plot(r/1.0e5,V1*V1,label=r'$\bar{v_{\theta}}\bar{v_{\theta}}/(\bar{\rho})^2$')
plt.xlabel(r'$R$ [km]')
plt.xlim(minR/1.0e5, maxR/1.0e5)
plt.legend()
plt.ylabel(r'[erg g$^{-1}$]')
plt.savefig(name3+'.png')
#plt.show()

name4='rhobar_2D'
plt.figure(6)
plt.plot(r/1.0e5,Rhobar)
plt.xlabel(r'$R$ [km]')
plt.ylabel(r'$\bar{\rho}$ [g cm$^{-3}$]')
plt.xscale('log')
plt.yscale('log')
plt.xlim(minR/1.0e5, maxR/1.0e5)
plt.savefig(name4+'.png')
#plt.show()

name5='xi_terms'
plt.figure(7)
plt.plot(r/1.0e5,Xi_th,label=r'$\xi_{\mathrm{thermal}}$')
plt.plot(r/1.0e5,Xi_turb,label=r'$\xi_{\mathrm{turb}}$')
plt.plot(r/1.0e5,Xi_eff,label=r'$\xi_{\mathrm{eff}}$')
plt.plot(r/1.0e5,critical,color='black',linestyle='dashed',label=r'$3/16 \gamma\left(r_{\mathrm{gain}}\right)$')
plt.legend()
plt.xlim(minR/1.0e5, maxR/1.0e5)
plt.xlabel(r'$R$ [km]')
plt.savefig(name5+'.png')
#plt.show()

name6='alpha'
plt.figure(8)
plt.plot(r/1.0e5,Alpha)
plt.xlabel(r'$R$ [km]')
plt.ylabel(r'$\alpha$')
plt.xlim(minR/1.0e5, maxR/1.0e5)
plt.savefig(name6+'.png')
#plt.show()

name7='K_energy'
plt.figure(9)
plt.plot(r/1.0e5,KE)
plt.xlabel(r'$R$ [km]')
plt.ylabel(r'$K$ [erg/g]')
plt.xlim(minR/1.0e5, maxR/1.0e5)
plt.savefig(name7+'.png')
#plt.show()

name8='Gamma'
plt.figure(10)
plt.plot(r/1.0e5,Gamma)
plt.xlabel(r'$R$ [km]')
plt.ylabel(r'$\gamma$')
plt.xlim(minR/1.0e5, maxR/1.0e5)
plt.savefig(name8+'.png')
#plt.show()

# --- vr and cs plots at rmax ---
name='VrCs'
plt.figure(11)
plt.plot(Time,VrO,color='black',label=r'$v_{r}$')
plt.plot(Time,CsO,color='red',label=r'$C_{s}$')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'[cm s$^{-1}$]')
plt.legend()
#plt.show()
plt.savefig(name+'.png')

