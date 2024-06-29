#import h5py
#import matplotlib.pyplot as plt
import numpy as np
import math
import sys
pythonpath=str(sys.argv[8])
sys.path.insert(0, pythonpath)
import athena_read

def fermi(n,eta):
    if(n==0.0):
        fermi_f = np.log(1.0 + np.exp(eta))
    elif(n==1.0):
        a  = np.exp(-1.0 * abs(eta))
        s  = (eta**2) / 2.0 + 1.6449341
        ff = a - a**2 / 4.0 + pow(a,3) / 9.0 - pow(a,4) / 16.0 + pow(a,5) / 25.0 - pow(a,6) / 36.0 + pow(a,7) / 49.0
        if(eta < 0):
            fermi_f = ff
        elif(eta == 0):
            fermi_f = s - ff
        elif(eta > 0):
            fermi_f = s - ff
    elif(n==2.0):
        a  = np.exp(-1.0 * abs(eta))
        s  = pow(eta,3) / 3.0 + 3.2898681 * eta
        ff = 2.0 * (a - a**2/8.0 + pow(a,3)/27.0 - pow(a,4)/64.0 + pow(a,5)/125.0 - pow(a,6)/216.0)
        if(eta<0):
            fermi_f = ff
        elif(eta==0):
            fermi_f = s + ff
        elif(eta>0):
            fermi_f = s + ff
    elif(n==3.0):
        a  = np.exp(-1.0 * abs(eta))
        s  = pow(eta,4) / 4.0 + 4.9348022 * eta**2 + 11.351273
        ff = 6.0 * (a - a**2/16.0 + pow(a,3)/81.0 - pow(a,4)/256.0)
        if(eta<0):
            fermi_f = ff
        elif(eta==0):
            fermi_f = s - ff
        elif(eta>0):
            fermi_f = s - ff
    elif(n==4.0):
        a  = np.exp(-1.0 * abs(eta))
        s  = pow(eta,5) / 5.0 + 6.5797363 * pow(eta,3) + 45.457576 * eta
        ff = 24.0 * (a - a**2/32.0 + pow(a,3)/243.0)
        if(eta<0):
            fermi_f = ff
        elif(eta==0):
            fermi_f = s + ff
        elif(eta>0):
            fermi_f = s + ff
    elif(n==5.0):
        a  = np.exp(-1.0 * abs(eta))
        s  = pow(eta,6) / 6.0 + 8.2246703 * pow(eta,4) + 113.64394 * eta**2 + 236.53226
        ff = 120.0 * (a - a**2/64.0 + pow(a,3)/729.0)
        if(eta<0):
            fermi_f = ff;
        elif(eta==0):
            fermi_f = s - ff
        elif(eta>0):
            fermi_f = s - ff

    return fermi_f

def HeatingCooling_Scheck(T,R,Ye,Rho,Lnue,Lnueb,enue,enueb):
	kbol        = 1.3806504*pow(10,-16)   # erg / K
	c           = 2.99792458*pow(10,10)   # cm / 2
	hbar        = 1.05457266*pow(10,-27)  # erg s
	mn          = 1.6749286*pow(10,-24)   # g
	MeV2erg     = 1.6021773*pow(10,-6)
	melectron   = 9.1093897*pow(10,-28)   # g
	delta       = 1.2935*MeV2erg # erg
	alpha       = 1.254
	sigmaweak   = 1.76*pow(10,-44)*(1.0+3.0*(alpha)**2)/(4.0*(melectron*c*c)**2)   # weak cross section cm^2 / erg^2
	fermi2      = fermi(2.0,0.0)
	fermi4      = fermi(4.0,0.0)
	fermi5      = fermi(5.0,0.0)
	fermi3      = fermi(3.0,0.0)
	fermifactor = fermi2 * fermi5 / ( fermi3 * fermi4 )
	epsnu       = enue*MeV2erg   # neutrino energy in ergs
	epsnubar    = enueb*MeV2erg   # anineutrino energy in ergs
	Lnu         = Lnue*pow(10,51) # assuming luminosities are the same (Lnue=Lnueb)

	fnu = 0.5 * (1.0 + np.sqrt(1.0 - (R / R)**2));
	Tprime      = T * kbol / (hbar * c) # cm^-1
	Terg        = T * kbol # erg
	tfact       = pow(27.0 * Ye * Rho + \
			np.sqrt(3.0) * np.sqrt(4.0 * math.pi**2 * mn**2 * pow(Tprime,6) + \
			243.0 * Ye**2 * Rho**2),1.0/3.0) # g^(1/3) cm^-1
	eta_e       = pow(12.0 * math.pi, 1.0/3.0) * (pow(math.pi * mn, 1.0/3.0) * tfact**2 - \
			math.pi * Tprime**2 * mn * pow(12.0,1.0/3.0)) \
			/ (6.0 * pow(mn,2.0/3.0) * tfact * Tprime) # unitless
	sinkconst   = sigmaweak * 0.5 * (Lnu) / (4.0 * math.pi * R**2 * fnu) # erg^-1 s^-1
	sourceconst = sigmaweak * 0.5 * c * pow(Tprime,3.0) / (math.pi**2 * mn) # erg^-2 s^-1 g^-1
	q_abs       = sinkconst / mn * ((1.0 - Ye) * (epsnu**2 * fermifactor + \
			2.0 * delta * epsnu * np.sqrt(fermi(2,0) * fermi(4,0)) / fermi(3,0) + delta**2) + \
			Ye * (epsnubar**2 * fermifactor + \
			3.0 * delta * epsnubar * np.sqrt(fermi(2,0) * fermi(4,0)) / fermi(3,0) + 3.0*delta**2)) # erg s^-1 g^-1
	q_em        = sourceconst * (Ye * (pow(Terg,3.0) * fermi(5.0,eta_e - delta / Terg) + \
			2.0 * delta * Terg**2 * fermi(4.0,eta_e - delta / Terg) + \
			delta**2 * Terg * fermi(3.0,eta_e - delta / Terg)) + \
			(1.0 - Ye) * (pow(Terg,3.0) * fermi(5.0,-1.0*eta_e) + \
			3.0 * delta * Terg**2 * fermi(4.0,-1.0*eta_e) + \
			3.0 * delta**2 * Terg * fermi(3.0,-1.0*eta_e) + \
			pow(delta,3.0) * fermi(2.0,-1.0*eta_e))) # erg s^-1 g^-1
	q_tot       = q_abs - q_em

	return[q_em,q_abs,q_tot]

def NR(r,rho,T,ye,Lnue,Lnueb,enue,enueb):
	# --- NR parameters
	T0      = T
	Tnew    = T0
	DeltaT0 = 1.0e-3
	tol     = 1.0e-4
	mod     = 1.0e0
	maxC    = 20
	eps     = 1.0e-10

	# --- Physical constants
	kbol_MeV = 8.61733326e-11;               # Boltzmann constant in MeV / K

	count = 0
	while (count<=maxC):
		# --- Differential terms
		dT0 = DeltaT0 * T0
		dT  = (T0 + dT0) - T0;

		# --- Scale terms and original function
		[C0,H0,qdot0] = HeatingCooling_Scheck(T0,r,ye,rho,Lnue,Lnueb,enue,enueb)

		# --- Function with differential term T0+dT0
		[dC,dH,dqdot] = HeatingCooling_Scheck(T0+dT0,r,ye,rho,Lnue,Lnueb,enue,enueb)

		# --- Rescale functions
		f    = qdot0 / H0
		f_dT = dqdot / H0

		# --- Derivative terms
		dqdotdT = (f_dT - f) / (dT * kbol_MeV)

		# --- NR step
		T1 = T0 - mod * f * (1.0 / dqdotdT) / kbol_MeV
		if (np.abs(T1 - T0) / np.abs(T1)<tol):
			Tnew = T1
			[Cnew,Hnew,qdotNew] = HeatingCooling_Scheck(Tnew,r,ye,rho,Lnue,Lnueb,enue,enueb)
			print('(NR) Root found for T!')
			print('(NR) iterations: ' + str(count+1))
			print('(NR) qdot(Tnew) = ' + '{:.12e}'.format(qdotNew) + ' erg/s/g')
			break
		else:
			T0 = T1

		# --- Checking iterations
		if (count>=maxC):
			print('(NR) Maximum iterations reached. Returning original guess temperature.')
			Tnew = T
			break

		count += 1

	return Tnew

# Not using index = 0 because weird stuff happens with it
index_1 = 2
index_2 = 1

# inputs
vrIC      = float(sys.argv[3])
Lnue      = float(sys.argv[4])
Lnueb     = float(sys.argv[5])
enue      = float(sys.argv[6])
enueb     = float(sys.argv[7])
Mdot      = float(sys.argv[9])
rIC       = float(sys.argv[10])

FILE_name_f = "accretion.prim."+sys.argv[2]+".athdf"
FILE_name_g = "accretion.uov."+sys.argv[2]+".athdf"
IC_name     = 'Lnu_'+str(Lnu)+'e51_Mdot_'+str(Mdot)+'_RPNS_'+sys.argv[1]+'km_extrapolate.txt'
f_data = FILE_name_f
g_data = FILE_name_g

data_prim = athena_read.athdf(f_data)
data_uov  = athena_read.athdf(g_data)

r1        = data_prim['x1f']
vr1       = (data_prim['vel1'])[0]
P1        = (data_prim['press'])[0]
rho1      = (data_prim['rho'])[0]
T         = (data_uov['dt1'])[0]
elecfrac1 = (data_prim['r0'])[0]

rho  = rho1[0,:]
vr1  = vr1[0,:]
T    = T[0,:]
X    = r1[:]
P    = P1[0,:]
Ye   = elecfrac1[0,:]

rho_y1 = rho[index_1]
rho_y2 = rho[index_2]
T_y1   = T[index_1]
T_y2   = T[index_2]
Ye_y1  = Ye[index_1]
Ye_y2  = Ye[index_2]
v_y1   = vr1[index_1]
v_y2   = vr1[index_2]

x1     = X[index_1]
x2     = X[index_2]
x_f    = float(sys.argv[1])
xa     = x_f * pow(10,5)

rho_y = rho_y1 - (xa - x1) / (x1 - x2) * (rho_y2 - rho_y1)
Ye_y  = Ye_y1  - (xa - x1) / (x1 - x2) * (Ye_y2  - Ye_y1)
v_y   = (xa/rIC) * vrIC

# --- Use T_y2 as a guess temperature for NR
Tg    = T_y2
T_y   = NR(xa,rho_y,Tg,Ye_y,Lnue,Lnueb,enue,enueb)



print('rho(start)        = ' + '{:.6e}'.format(rho[0]) + " g/cm^3")
print('rho(extrapolated) = ' + '{:.6e}'.format(rho_y) + " g/cm^3")
print('\n')
print('T(start)          = ' + '{:.6e}'.format(T[0])+ " K")
print('T(extrapolated)   = ' + '{:.6e}'.format(T_y) + " K")
print('\n')
print('Ye(start)         = ' + '{:.6e}'.format(Ye[0]))
print('Ye(extrapolated)  = ' + '{:.6e}'.format(Ye_y))
print('\n')
print('v(start)          = ' + '{:.6e}'.format(vr1[0]) + " cm/s")
print('v(extrapolated)   = ' + '{:.6e}'.format(v_y) + " cm/s")

rho_array    = rho[1:np.size(rho)]
newarray_rho = np.insert(rho_array,0,rho_y,axis=0)

vr1_array    = vr1[1:np.size(vr1)]
newarray_vr1 = np.insert(vr1_array,0,v_y,axis=0)

T_array      = T[1:np.size(T)]
newarray_T   = np.insert(T_array,0,T_y,axis=0)

Ye_array     = Ye[1:np.size(Ye)]
newarray_Ye  = np.insert(Ye_array,0,Ye_y,axis=0)

X_array      = X[1:np.size(X)]
newarray_X   = np.insert(X_array,0,xa,axis=0)

c=[newarray_X,newarray_rho,newarray_vr1,newarray_T,newarray_Ye]
with open(IC_name, "w") as file:
    for x in zip(*c):
        file.write("{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\n".format(*x))


