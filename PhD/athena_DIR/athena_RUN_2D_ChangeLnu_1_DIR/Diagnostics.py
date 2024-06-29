import matplotlib
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import sys
import h5py
import math
import os
sys.path.insert(0, '/home/dpochik/mhd-winds-new-EOS/vis/python')
import athena_read
plt.rcParams.update({'font.size': 10})
from matplotlib import ticker, cm
from scipy.interpolate import interp2d
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cmath
try:
    import helmeos
except ImportError:
    try:
        from .. import helmeos
    except ImportError:
        import sys
        sys.path.append(os.path.dirname(os.path.abspath(__file__)))
        import helmeos

DIR_name = str(sys.argv[1])
gif_name = str(sys.argv[2])
fig, [(ax1, ax2, ax3), (ax4, ax5, ax6)] = plt.subplots(2,3)

n1   = DIR_name+"/accretion.prim.00000.athdf"
data = athena_read.athdf(n1)
r    = data['x1v']
nR   = np.size(r)

kmin    = int(sys.argv[3])
k       = kmin
kmax    = int(sys.argv[4])
iFrames = int(kmax - k + 1)
tmax    = float(sys.argv[5])
time    = np.linspace(0,tmax,iFrames).T
rmin    = float(sys.argv[6])/1.0e5
rmax    = float(sys.argv[7])/1.0e5
Mdot    = float(sys.argv[8])
MachE   = float(sys.argv[9])
lcolor  = 'b'
kcolor  = 'black'
iMod    = int(3)
Mdot_U  = 1.3 * Mdot
Mdot_L  = 0.7 * Mdot

line1,  = ax1.loglog([],   [], lw=1.75, color=lcolor)
line2,  = ax2.loglog([],   [], lw=1.75, color=lcolor)
line3,  = ax3.loglog([],   [], lw=1.75, color=lcolor)
line4,  = ax4.loglog([],   [], lw=1.75, color=lcolor)
line5,  = ax5.loglog([],   [], lw=1.75, color=lcolor)
line5a, = ax5.loglog([],   [], lw=1.75, color=kcolor, linestyle='dashed')
line6,  = ax6.semilogx([], [], lw=1.75, color=lcolor)

kb     = 8.61733326e-11  # MeV/K
Na     = 6.02214e23      # Avogadro's constant
mb     = 1.0 / Na        # Baryon mass in g
abar   = 1.0
kb_erg = 1.380649e-16    # Boltzmann constant in erg/K

# Density
ax1.set_xlim(rmin,rmax)
ax1.set_ylim(4.0e6,1.0e12)
ax1.set_xlabel(r'$r$ [km]')
ax1.set_ylabel(r'$\rho$ [g/cm$^{3}$]')

# Temperature
ax2.set_xlim(rmin,rmax)
ax2.set_ylim(0.1,5.0)
ax2.set_xlabel(r'$r$ [km]')
ax2.set_ylabel(r'$T$ [MeV]')

# Pressure
ax3.set_xlim(rmin,rmax)
ax3.set_ylim(5.0e24,1.0e30)
ax3.set_xlabel(r'$r$ [km]')
ax3.set_ylabel(r'$P$ [erg/cm$^{3}$]')

# Velocity
ax4.set_xlim(rmin,rmax)
ax4.set_ylim(2.0e4,8.0e9)
ax4.set_xlabel(r'$r$ [km]')
ax4.set_ylabel(r'$\left|v_{r}\right|$ [cm/s$^{1}$]')

# Mdot
ax5.set_xlim(rmin,rmax)
ax5.set_ylim(Mdot_L,Mdot_U)
ax5.set_xlabel(r'$r$ [km]')
ax5.set_ylabel(r'$\dot{M}$ [M$_{\odot}$ s$^{-1}$]')

# Total Heating
ax6.set_xlim(rmin,rmax)
ax6.set_ylim(-20.0,3.0)
ax6.set_xlabel(r'$r$ [km]')
ax6.set_ylabel(r'$\dot{q}$ [$10^{21}$ erg/g/cm]')

lines = [line1, line2, line3, line4, line5, line5a, line6]
arr   = []

while(int(k)<=int(kmax)):
	if(int(k)%1==0):
		arr.append(int(k))
	k=int(k)+1

rho_avg  = np.zeros(shape=(nR,int(iFrames)))
T_avg    = np.zeros(shape=(nR,int(iFrames)))
P_avg    = np.zeros(shape=(nR,int(iFrames)))
vr_avg   = np.zeros(shape=(nR,int(iFrames)))
cs_avg   = np.zeros(shape=(nR,int(iFrames)))
Ye_avg   = np.zeros(shape=(nR,int(iFrames)))
qdot_avg = np.zeros(shape=(nR,int(iFrames)))
S_avg    = np.zeros(shape=(nR,int(iFrames)))
Mdot_avg = np.zeros(shape=(nR,int(iFrames)))
Mdot_EXP = np.zeros(nR)
for i in range(int(kmin),int(kmax+1)):
	if(i<10):
		n1=DIR_name+"/accretion.prim.0000"+str(i)+".athdf"
		n2=DIR_name+"/accretion.uov.0000"+str(i)+".athdf"
	if(i>=10 and i<100):
		n1=DIR_name+"/accretion.prim.000"+str(i)+".athdf"
		n2=DIR_name+"/accretion.uov.000"+str(i)+".athdf"
	if(i>=100 and i<1000):
		n1=DIR_name+"/accretion.prim.00"+str(i)+".athdf"
		n2=DIR_name+"/accretion.uov.00"+str(i)+".athdf"

	print('(DataRead) i = ' + str(i))
	data1  = athena_read.athdf(n1)
	data2  = athena_read.athdf(n2)
	rho    = (data1['rho'])[0]
	r      = data1['x1v']
	theta  = data1['x2v']
	p      = (data1['press'])[0]
	v      = (data1['vel1'])[0]
	ye     = (data1['r0'])[0]
	T      = (data2['dt1'])[0]
	qdot   = (data2['dt2'])[0]
	cs     = (data2['dt3'])[0]
	nTheta = np.size(theta)
	Ssum   = np.zeros(shape=(nTheta,nR,iFrames))
	fn     = os.path.join(os.path.dirname(os.path.abspath(__file__)), "helm_table.dat")
	by_filename = helmeos.HelmTable(fn=fn, temp_n=201, dens_n=541)
	EOSData     = by_filename.eos_DT(rho, T, abar, ye)
	Ssum[:,:,i] = EOSData['stot'] * (mb / kb_erg) # kb/baryon

	for k in range(0,nR): #i
		for l in range(0,nTheta-1): #j
			dtheta       = theta[l+1] - theta[l]
			rho_avg[k,i]  += (1.0 / 4.0) * dtheta * (rho[l,k]    * np.sin(theta[l]) + rho[l+1,k]    * np.sin(theta[l+1]))
			T_avg[k,i]    += (1.0 / 4.0) * dtheta * (T[l,k]      * np.sin(theta[l]) + T[l+1,k]      * np.sin(theta[l+1]))
			P_avg[k,i]    += (1.0 / 4.0) * dtheta * (p[l,k]      * np.sin(theta[l]) + p[l+1,k]      * np.sin(theta[l+1]))
			vr_avg[k,i]   += (1.0 / 4.0) * dtheta * (v[l,k]      * np.sin(theta[l]) + v[l+1,k]      * np.sin(theta[l+1]))
			cs_avg[k,i]   += (1.0 / 4.0) * dtheta * (cs[l,k]     * np.sin(theta[l]) + cs[l+1,k]     * np.sin(theta[l+1]))
			Ye_avg[k,i]   += (1.0 / 4.0) * dtheta * (ye[l,k]     * np.sin(theta[l]) + ye[l+1,k]     * np.sin(theta[l+1]))
			qdot_avg[k,i] += (1.0 / 4.0) * dtheta * (qdot[l,k]   * np.sin(theta[l]) + qdot[l+1,k]   * np.sin(theta[l+1]))
			S_avg[k,i]    += (1.0 / 4.0) * dtheta * (Ssum[l,k,i] * np.sin(theta[l]) + Ssum[l+1,k,i] * np.sin(theta[l+1]))

		Mdot_avg[k,i] = 4.0 * np.pi * r[k]**2 * rho_avg[k,i] * np.abs(vr_avg[k,i]) / 2.0e33
		Mdot_EXP[k]   = Mdot

print(np.size(rho_avg))
print('rho_avg[0,1] = ' + str(rho_avg[0,1]))

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

	DEL_T = time[np.size(time)-1] - time[0]
	Chi_AVG = 0.0
	for j in range(0,int(np.size(time)-1)):
		del_T    = time[j+1] - time[j]
		Chi_AVG += 1.0 / DEL_T * 0.5 * del_T * (Chi[j] + Chi[j+1])


	return[Chi, Chi_AVG]

def Mach(r, vr, cs):
	nR      = np.size(r)
	nT      = np.size(time)
	MachNum = np.zeros(nT)
	for i in range(0,nT):
		MachNum[i] = np.abs(vr[nR-1,i] / cs[nR-1,i]) / MachE

	return [MachNum]

def init():
	for i in range(0,len(lines)):
		lines[i].set_data([],[])
	return lines

def animate(i):
	print('(animate) i = ' + str(i))

	tformat  = '{:.2f}'.format
	plt.suptitle('t = '+tformat(i*float(tmax)/float(kmax) * 1.0e3)+' ms')
	plt.tight_layout(pad=0.1)
	fig.subplots_adjust(top=0.88)

	for j in range(0,len(lines)):
		if(j==0):
			lines[j].set_data(r/1.0e5,rho_avg[:,i])
		if(j==1):
			lines[j].set_data(r/1.0e5,T_avg[:,i]*kb)
		if(j==2):
			lines[j].set_data(r/1.0e5,P_avg[:,i])
		if(j==3):
			lines[j].set_data(r/1.0e5,np.abs(vr_avg[:,i]))
		if(j==4):
			lines[j].set_data(r/1.0e5,Mdot_avg[:,i])
		if(j==5):
			lines[j].set_data(r/1.0e5,Mdot_EXP)
		if(j==6):
			lines[j].set_data(r/1.0e5,qdot_avg[:,i]/1.0e21)
	return lines

anim = FuncAnimation(fig, animate, init_func=init,
                               frames=arr, interval=100, blit=True,repeat_delay=200)


anim.save(gif_name+'.gif', writer='imagemagick')

[Chi_t, Chi_T_AVG] = BruntVaisaila_Chi_of_t(rho_avg, vr_avg, T_avg, P_avg, S_avg, qdot_avg, Ye_avg, time, r, iMod)
[MachNumber]       = Mach(r, vr_avg, cs_avg)
Chi_AVG_Array      = np.zeros(iFrames)
name = 'Chi_of_t'
for n in range(kmin,kmax+1):
	Chi_AVG_Array[n] = Chi_T_AVG
plt.figure(2)
plt.plot(time,Chi_t,linewidth=2,color='black')
plt.plot(time,Chi_AVG_Array,linewidth=2,color='red',linestyle='dashed',label=r'$\left< \chi \right>_{t}$')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$\chi$')
plt.title(r'$\left< \chi \right>_{t} \sim $' + str(Chi_T_AVG))
plt.legend()
plt.savefig(name+'.png')

plt.figure(3)
name='Mach'
plt.plot(time,MachNumber,color='black')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'${\mathcal{M}} / {\mathcal{M}_{\mathrm{expected}}}$ at $r_{\mathrm{max}}$')
plt.title(r'$\mathcal{M} = $' + str(MachE))
plt.savefig(name+'.png')
