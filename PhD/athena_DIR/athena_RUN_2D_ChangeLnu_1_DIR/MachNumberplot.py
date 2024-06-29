import matplotlib

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(0, '/home/dpochik/mhd-winds-new-EOS/vis/python')
import athena_read

DIR_name = str(sys.argv[1])
gif_name = str(sys.argv[2])
fig = plt.figure()
k=int(sys.argv[3])
kmax=int(sys.argv[4])
tmax=float(sys.argv[5])
rmin=float(sys.argv[6])
rmax=float(sys.argv[7])
MachExpected=float(sys.argv[8])
l1acolor='k'
l1bcolor='r'

ymin = MachExpected - 0.5*MachExpected
ymax = MachExpected + 0.5*MachExpected

ax     = plt.axes(xlim=(rmin/1.0e5,rmax/1.0e5), ylim=(ymin, ymax),xlabel=r"$R$ [km]",ylabel=r"Mach number")
line1a = ax.semilogx([],   [], lw=1.75, color=l1acolor, label=r'$V_{r}/C_{s}$ output')[0]
line1b = ax.semilogx([],   [], lw=1.75, color=l1bcolor, linestyle='dashed', label=r'$V_{r}/C_{s}$=' + str(MachExpected))[0]

lines = [line1a, line1b]
arr=[]

while(int(k)<=int(kmax)):
    if(int(k)%1==0):
        arr.append(int(k))
    k=int(k)+1

def init():
    for i in range(0,len(lines)):
        lines[i].set_data([],[])
    return lines

def animate(i):

    print(i)

    if(i<10):
        n1=DIR_name+"/accretion.prim.0000"+str(i)+".athdf"
        n2=DIR_name+"/accretion.uov.0000"+str(i)+".athdf"

    if(i>=10 and i<100):
        n1=DIR_name+"/accretion.prim.000"+str(i)+".athdf"
        n2=DIR_name+"/accretion.uov.000"+str(i)+".athdf"

    if(i>=100 and i<1000):
        n1=DIR_name+"/accretion.prim.00"+str(i)+".athdf"
        n2=DIR_name+"/accretion.uov.00"+str(i)+".athdf"

    kb = 8.61733326e-11 # MeV/K

    data1 = athena_read.athdf(n1)
    data2 = athena_read.athdf(n2)

    r     = data1['x1v']
    v     = (data1['vel1'])[0][0]
    cs    = (data2['dt3'])[0][0]
    Mach  = np.abs(v/cs)
    nR    = np.size(r)
    MachE = np.zeros(nR)
    for j in range(0,nR):
      MachE[j] = float(MachExpected)

    tformat  = '{:.2f}'.format
    plt.suptitle('t = '+tformat(i*float(tmax)/float(kmax) * 1.0e3)+' ms')
    plt.tight_layout(pad=0.1)
    fig.subplots_adjust(top=0.88)
    for i in range(0,len(lines)):
        if(i==0):
            lines[i].set_data(r/1.0e5,Mach)
        if(i==1):
            lines[i].set_data(r/1.0e5,MachE)
    legend = plt.legend(loc='upper right')
    return lines

anim = FuncAnimation(fig, animate, init_func=init,
                               frames=arr, interval=100, blit=True,repeat_delay=200)


anim.save(gif_name+'_MachNumber.gif', writer='imagemagick')
