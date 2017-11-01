from pylab import *
import numpy as np

fig = matplotlib.pyplot.figure()
ratio = 1.3
sizex = 10.0
fig.set_size_inches(sizex,ratio*sizex)

# Read in data
data1  = loadtxt('data1.dat')
x0      = data1[:,0]
y0      = data1[:,1]
z0      = data1[:,2]
er0     = data1[:,17]

x=x0[::20]
y=y0[::20]
z=z0[::20]
er=er0[::20]

# Time
t = loadtxt('time.dat')
# Compute analytical solution
chi=1.0e10
x0=0.5
y0=0.5
z0=0.5
r=np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
E0=1.0e5
ana2 = E0/(8.0*(chi*np.pi*t)**1.5)
ana = ana2*np.exp(-r**2/(4.0*t*chi)) + 1.0

# Radiative energy
erad = subplot(311)
erad.semilogy(r,er,'.',color='black',label='Simulation')
erad.semilogy(r,ana,'.',color='red',label='Analytical')
erad.set_xlabel('Distance (cm)')
erad.set_ylabel('Radiative energy')
erad.legend()

# Relative error
error = subplot(312)
error.semilogy(r,abs(er-ana)/ana,'.',color='black')
error.set_xlabel('Distance (cm)')
error.set_ylabel('Error')

# 2d map
data2  = loadtxt('data2.dat')
x      = data2[:,0]
z      = data2[:,2]
er     = data2[:,3]

subplot(313)
imshow(np.log10(er.reshape(128,5,128)[:,4,:]),vmin=1,vmax=8,extent=(min(x),max(x),min(z),max(z)))
xlabel('Distance x (cm)')
ylabel('Distance z (cm)')
title('Log(er)')
colorbar()

savefig('dirac3d.pdf',bbox_inches='tight')
