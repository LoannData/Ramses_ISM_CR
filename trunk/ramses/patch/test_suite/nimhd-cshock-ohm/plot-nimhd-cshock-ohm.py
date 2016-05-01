from pylab import *

fig = matplotlib.pyplot.figure()
ratio = 1.3
sizex = 10.0
fig.set_size_inches(sizex,ratio*sizex)

# Read in data
data1 = loadtxt('data1.dat')
data2 = loadtxt('data2.dat')
data3 = loadtxt('data3.dat')
data4 = loadtxt('data4.dat')
data5 = loadtxt('data5.dat')
data6 = loadtxt('nimhd-cshock-ohm-ana.dat')

# Density
density = subplot(321)
density.plot(data6[:,0],data6[:,5],color='red')
density.plot(data1[:,0],data1[:,3],'o',color='black',markerfacecolor='none')
density.set_xlabel('Distance (cm)')
density.set_ylabel('Density (g/cm3)')

# Vx
vx = subplot(322)
vx.plot(data6[:,0],data6[:,3],color='red')
vx.plot(data2[:,0],data2[:,3],'o',color='black',markerfacecolor='none')
vx.set_xlabel('Distance (cm)')
vx.set_ylabel('Vx (cm/s)')

# Vy
vy = subplot(323)
vy.plot(data6[:,0],data6[:,4],color='red')
vy.plot(data3[:,0],data3[:,3],'o',color='black',markerfacecolor='none')
vy.set_xlabel('Distance (cm)')
vy.set_ylabel('Vy (cm/s)')

# By
by = subplot(324)
by.plot(data6[:,0],data6[:,1],color='red')
by.plot(data4[:,0],data4[:,3],'o',color='black',markerfacecolor='none')
by.set_xlabel('Distance (cm)')
by.set_ylabel('By')

# Pressure
pressure = subplot(325)
pressure.plot(data6[:,0],data6[:,2],color='red')
pressure.plot(data5[:,0],data5[:,3],'o',color='black',markerfacecolor='none')
pressure.set_xlabel('Distance (cm)')
pressure.set_ylabel('Pressure (g/cm/s2)')

savefig('nimhd-cshock-ohm.pdf',bbox_inches='tight')
