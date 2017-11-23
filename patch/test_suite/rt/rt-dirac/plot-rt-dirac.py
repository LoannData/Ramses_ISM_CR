import matplotlib as mpl
mpl.use('Agg')
from pylab import *

fig = matplotlib.pyplot.figure()
ratio = 1.0
sizex = 18.0
fig.set_size_inches(sizex,ratio*sizex)

pc=3.08e18
t=loadtxt('time.dat')

subplot(211)
# Density
data = loadtxt('data1.dat')
x    = data[:,0]/pc
z    = data[:,2]/pc
rho  = data[:,3].reshape(128,128)

imshow(np.log10(rho),origin='lower',interpolation='None',cmap='jet',extent=(min(x),max(x),min(z),max(z)))
xlabel('Distance x (pc)')
ylabel('Distance z (pc)')
title('t=%4.3f (Myr)'%t)
colorbar(label='Log(n)',shrink=0.5)

subplot(212)
# Pressure
data = loadtxt('data2.dat')
x    = data[:,0]/pc
z    = data[:,2]/pc
P    = data[:,3].reshape(128,128)

imshow(np.log10(P),origin='lower',interpolation='None',cmap='jet',extent=(min(x),max(x),min(z),max(z)))
xlabel('Distance x (pc)')
ylabel('Distance z (pc)')
colorbar(label='Log (P)',shrink=0.5)

# subplot(232)
# # Density
# data = loadtxt('data3.dat')
# y    = data[:,1]/pc
# z    = data[:,2]/pc
# rho  = data[:,3].reshape(128,128)

# imshow(np.log10(rho),origin='lower',interpolation='None',cmap='jet',extent=(min(y),max(y),min(z),max(z)))
# xlabel('Distance y (pc)')
# ylabel('Distance z (pc)')
# title('t=%4.3f (Myr)'%t)
# colorbar(label='Log(n)',shrink=0.5)

# subplot(235)
# # Pressure
# data = loadtxt('data4.dat')
# y    = data[:,1]/pc
# z    = data[:,2]/pc
# P    = data[:,3].reshape(128,128)

# imshow(np.log10(P),origin='lower',interpolation='None',cmap='jet',extent=(min(y),max(y),min(z),max(z)))
# xlabel('Distance y (pc)')
# ylabel('Distance z (pc)')
# colorbar(label='Log (P)',shrink=0.5)

# subplot(233)
# # Density
# data = loadtxt('data5.dat')
# x    = data[:,0]/pc
# y    = data[:,1]/pc
# rho  = data[:,3].reshape(128,128)

# imshow(np.log10(rho),origin='lower',interpolation='None',cmap='jet',extent=(min(x),max(x),min(y),max(y)))
# xlabel('Distance x (pc)')
# ylabel('Distance y (pc)')
# title('t=%4.3f (Myr)'%t)
# colorbar(label='Log(n)',shrink=0.5)

# subplot(236)
# # Pressure
# data = loadtxt('data6.dat')
# x    = data[:,0]/pc
# y    = data[:,1]/pc
# P    = data[:,3].reshape(128,128)

# imshow(np.log10(P),origin='lower',interpolation='None',cmap='jet',extent=(min(x),max(x),min(y),max(y)))
# xlabel('Distance x (pc)')
# ylabel('Distance y (pc)')
# colorbar(label='Log (P)',shrink=0.5)

savefig('rt-dirac.pdf',bbox_inches='tight')
