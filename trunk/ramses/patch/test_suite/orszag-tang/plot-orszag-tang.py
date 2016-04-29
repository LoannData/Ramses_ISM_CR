from pylab import *

fig = matplotlib.pyplot.figure()
ratio = 1.0
sizex = 18.0
fig.set_size_inches(sizex,ratio*sizex)

# Read in data
data = loadtxt('data1.dat')
rho  = data[:,2].reshape(256,256)

# Density
density = subplot(211)
imshow(rho,origin='lower',extent=[0.0,1.0,0.0,1.0],interpolation='None',cmap='jet',vmin=0.05,vmax=0.5)
xlabel('Distance x (cm)')
ylabel('Distance y (cm)')
colorbar(label='Density')

savefig('orszag-tang.pdf',bbox_inches='tight')
