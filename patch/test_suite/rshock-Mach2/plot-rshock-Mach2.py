from pylab import *

# Figure properties
fig = matplotlib.pyplot.figure()
ratio = 1.5
sizex = 8.0
fig.set_size_inches(sizex,ratio*sizex)

# Define constants
ar     = 7.56591469318689378e-015
kb     = 1.3806200e-16
mh     = 1.66e-24
mu_gas = 1.0
xshift = 1013.67

# Read in data
data   = loadtxt('data.dat')
amrlev = data[:, 0]
x      = data[:, 1] - xshift
rho    = data[:, 2]
u      = data[:, 3]
p      = data[:, 6]
er     = data[:,10]
Tr     = (er/ar)**0.25
Tg     = p*mu_gas*mh/(rho*kb)

# Read analytical solution
data_ana = loadtxt('rshock-Mach2-ana.dat')
x_ana    = data_ana[:,0]
rho_ana  = data_ana[:,2]
Tr_ana   = data_ana[:,4]
Tg_ana   = data_ana[:,3]

xmin = -1000.0
xmax =  1000.0

# Temperature
temperature = subplot(211)
temperature.plot(x,Tg,'o',color='black',markerfacecolor='none')
temperature.plot(x,Tr,'s',color='black',markerfacecolor='none')
temperature.plot(x_ana,Tr_ana,color='red')
temperature.plot(x_ana,Tg_ana,color='blue')
temperature.set_xlabel('Distance (cm)')
temperature.set_ylabel('Temperature (K)')
levels = temperature.twinx()
majorLocatorY = MultipleLocator(1.0)
levels.yaxis.set_major_locator(majorLocatorY)
levels.plot(x,amrlev,color='black',ls='dotted')
levels.set_ylabel('AMR Level')
temperature.set_xlim([xmin,xmax])

# Density
density = subplot(212)
density.plot(x,rho,'o',color='black',markerfacecolor='none')
density.plot(x_ana,rho_ana,color='red')
density.set_xlabel('Distance (cm)')
density.set_ylabel('Density (g/cm3)')
density.set_xlim([xmin,xmax])

savefig('rshock-Mach2.pdf',bbox_inches='tight')
