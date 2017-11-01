from pylab import *

# Figure properties
fig = matplotlib.pyplot.figure()
ratio = 0.8
sizex = 10.0
fig.set_size_inches(sizex,ratio*sizex)

# Define constants
ar     = 7.56591469318689378e-015
kb     = 1.3806200e-16
mh     = 1.66e-24
mu_gas = 1.0
xshift = 4205.08

# Read in data
data   = loadtxt('data.dat')
amrlev = data[:, 0]
x      = data[:, 1] - xshift
rho    = data[:, 2]
p      = data[:, 6]
er     = data[:,10]
Tr     = (er/ar)**0.25
Tg     = p*mu_gas*mh/(rho*kb)

# Read analytical solution
data_ana = loadtxt('rshock-Mach5-ana.dat')
xTg_ana   = data_ana[:,2]*1.e+05
Tg_ana    = data_ana[:,3]/(8.57274781455567/855.7)
xTr_ana   = data_ana[:,0]
Tr_ana    = data_ana[:,1]

xmin = -3000.0
xmax =   500.0
ymin =   100.0
ymax =  1000.0

# Temperature
temperature = subplot(111)
temperature.plot(x,Tg,'o',color='black',markerfacecolor='none')
temperature.plot(x,Tr,'s',color='black',markerfacecolor='none')
temperature.plot(xTr_ana,Tr_ana,color='red')
temperature.plot(xTg_ana,Tg_ana,color='blue')
temperature.set_xlabel('Distance (cm)')
temperature.set_ylabel('Temperature (K)')
levels = temperature.twinx()
majorLocatorY = MultipleLocator(1.0)
levels.yaxis.set_major_locator(majorLocatorY)
levels.plot(x,amrlev,color='black',ls='dotted')
levels.set_ylim([7,11])
levels.set_ylabel('AMR Level')
temperature.set_xlim([xmin,xmax])
temperature.set_ylim([ymin,ymax])

savefig('rshock-Mach5.pdf',bbox_inches='tight')
