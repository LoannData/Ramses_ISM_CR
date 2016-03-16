from pylab import *

fig = matplotlib.pyplot.figure()
ratio = 1.0
sizex = 10.0
fig.set_size_inches(sizex,ratio*sizex)

# Read in data
data   = loadtxt('data.dat')
amrlev = data[:, 0]
x      = data[:, 1]
er     = data[:,10]

# Compute analytical solution
t = loadtxt('time.dat')
chi=1.0e10
x0=0.5
E0=1.0e5
ana2 = E0/(2.0*(chi*pi*t)**.5)
ana = 1.+ana2*exp(-(x-x0)**2/(4.0*t*chi))

# Radiative energy
erad = subplot(211)
erad.semilogy(x,er,'o',color='black',markerfacecolor='none')
erad.semilogy(x,ana,color='red')
erad.set_xlabel('Distance (cm)')
erad.set_ylabel('Radiative energy')
levels1 = erad.twinx()
majorLocatorY = MultipleLocator(1.0)
levels1.yaxis.set_major_locator(majorLocatorY)
levels1.plot(x,amrlev,color='black',ls='dotted')
levels1.set_ylabel('AMR Level')

# Relative error
error = subplot(212)
error.semilogy(x,abs(er-ana)/ana,color='red')
error.set_xlabel('Distance (cm)')
error.set_ylabel('Percentage relative error')
levels2 = error.twinx()
majorLocatorY = MultipleLocator(1.0)
levels2.yaxis.set_major_locator(majorLocatorY)
levels2.plot(x,amrlev,color='black',ls='dotted')
levels2.set_ylabel('AMR Level')

savefig('dirac.pdf',bbox_inches='tight')
