from pylab import *

fig = matplotlib.pyplot.figure()
ratio = 1.3
sizex = 10.0
fig.set_size_inches(sizex,ratio*sizex)



# Read in data
data1  = loadtxt('cube1.dat')
x      = data[:,0]
By     = data[:,6]

data2  = loadtxt('cube2.dat')
z      = data[:,2]
By     = data[:,6]

data3  = loadtxt('cube3.dat')
x      = data[:,0]
z      = data[:,2]
By     = data[:,6]



# Time
t = loadtxt('time.dat')
# Compute analytical solution
mu    = 2.0 # dimensionality of the problem
beta  = 2.0
dx    = 0.5**lmax
dxx   = 0.0
alpha = -mu / (2.0+mu*beta)
delta = 1.0 / (2.0+mu*beta)
eta   = ((dx**mu/pi)/((0.5*delta*beta)**(1.0/beta) *gamma(0.5*mu)*gamma(1.0/beta+1.0)/gamma(1.0/beta+1.0+0.5*mu)))**(1.0/(mu+2.0/beta))
A     = sqrt(0.5*delta*beta*eta**2)

ana1 = abs(x-(0.5+dxx/2)) < eta*t**delta ? A*t**(alpha)*(1.-((x-(0.5+dxx/2))/(eta*t**delta))**2)**(1.0/beta) : 0.0



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
