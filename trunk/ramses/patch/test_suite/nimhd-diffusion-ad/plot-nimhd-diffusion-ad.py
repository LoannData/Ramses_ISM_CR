from pylab import *
import numpy as np
import re
import math

fig = matplotlib.pyplot.figure()
ratio = 1.3
sizex = 10.0
fig.set_size_inches(sizex,ratio*sizex)

# Read in data
data1  = loadtxt('data1.dat')
x1      = data1[:,0]
By1     = data1[:,3]

data2  = loadtxt('data2.dat')
z2      = data2[:,2]
By2     = data2[:,3]

data3  = loadtxt('data3.dat')
x3      = data3[:,0]
z3      = data3[:,2]
By3     = data3[:,3]

# Time
t = loadtxt('time.dat')

for line in open('nimhd-diffusion-ad.nml'):
    if re.search('levelmax',line):
        lmax=int(line.split('=')[1])

# Compute analytical solution
mu    = 2.0 # dimensionality of the problem
beta  = 2.0
dx    = 0.5**lmax
dxx   = dx
alpha = -mu / (2.0+mu*beta)
delta = 1.0 / (2.0+mu*beta)
eta   = ((dx**mu/np.pi)/((0.5*delta*beta)**(1.0/beta) *math.gamma(0.5*mu)*math.gamma(1.0/beta+1.0)/math.gamma(1.0/beta+1.0+0.5*mu)))**(1.0/(mu+2.0/beta))
A     = sqrt(0.5*delta*beta*eta**2)

x=np.arange(0,1,0.01)
ana = np.zeros(size(x))
ana[abs(x-(0.5+dxx/2)) < eta*t**delta] = A*t**(alpha)*(1.-((x[abs(x-(0.5+dxx/2)) < eta*t**delta]-(0.5+dxx/2))/(eta*t**delta))**2)**(1.0/beta)

ana1 = np.zeros(size(x1))
ana1[abs(x1-(0.5+dxx/2)) < eta*t**delta] = A*t**(alpha)*(1.-((x1[abs(x1-(0.5+dxx/2)) < eta*t**delta]-(0.5+dxx/2))/(eta*t**delta))**2)**(1.0/beta)

ana2 = np.zeros(size(z2))
ana2[abs(z2-(0.5+dxx/2)) < eta*t**delta] = A*t**(alpha)*(1.-((z2[abs(z2-(0.5+dxx/2)) < eta*t**delta]-(0.5+dxx/2))/(eta*t**delta))**2)**(1.0/beta)

# By(x)
by1 = subplot(311)
by1.plot(x1,By1,'o',color='red',label='simulation')
by1.plot(x,ana,color='black',label='analytical')
by1.set_xlabel('Distance x (cm)')
by1.set_ylabel('By')
by1.legend()
levels1 = by1.twinx()
majorLocatorY = MultipleLocator(2e-8)
levels1.yaxis.set_major_locator(majorLocatorY)
levels1.plot(x1,(By1-ana1)**2,color='black',ls='dotted',label='error L2')
levels1.set_ylabel('error')
levels1.legend(loc='lower right')
levels1.set_ylim([0,1.4e-7])

# By(z)
by2 = subplot(312)
by2.plot(z2,By2,'o',color='red',label='simulation')
by2.plot(x,ana,color='black',label='analytical')
by2.set_xlabel('Distance z (cm)')
by2.set_ylabel('By')
by2.legend()
levels2 = by2.twinx()
majorLocatorY = MultipleLocator(2e-8)
levels2.yaxis.set_major_locator(majorLocatorY)
levels2.plot(z2,((By2-ana2)**2),color='black',ls='dotted',label='error L2')
levels2.set_ylabel('error')
levels2.legend(loc='lower right')
levels2.set_ylim([0,1.4e-7])

# By(x,z)
By3=By3.reshape(int(sqrt(shape(By3)[0])),int(sqrt(shape(By3)[0])))
by3 = subplot(313)
by3.contour(By3,extent=(0.,1.,0.2,0.8))
by3.set_xlabel('Distance x (cm)')
by3.set_ylabel('Distance z (cm)')
by3.set_aspect(1./0.6)

savefig('nimhd-diffusion-ad.pdf',bbox_inches='tight')
