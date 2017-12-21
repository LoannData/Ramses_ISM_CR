import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses

fig = plt.figure()
ratio = 0.8
sizex = 12.0
fig.set_size_inches(sizex,ratio*sizex)

ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)
ax5 = ax1.twinx()
ax6 = ax3.twinx()

# Load RAMSES output
data   = visu_ramses.load_snapshot(2)
order  = data["x"].argsort()
x      = data["x"][order]
amrlev = data["level"][order]
pr     = data["cosmic_rays_pressure_1"][order]
ux     = np.abs(data["velocity_x"][order])
p      = data["thermal_pressure"][order]

# Compute analytical solution
t = data["time"]
chi=1.0e10
x0=0.5
E0=1.e5
ana2 = E0/(2.0*(chi*np.pi*t)**.5)
ana = 1.+ana2*np.exp(-(x-x0)**2/(4.0*t*chi))

majorLocatorY = plt.MultipleLocator(1.0)

print min(p),max(p),np.mean(p)
print min(pr),max(pr),np.mean(pr)
print min(ana),max(ana),np.mean(ana)

# Cosmic rays pressure
ax1.semilogy(x,pr,'o',color='black',markerfacecolor='none')
ax1.semilogy(x,ana,color='red')
ax1.set_xlabel('Distance (cm)')
ax1.set_ylabel('P')
ax5.yaxis.set_major_locator(majorLocatorY)
ax5.plot(x,amrlev,color='black',ls='dotted')
ax5.set_ylabel('AMR Level')

# Velocity
ax2.semilogy(x,ux+1.,'o',color='black',markerfacecolor='none')
ax2.set_xlabel('Distance (cm)')
ax2.set_ylabel('Velocity (cm/s)')

# Relative error
ax3.semilogy(x,abs(pr-ana)/ana,color='red')
ax3.set_xlabel('Distance (cm)')
ax3.set_ylabel('Relative error')
ax6.yaxis.set_major_locator(majorLocatorY)
ax6.plot(x,amrlev,color='black',ls='dotted')
ax6.set_ylabel('AMR Level')

# Pressure ratio
ax4.semilogy(x,pr/p,'o',color='black',markerfacecolor='none')
ax4.set_xlabel('Distance (cm)')
ax4.set_ylabel('Cosmic rays Pressure/Thermal pressure')

fig.subplots_adjust(wspace=0.33)
fig.savefig('heaviside-cr.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data,'heaviside-cr',tolerance={"velocity_x":1.5e-12})
