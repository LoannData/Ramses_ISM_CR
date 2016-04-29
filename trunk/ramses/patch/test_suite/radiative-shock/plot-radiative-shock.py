from pylab import *

# Figure properties
fig = matplotlib.pyplot.figure()
ratio = 1.5
sizex = 8.0
fig.set_size_inches(sizex,ratio*sizex)

# Define constants
ar        = 7.56591469318689378e-015
kb        = 1.3806200e-16
mh        = 1.66e-24
mu_gas    = 1.0
scale_d   = 7.78e-10
scale_l   = 7.00e+10
scale_v   = 6.00e+05
rescale_d = 1.0e-10
rescale_v = 1.0e+05

# Load colormap
jet = cm = plt.get_cmap('jet')

lwidth2 = 4

# Read in data
data       = loadtxt('data.dat')
t          = loadtxt('time.dat')
time       = t*scale_l/scale_v
[nz,ncols] = shape(data)
ngr        = ncols - 11
amrlev     = data[:, 0]
x          = data[:, 1]*scale_l + scale_v*time
rho        = data[:, 2]*scale_d
u          = data[:, 3]*scale_v
p          = data[:, 6]
er         = data[:,11:11+ngr]
Tr         = (er*scale_d*scale_v**2/ar)**0.25
Tr_tot     = (er.sum(axis=1)*scale_d*scale_v**2/ar)**0.25
Tg         = p*scale_v**2*mu_gas*mh/(rho/scale_d*kb)

# Read reference solution
data_ana   = loadtxt('radiative-shock-ref.dat')
[nz,ncols] = shape(data_ana)
ngr_ana    = ncols - 11
x_ana      = data[:, 1]*scale_l + scale_v*time
rho_ana    = data[:, 2]*scale_d
u_ana      = data[:, 3]*scale_v
p_ana      = data[:, 6]
er_ana     = data[:,11:11+ngr_ana]
Tr_ana     = (er_ana*scale_d*scale_v**2/ar)**0.25
Tr_tot_ana = (er_ana.sum(axis=1)*scale_d*scale_v**2/ar)**0.25
Tg_ana     = p_ana*scale_v**2*mu_gas*mh/(rho_ana/scale_d*kb)

xmin = 2.0e+10
xmax = 7.0e+10

# Density
density = subplot(311)
density.plot(x_ana,rho_ana,color='red',lw=lwidth2,ls='dashed')
density.plot(x,rho,color='black')
density.set_xlabel('Distance (cm)')
density.set_ylabel('Density (g/cm3)')
density.set_xlim([xmin,xmax])

# Velocity
velocity = subplot(312)
velocity.plot(x_ana,u_ana/rescale_v,color='red',lw=lwidth2,ls='dashed')
velocity.plot(x,u/rescale_v,color='black')
velocity.set_xlabel('Distance (cm)')
velocity.set_ylabel('Velocity (km/s)')
velocity.set_xlim([xmin,xmax])

# Temperature
temperature = subplot(313)
# Scale colormap
cNorm = matplotlib.colors.Normalize(vmin=0,vmax=ngr_ana-1)
scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm,cmap=jet)
for ig in range(ngr_ana):
    colorVal = scalarMap.to_rgba(ig)
    temperature.plot(x_ana,Tr_ana[:,ig],color=colorVal,lw=lwidth2,ls='dashed')
temperature.plot(x_ana,Tr_tot_ana,color='magenta',lw=lwidth2,ls='dashed')
cNorm = matplotlib.colors.Normalize(vmin=0,vmax=ngr-1)
scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm,cmap=jet)
for ig in range(ngr):
    colorVal = scalarMap.to_rgba(ig)
    temperature.plot(x,Tr[:,ig],color=colorVal,label='Tr('+str(ig+1)+')')
temperature.plot(x,Tr_tot,color='magenta',label='Tr tot')
temperature.plot(x_ana,Tg_ana,color='black',lw=lwidth2,ls='dashed')
temperature.plot(x,Tg,color='black',label='Tgas')
temperature.set_xlabel('Distance (cm)')
temperature.set_ylabel('Temperature (K)')
levels = temperature.twinx()
majorLocatorY = MultipleLocator(1.0)
levels.yaxis.set_major_locator(majorLocatorY)
levels.plot(x,amrlev,color='black',ls='dotted',label='Level')
levels.set_ylabel('AMR Level')
temperature.set_xlim([xmin,xmax])
temperature.legend()
levels.legend(loc='lower right')

savefig('radiative-shock.pdf',bbox_inches='tight')
