import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from scipy.interpolate import griddata

fig = plt.figure()
ratio = 0.75
sizex = 12.0
fig.set_size_inches(sizex,ratio*sizex)

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x      = data["x"]
y      = data["y"]
dx     = data["dx"]
rho    = data["density"]
p      = data["thermal_pressure"]
lev    = data["level"]
u      = np.sqrt(data["velocity_x"]**2 + data["velocity_y"]**2)

xmin = np.amin(x-0.5*dx)
xmax = np.amax(x+0.5*dx)
ymin = np.amin(y-0.5*dx)
ymax = np.amax(y+0.5*dx)

nx  = 128
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([x,y])
z1 = griddata(points,rho,(grid_x,grid_y),method='nearest')
z2 = griddata(points,u  ,(grid_x,grid_y),method='nearest')
z3 = griddata(points,p  ,(grid_x,grid_y),method='nearest')
z4 = griddata(points,lev,(grid_x,grid_y),method='nearest')

nc=17
im1 = ax1.contourf(xpx,ypx,z1,levels=np.linspace(np.amin(rho),np.amax(rho),nc))
im2 = ax2.contourf(xpx,ypx,z2,levels=np.linspace(np.amin(u),np.amax(u),nc),cmap='cubehelix')
im3 = ax3.contourf(xpx,ypx,z3,levels=np.linspace(np.amin(p),np.amax(p),nc),cmap='hot')
im4 = ax4.contourf(xpx,ypx,z4,levels=range(int(np.amin(lev))-1,int(np.amax(lev))+1),cmap='gnuplot')

cb1 = plt.colorbar(im1,ax=ax1,label='Density')
cb2 = plt.colorbar(im2,ax=ax2,label='Velocity')
cb3 = plt.colorbar(im3,ax=ax3,label='Pressure')
cb4 = plt.colorbar(im4,ax=ax4,label='AMR Level')
cb1.ax.yaxis.set_label_coords(-1.1,0.5)
cb2.ax.yaxis.set_label_coords(-1.1,0.5)
cb3.ax.yaxis.set_label_coords(-1.1,0.5)
cb4.ax.yaxis.set_label_coords(-1.1,0.5)

ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax4.set_xlabel('x')
ax4.set_ylabel('y')
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')
ax1.set_xlim([xmin,xmax])
ax1.set_ylim([ymin,ymax])
ax2.set_xlim([xmin,xmax])
ax2.set_ylim([ymin,ymax])
ax3.set_xlim([xmin,xmax])
ax3.set_ylim([ymin,ymax])
ax4.set_xlim([xmin,xmax])
ax4.set_ylim([ymin,ymax])

fig.subplots_adjust(wspace=0.25)
fig.savefig('implosion.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data,'implosion')
