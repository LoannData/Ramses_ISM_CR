import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import visu_ramses

fig = plt.figure()
ratio = 0.78
sizex = 20.0
fig.set_size_inches(sizex,ratio*sizex)

ax1 = fig.add_subplot(331)
ax2 = fig.add_subplot(332)
ax3 = fig.add_subplot(333)
ax4 = fig.add_subplot(334)
ax5 = fig.add_subplot(335)
ax6 = fig.add_subplot(336)
ax7 = fig.add_subplot(337)
ax8 = fig.add_subplot(338)
ax9 = fig.add_subplot(339)

# Load RAMSES output
data = visu_ramses.load_snapshot(2)

au = 1.5e13
alpha = 0.7

scale_d = data["unit_d"]
scale_l = data["unit_l"]
scale_t = data["unit_t"]
scale_b = np.sqrt(4.0*np.pi*scale_d*(scale_l/scale_t)**2)
xraw    = data["x"]*scale_l/au
yraw    = data["y"]*scale_l/au
zraw    = data["z"]*scale_l/au
dx      = data["dx"]*scale_l/au
lev     = data["level"]
rho     = np.log10(data["density"]*scale_d)
ux      = data["velocity_x"]*scale_l/scale_t/1.0e5
uy      = data["velocity_y"]*scale_l/scale_t/1.0e5
uz      = data["velocity_z"]*scale_l/scale_t/1.0e5
bx      = 0.5*(data["B_left_x"]+data["B_right_x"])*scale_b
by      = 0.5*(data["B_left_y"]+data["B_right_y"])*scale_b
bz      = 0.5*(data["B_left_z"]+data["B_right_z"])*scale_b
B       = np.log10(np.sqrt(bx**2 + by**2 + bz**2))
T       = np.log10(data["temperature"])
u       = np.sqrt(ux**2 + uy**2 + uz**2)

dmin = -19.5
dmax = -10.0
tmin = 0.8
tmax = 2.7
bmin = -5.0
bmax = 0.0
umin = -0.1
umax = 4.0

nx = 129
# Construct some edge specifiers for the histogram2d function call
d_edges = np.linspace(dmin,dmax,nx)
t_edges = np.linspace(tmin,tmax,nx)
b_edges = np.linspace(bmin,bmax,nx)
u_edges = np.linspace(umin,umax,nx)
# Call the numpy histogram2d function
za, yedges1, xedges1 = np.histogram2d(T,rho,bins=(t_edges,d_edges))
zb, yedges1, xedges1 = np.histogram2d(B,rho,bins=(b_edges,d_edges))
zc, yedges1, xedges1 = np.histogram2d(u,rho,bins=(u_edges,d_edges))
with np.errstate(divide="ignore",invalid="ignore"):
    #z1 = np.ma.masked_where(za == 0.0, np.log10(za))
    #z2 = np.ma.masked_where(zb == 0.0, np.log10(zb))
    z1 = np.log10(za)
    z2 = np.log10(zb)
    z3 = np.log10(zc)
# In the contour plots, x and y are the centers of the cells, instead of the edges.
d_mesh = np.zeros([nx-1])
t_mesh = np.zeros([nx-1])
b_mesh = np.zeros([nx-1])
u_mesh = np.zeros([nx-1])
for i in range(nx-1):
    d_mesh[i] = 0.5*(d_edges[i]+d_edges[i+1])
    t_mesh[i] = 0.5*(t_edges[i]+t_edges[i+1])
    b_mesh[i] = 0.5*(b_edges[i]+b_edges[i+1])
    u_mesh[i] = 0.5*(u_edges[i]+u_edges[i+1])
# Plot histograms
nc = 21
cont1 = ax1.contourf(d_mesh,t_mesh,z1,nc,cmap='Reds')
cont2 = ax2.contourf(d_mesh,b_mesh,z2,nc,cmap='Blues')
cont3 = ax3.contourf(d_mesh,u_mesh,z3,nc,cmap='Greens')
cont4 = ax1.contour (d_mesh,t_mesh,za,colors='r',levels=[1.0])
cont5 = ax2.contour (d_mesh,b_mesh,zb,colors='b',levels=[1.0])
cont6 = ax3.contour (d_mesh,u_mesh,zc,colors='g',levels=[1.0])

ax1.set_xlabel('log(rho)')
ax1.set_ylabel('log(T)')
ax2.set_xlabel('log(rho)')
ax2.set_ylabel('log(B)')
ax3.set_xlabel('log(rho)')
ax3.set_ylabel('Velocity')


# SLICES =====================================

dx_im = 3000.0

# Re-centre coordinates
x = xraw-0.5*data["boxlen"]*scale_l/au
y = yraw-0.5*data["boxlen"]*scale_l/au
z = zraw-0.5*data["boxlen"]*scale_l/au

dist = np.sqrt(x**2+y**2+z**2) - np.sqrt(3.0)*0.5*dx

cube = np.where(np.logical_and(np.abs(z) <= 0.5000000001*dx,np.abs(dist) <= dx_im*0.5*np.sqrt(2.0)))
im_x = x[cube]
im_y = y[cube]

xmin = -0.5*dx_im
xmax =  0.5*dx_im
ymin = -0.5*dx_im
ymax =  0.5*dx_im

nx  = 128
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([im_x,im_y])
z1 = griddata(points,rho[cube],(grid_x,grid_y),method='nearest')
z2 = griddata(points,T[cube]  ,(grid_x,grid_y),method='nearest')
z3 = griddata(points,ux[cube] ,(grid_x,grid_y),method='nearest')
z4 = griddata(points,uy[cube] ,(grid_x,grid_y),method='nearest')
z5 = np.around(griddata(points,lev[cube],(grid_x,grid_y),method='nearest'))

nc=21
im1 = ax4.contourf(xpx,ypx,z1,nc,cmap='jet')
im2 = ax7.contourf(xpx,ypx,z2,nc,cmap='hot')

ctr = ax4.contour(xpx,ypx,z5,colors='w',levels=range(0,20))
ax4.clabel(ctr,inline=1,fmt="%i")
vskip = 6
vec = ax7.quiver(xpx[::vskip],ypx[::vskip],z3[::vskip,::vskip],z4[::vskip,::vskip],color="w")

cb1 = plt.colorbar(im1,ax=ax4,label='log(Density)')
cb2 = plt.colorbar(im2,ax=ax7,label='log(T)')
cb1.ax.yaxis.set_label_coords(-1.1,0.5)
cb2.ax.yaxis.set_label_coords(-1.1,0.5)

xsink1  = data["sink1"][3]*scale_l/au-0.5*data["boxlen"]*scale_l/au
ysink1  = data["sink1"][4]*scale_l/au-0.5*data["boxlen"]*scale_l/au
zsink1  = data["sink1"][5]*scale_l/au-0.5*data["boxlen"]*scale_l/au
xsink2  = data["sink2"][3]*scale_l/au-0.5*data["boxlen"]*scale_l/au
ysink2  = data["sink2"][4]*scale_l/au-0.5*data["boxlen"]*scale_l/au
zsink2  = data["sink2"][5]*scale_l/au-0.5*data["boxlen"]*scale_l/au

circle1 = plt.Circle((xsink1,ysink1), data["r_sink"]*data["boxlen"]*scale_l/au,facecolor='w',edgecolor="k",linewidth=2,alpha=alpha)
circle2 = plt.Circle((xsink1,ysink1), data["r_sink"]*data["boxlen"]*scale_l/au,facecolor='w',edgecolor="k",linewidth=2,alpha=alpha)
circle3 = plt.Circle((xsink2,ysink2), data["r_sink"]*data["boxlen"]*scale_l/au,facecolor='w',edgecolor="k",linewidth=2,alpha=alpha)
circle4 = plt.Circle((xsink2,ysink2), data["r_sink"]*data["boxlen"]*scale_l/au,facecolor='w',edgecolor="k",linewidth=2,alpha=alpha)
ax4.add_artist(circle1)
ax4.add_artist(circle3)
ax7.add_artist(circle2)
ax7.add_artist(circle4)

ax4.set_xlabel('x')
ax4.set_ylabel('y')
ax7.set_xlabel('x')
ax7.set_ylabel('y')
ax4.set_aspect('equal')
ax7.set_aspect('equal')
ax4.set_xlim([xmin,xmax])
ax4.set_ylim([ymin,ymax])
ax7.set_xlim([xmin,xmax])
ax7.set_ylim([ymin,ymax])


# Sink 1 ===========================

dx_im = 500.0
xsink = data["sink1"][3]*scale_l/au
ysink = data["sink1"][4]*scale_l/au
zsink = data["sink1"][5]*scale_l/au

# Re-centre coordinates
x = xraw-xsink
y = yraw-ysink
z = zraw-zsink

dist = np.sqrt(x**2+y**2+z**2) - np.sqrt(3.0)*0.5*dx

cube1 = np.where(np.logical_and(np.abs(z) <= 0.5000000001*dx,np.abs(dist) <= dx_im*0.5*np.sqrt(2.0)))
im_x1 = x[cube1]
im_y1 = y[cube1]
cube2 = np.where(np.logical_and(np.abs(y) <= 0.5000000001*dx,np.abs(dist) <= dx_im*0.5*np.sqrt(2.0)))
im_x2 = x[cube2]
im_y2 = z[cube2]

xmin = -0.5*dx_im
xmax =  0.5*dx_im
ymin = -0.5*dx_im
ymax =  0.5*dx_im

nx  = 128
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([im_x1,im_y1])
z1 = griddata(points,rho[cube1],(grid_x,grid_y),method='nearest')
z3 = griddata(points,ux[cube1] ,(grid_x,grid_y),method='nearest')
z4 = griddata(points,uy[cube1] ,(grid_x,grid_y),method='nearest')

points = np.transpose([im_x2,im_y2])
q1 = griddata(points,rho[cube2],(grid_x,grid_y),method='nearest')
q3 = griddata(points,bx[cube2] ,(grid_x,grid_y),method='nearest')
q4 = griddata(points,bz[cube2] ,(grid_x,grid_y),method='nearest')

nc=21
im1 = ax5.contourf(xpx,ypx,z1,nc,cmap='jet')
im2 = ax6.contourf(xpx,ypx,q1,nc,cmap='jet')
vskip = 6
vec = ax5.quiver(xpx[::vskip],ypx[::vskip],z3[::vskip,::vskip],z4[::vskip,::vskip],color="k")
stm = ax6.streamplot(xpx,ypx,q3,q4,color="w")

cb1 = plt.colorbar(im1,ax=ax5,label='log(Density)')
cb2 = plt.colorbar(im2,ax=ax6,label='log(Density)')
cb1.ax.yaxis.set_label_coords(-1.1,0.5)
cb2.ax.yaxis.set_label_coords(-1.1,0.5)


circle1 = plt.Circle((0, 0), data["r_sink"]*data["boxlen"]*scale_l/au,facecolor='w',edgecolor="k",linewidth=2,alpha=alpha)
ax5.add_artist(circle1)
circle2 = plt.Circle((0, 0), data["r_sink"]*data["boxlen"]*scale_l/au,facecolor='w',edgecolor="k",linewidth=2,alpha=alpha)
ax6.add_artist(circle2)

ax5.set_xlabel('x')
ax5.set_ylabel('y')
ax5.set_aspect('equal')
ax5.set_xlim([xmin,xmax])
ax5.set_ylim([ymin,ymax])
ax6.set_xlabel('x')
ax6.set_ylabel('z')
ax6.set_aspect('equal')
ax6.set_xlim([xmin,xmax])
ax6.set_ylim([ymin,ymax])

msink = data["sink1"][1]
asink = data["sink1"][15]
ax5.text(-230,220,'Sink 1: %.3f Msun' % msink, color='k',bbox=dict(facecolor='w', edgecolor='k'),ha='left',va='center')
ax6.text(-230,220,'Sink 1: %.1f yr' % asink, color='k',bbox=dict(facecolor='w', edgecolor='k'),ha='left',va='center')


# Sink 2 ===========================

xsink = data["sink2"][3]*scale_l/au
ysink = data["sink2"][4]*scale_l/au
zsink = data["sink2"][5]*scale_l/au

# Re-centre coordinates
x = xraw-xsink
y = yraw-ysink
z = zraw-zsink

dist = np.sqrt(x**2+y**2+z**2) - np.sqrt(3.0)*0.5*dx

cube1 = np.where(np.logical_and(np.abs(z) <= 0.5000000001*dx,np.abs(dist) <= dx_im*0.5*np.sqrt(2.0)))
im_x1 = x[cube1]
im_y1 = y[cube1]
cube2 = np.where(np.logical_and(np.abs(y) <= 0.5000000001*dx,np.abs(dist) <= dx_im*0.5*np.sqrt(2.0)))
im_x2 = x[cube2]
im_y2 = z[cube2]

xmin = -0.5*dx_im
xmax =  0.5*dx_im
ymin = -0.5*dx_im
ymax =  0.5*dx_im

nx  = 128
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([im_x1,im_y1])
z1 = griddata(points,rho[cube1],(grid_x,grid_y),method='nearest')
z3 = griddata(points,ux[cube1] ,(grid_x,grid_y),method='nearest')
z4 = griddata(points,uy[cube1] ,(grid_x,grid_y),method='nearest')

points = np.transpose([im_x2,im_y2])
q1 = griddata(points,rho[cube2],(grid_x,grid_y),method='nearest')
q3 = griddata(points,bx[cube2] ,(grid_x,grid_y),method='nearest')
q4 = griddata(points,bz[cube2] ,(grid_x,grid_y),method='nearest')

nc=21
im1 = ax8.contourf(xpx,ypx,z1,nc,cmap='jet')
im2 = ax9.contourf(xpx,ypx,q1,nc,cmap='jet')
vskip = 6
vec = ax8.quiver(xpx[::vskip],ypx[::vskip],z3[::vskip,::vskip],z4[::vskip,::vskip],color="k")
stm = ax9.streamplot(xpx,ypx,q3,q4,color="w")

cb1 = plt.colorbar(im1,ax=ax8,label='log(Density)')
cb2 = plt.colorbar(im2,ax=ax9,label='log(Density)')
cb1.ax.yaxis.set_label_coords(-1.1,0.5)
cb2.ax.yaxis.set_label_coords(-1.1,0.5)

# Add sinks
circle1 = plt.Circle((0, 0), data["r_sink"]*data["boxlen"]*scale_l/au,facecolor='w',edgecolor="k",linewidth=2,alpha=alpha)
ax8.add_artist(circle1)
circle2 = plt.Circle((0, 0), data["r_sink"]*data["boxlen"]*scale_l/au,facecolor='w',edgecolor="k",linewidth=2,alpha=alpha)
ax9.add_artist(circle2)

ax8.set_xlabel('x')
ax8.set_ylabel('y')
ax8.set_aspect('equal')
ax8.set_xlim([xmin,xmax])
ax8.set_ylim([ymin,ymax])
ax9.set_xlabel('x')
ax9.set_ylabel('z')
ax9.set_aspect('equal')
ax9.set_xlim([xmin,xmax])
ax9.set_ylim([ymin,ymax])

msink = data["sink1"][1]
asink = data["sink1"][15]
ax8.text(-230,220,'Sink 2: %.3f Msun' % msink, color='k',bbox=dict(facecolor='w', edgecolor='k'),ha='left',va='center')
ax9.text(-230,220,'Sink 2: %.1f yr' % asink, color='k',bbox=dict(facecolor='w', edgecolor='k'),ha='left',va='center')


#fig.subplots_adjust(wspace=0.25)
fig.savefig('collapse-sink.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data,'collapse-sink')
