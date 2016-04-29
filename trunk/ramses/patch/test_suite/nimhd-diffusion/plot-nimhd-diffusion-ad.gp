reset

file1 = 'data1.dat'
file2 = 'data2.dat'
file3 = 'data3.dat'

lmin = int(system(sprintf("grep levelmin nimhd-diffusion-ad.nml | cut -d '=' -f2")))
lmax = int(system(sprintf("grep levelmax nimhd-diffusion-ad.nml | cut -d '=' -f2")))

tt = system(sprintf("cat time.dat"))
t  = tt + 0.0

# Compute analytical solution
mu    = 2.0 # dimensionality of the problem
beta  = 2.0
dx    = 0.5**lmax
dxx   = dx
alpha = -mu / (2.0+mu*beta)
delta = 1.0 / (2.0+mu*beta)
eta   = ((dx**mu/pi)/((0.5*delta*beta)**(1.0/beta) *gamma(0.5*mu)*gamma(1.0/beta+1.0)/gamma(1.0/beta+1.0+0.5*mu)))**(1.0/(mu+2.0/beta))
A     = sqrt(0.5*delta*beta*eta**2)
f1(x) = abs(x-(0.5+dxx/2.0)) < eta*t**delta ? A*t**(alpha)*(1.-((x-(0.5+dxx/2.0))/(eta*t**delta))**2)**(1.0/beta) : 0.0

set term post enh color portrait
set output 'nimhd-diffusion-ad.ps'

nxpanel=1
nypanel=2
panlgap=0.07

marginx1=0.14
marginx2=0.85
marginy1=0.02
marginy2=0.99

dxpanel=(marginx2-marginx1-(nxpanel-1)*panlgap)/nxpanel
dypanel=0.22

set multiplot

set lmargin at screen marginx1
set rmargin at screen marginx2
set bmargin at screen (marginy2-dypanel)
set tmargin at screen marginy2

set xrange [0.0:1.0]
set y2tics
set y2label 'error' rotate by -90 offset -2.5
set key left reverse Left
set ytics nomirror

sum = 0
set xlabel 'Distance x (cm)'
set ylabel 'By' offset 2.0
plot file1 u 1:4 w p pt 4 ps 2 lw 2 tit 'Simulation', f1(x) lt -1 lw 2 tit 'Analytical', file1 u 1:(sum=sum+($4-f1($1))**2) axis x1y2 w l lt 0 lw 3 tit 'error L2'

set bmargin at screen (marginy2-dypanel-dypanel-panlgap)
set tmargin at screen (marginy2-dypanel-panlgap)

sum = 0
set xlabel 'Distance z (cm)'
set ylabel 'By' offset 2.0
plot  file2 u 3:4 w p pt 4 ps 2 lw 2 tit 'Simulation', f1(x) lt -1 lw 2 tit 'Analytical', file2 u 3:(sum=sum+($4-f1($3))**2) axis x1y2 w l lt 0 lw 3 tit 'error L2'

xmin=0.0
xmax=1.0
ymin=0.2
ymax=0.8
xxyy=(ymax-ymin)/(xmax-xmin)

set bmargin at screen marginy1
set tmargin at screen (marginy1+xxyy*dxpanel)

nc   = 20
bmax = f1(0.5)
inc  = bmax/nc
set xlabel 'Distance x (cm)'
set ylabel 'Distance z (cm)' offset 0.0
set view map
set xrange[xmin:xmax]
set yrange[ymin:ymax]
set size ratio xxyy
# set dgrid3d 64,64
set contour
set cntrparam levels incremental 0.0,inc,bmax
unset surface
unset clabel
unset key
unset y2tics
splot file3 u 1:3:4 w l lt 1

unset multiplot

unset output
