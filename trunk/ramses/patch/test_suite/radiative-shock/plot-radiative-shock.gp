reset

ar     = 7.56591469318689378e-015
kb     = 1.3806200e-16
mh     = 1.66e-24
mu_gas = 1.0

scale_d = 7.78e-10
scale_l = 7.00e+10
scale_v = 6.00e+05

rescale_d = 1.0e-10
rescale_v = 1.0e+05

lwidth1 = 2
lwidth2 = 6

file1 = 'data.dat'
file2 = 'radiative-shock-ref.dat'

lmin = int(system(sprintf("grep levelmin tube1d.nml | cut -d '=' -f2")))
lmax = int(system(sprintf("grep levelmax tube1d.nml | cut -d '=' -f2")))

t  = system(sprintf("cat time.dat"))
tt = t + 0.0
tt = tt*scale_l/scale_v

nvar  = int(system(sprintf("grep nvar log | cut -d '=' -f3")))
ngr   = int(system(sprintf("grep groups log | cut -d ':' -f2")))

imin  = 12
imax  = imin+ngr-1

title1(n) = sprintf("Tr(%d)",n-imin+1)

set term post enh color portrait
set output 'radiative-shock.ps'

set xlabel 'Distance x (cm)'
set xrange[2.0e10:7.0e10]

nxpanel=1
nypanel=3
panlgap=0.08

marginx1=0.10
marginx2=0.90
marginy1=0.07
marginy2=0.99

dxpanel=(marginx2-marginx1-(nxpanel-1)*panlgap)/nxpanel
dypanel=(marginy2-marginy1-(nypanel-1)*panlgap)/nypanel

set multiplot
unset key

set lmargin at screen marginx1
set rmargin at screen marginx2
set bmargin at screen (marginy2-dypanel)
set tmargin at screen marginy2

set ylabel sprintf("Density (%2.1e g/cm3)",rescale_d) offset 2.0

plot file2 u (($2*scale_l)+scale_v*tt):($3*scale_d/rescale_d) w l lw lwidth2 lt 2 lc 1, file1 u (($2*scale_l)+scale_v*tt):($3*scale_d/rescale_d) w l lw lwidth1 lt 1 lc 0

set bmargin at screen (marginy1+dypanel+panlgap)
set tmargin at screen (marginy2-dypanel-panlgap)

set ylabel 'Velocity (km/s)' offset 2.0
plot file2 u (($2*scale_l)+scale_v*tt):($4*scale_v/rescale_v) w l lw lwidth2 lt 2 lc 1, file1 u (($2*scale_l)+scale_v*tt):($4*scale_v/rescale_v) w l lw lwidth1 lt 1 lc 0

set key
set ylabel 'Temperature (K)' offset 2.0

set bmargin at screen marginy1
set tmargin at screen (marginy1+dypanel)

set y2label 'AMR level'
set ytics nomirror
set y2tics
set autoscale  y
set y2tics lmin,1,lmax

plot for [i=imin:imax] file2 u (($2*scale_l)+scale_v*tt):(column(i)*scale_d*scale_v**2/ar)**0.25 w l lw lwidth2 lt 2 lc (i-imin+1) notitle axes x1y1, file2 using (($2*scale_l)+scale_v*tt):((sum [col=imin:imax] column(col)*scale_d*scale_v**2)/ar)**0.25 w l lw lwidth2 lt 2 lc (imax-imin+2) notitle axes x1y1, file2 u (($2*scale_l)+scale_v*tt):($7/$3*scale_v**2*(mh*mu_gas)/kb) w l lw lwidth2 lt 2 lc 0 notitle axes x1y1, for [i=imin:imax] file1 u (($2*scale_l)+scale_v*tt):(column(i)*scale_d*scale_v**2/ar)**0.25 w l lw lwidth1 lt 1 lc (i-imin+1) title title1(i) axes x1y1, file1 using (($2*scale_l)+scale_v*tt):((sum [col=imin:imax] column(col)*scale_d*scale_v**2)/ar)**0.25 w l lw lwidth1 lt 1 lc (imax-imin+2) title 'Tr tot' axes x1y1, file1 u (($2*scale_l)+scale_v*tt):($7/$3*scale_v**2*(mh*mu_gas)/kb) w l lw lwidth1 lt 1 lc 0 title 'Tgas' axes x1y1, file1 u (($2*scale_l)+scale_v*tt):1 w l lt 0 lc 0 title 'Level' axes x1y2

unset multiplot

unset output
