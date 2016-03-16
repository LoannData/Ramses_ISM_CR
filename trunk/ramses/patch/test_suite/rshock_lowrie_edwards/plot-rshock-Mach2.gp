reset

ar     = 7.56591469318689378e-015
kb     = 1.3806200e-16
mh     = 1.66e-24
mu_gas = 1.0

scale_d = 1.0e-13
scale_l = 1.0
scale_v = 1.0

xshift = 1013.67

lwidth = 2

file1 = 'data.dat'
file2 = 'rshock-Mach2-ana.dat'

lmin = int(system(sprintf("grep levelmin Mach2.nml | cut -d '=' -f2")))
lmax = int(system(sprintf("grep levelmax Mach2.nml | cut -d '=' -f2")))

set term post enh color portrait
set output 'rshock-Mach2.ps'

set xlabel 'Distance x (cm)'
set xrange[-1000:1000]

nxpanel=1
nypanel=2
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

set y2label 'AMR level' offset -1.5
set ytics nomirror
set autoscale  y
set y2tics lmin,1,lmax

set ylabel 'Temperature (K)'
plot file1 u ($2-xshift):((($11)/ar)**0.25) lc 0 lt 6, file1 u ($2-xshift):(($7)*mu_gas*mh/(($3)*kb)) lc 0 lt 4, file2 u 1:5 w l lw lwidth lt 1 lc 1, file2 u 1:4 w l lw lwidth lt 1 lc 3, file1 u ($2-xshift):1 w l lt 0 lc 0 title 'Level' axes x1y2

unset y2label
unset y2tics

set bmargin at screen marginy1
set tmargin at screen (marginy1+dypanel)

set ylabel sprintf("Density (%2.1e g/cm3)",scale_d)
plot file1 u ($2-xshift):($3/scale_d) lc 0 lt 6, file2 u 1:($3/scale_d) w l lw lwidth lt 1 lc 1

unset multiplot

unset output
