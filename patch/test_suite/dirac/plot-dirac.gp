reset

ar     = 7.56591469318689378e-015
kb     = 1.3806200e-16
mh     = 1.66e-24
mu_gas = 1.0

scale_d = 1.0
scale_l = 1.0
scale_v = 1.0

lwidth = 2

file1 = 'data.dat'

lmin = int(system(sprintf("grep levelmin dirac.nml | cut -d '=' -f2")))
lmax = int(system(sprintf("grep levelmax dirac.nml | cut -d '=' -f2")))

tt = system(sprintf("cat time.dat"))
t  = tt + 0.0 + 1.0e-20

chi=1.0e10
x0=0.5
E0=1.0e5
ana2 = E0/(2.0*(chi*pi*t)**.5)
ana(x) = 1.+ana2*exp(-(x-x0)**2/(4.0*t*chi))

set term post enh color portrait
set output 'dirac.ps'

set xlabel 'Distance x (cm)'

set multiplot layout 2, 1
unset key

set y2label 'AMR level'
set ytics nomirror
set y2tics
set autoscale  y
set y2tics lmin,1,lmax

set ylabel 'Radiative energy'
set logscale y
plot file1 u 2:11 lc 0 lt 6, ana(x) w l lw lwidth lt 1 lc 1 axes x1y1, file1 u 2:1 w l lt 0 lc 0 axes x1y2

set ylabel 'Percentage relative error'
plot file1 u 2:(abs(($11)-1.0-ana2*exp(-(($2)-x0)**2/(4.0*tt*chi)))/(1.0+ana2*exp(-(($2)-x0)**2/(4.0*tt*chi)))) w l lw lwidth axes x1y1, file1 u 2:1 w l lt 0 lc 0 axes x1y2

unset multiplot

unset output
