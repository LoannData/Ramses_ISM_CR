reset

ar     = 7.56591469318689378e-015
kb     = 1.3806200e-16
mh     = 1.66e-24
mu_gas = 1.0

scale_d = 1.0e-13
scale_l = 1.0
scale_v = 1.0

xshift = 4205.08

lwidth = 2

file1 = 'data.dat'
file2 = 'rshock-Mach5-ana.dat'

lmin = int(system(sprintf("grep levelmin rshock-Mach5.nml | cut -d '=' -f2")))
lmax = int(system(sprintf("grep levelmax rshock-Mach5.nml | cut -d '=' -f2")))

set term post enh color
set output 'rshock-Mach5.ps'

set xlabel 'Distance x (cm)'
set xrange[-3000:500]

unset key
set ylabel 'Temperature (K)'

set y2label 'AMR level'
set ytics nomirror
set autoscale  y
set y2tics lmin,1,lmax

plot file1 u ($2-xshift):((($11)/ar)**0.25) lc 0 lt 6, file1 u ($2-xshift):(($7)*mu_gas*mh/(($3)*kb)) lc 0 lt 4, file2 u (($3)*1.e5):(($4)/(8.57274781455567/855.7)) w l lw lwidth lt 1 lc 3, file2 u 1:2 w l lw lwidth lt 1 lc 1, file1 u ($2-xshift):1 w l lt 0 lc 0 title 'Level' axes x1y2

unset output
