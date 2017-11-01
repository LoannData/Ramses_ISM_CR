reset

file1 = 'data1.dat'
file2 = 'data2.dat'
file3 = 'data3.dat'
file4 = 'data4.dat'
file5 = 'data5.dat'
file6 = 'nimhd-cshock-ad-ana.dat'

lmin = int(system(sprintf("grep levelmin nimhd-cshock-ad.nml | cut -d '=' -f2")))
lmax = int(system(sprintf("grep levelmax nimhd-cshock-ad.nml | cut -d '=' -f2")))

tt = system(sprintf("cat time.dat"))
t  = tt + 0.0

set term post enh color portrait
set output 'nimhd-cshock-ad.ps'

set multiplot layout 3,2

set xrange [0.0:1.0]
unset key

set xlabel 'Distance x (cm)'

set ylabel 'Density (g/cm3)' offset 2.0
plot file1 u 1:4 w p pt 4 ps 1 lw 2, file6 u 1:6 w l lt -1 lw 2

set ylabel 'Vx (cm/s)' offset 2.0
plot file2 u 1:4 w p pt 4 ps 1 lw 2, file6 u 1:4 w l lt -1 lw 2

set ylabel 'Vy (cm/s)' offset 2.0
plot file3 u 1:4 w p pt 4 ps 1 lw 2, file6 u 1:5 w l lt -1 lw 2

set ylabel 'By' offset 2.0
plot file4 u 1:4 w p pt 4 ps 1 lw 2, file6 u 1:2 w l lt -1 lw 2

set ylabel 'Pressure (g/cm/s2)' offset 2.0
plot file5 u 1:4 w p pt 4 ps 1 lw 2, file6 u 1:3 w l lt -1 lw 2

unset multiplot

unset output
