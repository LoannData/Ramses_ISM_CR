reset

ar = 7.56591469318689378e-015
kb = 1.3806200e-16
mh = 1.66e-24

file1 = 'data1.dat'
file2 = 'data2.dat'
file3 = 'data3.dat'
file4 = 'data4.dat'
file5 = 'data5.dat'

tt = system(sprintf("cat time.dat"))
t  = tt + 0.0

unit_l = system(sprintf("grep unit_l output_00001/info_00001.txt | cut -d '=' -f2"))
unit_d = system(sprintf("grep unit_d output_00001/info_00001.txt | cut -d '=' -f2"))
unit_t = system(sprintf("grep unit_t output_00001/info_00001.txt | cut -d '=' -f2"))

scale_l = unit_l + 0.0
scale_d = unit_d + 0.0
scale_t = unit_t + 0.0

set term post enh color portrait
set output 'collapse-rhd.ps'

load "< awk '/t/ {print $0}' data2.dat"
load "< awk '/hmax/ {print $0}' data5.dat"

set multiplot

########### Temperature/Bfield/Density ############

reset

unset key
set lmargin at screen 0.10
set rmargin at screen 0.90
set bmargin at screen 0.70
set tmargin at screen 0.98

set label "t= %4.3f",t," Kyears" at screen 0.15, screen 0.95

set logscale x
set logscale y
set xlabel 'Density (g/cm^3)'
set ylabel 'Temperature (K)' tc lt 1
set y2label 'Magnetic field (G)' tc lt 3 offset -1.5
set format x "10^{%L}"
set format y "10^{%L}"
set ytics offset 0.6
set xrange[1.0e-19:5.0e-10]
set yrange[8.0:1.0e3]
set ytics nomirror
set y2tics offset -0.8
set logscale y2
set format y2 "10^{%L}"
plot file1 u (($7)*scale_d):((($18)*scale_d*((scale_l/scale_t)**2)/ar)**0.25) every 20 w d axes x1y1, file1 u (($7)*scale_d):(sqrt(0.25*((($11)+($14))**2+(($12)+($15))**2+(($13)+($16))**2)*4.0*pi*scale_d*((scale_l/scale_t)**2))) every 20 w d lc 3 axes x1y2

########### Column density ############

reset

unset key

set palette defined (0 "#000090",1 "#000fff",2 "#0090ff",3 "#0fffee",4 "#90ff70",5 "#ffee00",6 "#ff7000",7 "#ee0000",8 "#7f0000")

set lmargin at screen 0.45
set rmargin at screen 0.90
set bmargin at screen 0.10
set tmargin at screen 0.43

set label 'Column density (g/cm^2)' at screen 0.55, screen 0.60

set view map
set xlabel 'Distance (AU)'
set ylabel ''
set format y ''
set format cb "10^{%L}"
set logscale cb
set xtics 200
set xrange[-500:500]
set cbrange[1.0e-01:1.0e+04]
splot file2 index 1 u 2:3:4 w image

set lmargin at screen 0.45
set rmargin at screen 0.675
set bmargin at screen 0.43
set tmargin at screen 0.58

set xlabel ''
unset colorbox
set format x ''
splot file3 index 1 u 2:3:4 w image

set lmargin at screen 0.675
set rmargin at screen 0.90

splot file4 index 1 u 2:3:4 w image

########### Density azimuthal average ############

reset

set palette model RGB defined (0.2 "white", 0.45 "orange", 1 "red", 1.1 "black")

set view map
unset key

set xlabel 'Radius (AU)'
set ylabel 'Distance h (AU)' offset 2.0

range=hmax

set yrange [-range:range]
set xrange [0:range]
set palette maxcolors 100

set lmargin at screen 0.10
set rmargin at screen 0.35
set bmargin at screen 0.10
set tmargin at screen 0.58

set logscale cb
set format cb "10^{%L}"
set cbtics offset -0.8
set ytics rotate left
set xtics 500
splot file5 index 1 u 2:1:(($4/$3)) w image

set xlabel ''
set ylabel ''
set format x ''
set format y ''

set cbrange [0:*]
set contour
set cntrparam levels incremental -20,0.5,-6
unset surface
unset clabel
splot file5 index 1 u 2:1:(log10(abs($4/$3))) w l lt 1 lc 7 lw 2

set surface
unset contour

set cbtics offset 0,2.2
set cblabel "velocity in the rOz plane (km/s)" offset 0.0,5.0
set format cb "%1.1f"
set palette model RGB defined (0 "violet", 0.2 "blue", 0.6 "black")
set palette model RGB defined (0.1 "white", 0.2 "blue", 0.6 "green")
set colorbox horiz user origin 0.1,0.585 size 0.25,0.01
unset logscale cb
norm=range/18
splot file5 every 9:9 index 1 u 2:1:($7/$3/1000.):($5/$7*norm):($6/$7*norm):(0) w vectors filled lw 1.5 palette

unset multiplot

unset output
