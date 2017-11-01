reset

ar     = 7.56591469318689378e-015
kb     = 1.3806200e-16
mh     = 1.66e-24
mu_gas = 1.0

file1 = 'data1.dat'
file2 = 'data2.dat'

set palette defined (0 "#000090",1 "#000fff",2 "#0090ff",3 "#0fffee",4 "#90ff70",5 "#ffee00",6 "#ff7000",7 "#ee0000",8 "#7f0000")

lmin = int(system(sprintf("grep levelmin dirac3d.nml | cut -d '=' -f2")))
lmax = int(system(sprintf("grep levelmax dirac3d.nml | cut -d '=' -f2")))

tt = system(sprintf("cat time.dat"))
t  = tt + 0.0 + 1.0e-20

# Compute analytical solution
chi=1.0e10
x0=0.5
E0=1.0e5
ana2 = E0/(8.0*(chi*pi*t)**1.5)
ana(x) = ana2*exp(-(x)**2/(4.0*t*chi)) + 1.0

set term post enh color portrait
set output 'dirac3d.ps'

marginx1=0.14
marginx2=0.83

set multiplot

set lmargin at screen marginx1
set rmargin at screen marginx2
set bmargin at screen 0.55
set tmargin at screen 0.72

unset key
set logscale y
set xlabel 'Distance r (cm)'
set ylabel 'Error' offset 2.5
set xrange[0.0:0.7]
set yrange[1.0e-06:9.0]
plot  file1 u (sqrt(($1-0.5)*($1-0.5)+($2-0.5)*($2-0.5)+($3-0.5)*($3-0.5))):(abs((($18)-ana(sqrt(($1-0.5)*($1-0.5)+($2-0.5)*($2-0.5)+($3-0.5)*($3-0.5))))/(ana(sqrt(($1-0.5)*($1-0.5)+($2-0.5)*($2-0.5)+($3-0.5)*($3-0.5)))))) every 20 w d lc 0 lt 1

xmin=0.01
xmax=0.99
ymin=0.09
ymax=0.91
xxyy=(ymax-ymin)/(xmax-xmin)

set bmargin at screen 0.08
set tmargin at screen 0.47

unset logscale y
set logscale z
set logscale cb
unset cblabel
set xlabel 'Distance x (cm)'
set ylabel 'Distance z (cm)' offset 0.0
set view map
set xrange[xmin:xmax]
set yrange[ymin:ymax]
unset key
splot file2 u 1:3:4 w pm3d

set bmargin at screen 0.72
set tmargin at screen 0.99

unset xlabel
set xtics format " "
set ylabel 'Er' offset 2.5
set autoscale
set xrange[0.0:0.7]
set logscale y
set key
plot file1 u (sqrt(($1-0.5)*($1-0.5)+($2-0.5)*($2-0.5)+($3-0.5)*($3-0.5))):18 every 20 lc 0 lt 6 w d title 'Simulation', ana(x) lc 1 lt 1 lw 2 title 'Analytical'


unset multiplot

unset output
