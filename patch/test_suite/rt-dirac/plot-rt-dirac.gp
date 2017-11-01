reset

file1 = 'data1.dat'
file2 = 'data2.dat'
#file3 = 'data3.dat'
#file4 = 'data4.dat'
#file5 = 'data5.dat'
#file6 = 'data6.dat'

tt = system(sprintf("cat time.dat"))
t  = tt + 0.0

boxlen = system(sprintf("grep boxlen output_00001/info_00001.txt | cut -d '=' -f2"))
unit_l = system(sprintf("grep unit_l output_00001/info_00001.txt | cut -d '=' -f2"))
unit_d = system(sprintf("grep unit_d output_00001/info_00001.txt | cut -d '=' -f2"))
unit_t = system(sprintf("grep unit_t output_00001/info_00001.txt | cut -d '=' -f2"))
#mu_gas = system(sprintf("grep mu_gas output_00001/info_00001.txt | cut -d '=' -f2"))

scale_l = unit_l + 0.0
scale_d = unit_d + 0.0
scale_t = unit_t + 0.0

#units in cgs
Myr=1e6*365.25e0*86400e0
pc=3.08e18
mh=1.660531e-24

set palette defined (0 "#000090",1 "#000fff",2 "#0090ff",3 "#0fffee",4 "#90ff70",5 "#ffee00",6 "#ff7000",7 "#ee0000",8 "#7f0000")

set term post enh color portrait
set output 'rt-dirac.ps'

set multiplot

nxpanel=1
nypanel=2
panlgap=0.08

marginx1=0.10
marginx2=0.90
marginy1=0.07
marginy2=0.99

dxpanel=(marginx2-marginx1-(nxpanel-1)*panlgap)/nxpanel
dypanel=(marginy2-marginy1-(nypanel-1)*panlgap)/nypanel

unset key

set lmargin at screen marginx1
set rmargin at screen (marginx1+dxpanel)
set bmargin at screen (marginy2-dypanel)
set tmargin at screen marginy2

set size square
set view map
set xlabel 'Distance x (pc)'
set ylabel 'Distance z (pc)'
set cblabel 'Log (n) (cm-3)'
set title sprintf("t = %4.3f (Myr)",t*scale_t/Myr)
unset key
splot file1 u ($1/pc):($3/pc):(log10($4)) w pm3d

# unset y2label
# unset y2tics

# set lmargin at screen (marginx1+dxpanel+panlgap)
# set rmargin at screen (marginx2-dxpanel-panlgap)
# set xlabel 'Distance y (pc)'
# set ylabel 'Distance z (pc)'
# set cblabel 'Log (n) (cm-3)'
# set title sprintf("t = %4.3f (Myr)",t*scale_t/Myr)
# splot file3 u ($2/pc):($3/pc):(log10($4)) w pm3d

# set lmargin at screen (marginx2-dxpanel)
# set rmargin at screen marginx2
# set xlabel 'Distance x (pc)'
# set ylabel 'Distance y (pc)'
# set cblabel 'Log (n) (cm-3)'
# set title sprintf("t = %4.3f (Myr)",t*scale_t/Myr)
# splot file5 u ($1/pc):($2/pc):(log10($4)) w pm3d

set lmargin at screen marginx1
set rmargin at screen (marginx1+dxpanel)
set bmargin at screen marginy1
set tmargin at screen (marginy1+dypanel)
unset title
set cblabel 'Log(P) '
unset key
splot file2 u ($1/pc):($3/pc):(log10($4)) w pm3d

# set lmargin at screen (marginx1+dxpanel+panlgap)
# set rmargin at screen (marginx2-dxpanel-panlgap)
# set xlabel 'Distance y (pc)'
# set ylabel 'Distance z (pc)'
# set cblabel 'Log (P)'
# splot file4 u ($2/pc):($3/pc):(log10($4)) w pm3d

# set lmargin at screen (marginx2-dxpanel)
# set rmargin at screen marginx2
# set xlabel 'Distance x (pc)'
# set ylabel 'Distance y (pc)'
# set cblabel 'Log (P)'
# splot file6 u ($1/pc):($2/pc):(log10($4)) w pm3d

unset multiplot

unset output
