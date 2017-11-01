#!/bin/bash

testname="nimhd-diffusion-ad";

USEPYTHON=$1;

../visu/amr2cube -inp output_00002 -out data1.dat -typ 6 -ymi 0.5 -yma 0.52 -zmi 0.5 -zma 0.52 -fil txt;
../visu/amr2cube -inp output_00002 -out data2.dat -typ 6 -xmi 0.5 -xma 0.52 -ymi 0.5 -yma 0.52 -fil txt;
../visu/amr2cube -inp output_00002 -out data3.dat -typ 6 -ymi 0.5 -yma 0.52 -fil txt -bly T;

cat data1.dat data2.dat data3.dat > data.dat

if [ ${USEPYTHON} -eq 1 ] ; then

    python plot-${testname}.py;

else

    gnuplot plot-${testname}.gp;
    ps2pdf ${testname}.ps;
    rm ${testname}.ps;

fi

exit;
