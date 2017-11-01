#!/bin/bash

testname="nimhd-cshock-ad";

USEPYTHON=$1;

../visu/amr2cube -inp output_00002 -out data1.dat -typ  1 -ymi 0.5 -yma 0.54 -zmi 0.5 -zma 0.54 -fil txt;
../visu/amr2cube -inp output_00002 -out data2.dat -typ  2 -ymi 0.5 -yma 0.54 -zmi 0.5 -zma 0.54 -fil txt;
../visu/amr2cube -inp output_00002 -out data3.dat -typ  3 -ymi 0.5 -yma 0.54 -zmi 0.5 -zma 0.54 -fil txt;
../visu/amr2cube -inp output_00002 -out data4.dat -typ  6 -ymi 0.5 -yma 0.54 -zmi 0.5 -zma 0.54 -fil txt;
../visu/amr2cube -inp output_00002 -out data5.dat -typ 11 -ymi 0.5 -yma 0.54 -zmi 0.5 -zma 0.54 -fil txt;

cat data1.dat data2.dat data3.dat data4.dat data5.dat > data.dat

if [ ${USEPYTHON} -eq 1 ] ; then

    python plot-${testname}.py;

else

    gnuplot plot-${testname}.gp;
    ps2pdf ${testname}.ps;
    rm ${testname}.ps;

fi

exit;
