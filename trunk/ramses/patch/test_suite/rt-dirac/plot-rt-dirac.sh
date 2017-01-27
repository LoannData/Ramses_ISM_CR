#!/bin/bash

testname="rt-dirac";

USEPYTHON=$1;

../visu/amr2cube -inp output_00002 -out data1.dat -typ 1 -ymi 0.49 -yma 0.51 -fil txt -bly T;
../visu/amr2cube -inp output_00002 -out data2.dat -typ 11 -ymi 0.49 -yma 0.51 -fil txt -bly T;

../visu/amr2cube -inp output_00002 -out data3.dat -typ 1 -xmi 0.49 -xma 0.51 -fil txt -blz T;
../visu/amr2cube -inp output_00002 -out data4.dat -typ 11 -xmi 0.49 -xma 0.51 -fil txt -blz T;

../visu/amr2cube -inp output_00002 -out data5.dat -typ 1 -zmi 0.49 -zma 0.51 -fil txt -bly T;
../visu/amr2cube -inp output_00002 -out data6.dat -typ 11 -zmi 0.49 -zma 0.51 -fil txt -bly T;

cat data1.dat data2.dat data3.dat data4.dat data5.dat data6.dat > data7.dat
echo $(md5sum data7.dat | cut -d ' ' -f 1) > data.dat

if [ ${USEPYTHON} -eq 1 ] ; then

    python plot-${testname}.py;

else

    gnuplot plot-${testname}.gp;
    ps2pdf ${testname}.ps;
    rm ${testname}.ps;

fi

exit;
