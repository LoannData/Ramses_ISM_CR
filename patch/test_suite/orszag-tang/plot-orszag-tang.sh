#!/bin/bash

testname="orszag-tang";

USEPYTHON=$1;

../visu/amr2map -inp output_00002 -out data1.dat -typ 1 -fil ascii;

echo $(md5sum data1.dat | cut -d ' ' -f 1) > data.dat

if [ ${USEPYTHON} -eq 1 ] ; then

    python plot-${testname}.py;

else

    gnuplot plot-${testname}.gp;
    ps2pdf ${testname}.ps;
    rm ${testname}.ps;

fi

exit;
