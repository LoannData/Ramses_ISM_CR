#!/bin/bash

../visu/amr2cell -inp output_00002 -out data1.dat;
../visu/amr2cube -inp output_00002 -out data2.dat -typ 12 -ymi 0.48 -yma 0.52 -fil txt -bly T;

cat data1.dat data2.dat > data3.dat;
echo $(md5sum data3.dat | cut -d ' ' -f 1) > data.dat

exit;
