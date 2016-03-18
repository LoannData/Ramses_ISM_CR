#!/bin/bash

../visu/amr2cube -inp output_00002 -out data1.dat -typ 6 -ymi 0.5 -yma 0.52 -zmi 0.5 -zma 0.52 -fil txt;
../visu/amr2cube -inp output_00002 -out data2.dat -typ 6 -xmi 0.5 -xma 0.52 -ymi 0.5 -yma 0.52 -fil txt;
../visu/amr2cube -inp output_00002 -out data3.dat -typ 6 -ymi 0.5 -yma 0.52 -fil txt -bly T;

cat data1.dat data2.dat data3.dat > data.dat

exit;
