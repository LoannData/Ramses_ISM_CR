#!/bin/bash

../visu/amr2cell -inp output_00002 -out data1.dat;

../visu/amr2cold -inp output_00002 -out data2.dat -typ 1 -fil vtk -dir 1 -fen 3.0 -lma 12;
../visu/amr2cold -inp output_00002 -out data3.dat -typ 1 -fil vtk -dir 2 -fen 3.0 -lma 12;
../visu/amr2cold -inp output_00002 -out data4.dat -typ 1 -fil vtk -dir 3 -fen 3.0 -lma 12;

../visu/amr2azav -inp output_00002 -out data5.dat -typ 18 -fil vtk -dir 4 -fen 1.0 -lma 13 -nor 300;

cat data1.dat data2.dat data3.dat data4.dat data5.dat > data6.dat;
echo $(md5sum data6.dat | cut -d ' ' -f 1) > data.dat

rm barotropic_eos.dat Hosokawa_track.dat vaytet_grey_opacities*.bin groups.dat init_turb.data tab_eos.dat res*.dat

exit;
