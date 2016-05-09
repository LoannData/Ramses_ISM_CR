#!/bin/bash

#######################################################################
#
# Script to run the RAMSES test suite
#
# Usage:
#   ./run_test_suite.sh
#
# Options:
#   - Run the suite in parallel (on 4 cpus):
#       ./run_test_suite.sh -p 4
#   - Run using python visualization instead of gnuplot:
#       ./run_test_suite -y
#   - Run full test suite (including additional tests):
#       ./run_test_suite -f
#   - Do not delete results data:
#       ./run_test_suite -d
#   - Run in verbose mode:
#       ./run_test_suite -v
#   - Select test number (for test number 3 to 5, and 8):
#       ./run_test_suite -t 3-5,8
#   - Use batch submission :
#       ./run_test_suite -q "qsub"
#
#######################################################################

#######################################################################
# Determine the parameters for running the test suite
#######################################################################
MPI=0;
NCPU=1;
USE_PYTHON=0;
RUN_FULL_SUITE=false;
VERBOSE=false;
DELDATA=true;
SELECTTEST=false;
CLEAN_ALL=false;
QUEUE=false;
while getopts "cdfp:q:t:vy" OPTION; do
   case $OPTION in
      c)
         CLEAN_ALL=true;
      ;;
      d)
         DELDATA=false;
      ;;
      f)
         RUN_FULL_SUITE=true;
      ;;
      p)
         MPI=1;
         NCPU=$OPTARG;
      ;;
      q)
         QUEUE=true;
         SUBMITSTRING=$OPTARG;
      ;;
      t)
         SELECTTEST=true;
         TESTNUMBER=$OPTARG;
      ;;
      v)
         VERBOSE=true;
      ;;
      y)
         USE_PYTHON=1;
      ;;
   esac
done

MYLATEXPATH="/gpfs/data1/nvaytet/texlive2015/bin/x86_64-linux/";
MYMODULES="gnuplot anaconda openmpi/1.8.3-gnu4.8.4";

#######################################################################
# Setup paths and commands
#######################################################################
TEST_DIRECTORY=$(pwd);

length=${#TEST_DIRECTORY};
icut=$(($length - 17));

BASE_DIRECTORY="${TEST_DIRECTORY:0:${icut}}";
BIN_DIRECTORY="${BASE_DIRECTORY}/bin";
PATCH_DIRECTORY="${BASE_DIRECTORY}/patch/rhd";

VISU_DIR="${TEST_DIRECTORY}/visu";
DELETE_RESULTS="rm -rf output_* data*.dat time.dat";
DELETE_SOURCES="rm -f units.o condinit.o";
RETURN_TO_BIN="cd ${BIN_DIRECTORY}";
EXECNAME="test_exe_";
LOGFILE="${TEST_DIRECTORY}/test_suite.log";
echo -n > $LOGFILE;
COMPLETEDTESTS="${TEST_DIRECTORY}/completed_tests.txt";
if $QUEUE; then
   echo -n > $COMPLETEDTESTS;
fi
nbatch=5;
#COMP=$(grep COMP ${BIN_DIRECTORY}/Makefile | grep "=" | cut -d '"' -f2);

if [ ${MPI} -eq 1 ]; then
#    MPIFLAGS="--default-hostfile none"
   RUN_TEST_BASE="mpirun -np ${NCPU} ${MPIFLAGS}";
else
   RUN_TEST_BASE="./";
fi

TOTALSTARTTIME=$(date +%s);

#######################################################################
# Clean all directories and exit
#######################################################################
if $CLEAN_ALL ; then
   cd ${TEST_DIRECTORY};
   rm -r */output_* */*.log */data*.dat */time.dat */resdiff* */*.pdf */*.ps */*.o */*.mod */build* */ramses_* completed_tests.txt test_results.* test_suite.log */barotropic_eos.dat */Hosokawa_track.dat */vaytet_grey_opacities*.bin */groups.dat */init_turb.data */tab_eos.dat */res*.dat */submit*.sh;
   cd ${BIN_DIRECTORY};
   make clean;
   rm ${EXECNAME}*;
   cd ${VISU_DIR};
   make clean;
   exit;
fi

#######################################################################
# Welcome message
#######################################################################
echo "############################################";
echo "############################################" >> $LOGFILE;
if $RUN_FULL_SUITE ; then
   echo "Running extended RAMSES automatic test suite";
   echo "Running extended RAMSES automatic test suite" >> $LOGFILE;
else
   echo "Running standard RAMSES automatic test suite";
   echo "Running standard RAMSES automatic test suite" >> $LOGFILE;
fi
echo "############################################";
echo "############################################" >> $LOGFILE;

#######################################################################
# List of tests
#######################################################################

itest=0; # Test 1
testname[${itest}]="sod-tube";
testpatch[${itest}]="";
ndim[${itest}]=1;
flags[${itest}]="";

itest=$((itest + 1)); # Test 2
testname[${itest}]="imhd-tube";
testpatch[${itest}]="";
ndim[${itest}]=1;
flags[${itest}]="";

itest=$((itest + 1)); # Test 3
testname[${itest}]="nimhd-diffusion-ad";
testpatch[${itest}]="../patch/test_suite/nimhd-diffusion-ad";
ndim[${itest}]=3;
flags[${itest}]="NIMHD=1";

itest=$((itest + 1)); # Test 4
testname[${itest}]="nimhd-diffusion-ohm";
testpatch[${itest}]="../patch/test_suite/nimhd-diffusion-ohm";
ndim[${itest}]=3;
flags[${itest}]="NIMHD=1";

itest=$((itest + 1)); # Test 5
testname[${itest}]="dirac";
testpatch[${itest}]="../patch/test_suite/dirac";
ndim[${itest}]=1;
flags[${itest}]="NGRP=1 USE_FLD=1";

itest=$((itest + 1)); # Test 6
testname[${itest}]="rshock-Mach2";
testpatch[${itest}]="";
ndim[${itest}]=1;
flags[${itest}]="NGRP=1 USE_FLD=1";

itest=$((itest + 1)); # Test 7
testname[${itest}]="rshock-Mach5";
testpatch[${itest}]="";
ndim[${itest}]=1;
flags[${itest}]="NGRP=1 USE_FLD=1";

itest=$((itest + 1)); # Test 8
testname[${itest}]="radiative-shock";
testpatch[${itest}]="";
ndim[${itest}]=1;
flags[${itest}]="NGRP=4 USE_FLD=1";

itest=$((itest + 1)); # Test 9
testname[${itest}]="orszag-tang";
testpatch[${itest}]="../patch/test_suite/orszag-tang";
ndim[${itest}]=2;
flags[${itest}]="";

# Store number of standard tests
ntestsstandard=${#testname[@]};

# Additional tests ==========================

itest=$((itest + 1)); # Test 10
testname[${itest}]="dirac3d";
testpatch[${itest}]="../patch/test_suite/dirac3d";
ndim[${itest}]=3;
flags[${itest}]="NGRP=1 USE_FLD=1";

itest=$((itest + 1)); # Test 11
testname[${itest}]="nimhd-cshock-ad";
testpatch[${itest}]="";
ndim[${itest}]=3;
flags[${itest}]="NIMHD=1";

itest=$((itest + 1)); # Test 12
testname[${itest}]="nimhd-cshock-ohm";
testpatch[${itest}]="";
ndim[${itest}]=3;
flags[${itest}]="NIMHD=1";

itest=$((itest + 1)); # Test 13
testname[${itest}]="collapse-rhd";
testpatch[${itest}]="../patch/collapse";
ndim[${itest}]=3;
flags[${itest}]="NGRP=1 USE_FLD=1";

itest=$((itest + 1)); # Test 14
testname[${itest}]="collapse-baro";
testpatch[${itest}]="../patch/collapse";
ndim[${itest}]=3;
flags[${itest}]="NIMHD=1";

itest=$((itest + 1)); # Test 15
testname[${itest}]="collapse-ohm";
testpatch[${itest}]="../patch/collapse";
ndim[${itest}]=3;
flags[${itest}]="NGRP=1 USE_FLD=1 NIMHD=1";

ntestsfull=${#testname[@]};

# Count number of tests
if ${RUN_FULL_SUITE} ; then
   ntestsall=${ntestsfull};
else
   ntestsall=${ntestsstandard};
fi

ntests=$ntestsall;
all_tests_ok=true;

#######################################################################
# Select particular test if this was asked by user
#######################################################################
if $SELECTTEST ; then

   # Split test selection with commas
   s1=$(echo $TESTNUMBER | sed 's/,/ /'g);
   testsegs=( $s1 );
   nseg=${#testsegs[@]};
   
   # Search for dashes in individual segments
   ntests=0;
   for ((n=0;n<$nseg;n++)); do
      dashsearch=$(echo ${testsegs[n]} | grep '-');
      if [ ${#dashsearch} -gt 0 ] ; then
         istart=$(echo ${testsegs[n]} | cut -d '-' -f1);
         iend=$(echo ${testsegs[n]} | cut -d '-' -f2);
         is=$((istart - 1));
         ie=$((iend - 1));
         iep1=$(($ie + 1));
         for ((j=$is;j<$iep1;j++)); do
            if [ ${j} -ge 0 ] && [ ${j} -lt $ntestsfull ] ; then
               testnum[${ntests}]=$j;
               ntests=$((ntests + 1));
            else
               echo "Selected test ${j} does not exist! Ignoring test";
               echo "Selected test ${j} does not exist! Ignoring test" >> $LOGFILE;
            fi
         done
      else
         # No dash, just include test in list
         testnum[${ntests}]=$((${testsegs[n]} - 1));
         if [ ${testnum[${ntests}]} -gt $ntestsfull ] ; then
            echo "Selected test does not exist!";
            echo "Selected test does not exist!" >> $LOGFILE;
            exit;
         fi
         ntests=$((ntests + 1));
      fi
   done

else

   # Include all tests by default
   for ((n=0;n<$ntests;n++)); do
      testnum[n]=$n;
   done
   
fi

#######################################################################
# Write list of tests
#######################################################################
echo "Will perform the following tests:";
echo "Will perform the following tests:" >> $LOGFILE;
for ((i=0;i<$ntests;i++)); do
   n=${testnum[i]};
   j=$(($i + 1));
   if [ $ntests -gt 9 ] && [ $j -lt 10 ] ; then
      echo " [ ${j}] ${testname[n]}";
      echo " [ ${j}] ${testname[n]}" >> $LOGFILE;
   else
      echo " [${j}] ${testname[n]}";
      echo " [${j}] ${testname[n]}" >> $LOGFILE;
   fi
done
echo "--------------------------------------------";
echo "--------------------------------------------" >> $LOGFILE;

#######################################################################
# Prepare visualization software
#######################################################################
echo "Compiling visualization software";
echo "Compiling visualization software" >> $LOGFILE;
cd ${VISU_DIR};
if $VERBOSE ; then
   make;
else
   make >> $LOGFILE 2>&1;
fi
echo "--------------------------------------------";
echo "--------------------------------------------" >> $LOGFILE;

#######################################################################
# Loop through all tests
#######################################################################
itest=0;
count1=0;
count2=1;
for ((i=0;i<$ntests;i++)); do

   n=${testnum[i]};
   itest=$(($itest + 1));
   
   THISTESTLOG="${TEST_DIRECTORY}/${testname[n]}/${testname[n]}.log";
   echo -n > ${THISTESTLOG};

   # Check if we are using batch queue or not
   if $QUEUE; then

      # Pause if we are submitting in batches
      if [ $nbatch -gt 0 ] && [ $count1 -ge $nbatch ] ; then
         waitingforbatch=true;
         while $waitingforbatch ; do
            ncomp=$(wc -l $COMPLETEDTESTS | cut -d ' ' -f1);
            if [ $ncomp -ge $count2 ]; then
               waitingforbatch=false;
#               count2=$(($count2 + $nbatch));
               count2=$(($count2 + 1));
            fi
            sleep 1;
         done
#         count1=0;
         #count1=$(($count1 - 1));
      fi

      echo "Test ${itest}/${ntests}: ${testname[n]}";

      count1=$(($count1 + 1));
   
      SUBFILE="submit_${testname[n]}.sh";
      
      cd ${TEST_DIRECTORY}/${testname[n]};
      
      # Write submission script
      echo "#!/bin/bash" > $SUBFILE;
      echo "#$ -N ${testname[n]}" >> $SUBFILE;
      echo "#$ -o ${THISTESTLOG}" >> $SUBFILE;
      echo "#$ -j y" >> $SUBFILE;
      echo "#$ -cwd" >> $SUBFILE;
      echo "#$ -pe mpi ${NCPU}" >> $SUBFILE;
      echo "STARTTIME=\$(date +%s);" >> $SUBFILE;
      echo "source /etc/profile.d/modules.sh" >> $SUBFILE;
      echo "module purge" >> $SUBFILE;
      echo "module load ${MYMODULES}" >> $SUBFILE;
#      if [ $COMP == "GNU" ] ; then
#         echo "module load gnuplot" >> $SUBFILE;
#         echo "module load gnuplot openmpi/1.8.3-gnu4.8.4" >> $SUBFILE;
#       echo "module load "
#      elif [ $COMP == "INTEL" ] ; then
#         echo "module load openmpi/1.5.4-intel12.0.0" >> $SUBFILE;
#      else
#         echo "WARNING: unknown compiler";
#      fi
      echo "echo \"Test ${itest}/${ntests}: ${testname[n]}\" >> ${THISTESTLOG};" >> $SUBFILE;
      
      # Initial cleanup
      echo "mkdir build-${testname[n]};" >> $SUBFILE;
      echo "cd build-${testname[n]};" >> $SUBFILE;
      echo "cp ${BIN_DIRECTORY}/Makefile .;" >> $SUBFILE;
      BINSTRING=$(echo ${BIN_DIRECTORY} | sed 's/\//\\\//g');
      echo "sed -i 's/BINDIR = \./BINDIR = ${BINSTRING}/g' Makefile;" >> $SUBFILE;
      echo "echo \"Cleanup\";" >> $SUBFILE;
      echo "make clean;" >> $SUBFILE;

      # Compile source
      echo "echo \"Compiling source\";" >> $SUBFILE;
      echo "make EXEC=${EXECNAME} PATCH=${testpatch[n]} MPI=${MPI} NDIM=${ndim[n]} ${flags[n]};" >> $SUBFILE;
      echo "cd ../;" >> $SUBFILE;
      
      # Run tests
      echo "$DELETE_RESULTS;" >> $SUBFILE;
      echo "./prepare-${testname[n]}.sh;" >> $SUBFILE;
      echo "mv build-${testname[n]}/${EXECNAME}${ndim[n]}d ramses_${testname[n]};" >> $SUBFILE;
      
      echo "echo \"Running test\";" >> $SUBFILE;
      
#      RUN_TEST="${RUN_TEST_BASE}ramses_${testname[n]} ${testlist[n]}";     
#      echo "{ time { ${RUN_TEST}; } 2>> ${THISTESTLOG}; } 2> ${testname[n]}_stats.txt;" >> $SUBFILE;
      echo "${RUN_TEST_BASE}ramses_${testname[n]} ${testname[n]}.nml;" >> $SUBFILE;
      
      # Plot results
      echo "echo \"Plotting results\";" >> $SUBFILE;
      echo "./plot-${testname[n]}.sh ${USE_PYTHON};" >> $SUBFILE;
  #    if ${USE_GNUPLOT} ; then
  #       echo "module load gnuplot" >> $SUBFILE;
  #       echo "gnuplot plot-${testname[n]}.gp;" >> $SUBFILE;
  #       echo "ps2pdf ${testname[n]}.ps;" >> $SUBFILE;
  #       echo "rm ${testname[n]}.ps;" >> $SUBFILE;
  #    fi
  #    if ${USE_PYTHON} ; then
  #       echo "module load anaconda" >> $SUBFILE;
  #       echo "python plot-${testname[n]}.py;" >> $SUBFILE;
  #    fi
      
      # Check for differences in results
      echo "echo \"Analysing results\";" >> $SUBFILE;
      echo "difffile=\"resdiff-${testname[n]}\";" >> $SUBFILE;
      echo "diff data.dat ${testname[n]}\"-ref.dat\" >> \${difffile};" >> $SUBFILE;
      # Size of diff file?
      echo "diffoutput=\$(cat \${difffile});" >> $SUBFILE;
      echo "diffsize=\${#diffoutput};" >> $SUBFILE;
      echo "if [ \${diffsize} -gt 0 ]; then" >> $SUBFILE;
      echo "   echo \"Test failed!                          [FAIL]\";" >> $SUBFILE;
      echo "   all_tests_ok=false;" >> $SUBFILE;
      echo "else" >> $SUBFILE;
      echo "   echo \"Test passed                           [ OK ]\";" >> $SUBFILE;
      echo "fi" >> $SUBFILE;
      
      # Let script know that test has finished
      echo "echo \"${testname[n]}\" >> $COMPLETEDTESTS;" >> $SUBFILE;
      
      # Report elapsed time
      echo "ENDTIME=\$(date +%s);" >> $SUBFILE;
      echo "seconds=\$((\$ENDTIME - \$STARTTIME));" >> $SUBFILE;
      echo "echo \"Test duration: \$seconds s\" >> ${THISTESTLOG};" >> $SUBFILE;

      # Submit script in batch system
      ${SUBMITSTRING} $SUBFILE;

   else
   
      echo "Test ${itest}/${ntests}: ${testname[n]}";
      echo "Test ${itest}/${ntests}: ${testname[n]}" >> ${THISTESTLOG};

      STARTTIME=$(date +%s);

      # Initial cleanup
      $RETURN_TO_BIN;
      echo "Cleanup";
      echo "Cleanup" >> ${THISTESTLOG};
      if $VERBOSE ; then
         make clean 2>&1 | tee -a ${THISTESTLOG};
      else
         make clean >> ${THISTESTLOG} 2>&1;
      fi
      
      # Compile source
      echo "Compiling source";
      echo "Compiling source" >> ${THISTESTLOG};
      MAKESTRING="make EXEC=${EXECNAME} PATCH=${testpatch[n]} MPI=${MPI} NDIM=${ndim[n]} ${flags[n]}";
      if $VERBOSE ; then
         $MAKESTRING 2>&1 | tee -a ${THISTESTLOG};
      else
         $MAKESTRING >> ${THISTESTLOG} 2>&1;
      fi
      
      # Run tests
      cd ${TEST_DIRECTORY}/${testname[n]};
      $DELETE_RESULTS;
      if $VERBOSE ; then
         ./prepare-${testname[n]}.sh 2>&1 | tee -a ${THISTESTLOG};
      else
         ./prepare-${testname[n]}.sh >> ${THISTESTLOG} 2>&1;
      fi
      mv ${BIN_DIRECTORY}/${EXECNAME}${ndim[n]}d ramses_${testname[n]};
      
      RUN_TEST="${RUN_TEST_BASE}ramses_${testname[n]} ${testname[n]}.nml";
      echo "Running test";
      echo "Running test" >> ${THISTESTLOG};
   
      if $VERBOSE ; then
         ${RUN_TEST} 2>&1 | tee -a ${THISTESTLOG};
      else
         ${RUN_TEST} >> ${THISTESTLOG} 2>&1;
      fi

      # Plot results
      echo "Plotting results";
      echo "Plotting results" >> ${THISTESTLOG};
      if $VERBOSE ; then
         ./plot-${testname[n]}.sh ${USE_PYTHON} 2>&1 | tee -a ${THISTESTLOG};
      else
         ./plot-${testname[n]}.sh ${USE_PYTHON} >> ${THISTESTLOG} 2>&1;
      fi
#      if ${USE_GNUPLOT} ; then
#         gnuplot plot-${testname[n]}.gp;
#         ps2pdf ${testname[n]}.ps;
#         rm ${testname[n]}.ps;
#      fi
#      if ${USE_PYTHON} ; then
#         python plot-${testname[n]}.py;
#      fi
      
      # Check for differences in results
      echo "Analysing results";
      echo "Analysing results" >> ${THISTESTLOG};
      difffile="resdiff-${testname[n]}";
      diff data.dat ${testname[n]}"-ref.dat" &> ${difffile};
      # Size of diff file?
      diffoutput=$(cat ${difffile});
      diffsize=${#diffoutput};
      if [ ${diffsize} -gt 0 ]; then
         echo "Test failed!                          [FAIL]";
         echo "Test failed!                          [FAIL]" >> ${THISTESTLOG};
      else
         echo "Test passed                           [ OK ]";
         echo "Test passed                           [ OK ]" >> ${THISTESTLOG};
      fi
      
      ENDTIME=$(date +%s);
      seconds=$(($ENDTIME - $STARTTIME));
      echo "Test duration: $seconds s" >> ${THISTESTLOG};

      echo "--------------------------------------------";
      
   fi
   
#    # Check for differences in log files
#    echo "Analysing logs"; echo "Analyzing logs" >> $LOGFILE;
#    # Remove grids and memory from log file as they change with MPI
#    sed -i 's/has.*$//' log;
#    sed -i 's/mem.*$//' log;
#    # Apply difference (ignoring elapsed times)
#    diff log ${testname[n]}.log -I elapsed > ${testname[n]}"_logdiff.tex";
#    # Size of diff file?
#    diffoutput=$(cat ${testname[n]}"_logdiff.tex");
#    diffsize=${#diffoutput};
#    if [ ${diffsize} -gt 0 ]; then
#       diff_not_empty[n]=true;
#       echo "Test failed!"; echo "Test failed!" >> $LOGFILE;
#       all_tests_ok=false;
#       # Format diff output
#       sed -i 's/</\$<\$/g' ${testname[n]}"_logdiff.tex";
#       sed -i 's/>/\$>\$/g' ${testname[n]}"_logdiff.tex";
#       sed -i 's/%/\\%/g' ${testname[n]}"_logdiff.tex";
#       sed -i 's/$/ \\\\/g' ${testname[n]}"_logdiff.tex";
#    else
#       diff_not_empty[n]=false;
#       echo "Test passed"; echo "Test passed" >> $LOGFILE;
#    fi
   
#    echo "--------------------------------------------";
#    echo "--------------------------------------------" >> $LOGFILE;

done

#######################################################################
cd ${TEST_DIRECTORY};

if $QUEUE; then
   notcompleted=true;
   while ${notcompleted} ; do
      ncomp=$(wc -l $COMPLETEDTESTS | cut -d ' ' -f1);
      if [ ${ncomp} -eq ${ntests} ]; then
         notcompleted=false;
      fi
      sleep 5;
   done
fi

TOTALENDTIME=$(date +%s);
totalseconds=$((${TOTALENDTIME} - ${TOTALSTARTTIME}));
totalhours=$((${totalseconds} / 3600));
totalseconds=$((${totalseconds} % 3600));
totalminutes=$((${totalseconds} / 60));
totalseconds=$((${totalseconds} % 60));

for ((i=0;i<${ntests};i++)); do
   n=${testnum[i]};
   cat ${TEST_DIRECTORY}/${testname[n]}/${testname[n]}.log >> $LOGFILE;
   echo "--------------------------------------------" >> $LOGFILE;
done
#######################################################################

#######################################################################
# Generate pdf document with test results
#######################################################################
echo "Generating pdf document with test results";
echo "Generating pdf document with test results" >> $LOGFILE;
latexfile="test_results.tex";
echo "\documentclass[12pt]{article}" > $latexfile;
echo "\usepackage{graphicx,color}" >> $latexfile;
echo "\usepackage[colorlinks=true,linkcolor=blue]{hyperref}" >> $latexfile;
echo "\topmargin -1.3in" >> $latexfile;
echo "\textheight 10.1in" >> $latexfile;
echo "\oddsidemargin -0.7in" >> $latexfile;
echo "\evensidemargin -0.7in" >> $latexfile;
echo "\textwidth 7.7in" >> $latexfile;
echo >> $latexfile;
echo "\title{RAMSES test suite results}" >> $latexfile;
echo "\date{\today}" >> $latexfile;
echo "\author{${USER}}" >> $latexfile;
echo >> $latexfile;
echo "\begin{document}" >> $latexfile;
echo >> $latexfile;
echo "\maketitle" >> $latexfile;
echo >> $latexfile;
echo "\begin{table}[ht]" >> $latexfile;
echo "\centering" >> $latexfile;
echo "\caption{Test run summary using ${NCPU} processor(s)}" >> $latexfile;
echo "\begin{tabular}{|r|l|l|l|}" >> $latexfile;
echo "\hline" >> $latexfile;
echo "~ & Test name & Duration & Status\\\\" >> $latexfile;
echo "\hline" >> $latexfile;
for ((i=0;i<$ntests;i++)); do
   n=${testnum[i]};
   logfile="${TEST_DIRECTORY}/${testname[n]}/${testname[n]}.log";
   itest=$(($i + 1));
   checkfail=$(grep "\[FAIL\]" ${logfile});
   if [ ${#checkfail} -gt 0 ]; then
      status="\hyperref[fig-${testname[n]}]{\textcolor{red}{failed}}";
      all_tests_ok=false;
   else
      status="\hyperref[fig-${testname[n]}]{\textcolor{green}{passed}}";
   fi
   seconds=$(grep "Test duration:" ${logfile} | cut -d ' ' -f3);
   hours=$(($seconds / 3600));
   seconds=$(($seconds % 3600));
   minutes=$(($seconds / 60));
   seconds=$(($seconds % 60));
   echo $itest "& \hyperref[fig-${testname[n]}]{${testname[n]}} & ${hours}h${minutes}m${seconds}s & ${status} \\\\" >> $latexfile;
done
echo "\hline" >> $latexfile;
echo "\end{tabular}" >> $latexfile;
echo "\end{table}" >> $latexfile;
echo "\begin{center}" >> $latexfile;
echo "Total run time (including compilations): ${totalhours}h${totalminutes}m${totalseconds}s" >> $latexfile;
echo "\end{center}" >> $latexfile;
echo "\clearpage" >> $latexfile;
echo >> $latexfile;

for ((i=0;i<$ntests;i++)); do
   n=${testnum[i]};
   echo "\begin{figure}" >> $latexfile;
   echo "\centering" >> $latexfile;
   echo "\includegraphics[scale=0.7]{${TEST_DIRECTORY}/${testname[n]}/${testname[n]}.pdf}" >> $latexfile;
   echo "\caption{${testname[n]} test}" >> $latexfile;
   echo "\label{fig-${testname[n]}}" >> $latexfile;
   echo "\end{figure}" >> $latexfile;
   echo "\clearpage" >> $latexfile;
   echo >> $latexfile;
#    if ${diff_not_empty[n]} ; then
#       echo "{\bf Differences in log for test ${testname[n]}:}\\\\" >> $latexfile;
#       echo "\input{${TEST_DIRECTORY}/${testname[n]}/${testname[n]}_logdiff.tex}" >> $latexfile;
#       echo "\clearpage" >> $latexfile;
#       echo >> $latexfile;
#    fi
done
echo "\end{document}" >> $latexfile;
if $VERBOSE ; then
   ${MYLATEXPATH}pdflatex $latexfile 2>&1 | tee -a $LOGFILE;
   ${MYLATEXPATH}pdflatex $latexfile 2>&1 | tee -a $LOGFILE;
else
   ${MYLATEXPATH}pdflatex $latexfile >> $LOGFILE 2>&1;
   ${MYLATEXPATH}pdflatex $latexfile >> $LOGFILE 2>&1;
fi
if ${DELDATA} ; then
   rm -f ${latexfile/.tex/.log};
   rm -f ${latexfile/.tex/.aux};
   rm -f ${latexfile/.tex/.out};
   rm -f $latexfile;
fi

if $all_tests_ok ; then
   echo "All tests were completed successfully";
   echo "All tests were completed successfully" >> $LOGFILE;
else
   echo "There were some failed tests";
   echo "There were some failed tests" >> $LOGFILE;
fi

#######################################################################
# Clean up
#######################################################################

if ${DELDATA} ; then
   for ((i=0;i<$ntests;i++)); do
      n=${testnum[i]};
      cd ${TEST_DIRECTORY}/${testname[n]};
      $DELETE_RESULTS;
      rm -rf ${testname[n]}_stats.txt ${testname[n]}.pdf ${testname[n]}.log build-${testname[n]} ramses_${testname[n]} submit_${testname[n]}.sh;
   done
   rm -f $COMPLETEDTESTS;
   if $VERBOSE ; then
      cd ${VISU_DIR}
      make clean 2>&1 | tee -a ${LOGFILE};
      if [ !${QUEUE} ] ; then
         $RETURN_TO_BIN;
         make clean 2>&1 | tee -a ${LOGFILE};
      fi
   else
      cd ${VISU_DIR}
      make clean >> $LOGFILE 2>&1;
      if [ !${QUEUE} ] ; then
         $RETURN_TO_BIN;
         make clean >> $LOGFILE 2>&1;
      fi
   fi
   rm -f ${EXECNAME}*d;
fi

exit;
