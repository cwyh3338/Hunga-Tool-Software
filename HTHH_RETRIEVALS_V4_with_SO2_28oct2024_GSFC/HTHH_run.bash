#!/bin/bash

# This script runs the serial version of the HTHH retrieval tool

# Sample commands:
# HTHH_run.bash gfortran
# HTHH_run.bash ifort


# Do some set up

ulimit -s unlimited         # Maximum stack size of the main thread
#ulimit -c unlimited         # Maximum size of core file created if a problem occurs

# Run HTHH retrieval

echo
echo "making serial version of HTHH retrieval tool ..."
echo

# 2-band retrieval
#make Retrieval_MeasRatio_Mk2_with_SO2_V4.exe FC=$1 DEBUG=t  # serial dbg compile
#make Retrieval_MeasRatio_Mk2_with_SO2_V4.exe FC=$1          # serial normal compile
make Retrieval_MeasRatio_Mk2_with_SO2_V4.exe FC=$1 OPT=t    # serial opt compile

#echo
#echo 'at HTHH_run.bash hold/exit'
#echo
#exit 0

echo
echo "running Retrieval_MeasRatio_Mk2_with_SO2_V4.exe ..."
echo
./Retrieval_MeasRatio_Mk2_with_SO2_V4.exe
echo

#rm Retrieval_MeasRatio_Mk2_with_SO2_V4.exe
#rm makefile

echo
echo 'done'
echo

exit 0
