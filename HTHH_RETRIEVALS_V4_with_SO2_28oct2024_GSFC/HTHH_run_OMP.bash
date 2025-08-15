#!/bin/bash

# This script runs the HTHH retrieval tool which has been set up to run using OpenMP
# (see http://openmp.org for further information about OpenMP)

# Sample commands:
# HTHH_run_OMP.bash gfortran
# HTHH_run_OMP.bash ifort


# Do some set up

ulimit -s unlimited         # Maximum stack size of the main thread
#ulimit -c unlimited         # Maximum size of core file created if a problem occurs

#export OMP_STACKSIZE=50M    # Maximum stack size of OpenMP-spawned threads
#export OMP_STACKSIZE=100M
#export OMP_STACKSIZE=150M
#export OMP_STACKSIZE=200M
#export OMP_STACKSIZE=512M #For runs runs using <=139 wvls
export OMP_STACKSIZE=600M

echo
echo 'OMP_STACKSIZE = '$OMP_STACKSIZE
echo
#exit 0

# Run HTHH retrieval using OpenMP

echo
echo "making OpenMP version of HTHH retrieval tool ..."
echo

# 2-band retrieval
#make Retrieval_MeasRatio_Mk2_with_SO2_V4_OMP.exe FC=$1 COMPILE_OMP=t DEBUG=t         # serial dbg compile
#make Retrieval_MeasRatio_Mk2_with_SO2_V4_OMP.exe FC=$1 COMPILE_OMP=t OPT=t           # serial opt compile

#make Retrieval_MeasRatio_Mk2_with_SO2_V4_OMP.exe FC=$1 COMPILE_OMP=t OPENMP=t        # OpenMP norm compile
make Retrieval_MeasRatio_Mk2_with_SO2_V4_OMP.exe FC=$1 COMPILE_OMP=t OPENMP=t OPT=t  # OpenMP opt compile

#echo
#echo 'at HTHH_run_OMP.bash hold/exit'
#echo
#exit

echo
echo "running Retrieval_MeasRatio_Mk2_with_SO2_V4_OMP.exe ..."
echo
./Retrieval_MeasRatio_Mk2_with_SO2_V4_OMP.exe
echo

#rm Retrieval_MeasRatio_Mk2_with_SO2_V4_OMP.exe
#rm makefile

echo
echo 'done'
echo

exit 0
