rm *.o
rm *.mod
exit

gfortran -c -Wall RTSMie_Parameters.f90 
gfortran -c -Wall RTSMie_Inputs_Def.f90
gfortran -c -Wall RTSMie_Outputs_Def.f90
gfortran -c -Wall RTSMie_Lin_Outputs_Def.f90
gfortran -c -Wall RTSMie_Distribution.f90
gfortran -c -Wall RTSMie_IO_Readwrite.f90
gfortran -c -Wall RTSMie_Master.f90
gfortran -c -Wall RTSMie_Master_Bimodal.f90
gfortran -c -Wall RTSMie_Master_PLUS.f90
gfortran -c -Wall RTSMie_Master_Bimodal_PLUS.f90
exit

