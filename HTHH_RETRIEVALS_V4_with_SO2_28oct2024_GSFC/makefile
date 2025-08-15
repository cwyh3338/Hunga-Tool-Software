
#   MAKEFILE for the HUNGA-TONGA Retrievals, March-June 2022
#      R. Spurr, RT Solutions Inc.

#  -- Uses VLIDORT in LCS Mode, Version 2.8.3, but with Local PARS file for dimensioning
#  -- Uses RTSMIE linearized Mie Code, Version 1.7 (I/O type structures, otherwise Verison 1.5)
#  -- Dedicated subdirectory of TONGA-related Sourcecode.

#  VLIDORT UPGRADES
#  ================

# UPGRADE TO VLIDORT VERSION 2.7, May   2014
# UPGRADE TO VLIDORT VERSION 2.8, March 2017, together with FO code Version 1.5
# UPGRADE TO VLIDORT VERSION 2.8.1, August 2019
#      4/9/19.   Removal of TRANSFLUX MASTER Routine
#      4/19/19.  Introduction of mediaprops routines (3 in all)
#      5/1/19.   Supplement and FO codes no compilation changes
#      6/3/19.   Introduce 2 new special tests (Planetary/LWCoupling)
# UPGRADE TO VLIDORT VERSION 2.8.3, 31 March 2021

#  Define some variables
#  =====================

#  VLIDORT paths

#VPATH           = ../../../../RTSMODELS_2022/VLIDORT/Version_2.8.3_FINAL/RTS_VLIDORT2p8p3_f90pack_18aug2021
VPATH           = RTS_VLIDORT2p8p3_f90pack_31oct2022

VSUP_PATH       = $(VPATH)/vsup
FO_MAIN_PATH    = $(VPATH)/vlidort_focode
VLID_DEF_PATH   = $(VPATH)/vlidort_def
VLID_MAIN_PATH1 = $(VPATH)/vlidort_main/regular
VLID_MAIN_PATH2 = $(VPATH)/vlidort_main/linearized

#  Mie Path

# MIE_PATH  = ../../../../RTSMODELS_2022/MIE_TMAT/RAY_MIE_TMAT_2021/RTSMie_1p7_2022/RTSMie_sourcecode/
MIE_PATH  = RTSMie_sourcecode_V1p7

# VLIDORT test path

VLID_TEST_PATH = $(VPATH)/vlidort_s_test

# TONGA paths

TONGA_TEST_PATH    = .
TONGA_SOURCES_PATH = TONGA_Sourcecode

# Define module & object subdirectories

ifeq ($(COMPILE_OMP), t)
  MOD_PATH = mod_omp
  OBJ_PATH = obj_omp
else
  MOD_PATH = mod
  OBJ_PATH = obj
endif

MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

# Define default shell make will use

SHELL = /bin/bash

# Define FORTRAN90 compiler to use (can be defined here, but usually defined on makefile command line)

#  Intel
# FC = ifort

#  gfortran
#FC = gfortran

#  g95
# FC = g95

#  NAG
# FC = f95
# FC = nagfor

#
# Define FORTRAN90 compiler flags
#

FFLAGS = -c

# Additional flags for Intel
ifeq ($(FC), ifort)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -module $(MOD_PATH)
	FFLAGS_DEBUG = -g -warn all -check all -traceback
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -fbacktrace
	FFLAGS_DEBUG = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -fopenmp -frecursive
endif

# Additional flags for g95
ifeq ($(FC), g95)
#      older g95
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check
	FFLAGS_DEBUG = -g -Wall -Wuninitialized -fbounds-check
	FFLAGS_OPT = -O3
#      g95 v0.92
#	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -ftrace=full
endif

# Additional flags for NAG
#ifeq ($(FC), f95)
ifeq ($(FC), $(filter $(FC), f95 nagfor))
	FFLAGS := $(FFLAGS) -mdir $(MOD_PATH) -I$(MOD_PATH)
#	FFLAGS_DEBUG = -w=obs -w=unused -C=array -C=undefined -gline
	FFLAGS_DEBUG = -w=obs -w=unused -C=all -C=undefined -gline
#	FFLAGS_DEBUG = -g -C=all -C=undefined -gline -mtrace=all -nan
#	FFLAGS_DEBUG = -w=obs -w=unused -g -C=all -C=undefined -gline -mtrace=all -nan
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

#  for compilation detailed
#	FFLAGS := $(FFLAGS) $(FFLAGS_Catch)

# For debug build, use "make DEBUG=t"
ifeq ($(DEBUG), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG)
endif

# For optimized build, use "make OPT=t"
ifeq ($(OPT), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPT)
endif

# For parallel build using OpenMP, use "make OPENMP=t"
ifeq ($(OPENMP), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif

.SUFFIXES:

# Define list of source files
# (Note: ordering is important because of dependencies)

BASE_SOURCES =
SOURCES =
L_SOURCES =
LCS_SOURCES =

BASE_SOURCES +=   \
    vlidort_pars_local.f90

#BASE_SOURCES +=   \
#    vlidort_testpts.f90

SOURCES +=   \
   $(BASE_SOURCES) \
   $(VLID_MAIN_PATH1)/lapack_tools.f90           \
   $(VLID_DEF_PATH)/vlidort_inputs_def.f90       \
   $(VLID_DEF_PATH)/vlidort_sup_brdf_def.f90     \
   $(VLID_DEF_PATH)/vlidort_sup_ss_def.f90       \
   $(VLID_DEF_PATH)/vlidort_sup_sleave_def.f90   \
   $(VLID_DEF_PATH)/vlidort_sup_def.f90          \
   $(VLID_DEF_PATH)/vlidort_outputs_def.f90      \
   $(VLID_DEF_PATH)/vlidort_io_defs.f90          \
   $(VLID_DEF_PATH)/vlidort_work_def.f90         \
   $(VLID_MAIN_PATH1)/vlidort_aux.f90            \
   $(VLID_MAIN_PATH1)/vlidort_getplanck.f90      \
   $(VLID_MAIN_PATH1)/vlidort_geometry.f90       \
   $(VLID_MAIN_PATH1)/vlidort_Taylor.f90         \
   $(VLID_MAIN_PATH1)/vlidort_inputs.f90         \
   $(VLID_MAIN_PATH1)/vlidort_miscsetups.f90     \
   $(VLID_MAIN_PATH1)/vlidort_vfo_interface.f90  \
   $(VLID_MAIN_PATH1)/vlidort_thermalsup.f90     \
   $(VLID_MAIN_PATH1)/vlidort_solutions.f90      \
   $(VLID_MAIN_PATH1)/vlidort_bvproblem.f90      \
   $(VLID_MAIN_PATH1)/vlidort_PostProcessing.f90 \
   $(VLID_MAIN_PATH1)/vlidort_intensity.f90      \
   $(VLID_MAIN_PATH1)/vlidort_converge.f90       \
   $(VLID_MAIN_PATH1)/vlidort_writemodules.f90   \
   $(VLID_MAIN_PATH1)/vlidort_pack.f90           \
   $(VLID_MAIN_PATH1)/vlidort_unpack.f90         \
   $(VLID_MAIN_PATH1)/vlidort_mediaprops.f90     \
   $(VLID_MAIN_PATH1)/vlidort_masters_V2p8p3.f90

L_SOURCES += \
   $(VLID_DEF_PATH)/vlidort_lin_inputs_def.f90    \
   $(VLID_DEF_PATH)/vlidort_lin_sup_brdf_def.f90  \
   $(VLID_DEF_PATH)/vlidort_lin_sup_ss_def.f90    \
   $(VLID_DEF_PATH)/vlidort_lin_sup_sleave_def.f90\
   $(VLID_DEF_PATH)/vlidort_lin_sup_def.f90       \
   $(VLID_DEF_PATH)/vlidort_lin_outputs_def.f90   \
   $(VLID_DEF_PATH)/vlidort_lin_io_defs.f90       \
   $(VLID_DEF_PATH)/vlidort_lin_work_def.f90      \
   $(VLID_MAIN_PATH2)/vlidort_l_inputs.f90        \
   $(VLID_MAIN_PATH2)/vlidort_l_writemodules.f90  \
   $(VLID_MAIN_PATH2)/vlidort_la_miscsetups.f90   \
   $(VLID_MAIN_PATH2)/vlidort_l_pack.f90          \
   $(VLID_MAIN_PATH2)/vlidort_l_unpack.f90        \
   $(VLID_MAIN_PATH2)/vlidort_l_thermalsup.f90    \
   $(VLID_MAIN_PATH2)/vlidort_lpc_solutions.f90   \
   $(VLID_MAIN_PATH2)/vlidort_lpc_bvproblem.f90   \
   $(VLID_MAIN_PATH2)/vlidort_lbbf_jacobians_vector.f90 \
   $(VLID_MAIN_PATH2)/vlidort_ls_wfsurface.f90    \
   $(VLID_MAIN_PATH2)/vlidort_ls_wfsleave.f90

LCS_SOURCES += \
   $(VLID_MAIN_PATH2)/vlidort_lc_miscsetups.f90       \
   $(VLID_MAIN_PATH2)/vlidort_vfo_lcs_interface.f90   \
   $(VLID_MAIN_PATH2)/vlidort_lc_pack.f90             \
   $(VLID_MAIN_PATH2)/vlidort_lc_unpack.f90           \
   $(VLID_MAIN_PATH2)/vlidort_lc_mediaprops.f90       \
   $(VLID_MAIN_PATH2)/vlidort_lc_solutions.f90        \
   $(VLID_MAIN_PATH2)/vlidort_lc_bvproblem.f90        \
   $(VLID_MAIN_PATH2)/vlidort_lc_PostProcessing.f90   \
   $(VLID_MAIN_PATH2)/vlidort_lc_wfatmos.f90          \
   $(VLID_MAIN_PATH2)/vlidort_lcs_converge.f90        \
   $(VLID_MAIN_PATH2)/vlidort_lcs_masters_V2p8p3.f90

LPS_SOURCES += \
   $(VLID_MAIN_PATH2)/vlidort_lp_miscsetups.f90      \
   $(VLID_MAIN_PATH2)/vlidort_vfo_lps_interface.f90  \
   $(VLID_MAIN_PATH2)/vlidort_lp_pack.f90            \
   $(VLID_MAIN_PATH2)/vlidort_lp_unpack.f90          \
   $(VLID_MAIN_PATH2)/vlidort_lp_mediaprops.f90      \
   $(VLID_MAIN_PATH2)/vlidort_lp_solutions.f90       \
   $(VLID_MAIN_PATH2)/vlidort_lp_bvproblem.f90       \
   $(VLID_MAIN_PATH2)/vlidort_lp_PostProcessing.f90  \
   $(VLID_MAIN_PATH2)/vlidort_lp_wfatmos.f90         \
   $(VLID_MAIN_PATH2)/vlidort_lps_converge.f90       \
   $(VLID_MAIN_PATH2)/vlidort_lps_masters_V2p8p3.f90

# (Include vector supplement source files)
include $(VSUP_PATH)/makefile.vsup

# (Include first-order source files)
include $(FO_MAIN_PATH)/makefile.fo_vlidort

SOURCES_MIE = 
SOURCES_MIE += \
   $(MIE_PATH)/RTSMie_Parameters.f90      \
   $(MIE_PATH)/RTSMie_Inputs_Def.f90      \
   $(MIE_PATH)/RTSMie_Outputs_Def.f90     \
   $(MIE_PATH)/RTSMie_Lin_Outputs_Def.f90 \
   $(MIE_PATH)/RTSMie_Distribution.f90    \
   $(MIE_PATH)/RTSMie_IO_Readwrite.f90    \
   $(MIE_PATH)/RTSMie_Master.f90          \
   $(MIE_PATH)/RTSMie_Master_PLUS.f90     \
   $(MIE_PATH)/RTSMie_Master_Bimodal.f90  \
   $(MIE_PATH)/RTSMie_Master_Bimodal_PLUS.f90

SOURCES_MIE_OMP = 
SOURCES_MIE_OMP += \
   $(MIE_PATH)/RTSMie_Parameters.f90      \
   $(MIE_PATH)/RTSMie_Inputs_Def.f90      \
   $(MIE_PATH)/RTSMie_Outputs_Def.f90     \
   $(MIE_PATH)/RTSMie_Lin_Outputs_Def.f90 \
   $(MIE_PATH)/RTSMie_Distribution.f90    \
   $(MIE_PATH)/RTSMie_IO_Readwrite.f90    \
   $(MIE_PATH)/RTSMie_Master.f90          \
   $(MIE_PATH)/RTSMie_Master_PLUS.f90     \
   $(MIE_PATH)/RTSMie_Master_Bimodal.f90  \
   $(MIE_PATH)/RTSMie_Master_Bimodal_PLUS.f90 \
   $(MIE_PATH)/RTSMie_Master_Bimodal_OMP.f90  \
   $(MIE_PATH)/RTSMie_Master_Bimodal_PLUS_OMP.f90

SOURCES_TONGA =
SOURCES_TONGA += \
   $(TONGA_SOURCES_PATH)/TONGA_Pars.f90                \
   $(TONGA_SOURCES_PATH)/TONGA_Types_HPTORA_V5.f90     \
   $(TONGA_SOURCES_PATH)/TONGA_Create_HPTORA_V5.f90    \
   $(TONGA_SOURCES_PATH)/Loading_General_Routines.f90  \
   $(TONGA_SOURCES_PATH)/TONGA_Aerosol_Loading.f90     \
   $(TONGA_SOURCES_PATH)/TONGA_SO2_Loading.f90         \
   $(TONGA_SOURCES_PATH)/TONGA_Mie_Initialize.f90      \
   $(TONGA_SOURCES_PATH)/TONGA_RT_Initialize_V2.f90    \
   $(TONGA_SOURCES_PATH)/TONGA_Aux.f90                 \
   $(TONGA_SOURCES_PATH)/TONGA_Forward_Model_V4.f90    \
   $(TONGA_SOURCES_PATH)/TONGA_Inverse_Model_V4.f90

SOURCES_TONGA_OMP =
SOURCES_TONGA_OMP += \
   $(TONGA_SOURCES_PATH)/TONGA_Pars.f90                \
   $(TONGA_SOURCES_PATH)/TONGA_Types_HPTORA_V5.f90     \
   $(TONGA_SOURCES_PATH)/TONGA_Create_HPTORA_V5.f90    \
   $(TONGA_SOURCES_PATH)/Loading_General_Routines.f90  \
   $(TONGA_SOURCES_PATH)/TONGA_Aerosol_Loading.f90     \
   $(TONGA_SOURCES_PATH)/TONGA_SO2_Loading.f90         \
   $(TONGA_SOURCES_PATH)/TONGA_Mie_Initialize.f90      \
   $(TONGA_SOURCES_PATH)/TONGA_RT_Initialize_V2.f90    \
   $(TONGA_SOURCES_PATH)/TONGA_Aux.f90                 \
   $(TONGA_SOURCES_PATH)/TONGA_Forward_Model_V4.f90    \
   $(TONGA_SOURCES_PATH)/TONGA_Inverse_Model_V4_OMP.f90

SOURCES_SOLAR_LPCS = $(FO_SOURCES_L_Vector) + \
   $(SOURCES) \
   $(L_SOURCES) \
   $(LPS_SOURCES) \
   $(LCS_SOURCES) \
   $(VLID_TEST_PATH)/2p8p3_solar_lpcs_tester.f90

SOURCES_RETRIEVAL_V4 = $(FO_SOURCES_L_Vector) + \
   $(SOURCES) \
   $(L_SOURCES) \
   $(LCS_SOURCES) \
   $(SOURCES_MIE) \
   $(SOURCES_TONGA) \
   $(TONGA_TEST_PATH)/Retrieval_MeasRatio_Mk2_with_SO2_V4.f90

SOURCES_RETRIEVAL_V4_OMP = $(FO_SOURCES_L_Vector) + \
   $(SOURCES) \
   $(L_SOURCES) \
   $(LCS_SOURCES) \
   $(SOURCES_MIE_OMP) \
   $(SOURCES_TONGA_OMP) \
   $(TONGA_TEST_PATH)/Retrieval_MeasRatio_Mk2_with_SO2_V4_OMP.f90

# Define pattern rules for creating object files:

# For vlidort main source files
$(OBJ_PATH)/%.o : $(VLID_DEF_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_MAIN_PATH1)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_MAIN_PATH2)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLID_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(MIE_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(TONGA_SOURCES_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(TONGA_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# Define object files

F90SOURCES_SOLAR_LPCS := $(notdir $(filter %.f90, $(SOURCES_SOLAR_LPCS)))
F90OBJECTS_SOLAR_LPCS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_SOLAR_LPCS)))

F90SOURCES_TONGA_RETRIEVAL_V4 := $(notdir $(filter %.f90, $(SOURCES_RETRIEVAL_V4)))
F90OBJECTS_TONGA_RETRIEVAL_V4 := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_TONGA_RETRIEVAL_V4)))

F90SOURCES_TONGA_RETRIEVAL_V4_OMP := $(notdir $(filter %.f90, $(SOURCES_RETRIEVAL_V4_OMP)))
F90OBJECTS_TONGA_RETRIEVAL_V4_OMP := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_TONGA_RETRIEVAL_V4_OMP)))

# Define desired targets

s2p8p3_solar_lpcs_tester.exe: $(F90OBJECTS_SOLAR_LPCS)
	$(FC) $^ -o $@


#2-band versions (serial & OpenMP)
Retrieval_MeasRatio_Mk2_with_SO2_V4.exe: $(F90OBJECTS_TONGA_RETRIEVAL_V4)
	$(FC) $^ -o $@

Retrieval_MeasRatio_Mk2_with_SO2_V4_OMP.exe: $(F90OBJECTS_TONGA_RETRIEVAL_V4_OMP)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@
#Retrieval_MeasRatio_Mk2_with_SO2_V4_OMP.exe: $(F90OBJECTS_TONGA_RETRIEVAL_V4_OMP)
#	$(FC)  $^ -o $@


.PHONY: list_solar_test_files
list_solar_test_files:

	@echo
	@echo 
	@echo "The native solar test source set looks like:"
	@echo
	@echo "$(SOURCES_SOLAR_LPCS)"

	@echo
	@echo 
	@echo "The solar test source set looks like:"
	@echo
	@echo "$(F90SOURCES_SOLAR_LPCS)"

	@echo
	@echo 
	@echo "The solar test object set looks like:"
	@echo
	@echo "$(F90OBJECTS_SOLAR_LPCS)"
	@echo
	@echo

# for Tonga 2-band serial code
.PHONY: list_tonga_ret_files_2
list_tonga_ret_files_2:
	@echo
	@echo 
	@echo "The tonga ret 2-band source set looks like:"
	@echo
	@echo "$(F90SOURCES_TONGA_RETRIEVAL_V4)"

	@echo
	@echo 
	@echo "The tonga ret 2-band object set looks like:"
	@echo
	@echo "$(F90OBJECTS_TONGA_RETRIEVAL_V4)"
	@echo
	@echo

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

