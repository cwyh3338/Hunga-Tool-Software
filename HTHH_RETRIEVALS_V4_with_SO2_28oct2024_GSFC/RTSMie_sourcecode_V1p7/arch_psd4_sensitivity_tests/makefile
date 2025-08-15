
#  Mie Path

MIE_PATH  = ../RTSMie_sourcecode_V1p7/

#  mod/obj

MOD_PATH = mod
OBJ_PATH = obj
MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

# Define default shell make will use

SHELL = /bin/bash

# Define FORTRAN90 compiler to use (can be defined here, but usually defined on makefile command line)

#  Intel
#FC = ifort

#  gfortran
FC = gfortran

#  g95
#FC = g95

#  NAG
#FC = f95

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


SOURCES_MieDists = 
SOURCES_MieDists += \
   $(MIE_PATH)/RTSMie_Parameters.f90 \
   $(MIE_PATH)/RTSMie_Inputs_Def.f90 \
   $(MIE_PATH)/RTSMie_Outputs_Def.f90 \
   $(MIE_PATH)/RTSMie_Distribution.f90 \
   $(MIE_PATH)/RTSMie_Master.f90       \
   $(MIE_PATH)/MieDists_2pars_Sensitivity.f90

# Define pattern rules for creating object files:

$(OBJ_PATH)/%.o : $(MIE_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# Define object files

F90SOURCES_MieDists := $(notdir $(filter %.f90, $(SOURCES_MieDists)))
F90OBJECTS_MieDists  := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_MieDists)))

# Define desired target

All: MieDists_2pars_Sensitivity.exe
MieDists_2pars_Sensitivity.exe: $(F90OBJECTS_MieDists)
	$(FC) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

