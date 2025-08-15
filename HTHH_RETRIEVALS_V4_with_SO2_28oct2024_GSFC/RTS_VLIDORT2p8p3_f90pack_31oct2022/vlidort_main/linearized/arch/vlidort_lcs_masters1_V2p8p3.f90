
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2021.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VLIDORT_LCS_MASTER (master)                      #
! #                                                             #
! ###############################################################

!  1/31/21. Version 2.8.3. SEVERAL INNOVATIONS and CHANGES

!  1. Major Change for the Green's function implementation
!    -- Controlled by flag DO_CLASSICAL_SOLUTION. Only works for NSTOKES = 1 or 3
!    -- In vlidort_solutions, add GBEAM and GUSER subroutines using Green's function methods
!    -- Whole/part Sourceterm routines completely rewritten (Green's function included)
!    -- QuadIntens partial-layer subroutines completely rewritten (Green's function included)
!    -- Need additional Taylor routines to be used, for Green's function post-processing

!  2. Other changes
!    -- Converge routines have been moved to their own module (vlidort_converge.f90)
!    -- MSST option is now included, generates output for sphericity correction
!    -- DO_DOUBLET_GEOMETRY post-processing operation is in force
!    -- LATTICE/DOUBLET OFFSETS are created once and for all here in the main routine
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier component
!    -- Output type structures are filled directly in the converge routines
!    -- Use of do_debug_input dump flag, now controlled from outside the main routine
!    -- Use of the input quantity "TOLERANCE" to ASYMTX, help to avoid eigenproblem non-convergence
!    -- (2/16/21). Use of the FOURIER0_NSTOKES2 option 

      MODULE vlidort_lcs_masters1_m

      PUBLIC  :: VLIDORT_LCS_MASTER

      CONTAINS

      SUBROUTINE VLIDORT_LCS_MASTER ( do_debug_input, &
        VLIDORT_FixIn, &
        VLIDORT_ModIn, &
        VLIDORT_Sup, &
        VLIDORT_Out, &
        VLIDORT_LinFixIn, &
        VLIDORT_LinModIn, &
        VLIDORT_LinSup, &
        VLIDORT_LinOut )

!  parameter file

      USE VLIDORT_PARS_m

!  I/O Type structures

      USE VLIDORT_Inputs_def_m
      USE VLIDORT_Sup_InOut_def_m
      USE VLIDORT_Outputs_def_m

      USE VLIDORT_LinInputs_def_m
      USE VLIDORT_LinSup_InOut_def_m
      USE VLIDORT_LinOutputs_def_m

!  internal Type structure

      USE VLIDORT_Work_def_m
      USE VLIDORT_LinWork_def_m

!  Standard Modules
!  ----------------

!  Modules. Version 2.8, only use selections
!    -- 1/31/21. Version 2.8.3. EMULT_MASTER, EMULT_MASTER_OBSGEO now in MiscSetups

      USE VLIDORT_INPUTS_m      , only : VLIDORT_CHECK_INPUT_DIMS, VLIDORT_CHECK_INPUT, VLIDORT_DERIVE_INPUT
      USE VLIDORT_GEOMETRY_m    , only : VLIDORT_CHAPMAN
      USE VLIDORT_MISCSETUPS_m  , only : VLIDORT_MISCSETUPS, EMULT_MASTER, EMULT_MASTER_OBSGEO
      USE VLIDORT_THERMALSUP_m  , only : THERMAL_SETUP

!mick debug - VLIDORT_LCS_MASTER and VLIDORT_LCS_FOURIER temporarily put in two
!             separate modules; thus, we must have a USE statement
      USE VLIDORT_LCS_FOURIER1_m , only : VLIDORT_LCS_FOURIER

!  1/31/21. Version 2.8.3. Use VLIDORT_CONVERGE_m for the three convergenece routines

      USE VLIDORT_CONVERGE_m    , only : VLIDORT_CONVERGE, VLIDORT_CONVERGE_OBSGEO, VLIDORT_CONVERGE_DOUBLET

      USE VLIDORT_PACK_m
      USE VLIDORT_WRITEMODULES_m

!  7/7/16, RT Solutions. Version 2.8, only require the FO interface. CORRECTIONS removed
!      USE VLIDORT_CORRECTIONS         !  removed.

      USE VLIDORT_VFO_LCS_INTERFACE_m

!  4/9/19. TRANSFLUX superceded, replaced by internal "Adjusted_Backsub" routine
!  VLIDORT 2.8, 9/25/15. RT Solutions, new "transflux" module (Mark1, Mark2)
!  VLIDORT 2.8, 2/3/16 . RT Solutions, new "transflux" module (Mark 3)
!    USE Vlidort_transflux_Master_m
!      USE vlidort_transflux_m

!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination, and flux
!  --- THIS IS SEPARATE from the MEDIA_PROPERTIES installation.

!  Linearized Modules
!  ------------------

!  general linearizations (applicable to LCS and LPS)

      USE VLIDORT_L_INPUTS_m,  Only : VLIDORT_L_CHECK_INPUT_DIMS, VLIDORT_L_CHECK_INPUT
      USE VLIDORT_L_PACK_m
      USE VLIDORT_L_WRITEMODULES_m

!  thermal linearizations.

      USE VLIDORT_L_THERMALSUP_m, Only : THERMAL_SETUP_PLUS

!  Specific LP linearization modules

      USE VLIDORT_LC_PACK_m
      USE VLIDORT_LC_MISCSETUPS_m, Only : VLIDORT_LAC_MISCSETUPS, LC_EMULT_MASTER, LC_EMULT_MASTER_OBSGEO

!  1/31/21. Version 2.8.3. Use VLIDORT_LCS_CONVERGE_m for the three convergence routines

      USE VLIDORT_LCS_CONVERGE_m

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      IMPLICIT NONE

!  1/31/21. Version 2.8.3. input argument for the debug dumping factility

      LOGICAL, INTENT (IN) :: DO_DEBUG_INPUT

!  VLIDORT input structures

      TYPE(VLIDORT_Fixed_Inputs), INTENT (IN)       :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), INTENT (INOUT) :: VLIDORT_ModIn

!  VLIDORT supplements structure

      TYPE(VLIDORT_Sup_InOut), INTENT (INOUT)       :: VLIDORT_Sup

!  VLIDORT output structure

      TYPE(VLIDORT_Outputs), INTENT (OUT)           :: VLIDORT_Out

!  VLIDORT linearized input structures

      TYPE(VLIDORT_Fixed_LinInputs), INTENT (IN)       :: VLIDORT_LinFixIn
      TYPE(VLIDORT_Modified_LinInputs), INTENT (INOUT) :: VLIDORT_LinModIn

!  VLIDORT linearized supplements structure

      TYPE(VLIDORT_LinSup_InOut), INTENT (INOUT)       :: VLIDORT_LinSup

!  VLIDORT linearized output structure

      TYPE(VLIDORT_LinOutputs), INTENT (OUT)           :: VLIDORT_LinOut

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

      LOGICAL ::            DO_FULLRAD_MODE

      LOGICAL ::            DO_THERMAL_EMISSION
      LOGICAL ::            DO_SURFACE_EMISSION
      LOGICAL ::            DO_PLANE_PARALLEL

      LOGICAL ::            DO_UPWELLING
      LOGICAL ::            DO_DNWELLING

!      LOGICAL ::            DO_QUAD_OUTPUT       ! removed, 7/7/16
      LOGICAL ::            DO_TOA_CONTRIBS

      LOGICAL ::            DO_LAMBERTIAN_SURFACE

      LOGICAL ::            DO_SPECIALIST_OPTION_1
      LOGICAL ::            DO_SPECIALIST_OPTION_2
      LOGICAL ::            DO_SPECIALIST_OPTION_3

!  New 17 May 2012, surface leaving flags
      LOGICAL ::            DO_SURFACE_LEAVING
      LOGICAL ::            DO_SL_ISOTROPIC
      LOGICAL ::            DO_WATER_LEAVING       ! New 10/28/15
      LOGICAL ::            DO_FLUORESCENCE        ! New 10/28/15
      LOGICAL ::            DO_TF_ITERATION        ! New 7/7/16

!  Water-leaving output flag. 3/18/19 for Version 2.8.1
      
      LOGICAL ::            DO_WLADJUSTED_OUTPUT 

!   4/26/19 Added control for the media problem. Version 2.8.1
!     Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!       1 = Isotropic illumination from TOA, 2 = Isotropic illumination from BOA
      
      LOGICAL ::            DO_ALBTRN_MEDIA(2)

!   4/28/19 Added control for the planetary problem. Version 2.8.1
      !     Linked to the Media-problem, requires some flux and tranmsittance output (BOA unit illumination)

      LOGICAL ::            DO_PLANETARY_PROBLEM

!  TOA and BOA Illumination flags. 3/23/19 for Version 2.8.1
!   Airglow and Nighttime-viewing scenarios
      
      LOGICAL ::            DO_TOAFLUX
      LOGICAL ::            DO_BOAFLUX

!  1/31/21. Version 2.8.3. (RTS 2/16/21). Flag for using an NSTOKES = 2 calculation for Fourier 0
!    -- This must be set by hand, it is not a configuration file read

      LOGICAL ::            DO_FOURIER0_NSTOKES2

!  Order of Taylor series (including terms up to EPS^n).
!        Introduced 2/19/14 for Version 2.7
      
      INTEGER ::            TAYLOR_ORDER

!  Basic control numbers

      INTEGER ::            NSTOKES
      INTEGER ::            NSTREAMS
      INTEGER ::            NLAYERS
      INTEGER ::            NFINELAYERS
      INTEGER ::            N_THERMAL_COEFFS

!  Accuracy

      DOUBLE PRECISION ::   VLIDORT_ACCURACY

!  Special inputs (RT Solutions use only)

      INTEGER ::            NLAYERS_NOMS
      INTEGER ::            NLAYERS_CUTOFF

!  Water-leaving: Control for iterative calculation of transmittanaces
!    Variables added for Version 2.8, 7/6/16

      INTEGER          ::   TF_MAXITER
      DOUBLE PRECISION ::   TF_CRITERION

!  TOA Illumination, Flux value. 3/23/19 for Version 2.8.1
!    Must be solar-flux normalized. Designed for Airglow Studies
      
      DOUBLE PRECISION ::   TOAFLUX
      
!  BOA Illumination, Flux value. 3/23/19 for Version 2.8.1
!    Must be solar-flux normalized. Designed for Nighttime Studies.
      
      DOUBLE PRECISION ::   BOAFLUX

!  1/31/21. Version 2.8.3. This is the input tolerance variable

      DOUBLE PRECISION ::   ASYMTX_TOLERANCE

!  Flux factor

      DOUBLE PRECISION ::   FLUX_FACTOR

!  User levels

      INTEGER ::            N_USER_LEVELS

!  PTH and fine grid

      DOUBLE PRECISION ::   HEIGHT_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION ::   PRESSURE_GRID ( 0:MAXLAYERS )
      DOUBLE PRECISION ::   TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER ::            FINEGRID ( MAXLAYERS )
      DOUBLE PRECISION ::   RFINDEX_PARAMETER

!  Optical properties

      DOUBLE PRECISION ::   DELTAU_VERT_INPUT ( MAXLAYERS )
      DOUBLE PRECISION ::   GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  F-matrix optical properties (New for Version 2.8, 7/7/16)
      DOUBLE PRECISION ::   FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 )
      DOUBLE PRECISION ::   FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 )

      DOUBLE PRECISION ::   ALBEDO
      DOUBLE PRECISION ::   THERMAL_BB_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION ::   SURFBB

!  These "LTE" inputs were in Version 2.6.
!      DOUBLE PRECISION ::   LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
!      DOUBLE PRECISION ::   LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )

      DOUBLE PRECISION ::   ATMOS_WAVELENGTH

      LOGICAL ::            DO_DEBUG_WRITE
      LOGICAL ::            DO_WRITE_INPUT
      LOGICAL ::            DO_WRITE_SCENARIO
      LOGICAL ::            DO_WRITE_FOURIER
      LOGICAL ::            DO_WRITE_RESULTS

      CHARACTER (LEN=60) :: INPUT_WRITE_FILENAME
      CHARACTER (LEN=60) :: SCENARIO_WRITE_FILENAME
      CHARACTER (LEN=60) :: FOURIER_WRITE_FILENAME
      CHARACTER (LEN=60) :: RESULTS_WRITE_FILENAME

!  --------------------------
!  Standard Inputs - Modified
!  --------------------------

!  FOCorr and SSCorr Booleans. Completely revised for Version 2.8, 3/1/17.

      LOGICAL ::          DO_FOCORR                ! New 02 Jul 2013, used to be DO_FO_CALC
      LOGICAL ::          DO_FOCORR_EXTERNAL       ! renamed 2.8
      LOGICAL ::          DO_FOCORR_ALONE          ! Used to be DO_SSFULL - internal variable now (Version 2.8)
      LOGICAL ::          DO_FOCORR_NADIR
      LOGICAL ::          DO_FOCORR_OUTGOING

      LOGICAL ::          DO_SSCORR_USEFMAT        ! New, 2.8

!  1/31/21. Version 2.8.3. Variable has been removed
!      LOGICAL ::          DO_SSCORR_TRUNCATION

!  Additional Control for Externalized input (SLEAVE). Introduced 3/18/19 for Version 2.8.1
      
      LOGICAL ::          DO_EXTERNAL_WLEAVE

!  Solar

      LOGICAL ::          DO_SOLAR_SOURCES
      LOGICAL ::          DO_REFRACTIVE_GEOMETRY
      LOGICAL ::          DO_CHAPMAN_FUNCTION

!  1/31/21. Version 2.8.3. Add DO_MSSTS flag for source term output

      LOGICAL ::          DO_MSSTS

!  1/31/21. Version 2.8.3. Add DO_DOUBLET_GEOMETRY flag

      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  1/31/21. Version 2.8.3. new variable for RT Solution method
!            TRUE = Classical solution of PI, FALSE = Green's function solution

      LOGICAL ::          DO_CLASSICAL_SOLUTION

!  Performance flags

      LOGICAL ::          DO_DOUBLE_CONVTEST
      LOGICAL ::          DO_RAYLEIGH_ONLY
      LOGICAL ::          DO_DELTAM_SCALING
      LOGICAL ::          DO_SOLUTION_SAVING
      LOGICAL ::          DO_BVP_TELESCOPING

!  Other RT model flags

      LOGICAL ::          DO_USER_VZANGLES
      LOGICAL ::          DO_ADDITIONAL_MVOUT
      LOGICAL ::          DO_MVOUT_ONLY
      LOGICAL ::          DO_THERMAL_TRANSONLY
      LOGICAL ::          DO_OBSERVATION_GEOMETRY

!  Integers

      INTEGER ::          NGREEK_MOMENTS_INPUT

!  Geometry integers
!mick mod 1/5/2021 - added N_USER_DOUBLETS
!                  - modified floating point angle & level declarations
!                    as in "vlidort_masters"

      INTEGER ::          N_SZANGLES
      INTEGER ::          N_USER_RELAZMS
      INTEGER ::          N_USER_VZANGLES
      INTEGER ::          N_USER_OBSGEOMS
      INTEGER ::          N_USER_DOUBLETS

!  BOA solar zenith angles (degrees)

      DOUBLE PRECISION :: SZANGLES      ( MAX_SZANGLES )

!  User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      DOUBLE PRECISION :: USER_LEVELS   ( MAX_USER_LEVELS )

!  User-defined viewing zenith and azimuth angles input (degrees) 

      DOUBLE PRECISION :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      DOUBLE PRECISION :: USER_VZANGLES ( MAX_USER_VZANGLES )

!  User-defined Observation Geometry angle input
!   New variable, 25 OCtober 2012, for Observational Geometry input
!     -- 1/31/21. Version 2.8.3. (2/16/21). Doublet-Geometry angle input

      DOUBLE PRECISION :: USER_DOUBLETS ( MAX_USER_VZANGLES, 2 )
      DOUBLE PRECISION :: USER_OBSGEOMS ( MAX_USER_OBSGEOMS, 3 )

!  Spec height and earth radius

      DOUBLE PRECISION :: GEOMETRY_SPECHEIGHT
      DOUBLE PRECISION :: EARTH_RADIUS

!  Single scattering albedo

      DOUBLE PRECISION :: OMEGA_TOTAL_INPUT ( MAXLAYERS )

!  -----------------------
!  Standard Supplement I/O
!  -----------------------

!  BRDF Inputs
!  -----------

!  1/31/21. Version 2.8.3. SLEAVE Fourier terms defined locally (Drop MAXMOMENTS dimension)

      DOUBLE PRECISION ::   EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION ::   BRDF_F_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   BRDF_F   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_BRDF_F_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_BRDF_F   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION ::   EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION ::   USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  SS and DB I/O
!  -------------

!  SS/DB Intensity Results at all angles and optical depths
!   1/31/21. Version 2.8.3. Not used except for debug input and output write.
!   1/31/21. Version 2.8.3. (2/16/21). Not used at all now,
!      DOUBLE PRECISION ::  STOKES_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
!      DOUBLE PRECISION ::  STOKES_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Contribution function (TOA Upwelling only)
!  SS component of Diffuse Field
!      DOUBLE PRECISION ::  SS_CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  New Surface-Leaving Inputs, 17 May 12
!  -------------------------------------

      DOUBLE PRECISION ::   SLTERM_ISOTROPIC  ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION ::   SLTERM_USERANGLES ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

!  1/31/21. Version 2.8.3. SLEAVE Fourier terms defined locally.

      DOUBLE PRECISION ::   SLTERM_F_0      ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   USER_SLTERM_F_0 ( MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  ----------------
!  Standard Outputs
!  ----------------

!  Fourier-summed values
!   1/31/21. Version 2.8.3. (2/16/21). Not used at all now,
!      DOUBLE PRECISION :: STOKES ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  Results for mean-value output
!  -----------------------------

!mick mod 9/19/2017 - renamed fluxes for separating diffuse and direct

!  Complete Actinic and Regular Fluxes (including Direct terms)

      DOUBLE PRECISION :: MEANST_DIFFUSE ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FLUX_DIFFUSE   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )

!  Direct Fluxes only

      DOUBLE PRECISION :: DNMEANST_DIRECT ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION :: DNFLUX_DIRECT   ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  MEDIA-PROPERTY Output
!  ---------------------
      
!  4/26-28/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Output for User-angles and fluxes, Transmittances for solar beam

      DOUBLE PRECISION :: ALBMED_USER ( MAXSTOKES, MAX_USER_STREAMS ), ALBMED_FLUXES(MAXSTOKES,2)    !  TOA illumination
      DOUBLE PRECISION :: TRNMED_USER ( MAXSTOKES, MAX_USER_STREAMS ), TRNMED_FLUXES(MAXSTOKES,2)    !  BOA illumination
      DOUBLE PRECISION :: TRANSBEAM   ( MAXSTOKES, MAXBEAMS )                                        !  Planetary problem

!  TF Output: Complete/Direct Actinic and Regular Fluxes
!     RT Solutions, 9/25/15. removed 2/3/16 (superceded by Mark 3 output)
!      DOUBLE PRECISION :: MEANST_DIFFUSE_TF ( MAX_SZANGLES, MAXSTOKES )
!      DOUBLE PRECISION :: FLUX_DIFFUSE_TF ( MAX_SZANGLES, MAXSTOKES )
!      DOUBLE PRECISION :: DNMEANST_DIRECT_TF ( MAX_SZANGLES, MAXSTOKES )
!      DOUBLE PRECISION :: DNFLUX_DIRECT_TF ( MAX_SZANGLES, MAXSTOKES )

!  1/31/21. Version 2.8.3. Installed final version of DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier input, Converged_I output(upwelling case)
!    ==> Final MSST values are filled out directly in the Converge_Obsgeo routine 
!    ==> MSST situation expanded to include Downwelling case (as an alternative)

      DOUBLE PRECISION :: LAYER_MSSTS_F  ( MAXBEAMS, MAXSTOKES, MAXLAYERS  )
      DOUBLE PRECISION :: SURF_MSSTS_F   ( MAXBEAMS, MAXSTOKES  )

!  1/31/21. Version 2.8.3. Installed MSST LC and LS linearizations
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)

      DOUBLE PRECISION :: LC_SURF_MSSTS_F  ( MAXBEAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAXSTOKES, MAX_SURFACEWFS )
      DOUBLE PRECISION :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_SURFACEWFS )

!  Contribution functions (TOA Upwelling only)
!  -------------------------------------------

!  Fourier component of Diffuse Field
!      DOUBLE PRECISION ::  MS_CONTRIBS_F ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
!  Fourier-summed values
!      DOUBLE PRECISION ::  CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  Ancillary Output
!  ----------------

!  Fourier numbers used

      INTEGER ::          FOURIER_SAVED ( MAX_SZANGLES )

!  Number of geometries computed

      INTEGER ::          N_GEOMETRIES

!  Offsets for geometry indexing
!   -- 1/31/21. Version 2.8.3. Add Doublet geometry offsets (SZD_OFFSETS) 

      INTEGER         :: SZA_OFFSETS(MAX_SZANGLES)
      INTEGER         :: VZA_OFFSETS(MAX_SZANGLES,MAX_USER_VZANGLES)
      INTEGER         :: SZD_OFFSETS(MAX_SZANGLES)

!  -------------------------
!  Linearized Inputs - Fixed
!  -------------------------

!  Linearization control

      LOGICAL ::            DO_SIMULATION_ONLY
      LOGICAL ::            LAYER_VARY_FLAG ( MAXLAYERS )
      INTEGER ::            LAYER_VARY_NUMBER ( MAXLAYERS )

!  Number of Jacobians. Added TotalSurface 9/18/16, N_TOTALSURFACE = N_SURFACE + N_SLEAVE

      INTEGER ::            N_TOTALCOLUMN_WFS
      INTEGER ::            N_TOTALPROFILE_WFS
      INTEGER ::            N_TOTALSURFACE_WFS
      INTEGER ::            N_SURFACE_WFS
      INTEGER ::            N_SLEAVE_WFS

!  names

      CHARACTER (LEN=31) :: COLUMNWF_NAMES ( MAX_ATMOSWFS )
      CHARACTER (LEN=31) :: PROFILEWF_NAMES ( MAX_ATMOSWFS )

!  Linearized optical properties

      DOUBLE PRECISION ::   L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION ::   L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION ::   L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Linearized Fmatrix inputs for FO calculations
!mick fix 9/19/2017 - swapped layer & geo indices

      !DOUBLE PRECISION ::   L_FMATRIX_UP ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, 6 )
      !DOUBLE PRECISION ::   L_FMATRIX_DN ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS, 6 )
      DOUBLE PRECISION ::   L_FMATRIX_UP ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES,  6 )
      DOUBLE PRECISION ::   L_FMATRIX_DN ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES,  6 )

!  ----------------------------
!  Linearized Inputs - Variable
!  ----------------------------

!  Linearization flags
!    -- 5/24/21. Version 2.8.3. Add DO_INCLUDE_SLEAVEWFS

      LOGICAL ::         DO_COLUMN_LINEARIZATION
      LOGICAL ::         DO_PROFILE_LINEARIZATION
      LOGICAL ::         DO_ATMOS_LINEARIZATION
      LOGICAL ::         DO_SURFACE_LINEARIZATION
      LOGICAL ::         DO_LINEARIZATION
      LOGICAL ::         DO_SLEAVE_WFS, DO_INCLUDE_SLEAVEWFS

!  Control for  Blackbody Jacobians, New 28 March 2014
!   Replaces the two commented out flags
!      LOGICAL ::            DO_LTE_LINEARIZATION
!      LOGICAL ::            DO_SURFBB_LINEARIZATION

      LOGICAL :: DO_ATMOS_LBBF, DO_SURFACE_LBBF

!  -------------------------
!  Linearized Supplement I/O
!  -------------------------

!  BRDF Inputs
!  -----------

      DOUBLE PRECISION ::   LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Fourier-component BRDFs are defined locally
!     -- Drop the MAXMOMENTS dimension

      DOUBLE PRECISION ::   LS_BRDF_F_0 ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   LS_BRDF_F   ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

      DOUBLE PRECISION ::   LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION ::   LS_USER_BRDF_F   ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )

      DOUBLE PRECISION ::   LS_EMISSIVITY      ( MAX_SURFACEWFS, MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION ::   LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAXSTOKES, MAX_USER_STREAMS )

!  SS and DB I/O
!  -------------

!mick note 1/5/2021 - Version 2.8.3. Not used except for debug input and output write.
!   1/31/21. Version 2.8.3. (2/16/21). Not used at all now,
!    -- These DO NOT NEED to be left in, as they NOT NEEDED to be dummy-initialized for output.

!      DOUBLE PRECISION :: COLUMNWF_SS  ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
!      DOUBLE PRECISION :: COLUMNWF_DB  ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
!      DOUBLE PRECISION :: PROFILEWF_SS ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
!      DOUBLE PRECISION :: PROFILEWF_DB ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
!      DOUBLE PRECISION :: SURFACEWF_DB ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  New Surface-Leaving Inputs, 22 Aug 12
!  -------------------------------------

!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  1/31/21. Version 2.8.3. Fourier components are defined locally. Drop MAXMOMENTS Dimension

      DOUBLE PRECISION ::   LSSL_SLTERM_ISOTROPIC ( MAX_SLEAVEWFS, MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION ::   LSSL_SLTERM_USERANGLES( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

      DOUBLE PRECISION ::   LSSL_SLTERM_F_0      ( MAX_SLEAVEWFS, MAXSTOKES, MAXSTREAMS,      MAXBEAMS )
      DOUBLE PRECISION ::   LSSL_USER_SLTERM_F_0 ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  ------------------
!  Linearized Outputs
!  ------------------

!mick mod 9/19/2017 - renamed linearized fluxes for separating diffuse and direct
!   1/31/21. Version 2.8.3. (2/16/21). Fourier-summed values Not used at all now,
!      DOUBLE PRECISION :: COLUMNWF  ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
!      DOUBLE PRECISION :: SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  Mean intensity (actinic flux) and regular flux weighting functions. TOTAL

      DOUBLE PRECISION :: MEANST_DIFFUSE_COLWF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: FLUX_DIFFUSE_COLWF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Mean intensity (actinic flux) and regular flux weighting functions.
!  Direct beam only.

      DOUBLE PRECISION :: DNMEANST_DIRECT_COLWF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION :: DNFLUX_DIRECT_COLWF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  Surface weighting functions

      DOUBLE PRECISION :: MEANST_DIFFUSE_SURFWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

      DOUBLE PRECISION :: FLUX_DIFFUSE_SURFWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  4/26-29/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Linearized Output developed for Column Jacobians.

      DOUBLE PRECISION :: LC_ALBMED_USER  ( MAXSTOKES, MAX_USER_STREAMS, MAX_ATMOSWFS ) !  TOA illumination
      DOUBLE PRECISION :: LC_TRNMED_USER  ( MAXSTOKES, MAX_USER_STREAMS, MAX_ATMOSWFS ) !  BOA illumination
      
      DOUBLE PRECISION :: LC_ALBMED_FLUXES ( MAXSTOKES, 2, MAX_ATMOSWFS )    !  TOA illumination
      DOUBLE PRECISION :: LC_TRNMED_FLUXES ( MAXSTOKES, 2, MAX_ATMOSWFS )    !  BOA illumination
      
      DOUBLE PRECISION :: LC_TRANSBEAM ( MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS ) !  Planetary problem

!  BLACKBODY Jacobians, New 18 March 2014. Version 2.7
!  ===================================================

!  Enabled in Version 2.8

!  Postprocessed Jacobians.

      DOUBLE PRECISION :: ABBWFS_JACOBIANS &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS)
      DOUBLE PRECISION :: SBBWFS_JACOBIANS &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXSTOKES, MAX_DIRECTIONS)

!  Flux Jacobians.

      DOUBLE PRECISION :: ABBWFS_FLUXES &
          ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: SBBWFS_FLUXES &
          ( MAX_USER_LEVELS, 2, MAXSTOKES, MAX_DIRECTIONS)

!  RT Solutions Inc. RJD Spurr.  11 September 2009.
!    LTE linearization: Introduction of T-Jacobians for BB functions
!     Only works with pure thermal emission (no scattering)
!      Introduced for the GEOCAPE study. Version 2.6 Superseded
!      DOUBLE PRECISION :: LTE_ATMOSWF &
!          ( 0:MAXLAYERS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_DIRECTIONS )

!  Error handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NCHECKMESSAGES
      CHARACTER (LEN=120) :: CHECKMESSAGES (0:MAX_MESSAGES)
      CHARACTER (LEN=120) :: ACTIONS (0:MAX_MESSAGES)
      INTEGER ::             STATUS_CALCULATION
      CHARACTER (LEN=120) :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

! ######################################################################

!                       Local Arguments
!                       +++++++++++++++

! ######################################################################

!  Local bookkeeping
!  ----------------

!               Intent(In) To the Fourier  routine
!               Intent(In) To the Converge routine

!  Mode flags (Robfix, 13 January 2012. Add MSMODE_THERMAL flag)

      LOGICAL         :: DO_MSMODE_VLIDORT
      LOGICAL         :: DO_MSMODE_THERMAL

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER         :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER         :: NSTREAMS_2
      INTEGER         :: NTOTAL
      INTEGER         :: N_SUBDIAG, N_SUPDIAG

!  Quadrature weights and abscissae, and product

      DOUBLE PRECISION :: QUAD_STREAMS (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_WEIGHTS (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_STRMWTS (MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_HALFWTS ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_SINES   ( MAXSTREAMS )
      DOUBLE PRECISION :: QUAD_ANGLES  ( MAXSTREAMS )

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      !DOUBLE PRECISION :: USER_ANGLES   (MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_STREAMS  (MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES    (MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SECANTS  (MAX_USER_STREAMS)

!  Output optical depth masks and indices

      LOGICAL          :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER          :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER          :: LEVELMASK_UP  (MAX_USER_LEVELS)
      INTEGER          :: LEVELMASK_DN  (MAX_USER_LEVELS)

!  Off-grid optical depths (values, masks, indices)

      INTEGER          :: N_PARTLAYERS
      INTEGER          :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)
      DOUBLE PRECISION :: PARTLAYERS_VALUES       (MAX_PARTLAYERS)

!  Layer masks for doing integrated source terms

      LOGICAL         :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL         :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Indexing numbers

      INTEGER         :: N_VIEWING

!  1/31/21. Version 2.8.3. Post-processing masks.

      INTEGER         :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Local input solar zenith angles Cosines
!  ( Only required for refractive geometry attenuation of the solar beam)

      DOUBLE PRECISION :: SUN_SZA_COSINES(MAXLAYERS,MAXBEAMS)

!  Local solar zenith angles Cosines (regular case)

      DOUBLE PRECISION :: COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION :: SIN_SZANGLES ( MAX_SZANGLES )

!  Solar beam flags (always internal)

      LOGICAL         :: DO_MULTIBEAM (MAXBEAMS,0:MAXFOURIER)

!  Number of directions (1 or 2) and directional array

      INTEGER         :: N_DIRECTIONS
      INTEGER         :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Number of convergence tests

      INTEGER         :: N_CONVTESTS
      INTEGER         :: N_CONV_STREAMS

!  Arrays for setups and Corrections
!  ---------------------------------

!               Intent(In) To the Fourier routine

!  Local flags for the solution saving option

      INTEGER          :: LAYER_MAXMOMENTS (MAXLAYERS)

!  Initial transmittances * (secants)

      DOUBLE PRECISION :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Saved arrays for truncation factor and Delta-M scaling

      DOUBLE PRECISION :: TRUNC_FACTOR(MAXLAYERS)
      DOUBLE PRECISION :: FAC1(MAXLAYERS)

!  Derived Solar-beam Transmittance at all levels
!mick fix 9/19/2017 - added these four SOLARTRANS variables to facilitate correction of direct flux
!                     (standard & linearized)

      DOUBLE PRECISION :: LEVELS_SOLARTRANS   ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: SOLARBEAM_BOATRANS ( MAXBEAMS )

!  Derived Solar-beam Transmittance at BOA (Diagnostic output). 
!  Rob Fix 10/24/14, 11/17/14.
!    SOLARBEAM_BOATRANS is computed with Unscaled optical depths       ==> derivative  LC_SOLARBEAM_BOATRANS
!    TRANS_SOLAR_BEAM with scaled ODs. [ Same if  no Deltam-scaling]   ==> derivative  LC_SOLARBEAM_ATRANS

      DOUBLE PRECISION :: LC_SOLARBEAM_BOATRANS  ( MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_SOLARBEAM_ATRANS    ( MAXBEAMS, MAX_ATMOSWFS )

!  Derived Slant optical thickness inputs
!mick fix 9/19/2017 - added DELTAU_SLANT_UNSCALED & PARTAU_SLANT_UNSCALED
!                     to facilitate correction of direct flux

      DOUBLE PRECISION :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: DELTAU_SLANT_UNSCALED ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: PARTAU_SLANT_UNSCALED ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Scaled SSAs and phase function moments

      DOUBLE PRECISION :: OMEGA_TOTAL    ( MAXLAYERS )
      DOUBLE PRECISION :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  Linearized Input optical depths after delta-M scaling

      DOUBLE PRECISION :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_GREEKMAT_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

!  Derived Slant optical thickness inputs

      DOUBLE PRECISION :: L_DELTAU_SLANT (MAX_ATMOSWFS,MAXLAYERS,MAXLAYERS,MAXBEAMS)

!  Linearized truncation factor

      DOUBLE PRECISION :: L_TRUNC_FACTOR(MAX_ATMOSWFS,MAXLAYERS)

!  L'Hopital's rule logical variables

      LOGICAL          :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Coefficient functions for user-defined angles

      DOUBLE PRECISION :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      DOUBLE PRECISION :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Fourier component output
!  ------------------------

!               Intent(Out) from the Fourier routine
!               Intent(in)  to   the Converge routine

!  Fourier comonents User-defined solutions

      DOUBLE PRECISION :: STOKES_F &
          ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Fourier-component column weighting functions at user angles

      DOUBLE PRECISION :: COLUMNWF_F &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Fourier-component surface weighting functions at user angles

      DOUBLE PRECISION :: SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, &
            MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  Help arrays from the SS/DB correction routines
!  ==============================================

!  THIS SECTION HAS BEEN REMOVED for VERSION 2.8

!  Saved Legendre polynomials
!      !? DOUBLE PRECISION :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
!      !? DOUBLE PRECISION :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
!  Saved TMS (Nakajima-Tanaka) factor
!      DOUBLE PRECISION :: TMS ( MAXLAYERS )
!  Local truncation factors for additional DELTAM scaling
!      DOUBLE PRECISION :: SSFDEL ( MAXLAYERS )
!  Exact Phase function calculations
      !DOUBLE PRECISION :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
      !DOUBLE PRECISION :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)
!      DOUBLE PRECISION :: ZMAT_UP ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!      DOUBLE PRECISION :: ZMAT_DN ( MAX_GEOMETRIES, MAXLAYERS, MAXSTOKES, MAXSTOKES )
!  Cumulative single scatter source terms
!      DOUBLE PRECISION :: SS_CUMSOURCE_UP ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
!      DOUBLE PRECISION :: SS_CUMSOURCE_DN ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
!  Atmospheric attenuation before reflection
!      DOUBLE PRECISION :: ATTN_DB_SAVE ( MAX_GEOMETRIES )
!  Exact direct beam source terms
!      DOUBLE PRECISION :: EXACTDB_SOURCE ( MAX_GEOMETRIES, MAXSTOKES )
!  Cumulative direct bounce source terms
!      DOUBLE PRECISION :: DB_CUMSOURCE ( MAX_GEOMETRIES, MAXSTOKES, 0:MAXLAYERS )
!  Solar beam attenuation to BOA (required for exact DB calculation)
!      DOUBLE PRECISION :: BOA_ATTN ( MAX_GEOMETRIES )
!  Outgoing sphericity stuff - Whole and part-layer LOS transmittance factors
!      DOUBLE PRECISION :: UP_LOSTRANS    ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_LOSTRANS    ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: UP_LOSTRANS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_LOSTRANS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!  Whole and part-layer multipliers
!      DOUBLE PRECISION :: UP_MULTIPLIERS    ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_MULTIPLIERS    ( MAXLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: UP_MULTIPLIERS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!      DOUBLE PRECISION :: DN_MULTIPLIERS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )
!    Linearized output. Output required for later on.
!      DOUBLE PRECISION :: L_UP_LOSTRANS   (MAXLAYERS,     MAX_ATMOSWFS,MAX_GEOMETRIES)
!      DOUBLE PRECISION :: L_UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)
!  Solar beam attenuation to BOA (required for exact DB calculation)
!      DOUBLE PRECISION :: LC_BOA_ATTN(MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Arrays required at the Top level
!  ================================

!               Intent(In) To the Fourier routine

!  Input optical properties after delta-M scaling

      DOUBLE PRECISION :: DELTAU_VERT  ( MAXLAYERS )
      DOUBLE PRECISION :: PARTAU_VERT  ( MAX_PARTLAYERS )
      DOUBLE PRECISION :: OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Linearized input optical properties after delta-M scaling

      LOGICAL          :: DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Local input solar zenith angles by levels
!  ( Only required for refractive geometry attenuation of the solar beam)
!  These will be set internally if the refraction flag is set.

      DOUBLE PRECISION :: SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAX_SZANGLES )

!  Last layer to include Particular integral solution

      INTEGER          :: BEAM_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      DOUBLE PRECISION :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: LOCAL_CSZA     ( 0:MAXLAYERS, MAXBEAMS )

!  Solar beam attenuation

      DOUBLE PRECISION :: TRANS_SOLAR_BEAM ( MAXBEAMS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      DOUBLE PRECISION :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      DOUBLE PRECISION :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Reflectance flags

      LOGICAL         :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Local flags for the solution saving option

      LOGICAL         :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Local flags,  BVP telescoping enhancement

      LOGICAL         :: BVP_REGULAR_FLAG (0:MAXMOMENTS)

!  Masking for regular case. Required again for linearization

      !? INTEGER         :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      !? INTEGER         :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Telescoping initial flag (modified argument), Layer bookkeeping

      LOGICAL         :: DO_BVTEL_INITIAL
      INTEGER         :: BVTEL_FOURIER

!  Number of telescoped layers, active layers,  Size of BVP matrix 
      !INTEGER         :: NLAYERS_TEL, ACTIVE_LAYERS ( MAXLAYERS ), N_BVTELMATRIX_SIZE

!  Set up for band matrix compression
      !? INTEGER         :: BMAT_ROWMASK    ( MAXTOTAL, MAXTOTAL )
      !? INTEGER         :: BTELMAT_ROWMASK ( MAXTOTAL, MAXTOTAL )

!  Transmittance Setups
!  --------------------

!               Intent(In) To the Fourier routine

!  Discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      DOUBLE PRECISION :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      DOUBLE PRECISION :: T_UTUP_DISORDS(MAXSTREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION :: T_UTDN_DISORDS(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage

      DOUBLE PRECISION :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      DOUBLE PRECISION :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage 

      DOUBLE PRECISION :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Cumulative transmittance

      DOUBLE PRECISION :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!               Intent(In) To the Fourier routine

!  Forcing term multipliers (saved for whole atmosphere)

      DOUBLE PRECISION :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Partial layer multipliers

      DOUBLE PRECISION :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      DOUBLE PRECISION :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  LINEARIZED Arrays required at the Top level
!  ===========================================

!  Linearized Transmittance Setups
!  -------------------------------

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.
!  Version 2.8, dimensioning changed USER_LEVELS --> PARTLAYERS

      DOUBLE PRECISION :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTUP_DISORDS ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTDN_DISORDS ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0
!  Version 2.8, dimensioning changed USER_LEVELS --> PARTLAYERS

      DOUBLE PRECISION :: L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized transmittances, solar beam
!  Version 2.8, dimensioning changed USER_LEVELS --> PARTLAYERS

      DOUBLE PRECISION :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!mick fix 9/19/2017 - added LC_LEVELS_SOLARTRANS & LC_PARTIALS_SOLARTRANS

      DOUBLE PRECISION :: LC_LEVELS_SOLARTRANS   ( 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Beam multipliers
!  ---------------------------

!  Linearized whole layer multipliers

      DOUBLE PRECISION :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      DOUBLE PRECISION :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Thermal Setup outputs
!  ---------------------

!  Optical depth powers

      DOUBLE PRECISION :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      DOUBLE PRECISION :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!  Thermal coefficients, bookkeeping

      DOUBLE PRECISION :: THERMCOEFFS ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: TCOM1       ( MAXLAYERS, MAX_THERMAL_COEFFS )

!  Tranmsittance solutions

      DOUBLE PRECISION :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION :: T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Linearized optical depth powers

      DOUBLE PRECISION :: L_DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Linearized Thermal coefficients and bookkeeping

      DOUBLE PRECISION :: L_TCOM1       ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_THERMCOEFFS ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

!  Linearized Tranmsittance solutions

      DOUBLE PRECISION :: L_T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Vector FO setups
!  ----------------
!mick mod 9/19/2017 - added FO_STOKES_ATMOS, FO_STOKES_SURF, FO_COLUMNWF_ATMOS, FO_COLUMNWF_SURF (new output from vector FO code)

!  from VFO_LCS_MASTER_INTERFACE or VFO_MASTER_INTERFACE
!    SS/DB (solar), DTA/DTS (thermal) and total

      DOUBLE PRECISION :: FO_STOKES_SS    ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_DB    ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_STOKES_DTA   ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_DTS   ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION :: FO_STOKES_ATMOS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_STOKES_SURF  ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_STOKES       ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  from VFO_LCS_MASTER_INTERFACE: Column Weighting functions

      DOUBLE PRECISION :: FO_COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS,MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_COLUMNWF_DB &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION :: FO_COLUMNWF_DTA &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS,MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_COLUMNWF_DTS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION :: FO_COLUMNWF_ATMOS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION :: FO_COLUMNWF_SURF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_COLUMNWF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS,MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )

!  from VFO_LCS_MASTER_INTERFACE: Surface Jacobians (only DB and DTS)

      DOUBLE PRECISION :: FO_SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      DOUBLE PRECISION :: FO_SURFACEWF_DTS &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

      DOUBLE PRECISION :: FO_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  Additional output for the MSSTS requirements
!  --------------------------------------------

!  1/31/21. Version 2.8.3. Installed final version of DO_MSSTS code
!    ==> Additional Level SZA/VZA output (needed for the MSSTS)

      DOUBLE PRECISION :: FO_THETA_ALL ( 0:MAXLAYERS, MAX_GEOMETRIES )
      DOUBLE PRECISION :: FO_ALPHA     ( 0:MAXLAYERS, MAX_GEOMETRIES )

!  1/31/21. Version 2.8.3. Installed final version of DO_MSSTS code
!    ==> LOSTRANS_UP, needed for the MSST output, Upwelling option
!    ==> Added LOSTRANS_DN for the Downwelling Alternative.

      DOUBLE PRECISION :: FO_LOSTRANS_UP ( MAX_GEOMETRIES, MAXLAYERS )
      DOUBLE PRECISION :: FO_LOSTRANS_DN ( MAX_GEOMETRIES, MAXLAYERS )

      DOUBLE PRECISION :: FO_LC_LOSTRANS_UP ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: FO_LC_LOSTRANS_DN ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS )

!  Transmitted Fluxes for Adjusted Water-leaving
!  ---------------------------------------------

!  4/9/19. Additional output from the SFO Interface, for the sleave correction

      DOUBLE PRECISION :: FO_CUMTRANS    ( max_user_levels, MAX_GEOMETRIES )
      DOUBLE PRECISION :: FO_LC_CUMTRANS ( max_user_levels, MAX_GEOMETRIES, MAX_ATMOSWFS )

!  4/9/19. Surface leaving FO assignation (un-adjusted)

      DOUBLE PRECISION :: FO_SLTERM      (  MAXSTOKES, MAX_GEOMETRIES )
      DOUBLE PRECISION :: FO_LSSL_SLTERM (  MAXSTOKES, MAX_GEOMETRIES, MAX_SLEAVEWFS )

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  First Introduced 3/22/17 for LIDORT, based on VLIDORT code 
!    --  The LC Jacobian is an approximation.
!    --  Also depends on surface/sleave quantities, but these linearizations (LS/LSSL) are neglected

      DOUBLE PRECISION :: TRANS_ATMOS_FINAL      ( MAXBEAMS )
      DOUBLE PRECISION :: LC_TRANS_ATMOS_FINAL   ( MAXBEAMS, MAX_ATMOSWFS )
!      DOUBLE PRECISION :: LS_TRANS_ATMOS_FINAL   ( MAXBEAMS, MAX_SURFACEWFS )
!      DOUBLE PRECISION :: LSSL_TRANS_ATMOS_FINAL ( MAXBEAMS, MAX_SLEAVEWFS  )

!  Local variables
!  ---------------

!  Local quantities

      DOUBLE PRECISION :: CUMSOURCE_DB, SLTERM_LOCAL, LS_CUMSOURCE_DB, L_CUMSOURCE_DB, L_TRANS(MAX_ATMOSWFS)
      DOUBLE PRECISION :: LC_SLTERM_LOCAL(MAX_ATMOSWFS), LS_SLTERM_LOCAL(MAX_SLEAVEWFS), TFACTOR
      
!  Local flags for media problem and Planetary problem control
      
      LOGICAL          :: LOCAL_DO_ALBTRN_MEDIA(2)
      LOGICAL          :: LOCAL_DO_PLANETARY_PROBLEM
      DOUBLE PRECISION :: TRANS
      
!  Local flags

      LOGICAL ::          LOCAL_DO_NO_AZIMUTH
      LOGICAL ::          SAVE_DO_NO_AZIMUTH
      LOGICAL ::          LOCAL_ITERATION

!  inclusion flags

      LOGICAL ::          DO_INCLUDE_SURFACE
      LOGICAL ::          DO_INCLUDE_SURFEMISS
      LOGICAL ::          DO_INCLUDE_THERMEMISS

!  Fourier

      INTEGER ::          FOURIER
      INTEGER ::          N_FOURIERS

!  Misc. help

      INTEGER ::          OFF, O1, UA, UM, IB, G, G1, G2, L, Q, Q1, UTA, LUM, LUA, NMOMS, NWFS
      INTEGER ::          TESTCONV, LOCAL_N_USERAZM, STATUS_SUB
      INTEGER ::          IUNIT, SUNIT, FUNIT, RUNIT

!  For Convergence routines

      DOUBLE PRECISION :: AZM_ARGUMENT, DFC
      DOUBLE PRECISION :: AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXSTOKES )
      INTEGER ::          IBEAM_COUNT, IBEAM, NSOURCES
      LOGICAL ::          BEAM_ITERATION ( MAXBEAMS )
      INTEGER ::          BEAM_TESTCONV  ( MAXBEAMS )

!  Adjusted geometries. New, 2007.
!  -------------------------------

!  No longer required for Version 2.8, Have removed the Correction routines

!               Intent(Out) from the Adjust-geometry routine
!               Intent(In)  to   the Correction      routine
!      DOUBLE PRECISION :: USER_VZANGLES_ADJUST ( MAX_USER_VZANGLES )
!      DOUBLE PRECISION :: SZANGLES_ADJUST      ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
!      DOUBLE PRECISION :: USER_RELAZMS_ADJUST  ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
!      LOGICAL ::          ADJUST_SURFACE
!      DOUBLE PRECISION :: MODIFIED_ERADIUS

!  Transmittances and Fluxes for Water-leaving, Master routine (Mark 3)
!  --------------------------------------------------------------------
!    - Introduced 1/12/16, Validated 2/3/16. DISABLED Version 2.8.1.
!      DOUBLE PRECISION :: FLUX_DIFFUSE_FINAL  ( MAX_SZANGLES )

!  Helper variables
!  ----------------

!  flags

      LOGICAL ::          DO_NO_AZIMUTH
      LOGICAL ::          DO_ALL_FOURIER
      LOGICAL ::          DO_DIRECT_BEAM
!      LOGICAL ::          DO_CLASSICAL_SOLUTION  !!! removed for Version 2.8
      LOGICAL ::          DO_DBCORRECTION

!  Numbers

      INTEGER ::          NSTOKES_SQ
      INTEGER ::          NSTKS_NSTRMS
      INTEGER ::          NSTKS_NSTRMS_2
      INTEGER ::          NBEAMS
      INTEGER ::          NPARTICSOLS

!  Single scatter flux multipier

      DOUBLE PRECISION :: SS_FLUX_MULT

!  Solution bookkeeping

      INTEGER ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: DMAT ( MAXSTOKES, MAXSTOKES )
      INTEGER ::          GREEKMAT_INDEX ( 6 )
      LOGICAL ::          DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION :: FLUXVEC ( MAXSTOKES )
      DOUBLE PRECISION :: DFLUX ( MAXSTOKES )

!  output control
!    -- 1/31/21. Version 2.8.3. LOCAL_UM_START dropped

      INTEGER ::          N_OUT_STREAMS
      DOUBLE PRECISION :: OUT_ANGLES ( MAX_USER_STREAMS )

! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. SPECIALIST OUTPUT. Rayleigh Fourier output (TOA Upwelling only)
!      LOGICAL ::          SPECIAL_RAYF_OUTPUT
! #################### INNOVATIONS 5/5/20 ##############

!  level and partial control

      LOGICAL ::          DO_PARTLAYERS
      INTEGER ::          N_LAYERSOURCE_UP
      INTEGER ::          N_LAYERSOURCE_DN
      INTEGER ::          N_ALLLAYERS_UP
      INTEGER ::          N_ALLLAYERS_DN

!  Chapman factors
!mick fix 9/19/2017 - added PARTIAL_CHAPFACS

      DOUBLE PRECISION :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: PARTIAL_CHAPFACS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Local optical depths

      DOUBLE PRECISION :: TAUGRID_INPUT ( 0:MAXLAYERS )
      DOUBLE PRECISION :: TAUGRID ( 0:MAXLAYERS )

!  Contribution functions (TOA Upwelling only)
!  -------------------------------------------

!  Fourier component of Diffuse Field

      DOUBLE PRECISION ::  MS_CONTRIBS_F ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )

!  1/31/21. Version 2.8.3. Use type structure  variables directly.
!  Fourier-summed values
!      DOUBLE PRECISION ::  CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )
!  Single scatter
!      DOUBLE PRECISION :: SS_CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

!  Work type structures
!  --------------------

      TYPE(VLIDORT_Work_Miscellanous) :: Misc
      TYPE(VLIDORT_Work_Thermal)      :: Therm
      TYPE(VLIDORT_Work_Multiplier)   :: Mult
      !TYPE(VLIDORT_Work_Corrections)  :: Corr
      !TYPE(VLIDORT_Work_FirstOrder)   :: Fo

      TYPE(VLIDORT_LinWork_Miscellanous) :: LAC_Misc
      TYPE(VLIDORT_LinWork_Thermal)      :: L_Therm
      TYPE(VLIDORT_LinWork_Multiplier)   :: LC_Mult

!  Intermediate quantities

      DOUBLE PRECISION :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: PIMM_KM ( MAX_ALLSTRMS_P1 )

!  Local error handling

      LOGICAL ::          FAIL

!  Test variables
!      LOGICAL ::          DO_FDTEST =.FALSE.

!  This is the flag for a complete write-up of the inputs
!  1/31/21. Version 2.8.3. Now an input argument.
!      LOGICAL ::          DO_DEBUG_INPUT=.FALSE.
!      LOGICAL ::          DO_DEBUG_INPUT=.TRUE.

!  Initialize some variables
!  -------------------------

!mick debug
!write(*,*)
!write(*,*) 'in VLIDORT_LCS_MASTER'

!mick debug
!write(*,*) 'at 12'
!call system_mem_usage3(valueRSS)

!  Main status

      STATUS_CALCULATION = VLIDORT_SUCCESS
      STATUS_INPUTCHECK  = VLIDORT_SUCCESS

!mick fix 6/29/11 - initialize "Input checks"
!  Input checks

      NCHECKMESSAGES = 0
      CHECKMESSAGES  = ' '
      ACTIONS        = ' '

!  Model calculation

      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      TRACE_3 = ' '

!  Local user indices

      LUM = 1
      LUA = 1

!  Check input dimensions. New for Version 2.7
!  -------------------------------------------

!  regular

      CALL VLIDORT_CHECK_INPUT_DIMS &
      ( VLIDORT_FixIn, VLIDORT_ModIn, &
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  Linearized

      CALL VLIDORT_L_CHECK_INPUT_DIMS &
      ( VLIDORT_LinFixIn, VLIDORT_LinModIn, &
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  Fixed Boolean inputs
!  --------------------

      DO_FULLRAD_MODE        = VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE

      DO_THERMAL_EMISSION    = VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION
      DO_SURFACE_EMISSION    = VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION
      DO_PLANE_PARALLEL      = VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL

      DO_UPWELLING           = VLIDORT_FixIn%Bool%TS_DO_UPWELLING
      DO_DNWELLING           = VLIDORT_FixIn%Bool%TS_DO_DNWELLING

!      DO_QUAD_OUTPUT         = VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT  ! removed 7/7/16
      DO_TOA_CONTRIBS        = VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS

      DO_LAMBERTIAN_SURFACE  = VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE

!  1/31/21. Version 2.8.3.  Flag for calculating MSSTS output.

      DO_MSSTS               = VLIDORT_FixIn%Bool%TS_DO_MSSTS

      DO_SPECIALIST_OPTION_1 = VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_1
      DO_SPECIALIST_OPTION_2 = VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2
      DO_SPECIALIST_OPTION_3 = VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3

!  New 17 May 2012. Surface leaving
      DO_SURFACE_LEAVING     = VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC        = VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC

!  New Version 2.8. Water leaving flags.

      DO_WATER_LEAVING       = VLIDORT_FixIn%Bool%TS_DO_WATER_LEAVING   ! introduced 10/28/15
      DO_FLUORESCENCE        = VLIDORT_FixIn%Bool%TS_DO_FLUORESCENCE    ! introduced 10/28/15
      DO_TF_ITERATION        = VLIDORT_FixIn%Bool%TS_DO_TF_ITERATION    ! introduced 07/07/16

!  3/18/19. New for Version 2.8.1. Water leaving output flag.
      
      DO_WLADJUSTED_OUTPUT   = VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT

!  4/26/19. new for Version 2.8.1, introduced by R. Spurr
!    -- Control for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA

      DO_ALBTRN_MEDIA        = VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA
      
!  4/28/19. new for Version 2.8.1, introduced by R. Spurr
!    -- Control for Planetary problem.

      DO_PLANETARY_PROBLEM = VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM
      
!  TOA/BOA Illumination flags. 3/23/19 for Version 2.8.1
      
      DO_TOAFLUX    = VLIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION
      DO_BOAFLUX    = VLIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION

!  1/31/21. Version 2.8.3. (RTS 2/16/21). Flag for using an NSTOKES = 2 calculation for Fourier 0
!    -- This must be set by hand, it is not a configuration file read

      DO_FOURIER0_NSTOKES2 = VLIDORT_FixIn%Bool%TS_DO_FOURIER0_NSTOKES2

!  Fixed Control inputs
!  --------------------

!  Taylor parameter new, Version 2p7

      TAYLOR_ORDER     = VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER 

!  Numbers

      NSTOKES          = VLIDORT_FixIn%Cont%TS_NSTOKES
      NSTREAMS         = VLIDORT_FixIn%Cont%TS_NSTREAMS
      NMOMS            = 2*NSTREAMS - 1
      NLAYERS          = VLIDORT_FixIn%Cont%TS_NLAYERS
      NFINELAYERS      = VLIDORT_FixIn%Cont%TS_NFINELAYERS
      N_THERMAL_COEFFS = VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS

      VLIDORT_ACCURACY = VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY

      NLAYERS_NOMS     = VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS
      NLAYERS_CUTOFF   = VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF

!  1/31/21. Version 2.8.3. ASYMTX Tolerance variable

      ASYMTX_TOLERANCE = VLIDORT_FixIn%Cont%TS_ASYMTX_TOLERANCE

!  New for Version 2.8, Water-leaving iteration control. 7/7/16

      TF_MAXITER       = VLIDORT_FixIn%Cont%TS_TF_MAXITER
      TF_CRITERION     = VLIDORT_FixIn%Cont%TS_TF_CRITERION

!  TOA/BOA Illumination. 3/23/19 for Version 2.8.1
      
      TOAFLUX       = VLIDORT_FixIn%Cont%TS_TOA_ILLUMINATION
      BOAFLUX       = VLIDORT_FixIn%Cont%TS_BOA_ILLUMINATION

!  Fixed Beam inputs
!  -----------------

      FLUX_FACTOR      = VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR

!  Fixed User Value inputs
!  -----------------------

      N_USER_LEVELS    = VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS

!  Fixed Chapman Function inputs
!  -----------------------------

      HEIGHT_GRID(0:NLAYERS)      = VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0:NLAYERS)
      PRESSURE_GRID(0:NLAYERS)    = VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID(0:NLAYERS)
      TEMPERATURE_GRID(0:NLAYERS) = VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(0:NLAYERS)
      FINEGRID(1:NLAYERS)         = VLIDORT_FixIn%Chapman%TS_FINEGRID(1:NLAYERS)
      RFINDEX_PARAMETER           = VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER

!  Fixed Optical inputs
!  --------------------

      DELTAU_VERT_INPUT(1:NLAYERS) = VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(1:NLAYERS)
      GREEKMAT_TOTAL_INPUT(0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,1:NLAYERS,:) = &
               VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,1:NLAYERS,:)

!  New version 2.8, F-matrix input. 7/7/16
!mick fix 9/19/2017 - swapped layer & geo indices

      !FMATRIX_UP(:,1:NLAYERS,:) = VLIDORT_FixIn%Optical%TS_FMATRIX_UP(:,1:NLAYERS,:)
      !FMATRIX_DN(:,1:NLAYERS,:) = VLIDORT_FixIn%Optical%TS_FMATRIX_DN(:,1:NLAYERS,:)
      FMATRIX_UP(1:NLAYERS,:,:) = VLIDORT_FixIn%Optical%TS_FMATRIX_UP(1:NLAYERS,:,:)
      FMATRIX_DN(1:NLAYERS,:,:) = VLIDORT_FixIn%Optical%TS_FMATRIX_DN(1:NLAYERS,:,:)

      ALBEDO     = VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO

      THERMAL_BB_INPUT(0:NLAYERS)  = VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT(0:NLAYERS)
      SURFBB                = VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT

!  @@@ Rob Fix 1/31/11, TS_EMISSIVITY, TS_USER_EMISSIVITY are defined in
!                       in Structure VLIDORT_Sup%BRDF, so do not need to be copied here.
!      EMISSIVITY            = VLIDORT_FixIn%Optical%TS_EMISSIVITY
!      USER_EMISSIVITY       = VLIDORT_FixIn%Optical%TS_USER_EMISSIVITY

!  Replaced by LBBF facility. Version 2.7
!      LTE_DELTAU_VERT_INPUT(1:2,1:NLAYERS) = &
!        VLIDORT_FixIn%Optical%TS_LTE_DELTAU_VERT_INPUT(1:2,1:NLAYERS)
!      LTE_THERMAL_BB_INPUT(0:NLAYERS)      = &
!        VLIDORT_FixIn%Optical%TS_LTE_THERMAL_BB_INPUT(0:NLAYERS)

      ATMOS_WAVELENGTH        = VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH

!  Fixed Write inputs
!  ------------------

      DO_DEBUG_WRITE          = VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE

      DO_WRITE_INPUT          = VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT
      INPUT_WRITE_FILENAME    = VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME

      DO_WRITE_SCENARIO       = VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO
      SCENARIO_WRITE_FILENAME = VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME

      DO_WRITE_FOURIER        = VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER
      FOURIER_WRITE_FILENAME  = VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME

      DO_WRITE_RESULTS        = VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS
      RESULTS_WRITE_FILENAME  = VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME

!  Modified Boolean inputs
!  -----------------------

!  FOCORR and SSCORR Booleans. Completely reorganized for Version 2.8
!mick mod 9/19/2017 - DO_FOCORR_ALONE now defined internally
      !DO_FOCORR_ALONE         = VLIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE

!  FO flags

      DO_FOCORR               = VLIDORT_ModIn%MBool%TS_DO_FOCORR            !New 02 Jul 2013
      DO_FOCORR_EXTERNAL      = VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL
      DO_FOCORR_NADIR         = VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR
      DO_FOCORR_OUTGOING      = VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING
      DO_SSCORR_USEFMAT       = VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT

!  1/31/21. Version 2.8.3.  Flag removed
!      DO_SSCORR_TRUNCATION    = VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION

!  1/31/21. Version 2.8.3.  Add DO_CLASSICAL_SOLUTION flag

      DO_CLASSICAL_SOLUTION   = VLIDORT_ModIn%MBool%TS_DO_CLASSICAL_SOLUTION

!  Other Booleans

      DO_DOUBLE_CONVTEST      = VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST

      DO_SOLAR_SOURCES        = VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES
      DO_REFRACTIVE_GEOMETRY  = VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY
      DO_CHAPMAN_FUNCTION     = VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION

      DO_RAYLEIGH_ONLY        = VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY
      DO_DELTAM_SCALING       = VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING

      DO_SOLUTION_SAVING      = VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING
      DO_BVP_TELESCOPING      = VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING

      DO_USER_VZANGLES        = VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES
      DO_ADDITIONAL_MVOUT     = VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT
      DO_MVOUT_ONLY           = VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY
      DO_THERMAL_TRANSONLY    = VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY

!  1/31/21. Version 2.8.3.   Add DO_DOUBLET_GEOMETRY flag

      DO_OBSERVATION_GEOMETRY = VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY
      DO_DOUBLET_GEOMETRY     = VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY

!  Additional Control for Externalized Water-leaving input. Introduced 3/18/19 for Version 2.8.1

      DO_EXTERNAL_WLEAVE     = VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE

!  Modified Control inputs
!  -----------------------

      NGREEK_MOMENTS_INPUT = VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT

!  Modified Beam and User inputs
!  -----------------------------
!mick mod 1/5/2021 - added rewritten beam and user sections as in "vlidort_masters"

!  User levels & geometry specification height

      USER_LEVELS(1:N_USER_LEVELS)     = VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS)
      GEOMETRY_SPECHEIGHT              = VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT

!  Zero the local geometry numbers and arrays
!    -- 1/31/21. Version 2.8.3. (2/16/21). Newly added, with use of new doublet geometry

      N_USER_OBSGEOMS = 0 ; USER_OBSGEOMS     = ZERO
      N_USER_DOUBLETS = 0 ; USER_DOUBLETS     = ZERO
      N_SZANGLES      = 0 ; SZANGLES          = ZERO
      N_USER_VZANGLES = 0 ; USER_VZANGLES     = ZERO
      N_USER_RELAZMS  = 0 ; USER_RELAZMS      = ZERO

!  Geometries. Either Observational, Doublet or Lattice
!    * 1/31/21. Version 2.8.3. (2/16/21). section completely rewritten, new doublet goemetry option added

      IF ( DO_OBSERVATION_GEOMETRY ) THEN
        N_USER_OBSGEOMS                      = VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
        USER_OBSGEOMS(1:N_USER_OBSGEOMS,1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,1:3)
        N_SZANGLES                           = N_USER_OBSGEOMS
        N_USER_VZANGLES                      = N_USER_OBSGEOMS
        N_USER_RELAZMS                       = N_USER_OBSGEOMS
        SZANGLES     (1:N_USER_OBSGEOMS)     = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
        USER_VZANGLES(1:N_USER_OBSGEOMS)     = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        USER_RELAZMS (1:N_USER_OBSGEOMS)     = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
        N_SZANGLES                           = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
        SZANGLES(1:N_SZANGLES)               = VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:N_SZANGLES)
        N_USER_DOUBLETS                      = VLIDORT_ModIn%MUserVal%TS_N_USER_DOUBLETS
        USER_DOUBLETS(1:N_USER_DOUBLETS,1:2) = VLIDORT_ModIn%MUserVal%TS_USER_DOUBLETS(1:N_USER_DOUBLETS,1:2)
        N_USER_VZANGLES                      = N_USER_DOUBLETS
        N_USER_RELAZMS                       = N_USER_DOUBLETS
        USER_VZANGLES(1:N_USER_DOUBLETS)     = USER_DOUBLETS(1:N_USER_DOUBLETS,1)
        USER_RELAZMS (1:N_USER_DOUBLETS)     = USER_DOUBLETS(1:N_USER_DOUBLETS,2)
      ELSE
        N_SZANGLES                           = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
        SZANGLES(1:N_SZANGLES)               = VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:N_SZANGLES)
        N_USER_VZANGLES                      = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
        N_USER_RELAZMS                       = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
        USER_VZANGLES(1:N_USER_VZANGLES)     = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:N_USER_VZANGLES)
        USER_RELAZMS (1:N_USER_RELAZMS)      = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS)
      ENDIF

!      N_SZANGLES             = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
!      SZANGLES(1:N_SZANGLES) = VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:N_SZANGLES)
!      N_USER_RELAZMS                 = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
!      USER_RELAZMS(1:N_USER_RELAZMS) = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS)
!      N_USER_VZANGLES                  = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
!      USER_VZANGLES(1:N_USER_VZANGLES) = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:N_USER_VZANGLES)
!      USER_LEVELS(1:N_USER_LEVELS)     = VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS)
!      GEOMETRY_SPECHEIGHT              = VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT
!      N_USER_OBSGEOMS                      = VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
!      USER_OBSGEOMS(1:N_USER_OBSGEOMS,1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,1:3)

!  Modified Chapman Function inputs
!  --------------------------------

      !CHAPMAN_FACTORS     = VLIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS
      EARTH_RADIUS        = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS

!  Modified Optical inputs
!  -----------------------

      OMEGA_TOTAL_INPUT(1:NLAYERS) = VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS)

!  BRDF inputs
!  -----------
!mick mod 9/19/2017 - added IF conditions
!                   - note: emissivities left out of BRDF SURFACE block due to possibility
!                           of thermal being used in the Lambertian case

!  1/31/21. Version 2.8.3. BRDF copying only for Direct-bounce TERM.
!    -- Fourier BRDF copying now moved into Fourier loop.

      IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
        IF ( DO_USER_VZANGLES ) THEN
          EXACTDB_BRDFUNC (:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
            VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC(:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
        ENDIF
      ENDIF

!  @@@ Rob fix 1/31/11, Emissivities from earlier structure
!        Lambertian case was OK, as used internal definitions

      IF ( DO_SURFACE_EMISSION ) THEN
        EMISSIVITY(1:NSTOKES,1:NSTREAMS) = VLIDORT_Sup%BRDF%TS_EMISSIVITY(1:NSTOKES,1:NSTREAMS)
        IF ( DO_USER_VZANGLES ) THEN
          USER_EMISSIVITY(1:NSTOKES,1:N_USER_VZANGLES) = VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY(1:NSTOKES,1:N_USER_VZANGLES)
        ENDIF
      ENDIF

!  SLEAVE inputs  (This code introduced 17 May 2012)
!  -------------
!mick mod 9/19/2017 - added IF conditions

!  1/31/21. Version 2.8.3. SLEAVE copying only for Isotropic and direct TERM.
!    -- Fourier SLEAVE copying now moved into Fourier loop.

      IF ( DO_SURFACE_LEAVING ) THEN
        SLTERM_ISOTROPIC(1:NSTOKES,1:N_SZANGLES) = &
          VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(1:NSTOKES,1:N_SZANGLES)
        IF ( DO_USER_VZANGLES ) THEN
          SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
            VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
        ENDIF
      ENDIF

!  Fixed Linearized Control inputs
!  -------------------------------

      LAYER_VARY_FLAG  (1:NLAYERS) = VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG  (1:NLAYERS)
      LAYER_VARY_NUMBER(1:NLAYERS) = VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS)

      N_TOTALCOLUMN_WFS    = VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS
      N_TOTALPROFILE_WFS   = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
      N_SURFACE_WFS        = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS
      N_SLEAVE_WFS         = VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS

      N_TOTALSURFACE_WFS   = N_SURFACE_WFS + N_SLEAVE_WFS

      COLUMNWF_NAMES       = VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES
      PROFILEWF_NAMES      = VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES

!  Fixed Linearized Optical inputs
!  -------------------------------

      L_DELTAU_VERT_INPUT(1:N_TOTALCOLUMN_WFS,1:NLAYERS) = &
        VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(1:N_TOTALCOLUMN_WFS,1:NLAYERS)
      L_OMEGA_TOTAL_INPUT(1:N_TOTALCOLUMN_WFS,1:NLAYERS) = &
        VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(1:N_TOTALCOLUMN_WFS,1:NLAYERS)
      L_GREEKMAT_TOTAL_INPUT&
          (1:N_TOTALCOLUMN_WFS,0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,1:NLAYERS,:) = &
        VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT&
          (1:N_TOTALCOLUMN_WFS,0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,1:NLAYERS,:)

!  Version 2.8. Copy the linearized F-matrix inputs
!mick fix 9/19/2017 - swapped layer & geo indices

      !L_FMATRIX_UP ( 1:N_TOTALCOLUMN_WFS,:,1:NLAYERS,:) = &
      !  VLIDORT_LinFixIn%Optical%L_FMATRIX_UP(1:N_TOTALCOLUMN_WFS,:,1:NLAYERS,:)
      !L_FMATRIX_DN ( 1:N_TOTALCOLUMN_WFS,:,1:NLAYERS,:) = &
      !  VLIDORT_LinFixIn%Optical%L_FMATRIX_DN(1:N_TOTALCOLUMN_WFS,:,1:NLAYERS,:)
      L_FMATRIX_UP (1:N_TOTALCOLUMN_WFS,1:NLAYERS,:,:) = &
        VLIDORT_LinFixIn%Optical%TS_L_FMATRIX_UP(1:N_TOTALCOLUMN_WFS,1:NLAYERS,:,:)
      L_FMATRIX_DN (1:N_TOTALCOLUMN_WFS,1:NLAYERS,:,:) = &
        VLIDORT_LinFixIn%Optical%TS_L_FMATRIX_DN(1:N_TOTALCOLUMN_WFS,1:NLAYERS,:,:)

!  Modified Linearized Control inputs
!  ----------------------------------

      DO_COLUMN_LINEARIZATION  = VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION
      DO_PROFILE_LINEARIZATION = VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION
      DO_ATMOS_LINEARIZATION   = VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION

      DO_SURFACE_LINEARIZATION = VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION
      DO_LINEARIZATION         = VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION

      DO_SIMULATION_ONLY       = VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY

!  These 2.6 flags were superseded in Version 2.7 by LBBF flags
!      DO_LTE_LINEARIZATION    = VLIDORT_LinFixIn%Cont%TS_DO_LTE_LINEARIZATION
!      DO_SURFBB_LINEARIZATION = VLIDORT_LinFixIn%Cont%TS_DO_SURFBB_LINEARIZATION
      DO_ATMOS_LBBF   = VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF
      DO_SURFACE_LBBF = VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF

      DO_SLEAVE_WFS            = VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS

!  Linearized BRDF inputs
!  ----------------------
!mick mod 9/19/2017 - added IF conditions
!                   - note: emissivities left out of BRDF SURFACE block due to possibility
!                     of thermal being used in the Lambertian case

!  1/31/21. Version 2.8.3. BRDF copying only for Direct-bounce TERM.
!    -- Fourier BRDF copying now moved into Fourier loop.

      IF ( .NOT.DO_LAMBERTIAN_SURFACE .AND. DO_SURFACE_LINEARIZATION ) THEN
        IF ( DO_USER_VZANGLES ) THEN
          LS_EXACTDB_BRDFUNC(1:N_SURFACE_WFS,:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
            VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC&
                 (1:N_SURFACE_WFS,:,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
        ENDIF
      ENDIF

!  @@@ Rob fix 1/31/11, Linearized Emissivities from earlier structure
!                       were not copied --> wrong answers for BRDF cases
!        Lambertian case was OK, as used internal definitions

      IF ( DO_SURFACE_EMISSION ) THEN
        LS_EMISSIVITY(1:N_SURFACE_WFS,1:NSTOKES,1:NSTREAMS) = &
          VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY(1:N_SURFACE_WFS,1:NSTOKES,1:NSTREAMS)
        IF ( DO_USER_VZANGLES ) THEN
          LS_USER_EMISSIVITY(1:N_SURFACE_WFS,1:NSTOKES,1:N_USER_VZANGLES) = &
            VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY(1:N_SURFACE_WFS,1:NSTOKES,1:N_USER_VZANGLES)
        ENDIF
      ENDIF

!  Linearized SLEAVE inputs
!  ------------------------
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!mick mod 9/19/2017 - added IF conditions

!  1/31/21. Version 2.8.3. SLEAVE copying only for Isotropic and direct TERM.
!    -- Fourier SLEAVE copying now moved into Fourier loop.

      IF ( DO_SURFACE_LEAVING .AND. DO_SLEAVE_WFS ) THEN
        LSSL_SLTERM_ISOTROPIC(1:N_SLEAVE_WFS,1:NSTOKES,1:N_SZANGLES) = &
          VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC&
            (1:N_SLEAVE_WFS,1:NSTOKES,1:N_SZANGLES)
        IF ( DO_USER_VZANGLES ) THEN
          LSSL_SLTERM_USERANGLES(1:N_SLEAVE_WFS,1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
            VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES&
              (1:N_SLEAVE_WFS,1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
        ENDIF
      ENDIF

! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. Set Flag: Special Fourier-component output for Rayleigh + Planetary-problem TOA Upwelling situations
!      SPECIAL_RAYF_OUTPUT = &
!       ( DO_PLANETARY_PROBLEM .AND.DO_UPWELLING .AND..NOT.DO_FOCORR .AND.DO_RAYLEIGH_ONLY .AND. (USER_LEVELS(1).EQ.ZERO) )
! #################### INNOVATIONS 5/5/20 ##############

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  VLIDORT input debug

      IF (DO_DEBUG_INPUT) THEN
        CALL VLIDORT_DEBUG_INPUT_MASTER()
        IF (DO_COLUMN_LINEARIZATION .OR. DO_SURFACE_LINEARIZATION) THEN
          CALL VLIDORT_DEBUG_LIN_INPUT_MASTER()
        ENDIF
      ENDIF

!  initialize outputs
!  ------------------

!  Main outputs (Radiances and fluxes)
!   -- 1/31/21. Version 2.8.3. (2/16/21). Type-structure Stokes output now filled directly in Converge routines
!      STOKES(1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO

      MEANST_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO
      FLUX_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)   = ZERO
      DNMEANST_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)  = ZERO
      DNFLUX_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)    = ZERO

!  Special Media-property output. -- Introduced 4/26/19 R. Spurr. Pre-initialized here.
!     ** Output for User-angle streams, also fluxes. TRANSBEAM for the planetary problem.

      ALBMED_USER = ZERO ; ALBMED_FLUXES = ZERO
      TRNMED_USER = ZERO ; TRNMED_FLUXES = ZERO
      TRANSBEAM   = ZERO

!  4/28/19. Initialize the planetary problem outputs
!mick fix 8/20/2019 - switched to initializing the internal arrays

      VLIDORT_Out%Main%TS_PLANETARY_SBTERM    = ZERO
      VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM = ZERO

!  TF outputs, Rob fix, RT Solutions, 9/25/15, Superceded 2/3/16
!      MEANST_DIFFUSE_TF(1:N_SZANGLES,1:NSTOKES) = ZERO
!      FLUX_DIFFUSE_TF(1:N_SZANGLES,1:NSTOKES) = ZERO
!      DNMEANST_DIRECT_TF(1:N_SZANGLES,1:NSTOKES) = ZERO
!      DNFLUX_DIRECT_TF(1:N_SZANGLES,1:NSTOKES) = ZERO

!  Column weighting functions
!   -- 1/31/21. Version 2.8.3. (2/16/21). Type-structure Stokes output now filled directly in Converge routines
!      COLUMNWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO

!  Total mean intensity and flux weighting functions

      MEANST_DIFFUSE_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO
      FLUX_DIFFUSE_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO

!  Mean intensity and flux weighting functions.  Direct beam only.

      DNMEANST_DIRECT_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = ZERO
      DNFLUX_DIRECT_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = ZERO

!  Removed. 28 March 2014
!      !LTE column weighting functions
!      LTE_ATMOSWF(0:NLAYERS,1:N_USER_LEVELS,1:N_USER_VZANGLES,:) = ZERO

!  Surface weighting functions
!   -- 1/31/21. Version 2.8.3. (2/16/21). Type-structure Stokes output now filled directly in Converge routines
!      SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO

      MEANST_DIFFUSE_SURFWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO
      FLUX_DIFFUSE_SURFWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = ZERO

!mick fix 8/20/2019 - initialize linearized media-property & planetary problem output
 
!  Linearized Media-property output

      LC_ALBMED_USER(1:NSTOKES,1:N_USER_VZANGLES,1:MAX_ATMOSWFS) = ZERO
      LC_TRNMED_USER(1:NSTOKES,1:N_USER_VZANGLES,1:MAX_ATMOSWFS) = ZERO
      
      LC_ALBMED_FLUXES(1:NSTOKES,1:2,1:MAX_ATMOSWFS) = ZERO   
      LC_TRNMED_FLUXES(1:NSTOKES,1:2,1:MAX_ATMOSWFS) = ZERO   
      
      LC_TRANSBEAM(1:NSTOKES,1:N_SZANGLES,1:MAX_ATMOSWFS) = ZERO

!  Planetary problem output

      VLIDORT_LinOut%Col%TS_PLANETARY_TRANSTERM_COLWF(1:NSTOKES,:,1:MAX_ATMOSWFS) = ZERO
      VLIDORT_LinOut%Col%TS_PLANETARY_SBTERM_COLWF(1:MAX_ATMOSWFS) = ZERO

!  New 28 March 2014. BLACKBODY Linearization, Version 2.7

      ABBWFS_JACOBIANS(1:N_USER_LEVELS,:,:,:,:)  = ZERO
      ABBWFS_FLUXES   (1:N_USER_LEVELS,:,:,:,:)  = ZERO
      SBBWFS_JACOBIANS(1:N_USER_LEVELS,:,:,:)    = ZERO
      SBBWFS_FLUXES   (1:N_USER_LEVELS,:,:,:)    = ZERO

!  SS inputs
!  ---------

!  New 12 March 2012 --> IF SS results already available copy them. Modified flagging, Version 2.8

!  1/31/21. Version 2.8.3. (2/16/21). Local arrays no longer required. 
!    --Type structure arrays filled directly after VFO call, but zeroed here

      IF ( .NOT. DO_FOCORR_EXTERNAL ) THEN
         VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,:,1:NSTOKES,:) = ZERO
         VLIDORT_Sup%SS%TS_STOKES_DB(1:N_USER_LEVELS,:,1:NSTOKES)   = ZERO
      ENDIF

!  Linearized SS inputs
!  --------------------

      IF ( .not.DO_FOCORR_EXTERNAL ) THEN
         VLIDORT_LinSup%SS%Col%TS_COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,:,1:NSTOKES,:) = zero
         VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,:,1:NSTOKES)   = zero
         VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:,1:NSTOKES) = zero
      ENDIF

!  Define DO_FOCORR_ALONE flag before input checks
!mick fix 1/5/2021 - moved this line back to before VLIDORT_CHECK_INPUT

      DO_FOCORR_ALONE = ( .NOT.DO_FULLRAD_MODE .AND. DO_FOCORR )

!  Check input
!  -----------

!    Major revision of I/O output list, 25 October 2012.
!     ---- Observational Geometry control, New, 25 October 2012
!     ---- Automatic setting of NBEAMS, N_USER_VZANGLES, N_USER_RELAZMS, DO_USER_VZANGLES

!  %% DO_FULLRAD_MODE argument added (First line). R. Spurr, 05 March 2013
!  %%   Needed to ensure MS-only output in all cases when flagged
!  Add Flag DO_WLADJUSTED_OUTPUT (Water-leaving output). 3/18/19 for Version 2.8.1. 

!    Major revision of I/O output list, 07 july 2016, for Version 2.8.
!      Inputs are arranged by type.

!  1/31/21. Version 2.8.3. 
!   -- Add DO_CLASSICAL_SOLUTION to list arguments (Line 7). Also add NSTOKES (line 10)
!   -- Add DO_MSSTS flag to this list. Line 9, final Boolean. Arguments rearranged
!   -- Add DOUBLET_GEOMETRY to this list
!   -- Green's function only for solar sources with NSTOKES = 1 or 3. Add NSTOKES input.
!   -- Several restrictions on use of MSSTS option, need to be checked

      CALL VLIDORT_CHECK_INPUT ( &
           DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,                     & ! Input
           DO_THERMAL_TRANSONLY, DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_FLUORESCENCE,      & ! Input
           DO_WATER_LEAVING, DO_TF_ITERATION, DO_WLADJUSTED_OUTPUT, DO_EXTERNAL_WLEAVE,           & ! Boolean 
           DO_FULLRAD_MODE, DO_RAYLEIGH_ONLY, DO_SSCORR_USEFMAT, DO_DIRECT_BEAM,                  & ! Boolean
           DO_USER_VZANGLES, DO_DELTAM_SCALING, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY,     & ! Boolean
           DO_FOCORR, DO_FOCORR_EXTERNAL, DO_FOCORR_ALONE, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING,   & ! Boolean
           DO_REFRACTIVE_GEOMETRY, DO_CHAPMAN_FUNCTION, DO_PLANE_PARALLEL, DO_CLASSICAL_SOLUTION, & ! Boolean
           DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,            & ! Boolean
           DO_TOA_CONTRIBS, DO_ALL_FOURIER, DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3, DO_MSSTS, & ! Boolean
           TAYLOR_ORDER, NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, TF_MAXITER,                   & ! Integer
           N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS, N_USER_OBSGEOMS,      & ! Integer
           N_OUT_STREAMS, NLAYERS_NOMS, NLAYERS_CUTOFF, NGREEK_MOMENTS_INPUT, & ! Integer
           SZANGLES, USER_VZANGLES, USER_RELAZMS, USER_OBSGEOMS, USER_LEVELS, & ! Floating point
           OUT_ANGLES, TF_CRITERION, EARTH_RADIUS, GEOMETRY_SPECHEIGHT,       & ! Floating point
           HEIGHT_GRID, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,              & ! Floating point
           STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )                 ! Exception handling

!  Exception handling

      IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = VLIDORT_SERIOUS
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
        VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
        VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
        RETURN
      ELSE IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = VLIDORT_WARNING
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
!  ########## Rob Change, 3/28/2011 ################
!  Program will execute - these outputs are set at the end.
!          VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
!          VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
!          VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
!  ########## Rob Change, 3/28/2011 ################
      ENDIF

!  Extended input variables check

      IF ( DO_LINEARIZATION ) THEN

        CALL VLIDORT_L_CHECK_INPUT ( &
          DO_SIMULATION_ONLY, N_SURFACE_WFS,                 &
          DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION, &
          DO_ATMOS_LINEARIZATION, DO_SURFACE_LINEARIZATION,  &
          STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

        IF ( STATUS_SUB .EQ. VLIDORT_SERIOUS ) THEN
          STATUS_INPUTCHECK = VLIDORT_SERIOUS
          VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
          VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
          VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
          VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
          RETURN
        ELSE IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
          STATUS_INPUTCHECK = VLIDORT_WARNING
          VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
!  ########## Rob Change, 3/28/2011 ################
!  Program will execute - these outputs are set at the end.
!          VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
!          VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
!          VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
!  ########## Rob Change, 3/28/2011 ################
        ENDIF

      ENDIF

!  4/26/19. Additional checks on the Isotropic-illumination input.
!  ---------------------------------------------------------------

!  9/25/19. Bug: VLIDORT_Out%Status%TS_STATUS_INPUTCHECK was not set

!  5/5/20. Version 2.8.1 Upgrades 
!   ==> Relax condition on Rayleigh only, and on FOCORR_NADIR
!              ( Planetary problem works for Aerosols and FOCORR_NADIR )
!   ==> Now, only fails for (dark-surface, Lambertian case, no thermal)

!  Here is the older code..........(Pre 5/5/20)
!      IF (DO_PLANETARY_PROBLEM .OR. DO_ALBTRN_MEDIA(1) .OR. DO_ALBTRN_MEDIA(2) ) then
!         IF ( .NOT.DO_LAMBERTIAN_SURFACE .OR. DO_SURFACE_LEAVING .OR. DO_THERMAL_EMISSION &
!              .OR. DO_FOCORR .OR. (ALBEDO .NE.ZERO) ) then
!            STATUS_INPUTCHECK = VLIDORT_SERIOUS
!            VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK ! 9/25/19 this was not set
!            VLIDORT_Out%Status%TS_NCHECKMESSAGES = 1
!            VLIDORT_Out%Status%TS_CHECKMESSAGES(1) = 'Media_problem input check not valid'
!            VLIDORT_Out%Status%TS_ACTIONS(1)       = 'Check thermal/Rayleigh/Lambertian flags for this option'
!            RETURN
!         ENDIF
!      ENDIF

!  Here is the New Code..........(5/5/20 Upgrade)

      IF ( DO_PLANETARY_PROBLEM .OR. DO_ALBTRN_MEDIA(1) .OR. DO_ALBTRN_MEDIA(2) ) then
         IF ( (.NOT. DO_LAMBERTIAN_SURFACE) .OR. DO_SURFACE_LEAVING .OR. DO_THERMAL_EMISSION &
               .OR. (ALBEDO .NE. ZERO) .OR. (DO_FOCORR .AND. DO_FOCORR_OUTGOING) ) then
            STATUS_INPUTCHECK = VLIDORT_SERIOUS
            VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK ! 9/25/19 this was not set
            VLIDORT_Out%Status%TS_NCHECKMESSAGES   = 2
            VLIDORT_Out%Status%TS_CHECKMESSAGES(1:2) = 'Media/Planetary problem: input check not valid'
            VLIDORT_Out%Status%TS_ACTIONS(1)       = 'Either: Turn off Thermal_Emission/FOCORR_OUTGOING/Surface_Leaving flags'
            VLIDORT_Out%Status%TS_ACTIONS(2)       = 'Or    : Make sure Lambertian surface flag and set albedo to zero'
            RETURN
         ENDIF
      ENDIF

!  Bookkeeping and Preparation for Fourier call
!  --------------------------------------------

!   -- section moved here from before CHECK_INPUT call, 1/31/21. Version 2.8.3.

!  Proxy output

      FOURIER_SAVED(1:N_SZANGLES) = 0
      N_GEOMETRIES                = 0

!  1/31/21. Version 2.8.3. Zero All offsets

      SZA_OFFSETS(1:N_SZANGLES)   = 0
      VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES) = 0
      SZD_OFFSETS(1:N_SZANGLES) = 0

!  Single scatter correction: flux multiplier

      SS_FLUX_MULT = FLUX_FACTOR / PI4

!  Number of sources

      IF ( DO_SOLAR_SOURCES ) THEN
        NSOURCES = N_SZANGLES
        NBEAMS   = N_SZANGLES
      ELSE
        NSOURCES = 1
        NBEAMS   = 1
      ENDIF

!  1/31/21. Version 2.8.3. Added local post-processing control.

      PPSTREAM_MASK = 0 ; N_PPSTREAMS = N_USER_VZANGLES
      IF ( DO_OBSERVATION_GEOMETRY ) N_PPSTREAMS = 1
      DO IBEAM = 1, NBEAMS
        IF ( DO_OBSERVATION_GEOMETRY ) THEN
          PPSTREAM_MASK(1,IBEAM) = IBEAM
        ELSE
          DO UM = 1, N_PPSTREAMS
            PPSTREAM_MASK(UM,IBEAM) = UM
          ENDDO
        ENDIF
     ENDDO

!  Write input variables
!  ---------------------

!  Version 2.8 revision 7/8/16, remove QUAD_OUTPUT and CLASSICAL SOLUTION, revise input list.
!  open file, call standard input write, close file

      IF ( DO_WRITE_INPUT ) THEN
        IUNIT = VLIDORT_INUNIT
        OPEN(IUNIT,FILE=INPUT_WRITE_FILENAME,STATUS='REPLACE')

        CALL VLIDORT_WRITEINPUT ( IUNIT, &
           DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY, &
           DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_RAYLEIGH_ONLY,    &
           DO_NO_AZIMUTH, DO_CLASSICAL_SOLUTION, DO_MSSTS,                 &
           DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_LAMBERTIAN_SURFACE,      &
           DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_DIRECT_BEAM,       &
           NSTREAMS, NLAYERS, NBEAMS, N_USER_VZANGLES, N_USER_RELAZMS,     &
           N_USER_LEVELS, NGREEK_MOMENTS_INPUT, VLIDORT_ACCURACY, FLUX_FACTOR )

        IF ( DO_LINEARIZATION ) THEN
          CALL VLIDORT_L_WRITEINPUT ( &
            IUNIT, DO_SIMULATION_ONLY, DO_PROFILE_LINEARIZATION, &
            DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS,          &
            N_TOTALPROFILE_WFS, DO_SURFACE_LINEARIZATION,        &
            DO_ATMOS_LBBF, DO_SURFACE_LBBF, PROFILEWF_NAMES,     &
            COLUMNWF_NAMES )
        ENDIF

        CLOSE(IUNIT)
      ENDIF

!  Get derived inputs
!  ==================

!  Miscellaneous and layer input.
!  Version 2.8 revision 7/8/16, I/O list. USER_VZANGLES_ADJUST no longer present
!  Version 2.8 revision 3/1/17, I/O list for FOCORR flags
!mick fix 9/19/2017 - added DO_FOCORR_NADIR & DO_FOCORR_OUTGOING to input

!  1/31/21. Version 2.8.3. Introduce DO_DOUBLET_GEOMETRY input
!   -- re-ordered first 4 lines of input
!   -- DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY are not needed.

      CALL VLIDORT_DERIVE_INPUT ( &
        DO_FULLRAD_MODE, DO_UPWELLING, DO_DNWELLING, DO_FOCORR,                  & ! Input Boolean
        DO_FOCORR_ALONE, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_SSCORR_USEFMAT, & ! Input Boolean
        DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_RAYLEIGH_ONLY,              & ! Input Boolean
        DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY,          & ! Input Boolean
        DO_THERMAL_TRANSONLY, DO_SPECIALIST_OPTION_2, DO_SOLUTION_SAVING,        & ! Input Boolean
        DO_BVP_TELESCOPING, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_NO_AZIMUTH,   & ! Input Boolean
        DO_DOUBLE_CONVTEST, DO_ALL_FOURIER, DO_DIRECT_BEAM, DO_DBCORRECTION,     & ! Input Boolean
        NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SZANGLES,  NLAYERS_NOMS,    & ! Input Integer
        NGREEK_MOMENTS_INPUT, N_USER_RELAZMS, N_USER_VZANGLES,  N_OUT_STREAMS,   & ! Input Integer
        SZANGLES, USER_VZANGLES, USER_LEVELS,                                    & ! Input Floating point
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,              & ! Input Floating point
        COS_SZANGLES, SIN_SZANGLES, QUAD_STREAMS,                                & ! Output
        QUAD_WEIGHTS,QUAD_STRMWTS, QUAD_HALFWTS, QUAD_SINES, QUAD_ANGLES,        & ! Output
        DO_MSMODE_VLIDORT, NMOMENTS, NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,   & ! Output
        NSTKS_NSTRMS, NSTKS_NSTRMS_2, NBEAMS, NPARTICSOLS, NSTOKES_SQ,           & ! Output
        FLUXVEC, DMAT, MUELLER_INDEX, GREEKMAT_INDEX, DO_REAL_EIGENSOLVER,       & ! Output
        BVP_REGULAR_FLAG, LAYER_MAXMOMENTS, DO_LAYER_SCATTERING,                 & ! Output
        N_CONVTESTS, N_CONV_STREAMS, N_DIRECTIONS, WHICH_DIRECTIONS,             & ! Output
        USER_STREAMS, USER_SINES, USER_SECANTS, PARTLAYERS_OUTFLAG,              & ! Output
        PARTLAYERS_OUTINDEX, LEVELMASK_UP, LEVELMASK_DN,                         & ! Output
        DO_PARTLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,     & ! Output
        N_LAYERSOURCE_UP, N_LAYERSOURCE_DN, N_ALLLAYERS_UP, N_ALLLAYERS_DN,      & ! Output
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, TAUGRID_INPUT,            & ! Output
        STATUS_SUB, MESSAGE )

!  If there's no azimuth dependence, just do one value in azimuth loop
!mick fix 1/5/2021 - moved azimuth & offsets for indexing blocks after VLIDORT_DERIVE_INPUT

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_N_USERAZM = 1
      ELSE
        LOCAL_N_USERAZM = N_USER_RELAZMS
      ENDIF

!  Save some offsets for indexing geometries
!   This section revised for the Observational Geometry option

!  1/31/21. Version 2.8.3. Add Offsets for the Doublet Geometry option
!    -- rearrange code for better logic

      IF ( DO_OBSERVATION_GEOMETRY ) THEN
        N_VIEWING    = N_USER_OBSGEOMS
        N_GEOMETRIES = N_USER_OBSGEOMS
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
        N_VIEWING    = N_USER_VZANGLES
        N_GEOMETRIES = NSOURCES * N_VIEWING
        DO IBEAM = 1, NBEAMS
           SZD_OFFSETS(IBEAM) = N_VIEWING * ( IBEAM - 1 )
        END DO
      ELSE
        N_VIEWING    = N_USER_VZANGLES * LOCAL_N_USERAZM
        N_GEOMETRIES = NSOURCES * N_VIEWING
        DO IBEAM = 1, NBEAMS
          SZA_OFFSETS(IBEAM) = N_VIEWING * ( IBEAM - 1 )
          DO UM = 1, N_USER_VZANGLES
            VZA_OFFSETS(IBEAM,UM) = SZA_OFFSETS(IBEAM) + LOCAL_N_USERAZM * (UM - 1)
          END DO
        END DO
      ENDIF

!  Set thermal MS flag
!mick fix 9/19/2017 - changed def of DO_MSMODE_THERMAL flag based on implementation
!                     of new FO code.  When DO_FOCORR set, both solar AND THERMAL
!                     direct now come from the FO code.

      !DO_MSMODE_THERMAL = (.NOT.DO_FULLRAD_MODE) .AND. ( DO_SURFACE_EMISSION .AND.DO_THERMAL_EMISSION )
      DO_MSMODE_THERMAL = DO_MSMODE_VLIDORT .AND. ( DO_SURFACE_EMISSION .AND. DO_THERMAL_EMISSION )

!  Exception handling

!  Rob Fix 3/28/2011. Output from DERIVE_INPUTS are checks, not execution failures
!  Old code---------
!      IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
!        TRACE_1 = ''
!        TRACE_2 = 'Derive_Input Call in VLIDORT_LCS_MASTER'
!        TRACE_3 = ' ** VLIDORT_LCS_MASTER'
!        STATUS_INPUTCHECK = VLIDORT_WARNING
!        VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
!        VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
!        VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!        VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
!        VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!      ENDIF
!  New code-----------
      IF ( STATUS_SUB .EQ. VLIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = VLIDORT_WARNING
        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        NCHECKMESSAGES = NCHECKMESSAGES + 1
        CHECKMESSAGES(NCHECKMESSAGES)  = TRIM(ADJUSTL(MESSAGE))
        ACTIONS(NCHECKMESSAGES) = ' Action taken in VLIDORT_DERIVE_INPUT to set internal default'
      ENDIF

!  Rob Fix, 3/17/15. Serious Error check on use of FO calculation, must exit
!   --> Plan for Version 2.8 is to relax this condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IF (DO_FO_CALC .AND. DO_PARTLAYERS) THEN
!        STATUS_INPUTCHECK = VLIDORT_SERIOUS
!        VLIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
!        NCHECKMESSAGES = NCHECKMESSAGES + 1
!        CHECKMESSAGES(NCHECKMESSAGES) = ' Internal SS calculation using FO code: ONLY LAYER-BOUNDARY output (no partials)'
!        ACTIONS(NCHECKMESSAGES)       = ' Set USER-LEVEL output only for layer boundaries'
!        VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
!        VLIDORT_Out%Status%TS_CHECKMESSAGES  = CHECKMESSAGES
!        VLIDORT_Out%Status%TS_ACTIONS        = ACTIONS
!        RETURN
!      ENDIF

!  Geometry adjustment
!  -------------------

!  Section removed, Version 2.8
!mick mod 1/5/2021 - moved this geometry section (currently commented out) to
!                    before VLIDORT_CHAPMAN call

!  Adjust surface condition
!      ADJUST_SURFACE = .FALSE.
!      IF ( DO_SSCORR_OUTGOING ) THEN
!        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) ADJUST_SURFACE = .TRUE.
!      ENDIF
!  Perform adjustment
!      MODIFIED_ERADIUS = EARTH_RADIUS + GEOMETRY_SPECHEIGHT
!mick hold - 9/26/2012
!      IF ( DO_SOLAR_SOURCES ) THEN
!      IF ( .NOT. DO_OBSERVATION_GEOMETRY ) THEN
!        CALL MULTI_OUTGOING_ADJUSTGEOM                                &
!         ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS,         & ! Input
!           N_USER_VZANGLES,   N_SZANGLES,   N_USER_RELAZMS,           & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE,    & ! Input
!           USER_VZANGLES, SZANGLES, USER_RELAZMS,                     & ! Input
!           USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST,& ! Output
!           FAIL, MESSAGE, TRACE_1 )                                     ! Output
!      ELSE
!        CALL OBSGEOM_OUTGOING_ADJUSTGEOM                               &
!         ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS,          & ! Input
!           N_USER_VZANGLES,   N_SZANGLES,   N_USER_RELAZMS,            & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE,     & ! Input
!           USER_VZANGLES, SZANGLES, USER_RELAZMS,                      & ! Input
!           USER_VZANGLES_ADJUST, SZANGLES_ADJUST, USER_RELAZMS_ADJUST, & ! Output
!           FAIL, MESSAGE, TRACE_1 )                                      ! Output
!      ENDIF
!      ELSE
!        CALL LOSONLY_OUTGOING_ADJUSTGEOM                           &
!         ( MAX_USER_VZANGLES, N_USER_VZANGLES,                     & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE, & ! Input
!           USER_VZANGLES,                                          & ! Input
!           USER_VZANGLES_ADJUST,                                   & ! Output
!           FAIL, MESSAGE, TRACE_1 )                                  ! Output
!      ENDIF
!  Update exception handling. October 2010 2p4RTC
!      if ( fail ) then
!        TRACE_2 = ' Failure in multi_outgoing_adjustgeom'
!        TRACE_3 = ' ** VLIDORT_LCS_MASTER '
!        STATUS_CALCULATION = VLIDORT_SERIOUS
!        VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
!        VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
!        VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!        VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
!        VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!        RETURN
!      ENDIF

!  Chapman function calculation
!  ----------------------------

!  Calling statement revised I/O, Version 2.8 7/7/16.

!mick fix - comment out 1st IF
!mick fix 9/19/2017 - added input N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES
!                   - added output PARTIAL_CHAPFACS
!                   - moved the call to VLIDORT_CHAPMAN from before VLIDORT_DERIVE_INPUT
!                     to here to resolve an I/O contradiction

      !IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN

          CALL VLIDORT_CHAPMAN ( &
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,            & !  Input
            NLAYERS, N_SZANGLES, FINEGRID, SZANGLES,              & !  Input
            N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, & !  Input
            EARTH_RADIUS, RFINDEX_PARAMETER,                      & !  Input
            HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,         & !  Input
            CHAPMAN_FACTORS, PARTIAL_CHAPFACS, SZA_LOCAL_INPUT,   & !  Output
            SUN_SZA_COSINES, FAIL, MESSAGE, TRACE_1 )               !  Output

          IF (FAIL) THEN
            TRACE_2 = 'Direct call in VLIDORT_LCS_MASTER'
            TRACE_3 = ' ** VLIDORT_LCS_MASTER '
            STATUS_CALCULATION = VLIDORT_SERIOUS
            VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
            VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
          ENDIF

        ENDIF
      !ENDIF

!  #################
!  Set up operations
!  #################

!  Setups for Fourier = 0
!  ======================

!  Each call is foolowed by a packing routine

!  1a.  VLIDORT_MISCSETUPS      : Delta-M, average-secant formulation, transmittances
!  1b.  VLIDORT_LAC_MISCSETUPS  : Linearizations of Item 1a
!  2a.  THERMAL_SETUP           : Coefficients, direct multipliers
!  2b.  THERMAL_SETUP_PLUS      : Coefficients, direct multipliers AND LINEARIZATIONS
!  3a. LC_EMULT_MASTER /        : Beam source function multipliers.  Not required for the
!  3b. LC_EMULT_MASTER_OBSGEO     Full SS calculation in outgoing mode

!  1. VLIDORT MISCSETUPS
!  ---------------------

!  1/30/08. Telescoping setup is done in DERIVE_INPUTS (unlike LIDORT scalar code)
!  7/22/16. Revision of I/O lists for Version 2.8.

!  With linearization.

      IF ( DO_ATMOS_LINEARIZATION ) THEN

!mick fix 9/19/2017 - the following arguments added to facilitate correction of direct flux:
!                     to input  - DO_SOLAR_SOURCES, PARTIAL_CHAPFACS
!                     to output - DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED, 
!                                 LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS

        CALL VLIDORT_MISCSETUPS ( &
          DO_SOLAR_SOURCES, DO_DELTAM_SCALING, DO_PLANE_PARALLEL,        & ! Input flags
          DO_REFRACTIVE_GEOMETRY,                                        & ! Input flags
          DO_PARTLAYERS, DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY,      & ! Input flags
          DO_SOLUTION_SAVING, DO_SPECIALIST_OPTION_3, DO_TOA_CONTRIBS,   & ! Input flags
          NSTOKES, NLAYERS, NSTREAMS, N_USER_VZANGLES, N_USER_LEVELS,    & ! Input Numbers
          NBEAMS, NMOMENTS, N_PARTLAYERS, NLAYERS_CUTOFF,                & ! Input Numbers
          MUELLER_INDEX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,         & ! Input level/Stokes indices
          PARTLAYERS_OUTFLAG,  PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,   & ! Input Partlayers
          QUAD_STREAMS, USER_SECANTS, COS_SZANGLES, SUN_SZA_COSINES,     & ! Input streams, SZA cosines
          OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT,    & ! Input Optical
          TAUGRID_INPUT, CHAPMAN_FACTORS, PARTIAL_CHAPFACS,              & ! Input Chapman
          DELTAU_VERT, PARTAU_VERT, OMEGA_TOTAL, GREEKMAT_TOTAL,         & ! Output from DELTAMSCALE
          TAUGRID, DELTAU_SLANT, DELTAU_SLANT_UNSCALED,                  & ! Output from DELTAMSCALE
          PARTAU_SLANT_UNSCALED, LEVELS_SOLARTRANS,                      & ! Output from DELTAMSCALE
          PARTIALS_SOLARTRANS, SOLARBEAM_BOATRANS, TRUNC_FACTOR, FAC1,   & ! Output from DELTAMSCALE
          OMEGA_GREEK,                                                   & ! Output from SSALBINIT
          DO_REFLECTED_DIRECTBEAM, BEAM_CUTOFF, TRANS_SOLAR_BEAM,        & ! Output QSPREP
          INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA,                     & ! output QSPREP
          T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                & ! Output PREPTRANS (Discrete Ords.)
          T_DELT_MUBAR, T_UTUP_MUBAR, T_UTDN_MUBAR,                      & ! Output PREPTRANS (Solar beams)
          T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                      & ! Output PREPTRANS (User-vza)
          CUMTRANS, ITRANS_USERM )                                         ! Output PREPTRANS (auxiliary)       

!  1/31/21. Version 2.8.3. Need additional variable ITRANS_USERM to be packed
!    -- Argument list rearranged to follow more natural logic

        CALL VLIDORT_PACK_MISC ( &
          NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS,               & ! Input Numbers
          N_USER_VZANGLES, N_SZANGLES, NMOMENTS,                  & ! Input Numbers
          DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT,                 & ! packed Optical
          LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, OMEGA_GREEK,    & ! packed Optical
          T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,         & ! packed Trans D.O.
          T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, CUMTRANS,     & ! packed Trans User
          BEAM_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, & ! packed Solar
          INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,            & ! packed Average Secant
          T_UTDN_MUBAR, ITRANS_USERM, LOCAL_CSZA,                 & ! packed Trans Solar
          Misc )                                                    ! Output structure to be packed

!mick fix 9/19/2017 - DO_SOLAR_SOURCES NOT ADDED to input
        CALL VLIDORT_LAC_MISCSETUPS ( &
          DO_USER_VZANGLES, DO_DELTAM_SCALING, DO_PLANE_PARALLEL,           & ! Input flags
          DO_SOLUTION_SAVING, DO_ATMOS_LINEARIZATION,                       & ! Input flags 
          NSTOKES, NLAYERS, NSTREAMS, N_USER_VZANGLES, NBEAMS,              & ! Input numbers
          NMOMENTS, NSTOKES_SQ, N_PARTLAYERS, PARTLAYERS_LAYERIDX,          & ! Input numbers/partials
          OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                          & ! Input Optical props
          N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,            & ! Input Lin control 
          DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, DELTAU_SLANT_UNSCALED,    & ! Input derived properties
          PARTAU_SLANT_UNSCALED, LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,    & ! Input derived properties
          OMEGA_GREEK, TRUNC_FACTOR, FAC1,                                  & ! Input derived properties
          MUELLER_INDEX, QUAD_STREAMS, USER_SECANTS,                        & ! Input streams/Muller
          T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                   & ! Input Dis.Ord. Trans.
          BEAM_CUTOFF, T_DELT_MUBAR, T_UTDN_MUBAR, AVERAGE_SECANT,          & ! Input beam Trans.
          TRANS_SOLAR_BEAM, T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,       & ! Input User Trans.
          L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT, & ! Input Linearized Optical props.
          L_DELTAU_VERT, L_OMEGA_TOTAL, L_GREEKMAT_TOTAL, L_DELTAU_SLANT,        & ! Output Linearized scaled properties
          DO_SCATMAT_VARIATION, L_TRUNC_FACTOR, L_OMEGA_GREEK,                   & ! Output Linearized scaled properties
          LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, & ! Output linearized beam stuff
          LC_LEVELS_SOLARTRANS,  LC_PARTIALS_SOLARTRANS,                         & ! Output linearized beam stuff
          LC_SOLARBEAM_BOATRANS, LC_SOLARBEAM_ATRANS,                            & ! Output linearized beam stuff
          L_T_DELT_DISORDS, L_T_UTDN_DISORDS, L_T_UTUP_DISORDS,                  & ! output linearized DisOrd Trans.
          L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM )                         ! Output lineraized User Trans.

        CALL VLIDORT_PACK_LAC_MISC ( &
          NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS,                       & ! Input
          N_SZANGLES, N_USER_VZANGLES, NMOMENTS, N_TOTALCOLUMN_WFS,       & ! Input
          L_DELTAU_VERT, L_OMEGA_GREEK, LC_INITIAL_TRANS,                 & ! Input
          L_T_DELT_DISORDS, L_T_UTDN_DISORDS, L_T_UTUP_DISORDS,           & ! Input
          L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,                 & ! Input
          LC_AVERAGE_SECANT,  LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,           & ! Input
          LC_LEVELS_SOLARTRANS, LC_PARTIALS_SOLARTRANS,                   & ! Input
          LAC_Misc )                                                        ! Output

!  2. THERMAL SETUPS
!  -----------------

!  7/8/16, 7/22/16. Revision of I/O lists for Version 2.8.

        IF ( DO_THERMAL_EMISSION ) THEN

          CALL THERMAL_SETUP_PLUS ( &
            DO_USER_VZANGLES, DO_UPWELLING, DO_DNWELLING,                        & ! Flags
            DO_PARTLAYERS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,              & ! Flags
            DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,          & ! Linearization control
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES,            & ! Numbers basic
            THERMAL_BB_INPUT, USER_STREAMS,                                      & ! thermal input, streams
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,         & ! Level control
            OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT, L_OMEGA_TOTAL, L_DELTAU_VERT, & ! Input optical+linearized
            T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                            & ! Input transmittances
            L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,                      & ! Input linearized transmittances
            THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,                        & ! output thermal setups
            T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN,            & ! output thermal direct solutions
            L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_TCOM1,                & ! output Linearized thermal setups
            L_T_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_UP, L_T_UT_DIRECT_DN )     ! output Linearized thermal direct solutions

          CALL VLIDORT_PACK_THERM ( &
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, & ! Input Numbers
            THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Input Thermal setup
            T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN, & ! Input Thermal solutions
            Therm )                                                     ! Output

          CALL VLIDORT_PACK_L_THERM ( &
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, N_TOTALCOLUMN_WFS, & ! Input
            L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_TCOM1,                        & ! Input
            L_T_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_UP, L_T_UT_DIRECT_DN,            & ! Input
            L_Therm )                                                                      ! Output

        ENDIF

!  3a/b. EMULT_MASTERS.
!  --------------------

!  7/8/16, Revision of I/O lists for Version 2.8.

!mick fix 3/30/2015 - modified if condition
      !IF ( DO_SOLAR_SOURCES  ) THEN
!Rob fix 3/1/2017 - modified if condition
      !IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN

        IF ( DO_SOLAR_SOURCES .AND. DO_USER_VZANGLES) THEN
          IF (.NOT.DO_FOCORR_ALONE) THEN

!  Lattice or Doublet geometry
!  ---------------------------

            IF ( .NOT. DO_OBSERVATION_GEOMETRY ) THEN
              CALL EMULT_MASTER ( &
                DO_UPWELLING, DO_DNWELLING,                                   & ! Input flags
                NLAYERS, NBEAMS, N_USER_VZANGLES, TAYLOR_ORDER, N_PARTLAYERS, & ! Input numbers
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,  & ! Input Level control
                USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                       & ! Input optical, streams   
                BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,      & ! Input solar beam
                T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, ITRANS_USERM,       & ! Input User-stream Trans.
                EMULT_HOPRULE, SIGMA_M, SIGMA_P,                              & ! Output Multipliers
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )                  ! Output Multipliers

!  1/31/21. Version 2.8.3. Need additional variables SIGMA_M, SIGMA_P to be packed

              CALL VLIDORT_PACK_MULT ( &
                NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES,             & ! Input
                SIGMA_M, SIGMA_P, EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN, & ! Input
                Mult )                                                            ! Output

              CALL LC_EMULT_MASTER ( &
                DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, TAYLOR_ORDER,       & ! Input flags + Taylor
                NLAYERS, NBEAMS, N_USER_VZANGLES, N_PARTLAYERS, N_TOTALCOLUMN_WFS, & ! Input numbers
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,       & ! Input Level control
                USER_SECANTS, DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,             & ! Input optical, streams   
                BEAM_CUTOFF, T_DELT_MUBAR, EMULT_HOPRULE, SIGMA_M, SIGMA_P,        & ! Input solar beam + Multipliers
                T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, ITRANS_USERM,            & ! Input User-stream Trans.
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                      & ! Input Multipliers
                LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,              & ! Input linearized Beam stuff
                LC_T_UTDN_MUBAR, L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,   & ! Input linearized transmittanaces
                LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN )           ! Output linearized Mulitpliers

              CALL VLIDORT_PACK_LC_MULT ( &
                NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, N_TOTALCOLUMN_WFS, & ! Input
                LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN,              & ! Input
                LC_Mult )                                                                ! Output

!  Observational geometry
!  ----------------------

            ELSE

              CALL EMULT_MASTER_OBSGEO ( &
                DO_UPWELLING, DO_DNWELLING,                                  & ! Input flags
                NLAYERS, NBEAMS, TAYLOR_ORDER, N_PARTLAYERS,                 & ! Input numbers
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Input Level control
                USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                      & ! Input optical, streams   
                BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,     & ! Input solar beam
                T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, ITRANS_USERM,      & ! Input User-stream Trans.
                EMULT_HOPRULE, SIGMA_M, SIGMA_P,                             & ! Output Multipliers
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )                 ! Output Multipliers

!  1/31/21. Version 2.8.3. Need additional variables SIGMA_M, SIGMA_P to be packed

              CALL VLIDORT_PACK_MULT ( &
                NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES,             & ! Input
                SIGMA_M, SIGMA_P, EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN, & ! Input
                Mult )                                                            ! Output

              CALL LC_EMULT_MASTER_OBSGEO ( &
                DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, TAYLOR_ORDER,      & ! Input flags + Taylor
                NLAYERS, NBEAMS, N_PARTLAYERS, N_TOTALCOLUMN_WFS,                 & ! Input numbers
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input Level control
                USER_SECANTS, DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,            & ! Input optical, streams   
                BEAM_CUTOFF, T_DELT_MUBAR, EMULT_HOPRULE, SIGMA_M, SIGMA_P,       & ! Input solar beam + Multipliers
                T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, ITRANS_USERM,           & ! Input User-stream Trans.
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                     & ! Input Multipliers
                LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,             & ! Input linearized Beam stuff
                LC_T_UTDN_MUBAR, L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,  & ! Input linearized transmittanaces
                LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN )          ! Output linearized Mulitpliers

              CALL VLIDORT_PACK_LC_MULT ( &
                NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES, N_TOTALCOLUMN_WFS, & ! Input
                LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN,              & ! Input
                LC_Mult )                                                                ! Output

            ENDIF
          ENDIF
        ENDIF

!  No linearization
!  ================

      ELSE

!mick fix 9/19/2017 - the following arguments added to facilitate correction of direct flux:
!                     to input  - DO_SOLAR_SOURCES, PARTIAL_CHAPFACS
!                     to output - DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED, 
!                                 LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS

        CALL VLIDORT_MISCSETUPS ( &
          DO_SOLAR_SOURCES, DO_DELTAM_SCALING, DO_PLANE_PARALLEL,        & ! Input flags
          DO_REFRACTIVE_GEOMETRY,                                        & ! Input flags
          DO_PARTLAYERS, DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY,      & ! Input flags
          DO_SOLUTION_SAVING, DO_SPECIALIST_OPTION_3, DO_TOA_CONTRIBS,   & ! Input flags
          NSTOKES, NLAYERS, NSTREAMS, N_USER_VZANGLES, N_USER_LEVELS,    & ! Input Numbers
          NBEAMS, NMOMENTS, N_PARTLAYERS, NLAYERS_CUTOFF,                & ! Input Numbers
          MUELLER_INDEX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,         & ! Input level/Stokes indices
          PARTLAYERS_OUTFLAG,  PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,   & ! Input Partlayers
          QUAD_STREAMS, USER_SECANTS, COS_SZANGLES, SUN_SZA_COSINES,     & ! Input streams, SZA cosines
          OMEGA_TOTAL_INPUT, DELTAU_VERT_INPUT, GREEKMAT_TOTAL_INPUT,    & ! Input Optical
          TAUGRID_INPUT, CHAPMAN_FACTORS, PARTIAL_CHAPFACS,              & ! Input Chapman
          DELTAU_VERT, PARTAU_VERT, OMEGA_TOTAL, GREEKMAT_TOTAL,         & ! Output from DELTAMSCALE
          TAUGRID, DELTAU_SLANT, DELTAU_SLANT_UNSCALED,                  & ! Output from DELTAMSCALE
          PARTAU_SLANT_UNSCALED, LEVELS_SOLARTRANS,                      & ! Output from DELTAMSCALE
          PARTIALS_SOLARTRANS, SOLARBEAM_BOATRANS, TRUNC_FACTOR, FAC1,   & ! Output from DELTAMSCALE
          OMEGA_GREEK,                                                   & ! Output from SSALBINIT
          DO_REFLECTED_DIRECTBEAM, BEAM_CUTOFF, TRANS_SOLAR_BEAM,        & ! Output QSPREP
          INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA,                     & ! output QSPREP
          T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                & ! Output PREPTRANS (Discrete Ords.)
          T_DELT_MUBAR, T_UTUP_MUBAR, T_UTDN_MUBAR,                      & ! Output PREPTRANS (Solar beams)
          T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                      & ! Output PREPTRANS (User-vza)
          CUMTRANS, ITRANS_USERM )                                         ! Output PREPTRANS (auxiliary)       

!  1/31/21. Version 2.8.3. Need additional variable ITRANS_USERM to be packed
!    -- Argument list rearranged to follow more natural logic

      CALL VLIDORT_PACK_MISC ( &
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS,               & ! Input Numbers
        N_USER_VZANGLES, N_SZANGLES, NMOMENTS,                  & ! Input Numbers
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT,                 & ! packed Optical
        LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, OMEGA_GREEK,    & ! packed Optical
        T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,         & ! packed Trans D.O.
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, CUMTRANS,     & ! packed Trans User
        BEAM_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, & ! packed Solar
        INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,            & ! packed Average Secant
        T_UTDN_MUBAR, ITRANS_USERM, LOCAL_CSZA,                 & ! packed Trans Solar
        Misc )                                                    ! Output structure to be packed

!  2. THERMAL SETUPS
!  -----------------

!  7/8/16 . Revision of I/O lists for Version 2.8.

        IF ( DO_THERMAL_EMISSION ) THEN

          CALL THERMAL_SETUP ( &
            DO_USER_VZANGLES, DO_UPWELLING, DO_DNWELLING,                & ! Flags
            DO_PARTLAYERS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,      & ! Flags
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES,    & ! Numbers basic
            THERMAL_BB_INPUT, USER_STREAMS,                              & ! thermal input, streams
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Level control
            OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT,                       & ! Input optical
            T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                    & ! Input transmittances
            THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,                & ! output thermal setups
            T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN )     ! output thermal direct solutions

          CALL VLIDORT_PACK_THERM ( &
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_VZANGLES, & ! Input Numbers
            THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Input Thermal setup
            T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN, & ! Input Thermal solutions
            Therm )                                                     ! Output

        END IF

!  3a/b. EMULT_MASTERS.
!  --------------------

!  7/8/16 . Revision of I/O lists for Version 2.8.

!mick fix 3/30/2015 - modified if condition
      !IF ( DO_SOLAR_SOURCES ) THEN
!Rob fix 3/1/2017 - modified if condition
      !IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN

        IF ( DO_SOLAR_SOURCES .AND. DO_USER_VZANGLES) THEN
          IF (.NOT.DO_FOCORR_ALONE) THEN

!  Lattice or Doublet geometry
!  ---------------------------

            IF ( .NOT. DO_OBSERVATION_GEOMETRY ) THEN
              CALL EMULT_MASTER ( &
                DO_UPWELLING, DO_DNWELLING,                                   & ! Input flags
                NLAYERS, NBEAMS, N_USER_VZANGLES, TAYLOR_ORDER, N_PARTLAYERS, & ! Input numbers
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,  & ! Input Level control
                USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                       & ! Input optical, streams   
                BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,      & ! Input solar beam
                T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, ITRANS_USERM,       & ! Input User-stream Trans.
                EMULT_HOPRULE, SIGMA_M, SIGMA_P,                              & ! Output Multipliers
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )                  ! Output Multipliers

!  Observational geometry
!  ----------------------

            ELSE
              CALL EMULT_MASTER_OBSGEO ( &
                DO_UPWELLING, DO_DNWELLING,                                  & ! Input flags
                NLAYERS, NBEAMS, TAYLOR_ORDER, N_PARTLAYERS,                 & ! Input numbers
                PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Input Level control
                USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                      & ! Input optical, streams   
                BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,     & ! Input solar beam
                T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, ITRANS_USERM,      & ! Input User-stream Trans.
                EMULT_HOPRULE, SIGMA_M, SIGMA_P,                             & ! Output Multipliers
                EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )                 ! Output Multipliers
            ENDIF

!  1/31/21. Version 2.8.3. Need additional variables SIGMA_M, SIGMA_P to be packed

            CALL VLIDORT_PACK_MULT ( &
              NLAYERS, N_PARTLAYERS, N_SZANGLES, N_USER_VZANGLES,             & ! Input
              SIGMA_M, SIGMA_P, EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN, & ! Input
              Mult )                                                            ! Output

          ENDIF
        ENDIF

!  End clause linearization

      ENDIF

!  WATER-LEAVING TRANSMITTANCE CALCULATION
!  =======================================

!  - 04/09/19. Radical new departure. Now moved to the BVPROBLEM as an adjusted Backsub calculation
!              Scaling of SLTERMS by TRANS_ATMOS_FINAL is now done either inside the Fourier routine (MS)
!              or below for the direct term contribution (After Fourier call for Fourier 0)
!              The  LIDORT_TRANSFLUX_MASTER routine has been removed.

!  New section for Version 2.7a and 2.8
!  - 09/25/15. First programmed by R. Spurr for Version 2.7a, RT Solutions Inc.
!  - 12/24/15. Drop SL_Isotropic Constraint (if non-Isotropy, additional SL terms need transmittance scaling)
!  - 02/03/16. Mark 1 and Mark2 codes (Dark Surface Result). COMMENTED OUT IN FAVOR of Mark 3
!  - 07/08/16. Mark 3 code given iteration control (3 new inputs).

!  HERE IS THE OLD CODE.....
!      IF ( DO_SURFACE_LEAVING .AND. DO_WATER_LEAVING .AND. DO_SL_ISOTROPIC ) THEN
!      IF ( DO_SURFACE_LEAVING .AND. DO_WATER_LEAVING ) THEN
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         Mark 1, Mark2 code, Dark Surface Result. COMMENTED OUT 2/3/16
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         CALL VLIDORT_TRANSFLUX & 
!           ( NSTOKES, NSTREAMS, NLAYERS, NBEAMS, NMOMENTS, NSTREAMS_2,           & ! Input
!             NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,         & ! Input
!             FLUX_FACTOR, FLUXVEC, COS_SZANGLES, MUELLER_INDEX, DMAT, DFLUX,     & ! Input
!             QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS, LOCAL_CSZA, & ! Input
!             DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, DELTAU_VERT, OMEGA_GREEK, & ! Input
!             BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT,      & ! Input
!             MEANST_DIFFUSE_TF, FLUX_DIFFUSE_TF, DNMEANST_DIRECT_TF, DNFLUX_DIRECT_TF,     & ! InOut
!             STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                               ! Output
!       write(*,*)'MK2 Result (Dark)   = ',0,FLUX_DIFFUSE_TF(1,1),FLUX_DIFFUSE_TF(1,1) / flux_factor / local_csza(nlayers,1)
!  debug output. 17 december 2015.
!     It is the last entry that is correct.....
!         write(779,'(A,A)')'   SZA    COSSZA   F0      Fdirect      Tdirect      Fdiffuse   ',&
!                                                    '  TDiffuse     Fcomplete    Tcomplete  [T* = F*/F0/cossza]'
!         do ibeam = 1, nbeams
!!            write(*,'(a,i2,2(1p3e15.6,2x))')'Transflux',ibeam,MEANST_DIFFUSE_TF(ibeam,1:3),FLUX_DIFFUSE_TF(ibeam,1:3)
!             write(779,'(F8.3,F9.5,F6.2,2x,1p6e13.5)')szangles(ibeam),local_csza(nlayers,ibeam),flux_factor,&
!              DNFLUX_DIRECT_TF(ibeam,1),DNFLUX_DIRECT_TF(ibeam,1)/flux_factor/local_csza(nlayers,ibeam),&
!              FLUX_DIFFUSE_TF(ibeam,1)-DNFLUX_DIRECT_TF(ibeam,1),&
!              (FLUX_DIFFUSE_TF(ibeam,1)-DNFLUX_DIRECT_TF(ibeam,1))/flux_factor/local_csza(nlayers,ibeam),&
!              FLUX_DIFFUSE_TF(ibeam,1),FLUX_DIFFUSE_TF(ibeam,1)/flux_factor/local_csza(nlayers,ibeam)
!         enddo
!  error handling
!         IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
!            TRACE_3 = ' Called by VLIDORT_MASTER ' ; STATUS_CALCULATION = VLIDORT_SERIOUS
!            VLIDORT_Out%Status%TS_MESSAGE = MESSAGE ; VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!            VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2 ; VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!            VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION ; RETURN
!         ENDIF
!  Scale = Flux/F0/mu0 = Atmospheric Transmittance (Dark Result)
!         TRANS_ATMOS(1:nbeams) = FLUX_DIFFUSE_TF(1:nbeams,1) / flux_factor / local_csza(nlayers,1:nbeams)
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         Mark 3 code, Self-Consistent Result, by iteration.   2/3/16
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Very important to set the surface flag
!         DO_INCLUDE_SURFACE = .true.
!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination
!  --- Also need the discrete ordinate transmittances  
!         CALL VLIDORT_TRANSFLUX_MASTER &
!        ( DO_TF_ITERATION, TF_MAXITER, TF_CRITERION, DO_TOAFLUX, TOAFLUX,     & ! Input TF control, TOAISO (3/23/19 new)
!          NSTOKES, NSTREAMS, NLAYERS, NBEAMS, NMOMENTS, NSTREAMS_2,                 & ! Input Numbers
!          NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,               & ! Input Numbers
!          FLUX_FACTOR, FLUXVEC, COS_SZANGLES, MUELLER_INDEX, DMAT, DFLUX,           & ! Input Bookkeeping
!          DO_INCLUDE_SURFACE, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,            & ! Input surface 1/8/16 
!          DO_LAMBERTIAN_SURFACE, ALBEDO, BRDF_F, BRDF_F_0,               & ! Input surface 1/8/16
!          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_F_0,        & ! Input surface 1/8/16
!          QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS, LOCAL_CSZA,       & ! Input Quadrature
!          DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, DELTAU_VERT, OMEGA_GREEK,       & ! Input Optical
!          T_DELT_DISORDS, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT, & ! Input solar beam, 3/23/19 new
!          TRANS_ATMOS_FINAL, FLUX_DIFFUSE_FINAL,                                    & ! Output
!          STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                                     ! Output
!  Exception handling
!            IF ( STATUS_SUB == VLIDORT_SERIOUS ) THEN
!               TRACE_3 = 'VLIDORT_TRANSFLUX_MASTER failed, VLIDORT_MASTER'
!               STATUS_CALCULATION = VLIDORT_SERIOUS
!               VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
!               VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
!               VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!               VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
!               VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!               RETURN
!            ENDIF
!  Scale the Isotropic term
!  Add the same factor for the linearized terms (if present). 12/18/15, 2/3/16
!         do ibeam = 1, nbeams
!           TFACTOR = TRANS_ATMOS_FINAL(ibeam)
!           SLTERM_ISOTROPIC(1,IBEAM) = SLTERM_ISOTROPIC(1,IBEAM) * TFACTOR
!           if ( DO_SLEAVE_WFS ) THEN
!             DO Q = 1, N_SLEAVE_WFS
!               LSSL_SLTERM_ISOTROPIC(Q,1,IBEAM) =  LSSL_SLTERM_ISOTROPIC(Q,1,IBEAM) * TFACTOR
!             enddo
!           endif
!         enddo
!  Scale the Non-isotropic terms, 24 December 2015, if flagged
!  Add the same factor for the linearized terms (if present). 12/18/15, 2/3/16
!         if ( .NOT. DO_SL_ISOTROPIC ) THEN
!           do ibeam = 1, nbeams
!             TFACTOR = TRANS_ATMOS_FINAL(ibeam)
!             SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IBEAM) = &
!             SLTERM_USERANGLES(1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IBEAM) * TFACTOR
!             SLTERM_F_0(0:NMOMS,1:NSTOKES,1:NSTREAMS,IBEAM) = &
!             SLTERM_F_0(0:NMOMS,1:NSTOKES,1:NSTREAMS,IBEAM) * TFACTOR
!             USER_SLTERM_F_0(0:NMOMS,1:NSTOKES,1:N_USER_VZANGLES,IBEAM) = &
!             USER_SLTERM_F_0(0:NMOMS,1:NSTOKES,1:N_USER_VZANGLES,IBEAM) * TFACTOR
!             if ( DO_SLEAVE_WFS ) THEN
!               DO Q = 1, N_SLEAVE_WFS
!                 LSSL_SLTERM_USERANGLES(Q,1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IBEAM) = &
!                 LSSL_SLTERM_USERANGLES(Q,1:NSTOKES,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IBEAM) * TFACTOR
!                 LSSL_SLTERM_F_0(Q,0:NMOMS,1:NSTOKES,1:NSTREAMS,IBEAM) = &
!                 LSSL_SLTERM_F_0(Q,0:NMOMS,1:NSTOKES,1:NSTREAMS,IBEAM) * TFACTOR
!                 LSSL_USER_SLTERM_F_0(Q,0:NMOMS,1:NSTOKES,1:N_USER_VZANGLES,IBEAM) = &
!                 LSSL_USER_SLTERM_F_0(Q,0:NMOMS,1:NSTOKES,1:N_USER_VZANGLES,IBEAM) * TFACTOR
!               ENDDO
!             endif
!           enddo
!         endif
!  End of Transmittance calculation
!      ENDIF

!  SINGLE-SCATTER & DIRECT-BOUNCE CALCULATIONS
!  ===========================================

!  Version 2.8, Major revision, 7/8/16

!    - SSCORR and DBCORRECTION routines completely removed
!    - Calculations only done using the FO code, Version 1.5 (replaces FO Version 1.4 added 7/2/13)
!    - Rob Fix 3/17/15. Exception handling updated for serious error.
!    - Rob Fix 7/08/16. Master interface updated for FO 1.5, includes FMATRIX inputs.

!  Not required if no solar sources, and no user-angles
!mick fix 9/19/2017 - deactivated the DO_SOLAR_SOURCES IF condition due to the
!                     enhanced nature of the new internal FO code
!mick mod 9/19/2017 - added FO_STOKES_ATMOS, FO_STOKES_SURF, FO_COLUMNWF_ATMOS,
!                     FO_COLUMNWF_SURF (new output from vector FO code)
!mick fix 1/5/2021  - added DO_FOCORR_EXTERNAL to IF condition

      !IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_USER_VZANGLES ) THEN
          IF ( DO_FOCORR .AND. .NOT.DO_FOCORR_EXTERNAL ) THEN

!   -- 4/9/19. Added output FO Surface-leaving assignation + cumulative transmittances,, and linearizations.
!              Added input, water-leaving control
             
!  1/31/21. Version 2.8.3. DO_MSSTS option final installation
!    ==> MSST situations: Add LOSTRANS_UP/DN, LC_LOSTRANS_UP/DN, THETA_ALL, ALPHA to output list
!    ==> Last 3 lines reorganized output list.

!  1/31/21. Version 2.8.3. Add Doublet geometry flag and include offsets

!mick debug
!write(*,*) 'at 13'
!call system_mem_usage3(valueRSS)

             CALL VFO_LCS_MASTER_INTERFACE ( &
                DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                         & ! Input Sources flags
                DO_PLANE_PARALLEL, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_SSCORR_USEFMAT,          & ! Input SS control flags
                DO_DELTAM_SCALING, DO_UPWELLING, DO_DNWELLING, DO_PARTLAYERS,                       & ! Input RT Control flags
                DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY, DO_LAMBERTIAN_SURFACE,                & ! input RT Control flags
                DO_SURFACE_LEAVING, DO_WATER_LEAVING, DO_SL_ISOTROPIC,                              & ! Input Surface flags
                DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS,                   & ! Input Jacobian flags
                NSTOKES, NLAYERS, NFINELAYERS, NGREEK_MOMENTS_INPUT, SZD_OFFSETS, SZA_OFFSETS, VZA_OFFSETS, & ! Input nums/offsets
                N_SZANGLES, SZANGLES, N_USER_VZANGLES, USER_VZANGLES, N_USER_RELAZMS, USER_RELAZMS, & ! Input geometry
                N_USER_LEVELS, LEVELMASK_UP, LEVELMASK_DN, N_PARTLAYERS,                            & ! Input levels  control
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,    & ! Input partial control
                N_TOTALCOLUMN_WFS, N_SLEAVE_WFS, N_SURFACE_WFS, N_TOTALSURFACE_WFS,                 & ! Input numbers (Jacobians)
                EARTH_RADIUS, HEIGHT_GRID, SS_FLUX_MULT, FLUXVEC,                                   & ! Input Flux/Heights
                DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,                         & ! Inputs (Optical - Regular)
                DELTAU_VERT, FMATRIX_UP, FMATRIX_DN, TRUNC_FACTOR, THERMAL_BB_INPUT,                & ! Inputs (Optical - Regular)
                ALBEDO, EXACTDB_BRDFUNC, SURFBB, USER_EMISSIVITY,                                   & ! Inputs (Optical - Surface)
                SLTERM_ISOTROPIC, SLTERM_USERANGLES,                                                & ! Inputs (Optical - Surface)
                L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,                   & ! Inputs (Optical - Lin Atmos)
                L_DELTAU_VERT, L_TRUNC_FACTOR, L_FMATRIX_UP, L_FMATRIX_DN,                          & ! Inputs (Optical - Lin Atmos)
                LS_EXACTDB_BRDFUNC, LS_USER_EMISSIVITY,                                             & ! Inputs (Optical - Lin Surf)
                LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES,                                      & ! Inputs (Optical - Lin Surf)
                FO_STOKES_SS, FO_STOKES_DB, FO_STOKES_DTA, FO_STOKES_DTS,                           & ! Output - Stokes vectors
                FO_COLUMNWF_SS,  FO_COLUMNWF_DB, FO_COLUMNWF_DTA, FO_COLUMNWF_DTS,                  & ! Output - Column Jacobians
                FO_SURFACEWF_DB, FO_SURFACEWF_DTS,                                                  & ! Output - Surface Jacobians
                FO_STOKES_ATMOS, FO_STOKES_SURF, FO_STOKES,                                         & ! Output - Stokes composites
                FO_COLUMNWF_ATMOS, FO_COLUMNWF_SURF, FO_COLUMNWF, FO_SURFACEWF,                     & ! Output - Jacobian composites
                FO_CUMTRANS, FO_LOSTRANS_UP, FO_LOSTRANS_DN, FO_THETA_ALL, FO_ALPHA, FO_SLTERM,     & ! Output Auxiliary
                FO_LC_CUMTRANS, FO_LSSL_SLTERM, FO_LC_LOSTRANS_UP, FO_LC_LOSTRANS_DN,               & ! Output - Auxiliary
                FAIL, MESSAGE, TRACE_1, TRACE_2 )                                                     ! Output

!  Exception handling

            IF ( FAIL ) THEN
               TRACE_3 = 'VFO_LCS_MASTER_INTERFACE failed, VLIDORT_LCS_MASTER'
               STATUS_CALCULATION = VLIDORT_SERIOUS
               VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
               VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
               VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
               VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
               VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
               RETURN
            ENDIF

!  Copy FO results to VLIDORT arrays
!mick fix 9/19/2017  - modified argument passing for SURFACEWF_DB
!                      from FO_SURFACEWF_DB to FO_SURFACEWF
!                    - changed 1st dim of SURFACEWF_DB & FO_SURFACEWF
!                      from N_SURFACE_WFS to N_TOTALSURFACE_WFS
!mick mod 9/19/2017  - added IF conditions 
!mick note 3/22/2017 - Important! STOKES_SS, STOKES_DB, COLUMNWF_SS, COLUMNWF_DB, &
!                      SURFACEWF_DB contain BOTH solar AND thermal direct terms when
!                      computing in the crossover region!
!mick fix 1/5/2021   - added defining of VLIDORT_Sup STOKES_SS & STOKES_DB and VLIDORT_LinSup COLUMNWF_SS, COLUMNWF_DB,
!                      & SURFACEWF_DB

            IF ( DO_UPWELLING ) THEN
               VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,UPIDX) = &
                           FO_STOKES_ATMOS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,UPIDX)
               VLIDORT_Sup%SS%TS_STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)       = &
                           FO_STOKES_SURF (1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)

               IF ( DO_COLUMN_LINEARIZATION ) THEN
                  VLIDORT_LinSup%SS%Col%TS_COLUMNWF_SS&
                                           (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,UPIDX) = &
                          FO_COLUMNWF_ATMOS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,UPIDX)
                  VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB&
                                           (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
                          FO_COLUMNWF_SURF (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
               ENDIF

               IF ( DO_SURFACE_LINEARIZATION ) THEN
                 VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
                   FO_SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
               ENDIF
            ENDIF

            IF ( DO_DNWELLING ) THEN
               VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,DNIDX) = &
                           FO_STOKES_ATMOS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,DNIDX)

               IF ( DO_COLUMN_LINEARIZATION ) THEN
                  VLIDORT_LinSup%SS%Col%TS_COLUMNWF_SS&
                                           (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,DNIDX) = &
                          FO_COLUMNWF_ATMOS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,DNIDX)
               ENDIF
            ENDIF

!  1/31/21. Version 2.8.3. Old code commented out (Copy FO results to VLIDORT local arrays)
!            STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)   = &
!              FO_STOKES_ATMOS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)
!            STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)     = &
!              FO_STOKES_SURF(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
!            IF ( DO_COLUMN_LINEARIZATION ) THEN
!              COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:) = &
!                FO_COLUMNWF_ATMOS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)
!              COLUMNWF_DB(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
!                FO_COLUMNWF_SURF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
!            ENDIF
!            IF ( DO_SURFACE_LINEARIZATION ) THEN
!              !SURFACEWF_DB(1:N_SURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
!              !  FO_SURFACEWF_DB(1:N_SURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
!              SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
!                FO_SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
!            ENDIF

          ENDIF

!  End user-angle and solar-source if blocks

        ENDIF
      !ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@ Addition of SLEAVE WFs to Corrected Directbeam @@@@@@@@@@
!    R. Spurr, 22 August 2012
!              IF ( DO_SURFACE_LEAVING .AND. DO_SLEAVE_WFS ) then
!                CALL VLIDORT_LSSL_DBCORRECTION ( &
!                  DO_SSCORR_OUTGOING, DO_UPWELLING, DO_SL_ISOTROPIC,       &
!                  DO_OBSERVATION_GEOMETRY,                                 &
!                  NSTOKES, NLAYERS, NBEAMS, N_SLEAVE_WFS, N_SURFACE_WFS,   &
!                  N_USER_VZANGLES, N_USER_RELAZMS, N_USER_LEVELS,           &
!                  PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 &
!                  PARTLAYERS_LAYERIDX, LEVELMASK_UP,                 &
!                  N_GEOMETRIES, VZA_OFFSETS, DO_REFLECTED_DIRECTBEAM,      &
!                  SS_FLUX_MULT, LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES, &
!                  T_DELT_USERM, T_UTUP_USERM, UP_LOSTRANS, UP_LOSTRANS_UT, &
!                  SURFACEWF_DB )
!              ENDIF
!@@@@@@@@@@ END Addition of SLEAVE Corrected Directbeam @@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  1/31/21. Version 2.8.3. DO_MSSTS option final installation
!    ==> MSSTS Auxiliary output from VFO is copied here. EITHER UPWELLING or DOWNWELLING
!    ==> [output for Multiple scatter Sphericity corrections, filled directly in Converge routines]

      IF ( DO_MSSTS ) THEN
         VLIDORT_Out%Main%TS_PATHGEOMS  (1,0:NLAYERS)                      = FO_THETA_ALL  (0:NLAYERS,1)
         VLIDORT_Out%Main%TS_PATHGEOMS  (2,0:NLAYERS)                      = FO_ALPHA      (0:NLAYERS,1)
         IF ( DO_UPWELLING ) THEN
            VLIDORT_Out%Main%TS_LOSTRANS   (1:N_SZANGLES,1:NLAYERS) = FO_LOSTRANS_UP(1:N_SZANGLES,1:NLAYERS)
            IF ( do_COLUMN_LINEARIZATION ) THEN
               DO Q = 1, N_TOTALCOLUMN_WFS
                  VLIDORT_LinOut%Col%TS_LC_LOSTRANS(Q,1:N_SZANGLES,1:NLAYERS) = FO_LC_LOSTRANS_UP(1:N_SZANGLES,1:NLAYERS,Q)
               ENDDO
            ENDIF
         ELSE IF ( DO_DNWELLING ) THEN
            VLIDORT_Out%Main%TS_LOSTRANS   (1:N_SZANGLES,1:NLAYERS) = FO_LOSTRANS_DN(1:N_SZANGLES,1:NLAYERS)
            IF ( do_COLUMN_LINEARIZATION ) THEN
               DO Q = 1, N_TOTALCOLUMN_WFS
                  VLIDORT_LinOut%Col%TS_LC_LOSTRANS(Q,1:N_SZANGLES,1:NLAYERS) = FO_LC_LOSTRANS_DN(1:N_SZANGLES,1:NLAYERS,Q)
               ENDDO
            ENDIF
         ENDIF
      ENDIF

!  ####################
!   MAIN FOURIER LOOP
!  ####################

!  Initialise Fourier loop
!  =======================

!  Set Number of Fourier terms (NMOMENTS = Maximum).
!    ( Starting from 0 = Fundamental )

      SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
      LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

!  No azimuth dependency for following cases
!   Other cases now disabled, 17 Janaury 2006

!      IF ( DO_TRANSMITTANCE_ONLY ) THEN
!        LOCAL_DO_NO_AZIMUTH = .TRUE.
!      ENDIF

!mick mod 3/30/2015 - consolidated if conditions
      IF ( .NOT.DO_SOLAR_SOURCES .OR. DO_MVOUT_ONLY ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.    
      ENDIF

!  set Fourier number (2 for Rayleigh only)

      IF ( LOCAL_DO_NO_AZIMUTH .OR. DO_FOCORR_ALONE ) THEN
        N_FOURIERS = 0
      ELSE
        IF ( DO_RAYLEIGH_ONLY  ) THEN
          N_FOURIERS = 2
        ELSE
          N_FOURIERS = NMOMENTS
        ENDIF
      ENDIF

!  re-set no-azimuth flag

      DO_NO_AZIMUTH = LOCAL_DO_NO_AZIMUTH

!  Initialise BVP telescoping. Important.

      DO_BVTEL_INITIAL = DO_BVP_TELESCOPING

!  Fourier loop
!  ============

!  Initialize

      LOCAL_ITERATION   = .TRUE.
      FOURIER = -1
      TESTCONV          = 0

!  set up solar beam flags. Required in all cases.
!   ---Even required for the thermal.....

      DO IBEAM = 1, NBEAMS
        BEAM_TESTCONV  ( IBEAM )  = 0
        BEAM_ITERATION ( IBEAM ) = .TRUE.
        DO L = 0, MAXFOURIER
          DO_MULTIBEAM   ( IBEAM, L ) = .TRUE.
        ENDDO
      ENDDO

!  Start Fourier loop
!  ------------------

      DO WHILE ( LOCAL_ITERATION .AND. FOURIER.LT.N_FOURIERS )

!  Fourier counter

        FOURIER = FOURIER + 1

!  Local start of user-defined streams. Should always be 1. No zenith tolerance now.
!    -- 1/31/21. Version 2.8.3. This has been generally disabled in favor of masking
!         LOCAL_UM_START = 1

!  4/28/19. Local media problem and planetary-problem flags
!    -- if the Planetary problem is set, then must have LOCAL_ALBTRN for BOA Unit illumination #2
!       (regardless of the input values of DO_ALBTRN_MEDIA)        
        
        LOCAL_DO_ALBTRN_MEDIA = .false.
        IF ( DO_ALBTRN_MEDIA(1) ) LOCAL_DO_ALBTRN_MEDIA(1) = ( FOURIER == 0 )
        IF ( DO_ALBTRN_MEDIA(2) ) LOCAL_DO_ALBTRN_MEDIA(2) = ( FOURIER == 0 )

        LOCAL_DO_PLANETARY_PROBLEM = .false.
        IF ( DO_PLANETARY_PROBLEM  ) then
           LOCAL_DO_PLANETARY_PROBLEM = ( FOURIER == 0 )
           LOCAL_DO_ALBTRN_MEDIA(2)   = ( FOURIER == 0 )
        ENDIF

!  Local copying of BRDF/SLEAVE inputs
!  -----------------------------------
  
!  1/31/21. Version 2.8.3. Copy Local BRDF Fourier-component Input (only what you need)

        IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
           BRDF_F_0(:,1:NSTREAMS,1:N_SZANGLES) = &
              VLIDORT_Sup%BRDF%TS_BRDF_F_0(FOURIER,:,1:NSTREAMS,1:N_SZANGLES)
           BRDF_F  (:,1:NSTREAMS,1:NSTREAMS)   = &
              VLIDORT_Sup%BRDF%TS_BRDF_F  (FOURIER,:,1:NSTREAMS,1:NSTREAMS)
           IF ( DO_USER_VZANGLES ) THEN
              USER_BRDF_F_0(:,1:N_USER_VZANGLES,1:N_SZANGLES) = &
                 VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0(FOURIER,:,1:N_USER_VZANGLES,1:N_SZANGLES)
              USER_BRDF_F(:,1:N_USER_VZANGLES,1:NSTREAMS)      = &
                 VLIDORT_Sup%BRDF%TS_USER_BRDF_F(FOURIER,:,1:N_USER_VZANGLES,1:NSTREAMS)
           ENDIF
        ENDIF

!  1/31/21. Version 2.8.3. Copy Local Linearized BRDF Fourier-component Input (only what you need)

        IF ( .NOT.DO_LAMBERTIAN_SURFACE .AND. DO_SURFACE_LINEARIZATION ) THEN
           LS_BRDF_F_0(1:N_SURFACE_WFS,:,1:NSTREAMS,1:N_SZANGLES) = &
                  VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0 (1:N_SURFACE_WFS,FOURIER,:,1:NSTREAMS,1:N_SZANGLES)
           LS_BRDF_F(1:N_SURFACE_WFS,:,1:NSTREAMS,1:NSTREAMS)     = &
                  VLIDORT_LinSup%BRDF%TS_LS_BRDF_F   (1:N_SURFACE_WFS,FOURIER,:,1:NSTREAMS,1:NSTREAMS)
           IF ( DO_USER_VZANGLES ) THEN
              LS_USER_BRDF_F_0(1:N_SURFACE_WFS,:,1:N_USER_VZANGLES,1:N_SZANGLES) = &
                  VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0(1:N_SURFACE_WFS,FOURIER,:,1:N_USER_VZANGLES,1:N_SZANGLES)
              LS_USER_BRDF_F(1:N_SURFACE_WFS,:,1:N_USER_VZANGLES,1:NSTREAMS)     = &
                  VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F  (1:N_SURFACE_WFS,FOURIER,:,1:N_USER_VZANGLES,1:NSTREAMS)
           ENDIF
        ENDIF

!  1/31/21. Version 2.8.3. Copy Local SLEAVE Fourier-component Input (only what you need)

        IF ( DO_SURFACE_LEAVING ) THEN
           SLTERM_F_0(1:NSTOKES,1:NSTREAMS,1:N_SZANGLES) = &
               VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0(FOURIER,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES)
           IF ( DO_USER_VZANGLES ) THEN
              USER_SLTERM_F_0(1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES) = &
                  VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(FOURIER,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES)
           ENDIF
        ENDIF

!  1/31/21. Version 2.8.3. Fourier SLEAVE copying now moved into Fourier loop.

        IF ( DO_SURFACE_LEAVING .AND. DO_SLEAVE_WFS ) THEN
           LSSL_SLTERM_F_0(1:N_SLEAVE_WFS,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES) = &
               VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0(1:N_SLEAVE_WFS,FOURIER,1:NSTOKES,1:NSTREAMS,1:N_SZANGLES)
           IF ( DO_USER_VZANGLES ) THEN
              LSSL_USER_SLTERM_F_0(1:N_SLEAVE_WFS,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES) = &
               VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0(1:N_SLEAVE_WFS,FOURIER,1:NSTOKES,1:N_USER_VZANGLES,1:N_SZANGLES)
           ENDIF
        ENDIF

!  Main call to VLidort Fourier module
!  -----------------------------------

!  7/22/16. Version 2.8 revision. Cleanup I/O listings, remove "Classical_Solution, Quad_Output"

!        write(*,*)' ..calculating fourier component',FOURIER

!   Version 2.8.1, Control for TOA/BOA isotropic illumination added, 3/23/19

!  4/9/19. Added Inputs, Water-leaving control, SOLARBEAM_BOATRANS and linearization.
!          Added Output, TRANS_ATMOS_FINAL amd LC Jacobian (for water-leaving self-consistency)    

!  4/28/19 Module for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!    -- introduced by R. Spurr 4/26/19. Controlled by flags LOCAL_DO_ALBTRN_MEDIA, LOCAL_DO_PLANETARY_PROBLEM
!    -- Associated outputs are TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES   

!  1/31/21. Version 2.8.3. Several additional changes to this argument list
!    -- Include flags DO_CLASSICAL_SOLUTION (Line 3), DO_MSSTS (line 8)
!    -- Include MSSTS outputs LAYER_MSSTS_F, SURF_MSSTS_F
!    -- Include ASYMTX Tolerance variable to the list (Line 8)
!    -- Use PPSTREAM masking system; replaces use of Observation/Doublet/Lattice
!    -- (RTS 2/16/21). Introduce Flag for using an NSTOKES = 2 calculation for Fourier 0 (set by hand)

!  5/24/21. Version 2.8.3. Add DO_INCLUDE_SLEAVEWFS as an output

!mick debug
!write(*,*)
!write(*,*) 'entering VLIDORT_LCS_FOURIER'

!mick debug
!write(*,*) 'at 15'
!call system_mem_usage3(valueRSS)

        CALL VLIDORT_LCS_FOURIER ( FOURIER, DO_FOURIER0_NSTOKES2, &
            DO_UPWELLING, DO_DNWELLING, DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY,                  & !Input flags (RT operation)
            DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,           & !Input flags (RT operation)
            DO_PLANE_PARALLEL, DO_LAYER_SCATTERING, DO_SOLUTION_SAVING, DO_BVTEL_INITIAL,           & !Input flags (performance)
            DO_CLASSICAL_SOLUTION,                                                                  & !Input flags
            DO_MSMODE_VLIDORT, DO_MULTIBEAM, LOCAL_DO_ALBTRN_MEDIA, LOCAL_DO_PLANETARY_PROBLEM,     & !Input flags (Beam/Planetary)
            DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,                             & !Input flags (Surface)
            DO_MSMODE_THERMAL, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_THERMAL_TRANSONLY,      & !Input flags (thermal)
            DO_FOCORR_ALONE, DO_DBCORRECTION, DO_TOA_CONTRIBS, DO_PARTLAYERS, DO_REAL_EIGENSOLVER,  & !Input Bookkeeping
            DO_MSSTS, DO_TOAFLUX, TOAFLUX, DO_BOAFLUX, BOAFLUX, ASYMTX_TOLERANCE, DO_DEBUG_WRITE,   & !Input MSSTS & TOA/BOA fluxes
            DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION, TF_MAXITER, TF_CRITERION,        & !Input Water-leaving control
            DO_ATMOS_LINEARIZATION, DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION,              & !Input Linearization flags
            DO_SIMULATION_ONLY, DO_SLEAVE_WFS, DO_ATMOS_LBBF, DO_SURFACE_LBBF,                      & !Input Linearization flags
            NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_VZANGLES, N_USER_LEVELS, N_THERMAL_COEFFS,   & !Input Basic integers
            NMOMENTS, NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,       & !Input Derived integers
            N_PPSTREAMS, PPSTREAM_MASK, N_PARTLAYERS, N_DIRECTIONS, TAYLOR_ORDER,                   & !Input Derived integers
            N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,                     & !In Num(Lin)
            N_ALLLAYERS_UP, N_ALLLAYERS_DN, FLUX_FACTOR, FLUXVEC, COS_SZANGLES, SZA_LOCAL_INPUT, SUN_SZA_COSINES,   & !In SZAs/FLUX
            USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,                     & !In Streams
            MUELLER_INDEX, DMAT, BVP_REGULAR_FLAG, WHICH_DIRECTIONS, SOLARBEAM_BOATRANS,                            & !In Bkkeep
            LEVELMASK_UP, LEVELMASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,               & !In Bkkeep
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, LC_SOLARBEAM_BOATRANS, LC_SOLARBEAM_ATRANS,              & !In Bkkeep
            Misc, Therm, Mult, LAC_Misc, L_Therm, LC_Mult,                                                          & !In Packing
            ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0,                                                   & !In ALB/BRDF
            SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0, SURFBB, EMISSIVITY, USER_EMISSIVITY,                     & !In SL/EMISS
            LS_BRDF_F_0, LS_BRDF_F,  LS_USER_BRDF_F_0, LS_USER_BRDF_F, LS_EMISSIVITY, LS_USER_EMISSIVITY,           & !In LS
            LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,                                           & !In LSSL
            PIMM_11, PIMM_KM, BVTEL_FOURIER, DO_INCLUDE_THERMEMISS, DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS,       & !Out Bkkeep
            STOKES_F, MEANST_DIFFUSE, FLUX_DIFFUSE, DNMEANST_DIRECT, DNFLUX_DIRECT, MS_CONTRIBS_F,                  & !Out STOKES
            COLUMNWF_F,  MEANST_DIFFUSE_COLWF,  FLUX_DIFFUSE_COLWF,  DNMEANST_DIRECT_COLWF, DNFLUX_DIRECT_COLWF,    & !Out COL JACS
            SURFACEWF_F, MEANST_DIFFUSE_SURFWF, FLUX_DIFFUSE_SURFWF, DO_INCLUDE_SLEAVEWFS,                          & !Out SURF JACS
            ABBWFS_JACOBIANS, ABBWFS_FLUXES, SBBWFS_JACOBIANS, SBBWFS_FLUXES,                                       & !Out LBBF JACS
            TRANS_ATMOS_FINAL, TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,                   & !Out 4/26/19
            LC_TRANS_ATMOS_FINAL, LC_TRANSBEAM, LC_ALBMED_USER, LC_ALBMED_FLUXES, LC_TRNMED_USER, LC_TRNMED_FLUXES, & !Out 4/26/19
            LAYER_MSSTS_F, SURF_MSSTS_F, LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F,      & !Output SPECIAL
            STATUS_SUB, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )                                                          !Out Status

!  error handling

        IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
           TRACE_3 = ' Called by VLIDORT_LCS_MASTER '
           STATUS_CALCULATION = VLIDORT_SERIOUS
           VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
           VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
           VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
           VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
           VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3
           RETURN
        ENDIF

!  4/9/19. FO-DB adjustment for Fourier 0, using self-consistent water-leaving term

!  1/31/21. Version 2.8.3. Use Type structure variables directly
!      ==> Replace arrays STOKES_DB with VLIDORT_Sup%SS%TS_STOKES_DB. Similarly, linearization
!  1/31/21. Version 2.8.3. Add Doublet geometry option. Add FOCORR_EXTERNAL to if clause

!  5/19/21. Version 2.8.3. Mick Fix - added multiplicative factor FO_CUMTRANS(UTA,G) in defining LS_CUMSOURCE_DB
!                     in SLEAVE WF portion of each geometry section below

!  5/28/21. Version 2.8.3. Careful with surface wf indexing.....use Q1 instead of Q

        if ( FOURIER .eq. 0 .and. DO_FOCORR .and. .NOT.DO_FOCORR_EXTERNAL .and. DO_WATER_LEAVING ) then

!  Obsgeom

           O1 = 1
           IF ( DO_OBSERVATION_GEOMETRY ) THEN
              DO IB = 1, NBEAMS

                 SLTERM_LOCAL = FO_SLTERM(O1,IB) * TRANS_ATMOS_FINAL(IB)
                 DO UTA = 1, N_USER_LEVELS
                     CUMSOURCE_DB = FO_CUMTRANS(UTA,IB) * SLTERM_LOCAL
                     VLIDORT_Sup%SS%TS_STOKES_DB(UTA,IB,O1) = VLIDORT_Sup%SS%TS_STOKES_DB(UTA,IB,O1) + CUMSOURCE_DB
                 ENDDO

                 IF ( DO_COLUMN_LINEARIZATION ) THEN
                    NWFS = N_TOTALCOLUMN_WFS
                    LC_SLTERM_LOCAL(1:NWFS) = FO_SLTERM(O1,IB) * LC_TRANS_ATMOS_FINAL(IB,1:NWFS)
                    DO UTA = 1, N_USER_LEVELS
                       DO Q = 1, NWFS
                          L_CUMSOURCE_DB =    FO_CUMTRANS(UTA,IB)   * LC_SLTERM_LOCAL(Q) &
                                         + FO_LC_CUMTRANS(UTA,IB,Q) *    SLTERM_LOCAL
                          VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB(Q,UTA,IB,O1) = &
                             VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB(Q,UTA,IB,O1) + L_CUMSOURCE_DB
                       ENDDO
                    ENDDO
                 ENDIF

                 IF ( DO_WATER_LEAVING .AND. DO_SLEAVE_WFS ) THEN
                    NWFS = N_SLEAVE_WFS
                    LS_SLTERM_LOCAL(1:NWFS) = FO_LSSL_SLTERM(O1,IB,1:NWFS) * TRANS_ATMOS_FINAL(IB)
                    DO UTA = 1, N_USER_LEVELS
                       DO Q = 1, NWFS
                          LS_CUMSOURCE_DB = FO_CUMTRANS(UTA,IB) * LS_SLTERM_LOCAL(Q) ; Q1 = Q + N_SURFACE_WFS
                          VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,IB,O1) = &
                             VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,IB,O1) + LS_CUMSOURCE_DB
                       ENDDO
                    ENDDO
                 ENDIF
              ENDDO

!  Doublet

           ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN

!  5/24/21. Version 2.8.3. Rob Fix - Doublet
!       -- SLTERM_LOCAL    must go inside the geometry loop, and FO_SLTERM must use index g (not g1)
!       -- LC_SLTERM_LOCAL must go inside the geometry loop, and FO_SLTERM must use index g (not g1)
!       -- LS_SLTERM_LOCAL must go inside the geometry loop, and FO_LSSL_SLTERM must use index g (not g1)

              DO IB = 1, NBEAMS

                 G1 = SZD_OFFSETS(IB) + 1 ; G2 = SZD_OFFSETS(IB) + N_USER_VZANGLES
!                 SLTERM_LOCAL = FO_SLTERM(O1,G1) * TRANS_ATMOS_FINAL(IB)
                 DO UTA = 1, N_USER_LEVELS ; DO G = G1, G2
                    SLTERM_LOCAL = FO_SLTERM(O1,G) * TRANS_ATMOS_FINAL(IB)
                    CUMSOURCE_DB = FO_CUMTRANS(UTA,G) * SLTERM_LOCAL
                    VLIDORT_Sup%SS%TS_STOKES_DB(UTA,G,O1) = VLIDORT_Sup%SS%TS_STOKES_DB(UTA,G,O1) + CUMSOURCE_DB
                 ENDDO ; ENDDO

                 IF ( DO_COLUMN_LINEARIZATION ) THEN
                    NWFS = N_TOTALCOLUMN_WFS
!                    LC_SLTERM_LOCAL(1:NWFS) = FO_SLTERM(O1,G1) * LC_TRANS_ATMOS_FINAL(IB,1:NWFS)
                    DO UTA = 1, N_USER_LEVELS ; DO G = G1, G2
                       LC_SLTERM_LOCAL(1:NWFS) = FO_SLTERM(O1,G) * LC_TRANS_ATMOS_FINAL(IB,1:NWFS)
                       DO Q = 1, NWFS
                          L_CUMSOURCE_DB =    FO_CUMTRANS(UTA,G)   * LC_SLTERM_LOCAL(Q) &
                                         + FO_LC_CUMTRANS(UTA,G,Q) *    SLTERM_LOCAL
                          VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB(Q,UTA,G,O1) = &
                             VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB(Q,UTA,G,O1) + L_CUMSOURCE_DB
                       ENDDO
                    ENDDO ; ENDDO
                 ENDIF

                 IF ( DO_WATER_LEAVING .AND. DO_SLEAVE_WFS ) THEN
                    NWFS = N_SLEAVE_WFS
!                    LS_SLTERM_LOCAL(1:NWFS) = FO_LSSL_SLTERM(O1,G1,1:NWFS) * TRANS_ATMOS_FINAL(IB)
                    DO UTA = 1, N_USER_LEVELS ; DO G = G1, G2
                       LS_SLTERM_LOCAL(1:NWFS) = FO_LSSL_SLTERM(O1,G,1:NWFS) * TRANS_ATMOS_FINAL(IB)
                       DO Q = 1, NWFS
                          LS_CUMSOURCE_DB =  LS_SLTERM_LOCAL(Q) ; Q1 = Q + N_SURFACE_WFS
                          VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,G,O1) = &
                             VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,G,O1) + LS_CUMSOURCE_DB
                       ENDDO
                    ENDDO ; ENDDO
                 ENDIF
 
              ENDDO

!  Lattice

           ELSE

!  5/24/21. Version 2.8.3. Rob Fix - Lattice
!       -- SLTERM_LOCAL    must go inside the geometry loop, and FO_SLTERM must use index g (not g1)
!       -- LC_SLTERM_LOCAL must go inside the geometry loop, and FO_SLTERM must use index g (not g1)
!       -- LS_SLTERM_LOCAL must go inside the geometry loop, and FO_LSSL_SLTERM must use index g (not g1)

              DO IB = 1, NBEAMS ; DO UM = 1, N_USER_VZANGLES
                 G1 = VZA_OFFSETS(IB,UM) + 1 ; G2 = VZA_OFFSETS(IB,UM) + N_USER_RELAZMS
!                 SLTERM_LOCAL = FO_SLTERM(O1,G1) * TRANS_ATMOS_FINAL(IB)
                 DO UTA = 1, N_USER_LEVELS ; DO G = G1, G2
                    SLTERM_LOCAL = FO_SLTERM(O1,G) * TRANS_ATMOS_FINAL(IB)
                    CUMSOURCE_DB = FO_CUMTRANS(UTA,G) * SLTERM_LOCAL
                    VLIDORT_Sup%SS%TS_STOKES_DB(UTA,G,O1) = VLIDORT_Sup%SS%TS_STOKES_DB(UTA,G,O1) + CUMSOURCE_DB
                 ENDDO ; ENDDO

                 IF ( DO_COLUMN_LINEARIZATION ) THEN
                    NWFS = N_TOTALCOLUMN_WFS
!                   LC_SLTERM_LOCAL(1:NWFS) = FO_SLTERM(O1,G1) * LC_TRANS_ATMOS_FINAL(IB,1:NWFS)
                    DO UTA = 1, N_USER_LEVELS ; DO G = G1, G2
                       LC_SLTERM_LOCAL(1:NWFS) = FO_SLTERM(O1,G) * LC_TRANS_ATMOS_FINAL(IB,1:NWFS)
                       DO Q = 1, NWFS
                          L_CUMSOURCE_DB =    FO_CUMTRANS(UTA,G)   * LC_SLTERM_LOCAL(Q) &
                                         + FO_LC_CUMTRANS(UTA,G,Q) *    SLTERM_LOCAL
                          VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB(Q,UTA,G,O1) = &
                             VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB(Q,UTA,G,O1) + L_CUMSOURCE_DB
                       ENDDO
                    ENDDO ; ENDDO
                 ENDIF

                 IF ( DO_WATER_LEAVING .AND. DO_SLEAVE_WFS ) THEN
                    NWFS = N_SLEAVE_WFS
!                    LS_SLTERM_LOCAL(1:NWFS) = FO_LSSL_SLTERM(O1,G1,1:NWFS) * TRANS_ATMOS_FINAL(IB)
                    DO UTA = 1, N_USER_LEVELS ; DO G = G1, G2
                       LS_SLTERM_LOCAL(1:NWFS) = FO_LSSL_SLTERM(O1,G,1:NWFS) * TRANS_ATMOS_FINAL(IB)
                       DO Q = 1, NWFS
                          LS_CUMSOURCE_DB = FO_CUMTRANS(UTA,G) * LS_SLTERM_LOCAL(Q) ; Q1 = Q + N_SURFACE_WFS
                          VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,G,O1) = &
                            VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,G,O1) + LS_CUMSOURCE_DB
                       ENDDO
                    ENDDO ; ENDDO
                 ENDIF
              ENDDO ; ENDDO
           ENDIF

!  END WATER-LEAVING CORRECTION

        ENDIF

!  Output for WLADJUSTED Water-Leaving . Introduced 4/22/19 for Version 2.8.1
      ! Sleave Results need to be modified from their original inputs.
      ! Debug results - USE WITH CAUTION. Note the preliminary zeroing to avoid unassigned arrays.
!  -- 1/31/21. Version 2.8.3, SLTERM_F arrays defined locally for each Fourier

        IF ( FOURIER.EQ.0 .AND. DO_WATER_LEAVING .AND. DO_WLADJUSTED_OUTPUT ) THEN
           VLIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC = ZERO
           VLIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT    = ZERO
           VLIDORT_Out%WLOut%TS_WLADJUSTED_F_Ords_0  = ZERO
           VLIDORT_Out%WLOut%TS_WLADJUSTED_F_User_0  = ZERO
           DO IB = 1, NBEAMS
              TFACTOR = TRANS_ATMOS_FINAL(ib)
              VLIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC(O1,IB)                   = TFACTOR * SLTERM_ISOTROPIC(O1,IB)
              VLIDORT_Out%WLOut%TS_WLADJUSTED_F_Ords_0(FOURIER,O1,1:NSTREAMS,IB) = TFACTOR * SLTERM_F_0(O1,1:NSTREAMS,IB)
              IF ( DO_USER_VZANGLES ) THEN
                 VLIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT(O1,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IB) = &
                                         TFACTOR * SLTERM_USERANGLES(O1,1:N_USER_VZANGLES,1:N_USER_RELAZMS,IB)
                 VLIDORT_Out%WLOut%TS_WLADJUSTED_F_User_0(FOURIER,O1,1:N_USER_VZANGLES,IB) = &
                                         TFACTOR * USER_SLTERM_F_0(O1,1:N_USER_VZANGLES,IB)
              ENDIF
           ENDDO
        ENDIF

! #################### INNOVATIONS 5/5/20 ##############
!  4/29/20. Special Fourier-component output for Rayleigh + Planetary-problem TOA Upwelling situations
!    ==> Any number of profile wfs allowed for Fourier 0, 1, 2
!    ==> Only 1 SURFACEWF w.r.t ALBEDO     for Fourier 0
!         IF ( SPECIAL_RAYF_OUTPUT ) then
!           UTA = 1 ; M =  FOURIER
!           IF ( DO_OBSERVATION_GEOMETRY ) THEN
!             DO IB =1, NBEAMS ; DO O1 = 1, NSTOKES ; LUM = 1
!               VLIDORT_Out%Main%TS_TOAUP_RAYSTOKES_FOURIER(LUM,IB,O1,M) =  STOKES_F(UTA,LUM,IB,O1,UPIDX)
!             ENDDO ; ENDDO
!             if ( do_column_linearization ) then
!               NWFS = N_TOTALCOLUMN_WFS
!               DO IB =1, NBEAMS ; DO O1 = 1, NSTOKES ; LUM = 1
!                 VLIDORT_LinOut%Col%TS_TOAUP_RAYCOLWF_FOURIER(1:NWFS,LUM,IB,O1,M) = COLUMNWF_F(1:NWFS,UTA,LUM,IB,O1,UPIDX)
!               ENDDO ; ENDDO
!             ENDIF
!             if ( do_surface_linearization .AND. Fourier.EQ.0  ) then
!               DO IB =1, NBEAMS ; DO O1 = 1, NSTOKES ; LUM = 1
!                 VLIDORT_LinOut%Surf%TS_TOAUP_RAYSURFWF_FOURIER(LUM,IB,O1) = SURFACEWF_F(1,UTA,LUM,IB,O1,UPIDX)
!               ENDDO ; ENDDO
!             ENDIF
!           ELSE
!             DO IB =1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
!               VLIDORT_Out%Main%TS_TOAUP_RAYSTOKES_FOURIER(UM,IB,O1,M) =  STOKES_F(UTA,UM,IB,O1,UPIDX)
!             ENDDO ; ENDDO ; ENDDO
!             if ( do_column_linearization ) then
!               NWFS = N_TOTALCOLUMN_WFS
!               DO IB =1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
!                 VLIDORT_LinOut%Col%TS_TOAUP_RAYCOLWF_FOURIER(1:NWFS,UM,IB,O1,M) = COLUMNWF_F(1:NWFS,UTA,UM,IB,O1,UPIDX)
!               ENDDO ; ENDDO ; ENDDO
!             ENDIF
!             if ( do_surface_linearization .AND. Fourier.EQ.0  ) then
!               DO IB =1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
!                 VLIDORT_LinOut%Surf%TS_TOAUP_RAYSURFWF_FOURIER(UM,IB,O1) = SURFACEWF_F(1,UTA,UM,IB,O1,UPIDX)
!               ENDDO ; ENDDO ; ENDDO
!             ENDIF
!           ENDIF
!         ENDIF
! #################### INNOVATIONS 5/5/20 ##############

!  Fourier summation and Convergence examination
!  ---------------------------------------------

!mick fix 3/30/2015 - added if condition and moved azimuth block from
!                     before call to VLIDORT_LCS_FOURIER to here

!  Begin convergence if block

        IF ( .NOT.DO_MVOUT_ONLY ) THEN

!  Azimuth cosine/sine factors
!    - Use of adjusted geometries retired for Version 2.8
!mick mod 1/5/2021 - added rewritten section with doublet from "vlidort_masters"
!mick fix 1/5/2021 - added DO_DOUBLET_GEOMETRY to 1st IF condition

          IF ( FOURIER .GT. 0 ) THEN
            DFC = DBLE(FOURIER)
            IF ( .NOT.DO_OBSERVATION_GEOMETRY .AND. .NOT.DO_DOUBLET_GEOMETRY ) THEN
              DO UA = 1, LOCAL_N_USERAZM
                AZM_ARGUMENT = USER_RELAZMS(UA) * DFC * DEG_TO_RAD
                AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,1)  = COS(AZM_ARGUMENT)
                AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,3)  = SIN(AZM_ARGUMENT)
                AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,2)  = AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,1)
                AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,4)  = AZMFAC(1:N_USER_VZANGLES,1:NBEAMS,UA,3)
              ENDDO
            ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
              DO UM = 1, N_USER_VZANGLES
                AZM_ARGUMENT = USER_RELAZMS(UM) * DFC * DEG_TO_RAD
                AZMFAC(UM,1:N_SZANGLES,LUA,1)  = COS(AZM_ARGUMENT)
                AZMFAC(UM,1:N_SZANGLES,LUA,3)  = SIN(AZM_ARGUMENT)
                AZMFAC(UM,1:N_SZANGLES,LUA,2)  = AZMFAC(UM,1:N_SZANGLES,LUA,1)
                AZMFAC(UM,1:N_SZANGLES,LUA,4)  = AZMFAC(UM,1:N_SZANGLES,LUA,3)
              ENDDO
            ELSE
              DO IB = 1, NSOURCES
                AZM_ARGUMENT = USER_RELAZMS(IB) * DFC * DEG_TO_RAD
                AZMFAC(LUM,IB,LUA,1)  = COS(AZM_ARGUMENT)
                AZMFAC(LUM,IB,LUA,3)  = SIN(AZM_ARGUMENT)
                AZMFAC(LUM,IB,LUA,2)  = AZMFAC(LUM,IB,LUA,1)
                AZMFAC(LUM,IB,LUA,4)  = AZMFAC(LUM,IB,LUA,3)
              ENDDO
            ENDIF
          ENDIF

!   -- only done for beams which are still not converged.
!      This is controlled by flag DO_MULTIBEAM

!   -- new criterion, SS is added for Fourier = 0, as this means that
!      higher-order terms will be relatively smaller, which implies
!      faster convergence in some circumstances (generally not often).

          IBEAM_COUNT = 0
          DO IBEAM = 1, NBEAMS
            IF ( DO_MULTIBEAM ( IBEAM, FOURIER ) ) THEN

!  Convergence and radiance summation. Version 2.8, Clean up I/O listings, 7/8/16

!  Fourier summation of linearization quantities
!     mick fix 9/6/2012 - added N_SLEAVE_WFS to call

!  1/31/21. Version 2.8.3. Several changes, including use of MSST output.
!    -- Convergence subroutines now have their own module.
!    -- Add Doublet geometry option (new convergence routine). Version 2.8.2 Feature
!    -- Add DO_MSSTS (input) and LAYER_MSSTS_F, SURF_MSSTS_F (outputs) for the OBSGEO routine
!    -- Use/Fill Regular    Type structure variables directly (replaces STOKES, STOKES_SS, STOKES_DB)
!    -- Use/Fill linearized type structure variables directly (replaces COLUMNWF, SURFACEWF, associated FO inputs)
!    -- Drop LOCAL_UM_START from the Lattice converge routines.

              IF ( DO_OBSERVATION_GEOMETRY ) THEN

!  Fourier summation of regular quantities

                CALL VLIDORT_CONVERGE_OBSGEO ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_DBCORRECTION, DO_UPWELLING, DO_NO_AZIMUTH,     & ! Input flags
                  DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST, DO_TOA_CONTRIBS,        & ! Input flags
                  DO_MSSTS, NSTOKES, NSTREAMS, NLAYERS, N_USER_LEVELS, IBEAM, FOURIER,          & ! Input numbers
                  N_CONVTESTS, VLIDORT_ACCURACY, AZMFAC, N_DIRECTIONS, WHICH_DIRECTIONS,        & ! Input Bookkeep, Conv.
                  STOKES_F, MS_CONTRIBS_F, LAYER_MSSTS_F, SURF_MSSTS_F, VLIDORT_Sup%SS,         & ! Input/Output fields
                  VLIDORT_Out%Main, FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )  ! Output diagnostics

!  Fourier summation of linearization quantities
!   -- 5/24/21. Version 2.8.3. Add DO_INCLUDE_SLEAVEWFS (for surface-leaving weighting functions)

                IF ( DO_LINEARIZATION ) THEN
                  CALL VLIDORT_LCS_CONVERGE_OBSGEO ( &
                    DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH, DO_DBCORRECTION,  & ! Input flags
                    DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_LAMBERTIAN_SURFACE,  & ! Input flags
                    DO_MSSTS, IBEAM, FOURIER, NSTOKES, NLAYERS, N_USER_LEVELS, N_DIRECTIONS,   & ! Input numbers/indices
                    N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS, DO_INCLUDE_SLEAVEWFS,      & ! Input linearization control
                    AZMFAC, WHICH_DIRECTIONS, COLUMNWF_F, SURFACEWF_F,                         & ! Input Bookkeeping/Fourier Jacs 
                    LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F,      & ! Input Jacobian MSST inputs
                    VLIDORT_LinSup%SS, VLIDORT_LinOut )                                          ! Input/Output fields
                ENDIF

              ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN

!  Fourier summation of regular quantities

                CALL VLIDORT_CONVERGE_DOUBLET ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_DBCORRECTION, DO_UPWELLING,                          & ! Input flags
                  DO_NO_AZIMUTH, DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST,                & ! Input flags
                  NSTOKES, NSTREAMS, N_OUT_STREAMS, N_USER_LEVELS, IBEAM, FOURIER,                    & ! Input numbers
                  N_CONVTESTS, VLIDORT_ACCURACY, SZD_OFFSETS, AZMFAC, N_DIRECTIONS, WHICH_DIRECTIONS, & ! Input Bookkeep, Conv.
                  STOKES_F, VLIDORT_Sup%SS, VLIDORT_Out%Main,                                         & ! Input/Output fields
                  FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM)  )                         ! Output Convergence

!  Fourier summation of linearization quantities
!   -- 5/24/21. Version 2.8.3. Add DO_INCLUDE_SLEAVEWFS (for surface-leaving weighting functions)

                IF ( DO_LINEARIZATION ) THEN
                  CALL VLIDORT_LCS_CONVERGE_DOUBLET ( &
                    DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH, DO_DBCORRECTION, & ! Input flags
                    DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_LAMBERTIAN_SURFACE, & ! Input flags
                    IBEAM, FOURIER, NSTOKES, N_OUT_STREAMS, N_USER_LEVELS, N_DIRECTIONS,      & ! Input control numbers
                    N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS, DO_INCLUDE_SLEAVEWFS,     & ! Input linearization control
                    SZD_OFFSETS, WHICH_DIRECTIONS, AZMFAC,                                    & ! Input bookkeeping
                    COLUMNWF_F, SURFACEWF_F, VLIDORT_LinSup%SS, VLIDORT_LinOut )                ! Input Azm, fields
                ENDIF

              ELSE  !LATTICE

!  Fourier summation of regular quantities

                CALL VLIDORT_CONVERGE ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_DBCORRECTION, DO_UPWELLING, DO_NO_AZIMUTH, & ! Input flags
                  DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST, DO_TOA_CONTRIBS,    & ! Input flags
                  NSTOKES, NSTREAMS, NLAYERS, N_OUT_STREAMS, N_USER_RELAZMS, N_USER_LEVELS, & ! Input control numbers
                  IBEAM, FOURIER, N_CONVTESTS, VLIDORT_ACCURACY, VZA_OFFSETS, AZMFAC,       & ! Input numbers, convergence 
                  N_DIRECTIONS, WHICH_DIRECTIONS, LOCAL_N_USERAZM,                          & ! Input bookkeeping
                  STOKES_F, MS_CONTRIBS_F, VLIDORT_Sup%SS, VLIDORT_Out%Main,                & ! Input and output fields
                  FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )                ! Output diagnostics

!  Fourier summation of linearization quantities
!   -- 5/24/21. Version 2.8.3. Add DO_INCLUDE_SLEAVEWFS (for surface-leaving weighting functions)

                IF ( DO_LINEARIZATION ) THEN
                  CALL VLIDORT_LCS_CONVERGE ( &
                    DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH, DO_DBCORRECTION, & ! Input flags
                    DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_LAMBERTIAN_SURFACE, & ! Input flags
                    IBEAM, FOURIER, NSTOKES, N_OUT_STREAMS, N_USER_LEVELS, N_DIRECTIONS,      & ! Input control numbers
                    N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS, DO_INCLUDE_SLEAVEWFS,     & ! Input linearization control
                    VZA_OFFSETS, WHICH_DIRECTIONS, LOCAL_N_USERAZM, AZMFAC,                   & ! Input bookkeeping
                    COLUMNWF_F, SURFACEWF_F, VLIDORT_LinSup%SS, VLIDORT_LinOut )                ! Input/Output fields
                ENDIF

              ENDIF

!  Check number of beams already converged

              IF ( BEAM_ITERATION(IBEAM) ) THEN
                IBEAM_COUNT = IBEAM_COUNT + 1
              ELSE
                DO L = FOURIER+1,MAXFOURIER
                  DO_MULTIBEAM (IBEAM,L) = .FALSE.
                ENDDO
              ENDIF

!  end beam count loop

            ENDIF
          ENDDO

!  If all beams have converged, stop iteration

          IF ( IBEAM_COUNT .EQ. 0 ) LOCAL_ITERATION = .FALSE.

!  End convergence if block

        END IF

!  Fourier output
!  --------------

!  Open file if Fourier = 0
!  Write Standard   Fourier output (radiances)
!  Write additional Fourier output (linearizations)
!  Close file if iteration has finished
!  New comment:
!    If the SS correction is set, Fourier=0 will include SS field
!mick mod 1/5/2021 - re-introduced updated subroutine VLIDORT_WRITEFOURIER

        IF ( DO_WRITE_FOURIER ) THEN
          FUNIT = VLIDORT_FUNIT
          IF ( FOURIER .EQ. 0 ) OPEN(FUNIT, FILE = FOURIER_WRITE_FILENAME, STATUS='REPLACE')
          CALL VLIDORT_WRITEFOURIER ( FUNIT, FOURIER, N_USER_LEVELS, N_OUT_STREAMS, LOCAL_N_USERAZM, &
                                      NBEAMS, NSTOKES, N_DIRECTIONS, WHICH_DIRECTIONS, STOKES_F)
          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
        ENDIF

!  1/31/21. Version 2.8.3.  Copying for Sleave Fourier components must be done here
!    -- Additional Control for Externalized input (SLEAVE). Introduced 4/22/19 for Version 2.8.1
!    -- Sleave Results May have been modified from their original inputs.
!    -- Allows you to come out with modified SLEAVE, so you can use it!!!!!!!!!!!!

        IF ( DO_EXTERNAL_WLEAVE ) THEN
           O1 = 1
           IF ( FOURIER .EQ. 0 ) THEN
              VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(O1,1:N_SZANGLES) = SLTERM_ISOTROPIC(O1,1:N_SZANGLES)
              IF ( DO_USER_VZANGLES ) THEN
                 VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(O1,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES) = &
                                       SLTERM_USERANGLES(O1,1:N_USER_VZANGLES,1:N_USER_RELAZMS,1:N_SZANGLES)
              ENDIF
           ENDIF
           VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0(FOURIER,O1,1:NSTREAMS,1:N_SZANGLES) = SLTERM_F_0(O1,1:NSTREAMS,1:N_SZANGLES)
           VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(FOURIER,O1,1:N_USER_VZANGLES,1:N_SZANGLES) = &
                                 USER_SLTERM_F_0(O1,1:N_USER_VZANGLES,1:N_SZANGLES)
        ENDIF

!  End Fourier loop

      ENDDO

!  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  Geophysical input (scenario) write
!  ----------------------------------

      IF ( DO_WRITE_SCENARIO ) THEN
        SUNIT = VLIDORT_SCENUNIT
        OPEN(SUNIT,FILE=SCENARIO_WRITE_FILENAME,STATUS='REPLACE')

        CALL VLIDORT_WRITESCEN ( &
          SUNIT, DO_DELTAM_SCALING, NSTREAMS, NLAYERS, NGREEK_MOMENTS_INPUT,      &
          N_SZANGLES, SZANGLES, N_USER_RELAZMS, USER_RELAZMS,                     &
          N_USER_VZANGLES, USER_VZANGLES, N_USER_LEVELS, USER_LEVELS,             &
          OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, DO_LAMBERTIAN_SURFACE, ALBEDO, &
          DO_THERMAL_EMISSION, N_THERMAL_COEFFS,  DO_SURFACE_EMISSION, SURFBB,    &
          QUAD_STREAMS, QUAD_WEIGHTS, QUAD_ANGLES, DO_NO_AZIMUTH,                 &
          NMOMENTS, GREEKMAT_INDEX, TAUGRID_INPUT, OMEGA_TOTAL, GREEKMAT_TOTAL, TAUGRID )

        IF ( DO_LINEARIZATION ) THEN
          CALL VLIDORT_L_WRITESCEN ( &
            SUNIT, NLAYERS, DO_ATMOS_LINEARIZATION, &
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
            L_OMEGA_TOTAL_INPUT, L_DELTAU_VERT_INPUT )
        ENDIF

        CLOSE(SUNIT)
      ENDIF

!  ========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (IN/OUT variables)
!  ========================================================

!  FOCORR Booleans reorganized for Version 2.8
!mick mod 9/19/2017 - reordered FO variables to conform to newly modified input type structure
!                   - DO_FOCORR_ALONE now defined internally

!  copy FO flags

      VLIDORT_ModIn%MBool%TS_DO_FOCORR              = DO_FOCORR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL     = DO_FOCORR_EXTERNAL
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR        = DO_FOCORR_NADIR
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING     = DO_FOCORR_OUTGOING
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT      = DO_SSCORR_USEFMAT

!  1/31/21. Version 2.8.3. Copy two additional flags (DO_DOUBLET, CLASSICAL)

      VLIDORT_ModIn%MBool%TS_DO_CLASSICAL_SOLUTION  = DO_CLASSICAL_SOLUTION
      VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY    = DO_DOUBLET_GEOMETRY

!  1/31/21. Version 2.8.3. Drop the DO_SSCORR_TRUNCATION
!      VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION   = DO_SSCORR_TRUNCATION

!  Additional Control for Externalized input (SLEAVE). Introduced 3/18/19 for Version 2.8.1

      VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE     = DO_EXTERNAL_WLEAVE

!  Solar control

      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = DO_SOLAR_SOURCES
      VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = DO_REFRACTIVE_GEOMETRY
      VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = DO_CHAPMAN_FUNCTION

      !VLIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY      = DO_ISOTROPIC_ONLY
      !VLIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH          = DO_NO_AZIMUTH
      !VLIDORT_ModIn%MBool%TS_DO_ALL_FOURIER         = DO_ALL_FOURIER

!  Performance control

      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY       = DO_RAYLEIGH_ONLY
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING      = DO_DELTAM_SCALING
      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING     = DO_SOLUTION_SAVING
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING     = DO_BVP_TELESCOPING

!  RT Model control

      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES       = DO_USER_VZANGLES
      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = DO_ADDITIONAL_MVOUT
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = DO_MVOUT_ONLY
      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = DO_THERMAL_TRANSONLY
      VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY= DO_OBSERVATION_GEOMETRY

!  Modified control inputs

      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT   = NGREEK_MOMENTS_INPUT

!  Modified beam inputs

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES             = N_SZANGLES
      VLIDORT_ModIn%MSunRays%TS_SZANGLES(1:N_SZANGLES) = SZANGLES(1:N_SZANGLES)

!  Modified user value inputs

      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS         = N_USER_RELAZMS
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS) = USER_RELAZMS(1:N_USER_RELAZMS)

!mick fix 9/19/2017 - added IF condition
      IF ( DO_USER_VZANGLES ) THEN
        VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES      = N_USER_VZANGLES
        VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:N_USER_VZANGLES) = USER_VZANGLES(1:N_USER_VZANGLES)
      ENDIF

      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS)             = USER_LEVELS(1:N_USER_LEVELS)

      VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS     = N_USER_OBSGEOMS
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,:) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,:)

!  Modified Chapman function inputs

      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS        = EARTH_RADIUS

!  Modified optical inputs

      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS) = OMEGA_TOTAL_INPUT(1:NLAYERS)

!  Modified linearized control variables

!      First three are additional variables, Version 2.7

      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY = DO_SIMULATION_ONLY
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF      = DO_ATMOS_LBBF
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF    = DO_SURFACE_LBBF

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = DO_COLUMN_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = DO_PROFILE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = DO_ATMOS_LINEARIZATION

      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = DO_SURFACE_LINEARIZATION
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = DO_LINEARIZATION

!  ==========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (pure OUT variables)
!  ==========================================================

!  Radiances and fluxes
!    Output direct Mean/Flux added 17 May 2012
!mick fix 1/5/2021 - shut off passing of STOKES as in 2.8.3 "vlidort_master"
!                    (output type structure now filled directly)

!      VLIDORT_Out%Main%TS_STOKES(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)    = &
!        STOKES(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)

!  1/31/21. Version 2.8.3. Integrated values - Only  copy if necessary.

      VLIDORT_Out%Main%TS_MEANST_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        MEANST_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)
      VLIDORT_Out%Main%TS_FLUX_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        FLUX_DIFFUSE(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)

      VLIDORT_Out%Main%TS_DNMEANST_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = &
        DNMEANST_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)
      VLIDORT_Out%Main%TS_DNFLUX_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = &
        DNFLUX_DIRECT(1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)

!  4/26/19. Media properties output.
!mick mod 1/5/2021 - added IF conditions to only fill what you need (as in 2.8.3 "vlidort_master")

      IF ( DO_ALBTRN_MEDIA(1) ) THEN
         VLIDORT_Out%Main%TS_ALBMED_USER(1:NSTOKES,1:N_USER_VZANGLES) = ALBMED_USER(1:NSTOKES,1:N_USER_VZANGLES)
         VLIDORT_Out%Main%TS_ALBMED_FLUXES(1:NSTOKES,1:2)             = ALBMED_FLUXES(1:NSTOKES,1:2)
      ENDIF
      IF ( DO_ALBTRN_MEDIA(2) ) THEN
         VLIDORT_Out%Main%TS_TRNMED_USER(1:NSTOKES,1:N_USER_VZANGLES) = TRNMED_USER(1:NSTOKES,1:N_USER_VZANGLES)
         VLIDORT_Out%Main%TS_TRNMED_FLUXES(1:NSTOKES,1:2)             = TRNMED_FLUXES(1:NSTOKES,1:2)
      ENDIF

!  1/31/21. Version 2.8.3. DO_MSSTS option final installation
!    ==> MSSTS auxiliary output is copied from VFO (see above)
!    ==> MSSTS main output is filled directly in Converge routines

!  4/28/19. Planetary problem output
!  ---------------------------------

!  1/31/21. Version 2.8.3. ( 5/5/20. Version 2.8.1  Upgrade )
!    ==> Revised code to use only the first index of TRANSBEAM ==> makes the Q-problem valid

!  Here is the Old code (pre 5/5/20 Upgrade) 
!       IF ( DO_PLANETARY_PROBLEM ) THEN
!         VLIDORT_Out%Main%TS_PLANETARY_SBTERM = TRNMED_FLUXES(1,2)
!         IF ( DO_OBSERVATION_GEOMETRY ) THEN
!            DO IB =1, N_GEOMETRIES ; DO O1 = 1, NSTOKES
!               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,IB) = TRANSBEAM(O1,IB) * TRNMED_USER(O1,IB) / PIE
!            ENDDO ; ENDDO
!         ELSE
!            DO IB =1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
!               OFF = VZA_OFFSETS(IB,UM) ; TRANS = TRANSBEAM(O1,IB) * TRNMED_USER(O1,UM) / PIE
!               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,OFF+1:OFF+N_USER_RELAZMS) = TRANS
!            ENDDO ; ENDDO ; ENDDO
!         ENDIF   
!      ENDIF

!  1/31/21. Version 2.8.3. ( 5/5/20. Version 2.8.1  Upgrade )
!   -- 1/31/21. Version 2.8.3. (2/16/21).  Add the Doublet offset settings

      IF ( DO_PLANETARY_PROBLEM ) THEN
         VLIDORT_Out%Main%TS_PLANETARY_SBTERM = TRNMED_FLUXES(1,2)
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            DO IB =1, N_GEOMETRIES ; DO O1 = 1, NSTOKES
               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,IB) = TRANSBEAM(1,IB) * TRNMED_USER(O1,IB) / PIE
            ENDDO ; ENDDO
         ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
            DO IB = 1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
               OFF = SZD_OFFSETS(IB) ; TRANS = TRANSBEAM(1,IB) * TRNMED_USER(O1,UM) / PIE
               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,OFF+UM) = TRANS
            ENDDO ; ENDDO ; ENDDO
         ELSE
            DO IB =1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
               OFF = VZA_OFFSETS(IB,UM) ; TRANS = TRANSBEAM(1,IB) * TRNMED_USER(O1,UM) / PIE
               VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(O1,OFF+1:OFF+N_USER_RELAZMS) = TRANS
            ENDDO ; ENDDO ; ENDDO
         ENDIF   
      ENDIF
      
!  Additional Control for Externalized input (SLEAVE). Introduced 3/18/19 for Version 2.8.1
!   -- 1/31/21. Version 2.8.3. Original code now moved to within Fouirer loop

!  new 12 March 2012
!   IF SS results already available, no need to copy them !
!mick fix 1/5/2021 - shut off passing of STOKES_SS & STOKES_DB - output type structures now filled directly

      !IF ( .NOT. DO_FOCORR_EXTERNAL ) THEN
      !   VLIDORT_Sup%SS%TS_STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:) = &
      !     STOKES_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)
      !   VLIDORT_Sup%SS%TS_STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
      !     STOKES_DB(1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
      !ENDIF

!  Bookkeeping
!   -- 1/31/21. Version 2.8.3. (2/16/21).  Add the Doublet offset settings

      VLIDORT_Out%Main%TS_FOURIER_SAVED(1:N_SZANGLES) = FOURIER_SAVED(1:N_SZANGLES)
      VLIDORT_Out%Main%TS_N_GEOMETRIES                = N_GEOMETRIES

      VLIDORT_Out%Main%TS_SZD_OFFSETS(1:N_SZANGLES)   = SZD_OFFSETS(1:N_SZANGLES)
      VLIDORT_Out%Main%TS_SZA_OFFSETS(1:N_SZANGLES)   = SZA_OFFSETS(1:N_SZANGLES)
      VLIDORT_Out%Main%TS_VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES) = VZA_OFFSETS(1:N_SZANGLES,1:N_USER_VZANGLES)

!  Solar Beam Transmittance to BOA
!  rob fix 11/17/2014, for diagnostic use only

      VLIDORT_Out%Main%TS_SOLARBEAM_BOATRANS (1:NBEAMS) = SOLARBEAM_BOATRANS (1:NBEAMS)

!  Column weighting functions
!    Output direct Mean/Flux added 17 May 2012
!mick fix 1/5/2021 - shut off passing of COLUMNWF - output type structure now filled directly

      !VLIDORT_LinOut%Col%TS_COLUMNWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:) = &
      !  COLUMNWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)

      VLIDORT_LinOut%Col%TS_MEANST_DIFFUSE_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        MEANST_DIFFUSE_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)

      VLIDORT_LinOut%Col%TS_FLUX_DIFFUSE_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        FLUX_DIFFUSE_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)

      VLIDORT_LinOut%Col%TS_DNMEANST_DIRECT_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = &
        DNMEANST_DIRECT_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)

      VLIDORT_LinOut%Col%TS_DNFLUX_DIRECT_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES) = &
        DNFLUX_DIRECT_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES)

!  4/26/19. Media properties output.
!mick mod 1/5/2021 - added IF conditions to only fill what you need (as in 2.8.3 "vlidort_master")

      IF ( DO_ALBTRN_MEDIA(1) ) THEN
        VLIDORT_LinOut%Col%TS_ALBMED_USER_COLWF(1:NSTOKES,1:N_USER_VZANGLES,:) = &
                                 LC_ALBMED_USER(1:NSTOKES,1:N_USER_VZANGLES,:)
        VLIDORT_LinOut%Col%TS_ALBMED_FLUXES_COLWF(1:NSTOKES,:,:)               = &
                                 LC_ALBMED_FLUXES(1:NSTOKES,:,:)
      ENDIF
      IF ( DO_ALBTRN_MEDIA(2) ) THEN
        VLIDORT_LinOut%Col%TS_TRNMED_USER_COLWF(1:NSTOKES,1:N_USER_VZANGLES,:) = &
                                 LC_TRNMED_USER(1:NSTOKES,1:N_USER_VZANGLES,:)
        VLIDORT_LinOut%Col%TS_TRNMED_FLUXES_COLWF(1:NSTOKES,:,:)               = &
                                 LC_TRNMED_FLUXES(1:NSTOKES,:,:)
      ENDIF

!  4/28/19. Planetary problem output
!  ---------------------------------

!  5/5/20. Version 2.8.1 Upgrades 
!    ==> Revised code to use only the first index of TRANSBEAM ==> makes the Q-problem valid
!    -- 1/31/21. Version 2.8.3. (2/16/21).  Add the Doublet offset settings

      IF ( DO_PLANETARY_PROBLEM .AND. DO_COLUMN_LINEARIZATION ) THEN
         VLIDORT_LinOut%Col%TS_PLANETARY_SBTERM_COLWF(:) = LC_TRNMED_FLUXES(1,2,:)
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            DO IB = 1, N_GEOMETRIES ; DO O1 = 1, NSTOKES
               DO Q = 1, N_TOTALCOLUMN_WFS
                  VLIDORT_LinOut%Col%TS_PLANETARY_TRANSTERM_COLWF(O1,IB,Q) = &
                          ( LC_TRANSBEAM(1,IB,Q) * TRNMED_USER(O1,IB) + TRANSBEAM(1,IB) * LC_TRNMED_USER(O1,IB,Q) ) / PIE
               ENDDO
            ENDDO ; ENDDO
         ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
            DO IB =1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
               OFF = SZD_OFFSETS(IB) ; TRANS = TRANSBEAM(1,IB) * TRNMED_USER(O1,UM) / PIE
               DO Q = 1, N_TOTALCOLUMN_WFS
                  L_TRANS(Q) =  ( LC_TRANSBEAM(1,IB,Q) * TRNMED_USER(O1,UM) + TRANSBEAM(1,IB) * LC_TRNMED_USER(O1,UM,Q) ) / PIE
                  VLIDORT_LinOut%Col%TS_PLANETARY_TRANSTERM_COLWF(O1,OFF+UM,Q) = L_TRANS(Q)
               ENDDO
            ENDDO ; ENDDO ; ENDDO
         ELSE
            DO IB =1, NBEAMS ; DO UM = 1, N_USER_VZANGLES ; DO O1 = 1, NSTOKES
               OFF = VZA_OFFSETS(IB,UM) ; TRANS = TRANSBEAM(1,IB) * TRNMED_USER(O1,UM) / PIE
               DO Q = 1, N_TOTALCOLUMN_WFS
                  L_TRANS(Q) =  ( LC_TRANSBEAM(1,IB,Q) * TRNMED_USER(O1,UM) + TRANSBEAM(1,IB) * LC_TRNMED_USER(O1,UM,Q) ) / PIE
                  VLIDORT_LinOut%Col%TS_PLANETARY_TRANSTERM_COLWF(O1,OFF+1:OFF+N_USER_RELAZMS,Q) = L_TRANS(Q)
               ENDDO
            ENDDO ; ENDDO ; ENDDO
         ENDIF   
      ENDIF
      
!  Superseded, Version 2.7
!      VLIDORT_LinOut%Atmos%TS_LTE_ATMOSWF(0:NLAYERS,1:N_USER_LEVELS,1:N_USER_VZANGLES,:) = &
!        LTE_ATMOSWF(0:NLAYERS,1:N_USER_LEVELS,1:N_USER_VZANGLES,:)

!  Surface weighting functions
!mick fix 1/5/2021 - shut off passing of SURFACEWF - output type structure now filled directly

      !VLIDORT_LinOut%Surf%TS_SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:) = &
      !  SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)

      VLIDORT_LinOut%Surf%TS_MEANST_DIFFUSE_SURFWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        MEANST_DIFFUSE_SURFWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)

      VLIDORT_LinOut%Surf%TS_FLUX_DIFFUSE_SURFWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:) = &
        FLUX_DIFFUSE_SURFWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_SZANGLES,1:NSTOKES,:)

!  SS weighting functions
!mick fix 1/5/2021 - shut off passing of COLUMNWF_SS, COLUMNWF_DB, & SURFACEWF_DB - output type structures now filled directly

      !IF ( .NOT. DO_FOCORR_EXTERNAL ) THEN
      !   VLIDORT_LinSup%SS%Col%TS_COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:) = &
      !     COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES,:)
      !   VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)   = &
      !     COLUMNWF_DB(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
      !   VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES) = &
      !     SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,1:NSTOKES)
      !ENDIF

!  BLACKBODY JACOBIANS, New 28 March 2014. 

      VLIDORT_LinOut%Atmos%TS_ABBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_VZANGLES,0:NLAYERS,1:NSTOKES,1:2) = &
          ABBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_VZANGLES,0:NLAYERS,1:NSTOKES,1:2)
      VLIDORT_LinOut%Atmos%TS_ABBWFS_FLUXES(1:N_USER_LEVELS,1:2,0:NLAYERS,1:NSTOKES,1:2) = &
          ABBWFS_FLUXES(1:N_USER_LEVELS,1:2,0:NLAYERS,1:NSTOKES,1:2)

      VLIDORT_LinOut%Surf%TS_SBBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_VZANGLES,1:NSTOKES,1:2) = &
          SBBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_VZANGLES,1:NSTOKES,1:2)
      VLIDORT_LinOut%Surf%TS_SBBWFS_FLUXES(1:N_USER_LEVELS,1:2,1:NSTOKES,1:2) = &
          SBBWFS_FLUXES(1:N_USER_LEVELS,1:2,1:NSTOKES,1:2)

!  Exception handling

      VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  = STATUS_INPUTCHECK
      VLIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION

      VLIDORT_Out%Status%TS_NCHECKMESSAGES = NCHECKMESSAGES
      VLIDORT_Out%Status%TS_CHECKMESSAGES(0:NCHECKMESSAGES)  = CHECKMESSAGES(0:NCHECKMESSAGES)
      VLIDORT_Out%Status%TS_ACTIONS(0:NCHECKMESSAGES)        = ACTIONS(0:NCHECKMESSAGES)

      VLIDORT_Out%Status%TS_MESSAGE = MESSAGE
      VLIDORT_Out%Status%TS_TRACE_1 = TRACE_1
      VLIDORT_Out%Status%TS_TRACE_2 = TRACE_2
      VLIDORT_Out%Status%TS_TRACE_3 = TRACE_3

!  ===================================
!  END COPY LOCAL VARIABLES TO OUTPUTS
!  ===================================

!  Major result output
!  ===================

!  Revised Version 2.8, 7/8/16. Cleaned up I/O listings (but not for the linearized routine....)
!mick mod 1/5/2021 - reactivated VLIDORT_WRITERESULTS and VLIDORT_L_WRITERESULTS:
!                      * moved subroutines from before output copy to here
!                      * updated call statements
!                      * switched to using output type structures VLIDORT_Out and VLIDORT_LinOut

      IF ( DO_WRITE_RESULTS ) THEN
        RUNIT = VLIDORT_RESUNIT
        OPEN(RUNIT,FILE=RESULTS_WRITE_FILENAME,STATUS='REPLACE')

        CALL VLIDORT_WRITERESULTS ( &
          RUNIT, DO_FULLRAD_MODE, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_DOUBLE_CONVTEST, &
          DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,   &
          NSTOKES, VLIDORT_ACCURACY, SZANGLES, N_USER_RELAZMS,                             &
          USER_RELAZMS, N_USER_LEVELS, USER_LEVELS, HEIGHT_GRID,  DELTAU_VERT_INPUT,       &
          DO_NO_AZIMUTH, NBEAMS, N_DIRECTIONS, WHICH_DIRECTIONS, N_OUT_STREAMS,            &
          OUT_ANGLES, PARTLAYERS_OUTFLAG, VZA_OFFSETS, TAUGRID_INPUT, DO_MULTIBEAM,        &
          VLIDORT_Out%Main, MEANST_DIFFUSE, FLUX_DIFFUSE, FOURIER_SAVED )

        IF ( DO_LINEARIZATION ) THEN
          CALL VLIDORT_L_WRITERESULTS ( &
            RUNIT, DO_FULLRAD_MODE, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_DOUBLE_CONVTEST, &
            DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,   &
            DO_LAMBERTIAN_SURFACE, DO_NO_AZIMUTH, DO_MULTIBEAM,                              &
            NSTOKES, NLAYERS, VLIDORT_ACCURACY, SZANGLES, N_USER_RELAZMS, USER_RELAZMS,      &
            N_USER_LEVELS, USER_LEVELS, NBEAMS, N_DIRECTIONS, WHICH_DIRECTIONS,              &
            N_OUT_STREAMS, OUT_ANGLES, VZA_OFFSETS, FOURIER_SAVED, &
            DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION,     &
            N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, & 
            DO_SURFACE_LINEARIZATION, N_SURFACE_WFS, &
            PROFILEWF_NAMES, COLUMNWF_NAMES, VLIDORT_LinOut )
        ENDIF

        CLOSE(RUNIT)
      ENDIF

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE VLIDORT_DEBUG_INPUT_MASTER()

!  3/28/14. Changes for Version 2.7. remove LTE linearization references
!  7/8/16. Removed DO_QUAD_OUTPUT for Version 2.8. I/O list needs a clean-up
!mick mod 9/19/2017 - added DO_FLUORESCENCE, DO_TF_ITERATION, TAYLOR_ORDER, FMATRIX_UP, FMATRIX_DN,
!                           ATMOS_WAVELENGTH, TF_MAXITER, TF_CRITERION to argument list

!  Additional Control for SLEAVE (DO_WLADJUSTED_OUTPUT,DO_EXTERNAL_WLEAVE).
!     Introduced 3/18/19 for Version 2.8.1
         
!  4/26/19. Record the Media-problem inputs     (DO_ALBTRN_MEDIA)
!  4/28/19. Record the Planetary-problem inputs (DO_PLANETARY_PROBLEM)
!  3/23/19, Version 2.8.1, Record Control for TOA/BOA illumination
        
!  1/31/21. Version 2.8.3. Add 3 new input flags. 
!    -- 3 new flags are: DO_CLASSICAL_SOLUTION, DO_MSSTS, DO_DOUBLET_GEOMETRY
!    -- re-order some of the Boolean input arguments (first 8 lines)
!    -- Drop SSCORR_TRUNCATION (disabled). Add ASYMTX_TOLERANCE
!    -- (RTS 2/16/21). Add DO_FOURIER0_NSTOKES2 flag

      CALL VLIDORT_WRITE_STD_INPUT ( DO_DEBUG_WRITE, &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_THERMAL_TRANSONLY, DO_SURFACE_EMISSION,                   & ! Sources Main
        DO_TOAFLUX, DO_BOAFLUX, DO_ALBTRN_MEDIA, DO_PLANETARY_PROBLEM,                                      & ! Sources Other
        DO_FULLRAD_MODE, DO_UPWELLING, DO_DNWELLING, DO_CLASSICAL_SOLUTION, DO_MSSTS, DO_FOURIER0_NSTOKES2, & ! Main RT control
        DO_RAYLEIGH_ONLY, DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, DO_CHAPMAN_FUNCTION,                   & ! RT Control
        DO_FOCORR, DO_FOCORR_EXTERNAL, DO_FOCORR_NADIR,  DO_FOCORR_OUTGOING, DO_SSCORR_USEFMAT,             & ! FO (first-order)
        DO_DELTAM_SCALING, DO_DOUBLE_CONVTEST, DO_SOLUTION_SAVING,  DO_BVP_TELESCOPING,                     & ! RT Performance
        DO_USER_VZANGLES, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, & ! RT Post-processing
        DO_TOA_CONTRIBS,  DO_SPECIALIST_OPTION_1,  DO_SPECIALIST_OPTION_2, DO_SPECIALIST_OPTION_3,          & ! RT Specialist
        DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_EXTERNAL_WLEAVE,                     & ! Surface Control
        DO_WATER_LEAVING, DO_FLUORESCENCE, DO_TF_ITERATION, DO_WLADJUSTED_OUTPUT, TF_MAXITER, TF_CRITERION, & ! Lw Control
        TAYLOR_ORDER, NSTOKES, NSTREAMS, NLAYERS, NFINELAYERS, N_THERMAL_COEFFS, NGREEK_MOMENTS_INPUT, & ! Main numbers
        NLAYERS_NOMS, NLAYERS_CUTOFF, VLIDORT_ACCURACY, FLUX_FACTOR, EARTH_RADIUS, RFINDEX_PARAMETER,  & ! Flux/Acc/Radius
        ASYMTX_TOLERANCE, N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES, N_USER_OBSGEOMS, N_USER_LEVELS, & ! Geometry/Output control
        SZANGLES,   USER_RELAZMS,   USER_VZANGLES,   USER_OBSGEOMS,   USER_LEVELS,                     & ! Geometry/Output control
        HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, FINEGRID, GEOMETRY_SPECHEIGHT, ATMOS_WAVELENGTH, & ! Atmospheric inputs
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, FMATRIX_UP, FMATRIX_DN,            & ! Optical Inputs
        THERMAL_BB_INPUT, ALBEDO, SURFBB, TOAFLUX, BOAFLUX,                                            & ! Optical Inputs
        DO_WRITE_INPUT,    INPUT_WRITE_FILENAME,    DO_WRITE_FOURIER,       DO_WRITE_RESULTS,          & ! Debug control
        DO_WRITE_SCENARIO, SCENARIO_WRITE_FILENAME, FOURIER_WRITE_FILENAME, RESULTS_WRITE_FILENAME )     ! Debug control

!  1/31/21. Version 2.8.3. (2/16/21). BRDF arrays input directly from type structure. Restore Full Fourier loops
!mick mod 9/19/2017 - added NMOMENTS to input. No longer required

      IF (.NOT. DO_LAMBERTIAN_SURFACE) THEN
        CALL VLIDORT_WRITE_SUP_BRDF_INPUT ( &
          DO_USER_VZANGLES, DO_SURFACE_EMISSION, NSTOKES, NSTREAMS, N_SZANGLES,   &
          N_USER_VZANGLES, N_USER_RELAZMS,   VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC, &
          VLIDORT_Sup%BRDF%TS_BRDF_F_0,      VLIDORT_Sup%BRDF%TS_BRDF_F,          &
          VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0, VLIDORT_Sup%BRDF%TS_USER_BRDF_F,     &
          VLIDORT_Sup%BRDF%TS_EMISSIVITY,    VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY )
      END IF

!  1/31/21. Version 2.8.3. (2/16/21). SS arrays input directly from type structure. 

      IF (DO_FOCORR_EXTERNAL) THEN
        CALL VLIDORT_WRITE_SUP_SS_INPUT ( &
          NSTOKES, N_USER_LEVELS, VLIDORT_Sup%SS%TS_STOKES_SS, VLIDORT_Sup%SS%TS_STOKES_DB )
      END IF

!  1/31/21. Version 2.8.3. (2/16/21). SLEAVE arrays input directly from type structure. Restore Full Fourier loops

      IF (DO_SURFACE_LEAVING) THEN
        CALL VLIDORT_WRITE_SUP_SLEAVE_INPUT ( &
          DO_USER_VZANGLES, NSTOKES, NSTREAMS, N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS, &
          VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC, VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES,  &
          VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0,       VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0 )
      END IF

      END SUBROUTINE VLIDORT_DEBUG_INPUT_MASTER



      SUBROUTINE VLIDORT_DEBUG_LIN_INPUT_MASTER()

!  3/28/14. Changes for Version 2.7. remove LTE linearization references. Add LBBF

      CALL VLIDORT_WRITE_LIN_INPUT ( &
        NSTOKES, NLAYERS, NGREEK_MOMENTS_INPUT, N_SZANGLES, N_USER_RELAZMS, &
        N_USER_VZANGLES, DO_OBSERVATION_GEOMETRY, DO_SIMULATION_ONLY,       &
        N_TOTALCOLUMN_WFS, N_TOTALPROFILE_WFS, N_SURFACE_WFS, N_SLEAVE_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                                 &
        COLUMNWF_NAMES, PROFILEWF_NAMES, DO_ATMOS_LBBF, DO_SURFACE_LBBF,    &
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,   &
        L_FMATRIX_UP, L_FMATRIX_DN,                                         &
        DO_COLUMN_LINEARIZATION, DO_PROFILE_LINEARIZATION, DO_ATMOS_LINEARIZATION,&
        DO_SURFACE_LINEARIZATION, DO_LINEARIZATION, DO_SLEAVE_WFS )

!  1/31/21. Version 2.8.3. (2/16/21). BRDF arrays input directly from type structure. Restore Full Fourier loops
!mick mod 9/19/2017 - added NMOMENTS to input. No longer required

      IF (.NOT. DO_LAMBERTIAN_SURFACE .AND. DO_SURFACE_LINEARIZATION) THEN
        CALL VLIDORT_WRITE_LIN_SUP_BRDF_INPUT ( &
          DO_USER_VZANGLES, DO_SURFACE_EMISSION, NSTOKES, NSTREAMS, NBEAMS, N_USER_VZANGLES,  &
          N_USER_RELAZMS, N_SURFACE_WFS,           VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC, &
          VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0,      VLIDORT_LinSup%BRDF%TS_LS_BRDF_F,          &
          VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0, VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F,     &
          VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY,    VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY )
      END IF

!  1/31/21. Version 2.8.3. (2/16/21). SS arrays input directly from type structure.
!    -- Follow the LIDORT pattern, dedicated subroutine only for profile output.

      IF (DO_FOCORR_EXTERNAL) THEN
        CALL VLIDORT_WRITE_LCS_SUP_SS_INPUT ( &
          DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION,                      &
          NSTOKES, NLAYERS, N_USER_LEVELS, N_TOTALCOLUMN_WFS, N_TOTALSURFACE_WFS, &
          VLIDORT_LinSup%SS%Col%TS_COLUMNWF_SS, &
          VLIDORT_LinSup%SS%Col%TS_COLUMNWF_DB, &
          VLIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB )
      END IF

!  1/31/21. Version 2.8.3. (2/16/21). SLEAVE arrays input directly from type structure. Restore Full Fourier loops

      IF (DO_SURFACE_LEAVING .AND. DO_SURFACE_LINEARIZATION .AND. DO_SLEAVE_WFS) THEN
        CALL VLIDORT_WRITE_LIN_SUP_SLEAVE_INPUT ( &
          DO_USER_VZANGLES, NSTOKES, NSTREAMS, N_SZANGLES, N_USER_VZANGLES, N_USER_RELAZMS, N_SLEAVE_WFS,  &
          VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC, VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES, &
          VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0,       VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0 )
      END IF

      END SUBROUTINE VLIDORT_DEBUG_LIN_INPUT_MASTER

      END SUBROUTINE VLIDORT_LCS_MASTER

!  End Module

      END MODULE vlidort_lcs_masters1_m

