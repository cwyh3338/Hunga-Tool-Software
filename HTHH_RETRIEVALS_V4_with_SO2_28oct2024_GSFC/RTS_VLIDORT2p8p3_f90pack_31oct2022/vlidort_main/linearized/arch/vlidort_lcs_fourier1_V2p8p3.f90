
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
! #            VLIDORT_LCS_FOURIER (master)                     #
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

      MODULE vlidort_lcs_fourier1_m

      PUBLIC  :: VLIDORT_LCS_FOURIER

      CONTAINS
!
!  1/31/21. Version 2.8.3. Several additional changes to this argument list
!    -- Include flags DO_CLASSICAL_SOLUTION (Line 3), DO_MSSTS (line 8)
!    -- Include MSSTS outputs LAYER_MSSTS_F, SURF_MSSTS_F
!    -- Include ASYMTX Tolerance variable to the list (Line 8)
!    -- Use PPSTREAM masking system, replaces LOCAL_UM_START.
!    - (5/24/21). Add DO_INCLUDE_SLEAVEWFS

      SUBROUTINE VLIDORT_LCS_FOURIER ( FOURIER, DO_FOURIER0_NSTOKES2, &
            DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY,                   & !Input flags (RT operation)
            DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,           & !Input flags (RT operation)
            DO_PLANE_PARALLEL, DO_LAYER_SCATTERING, DO_SOLUTION_SAVING, DO_BVTEL_INITIAL,           & !Input flags (performance)
            DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT, DO_MULTIBEAM, DO_ALBTRN_MEDIA,                & !Input flags (Beam/Planetary)
            DO_PLANETARY_PROBLEM, DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,       & !Input flags (Planetary/Surf)
            DO_MSMODE_THERMAL, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_THERMAL_TRANSONLY,      & !Input flags (thermal)
            DO_FOCORR_ALONE, DO_DBCORRECTION, DO_TOA_CONTRIBS, DO_PARTLAYERS, DO_REAL_EIGENSOLVER,  & !Input Bookkeeping
            DO_MSSTS, DO_TOAFLUX, TOAFLUX, DO_BOAFLUX, BOAFLUX, TOLERANCE, DO_DEBUG_WRITE,          & !Input MSSTS & TOA/BOA fluxes
            DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION, TF_MAXITER, TF_CRITERION,        & !Input Water-leaving control
            DO_ATMOS_LINEARIZATION, DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION,              & !Input Linearization flags
            DO_SIMULATION_ONLY, DO_SLEAVE_WFS, DO_ATMOS_LBBF, DO_SURFACE_LBBF,                      & !Input Linearization flags
            NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_LEVELS, N_THERMAL_COEFFS,    & !Input Basic integers
            NMOMENTS, NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,       & !Input Derived integers
            N_PPSTREAMS, PPSTREAM_MASK, N_PARTLAYERS, N_DIRECTIONS, TAYLOR_ORDER,                   & !Input Derived integers
            N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,                    & !In Num(Lin)
            N_ALLLAYERS_UP, N_ALLLAYERS_DN, FLUX_FACTOR, FLUXVEC, COS_SZANGLES, SZA_LOCAL_INPUT, SUN_SZA_COSINES,  & !In SZAs/FLUX
            USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,                    & !In Streams
            MUELLER_INDEX, DMAT, BVP_REGULAR_FLAG, WHICH_DIRECTIONS, SOLARBEAM_BOATRANS,                           & !In Bkkeep
            LEVELMASK_UP, LEVELMASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,              & !In Bkkeep
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, LC_SOLARBEAM_BOATRANS, LC_SOLARBEAM_ATRANS,             & !In Bkkeep
            Misc, Therm, Mult, LAC_Misc, L_Therm, LC_Mult,                                                         & !In Packing
            ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0,                                                  & !In ALB/BRDF
            SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0, SURFBB, EMISSIVITY, USER_EMISSIVITY,                    & !In SL/EMISS
            LS_BRDF_F_0, LS_BRDF_F,  LS_USER_BRDF_F_0, LS_USER_BRDF_F, LS_EMISSIVITY, LS_USER_EMISSIVITY,          & !In LS
            LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,                                          & !In LSSL
            PIMM_11, PIMM_KM, BVTEL_FOURIER, DO_INCLUDE_THERMEMISS, DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS,      & !Out Bkkeep
            STOKES_F, MEANST_DIFFUSE, FLUX_DIFFUSE, DNMEANST_DIRECT, DNFLUX_DIRECT, MS_CONTRIBS_F,                 & !Out STOKES
            COLUMNWF_F,  MEANST_DIFFUSE_COLWF,  FLUX_DIFFUSE_COLWF,  DNMEANST_DIRECT_COLWF, DNFLUX_DIRECT_COLWF,   & !Out COL JACS
            SURFACEWF_F, MEANST_DIFFUSE_SURFWF, FLUX_DIFFUSE_SURFWF, DO_INCLUDE_SLEAVEWFS,                         & !Out SURF JACS
            ABBWFS_JACOBIANS, ABBWFS_FLUXES, SBBWFS_JACOBIANS, SBBWFS_FLUXES,                                      & !Out LBBF JACS
            TRANS_ATMOS_FINAL, TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,                  & !Out 4/26/19
            LC_TRANS_ATMOS_FINAL, LC_TRANSBEAM, LC_ALBMED_USER, LC_ALBMED_FLUXES, LC_TRNMED_USER, LC_TRNMED_FLUXES,& !Out 4/26/19
            LAYER_MSSTS_F, SURF_MSSTS_F, LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F,     & !Output SPECIAL
            STATUS, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )                                                             !Out Status

!  This line was removed, 7/8/16.
!        DO_SPECIALIST_OPTION_2, DO_DEBUG_WRITE, DO_FDTEST,                                         & !Input

!   Version 2.8.1, Control for TOA/BOA isotropic illumination added, 3/23/19
!      --- NOT THE SAME AS THE ALBTRN_MEDIA condition.         
        
!  Argument list revised for Version 2.8.1, 4/9/19, 3/23/19
!    1.    Water-leaving adjustment control added (DO_WATER_LEAVING, DO_TF_ITERATION, TF_MAXITER, TF_CRITERION)
!    2.    TOA/BOA isotropic illumination control, DO_TOAFLUX, TOAFLUX, DO_BOAFLUX, BOAFLUX    
!    3.    TRANS_ATMOS_FINAL, added to the output list (Self-adjusting water-leaving formulation)
!    4.    SOLARBEAM_BOATRANS, TRANS_SOLAR_BEAM distinguished.
      
!  Argument list revised for Version 2.8.1, 4/26/19, 4/28/19
!  4/26/19 Module for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!    -- introduced by R. Spurr 4/26/19. Controlled by flags DO_ALBTRN_MEDIA, DO_PLANETARY_PROBLEM
     
!  Complete Fourier component calculation for the Extended Code
!    Stokes vector and Weighting function computations

!  1/31/21. Version 2.8.3. Changes include
!    - introduction of flags DO_CLASSICAL_SOLUTION, DO_MSSTS (control Green's function, MSST output)
!    - Introduction of TOLERANCE variable for the ASYMTX ekgenroutine
!    - BRDF/SLEAVE Fourier arrays defined locally for each Fourier component
!    - MSSTS outputs LAYER_MSSTS_F, SURF_MSSTS_F now added to argument list.
!    - Additional unpacking required (SIGMA_M, SIGMA_P, ITRANS_USERM) for Green's function solution
!    - Complete set of internal arrays defined for the Green's function solution
!    - Use of post-processing mask. HMULT_MASTER now in solutions module
!    - (RTS 2/16/21). Introduce Fourier0_NSTOKES2 flag for fast calculation
!    - (5/24/21). Add DO_INCLUDE_SLEAVEWFS

      USE VLIDORT_PARS_m

!  Internal types

      USE VLIDORT_Work_def_m
      USE VLIDORT_LinWork_def_m

!  Standard code modules, restricted use
!   1/31/21. Version 2.8.3.  MULTIPLIERS module has been dispersed.

      USE VLIDORT_MISCSETUPS_m , Only : VLIDORT_DIRECTRADIANCE, VLIDORT_PIMATRIX_SETUP_OMP
      USE VLIDORT_THERMALSUP_m , Only : THERMAL_CLSOLUTION, THERMAL_STERMS_UP, THERMAL_STERMS_DN
      USE VLIDORT_SOLUTIONS_m
      USE VLIDORT_BVPROBLEM_m  , Only : BVP_MATRIXSETUP_MASTER,    BVP_SOLUTION_MASTER, &
                                        BVPTEL_MATRIXSETUP_MASTER, BVPTEL_SOLUTION_MASTER
      USE VLIDORT_INTENSITY_m  , Only : VLIDORT_UPUSER_INTENSITY, VLIDORT_DNUSER_INTENSITY, VLIDORT_INTEGRATED_OUTPUT

!  4/26/19. Media-properties routines, Mark II, 4/28/19.
   
      USE VLIDORT_MediaProps_m
      USE VLIDORT_LC_MediaProps_m

!  Valid for both LCS and LPS (this subroutine)

      USE VLIDORT_LPC_SOLUTIONS_m
      USE VLIDORT_L_THERMALSUP_m, Only : THERMAL_CLSOLUTION_PLUS, THERMAL_STERMS_UP_PLUS, THERMAL_STERMS_DN_PLUS

!  valid only for LC linearization (this subroutine)

      USE VLIDORT_LC_SOLUTIONS_m
      USE VLIDORT_LC_BVPROBLEM_m , Only : LC_BVP_SOLUTION_MASTER, LC_BVPTEL_SOLUTION_MASTER
      USE VLIDORT_LC_WFATMOS_m   , Only : UPUSER_COLUMNWF, DNUSER_COLUMNWF, MIFLUX_COLUMNWF

!  Addition of LBBF Jacobians, Version 2.7, 3/28/14. Enabled for Version 2.8, April 2016

      use vlidort_lbbf_jacobians_m

!  Surface linearization (present in LCS or LPS)
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@

      USE VLIDORT_LS_WFSURFACE_m, Only : SURFACEWF_BVP_SOLUTION, SURFACEWF_BVPTEL_SOLUTION, SURFACEWF_POSTPROCESS_MASTER
      USE VLIDORT_LS_WFSLEAVE_m

!  Unpack

      USE VLIDORT_UNPACK_m
      USE VLIDORT_L_UNPACK_m
      USE VLIDORT_LC_UNPACK_m

      IMPLICIT NONE

!  INPUT, INTENT(IN) ONLY
!  ======================

!  Input Fourier component number

      INTEGER, intent(in)  :: FOURIER

!  1/31/21. Version 2.8.3. (RTS 2/16/21). 
!    -- Introduce Flag for using NSTOKES = 2 for Fourier 0 (set by hand)

      LOGICAL, INTENT(IN) ::           DO_FOURIER0_NSTOKES2

!  Flags
!  -----

!  SS control

!      LOGICAL, INTENT (IN) ::          DO_SSCORR_NADIR  ! Dropped 2.8
      LOGICAL, INTENT (IN) ::          DO_FOCORR_ALONE   ! Renamed 2.8

!  Solar control

      LOGICAL, INTENT (IN) ::          DO_SOLAR_SOURCES
      LOGICAL, INTENT (IN) ::          DO_REFRACTIVE_GEOMETRY

!  RT control

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::          DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT (IN) ::          DO_MVOUT_ONLY
      LOGICAL, INTENT (IN) ::          DO_SOLUTION_SAVING
      LOGICAL, INTENT (IN) ::          DO_MSMODE_VLIDORT
      LOGICAL, INTENT (IN) ::          DO_PLANE_PARALLEL

!  1/31/21. Version 2.8.3. Add DO_CLASSICAL_SOLUTION to list of arguments (Line 3)

      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION

!  Convergence tracker

      LOGICAL, INTENT (IN) ::          DO_MULTIBEAM ( MAXBEAMS, 0:MAXFOURIER )

!  4/26/19. Control flags for the isotropic-illuminated Media calculations

      LOGICAL, intent(in)  ::          DO_ALBTRN_MEDIA(2)
      
!  4/28/19  Added control for the planetary problem      

      LOGICAL, intent(in)  ::          DO_PLANETARY_PROBLEM
      
!  4/9/19. Additional Water-leaving control

      LOGICAL, intent(in)  ::          DO_WATER_LEAVING
      LOGICAL, intent(in)  ::          DO_EXTERNAL_WLEAVE
      LOGICAL, intent(in)  ::          DO_TF_ITERATION
      INTEGER, INTENT (IN) ::          TF_MAXITER
      DOUBLE PRECISION, INTENT (IN) :: TF_CRITERION

!  surface

      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE
!      LOGICAL, INTENT (IN) ::          DO_DIRECT_BEAM
      LOGICAL, INTENT (IN) ::          DO_DBCORRECTION
      LOGICAL, INTENT (IN) ::          DO_SURFACE_EMISSION

!  surface leaving

      LOGICAL, INTENT (IN) ::          DO_SURFACE_LEAVING
      LOGICAL, INTENT (IN) ::          DO_SL_ISOTROPIC

!  Output control

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_PARTLAYERS

!  thermal

      LOGICAL, INTENT (IN) ::          DO_THERMAL_TRANSONLY
      LOGICAL, INTENT (IN) ::          DO_MSMODE_THERMAL
      LOGICAL, INTENT (IN) ::          DO_THERMAL_EMISSION

!  Bookkeeping

      LOGICAL, INTENT (IN) ::          DO_REAL_EIGENSOLVER  ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING  ( 0:MAXMOMENTS, MAXLAYERS )

!  Other flags

      LOGICAL, INTENT (IN) ::          DO_TOA_CONTRIBS

!  1/31/21. Version 2.8.3. Flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::          DO_MSSTS

!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19

      LOGICAL         , INTENT(IN) :: DO_TOAFLUX, DO_BOAFLUX 
      DOUBLE PRECISION, INTENT(IN) :: TOAFLUX, BOAFLUX

!  debug write flag

      LOGICAL, INTENT (IN) ::          DO_DEBUG_WRITE

!  Disabled flags, removed, 7/8/16 for Version 2.8

      !LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      !LOGICAL, INTENT (IN) ::          DO_FO_CALC !New 02 Jul 2013
      !LOGICAL, INTENT (IN) ::          DO_SSCORR_TRUNCATION
      !LOGICAL, INTENT (IN) ::          DO_DELTAM_SCALING
      !LOGICAL, INTENT (IN) ::          DO_BVP_TELESCOPING
!      LOGICAL, INTENT (IN) ::          DO_FDTEST
!      LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_2
!      LOGICAL, INTENT (IN) ::          DO_SPECIALIST_OPTION_3
!      LOGICAL, INTENT (IN) ::          DO_QUAD_OUTPUT             ! removed, Version 2.8 7/8/16
!      LOGICAL, INTENT (IN) ::          DO_CLASSICAL_SOLUTION      ! removed, Version 2.8 7/8/16

!  Integers
!  --------

!  Basic

      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_THERMAL_COEFFS

!  derived

!      INTEGER, INTENT (IN) ::          NSTOKES_SQ ! removed 2.8
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL

      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2
      INTEGER, INTENT (IN) ::          N_PARTLAYERS
      INTEGER, INTENT (IN) ::          N_DIRECTIONS

      INTEGER, INTENT (IN) ::          N_ALLLAYERS_UP
      INTEGER, INTENT (IN) ::          N_ALLLAYERS_DN

!  Version 2p7 input, 2/19/14
      INTEGER, INTENT (IN) ::          TAYLOR_ORDER

!  1/31/21. Version 2.8.3. Post-processing masks. 

      INTEGER, INTENT (IN) ::           N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Disabled

      !INTEGER, INTENT (IN) ::          NGREEK_MOMENTS_INPUT
      !INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      !INTEGER, INTENT (IN) ::          NFINELAYERS
      !INTEGER, INTENT (IN) ::          N_GEOMETRIES
      !INTEGER, INTENT (IN) ::          NLAYERS_CUTOFF

!  Floating point and Misc.
!  ------------------------

!  1/31/21. Version 2.8.3.  Add the tolerance variable

      DOUBLE PRECISION, INTENT (IN) :: TOLERANCE

!  FLux control

      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      DOUBLE PRECISION, INTENT (IN) :: FLUXVEC ( MAXSTOKES )

!  SZA control

      DOUBLE PRECISION, INTENT (IN) :: COS_SZANGLES ( MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAX_SZANGLES )
      DOUBLE PRECISION, INTENT (IN) :: SUN_SZA_COSINES ( MAXLAYERS, MAX_SZANGLES )

!  User-angle control

      DOUBLE PRECISION, INTENT (IN) :: USER_STREAMS  ( MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_SECANTS  ( MAX_USER_STREAMS )

!  Quadrature (discrete ordinates)

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )

!  solution bookkeeping
!   -- 1/31/21. Version 2.8.3. Now using Lattice/Doublet/Obsgeom post-processing mask. Removed LOCAL_UM_START

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DMAT ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DFLUX ( MAXSTOKES )
      LOGICAL, INTENT (IN) ::          BVP_REGULAR_FLAG ( 0:MAXMOMENTS )
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Rob fix 11/27/14. Proxy for new output.
!     NOTE. SOLARBEAM_BOATRANS is computed with Unscaled optical depths, TRANS_SOLAR_BEAM with scaled ODs.

      DOUBLE PRECISION, INTENT (IN) :: SOLARBEAM_BOATRANS ( MAXBEAMS )

!  Level output control

      INTEGER, INTENT (IN) ::          LEVELMASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          LEVELMASK_DN  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )

!  Disabled

      !DOUBLE PRECISION, INTENT (IN) :: SS_FLUX_MULT
      !DOUBLE PRECISION, INTENT (IN) :: SZANGLES ( MAX_SZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: SZANGLES_ADJUST ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      !DOUBLE PRECISION, INTENT (IN) :: SIN_SZANGLES ( MAX_SZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES ( MAX_USER_VZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: USER_VZANGLES_ADJUST ( MAX_USER_VZANGLES )
      !INTEGER, INTENT (IN) ::          VZA_OFFSETS  ( MAX_SZANGLES, MAX_USER_VZANGLES )
      !DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      !DOUBLE PRECISION, INTENT (IN) :: USER_RELAZMS_ADJUST ( MAX_USER_VZANGLES, MAX_SZANGLES, MAX_USER_RELAZMS )
      !DOUBLE PRECISION, INTENT (IN) :: USER_LEVELS   ( MAX_USER_LEVELS )
      !INTEGER, INTENT (IN) ::          LAYER_MAXMOMENTS ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: PARTLAYERS_VALUES ( MAX_PARTLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: CHAPMAN_FACTORS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: EARTH_RADIUS
      !DOUBLE PRECISION, INTENT (IN) :: HEIGHT_GRID ( 0:MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: GREEKMAT_TOTAL_INPUT( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
      !DOUBLE PRECISION, INTENT (IN) :: TAUGRID_INPUT ( 0:MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  surface

      DOUBLE PRECISION, INTENT (IN) :: SURFBB

!  Work type structures

      TYPE(VLIDORT_Work_Miscellanous), INTENT (IN) :: Misc
      TYPE(VLIDORT_Work_Thermal),      INTENT (IN) :: Therm
      TYPE(VLIDORT_Work_Multiplier),   INTENT (IN) :: Mult

!  Linearized input
!  ----------------

!  flags and numbers
      
      LOGICAL, INTENT (INOUT) ::       DO_SIMULATION_ONLY        ! Intent changed, Version 2.7
      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT (IN) ::          DO_SLEAVE_WFS
      LOGICAL, INTENT (IN) ::          LAYER_VARY_FLAG   ( MAXLAYERS )
      
      INTEGER, INTENT (IN) ::          LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT (IN) ::          N_TOTALCOLUMN_WFS
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS

!  Control for  Blackbody Jacobians, New 28 March 2014. Version 2.7

      LOGICAL, INTENT (INOUT)  :: DO_ATMOS_LBBF, DO_SURFACE_LBBF

!  Work type structures

      TYPE(VLIDORT_LinWork_Miscellanous), INTENT (IN) :: LAC_Misc
      TYPE(VLIDORT_LinWork_Thermal),      INTENT (IN) :: L_Therm
      TYPE(VLIDORT_LinWork_Multiplier),   INTENT (IN) :: LC_Mult

!  Disabled

      !DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      !DOUBLE PRECISION, INTENT (IN) :: L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
!      DOUBLE PRECISION, INTENT (IN) :: LTE_DELTAU_VERT_INPUT ( 2, MAXLAYERS )
!      DOUBLE PRECISION, INTENT (IN) :: LTE_THERMAL_BB_INPUT ( 0:MAXLAYERS )
!      LOGICAL ::            DO_LTE_LINEARIZATION
!      LOGICAL ::            DO_SURFBB_LINEARIZATION

!  Surface, BRDF and SLEAVE
!  ------------------------

!  Albedo

      DOUBLE PRECISION, INTENT (IN) :: ALBEDO

!  BRDF arrays
!    -- 1/31/21. Version 2.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

      DOUBLE PRECISION, INTENT (IN) :: BRDF_F        ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F_0      ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_BRDF_F_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivity

      DOUBLE PRECISION, INTENT (IN) :: EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  New surface-leaving stuff 17 May 2012
!    -- 1/31/21. Version 2.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_ISOTROPIC ( MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   SLTERM_F_0       ( MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   USER_SLTERM_F_0  ( MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!  Disabled

      !DOUBLE PRECISION, INTENT (IN) :: EXACTDB_BRDFUNC    ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) :: SLTERM_USERANGLES  ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

!  Linearized BRDF and SLEAVE
!  --------------------------

!  Linearized BRDF stuff
!    -- 1/31/21. Version 2.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F        ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_BRDF_F_0      ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS )

!  Linearized Emissivity

      DOUBLE PRECISION, INTENT (IN) :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAXSTOKES,MAX_USER_STREAMS )
      DOUBLE PRECISION, INTENT (IN) :: LS_EMISSIVITY      ( MAX_SURFACEWFS, MAXSTOKES,MAXSTREAMS )

!  Addition of Linearized SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@
!    -- 1/31/21. Version 2.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_ISOTROPIC ( MAX_SLEAVEWFS, MAXSTOKES, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_F_0       ( MAX_SLEAVEWFS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::   LSSL_USER_SLTERM_F_0  ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS )

!Rob  fix 4/9/19    - added LC_SOLARBEAM_BOATRANS, LC_SOLARBEAM_ATRANS

      DOUBLE PRECISION, intent(in)  :: LC_SOLARBEAM_BOATRANS  ( MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION, intent(in)  :: LC_SOLARBEAM_ATRANS    ( MAXBEAMS, MAX_ATMOSWFS )

!  Disabled

      !DOUBLE PRECISION, INTENT (IN) :: LS_EXACTDB_BRDFUNC &
      !    ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      !DOUBLE PRECISION, INTENT (IN) ::   LSSL_SLTERM_USERANGLES &
      !    ( MAX_SLEAVEWFS, MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

!  OUTPUTS, Intent(INOUT)
!  ======================

!  Flags

      LOGICAL, INTENT (INOUT) ::          DO_BVTEL_INITIAL
      INTEGER, INTENT (INOUT) ::          BVTEL_FOURIER

!  Help arrays for the PI MAtrix Setup to make if thread-safe

      DOUBLE PRECISION, INTENT (INOUT) :: PIMM_11 ( MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION, INTENT (INOUT) :: PIMM_KM ( MAX_ALLSTRMS_P1 )

!  Regular output
!  --------------
      
!  Mean-value results will be done for Fourier = 0, then input again.
 
      DOUBLE PRECISION, INTENT (INOUT) :: MEANST_DIFFUSE ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: FLUX_DIFFUSE ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) :: DNMEANST_DIRECT ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) :: DNFLUX_DIRECT ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  First Introduced 3/22/17 for LIDORT, based on VLIDORT code 

      DOUBLE PRECISION, INTENT (INOUT) :: TRANS_ATMOS_FINAL  ( MAX_SZANGLES )

!  4/26/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Output for User-angles and fluxes. Output of Beam transmittance for Planetary problem.

      DOUBLE PRECISION, INTENT (INOUT) :: ALBMED_USER ( MAXSTOKES, MAX_USER_STREAMS ), &
                                          ALBMED_FLUXES( MAXSTOKES, 2 ) !  TOA illumination
      DOUBLE PRECISION, INTENT (INOUT) :: TRNMED_USER ( MAXSTOKES, MAX_USER_STREAMS ), &
                                          TRNMED_FLUXES( MAXSTOKES, 2 ) !  BOA illumination

!  4/28/19. Special Output of Beam transmittance for Planetary problem.

      DOUBLE PRECISION, intent(inout)  :: TRANSBEAM   ( MAXSTOKES, MAXBEAMS )

!  Linearized results
!  ------------------

!  Linearized Flux outputs
      
      DOUBLE PRECISION, INTENT (INOUT) ::  MEANST_DIFFUSE_COLWF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  DNMEANST_DIRECT_COLWF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_DIFFUSE_COLWF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  DNFLUX_DIRECT_COLWF &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) ::  MEANST_DIFFUSE_SURFWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (INOUT) ::  FLUX_DIFFUSE_SURFWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  5/24/21. Version 2.8.3. Must output the include SLEAVEWFS flag

      LOGICAL  , intent(inout) :: DO_INCLUDE_SLEAVEWFS

!  New Code Version 2.7. LBBF outputs. Superceded LTE_ATMOSWF output (now dropped)
!     Outputs are all Pre-zeroed in the calling Masters
!     Postprocessed and Flux Jacobians.

      DOUBLE PRECISION, INTENT(INOUT) :: ABBWFS_JACOBIANS &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS)
      DOUBLE PRECISION, INTENT(INOUT) :: SBBWFS_JACOBIANS &
          ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXSTOKES, MAX_DIRECTIONS)

      DOUBLE PRECISION, INTENT(INOUT) :: ABBWFS_FLUXES &
          ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT(INOUT) :: SBBWFS_FLUXES &
          ( MAX_USER_LEVELS, 2, MAXSTOKES, MAX_DIRECTIONS)

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  The LC Jacobian is an approximation.
!    --  Also depends on surface/sleave quantities, but these linearizations (LS/LSSL) are neglected

      DOUBLE PRECISION, INTENT (INOUT) :: LC_TRANS_ATMOS_FINAL ( MAXBEAMS, MAX_ATMOSWFS )

!  4/26-29/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Linearized Output developed for Column Jacobians.

      DOUBLE PRECISION, INTENT (INOUT) :: LC_ALBMED_USER  ( MAXSTOKES, MAX_USER_STREAMS, MAX_ATMOSWFS ) !  TOA illumination
      DOUBLE PRECISION, INTENT (INOUT) :: LC_TRNMED_USER  ( MAXSTOKES, MAX_USER_STREAMS, MAX_ATMOSWFS ) !  BOA illumination
      
      DOUBLE PRECISION, INTENT (INOUT) :: LC_ALBMED_FLUXES ( MAXSTOKES, 2, MAX_ATMOSWFS )    !  TOA illumination
      DOUBLE PRECISION, INTENT (INOUT) :: LC_TRNMED_FLUXES ( MAXSTOKES, 2, MAX_ATMOSWFS )    !  BOA illumination
      
      DOUBLE PRECISION, INTENT (INOUT) :: LC_TRANSBEAM ( MAXSTOKES, MAXBEAMS, MAX_ATMOSWFS ) !  Planetary problem

!  Disabled

      !DOUBLE PRECISION, INTENT (INOUT) :: TAUGRID ( 0:MAXLAYERS )
      !DOUBLE PRECISION, INTENT (INOUT) :: OMEGA_TOTAL ( MAXLAYERS )
      !DOUBLE PRECISION, INTENT (INOUT) :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )

      !DOUBLE PRECISION, INTENT (INOUT) :: STOKES_SS   ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      !DOUBLE PRECISION, INTENT (INOUT) :: STOKES_DB   ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      !DOUBLE PRECISION, INTENT (INOUT) :: SS_CONTRIBS ( MAX_GEOMETRIES, MAXSTOKES, MAXLAYERS )

      !DOUBLE PRECISION, INTENT (INOUT) :: COLUMNWF_SS  ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS )
      !DOUBLE PRECISION, INTENT (INOUT) :: COLUMNWF_DB  ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )
      !DOUBLE PRECISION, INTENT (INOUT) :: SURFACEWF_DB ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES )

!  2p6 Variables replaced 
!      DOUBLE PRECISION, INTENT (INOUT) ::  LTE_ATMOSWF &
!          ( 0:MAXLAYERS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_DIRECTIONS )

!  OUTPUTS, Intent(OUT)
!  ====================

!  Inclusion flags

      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFACE
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_THERMEMISS

!  Fourier components

      DOUBLE PRECISION, INTENT (OUT) ::  MS_CONTRIBS_F ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (OUT) ::  STOKES_F ( MAX_USER_LEVELS, MAX_USER_VZANGLES, &
                                                    MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

! 1/31/21. Version 2.8.3. . DO_MSSTS option final installation.
!    -- Additional layer_mssts and surf_mssts, Fourier component output (upwelling case)
!    -- MSST situation expanded to include Downwelling case (as an alternative, not both!)

      DOUBLE PRECISION, INTENT (OUT) :: LAYER_MSSTS_F  ( MAX_SZANGLES, MAXSTOKES, MAXLAYERS  )
      DOUBLE PRECISION, INTENT (OUT) :: SURF_MSSTS_F   ( MAX_SZANGLES, MAXSTOKES  )

!  Jacobian Fourier components

      DOUBLE PRECISION, INTENT (OUT) ::  COLUMNWF_F &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )
      DOUBLE PRECISION, INTENT (OUT) ::  SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_VZANGLES, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS )

!  1/31/21. Version 2.8.3. Installed MSST linearizations
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)

      DOUBLE PRECISION, INTENT (OUT) :: LC_SURF_MSSTS_F  ( MAXBEAMS, MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LC_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (OUT) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAXSTOKES, MAX_SURFACEWFS )
      DOUBLE PRECISION, INTENT (OUT) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXSTOKES, MAXLAYERS, MAX_SURFACEWFS )

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  Local variables
!  ================

!  From VLIDORT_MISCSETUPS
!  -----------------------

!  deltam-scaled optical
!mick fix 9/19/2017 - added LEVELS_SOLARTRANS & PARTIALS_SOLARTRANS to facilitate correction of direct flux

      DOUBLE PRECISION :: DELTAU_VERT ( MAXLAYERS )
      DOUBLE PRECISION :: PARTAU_VERT ( MAX_PARTLAYERS )
      DOUBLE PRECISION :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Transmittances

      DOUBLE PRECISION :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_UTUP_DISORDS ( MAXSTREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UTDN_DISORDS ( MAXSTREAMS, MAX_PARTLAYERS )

      DOUBLE PRECISION :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      DOUBLE PRECISION :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

      DOUBLE PRECISION :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      DOUBLE PRECISION :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )

!  Beam pseudo-spherical

      INTEGER ::          BEAM_CUTOFF ( MAXBEAMS )
      DOUBLE PRECISION :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

      DOUBLE PRECISION :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: LOCAL_CSZA ( 0:MAXLAYERS, MAXBEAMS )

!  Disabled

      !DOUBLE PRECISION :: TRUNC_FACTOR ( MAXLAYERS )
      !DOUBLE PRECISION :: FAC1 ( MAXLAYERS )
      !DOUBLE PRECISION :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      !DOUBLE PRECISION :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Need the following variables now (formerly disabled)

      DOUBLE PRECISION :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  From THERMAL_SETUP
!  ------------------

!  Bookkeeping

      DOUBLE PRECISION :: THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS )
      DOUBLE PRECISION :: TCOM1        ( MAXLAYERS, MAX_THERMAL_COEFFS )

!  Driect thermal solutions

      DOUBLE PRECISION :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_UP  ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UT_DIRECT_DN  ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  From EMULT_MASTER
!  -----------------

!  Multipliers

      DOUBLE PRECISION :: EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )
      DOUBLE PRECISION :: UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS )

!  1/31/21. Version 2.8.3. Need the following variables now (formerly disabled)

      DOUBLE PRECISION :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      DOUBLE PRECISION :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Disabled

      !LOGICAL ::          EMULT_HOPRULE ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      !DOUBLE PRECISION :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      !DOUBLE PRECISION :: SIGMA_P  ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  FROM VLIDORT_L_MISCSETUPS
!  -------------------------
!mick fix 9/19/2017 - added LC_LEVELS_SOLARTRANS & LC_PARTIALS_SOLARTRANS to facilitate
!                     correction of linearized direct flux

      DOUBLE PRECISION :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      DOUBLE PRECISION :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTDN_DISORDS ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTUP_DISORDS ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      DOUBLE PRECISION :: LC_LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Disabled

      !DOUBLE PRECISION :: L_OMEGA_TOTAL ( MAX_ATMOSWFS, MAXLAYERS )
      !DOUBLE PRECISION :: L_GREEKMAT_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
      !DOUBLE PRECISION :: L_DELTAU_SLANT   ( MAX_ATMOSWFS, MAXLAYERS, MAXLAYERS, MAXBEAMS )
      !LOGICAL ::          DO_SCATMAT_VARIATION  ( MAXLAYERS, MAX_ATMOSWFS )
      !DOUBLE PRECISION :: L_TRUNC_FACTOR  ( MAX_ATMOSWFS, MAXLAYERS )

!  FROM THERMAL_SETUP_PLUS
!  -----------------------

      DOUBLE PRECISION :: L_THERMCOEFFS  ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_DELTAU_POWER ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_XTAU_POWER   ( MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_TCOM1        ( MAXLAYERS, MAX_THERMAL_COEFFS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  FROM L_EMULT_MASTER
!  -------------------

      DOUBLE PRECISION :: LC_EMULT_UP &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_EMULT_DN &
          ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_UT_EMULT_UP &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_UT_EMULT_DN &
          ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  SOLUTION VARIABLES
!  ==================

!  At quadrature angles

      DOUBLE PRECISION :: PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQM ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION :: PI_XQM_POST ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQM_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQP_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  At user angles

      DOUBLE PRECISION :: PI_XUP ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUM ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION :: PI_XUM_POST ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUP_PRE  ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

!  At solar angles

      DOUBLE PRECISION :: PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  Thermal solutions of the RTE
!  ----------------------------

!  Discrete ordinate solutions

      DOUBLE PRECISION :: T_WUPPER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: T_WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      DOUBLE PRECISION :: UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS )

!  User-solutions

      DOUBLE PRECISION :: U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS )

!  Post processed solutions

      DOUBLE PRECISION :: LAYER_TSUP_UP ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: LAYER_TSUP_DN ( MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Homogeneous RTE
!  ---------------

!  Eigenmatrices

      DOUBLE PRECISION :: SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )

      DOUBLE PRECISION :: EIGENMAT_SAVE ( MAXEVALUES, MAXEVALUES, MAXLAYERS )

!  Eigensolutions
!    -- 1/31/21. Version 2.8.3. KEIGEN_CSQ now has extra MAXLAYERS dimension

      DOUBLE PRECISION :: REAL_KSQ   ( MAXSTRMSTKS )
      DOUBLE PRECISION :: IMAG_KSQ   ( MAXSTRMSTKS )
      DOUBLE PRECISION :: KEIGEN_CSQ ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: LEFT_EVEC  ( MAXSTRMSTKS, MAXSTRMSTKS )
      DOUBLE PRECISION :: RITE_EVEC  ( MAXSTRMSTKS, MAXSTRMSTKS )

      LOGICAL ::          EIGENDEGEN  ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER ::          EIGENMASK_R ( MAXEVALUES )
      INTEGER ::          EIGENMASK_C ( MAXEVALUES )

!  auxiliary solutions

      DOUBLE PRECISION :: FWD_SUMVEC  ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: FWD_DIFVEC  ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

!  Main solution variables

      INTEGER ::          K_REAL     ( MAXLAYERS )
      INTEGER ::          K_COMPLEX  ( MAXLAYERS )
      DOUBLE PRECISION :: KEIGEN     ( MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION :: SOLA_XPOS  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: SOLB_XNEG  ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  1/31/21. Version 2.8.3. Saved quantities for the Green function solution (Normalizations)
!    --INDICES REVERSED FROM LIDORT values

      DOUBLE PRECISION :: NORM_SAVED ( MAXEVALUES, MAXLAYERS)

!  Eigenstream transmittances

      DOUBLE PRECISION :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )
      DOUBLE PRECISION :: T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS )

!  User homogeneous solutions

      DOUBLE PRECISION :: HELPSTOKES ( 0:MAXMOMENTS, MAXEVALUES, MAXSTOKES )
      DOUBLE PRECISION :: ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )

!  Homogeneous solution multipliers

      LOGICAL ::          HSINGO  ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
      DOUBLE PRECISION :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )

      DOUBLE PRECISION :: UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )
      DOUBLE PRECISION :: UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Solar beam particular integral
!  ------------------------------

!  associated arrays for solving 

      DOUBLE PRECISION :: QSUMVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: QDIFVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: QVEC_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION :: QDIF_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION :: QMAT_SAVE ( MAXSTRMSTKS, MAXSTRMSTKS )
      INTEGER ::          QPIVOT ( MAXSTRMSTKS )

!  Quadrature Solutions

      DOUBLE PRECISION :: BVEC   ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  User solutions

      DOUBLE PRECISION :: HELPSTOKES_BEAM ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: UPAR_DN_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  BVProblem arrays
!  ----------------

!  Regular BVP

      DOUBLE PRECISION :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER ::          IPIVOT   ( MAXTOTAL )

!  One-layer BVP

      DOUBLE PRECISION :: SMAT2    ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER ::          SIPIVOT  ( MAXSTRMSTKS_2 )

!  Column vectors [mick fix 7/29/2014 - added to make VLIDORT threadsafe]

      DOUBLE PRECISION :: COL2    ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION :: SCOL2   ( MAXSTRMSTKS_2, MAXBEAMS )

!  Integration constants (BVP results)

      DOUBLE PRECISION :: LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION :: MCON ( MAXSTRMSTKS, MAXLAYERS )

!  Telescoped BVP setup

      INTEGER ::          ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER ::          N_BVTELMATRIX_SIZE
      INTEGER ::          N_BVTELMATRIX_SUPDIAG
      INTEGER ::          N_BVTELMATRIX_SUBDIAG
      INTEGER ::          NLAYERS_TEL

!  Telescoped BVP

      DOUBLE PRECISION :: BANDTELMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER ::          IPIVOTTEL   ( MAXTOTAL )
      DOUBLE PRECISION :: COLTEL2     ( MAXTOTAL, MAXBEAMS ) ! [mick fix 7/29/2014 - added to make VLIDORT threadsafe]

!  Surface-reflected solutions
!  ---------------------------

!  Thermal BOA term

      DOUBLE PRECISION :: BOA_THTONLY_SOURCE ( MAXSTREAMS, MAXSTOKES )

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving. These arrays NOT ENABLED as of 4/9/19.
     
      DOUBLE PRECISION :: LS_TRANS_ATMOS_FINAL   ( MAXBEAMS, MAX_SURFACEWFS )
      DOUBLE PRECISION :: LSSL_TRANS_ATMOS_FINAL ( MAXBEAMS, MAX_SLEAVEWFS  )
      
!  Direct beam solutions, 4/9/19 renamed

      DOUBLE PRECISION :: ATMOS_ATTN ( MAXBEAMS )
      DOUBLE PRECISION :: RF_DIRECT_BEAM      ( MAXSTREAMS,       MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION :: RF_USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
!      DOUBLE PRECISION :: DIRECT_BEAM      ( MAXSTREAMS, MAXBEAMS, MAXSTOKES )
!      DOUBLE PRECISION :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!  4/9/19 Surface-leaving contributions, added
      
      DOUBLE PRECISION :: SL_QUADTERM ( MAXSTREAMS,       MAXBEAMS, MAXSTOKES )
      DOUBLE PRECISION :: SL_USERTERM ( MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )
      
!  BRDF auxiliary array

      DOUBLE PRECISION :: AXBID_F  ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

!  reflected downwelling solutions

      DOUBLE PRECISION :: STOKES_DOWNSURF ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION :: R2_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: R2_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: R2_BEAM ( MAXSTREAMS, MAXSTOKES )

!  Cumulative discrete-ordinate transmittances to surface (Telescoped problem only)

      DOUBLE PRECISION :: CUMTRANSDOM(MAXSTREAMS)
      DOUBLE PRECISION :: CUMQUADDOM (MAXSTREAMS)

!  Cumulative source terms
!  -----------------------

!  Post-processed

      DOUBLE PRECISION :: CUMSOURCE_UP ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )
      DOUBLE PRECISION :: CUMSOURCE_DN ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )

!  Linearized Thermal Solutions
!  ----------------------------

      DOUBLE PRECISION :: L_T_WUPPER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_WLOWER ( MAXSTREAMS_2, MAXLAYERS, MAX_ATMOSWFS)

      DOUBLE PRECISION :: L_UT_T_PARTIC ( MAXSTREAMS_2, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_U_TPOS1 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_U_TNEG1 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_U_TPOS2 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_U_TNEG2 ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_LAYER_TSUP_UP ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_TSUP_DN ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized user solutions

      DOUBLE PRECISION :: L_UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_BOA_THTONLY_SOURCE (MAXSTREAMS, MAX_ATMOSWFS )

!  Linearized Multipliers

      DOUBLE PRECISION :: L_HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_UT_HMULT_UU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UT_HMULT_UD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UT_HMULT_DU ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UT_HMULT_DD ( MAXEVALUES, MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_ZETA_M ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_ZETA_P ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS,  MAX_ATMOSWFS )

!  Linearized eigenproblem

      DOUBLE PRECISION :: L_SAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_DAB &
          ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_EIGENMAT ( MAXEVALUES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_KEIGEN   ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTUP_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_T_UTDN_EIGEN ( MAXEVALUES, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Particular integral

      DOUBLE PRECISION :: LC_BVEC ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_UPAR_DN_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UPAR_UP_1  ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_UPAR_DN_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  linearized integration constants

      DOUBLE PRECISION :: NCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION :: PCON_SWF ( MAX_SURFACEWFS, MAXSTRMSTKS, MAXLAYERS )

      DOUBLE PRECISION :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!      DOUBLE PRECISION :: LSSL_DIRECT_BEAM      ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS, MAXSTOKES )
!      DOUBLE PRECISION :: LSSL_USER_DIRECT_BEAM ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@  NEW SECTION  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  1/31/21. Version 2.8.3. Introduce Green's function holding arrays
!  -----------------------------------------------------------------

!  Output from Green function solution, for use in the linearizations

      DOUBLE PRECISION :: DMI (MAXSTREAMS,MAXSTOKES), DPI (MAXSTREAMS,MAXSTOKES)

!  Saved quantities for the Green function solution

      DOUBLE PRECISION :: ATERM_SAVE (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION :: BTERM_SAVE (MAXEVALUES,MAXLAYERS)

!   Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      DOUBLE PRECISION :: L_ATERM_SAVE ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_BTERM_SAVE ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_NORM_SAVED ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Layer C and D functions, Multipliers GFUNC

      DOUBLE PRECISION :: CFUNC (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION :: DFUNC (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION :: GFUNC_UP (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION :: GFUNC_DN (MAXEVALUES,MAXLAYERS)

!  1/31/21. Version 2.8.3. Green's function multipliers for off-grid optical depths

      DOUBLE PRECISION :: UT_GMULT_UP(MAXEVALUES,MAX_PARTLAYERS)
      DOUBLE PRECISION :: UT_GMULT_DN(MAXEVALUES,MAX_PARTLAYERS)

!  Holding arrays for multiplier coefficients

      DOUBLE PRECISION :: GAMMA_M (MAXEVALUES,MAXLAYERS)
      DOUBLE PRECISION :: GAMMA_P (MAXEVALUES,MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      DOUBLE PRECISION :: PMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION :: PMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION :: PMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)
      DOUBLE PRECISION :: PMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAXLAYERS)

!  Post-processed Greens function multipliers (Partial layer)

      DOUBLE PRECISION :: UT_PMULT_DU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION :: UT_PMULT_DD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION :: UT_PMULT_UU(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)
      DOUBLE PRECISION :: UT_PMULT_UD(MAXEVALUES,MAX_USER_STREAMS,MAX_PARTLAYERS)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@   END NEW SECTION   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Other local variables
!  ---------------------

!mick mod 9/19/2017 - added and activated the variable DO_LOCALBEAM as in LIDORT

!  Local direct beam reflectance

      LOGICAL ::          DO_LOCALBEAM ( MAXBEAMS )

!  Indices and help for the Planetary problem

      INTEGER ::          I, IBEAM, IPARTIC, LAYER, N, O1, K, K0, K1, K2, KO1, LUM, UM, Q
      DOUBLE PRECISION :: TRANSQUAD(MAXSTREAMS,MAXSTOKES), TRANSDIRECT
      DOUBLE PRECISION :: L_TRANSQUAD(MAXSTREAMS,MAXSTOKES,MAX_ATMOSWFS), L_TRANS_DIFF, L_TRANS_DIRC
      DOUBLE PRECISION :: SHOM_R, SHOM_CR, SHOM, LXR, MXR, LXR_CR, LXR_CI, HOM1CR,  MXR_CR
      DOUBLE PRECISION :: HOM1, HOM2, HOM2CR, HOM3CR
      DOUBLE PRECISION :: NXR, PXR, NXR1, NXR2, PXR1
      DOUBLE PRECISION :: LXR1, LXR2, LLXR, MLXR, LLXR1, MLXR1, LLXR2

!  Local inclusion flags
!   Version 2.8.1, Control for TOA/BOA  illumination added, 3/23/19
!  5/24/21. Version 2.8.3. Add Include SURFACEWF flag

      LOGICAL ::          DO_INCLUDE_MVOUTPUT
      LOGICAL ::          DO_INCLUDE_DIRECTRF   ! 4/9/19 renamed
      LOGICAL ::          DO_INCLUDE_DIRECTSL   ! 4/9/19 New
!      LOGICAL ::          DO_INCLUDE_DIRECTBEAM  ! Replaced 4/9/19
      LOGICAL ::          DO_INCLUDE_TOAFLUX
      LOGICAL ::          DO_INCLUDE_BOAFLUX
      LOGICAL ::          DO_INCLUDE_SURFACEWF

!  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION :: FLUX_MULT
      DOUBLE PRECISION :: DELTA_FACTOR
      DOUBLE PRECISION :: SURFACE_FACTOR, SL, TFACTOR, RATIO

!  Weighting function indices

      LOGICAL ::          DO_RTSOL_VARY
      INTEGER ::          NPARAMS_VARY , VARIATION_INDEX

!  Error tracing variables

      CHARACTER (LEN=2) :: CF
      INTEGER ::           STATUS_SUB

!TESTING

!      INTEGER :: K,Q,UM,UTA,V
!      LOGICAL :: DO_DEBUG=.FALSE.
!      LOGICAL ::           FAIL

!  Progress

      LOGICAL, PARAMETER :: DO_WRITE_SCREEN = .FALSE.

!  1/31/21. Version 2.8.3. (RTS 2/16/21). Local NSTOKES and related control integers

      INTEGER :: LC_NSTK, LC_NSTKNSTRM, LC_NSTKNSTRM_2, LC_NTT, LC_NSUBD, LC_NSUPD

!  ##############
!  initialization
!  ##############

!mick debug
!write(*,*)
!write(*,*) 'in VLIDORT_LCS_FOURIER'

!mick debug
!write(*,*) 'at 16'
!call system_mem_usage3(valueRSS)

!  debug

      IF ( DO_DEBUG_WRITE ) THEN
        write(*,*) ; write(*,*) 'FOURIER = ' ,FOURIER
      END IF

!  module status and message initialization

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      TRACE_3 = ' '

!  1/31/21. Verson 2.8.3. This has been removed now.
!      LOCAL_UM_START = 1

!  Set local flags
!  ---------------

!  Local simulation only
!   Definition extended, 28 March 2014

      DO_SIMULATION_ONLY = ( .NOT. DO_COLUMN_LINEARIZATION .AND. .NOT. DO_SURFACE_LINEARIZATION .AND. &
                             .NOT. do_ATMOS_LBBF            .AND. .NOT. do_surface_LBBF )

!  inclusion of thermal surface emission term, only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE.
      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_SURFEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of thermal emission term, only for Fourier = 0

      DO_INCLUDE_THERMEMISS = .FALSE.
      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_THERMEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of TOA  illumination, only for Fourier = 0
!     Version 2.8.1, Control for TOA  illumination added, 3/23/19

      DO_INCLUDE_TOAFLUX = .FALSE.
      IF ( DO_TOAFLUX ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_TOAFLUX = .TRUE.
        ENDIF
      ENDIF

!  inclusion of BOA  illumination, only for Fourier = 0
!     Version 2.8.1, Control for BOA  illumination added, 3/23/19

      DO_INCLUDE_BOAFLUX = .FALSE.
      IF ( DO_BOAFLUX ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_BOAFLUX = .TRUE.
        ENDIF
      ENDIF

!  Surface flag (for inclusion of some kind of reflecting boundary)
!  9/9/19. Restore original: "DO_INCLUDE_SURFACE = .TRUE." needed for surfacewf even when ALBEDO = 0
!  5/24/21. Version 2.8.3. Add Include SURFACEWFS flag

      DO_INCLUDE_SURFACEWF = .TRUE.
      DO_INCLUDE_SURFACE   = .TRUE.
      IF ( DO_LAMBERTIAN_SURFACE ) THEN
         IF ( FOURIER .NE. 0 ) THEN
            DO_INCLUDE_SURFACE   = .FALSE.
            DO_INCLUDE_SURFACEWF = .FALSE.
!        ELSE
! ==> REMOVE THIS CONDITION, CLEARLY WRONG. Zero albedo case relaxed, if surface leaving. 12/17/15
!          IF ( ALBEDO .EQ. ZERO .AND..NOT.DO_SURFACE_LEAVING ) THEN
!            DO_INCLUDE_SURFACE = .FALSE.
!          ENDIF
         ENDIF
      ENDIF

!  5/24/21. Version 2.8.3. Add Include SLEAVEWFS flag
!   -- This is an output. Only True with EXTERNAL_WLEAVE and water-leaving, and fluorescence

      DO_INCLUDE_SLEAVEWFS = .FALSE.
      IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS) THEN
        IF ( DO_WATER_LEAVING ) THEN
          IF ( DO_EXTERNAL_WLEAVE ) then
            IF  ( ( DO_SL_ISOTROPIC .and. FOURIER.eq.0 ) .or..not.DO_SL_ISOTROPIC ) DO_INCLUDE_SLEAVEWFS = .true.
          ENDIF
        ELSE
          IF ( FOURIER.eq.0 ) DO_INCLUDE_SLEAVEWFS = .true.
        ENDIF
      ENDIF

!  Surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR = TWO
        DELTA_FACTOR   = ONE
      ELSE
        SURFACE_FACTOR = ONE
        DELTA_FACTOR   = TWO
      ENDIF

!  Flux multipliers

!       = 1 / 4.pi with beam sources,  = 1 for Thermal alone.    ! OLD
!       = 1                                                      ! NEW

!      FLUX_MULT   = SS_FLUX_MULT * DELTA_FACTOR     ! OLD
      FLUX_MULT = DELTA_FACTOR                             ! NEW

      IF ( DO_INCLUDE_THERMEMISS.AND..NOT.DO_SOLAR_SOURCES) THEN
        FLUX_MULT = DELTA_FACTOR
      ENDIF

!  Inclusion of mean value output

      DO_INCLUDE_MVOUTPUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_MVOUTPUT = .TRUE.
        ENDIF
      ENDIF

!  ##################
!  mick chg 7/17/2014 - moved misc setup, thermal setup, emult setup, sscor codes (2),
!                       db code, and fo code up to VLIDORT_LCS_MASTER.  There, the needed
!                       vars are packed into type structures to bypass f90 continuation
!                       line limits and passed down.  Here, they are unpacked.
!                     - (Later) Only the misc-setups, thermal-setups, emult setups are UNPACKED HERE
!  (Rev. 7/22/16, 2.8) - Cleaned up

!  1/31/21. Version 2.8.3. Need additional variable ITRANS_USERM to be unpacked.
!    -- Argument list rearranged to follow more natural logic 

      CALL VLIDORT_UNPACK_MISC ( Misc,                             & ! Input structure to be unpacked
        NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS,                  & ! Input Numbers
        N_USER_STREAMS, NBEAMS, NMOMENTS,                          & ! Input Numbers
        DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT,                    & ! unpacked Optical
        LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, OMEGA_GREEK,       & ! unpacked Optical
        T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,            & ! unpacked Trans D.O.
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM, CUMTRANS,        & ! unpacked Trans User
        BEAM_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,    & ! unpacked Solar
        INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR, & ! unpacked Average Secant
        ITRANS_USERM, LOCAL_CSZA )                                   ! unpacked Misc

      IF ( DO_THERMAL_EMISSION ) THEN
        CALL VLIDORT_UNPACK_THERM ( Therm,                          & !Input Sstructure to be unpacked
          NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,  & ! Input Numbers
          THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,             & ! Output
          T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN )  ! Output
      END IF

!  Flag changed 3/1/17.
!        IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN

!  1/31/21. Version 2.8.3. Need additional variables SIGMA_M, SIGMA_P to be unpacked.

      IF ( DO_SOLAR_SOURCES .AND. DO_USER_STREAMS ) THEN
        IF ( .NOT. DO_FOCORR_ALONE ) THEN
          CALL VLIDORT_UNPACK_MULT ( Mult,                  & ! Input structure to be unpacked
            NLAYERS, N_PARTLAYERS, NBEAMS, N_USER_STREAMS,  & ! Input numbers
            SIGMA_M, SIGMA_P, EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )    ! Output
        ENDIF
      ENDIF

!  Linearized

      IF ( DO_ATMOS_LINEARIZATION ) THEN

        CALL VLIDORT_UNPACK_LAC_MISC ( LAC_Misc,                & ! Input
          NSTOKES, NSTREAMS, NLAYERS, N_PARTLAYERS,             & ! Input
          NBEAMS, N_USER_STREAMS, NMOMENTS, N_TOTALCOLUMN_WFS,  & ! Input
          L_DELTAU_VERT, L_OMEGA_GREEK,                         & ! Output
          L_T_DELT_DISORDS, L_T_UTDN_DISORDS, L_T_UTUP_DISORDS, & ! Output
          L_T_DELT_USERM,   L_T_UTDN_USERM,   L_T_UTUP_USERM,   & ! Output
          LC_AVERAGE_SECANT, LC_INITIAL_TRANS,                  & ! Output
          LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,                     & ! Output
          LC_LEVELS_SOLARTRANS, LC_PARTIALS_SOLARTRANS )          ! Output

        IF ( DO_THERMAL_EMISSION ) THEN
          CALL VLIDORT_UNPACK_L_THERM ( L_Therm,                                        & ! Input
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS, N_TOTALCOLUMN_WFS, & ! Input
            L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_TCOM1,                       & ! Output
            L_T_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_UP, L_T_UT_DIRECT_DN )            ! Output
        END IF

!  Flag changed 3/1/17.
!        IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN

        IF ( DO_SOLAR_SOURCES .AND. DO_USER_STREAMS ) THEN
          IF ( .NOT. DO_FOCORR_ALONE ) THEN
            CALL VLIDORT_UNPACK_LC_MULT ( LC_Mult,                     & ! Input
              NLAYERS, N_PARTLAYERS, NBEAMS, N_USER_STREAMS,           & ! Input
              N_TOTALCOLUMN_WFS,                                       & ! Input
              LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN ) ! Output
          ENDIF
        ENDIF

      ENDIF

!  ##################

!  Direct beam flag (only if above albedo flag has been set)
!mick fix 7/11/2014 - this IF block moved from above and modified to accomodate move of
!                     "Fourier = 0" subroutines
!mick mod 9/19/2017 - commented out this IF block (DO_REFLECTED_DIRECTBEAM already defined
!                     in VLIDORT_LCS_MASTER & passed down VIA the "Misc" type structure)
!                   - activated DO_LOCALBEAM as in LIDORT

      !IF ( FOURIER .GT. 0 ) THEN
      !  IF ( DO_DIRECT_BEAM .AND. DO_SOLAR_SOURCES ) THEN
      !    IF ( DO_INCLUDE_SURFACE ) THEN
      !      DO IBEAM = 1, NBEAMS
      !        DO_REFLECTED_DIRECTBEAM(IBEAM) = .TRUE.
      !      ENDDO
      !    ELSE
      !      DO IBEAM = 1, NBEAMS
      !        DO_REFLECTED_DIRECTBEAM(IBEAM) = .FALSE.
      !      ENDDO
      !    ENDIF
      !  ENDIF
      !ENDIF

      IF ( DO_SOLAR_SOURCES .AND. DO_INCLUDE_SURFACE ) THEN
        DO_LOCALBEAM(1:NBEAMS) = DO_REFLECTED_DIRECTBEAM(1:NBEAMS)
      ELSE
        DO_LOCALBEAM(1:NBEAMS) = .FALSE.
      ENDIF

!  Surface direct beam (Not required if no solar sources)
!    - New Surface-Leaving arguments added, 17 May 2012
!    - Revised I/O listings for Version 2.8, 7/8/16
!mick mod 9/19/2017 - replaced DO_REFLECTED_DIRECTBEAM with DO_LOCALBEAM

!  4/9/19. Argument list refined to include Water-leaving control
!  4/9/19. Here, SL output ONLY for non water-leaving, or water-leaving external, otherwise zeroed

!  1/31/21. Version 2.8.3. BRDF/SLEAVE input arrays are defined locally for each Fourier.
!  1/31/21. Version 2.8.3. Use PostProcessing mask input, instead of N_USER_STREAMS, LOCAL_UM_START
!  5/24/21. Version 2.8.3. Add the TF_iteration flag
!  5/24/21. Version 2.8.3. Must introduce TRANS_ATMOS_FINAL as argument into this routine

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( .NOT.DO_FOCORR_ALONE ) THEN
          CALL VLIDORT_DIRECTRADIANCE ( DO_USER_STREAMS, DO_REFRACTIVE_GEOMETRY,           & ! Input flags (general)
                DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_SURFACE_LEAVING,             & ! Input flags (surface)
                DO_SL_ISOTROPIC, DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,    & ! input
                NSTOKES, NSTREAMS, NBEAMS, NLAYERS, N_PPSTREAMS, PPSTREAM_MASK,            & ! input
                FOURIER, MUELLER_INDEX, FLUX_FACTOR, FLUXVEC, DELTA_FACTOR,                & ! input
                SZA_LOCAL_INPUT, COS_SZANGLES, TRANS_SOLAR_BEAM, DO_LOCALBEAM,             & ! Input solar beam
                ALBEDO, BRDF_F_0, USER_BRDF_F_0,                                           & ! Input surface reflectance
                TRANS_ATMOS_FINAL, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,          & ! input surface leaving
                ATMOS_ATTN, RF_DIRECT_BEAM, RF_USER_DIRECT_BEAM, SL_QUADTERM, SL_USERTERM )  ! Output reflected direct beams.
        ENDIF

!  Useful debug
!      do IBEAM = 1, nbeams
!         if ( FOURIER.EQ.0 ) write(*,*)IBEAM,ATMOS_ATTN(IBEAM),RF_DIRECT_BEAM(1,IBEAM,1), RF_USER_DIRECT_BEAM(1,IBEAM,1)
!         if ( FOURIER.EQ.0 ) write(*,*)IBEAM,SL_QUADTERM(1,IBEAM,1), SL_USERTERM(1,IBEAM,1)
!      enddo

!  Addition of SLEAVE weighting function setups, R. Spurr, 22 August 2012
!  4/9/19. Only done here for non water-leaving, or water-leaving external, otherwise zeroed
!  4/9/19. No point in having this routine, can do what you need in the one master routine.         
!        IF ( DO_SURFACE_LEAVING .AND. DO_SLEAVE_WFS .AND.FOURIER .EQ. 0 ) then
!         CALL VLIDORT_LSSL_DBSETUPS ( &
!           DO_OBSERVATION_GEOMETRY, DO_SL_ISOTROPIC, DO_USER_STREAMS,    &
!           DO_REFLECTED_DIRECTBEAM, FOURIER,                             &
!           NSTOKES, NSTREAMS, NBEAMS, N_USER_STREAMS, N_SLEAVE_WFS,      &
!           FLUX_FACTOR, DELTA_FACTOR,                                    &
!           LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0, &
!           LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM )
!        ENDIF

!  End solar sources only clause for corrections

      ENDIF

!  ###################
!  Spherical functions
!  ###################

!  Get Pi Matrices for this Fourier component, Including all beam angles !!
!   - Version 2.7. Use "OMP" routine for thread-safety. Mick fix.
!   - Version 2.8. Clean up arguments 7/8/16

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        CALL VLIDORT_PIMATRIX_SETUP_OMP ( &
          DO_REFRACTIVE_GEOMETRY, DO_USER_STREAMS, FOURIER,             & ! Input flags and Fourier
          NSTOKES, NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS,           & ! Input numbers
          QUAD_STREAMS, USER_STREAMS, COS_SZANGLES, SUN_SZA_COSINES,    & ! Input streams and angles
          NMOMENTS, MUELLER_INDEX, DMAT,                                & ! auxiliary inputs
          PIMM_11, PIMM_KM, PI_XQP, PI_XQM, PI_XUP, PI_XUM, PI_X0P,     & ! Output Pi-Matrix stuff
          PI_XQM_POST, PI_XQM_PRE, PI_XQP_PRE, PI_XUM_POST, PI_XUP_PRE )  ! Output Pi-Matrix stuff
      ENDIF

!  ###################
!  LOCAL NSTOKES USAGE
!  ###################

!  1/31/21. Version 2.8.3. (RTS 2/16/21). 
!    -- Introduce control for Fast execution of Fourier-zero with LOCAL_NSTOKES = 2
!rob & mick fix 1/5/2021 - added Lambertian IF condition

      IF ( DO_Fourier0_NSTOKES2 .and. NSTOKES.eq.3 .and. FOURIER.eq.0 .and. DO_LAMBERTIAN_SURFACE ) THEN
        LC_NSTK = 2 ; LC_NSTKNSTRM = NSTREAMS * LC_NSTK ; LC_NSTKNSTRM_2 = LC_NSTKNSTRM * 2
        LC_NTT  = LC_NSTKNSTRM_2 * NLAYERS
        IF ( NLAYERS .EQ. 1 ) THEN
          LC_NSUBD = 2*LC_NSTKNSTRM - 1 ; LC_NSUPD = LC_NSUBD
        ELSE
          LC_NSUBD = 3*LC_NSTKNSTRM - 1 ; LC_NSUPD = LC_NSUBD
        ENDIF
      ELSE
        LC_NSTK = NSTOKES ; LC_NSTKNSTRM   = NSTKS_NSTRMS ; LC_NSTKNSTRM_2 = NSTKS_NSTRMS_2
        LC_NTT  = NTOTAL  ; LC_NSUBD       = N_SUBDIAG    ; LC_NSUPD       = N_SUPDIAG
      ENDIF

!  #########################################################
!  RT differential equation Eigensolutions + linearizations
!  #########################################################

!  Version 2.8, remove GOTO statements
!  Go to continuation point for thermal transmittance only
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 8899

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN

!  Start layer loop

        DO LAYER = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer
!   Version 2.8, Cleaned-up I/O presentation. 7/8/16.

!  1/31/21. Version 2.8.3. Add TOLERANCE variable.
!    -- (RTS 2/16/21). use local variables LC_NSTK, LC_NSTKNSTRM

          CALL VLIDORT_QHOM_SOLUTION ( &
            DO_SOLUTION_SAVING, LAYER, FOURIER, TOLERANCE,                        & ! Input Flag, Fourier/Layer
            LC_NSTK, NSTREAMS, N_USER_LEVELS, NMOMENTS, LC_NSTKNSTRM,             & ! Input numbers
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,         & ! Input partial control
            QUAD_STREAMS, QUAD_HALFWTS, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, & ! Quadrature and bookkeeping
            DELTAU_VERT, PARTAU_VERT, OMEGA_GREEK, PI_XQP, PI_XQM, PI_XQM_PRE,    & ! Input optical and PI matrices
            T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                       & ! Input transmittances (discrete Ords.)
            SAB, DAB, EIGENMAT_SAVE, REAL_KSQ, IMAG_KSQ,                          & ! Output Eigenproblem
            LEFT_EVEC, RITE_EVEC, EIGENDEGEN, EIGENMASK_R, EIGENMASK_C,           & ! Output Eigenproblem
            FWD_SUMVEC, FWD_DIFVEC, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,     & ! Output Homogeneous solutions
            K_REAL, K_COMPLEX, KEIGEN, KEIGEN_CSQ, SOLA_XPOS, SOLB_XNEG,          & ! Output Homogeneous solutions
            STATUS_SUB, MESSAGE, TRACE_1 )                                          ! Exception handling

!  .. error tracing

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER
            TRACE_2 = 'Called in VLIDORT_LCS_FOURIER, Fourier component '//CF
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

!  Get Post-processing ("user") solutions for this layer
!   Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    -- 1/31/21. Version 2.8.3.  LOCAL_UM_START (argument dropped)
!    -- 1/31/21. Version 2.8.3.  (RTS 2/16/21). use local variable LC_NSTK instead of NSTOKES

          IF  ( STERM_LAYERMASK_UP(LAYER) .OR. STERM_LAYERMASK_DN(LAYER) ) THEN
            IF ( DO_USER_STREAMS ) THEN
              CALL VLIDORT_UHOM_SOLUTION ( &
                DO_UPWELLING, DO_DNWELLING, LAYER, FOURIER, LC_NSTK,         & ! Input flags
                NSTREAMS, N_USER_STREAMS, NMOMENTS, DO_LAYER_SCATTERING,     & ! Input numbers
                QUAD_HALFWTS, USER_SECANTS, OMEGA_GREEK,                     & ! Input quadratures + optical
                PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE, PI_XUM_POST, PI_XUP_PRE, & ! Input PI Matrices
                K_REAL, K_COMPLEX, KEIGEN, KEIGEN_CSQ, SOLA_XPOS, SOLB_XNEG, & ! Input Homog. RTE solutions
                UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                  & ! Output user Homog
                HSINGO, HELPSTOKES, ZETA_M, ZETA_P )                           ! Output user Homog
            ENDIF
          ENDIF

!  Additional Solutions for Linearization
!  --------------------------------------

          IF ( DO_ATMOS_LINEARIZATION ) THEN

!  Parameter control

            IF ( DO_COLUMN_LINEARIZATION ) THEN
              DO_RTSOL_VARY = .TRUE.
              NPARAMS_VARY  = N_TOTALCOLUMN_WFS
            ENDIF

!  Linearizations of Discrete Ordinate solutions
!    -- 1/31/21. Version 2.8.3. KEIGEN_CSQ now has an extra "LAYERS" dimension
!    -- 1/31/21. Version 2.8.3.  (RTS 2/16/21). use local variable LC_NSTK instead of NSTOKES

!         if ( do_write_screen) write(*,*)'l_homsolution',layer
            CALL VLIDORT_L_QHOM_SOLUTION ( &
              DO_SOLUTION_SAVING, LAYER, FOURIER, LC_NSTK,                          & ! Input Flag and Indices
              NSTREAMS, N_PARTLAYERS, DO_RTSOL_VARY, NPARAMS_VARY, NMOMENTS,        & ! Input Numbers+Lin.Control
              NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2, PARTLAYERS_LAYERIDX,        & ! Input Numbers
              DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, QUAD_STREAMS, QUAD_HALFWTS, & ! Input bookkeeping and quadrature
              DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT, L_OMEGA_GREEK,               & ! Input Optical
              PI_XQP, PI_XQM_PRE, SAB, DAB, EIGENMAT_SAVE, REAL_KSQ, IMAG_KSQ,      & ! Input Homog solution
              LEFT_EVEC, RITE_EVEC, EIGENDEGEN, EIGENMASK_R, EIGENMASK_C,           & ! Input Homog solution
              K_REAL, K_COMPLEX, KEIGEN, KEIGEN_CSQ, FWD_SUMVEC, FWD_DIFVEC,        & ! Input Homog solution
              T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                             & ! Input Homog solution
              L_T_DELT_DISORDS, L_T_UTDN_DISORDS, L_T_UTUP_DISORDS,                 & ! Input linearized Homog Trans
              L_SAB, L_DAB, L_EIGENMAT, L_KEIGEN, L_SOLA_XPOS, L_SOLB_XNEG,         & ! Output lineaerized Homog. Sol.
              L_T_DELT_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,                       & ! Output linearized Homog Trans
              STATUS_SUB, MESSAGE, TRACE_1 )

!  .. error tracing

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_2 = 'Called in VLIDORT_LCS_FOURIER, Fourier component '//CF
              STATUS = VLIDORT_SERIOUS
              RETURN
            ENDIF

!  Get Linearizations of ("user") solutions for this layer
!    -- mick fix 3/30/2015 - added IF condition
!    -- 1/31/21. Version 2.8.3.  LOCAL_UM_START (argument dropped)
!    -- 1/31/21. Version 2.8.3.  (RTS 2/16/21). use local variable LC_NSTK instead of NSTOKES

            IF ( DO_USER_STREAMS ) THEN
              CALL VLIDORT_L_UHOM_SOLUTION ( &
                DO_UPWELLING, DO_DNWELLING, DO_DEBUG_WRITE,              & ! Input flags
                LAYER, FOURIER, LC_NSTK, NSTREAMS, N_USER_STREAMS,       & ! Input numbers
                DO_RTSOL_VARY, NPARAMS_VARY, NMOMENTS,                   & ! Input Lin-control, numbers
                QUAD_HALFWTS, USER_SECANTS, DO_LAYER_SCATTERING,         & ! Input streams/scat
                STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                  & ! Input bookkeeping
                OMEGA_GREEK, L_OMEGA_GREEK, PI_XQP, PI_XQM_PRE,          & ! Input optical + PI-quad
                PI_XUP, PI_XUM, PI_XUM_POST, PI_XUP_PRE,                 & ! Input PI-User
                K_REAL, K_COMPLEX, KEIGEN, HELPSTOKES,                   & ! Input eigensolution
                UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,              & ! Input User solutions
                ZETA_M, ZETA_P, L_KEIGEN, L_SOLA_XPOS, L_SOLB_XNEG,      & ! Input Zeta + linearized Homog
                L_UHOM_DNDN, L_UHOM_DNUP, L_UHOM_UPDN, L_UHOM_UPUP,      & ! Output Linearized user solutions
                L_ZETA_M, L_ZETA_P )                                       ! Output linearized Zeta
            ENDIF

!  End linearization control

          ENDIF

!  end layer loop

        ENDDO

!  1/31/21. Version 2.8.3. Greens function solution.
!    -- New subroutines --> Add Calculation of Normalization factors
!    -- 1/31/21. Version 2.8.3.  (RTS 2/16/21). use local variable LC_NSTK instead of NSTOKES

        IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
          CALL VLIDORT_QHOM_NORMS &
             ( FOURIER, LC_NSTK, NSTREAMS, NLAYERS, QUAD_STRMWTS, & ! Input numbers
               K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,           & ! Input solutions
               NORM_SAVED )                                         ! Output
          IF ( DO_ATMOS_LINEARIZATION ) THEN
            CALL VLIDORT_L_QHOM_NORMS &
              ( FOURIER, LC_NSTK, NSTREAMS, NLAYERS, QUAD_STRMWTS,      & ! Input
                LAYER_VARY_FLAG, LAYER_VARY_NUMBER, K_REAL, K_COMPLEX,  & ! Input
                SOLA_XPOS, SOLB_XNEG, L_SOLA_XPOS, L_SOLB_XNEG,         & ! Input
                L_NORM_SAVED )                                            ! Output
          ENDIF
        ENDIF

!  Prepare homogeneous solution multipliers
!    -- mick fix 3/30/2015 - added if DO_USER_STREAMS condition
!    -- Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    -- 1/31/21. Version 2.8.3.  LOCAL_UM_START (arguments dropped)

        IF ( DO_USER_STREAMS ) THEN

          CALL HMULT_MASTER ( &
            DO_UPWELLING, DO_DNWELLING,                                   & ! flags
            NLAYERS, N_USER_STREAMS, N_USER_LEVELS, TAYLOR_ORDER,         & ! Basic numbers
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                       & ! whole-layer control
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! partial-layer control
            USER_SECANTS, DELTAU_VERT, PARTAU_VERT,                       & ! secants, optical
            T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                     & ! User-stream transmittances
            T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                     & ! eigenstreamTransmittances
            K_REAL, K_COMPLEX, HSINGO, ZETA_M, ZETA_P,                    & ! RTE Eigen solution
            HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD ) ! Output Multipliers

          IF ( DO_COLUMN_LINEARIZATION ) THEN
            CALL L_HMULT_MASTER ( &
              DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER, NLAYERS, N_PARTLAYERS,    & ! Flags/Taylor/Numbers
              N_USER_STREAMS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                 & ! Numbers + Lin-control
              STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,  PARTLAYERS_LAYERIDX,       & ! Output control
              USER_SECANTS, DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,              & ! sstreams/Optical
              T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                           & ! Transmittances User
              T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                           & ! Transmittances Homog.
              K_REAL, K_COMPLEX, HSINGO, HMULT_1, HMULT_2, ZETA_M, ZETA_P,        & ! Multipliers and Zetas
              L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,                     & ! Linearized Transmittances User
              L_T_DELT_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,                     & ! LinearizedTransmittances Homog.
              L_KEIGEN, L_ZETA_M, L_ZETA_P,                                       & ! Linearized Zeta/eigenvalues
              L_HMULT_1, L_HMULT_2, L_UT_HMULT_UU,                                & ! Output - Linearized Multipliers 
              L_UT_HMULT_UD, L_UT_HMULT_DU, L_UT_HMULT_DD )                         ! Output - Linearized Multipliers 
          ENDIF

        ENDIF

!  ############################################
!   boundary value problem - MATRIX PREPARATION
!  ############################################

!  standard case using compression of band matrices, etc..
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    - 1/31/21. Version 2.8.3. BRDF_F argument defined locally for each Fourier.
!    - 1/31/21. Version 2.8.3. (RTS 2/16/21). Use local bookkeeping variables LC_NSTK etc....

        IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN

          CALL BVP_MATRIXSETUP_MASTER ( &
            DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,       & ! Input Flags
            FOURIER, LC_NSTK, NSTREAMS, NLAYERS,             & ! Input Numbers
            NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2,        & ! Input Numbers
            LC_NTT, LC_NSUBD, LC_NSUPD, SURFACE_FACTOR,      & ! Input Numbers
            QUAD_STRMWTS, ALBEDO, BRDF_F,                    & ! Input Surface stuff
            MUELLER_INDEX, K_REAL, K_COMPLEX,                & ! Input Bookkeeping
            SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,              & ! Input RTE stuff
            R2_HOMP, R2_HOMM, AXBID_F,                       & ! Output Surface reflection
            BANDMAT2, IPIVOT, SMAT2, SIPIVOT,                & ! Output BVP Matrices
            STATUS_SUB, MESSAGE, TRACE_1 )                     ! Exception handling

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER
            TRACE_2 = 'Error from BVP_MATRIXSETUP_MASTER, '// 'Called in VLIDORT_LPS_FOURIER, Fourier # '//CF
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

!  Telescoped case
!    - 1/31/21. Version 2.8.3. (RTS 2/16/21). Use local bookkeeping variables LC_NSTK etc....

        ELSE

          CALL BVPTEL_MATRIXSETUP_MASTER ( &
            DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,          & ! Input Flags
            DO_LAYER_SCATTERING, FOURIER, LC_NSTK, NSTREAMS,    & ! Input Numbers
            NLAYERS, NSTREAMS_2, LC_NSTKNSTRM_2, LC_NSTKNSTRM,  & ! Input Numbers
            DO_BVTEL_INITIAL, BVTEL_FOURIER,                    & ! BVP Tel Control
            N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,     & ! BVP Tel Control
            N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,       & ! BVP Tel Control
            SURFACE_FACTOR, MUELLER_INDEX, K_REAL, K_COMPLEX,   & ! Input Bookkeeping
            QUAD_STRMWTS, ALBEDO, BRDF_F,                       & ! Input Surface inputs
            SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_DELT_DISORDS, & ! Input RTE stuff
            R2_HOMP, R2_HOMM, AXBID_F, CUMTRANSDOM, CUMQUADDOM, & ! Output Surface reflection
            BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT,             & ! Output BVP Matrices
            STATUS_SUB, MESSAGE, TRACE_1 )                        ! Exception handling

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER
            TRACE_2 = 'Error from BVPTEL_MATRIXSETUP_MASTER, '//'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

        ENDIF

!  4/26/19. Add Call to Media properties subroutine (regular or LC_Jacobian)
!   -- Stand-alone output, but this is a necessary call for the planetary problem
!   -- 1/31/21. Version 2.8.3. (RTS 2/16/21). Use local bookkeeping variables LC_NSTK etc....

        IF ( DO_ALBTRN_MEDIA(1) .OR. DO_ALBTRN_MEDIA(2) .OR. DO_PLANETARY_PROBLEM ) THEN
           IF ( DO_COLUMN_LINEARIZATION ) THEN
              CALL VLIDORT_LC_MediaProps &
               ( DO_USER_STREAMS, DO_ALBTRN_MEDIA, DO_COLUMN_LINEARIZATION,          & ! Input flags
                 LC_NSTK, NLAYERS, NSTREAMS, N_USER_STREAMS, NSTREAMS_2,             & ! Input numbers
                 N_TOTALCOLUMN_WFS, LC_NSTKNSTRM, LC_NSTKNSTRM_2, LC_NTT,            & ! Input numbers
                 LC_NSUBD, LC_NSUPD, QUAD_STRMWTS, DELTAU_VERT,                      & ! Input quad/delta
                 K_REAL, K_COMPLEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG,              & ! Input Homog solutions
                 USER_STREAMS, T_DELT_USERM, UHOM_UPDN, UHOM_UPUP,                   & ! Input User solutions
                 HMULT_1, HMULT_2, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                 & ! Input Multipliers, BVP
                 L_DELTAU_VERT, L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG,            & ! Input Lin solutions
                 L_T_DELT_USERM, L_UHOM_UPDN, L_UHOM_UPUP, L_HMULT_1, L_HMULT_2,     & ! Input Lin solutions
                 ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,             & ! output Main
                 LC_ALBMED_USER, LC_ALBMED_FLUXES, LC_TRNMED_USER, LC_TRNMED_FLUXES, & ! output Linearized
                 STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                               ! Output
              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                 TRACE_3 = 'Error from VLIDORT_LC_MediaProps, '//'Called in VLIDORT_LCS_FOURIER'
                 STATUS = VLIDORT_SERIOUS ; RETURN
              ENDIF
           ELSE
              CALL VLIDORT_MediaProps &
               ( DO_USER_STREAMS, DO_ALBTRN_MEDIA, LC_NSTK,                       & ! Input
                 NLAYERS, NSTREAMS, N_USER_STREAMS, LC_NSTKNSTRM, LC_NSTKNSTRM_2, & ! Input
                 LC_NTT, LC_NSUBD, LC_NSUPD, QUAD_STRMWTS, DELTAU_VERT,           & ! Input
                 K_REAL, K_COMPLEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG,           & ! Input Homog solutions
                 USER_STREAMS, T_DELT_USERM, UHOM_UPDN, UHOM_UPUP,                & ! Input User solutions
                 HMULT_1, HMULT_2, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,              & ! Input Multipliers, BVP
                 ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,          & ! output
                 STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                            ! Output
              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                 TRACE_3 = 'Error from VLIDORT_MediaProps, '//'Called in VLIDORT_LCS_FOURIER'
                 STATUS = VLIDORT_SERIOUS ; RETURN
              ENDIF
           ENDIF
        ENDIF
        
!  End scattering calculation (replaces GOTO 8899)

      ENDIF

!  Continuation point for avoiding the scattering calculations
! 8899 continue. Removed, Version 2.8

!  ################
!  Thermal Solution (2 steps)
!  ################

!  THERMAL SOLUTIONS
!  =================

!  Separate calls if linearization is required

!  1. Find the Particular solution (also for transmittance only)
!  2. Compute thermal layer source terms. (Upwelling and Downwelling)
!    These will be scaled up by factor 4.pi if solar beams as well

!    ** REMARK. No Green's function treatment here
!    ** Version 2.8, Cleaned-up I/O presentations, all thermal routines. 7/8/16.

      IF ( DO_INCLUDE_THERMEMISS ) THEN

!  With linearization

        IF ( DO_ATMOS_LINEARIZATION ) THEN

!  discrete ordinate particular integral for thermal sources.
!    - mick fix 9/19/2017 - added DO_USER_STREAMS to input
!    - 1/31/21. Version 2.8.3. Drop LOCAL_UM_START argument. User-stream do-loops start with 1 now.

          CALL THERMAL_CLSOLUTION_PLUS ( &
            DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                 & ! input flags
            DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_THERMAL_TRANSONLY,    & ! Input flags
            DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,  & ! Input linearization control
            NSTREAMS, NLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,         & ! Input basic numbers
            NMOMENTS, NSTREAMS_2, N_PARTLAYERS,                          & ! Input numbers
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Input Level control
            QUAD_STREAMS, QUAD_HALFWTS, OMEGA_GREEK, SAB, DAB,           & ! Input optical and SAB/DAB
            T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,              & ! Input Discrete Ord. Trans.
            PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE,                          & ! Input PI matrices, 
            THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,                & ! Input thermal setups
            L_OMEGA_GREEK, L_SAB, L_DAB,                                 & ! Input Linearized optical and SAB/DAB
            L_T_DELT_DISORDS, L_T_UTDN_DISORDS, L_T_UTUP_DISORDS,        & ! Input Linearized Discrete Ord. Trans.
            L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER, L_TCOM1,        & ! Input Linearized thermal setups
            T_WUPPER, T_WLOWER, UT_T_PARTIC,                             & ! Output thermal solutions
            U_TPOS1, U_TNEG1, U_TPOS2, U_TNEG2,                          & ! Output User thermal solutions
            L_T_WUPPER, L_T_WLOWER, L_UT_T_PARTIC,                       & ! Output Linearized thermal solutions
            L_U_TPOS1, L_U_TNEG1, L_U_TPOS2, L_U_TNEG2,                  & ! Output Linearized User thermal solutions
            STATUS_SUB, MESSAGE, TRACE_1 )                                 ! Exception handling

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            TRACE_2 = 'Called in VLIDORT_LCS_FOURIER, Fourier 0'
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

!mick fix 3/30/2015 - modified if conditions for up & dn thermal source terms
!    - 1/31/21. Version 2.8.3. Drop LOCAL_UM_START argument. User-stream do-loops start with 1 now.

          !IF ( DO_UPWELLING ) THEN
          IF ( DO_UPWELLING .AND. DO_USER_STREAMS) THEN
            CALL THERMAL_STERMS_UP_PLUS ( &
              DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, DO_PARTLAYERS,                 & ! Input flags
              DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,            & ! Input Linearization control
              NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,               & ! Input numbers
              N_ALLLAYERS_UP, PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, USER_STREAMS, & ! Input level control, user-streams
              T_DELT_USERM, T_UTUP_USERM, L_T_DELT_USERM, L_T_UTUP_USERM,            & ! Input transmittances
              DELTAU_POWER, XTAU_POWER, L_DELTAU_POWER, L_XTAU_POWER,                & ! Input Thermal setups
              U_TPOS1, U_TPOS2, L_U_TPOS1, L_U_TPOS2,                                & ! Input Thermal solutions               
              T_DIRECT_UP, T_UT_DIRECT_UP, L_T_DIRECT_UP, L_T_UT_DIRECT_UP,          & ! Input thermal direct solutions
              LAYER_TSUP_UP, LAYER_TSUP_UTUP, L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP )     ! Output user RTE thermal
          ENDIF

          !IF ( DO_DNWELLING ) THEN
          IF ( DO_DNWELLING .AND. DO_USER_STREAMS ) THEN
            CALL THERMAL_STERMS_DN_PLUS ( &
              DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, DO_PARTLAYERS,                 & ! Input flags
              DO_ATMOS_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,            & ! Input Linearization control
              NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,               & ! Input numbers
              N_ALLLAYERS_DN, PARTLAYERS_LAYERIDX, STERM_LAYERMASK_DN, USER_STREAMS, & ! Input level control, user-streams
              T_DELT_USERM, T_UTDN_USERM, L_T_DELT_USERM, L_T_UTDN_USERM,            & ! Input transmittances
              DELTAU_POWER, XTAU_POWER, L_DELTAU_POWER, L_XTAU_POWER,                & ! Input Thermal setups
              U_TNEG1, U_TNEG2, L_U_TNEG1, L_U_TNEG2,                                & ! Input Thermal solutions               
              T_DIRECT_DN, T_UT_DIRECT_DN, L_T_DIRECT_DN, L_T_UT_DIRECT_DN,          & ! Input thermal direct solutions
              LAYER_TSUP_DN, LAYER_TSUP_UTDN, L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN )     ! Output user RTE thermal
          ENDIF

!  No linearization
!  @@@ Rob fix 1/31/11, This line was wrong (see below)
!          PI_XQP, PI_XQM, PI_XUP, PI_XQM_PRE, &

        ELSE

!  discrete ordinate particular integral for thermal sources.
!    - Rob fix 1/31/11, This line was wrong (see below)     PI_XQP, PI_XQM, PI_XUP, PI_XQM_PRE, &
!mick fix 9/19/2017 - added DO_USER_STREAMS to input
!    - 1/31/21. Version 2.8.3. Drop LOCAL_UM_START argument. User-stream do-loops start with 1 now.

          CALL THERMAL_CLSOLUTION ( &
            DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                 & ! input flags
            DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_THERMAL_TRANSONLY,    & ! Input flags
            NSTREAMS, NLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,         & ! Input basic numbers
            NMOMENTS, NSTREAMS_2, N_PARTLAYERS,                          & ! Input numbers
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Input level control
            QUAD_STREAMS, QUAD_HALFWTS, OMEGA_GREEK, SAB, DAB,           & ! Input optical and SAB/DAB
            T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,              & ! Input Discrete Ord. Trans.
            PI_XQP, PI_XUP, PI_XUM, PI_XQM_PRE,                          & ! Input PI matrices, Corrected 1/31/11
            THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,                & ! Input thermal setups
            T_WUPPER, T_WLOWER, UT_T_PARTIC,                             & ! Output thermal solutions
            U_TPOS1, U_TNEG1, U_TPOS2, U_TNEG2,                          & ! Output User thermal solutions
            STATUS_SUB, MESSAGE, TRACE_1 )                                 ! Exception handling

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            TRACE_2 = 'Called in VLIDORT_LCS_FOURIER, Fourier 0'
            STATUS = VLIDORT_SERIOUS
            RETURN
          ENDIF

!  User solutions and post processing.
!mick fix 3/30/2015 - modified if conditions for up & dn thermal source terms
!    - 1/31/21. Version 2.8.3. Drop LOCAL_UM_START argument. User-stream do-loops start with 1 now.

          IF ( DO_UPWELLING .AND. DO_USER_STREAMS) THEN
            CALL THERMAL_STERMS_UP ( &
              DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, DO_PARTLAYERS,    & ! Input flags
              NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,  & ! Input numbers
              N_ALLLAYERS_UP, PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP,  & ! Input level control
              USER_STREAMS, T_DELT_USERM, T_UTUP_USERM,                 & ! Input User streams and transmittances
              DELTAU_POWER, XTAU_POWER, U_TPOS1, U_TPOS2,               & ! Input Thermal setups/solutions
              T_DIRECT_UP, T_UT_DIRECT_UP,                              & ! Input thermal direct solutions
              LAYER_TSUP_UP, LAYER_TSUP_UTUP )                            ! Output user RTE thermal
          ENDIF

          IF ( DO_DNWELLING .AND. DO_USER_STREAMS ) THEN
            CALL THERMAL_STERMS_DN ( &
              DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, DO_PARTLAYERS,    & ! Input flags
              NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,  & ! Input numbers
              N_ALLLAYERS_DN, PARTLAYERS_LAYERIDX, STERM_LAYERMASK_DN,  & ! Input level control
              USER_STREAMS, T_DELT_USERM, T_UTDN_USERM,                 & ! Input User streams and transmittances
              DELTAU_POWER, XTAU_POWER, U_TNEG1, U_TNEG2,               & ! Input Thermal setups/solutions
              T_DIRECT_DN, T_UT_DIRECT_DN,                              & ! Input thermal direct solutions
              LAYER_TSUP_DN, LAYER_TSUP_UTDN )                            ! Output user RTE thermal
          ENDIF

        ENDIF

!  End include thermal emission

      ENDIF

!  Skip the thermal-only section if there are solar sources
!      IF ( DO_SOLAR_SOURCES ) GO TO 455

!  Version 2.8, Get rid of the GO TO 455 statement

      IF ( .NOT. DO_SOLAR_SOURCES ) THEN

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Only one solution, local direct_beam flag NOT set
!mick note 9/19/2017 - DO_INCLUDE_DIRECTBEAM in the thermal-only case here
!  covers the roles of both DO_LOCALBEAM(IBEAM) and DO_INCLUDE_DIRECTBEAM
!  in the solar case later

        IPARTIC = 1
!        DO_INCLUDE_DIRECTBEAM = .FALSE.       ! replaced 4/28/19
        DO_INCLUDE_DIRECTRF = .FALSE.
        DO_INCLUDE_DIRECTSL = .FALSE.

!  Version 2.8, Get rid of the GOTO 566 statement
!  Avoid the scattering solutions if not flagged
!      IF ( DO_THERMAL_TRANSONLY ) GO TO 566

        IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN

!  Thermal-only. Find the BVP solution and intensity field
!  -------------------------------------------------------

!  set the BVP PI solution at the lower/upper boundaries
!          O1 = 1
!          DO LAYER = 1, NLAYERS
!            DO I = 1, NSTREAMS_2
!              WUPPER(I,O1,LAYER) = T_WUPPER(I,LAYER)
!              WLOWER(I,O1,LAYER) = T_WLOWER(I,LAYER)
!            ENDDO
!          ENDDO
!mick fix 2/17/11 - include remaining stokes components as needed
!Rob Fix  April 2016. Original code was correct
!          IF (LC_NSTK .EQ. 1) THEN
!            DO LAYER = 1, NLAYERS
!              DO I = 1, NSTREAMS_2
!                WUPPER(I,1,LAYER) = T_WUPPER(I,LAYER) ; WLOWER(I,1,LAYER) = T_WLOWER(I,LAYER)
!              ENDDO
!            ENDDO
!          ELSE
!            DO LAYER = 1, NLAYERS
!              DO O1 = 1, LC_NSTK
!                DO I = 1, NSTREAMS_2
!                  WUPPER(I,O1,LAYER) = T_WUPPER(I,LAYER) ; WLOWER(I,O1,LAYER) = T_WLOWER(I,LAYER)
!                ENDDO
!              ENDDO
!            ENDDO
!          END IF

!mick fix 9/19/2017 - still bug with original code;
!                     new fix: initialize all elements of WUPPER and WLOWER
          WUPPER = ZERO
          WLOWER = ZERO
          O1 = 1
          WUPPER(1:NSTREAMS_2,O1,1:NLAYERS) = T_WUPPER(1:NSTREAMS_2,1:NLAYERS)
          WLOWER(1:NSTREAMS_2,O1,1:NLAYERS) = T_WLOWER(1:NSTREAMS_2,1:NLAYERS)

!  Solve the boundary value problem. No telescoping here
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
      
!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19
!      -- Rob fix  4/9/2019  - Major overhaul for adjusted BVP solution

!  1/31/21. Version 2.8.3.  No Change in this calling routine
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK etc....

          CALL BVP_SOLUTION_MASTER ( &
             DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_SURFEMISS,          & ! Input Surface Flags
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_INCLUDE_DIRECTRF,                 & ! Input Surface Flags
             DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX, FOURIER, IPARTIC,                 & ! Input illumination flags, indices
             DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,                    & ! Input Water-leaving flags
             LC_NSTK, NSTREAMS, NLAYERS, NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2,     & ! Numbers
             LC_NTT, LC_NSUBD, LC_NSUPD, K_REAL, K_COMPLEX, MUELLER_INDEX,             & ! Numbers, bookkeeping
             TF_MAXITER, TF_CRITERION, FLUX_FACTOR, FLUXVEC, SURFACE_FACTOR,           & ! Input factors/WL control
             QUAD_STRMWTS, ALBEDO, AXBID_F, SURFBB, EMISSIVITY,                        & ! Input surface terms
             SLTERM_ISOTROPIC, SLTERM_F_0, TOAFLUX, BOAFLUX,                           & ! Input Sleave/illumination
             SOLARBEAM_BOATRANS, LOCAL_CSZA, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, & ! Input Direct-flux
             RF_DIRECT_BEAM, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, WUPPER, WLOWER,       & ! Input RTE solutions
             BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                                         & ! Input BVP matrices/pivots             
             COL2, SCOL2, TRANS_ATMOS_FINAL, SL_QUADTERM, R2_BEAM, LCON, MCON,         & ! Modified input/output
             STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                                     ! Exception handling

!  Exception handling

          IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
            write(CF,'(I2)')FOURIER
            TRACE_3 = 'Error return from BVP_SOLUTION_MASTER (thermal-only). ' &
                       //'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
            STATUS = VLIDORT_SERIOUS ; RETURN
          ENDIF

!  Continuation point for avoiding thermal scattering
! 566  CONTINUE. removed, Versio 2.8

!  End thermal scattering clause

        ENDIF

!  Post-processing - upwelling thermal-only field
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    - Version 2.8.1, Control for BOA illumination added, 3/23/19
!    - Version 2.8.1, 4/9/19. DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL added
!    - Version 2.8.1, 4/9/19. RF_USER_DIRECT_BEAM, SL_USERTERM have been added

!  1/31/21. Version 2.8.3. SEVERAL INNOVATIONS and CHANGES
!    -- Input taylor Order TAYLOR_ORDER
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F, SURF_MSSTS_F
!    -- DO_CLASSICAL_SOLUTION flag added to control solution method (Greens vs. Classical)
!    -- Following additional arrays are Inputs and outputs for the Green's function
!           DELTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, ITRANS_USERM, & ! Input Greens function stuff
!           ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, SIGMA_P,                   & ! Input Greens function stuff
!           PMULT_UU, PMULT_UD,  UT_PMULT_UU, UT_PMULT_UD,                       & ! Output Greens function multipliers
!    -- Use post-processing masks N_PPSTREAMS,PPSTREAM_MASK
!    -- (RTS 2/16/21). Use local bookkeeping variable LC_NSTK, as well as NSTOKES

        IF ( DO_UPWELLING ) THEN
          CALL VLIDORT_UPUSER_INTENSITY ( &
            DO_USER_STREAMS,    DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,            & ! Input flags (RT mode)
            DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,           & ! Input flags (sources)
            DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,   DO_INCLUDE_SURFEMISS,           & ! Input flags (Surface)
            DO_DBCORRECTION,     DO_INCLUDE_DIRECTRF,     DO_INCLUDE_DIRECTSL,           & ! Input flags (Surface)
            DO_LAYER_SCATTERING, DO_MSMODE_VLIDORT,       DO_MSMODE_THERMAL,             & ! Input flags (RT mode)
            DO_TOA_CONTRIBS, DO_INCLUDE_BOAFLUX, DO_CLASSICAL_SOLUTION, DO_MSSTS,        & ! Input flags (RT mode)
            FOURIER, IPARTIC, NSTOKES, LC_NSTK, NSTREAMS, NLAYERS, N_USER_LEVELS,        & ! Input numbers (basic)
            N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_UP, MUELLER_INDEX, TAYLOR_ORDER,       & ! Input bookkeeping + levels
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                & ! Input partial-layer control
            FLUX_MULT, BOAFLUX, SURFACE_FACTOR, QUAD_WEIGHTS, QUAD_STRMWTS,              & ! Input Flux and quadrature
            T_DELT_DISORDS,  T_DELT_USERM, T_UTUP_USERM, CUMTRANS,                       & ! Input Transmittances
            ALBEDO, BRDF_F, USER_BRDF_F, SURFBB, EMISSIVITY, USER_EMISSIVITY,            & ! Input Surface BRDF/Emiss.
            K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,                       & ! Input Homog. RTE Soln.
            WLOWER, LCON, MCON, T_WLOWER, LAYER_TSUP_UP, LAYER_TSUP_UTUP,                & ! Input RTE PI and thermal
            DELTAU_VERT, PARTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,          & ! Input Beam for Greens 
            T_UTDN_MUBAR, ITRANS_USERM, ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,        & ! Input Greens function
            RF_USER_DIRECT_BEAM, SL_USERTERM, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2,& ! Input User solutions
            HMULT_1, HMULT_2, EMULT_UP, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, SIGMA_P,  & ! Input multipliers
            PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,                                & ! Output multipliers (Greens)
            STOKES_DOWNSURF, BOA_THTONLY_SOURCE, MS_CONTRIBS_F,                          & ! Output 1 (Auxiliary)
            STOKES_F, CUMSOURCE_UP, LAYER_MSSTS_F, SURF_MSSTS_F )                          ! Output 2 (Main)
        ENDIF

!  Post-processing - Downwelling thermal-only field
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    - Version 2.8.1, Control for TOA illumination added, 3/23/19

!  1/31/21. Version 2.8.3. SEVERAL INNOVATIONS and CHANGES
!    -- Input taylor Order TAYLOR_ORDER
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F
!    -- DO_CLASSICAL_SOLUTION flag added to control solution method (Greens vs. Classical)
!    -- Following additional arrays are Inputs and outputs for the Green's function
!           DELTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, ITRANS_USERM, & ! Input Greens function stuff
!           ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, SIGMA_M,                   & ! Input Greens function stuff
!           PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,                        & ! Output Greens function multipliers
!    -- Use N_PPSTREAMS and PPSTREAM_MASK, post-processing mask
!    -- (RTS 2/16/21). Use local bookkeeping variable LC_NSTK, as well as NSTOKES

        IF ( DO_DNWELLING ) THEN
          CALL VLIDORT_DNUSER_INTENSITY ( &
            DO_USER_STREAMS,     DO_OBSERVATION_GEOMETRY, DO_INCLUDE_TOAFLUX,        & ! Input flags (RT mode)
            DO_SOLAR_SOURCES,    DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,      & ! Input flags (sources)
            DO_LAYER_SCATTERING, DO_MSMODE_VLIDORT, DO_CLASSICAL_SOLUTION, DO_MSSTS, & ! Input flags (RT mode)
            FOURIER, IPARTIC, NSTOKES, LC_NSTK, NLAYERS, N_USER_LEVELS, N_PPSTREAMS, & ! Input numbers (basic)
            PPSTREAM_MASK, LEVELMASK_DN, TAYLOR_ORDER,                               & ! Input bookkeeping + levels
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,            & ! Input partial-layer control
            FLUX_MULT, TOAFLUX, T_DELT_USERM, T_UTDN_USERM,                          & ! Input Transmittances, Flux
            K_REAL, K_COMPLEX, LCON, MCON, LAYER_TSUP_DN, LAYER_TSUP_UTDN,           & ! Input RTE Sol + thermal
            DELTAU_VERT, PARTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,      & ! Input Beam for Greens
            T_UTDN_MUBAR, ITRANS_USERM, ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,    & ! Input Greens function 2   
            UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2, SIGMA_M,                     & ! Input User solutions
            HMULT_1, HMULT_2, EMULT_DN, UT_HMULT_DD, UT_HMULT_DU, UT_EMULT_DN,       & ! Input multipliers
            PMULT_DD, PMULT_DU, UT_PMULT_DU, UT_PMULT_DD,                            & ! Output mutlipliers (Greens)
            STOKES_F, CUMSOURCE_DN, LAYER_MSSTS_F  )                                   ! Main output
        ENDIF

!  mean value (integrated) output for thermal-only field
!    - Can also use this to get debug Quadrature output
!    - @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16. Removed QUAD_OUTPUT flag.
!mick mod 9/19/2017 - added LEVELS_SOLARTRANS & PARTIALS_SOLARTRANS to input arguments

!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19
        
!  1/31/21. Version 2.8.3.
!   -- Add flag DO_CLASSICAL_SOLUTION, controls use of Greens function
!   -- Additional inputs BVEC and INITIAL_TRANS, T_UTDN_MUBAR, TAYLOR ORDER, LC_NSTKNSTRM
!   -- Additional inputs ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P foro the Green's function
!   -- Additional outputs are the partial-layer Green's function multipliers UT_GMULT_UP, UT_GMULT_DN
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

        IF ( DO_INCLUDE_MVOUTPUT ) THEN
          CALL VLIDORT_INTEGRATED_OUTPUT ( &
             DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTRF,                         & ! Input flags
             DO_CLASSICAL_SOLUTION, DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX,    & ! Input flags
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,    & ! Input flags
             IPARTIC, LC_NSTK, NSTREAMS, NLAYERS, N_USER_LEVELS, N_DIRECTIONS, & ! Input numbers
             LC_NSTKNSTRM, LEVELMASK_UP, LEVELMASK_DN, WHICH_DIRECTIONS,    & ! Input Bookkeeping
             PARTLAYERS_LAYERIDX, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,     & ! Input partial layer control
             TAYLOR_ORDER, FLUX_MULT, TOAFLUX, BOAFLUX, FLUX_FACTOR,           & ! Flux inputs
             FLUXVEC, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, PARTAU_VERT,   & ! Quadrature inputs
             BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,           & ! Solar beam Param.
             LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, LOCAL_CSZA,               & ! Solar beam transmittances
             T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                   & ! Discrete ordinate transmittances
             ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P,                         & ! Input Greens function
             K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, LCON, MCON,              & ! Input RTE
             T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, BVEC, WUPPER, WLOWER,   & ! Input RTE
             T_WLOWER, T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,              & ! Input Thermal
             MEANST_DIFFUSE, FLUX_DIFFUSE, DNMEANST_DIRECT, DNFLUX_DIRECT,     & ! MAIN output
             UT_GMULT_UP, UT_GMULT_DN )                                          ! Auxiliary Output
        ENDIF

!  LTE linearization (new section, 14 September 2009)
!  This routine is now defunct in Verison 2.7. Replaced by LBBF stuff. 3/28/14
!      IF ( DO_THERMAL_TRANSONLY .AND. DO_LTE_LINEARIZATION ) THEN
!        CALL THERMAL_LTE_LINEARIZATION ( &
!          DO_INCLUDE_SURFACE, SURFACE_FACTOR,  FLUX_MULT, &
!          DO_UPWELLING, DO_DNWELLING, NSTREAMS, NLAYERS, N_USER_LEVELS, &
!          DELTAU_VERT_INPUT, ALBEDO, THERMAL_BB_INPUT, &
!          QUAD_STREAMS, QUAD_STRMWTS, N_USER_STREAMS, LOCAL_UM_START, USER_STREAMS, &
!          LEVELMASK_UP, LEVELMASK_DN, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
!          T_DELT_DISORDS, T_DELT_USERM, CUMSOURCE_UP, CUMSOURCE_DN, THERMCOEFFS, &
!          LTE_DELTAU_VERT_INPUT, LTE_THERMAL_BB_INPUT, LTE_ATMOSWF )
!      ENDIF

!  Thermal only. Avoid weighting functions all together
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        IF ( .NOT. DO_SIMULATION_ONLY ) THEN

!  Thermal Only: Atmospheric Bulk weighting functions
!  --------------------------------------------------

          IF ( DO_COLUMN_LINEARIZATION ) THEN

!   variation index = 0

            VARIATION_INDEX = 0

!  Avoidance of BVP problem for transmittance only
!   Version 2.8. remove GOTO 789 statement
!        IF ( DO_THERMAL_TRANSONLY ) GO TO 789

            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN

!  4/9/19. Additional code to set the LC derivative of TRANS_ATMOS_FINAL
!     Here in the thermal-regime, there is no water-leaving, so initialize to zero          

              LC_TRANS_ATMOS_FINAL(IBEAM,1:N_TOTALCOLUMN_WFS) = ZERO
          
!  Solve the Regular BVP linearization

              IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN

!  Get the linearized BVP solution for this beam component
!  4/9/19. Water-leaving control. Also need LC_TRANS_ATMOS_FINAL

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> Several New Inputs arguments to LC_BVP_SOLUTION_MASTER, LC_BVPTEL_SOLUTION_MASTER : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LC_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE
!           ** BRDF Fourier inputs are defined locally for each Fourier component
!mick fix 1/5/2021 - adjusted ordering of some input vars based on updated LC_BVP_SOLUTION_MASTER subroutine statement
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

                CALL LC_BVP_SOLUTION_MASTER ( &
                   DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTRF,         & ! Flags
                   DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,              & ! Flags
                   DO_CLASSICAL_SOLUTION, DO_LAYER_SCATTERING, FOURIER, IPARTIC,           & ! Flags/Indices
                   LC_NSTK, NSTREAMS, NLAYERS, N_TOTALCOLUMN_WFS, TAYLOR_ORDER,            & ! Numbers (Basic)
                   NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2, LC_NTT, LC_NSUBD, LC_NSUPD,   & ! Numbers (derived)
                   MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR, QUAD_STRMWTS,         & ! Bookkeeping
                   DELTAU_SLANT, DELTAU_VERT, L_DELTAU_VERT, RF_DIRECT_BEAM, SL_QUADTERM,      & ! Optical and direct beam
                   BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                   & ! Input, Beam parameterization
                   LC_TRANS_ATMOS_FINAL, LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Input, Linearized Beam param
                   CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
                   SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, BVEC,             & ! Input, Homogeneous/Classical
                   BANDMAT2, IPIVOT, SMAT2, SIPIVOT, ALBEDO, BRDF_F,                 & ! Input, BVP matrices, surface
                   L_T_DELT_EIGEN, L_KEIGEN, L_SOLA_XPOS, L_SOLB_XNEG,               & ! Input, Linearized Homog solutions
                   L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER, LC_BVEC,      & ! Input, Linearized Greens/Thermal/BVEC
                   NCON, PCON, L_WLOWER, L_WUPPER,                                   & ! output - Linearized Constants + PI
                   STATUS_SUB, MESSAGE, TRACE_1 )                                      ! Exception handling

!  Exception handling

                IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  TRACE_2 = 'Error return from LC_BVP_SOLUTION_MASTER, thermal-only solutions '
                  TRACE_3 =  'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
                  STATUS = VLIDORT_SERIOUS ; RETURN
                ENDIF

              ENDIF

!   Version 2.8. remove GOTO 789 statement
!  Continuation point for avoiding scattering solution
! 789    CONTINUE

!   End thermal scattering clause

            ENDIF

!  Post-processing for Column weighting functions (upwelling)
!  ----------------------------------------------------------

!  Streamlined for  Version 2.8. 7/22/16
!   4/9/19.  For the Column-Linearized BOA source term, need additional inputs
!             -- Used 2 flags to control direct radiances (reflected-beam, surface-leaving)
!             -- Note use of LC_TRANS_ATMOS_FINAL to go with adjusted water-leaving

!  1/31/21. Version 2.8.3.
!    ==> New input arguments to the post processing (Green's function) : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LC_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE, SIGMA_P
!           ** Green's function arrays ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_UU, PMULT_UD
!    ==> Other Input/Output arguments and other changes
!           ** for the multiple scatter source term linearization, use control flag DO_MSSTS
!           ** LInearizations of MSST functions (LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F), now output
!           ** BRDF Fourier inputs are defined locally for each Fourier component
!   -- (RTS 2/16/21). Use local bookkeeping variable LC_NSTK in addition to NSTOKES

            IF ( DO_UPWELLING ) THEN
              CALL UPUSER_COLUMNWF ( &
                DO_USER_STREAMS,    DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,                 & ! Input flags (RT mode)
                DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,                & ! Input flags (sources)
                DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT,    DO_INCLUDE_SURFACE,                  & ! Input flags (sources)
                DO_LAMBERTIAN_SURFACE, DO_DBCORRECTION, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,     & ! Input flags (Surface)
                DO_LAYER_SCATTERING, DO_MSSTS, FOURIER, IPARTIC, NSTOKES, LC_NSTK, NSTREAMS, NLAYERS, & ! Input flags/numbers (basic)
                N_PPSTREAMS, PPSTREAM_MASK, N_USER_LEVELS, N_TOTALCOLUMN_WFS,                     & ! Input numbers (basic)
                LEVELMASK_UP, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,       & ! Input Level output control
                TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,              & ! Input Taylor/Optical
                FLUX_MULT, QUAD_WEIGHTS, QUAD_STRMWTS, USER_SECANTS, MUELLER_INDEX,               & ! Input Flux/quad/Mueller
                BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_DISORDS,           & ! Input Beam param
                T_DELT_USERM, T_UTUP_USERM, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,          & ! Input UTrans/Surface.
                RF_USER_DIRECT_BEAM, SL_USERTERM, LC_TRANS_ATMOS_FINAL, CUMSOURCE_UP,             & ! Input surfradiances/Cumsource
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON,                & ! Input Homogeneous Solutions
                ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_UU, PMULT_UD, SIGMA_P,            & ! Input Green's function vars
                T_WLOWER, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, HMULT_1, HMULT_2,           & ! Input Thermal/Homog
                EMULT_UP, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, UT_PMULT_UU, UT_PMULT_UD,        & ! Input multipliers  
                LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, L_T_DELT_DISORDS,           & ! Lin BeamParm/DODTrans
                L_T_DELT_USERM, L_T_UTUP_USERM, L_T_WLOWER, L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP,   & ! Lin Trans/Thermal
                L_KEIGEN, L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG, L_UHOM_UPDN, L_UHOM_UPUP,     & ! Lin Homog solutions
                L_HMULT_1, L_HMULT_2, LC_EMULT_UP, L_UT_HMULT_UU, L_UT_HMULT_UD, LC_UT_EMULT_UP,  & ! Lin multipliers
                L_UPAR_UP_1, LC_UPAR_UP_2, L_ATERM_SAVE, L_BTERM_SAVE, L_WLOWER, NCON, PCON,      & ! Lin PI Solutions
                L_BOA_THTONLY_SOURCE, COLUMNWF_F, LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F )               ! Output
            ENDIF

!  Post-processing for Column weighting functions (Downwelling)
!  ------------------------------------------------------------

!  Streamlined for  Version 2.8. 7/22/16

!  1/31/21. Version 2.8.3.
!    ==> New input arguments to the post processing (Green's function) : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LC_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE, SIGMA_M
!           ** Green's function arrays ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_DU, PMULT_DD
!    ==> Other Input/Output arguments and other changes
!           ** for the multiple scatter source term linearization, use control flag DO_MSSTS
!           ** LInearizations of MSST functions (LC_LAYER_MSSTS_F), now output
!           ** BRDF Fourier inputs are defined locally for each Fourier component
!   -- (RTS 2/16/21). Use local bookkeeping variable LC_NSTK in addition to NSTOKES
!mick fix 1/5/2021 - added NLAYERS to inputs

            IF ( DO_DNWELLING ) THEN
              CALL DNUSER_COLUMNWF ( &
                DO_USER_STREAMS,   DO_OBSERVATION_GEOMETRY, DO_CLASSICAL_SOLUTION,                & ! Input flags (RT mode)
                DO_SOLAR_SOURCES,  DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,                 & ! Input flags (sources)
                DO_MSMODE_VLIDORT, DO_LAYER_SCATTERING,     DO_MSSTS,  FOURIER, IPARTIC,          & ! Input flags (sources)
                NSTOKES, LC_NSTK, NLAYERS, N_USER_LEVELS, N_TOTALCOLUMN_WFS,                      & ! Input flags/numbers (basic)
                N_PPSTREAMS, PPSTREAM_MASK,                                                       & ! Input flags/numbers (basic)
                LEVELMASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,       & ! Input Level output control
                TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT, FLUX_MULT, USER_SECANTS,   & ! Input Taylor/Opt/Secants
                BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,           & ! Input Beam param
                T_DELT_USERM, T_UTDN_USERM, CUMSOURCE_DN, K_REAL, K_COMPLEX, LCON, MCON,          & ! Input User/Cumsource/BVP
                ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_DU, PMULT_DD, SIGMA_M,            & ! Input Green's function vars
                UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2, HMULT_1, HMULT_2,                     & ! Homog
                EMULT_DN, UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, UT_PMULT_DU, UT_PMULT_DD,        & ! Input multipliers  
                L_T_DELT_USERM, L_T_UTDN_USERM, LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Lin BeamParm/Trans
                L_KEIGEN, L_UHOM_DNDN, L_UHOM_DNUP, L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN,               & ! Lin Homog/Thermal
                L_HMULT_1, L_HMULT_2, LC_EMULT_DN, L_UT_HMULT_DU, L_UT_HMULT_DD, LC_UT_EMULT_DN,      & ! Lin multipliers
                L_UPAR_DN_1, LC_UPAR_DN_2, L_ATERM_SAVE, L_BTERM_SAVE, NCON, PCON,                    & ! Lin PI Solutions
                COLUMNWF_F, LC_LAYER_MSSTS_F )                                                          ! Output
            ENDIF

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument
!    - Version 2.8, Cleaned-up I/O presentation. 7/22/16. Removed QUAD_OUTPUT flag.
!mick mod 9/19/2017 - added LC_LEVELS_SOLARTRANS & LC_PARTIALS_SOLARTRANS to input arguments

!  1/31/21. Version 2.8.3. Rewritten to incorporate new Linearized Green's function solution
!     -- New arguments include DO_CLASSICAL, TAYLOR_ORDER, PARTAU_VERT, L_DELTAU_VERT
!     -- New Greens function variables ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, UT_GMULT_UP, UT_GMULT_DN
!     -- Logic changed to perform additional Partial-layer Green function multipliers initially
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

            IF ( DO_INCLUDE_MVOUTPUT ) THEN
              CALL MIFLUX_COLUMNWF ( &
                DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTRF, DO_CLASSICAL_SOLUTION,           & ! Input flags
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,             & ! Input flags
                IPARTIC, N_TOTALCOLUMN_WFS, LC_NSTK, NSTREAMS, LC_NSTKNSTRM, NLAYERS,      & ! Input numbers
                N_USER_LEVELS, N_DIRECTIONS, WHICH_DIRECTIONS, LEVELMASK_UP, LEVELMASK_DN, & ! Level/Dir output control
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,              & ! Partial layer output control
                FLUX_MULT, FLUX_FACTOR, FLUXVEC, LOCAL_CSZA, TAYLOR_ORDER,                 & ! Input Flux/Angles/Taylor
                QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, PARTAU_VERT, L_DELTAU_VERT,      & ! Input quadrature/Optical
                BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! Beam parameterization
                T_DELT_DISORDS, T_UTDN_DISORDS, T_UTUP_DISORDS, K_REAL, K_COMPLEX,         & ! DODTrans/Bookkeep
                SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,            & ! Homog solution
                ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, UT_GMULT_UP, UT_GMULT_DN,        & ! Green's function solution
                WUPPER, LCON, MCON, T_WLOWER, T_WUPPER, BOA_THTONLY_SOURCE,                & ! Solar/BVP/Thermal
                LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,                      & ! Lin Beam parameterization
                LC_T_UTDN_MUBAR, LC_LEVELS_SOLARTRANS, LC_PARTIALS_SOLARTRANS,             & ! Lin Beam transmittances
                L_T_DELT_DISORDS, L_T_UTDN_DISORDS, L_T_UTUP_DISORDS, L_KEIGEN,            & ! Lin DODTrans/Keigen
                L_SOLA_XPOS, L_SOLB_XNEG, L_T_DELT_EIGEN, L_T_UTDN_EIGEN, L_T_UTUP_EIGEN,  & ! Lin Homog solution
                L_WLOWER, L_WUPPER, NCON, PCON, L_ATERM_SAVE, L_BTERM_SAVE,                & ! Lin PI/Greens Solutons
                L_UT_T_PARTIC, L_T_WLOWER, L_T_WUPPER, L_BOA_THTONLY_SOURCE,               & ! Lin Thermal solutions
                MEANST_DIFFUSE_COLWF, DNMEANST_DIRECT_COLWF,                               & ! Lin Output Actinic Fluxes
                FLUX_DIFFUSE_COLWF, DNFLUX_DIRECT_COLWF )                                    ! Lin Output Regular Fluxes
            ENDIF

!  End atmospheric column weighting functions

          ENDIF

!  Thermal-only: Surface Reflectance weighting functions
!  -----------------------------------------------------

!mick fix 9/6/2012 - added (N_SURFACE_WFS > 0) condition to IF

          IF ( DO_SURFACE_LINEARIZATION .AND. DO_INCLUDE_SURFACE .AND. (N_SURFACE_WFS > 0) ) THEN

!  4/9/19. Additional code to set the LS derivative of TRANS_ATMOS_FINAL
!     Here in the thermal-regime, there is no water-leaving, so initialize to zero          

            LS_TRANS_ATMOS_FINAL(IBEAM,1:N_SURFACE_WFS) = ZERO
          
!  Surface WFs; Solve boundary value problem
!mick fix 9/19/2017 - replaced DO_REFLECTED_DIRECTBEAM(IBEAM) with DO_LOCALBEAM(IBEAM)
!                     which in turn is replaced by DO_INCLUDE_DIRECTRF here in the thermal

!   4/9/19. Regular Solution. Introduce variability for adjusted water-leaving transmittance (handle only)

!  1/31/21. Version 2.8.3. BRDF arrays are defined locally for each Fourier
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

            CALL SURFACEWF_BVP_SOLUTION ( FOURIER, &
                DO_INCLUDE_DIRECTRF, DO_INCLUDE_SURFEMISS, DO_LAMBERTIAN_SURFACE,     & ! Input flags
                DO_WATER_LEAVING, IPARTIC, LC_NSTK, NSTREAMS, NLAYERS, N_SURFACE_WFS, & ! Input
                LC_NSTKNSTRM, LC_NSTKNSTRM_2, LC_NTT, LC_NSUPD, LC_NSUBD,             & ! Input
                MUELLER_INDEX, K_REAL, K_COMPLEX, QUAD_STRMWTS,                       & ! Input
                SURFACE_FACTOR, ATMOS_ATTN, LS_TRANS_ATMOS_FINAL, SL_QUADTERM,        & ! input surface
                SURFBB, LS_EMISSIVITY, LS_BRDF_F, LS_BRDF_F_0,                        & ! Input
                T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, WLOWER,                           & ! Input
                BANDMAT2, IPIVOT, SMAT2, SIPIVOT, LCON, MCON,                         & ! Input
                NCON_SWF, PCON_SWF, STATUS_SUB, MESSAGE, TRACE_1 )                      ! Output

!  Exception handling

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_2 = 'Error return from SURFACEWF_BVP_SOLUTION, '//'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
              STATUS = VLIDORT_SERIOUS ; RETURN
            ENDIF

!  Surface WFs: Postprocessing Master. [DO_INCLUDE_DIRECTRF = .false. (automatic here)]
!  4/9/19. Separate terms from surface (reflected-DB and Sleave). Revised I/O list.

!  1/31/21. Version 2.8.3. Several minor changes to this call.
!    -- Drop LOCAL_UM_START; all UM loops start with 1
!    -- Use the post processing mask system, don't need separate OBSGEOM calculation
!    -- Introduce DO_MSSTS flag, calculate linearized MSSTS output LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F
!    -- BRDF arrays USER_BRDF_F, LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0 all defined locally each Fourier
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others
!mick fix 1/5/2021 - added NSTOKES back to inputs (in addition to LC_NSTK)

            CALL SURFACEWF_POSTPROCESS_MASTER &
              ( DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_MSSTS,                       & ! Input flags (general)
                DO_DBCORRECTION, DO_INCLUDE_MVOUTPUT, DO_OBSERVATION_GEOMETRY,               & ! Input flags (general)
                DO_INCLUDE_SURFEMISS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,               & ! Input flags (thermal)
                DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,             & ! Input flags (surface)
                FOURIER, IPARTIC, NSTOKES, LC_NSTK, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS, & ! Input Control numbers
                N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_UP, LEVELMASK_DN,                      & ! Input Level output control
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                & ! Input partlayer control
                N_DIRECTIONS, WHICH_DIRECTIONS, MUELLER_INDEX,                               & ! Input bookkeeping
                FLUX_MULT, FLUXVEC, SURFACE_FACTOR, QUAD_WEIGHTS, QUAD_STRMWTS,              & ! Input flux/quads
                SL_USERTERM, ALBEDO, USER_BRDF_F, SURFBB, LS_USER_EMISSIVITY, LS_EMISSIVITY, & ! Input Surface stuff
                LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_TRANS_ATMOS_FINAL,           & ! Input Surface stuff
                T_DELT_DISORDS, T_UTUP_DISORDS, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,    & ! Input trans (User/Disords)
                K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,                       & ! Input Homog. solutions
                T_UTUP_EIGEN, T_UTDN_EIGEN, UHOM_UPDN, UHOM_UPUP, UHOM_DNDN, UHOM_DNUP,      & ! Input Homog. Solutions
                HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,  UT_HMULT_DU, UT_HMULT_DD,       & ! Input Homog.Multipliers
                ATMOS_ATTN, STOKES_DOWNSURF, NCON_SWF, PCON_SWF,                             & ! BOA downwelling, Lin BVP
                SURFACEWF_F, LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F, MEANST_DIFFUSE_SURFWF, FLUX_DIFFUSE_SURFWF) ! Output

!  end of surface weighting functions

          ENDIF

!  New, 28 March 2014. Linearization for BLACKBODY (Verson 2.7). Enabled for Version 2.8
!    -- 1/31/21. Version 2.8.3. BRDF arrays are defined locally for each Fourier
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

          IF ( ( DO_ATMOS_LBBF .OR. DO_SURFACE_LBBF ) .AND. Fourier.EQ.0 ) THEN
            IF ( N_PARTLAYERS .GT. 0 ) THEN
              CALL VLIDORT_LBBF_JACOBIANS_WPARTIALS &
               ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,      & ! Input
                 DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,              & ! Input
                 DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,   & ! input
                 DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,                 & ! input
                 LC_NSTK, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, & ! input
                 NMOMENTS, NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2,        & ! Input
                 LC_NTT, LC_NSUPD, LC_NSUBD,  MUELLER_INDEX,                & ! input
                 N_PARTLAYERS, PARTLAYERS_LAYERIDX,                         & ! Input
                 PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
                 LEVELMASK_UP, LEVELMASK_DN,                                & ! Input
                 USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input
                 QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                  & ! input
                 SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,               & ! input
                 EMISSIVITY, USER_EMISSIVITY, FLUX_MULT,                    & ! input
                 DELTAU_VERT, PARTAU_VERT, OMEGA_GREEK,                     & ! Input
                 T_DELT_DISORDS, T_UTDN_DISORDS, T_UTUP_DISORDS,            & ! input
                 PI_XQP, PI_XQP_PRE, PI_XUP, PI_XUM, SAB, DAB,              & ! input
                 K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,                   & ! input
                 T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,                  & ! input
                 T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                  & ! Input
                 BANDMAT2, IPIVOT, SMAT2, SIPIVOT, HMULT_1, HMULT_2,        & ! Input
                 UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,        & ! Input
                 UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                & ! input
                 ABBWFS_JACOBIANS, ABBWFS_FLUXES,                           & ! Output
                 SBBWFS_JACOBIANS, SBBWFS_FLUXES,                           & ! Output
                 STATUS_SUB, MESSAGE, TRACE_1 )                               ! Output
            ELSE
              CALL VLIDORT_LBBF_JACOBIANS_WHOLE &
               ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,      & ! Input
                 DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,              & ! Input
                 DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,   & ! input
                 DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,                 & ! input
                 LC_NSTK, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, & ! input
                 NMOMENTS, NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2,        & ! Input
                 LC_NTT, LC_NSUPD, LC_NSUBD,  MUELLER_INDEX,                & ! input
                 LEVELMASK_UP, LEVELMASK_DN,                                & ! Input
                 USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input
                 QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                  & ! input
                 SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,               & ! input
                 EMISSIVITY, USER_EMISSIVITY,                               & ! input
                 FLUX_MULT, DELTAU_VERT, OMEGA_GREEK,                       & ! Input
                 T_DELT_DISORDS, T_DELT_USERM,                              & ! input
                 PI_XQP, PI_XQP_PRE, PI_XUP, PI_XUM, SAB, DAB,              & ! input
                 K_REAL, K_COMPLEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG,     & ! input
                 BANDMAT2, IPIVOT, SMAT2, SIPIVOT, HMULT_1, HMULT_2,        & ! Input
                 UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                & ! input
                 ABBWFS_JACOBIANS, ABBWFS_FLUXES,                           & ! Output
                 SBBWFS_JACOBIANS, SBBWFS_FLUXES,                           & ! Output
                 STATUS_SUB, MESSAGE, TRACE_1 )                               ! Output
            ENDIF

!  Error handling

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
               TRACE_2 = 'Error return from vlidort_lbbf_jacobians, thermal only'// &
                         'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
              STATUS = VLIDORT_SERIOUS ; RETURN
            ENDIF

!  End of LBBF Jacobians

          ENDIF

!  End not simulation-only clause

        ENDIF

!  Finish Thermal only.

        RETURN

!  Continuation point. Removed for Version 2.8
! 455  CONTINUE

!  End of Thermal only calculation

      ENDIF

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!  Start loop over desired beam solutions

      DO IBEAM = 1, NBEAMS

        IF ( DO_MULTIBEAM(IBEAM,FOURIER) ) THEN

!  Step 1. Solar beam Particular solutions + linearizations
!  ========================================================

          DO LAYER = 1, NLAYERS

!  Parameter control for the linearization

            IF ( DO_ATMOS_LINEARIZATION ) THEN
              IF ( DO_COLUMN_LINEARIZATION ) THEN
                DO_RTSOL_VARY = .TRUE.
                NPARAMS_VARY  = N_TOTALCOLUMN_WFS
              ENDIF
            ENDIF

!  1/31/21. Version 2.8.3. Either Classical or Greens
!    -- Green's function option is new for this version (GBEAM solutions)

!  Classical
!  ---------

!    -- Version 2.8, Cleaned-up I/O presentations. 7/8/16.
!    -- Classical Subroutines renamed (QBEAM --> CBEAM). Argument list cleaned up.
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

            IF ( DO_CLASSICAL_SOLUTION ) THEN

!  regular

              CALL VLIDORT_CBEAM_SOLUTION &                                 ! CLASSICAL SOLUTION
                ( LAYER, FOURIER, IBEAM,                                  & ! Input indices
                  LC_NSTK, NSTREAMS, NMOMENTS, NSTREAMS_2, LC_NSTKNSTRM,  & ! Input numbers
                  FLUX_FACTOR, DFLUX, QUAD_STREAMS,                       & ! Input Flux and quadrature
                  DO_LAYER_SCATTERING, OMEGA_GREEK, BEAM_CUTOFF,          & ! Input bookkeeping + optical
                  T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT,            & ! Input beam attenuation
                  PI_XQP, PI_XQM, PI_X0P, PI_XQM_POST,                    & ! Input PI matrices
                  SAB, DAB, EIGENMAT_SAVE,                                & ! Input matrices from RTE
                  QSUMVEC_SAVE, QDIFVEC_SAVE, QVEC_SAVE, QDIF_SAVE,       & ! Output beam auxiliary vectors
                  QMAT_SAVE, QPIVOT, BVEC, WUPPER, WLOWER,                & ! Output matrix and Beam solutions
                  STATUS_SUB, MESSAGE, TRACE_1 )                            ! Exception handling

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER
                TRACE_2 = 'Error return from VLIDORT_CBEAM_SOLUTION, '//'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS ; RETURN
              ENDIF

!  Linearized

              IF ( DO_ATMOS_LINEARIZATION ) THEN

                CALL VLIDORT_LC_CBEAM_SOLUTION ( &
                  DO_PLANE_PARALLEL, DO_RTSOL_VARY, LAYER, FOURIER, IBEAM, & ! Input flags/numbers
                  LC_NSTK, NSTREAMS, NMOMENTS, NPARAMS_VARY,               & ! Input numbers
                  NSTREAMS_2, LC_NSTKNSTRM, DO_LAYER_SCATTERING,           & ! Input bookkeeping
                  QUAD_STREAMS, FLUX_FACTOR, BEAM_CUTOFF, AVERAGE_SECANT,  & ! Input quads/Beam
                  DFLUX, PI_XQP, PI_X0P, PI_XQM_POST, SAB, DAB,            & ! Input Eigenproblem
                  QSUMVEC_SAVE, QDIFVEC_SAVE, QVEC_SAVE, QDIF_SAVE,        & ! Input solution
                  QMAT_SAVE, QPIVOT, L_SAB, L_DAB, L_EIGENMAT,             & ! Input solution, lin Eigen.
                  L_OMEGA_GREEK, LC_AVERAGE_SECANT,                        & ! Input optical + beam
                  LC_BVEC, STATUS_SUB, MESSAGE, TRACE_1 )                    ! Output + status

                IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  TRACE_2 = 'Error from VLIDORT_LC_QBEAM_SOLUTION, '//'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
                  STATUS = VLIDORT_SERIOUS ; RETURN
                ENDIF

              ENDIF

!  Greens function
!  ---------------

!  New subroutines.
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

            ELSE

              CALL VLIDORT_GBEAM_SOLUTION &                                        ! GREENS FUNCTION
                ( FOURIER, IBEAM, LAYER, TAYLOR_ORDER,                                & ! Input
                  LC_NSTK, NSTREAMS, NSTREAMS_2, LC_NSTKNSTRM, NMOMENTS,              & ! input
                  FLUX_FACTOR, DFLUX, QUAD_WEIGHTS, DO_LAYER_SCATTERING, BEAM_CUTOFF, & ! Input
                  INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                        & ! Beam parameterization
                  DELTAU_VERT, OMEGA_GREEK, PI_XQP, PI_XQM, PI_X0P, PI_XQM_POST,      & ! GSF/Scat Optical
                  KEIGEN, KEIGEN_CSQ, K_REAL, K_COMPLEX,                              & ! Eigenvalues                      
                  SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, NORM_SAVED,                     & ! Eigensolutions
                  DPI, DMI, ATERM_SAVE, BTERM_SAVE, CFUNC, DFUNC,                     & ! Output
                  GAMMA_M, GAMMA_P, GFUNC_UP, GFUNC_DN, WUPPER, WLOWER )                ! Output

              IF ( DO_ATMOS_LINEARIZATION ) THEN
                CALL VLIDORT_L_GBEAM_SOLUTION &
                  ( FOURIER, IBEAM, LAYER, LC_NSTK, NSTREAMS,              & ! input Numbers
                    LC_NSTKNSTRM, NMOMENTS, DO_RTSOL_VARY, NPARAMS_VARY,   & ! input Numbers
                    DO_LAYER_SCATTERING, FLUX_FACTOR, DFLUX, QUAD_WEIGHTS, & ! Input Bookkeeping
                    BEAM_CUTOFF, PI_XQP, PI_XQM, PI_X0P, PI_XQM_POST,      & ! Input GSF
                    L_OMEGA_GREEK, K_REAL, SOLA_XPOS, SOLB_XNEG,           & ! Input Optical/Eigensolutions
                    NORM_SAVED, ATERM_SAVE, BTERM_SAVE, DMI, DPI,          & ! Input Green's function stuff
                    L_NORM_SAVED, L_SOLA_XPOS, L_SOLB_XNEG,                & ! Input linearized solutions
                    L_ATERM_SAVE, L_BTERM_SAVE )                             ! Output
              ENDIF

            ENDIF

!  user solutions
!  --------------

!  1/31/21. Version 2.8.3.
!    -- Either classical (CUSER, renamed) or Green's function (GUSER) solutions
!    -- Green's function option is new for this version.  
!    -- Drop LOCAL_UM_START Argument from both routines.
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

            IF ( STERM_LAYERMASK_UP(LAYER) .OR. STERM_LAYERMASK_DN(LAYER) ) THEN
              IF ( DO_USER_STREAMS ) THEN

!  Classical

                IF ( DO_CLASSICAL_SOLUTION ) THEN

                  CALL VLIDORT_CUSER_SOLUTION ( &                                         ! CLASSICAL SOLUTION
                    DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,                & ! Input flags
                    LAYER, FOURIER, IBEAM, LC_NSTK, NSTREAMS, N_USER_STREAMS,           & ! Input numbers
                    NMOMENTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                   & ! Input bookkeeping
                    DO_LAYER_SCATTERING, FLUX_FACTOR, DFLUX, BEAM_CUTOFF, QUAD_HALFWTS, & ! Input flux/quadrature/Bookkeeping
                    OMEGA_GREEK, PI_XQP, PI_XUP, PI_XUM, PI_X0P, PI_XQM_PRE, BVEC,      & ! Input optical/PI/Beam-PI
                    HELPSTOKES_BEAM, UPAR_DN_1, UPAR_DN_2, UPAR_UP_1, UPAR_UP_2 )         ! Output user solutions ( &

                  IF ( DO_ATMOS_LINEARIZATION ) THEN
                    CALL VLIDORT_LC_CUSER_SOLUTION ( &
                      DO_OBSERVATION_GEOMETRY, DO_UPWELLING, DO_DNWELLING,         & ! Input flags
                      DO_RTSOL_VARY, LAYER, FOURIER, IBEAM,                        & ! Input flags/indices
                      LC_NSTK, NSTREAMS, NMOMENTS, N_USER_STREAMS, NPARAMS_VARY,   & ! Input numbers
                      STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DO_LAYER_SCATTERING, & ! Bookkeeping
                      DFLUX, FLUX_FACTOR, QUAD_HALFWTS,                            & ! Bookkeeping and quads
                      PI_XQP, PI_XUP, PI_XUM, PI_X0P, PI_XQM_PRE, BEAM_CUTOFF,     & ! Pi Matrices and cutoff
                      OMEGA_GREEK, HELPSTOKES_BEAM, LC_BVEC, L_OMEGA_GREEK,        & ! Solutions and optical
                      L_UPAR_DN_1, L_UPAR_UP_1, LC_UPAR_DN_2, LC_UPAR_UP_2 )         ! Output
                  ENDIF

!  Greens (completely new subroutines)
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

                ELSE

                  CALL VLIDORT_GUSER_SOLUTION ( &
                    DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,      & ! Input flags
                    LAYER, FOURIER, IBEAM, LC_NSTK, N_USER_STREAMS, NMOMENTS, & ! Input numbers
                    DO_LAYER_SCATTERING, BEAM_CUTOFF, FLUX_FACTOR,            & ! Input Bookkeeping
                    DFLUX, OMEGA_GREEK, PI_XUP, PI_XUM, PI_X0P,               & ! Input optical/PI/flux
                    UPAR_DN_1, UPAR_UP_1 )                                      ! Output user solutions

                  IF ( DO_ATMOS_LINEARIZATION ) THEN
                    CALL VLIDORT_L_GUSER_SOLUTION ( &
                      DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,  & ! Input Flags
                      LAYER, FOURIER, IBEAM, LC_NSTK, N_USER_STREAMS,       & ! Input numbers
                      NMOMENTS, DO_RTSOL_VARY, NPARAMS_VARY,                & ! Input Numbers
                      DO_LAYER_SCATTERING, BEAM_CUTOFF, FLUX_FACTOR, DFLUX, & ! Input flux/Bookkeeping
                      L_OMEGA_GREEK, PI_XUP, PI_XUM, PI_X0P,                & ! Input optical/PI
                      L_UPAR_DN_1, L_UPAR_UP_1 )                              ! Output
                  ENDIF

                ENDIF

!  end post-processing clause

              END IF
            END IF

!  end layer loop

          END DO

!  Add thermal solutions if flagged
!  ---------------------------------

!    NO modulus on the thermal contribution (Bug fixed 26 January 2010)

          IF ( DO_INCLUDE_THERMEMISS ) THEN
            O1 = 1
            DO N = 1, NLAYERS
             DO I = 1, NSTREAMS_2
! 1/26/10       WUPPER(I,O1,N) = WUPPER(I,O1,N) + pi4*T_WUPPER(I,N)
! 1/26/10       WLOWER(I,O1,N) = WLOWER(I,O1,N) + pi4*T_WLOWER(I,N)
                WUPPER(I,O1,N) = WUPPER(I,O1,N) + T_WUPPER(I,N)
                WLOWER(I,O1,N) = WLOWER(I,O1,N) + T_WLOWER(I,N)
              ENDDO
            ENDDO
          ENDIF

!  Adding the linearized thermal solutions is done later.....

!  2A. Boundary Value problem
!  ==========================

!  Standard case using compressed-band matrices, etc..
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!mick mod 9/19/2017 - replaced DO_REFLECTED_DIRECTBEAM(IBEAM) with DO_LOCALBEAM(IBEAM) in
!                     both BVP_SOLUTION_MASTER & BVPTEL_SOLUTION_MASTER

!       if ( do_write_screen) write(*,*)'bvp solution',ibeam
          IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN

!  Get the BVP solution (regular case) for this beam component
!mick mod 9/19/2017 - replaced DO_REFLECTED_DIRECTBEAM(IBEAM) with DO_LOCALBEAM(IBEAM)

!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19
!      -- Rob fix  4/9/2019  - Major overhaul for adjusted BVP solution

!  1/31/21. Version 2.8.3.  No Change in this calling routine
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

            CALL BVP_SOLUTION_MASTER ( & 
               DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_INCLUDE_SURFEMISS,          & ! Input Surface Flags
               DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_LOCALBEAM(IBEAM),                 & ! Input Surface Flags
               DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX, FOURIER, IBEAM,                   & ! Input illumination flags, indices
               DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,                    & ! Input Water-leaving flags
               LC_NSTK, NSTREAMS, NLAYERS, NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2,     & ! Numbers
               LC_NTT, LC_NSUBD, LC_NSUPD, K_REAL, K_COMPLEX, MUELLER_INDEX,             & ! Numbers, bookkeeping
               TF_MAXITER, TF_CRITERION, FLUX_FACTOR, FLUXVEC, SURFACE_FACTOR,           & ! Input factors/WL control
               QUAD_STRMWTS, ALBEDO, AXBID_F, SURFBB, EMISSIVITY,                        & ! Input surface terms
               SLTERM_ISOTROPIC, SLTERM_F_0, TOAFLUX, BOAFLUX,                           & ! Input Sleave/illumination
               SOLARBEAM_BOATRANS, LOCAL_CSZA, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, & ! Input Direct-flux
               RF_DIRECT_BEAM, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, WUPPER, WLOWER,       & ! Input RTE solutions
               BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                                         & ! Input BVP matrices/pivots
               COL2, SCOL2, TRANS_ATMOS_FINAL, SL_QUADTERM, R2_BEAM, LCON, MCON,         & ! Modified input/output
               STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                                     ! Exception handling

!  Exception handling

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_3 = 'Error return from BVP_SOLUTION_MASTER (Beam solution), '//&
                         'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
              STATUS = VLIDORT_SERIOUS ; RETURN
            ENDIF

!  Telescoped case
!mick mod 9/19/2017 - replaced DO_REFLECTED_DIRECTBEAM(IBEAM) with DO_LOCALBEAM(IBEAM)
!    - Version 2.8, Substantial upgrade to this routine (BRDF with Telescoping)
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.

          ELSE

!  1/31/21. Version 2.8.3.  No Change in this calling routine
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

            CALL BVPTEL_SOLUTION_MASTER ( &
              DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,                   & ! Flags
              DO_LOCALBEAM(IBEAM), FOURIER, IBEAM, LC_NSTK,                & ! Numbers
              NSTREAMS, NLAYERS, NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2, & ! Numbers
              N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,              & ! BVP Tel Control
              N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,                & ! BVP Tel Control
              MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR,            & ! Input Bookkeeping
              ALBEDO, AXBID_F, CUMTRANSDOM, CUMQUADDOM,                    & ! Input Surface inputs
              WLOWER, WUPPER, RF_DIRECT_BEAM,                              & ! Input RTE stuff
              SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_DELT_DISORDS,          & ! Input RTE stuff
              BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT, COLTEL2, SCOL2,      & ! Input BVProblem
              R2_BEAM, LCON, MCON, STATUS_SUB, MESSAGE, TRACE_1 )            ! Output and Status

            IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_2 = 'Error return from BVPTEL_SOLUTION_MASTER, '//'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
              STATUS = VLIDORT_SERIOUS
              RETURN
            ENDIF

          ENDIF

!  New 4/28/19. Diffuse and Direct solar flux at the bottom of the atmosphere
!  here is where you save the solutions you need for the planetary problem
!     SBTERM is already available, TRANSTERM(IB,UM) = TRANSBEAM(IB) * TRNMED_USER(UM) / PIE
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

          IF ( DO_PLANETARY_PROBLEM ) THEN
             N = NLAYERS ;  KO1 = K_REAL(N) + 1
             DO I = 1, NSTREAMS
                DO O1 = 1, LC_NSTK
                   SHOM_R = ZERO ; SHOM_CR = ZERO
                   DO K = 1, K_REAL(N)
                      LXR  = LCON(K,N)*SOLA_XPOS(I,O1,K,N) * T_DELT_EIGEN(K,N)
                      MXR  = MCON(K,N)*SOLB_XNEG(I,O1,K,N)
                      SHOM_R = SHOM_R + LXR + MXR
                   ENDDO
                   DO K = 1, K_COMPLEX(N)
                      K0 = 2*K-2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                      LXR_CR =  LCON(K1,N) * SOLA_XPOS(I,O1,K1,N) - LCON(K2,N) * SOLA_XPOS(I,O1,K2,N)
                      LXR_CI =  LCON(K1,N) * SOLA_XPOS(I,O1,K2,N) + LCON(K2,N) * SOLA_XPOS(I,O1,K1,N)
                      MXR_CR =  MCON(K1,N) * SOLB_XNEG(I,O1,K1,N) - MCON(K2,N) * SOLB_XNEG(I,O1,K2,N)
                      HOM1CR = LXR_CR*T_DELT_EIGEN(K1,N) - LXR_CI*T_DELT_EIGEN(K2,N)
                      SHOM_CR = SHOM_CR + HOM1CR + MXR_CR
                   ENDDO
                   TRANSQUAD(I,O1) = SHOM_R + SHOM_CR + WLOWER(I,O1,N)
                ENDDO
             ENDDO
             TRANSDIRECT = LOCAL_CSZA(NLAYERS,IBEAM) * FLUX_FACTOR * TRANS_SOLAR_BEAM(IBEAM)      ! Direct
             DO O1 = 1, LC_NSTK
                TRANSBEAM(O1,IBEAM) = PI2 * DOT_PRODUCT(TRANSQUAD(1:NSTREAMS,O1),QUAD_STRMWTS(1:NSTREAMS)) ! Diffuse
                TRANSBEAM(O1,IBEAM) = TRANSBEAM(O1,IBEAM) + TRANSDIRECT * FLUXVEC(O1)
             ENDDO
          ENDIF            
           
!  2B. Stokes vector Post Processing
!  =================================

!  upwelling
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.

          IF ( DO_UPWELLING ) THEN

!  VLIDORT COMMENTS:  Direct beam inclusion flag:
!   This now has the DBCORRECTION option: if the DBCORRECTION Flag
!   is set, then we will be doing exact calculations of the reflected
!   directbeam, so we do not need to include it in the Post-processing.
!   However, the direct beam will need to be included in the basic RT
!   solution (the BVP), and this is controlled separately by the
!   DO_REFLECTED_DIRECTBEAM(IBEAM) flags.
!     R. Spurr, RT Solutions, Inc., 19 August 2005.

!Rob fix 4/12/12 - added DO_MSMODE_VLIDORT

!            DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND. &
!                 (DO_REFLECTED_DIRECTBEAM(IBEAM) .AND. .NOT.DO_DBCORRECTION) ) .AND. .NOT.DO_MSMODE_VLIDORT

!Rob Fix 5/27/19. DO_INCLUDE_DIRECTRF same as before (DO_INCLUDE_DIRECTBEAM)
!                 DO_INCLUDE_DIRECTSL similarly constructed, using DO_DBCORRECTION      

            DO_INCLUDE_DIRECTRF = &
             ( DO_UPWELLING .AND. ( DO_LOCALBEAM(IBEAM).AND..NOT.DO_DBCORRECTION) ) .AND. .NOT. DO_MSMODE_VLIDORT
            DO_INCLUDE_DIRECTSL = &
             ( FOURIER.EQ.0 .AND. ( DO_SURFACE_LEAVING .AND..NOT.DO_DBCORRECTION) ) .AND. .NOT. DO_MSMODE_VLIDORT

!  4/28/19. Version 2.8.1. Need to set the User-defined surface leaving, if not set.
!     -- Get adjusted User-term surface-leaving contribution

!  1/31/21. Version 2.8.3. USER_SLTERM_F_0 defined locally, drop FOURIER index.
!  1/31/21. Version 2.8.3. Use post-processing mask.
!  1/31/21. Version 2.8.3. Revised Code, Revision of this condition
!     --> User-term Water-leaving terms need to be calculated from scratch for Fourier zero.
!     --> User-term Water-leaving Terms always need adjustment, all Fourier

            IF ( DO_USER_STREAMS .AND. DO_INCLUDE_DIRECTSL .AND. (DO_WATER_LEAVING .AND..NOT.DO_EXTERNAL_WLEAVE) ) then
               O1 = 1 ; TFACTOR = TRANS_ATMOS_FINAL(IBEAM) * FLUX_FACTOR / DELTA_FACTOR
               IF ( FOURIER.EQ.0 ) then
                  IF ( DO_SL_ISOTROPIC ) THEN
                     SL = SLTERM_ISOTROPIC(O1,IBEAM) * TFACTOR
                     DO LUM = 1, N_PPSTREAMS
                        UM = PPSTREAM_MASK(LUM,IBEAM) ; SL_USERTERM(UM,IBEAM,O1) =  SL
                     ENDDO
                  ELSE
                     DO LUM = 1, N_PPSTREAMS
                        UM = PPSTREAM_MASK(LUM,IBEAM)
                        SL_USERTERM(UM,IBEAM,O1) = USER_SLTERM_F_0(O1,UM,IBEAM) * TFACTOR
                     ENDDO
                  ENDIF
               ELSE
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IBEAM)
                     SL_USERTERM(UM,IBEAM,O1) = SL_USERTERM(UM,IBEAM,O1) * TFACTOR
                  ENDDO
               ENDIF
            ENDIF

!  debug
!       write(*,*)'SLUSERTERM',IBEAM,DO_INCLUDE_DIRECTSL,FOURIER,SL_USERTERM(1,IBEAM,1)
            
!  Now, call the post-processing routine
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    - Version 2.8.1, Control for BOA illumination added, 3/23/19
!    - Version 2.8.1, 4/9/19. DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL added
!    - Version 2.8.1, 4/9/19. RF_USER_DIRECT_BEAM, SL_USERTERM have been added

!  1/31/21. Version 2.8.3. SEVERAL INNOVATIONS and CHANGES
!    -- Input taylor Order TAYLOR_ORDER
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F, SURF_MSSTS_F
!    -- DO_CLASSICAL_SOLUTION flag added to control solution method (Greens vs. Classical)
!    -- Following additional arrays are Inputs and outputs for the Green's function
!           DELTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, ITRANS_USERM, & ! Input Greens function stuff
!           ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, SIGMA_P,                   & ! Input Greens function stuff
!           PMULT_UU, PMULT_UD,  UT_PMULT_UU, UT_PMULT_UD,                       & ! Output Greens function multipliers
!    -- Use post-processing masks N_PPSTREAMS,PPSTREAM_MASK
!    -- (RTS 2/16/21). Use local bookkeeping variable LC_NSTK, as well as NSTOKES

            CALL VLIDORT_UPUSER_INTENSITY ( &
              DO_USER_STREAMS,    DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,            & ! Input flags (RT mode)
              DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,           & ! Input flags (sources)
              DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,   DO_INCLUDE_SURFEMISS,           & ! Input flags (Surface)
              DO_DBCORRECTION,     DO_INCLUDE_DIRECTRF,     DO_INCLUDE_DIRECTSL,           & ! Input flags (Surface)
              DO_LAYER_SCATTERING, DO_MSMODE_VLIDORT,       DO_MSMODE_THERMAL,             & ! Input flags (RT mode)
              DO_TOA_CONTRIBS, DO_INCLUDE_BOAFLUX, DO_CLASSICAL_SOLUTION, DO_MSSTS,        & ! Input flags (RT mode)
              FOURIER, IBEAM, NSTOKES, LC_NSTK, NSTREAMS, NLAYERS, N_USER_LEVELS,          & ! Input numbers (basic)
              N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_UP, MUELLER_INDEX, TAYLOR_ORDER,       & ! Input bookkeeping + levels
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                & ! Input partial-layer control
              FLUX_MULT, BOAFLUX, SURFACE_FACTOR, QUAD_WEIGHTS, QUAD_STRMWTS,              & ! Input Flux and quadrature
              T_DELT_DISORDS,  T_DELT_USERM, T_UTUP_USERM, CUMTRANS,                       & ! Input Transmittances
              ALBEDO, BRDF_F, USER_BRDF_F, SURFBB, EMISSIVITY, USER_EMISSIVITY,            & ! Input Surface BRDF/Emiss.
              K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,                       & ! Input Homog. RTE Soln.
              WLOWER, LCON, MCON, T_WLOWER, LAYER_TSUP_UP, LAYER_TSUP_UTUP,                & ! Input RTE PI and thermal
              DELTAU_VERT, PARTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,          & ! Input Beam for Greens 
              T_UTDN_MUBAR, ITRANS_USERM, ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,        & ! Input Greens function
              RF_USER_DIRECT_BEAM, SL_USERTERM, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2,& ! Input User solutions
              HMULT_1, HMULT_2, EMULT_UP, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, SIGMA_P,  & ! Input multipliers
              PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,                                & ! Output multipliers (Greens)
              STOKES_DOWNSURF, BOA_THTONLY_SOURCE, MS_CONTRIBS_F,                          & ! Output 1 (Auxiliary)
              STOKES_F, CUMSOURCE_UP, LAYER_MSSTS_F, SURF_MSSTS_F )                          ! Output 2 (Main)

!  End do upwelling

          ENDIF

!  Downwelling
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16.
!    - Version 2.8.1,  3/23/19. Introduce Control for including TOA illumination
        
!  1/31/21. Version 2.8.3. SEVERAL INNOVATIONS and CHANGES
!    -- Input taylor Order TAYLOR_ORDER
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F
!    -- DO_CLASSICAL_SOLUTION flag added to control solution method (Greens vs. Classical)
!    -- Following additional arrays are Inputs and outputs for the Green's function
!           DELTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, ITRANS_USERM, & ! Input Greens function stuff
!           ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M, SIGMA_M,                   & ! Input Greens function stuff
!           PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,                        & ! Output Greens function multipliers
!    -- Use N_PPSTREAMS and PPSTREAM_MASK, post-processing mask
!    -- (RTS 2/16/21). Use local bookkeeping variable LC_NSTK, as well as NSTOKES

          IF ( DO_DNWELLING ) THEN
            CALL VLIDORT_DNUSER_INTENSITY ( &
              DO_USER_STREAMS,     DO_OBSERVATION_GEOMETRY, DO_INCLUDE_TOAFLUX,        & ! Input flags (RT mode)
              DO_SOLAR_SOURCES,    DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,      & ! Input flags (sources)
              DO_LAYER_SCATTERING, DO_MSMODE_VLIDORT, DO_CLASSICAL_SOLUTION, DO_MSSTS, & ! Input flags (RT mode)
              FOURIER, IBEAM, NSTOKES, LC_NSTK, NLAYERS, N_USER_LEVELS, N_PPSTREAMS,   & ! Input numbers (basic)
              PPSTREAM_MASK, LEVELMASK_DN, TAYLOR_ORDER,                               & ! Input bookkeeping + levels
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,            & ! Input partial-layer control
              FLUX_MULT, TOAFLUX, T_DELT_USERM, T_UTDN_USERM,                          & ! Input Transmittances, Flux
              K_REAL, K_COMPLEX, LCON, MCON, LAYER_TSUP_DN, LAYER_TSUP_UTDN,           & ! Input RTE Sol + thermal
              DELTAU_VERT, PARTAU_VERT, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,      & ! Input Beam for Greens
              T_UTDN_MUBAR, ITRANS_USERM, ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,    & ! Input Greens function 2   
              UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2, SIGMA_M,                     & ! Input User solutions
              HMULT_1, HMULT_2, EMULT_DN, UT_HMULT_DD, UT_HMULT_DU, UT_EMULT_DN,       & ! Input multipliers
              PMULT_DD, PMULT_DU, UT_PMULT_DU, UT_PMULT_DD,                            & ! Output mutlipliers (Greens)
              STOKES_F, CUMSOURCE_DN, LAYER_MSSTS_F  )                                   ! Main output
          ENDIF

!  mean value (integrated) output
!    - Can also use this to get debug Quadrature output
!    - @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument
!    - Version 2.8, Cleaned-up I/O presentation. 7/8/16. Removed QUAD_OUTPUT flag.
!mick mod 9/19/2017 - replaced DO_INCLUDE_DIRECTBEAM with DO_LOCALBEAM(IBEAM)
!                   - added LEVELS_SOLARTRANS & PARTIALS_SOLARTRANS to input arguments

!   Version 2.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination
!   Version 2.8.1, Control for TOA/BOA illumination added, 3/23/19

!  1/31/21. Version 2.8.3.
!   -- Add flag DO_CLASSICAL_SOLUTION, controls use of Greens function
!   -- Additional inputs BVEC and INITIAL_TRANS, T_UTDN_MUBAR, TAYLOR ORDER, LC_NSTKNSTRM
!   -- Additional inputs ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P foro the Green's function
!   -- Additional outputs are the partial-layer Green's function multipliers UT_GMULT_UP, UT_GMULT_DN
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL VLIDORT_INTEGRATED_OUTPUT ( &
               DO_INCLUDE_MVOUTPUT, DO_LOCALBEAM(IBEAM),                       & ! Input flags
               DO_CLASSICAL_SOLUTION, DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX,  & ! Input flags
               DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,  & ! Input flags
               IBEAM, LC_NSTK, NSTREAMS, NLAYERS, N_USER_LEVELS, N_DIRECTIONS, & ! Input numbers
               LC_NSTKNSTRM, LEVELMASK_UP, LEVELMASK_DN, WHICH_DIRECTIONS,     & ! Input Bookkeeping
               PARTLAYERS_LAYERIDX, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,   & ! Input partial layer control
               TAYLOR_ORDER, FLUX_MULT, TOAFLUX, BOAFLUX, FLUX_FACTOR,         & ! Flux inputs
               FLUXVEC, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, PARTAU_VERT, & ! Quadrature inputs
               BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,         & ! Solar beam Param.
               LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, LOCAL_CSZA,             & ! Solar beam transmittances
               T_DELT_DISORDS, T_UTUP_DISORDS, T_UTDN_DISORDS,                 & ! Discrete ordinate transmittances
               ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P,                       & ! Input Greens function
               K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, LCON, MCON,            & ! Input RTE
               T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, BVEC, WUPPER, WLOWER, & ! Input RTE
               T_WLOWER, T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,            & ! Input Thermal
               MEANST_DIFFUSE, FLUX_DIFFUSE, DNMEANST_DIRECT, DNFLUX_DIRECT,   & ! MAIN output
               UT_GMULT_UP, UT_GMULT_DN )                                        ! Auxiliary Output
          ENDIF

!  Version 2.8 remove GOTO statement, combine with linearization flag
!  Finished this Beam solution, if only Stokes vector is required
!          IF ( DO_SIMULATION_ONLY ) GO TO 4000

!  Clause

          IF ( .NOT. DO_SIMULATION_ONLY ) THEN

!  Step 3. Atmospheric Bulk weighting functions
!  --------------------------------------------

            IF ( DO_COLUMN_LINEARIZATION ) THEN

!   variation index = 0

              VARIATION_INDEX = 0

!  4/9/19. Need to calculate LC_TRANS_ATMOS_FINAL for the water- leaving case
!          Assumes that (proportionally) LC derivatives are in the same ratio as those for the first-guess
!            This is an approximation for the iteration, but is exact for the Gordon result              

              IF ( FOURIER.EQ.0 .AND. DO_WATER_LEAVING ) THEN
                 RATIO = HALF * TRANS_ATMOS_FINAL(IBEAM) / SOLARBEAM_BOATRANS(IBEAM)
                 DO Q = 1, N_TOTALCOLUMN_WFS
                    LC_TRANS_ATMOS_FINAL(IBEAM,Q) = RATIO * LC_SOLARBEAM_BOATRANS(IBEAM,Q)
                 ENDDO
              ENDIF
              
!  3A. Solve the linearized BVP
!  ============================

!  (a) Regular BVP linearization

              IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN

!  Get the linearized BVP solution for this beam component
!    - Version 2.8, Cleaned-up I/O presentation. 7/23/16.
!mick mod 9/19/2017 - replaced DO_REFLECTED_DIRECTBEAM(IBEAM) with DO_LOCALBEAM(IBEAM) in
!                     both LC_BVP_SOLUTION_MASTER & LC_BVPTEL_SOLUTION_MASTER

!  4/9/19. Version 2.8.1, Water-leaving control. Also need LC_TRANS_ATMOS_FINAL

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> Several New Inputs arguments to LC_BVP_SOLUTION_MASTER, LC_BVPTEL_SOLUTION_MASTER : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LC_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE
!           ** BRDF Fourier inputs are defined locally for each Fourier component
!mick fix 1/5/2021 - adjusted ordering of some input vars based on updated LC_BVP_SOLUTION_MASTER subroutine statement
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

                CALL LC_BVP_SOLUTION_MASTER ( &
                   DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_LOCALBEAM(IBEAM),         & ! Flags
                   DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,              & ! Flags
                   DO_CLASSICAL_SOLUTION, DO_LAYER_SCATTERING, FOURIER, IBEAM,             & ! Flags/Indices
                   LC_NSTK, NSTREAMS, NLAYERS, N_TOTALCOLUMN_WFS, TAYLOR_ORDER,            & ! Numbers (Basic)
                   NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2, LC_NTT, LC_NSUBD, LC_NSUPD,   & ! Numbers (derived)
                   MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR, QUAD_STRMWTS,         & ! Bookkeeping
                   DELTAU_SLANT, DELTAU_VERT, L_DELTAU_VERT, RF_DIRECT_BEAM, SL_QUADTERM,      & ! Optical and direct beam
                   BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                   & ! Input, Beam parameterization
                   LC_TRANS_ATMOS_FINAL, LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Input, Linearized Beam param
                   CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, & ! Input, Greens Function
                   SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, BVEC,             & ! Input, Homogeneous/Classical
                   BANDMAT2, IPIVOT, SMAT2, SIPIVOT, ALBEDO, BRDF_F,                 & ! Input, BVP matrices, surface
                   L_T_DELT_EIGEN, L_KEIGEN, L_SOLA_XPOS, L_SOLB_XNEG,               & ! Input, Linearized Homog solutions
                   L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER, LC_BVEC,      & ! Input, Linearized Greens/Thermal/BVEC
                   NCON, PCON, L_WLOWER, L_WUPPER,                                   & ! Output - Linearized Constants + PI
                   STATUS_SUB, MESSAGE, TRACE_1 )                                      ! Exception handling

                IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  TRACE_2 = 'Error return from LC_BVP_SOLUTION_MASTER, Beam solutions '
                  TRACE_3 =  'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
                  STATUS = VLIDORT_SERIOUS ; RETURN
                ENDIF

!  (b) telescoped BVP linearization
!    - Version 2.8, Substantial upgrade to this routine (BRDF with Telescoping)
!    - Version 2.8, Cleaned-up I/O presentation. 7/23/16.

!  1/31/21. Version 2.8.3. Introduce Green's function alternative PI solution method.
!    ==> Several New Inputs arguments to LC_BVP_SOLUTION_MASTER, LC_BVPTEL_SOLUTION_MASTER : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LC_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE
!           ** Green's function arrays CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE
!           ** BRDF Fourier inputs are defined locally for each Fourier component
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

              ELSE

                CALL LC_BVPTEL_SOLUTION_MASTER ( &
                  DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, DO_LOCALBEAM(IBEAM),             & ! Flags
                  DO_CLASSICAL_SOLUTION, DO_LAYER_SCATTERING, FOURIER, IBEAM,                 & ! Flags/Indices
                  TAYLOR_ORDER, LC_NSTK, NSTREAMS, NLAYERS, N_TOTALCOLUMN_WFS,                & ! Basic Control Numbers
                  NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2, NLAYERS_TEL, ACTIVE_LAYERS,       & ! Other Numbers
                  N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,           & ! BVPTel Control
                  MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR, RF_DIRECT_BEAM,           & ! Bookkeeping/Surface
                  DELTAU_SLANT, DELTAU_VERT, L_DELTAU_VERT, T_DELT_DISORDS, L_T_DELT_DISORDS, & ! Optical.Direct/Disords
                  BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                   & ! Beam parameterization
                  SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON, BVEC, WLOWER,          & ! Homogeneous Solutions
                  QUAD_STRMWTS, ALBEDO, BRDF_F,                                          & ! Surface inputs
                  CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                    & ! Input, Greens Function
                  ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,                    & ! Input, Greens Function 
                  L_T_DELT_EIGEN, L_KEIGEN, L_SOLA_XPOS, L_SOLB_XNEG, LC_BVEC,           & ! Linearized Homogeneous/Classical
                  LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,                  & ! Linearized Beam parameterization
                  BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT,                                & ! BVPTel matrices
                  L_WLOWER, L_WUPPER, NCON, PCON,                                        & ! Output solutions
                  STATUS_SUB, MESSAGE, TRACE_1 )                                           ! Exception handling

                IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  TRACE_2 = 'Error from LC_BVPTEL_SOLUTION_MASTER, '//&
                       'Column Jacobians, Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
                  STATUS = VLIDORT_SERIOUS ; RETURN
                ENDIF

              ENDIF

!  New 4/28/19. Diffuse and Direct Linearized solar fluxes at the bottom of the atmosphere
!  here is where you save the solutions you need for the planetary problem
!     SBTERM is already available, TRANSTERM(IB,UM) = TRANSBEAM(IB) * TRNMED_USER(UM) / PIE
!     -- Linearization similar to that in (QUADINTENS_LEVEL_DN), lc_wfatmos routine

!  BugFix 8/12/19. First 2 indices of LC_TRANSBEAM needed to be reversed.
             
              IF ( DO_PLANETARY_PROBLEM ) THEN
                N = NLAYERS ;  KO1 = K_REAL(N) + 1
                DO I = 1, NSTREAMS
                  DO O1 = 1, LC_NSTK
                    DO Q = 1, N_TOTALCOLUMN_WFS
                      SHOM_R = ZERO ; SHOM_CR = ZERO
                      DO K = 1, K_REAL(N)
                        LXR  = LCON(K,N)   *   SOLA_XPOS(I,O1,K,N)
                        LLXR = LCON(K,N)   * L_SOLA_XPOS(I,O1,K,N,Q)
                        MLXR = MCON(K,N)   * L_SOLB_XNEG(I,O1,K,N,Q)
                        NXR  = NCON(K,N,Q) *   SOLA_XPOS(I,O1,K,N)
                        PXR  = PCON(K,N,Q) *   SOLB_XNEG(I,O1,K,N)
                        HOM1 = ( NXR + LLXR ) * T_DELT_EIGEN(K,N) + LXR * L_T_DELT_EIGEN(K,N,Q)
                        HOM2 = PXR + MLXR
                        SHOM_R = SHOM_R + HOM1 + HOM2
                      ENDDO
                      DO K = 1, K_COMPLEX(N)
                        K0 = 2 * K - 2 ; K1 = KO1 + K0 ; K2 = K1  + 1
                        NXR1  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K1,N) - NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K2,N)
                        NXR2  =   NCON(K1,N,Q) *   SOLA_XPOS(I,O1,K2,N) + NCON(K2,N,Q) *   SOLA_XPOS(I,O1,K1,N)
                        PXR1  =   PCON(K1,N,Q) *   SOLB_XNEG(I,O1,K1,N) - PCON(K2,N,Q) *   SOLB_XNEG(I,O1,K2,N)
                        LXR1  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K1,N) - LCON(K2,N) *   SOLA_XPOS(I,O1,K2,N)
                        LXR2  =   LCON(K1,N) *   SOLA_XPOS(I,O1,K2,N) + LCON(K2,N) *   SOLA_XPOS(I,O1,K1,N)
                        LLXR1  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K1,N,Q) - LCON(K2,N) * L_SOLA_XPOS(I,O1,K2,N,Q)
                        LLXR2  =   LCON(K1,N) * L_SOLA_XPOS(I,O1,K2,N,Q) + LCON(K2,N) * L_SOLA_XPOS(I,O1,K1,N,Q)
                        MLXR1  =   MCON(K1,N) * L_SOLB_XNEG(I,O1,K1,N,Q) - MCON(K2,N) * L_SOLB_XNEG(I,O1,K2,N,Q)
                        HOM1CR =   ( NXR1 + LLXR1 ) *   T_DELT_EIGEN(K1,N)   - ( NXR2 + LLXR2 ) *   T_DELT_EIGEN(K2,N)
                        HOM2CR =             LXR1   * L_T_DELT_EIGEN(K1,N,Q) -           LXR2   * L_T_DELT_EIGEN(K2,N,Q)
                        HOM3CR = PXR1 + MLXR1
                        SHOM_CR = SHOM_CR + HOM1CR + HOM2CR + HOM3CR
                      ENDDO
                      SHOM = SHOM_R + SHOM_CR
                      L_TRANSQUAD(I,O1,Q) = L_WLOWER(I,O1,N,Q) + SHOM
                    ENDDO
                  ENDDO
                ENDDO
                DO Q = 1, N_TOTALCOLUMN_WFS
                   L_TRANS_DIRC = LOCAL_CSZA(NLAYERS,IBEAM) * FLUX_FACTOR * LC_SOLARBEAM_ATRANS(IBEAM,Q)   ! Direct
                   DO O1 = 1, LC_NSTK
                      L_TRANS_DIFF = PI2 * DOT_PRODUCT(L_TRANSQUAD(1:NSTREAMS,O1,Q),QUAD_STRMWTS(1:NSTREAMS))    ! Diffuse
!                      LC_TRANSBEAM(IBEAM,O1,Q) = L_TRANS_DIFF + L_TRANS_DIRC * FLUXVEC(O1)   Bug 8/12/19
                      LC_TRANSBEAM(O1,IBEAM,Q) = L_TRANS_DIFF + L_TRANS_DIRC * FLUXVEC(O1)
                   ENDDO
                ENDDO   
              ENDIF

!  3B. Post-processing for the weighting functions
!  ===============================================

!  Streamlined for Version 2.8. 7/22/16
!mick fix 9/19/2017 - replaced IPARTIC with IBEAM in UPUSER_COLUMNWF, DNUSER_COLUMNWF,
!                     & MIFLUX_COLUMNWF in this section

!  Post-processing for Column weighting functions (Upwelling)
!  ----------------------------------------------------------

!  Streamlined for  Version 2.8. 7/22/16
!   4/9/19.  For the Column-Linearized BOA source term, need additional inputs
!             -- Used 2 flags to control direct radiances (reflected-beam, surface-leaving)
!             -- Note use of LC_TRANS_ATMOS_FINAL to go with adjusted water-leaving

!  1/31/21. Version 2.8.3.
!    ==> New input arguments to the post processing (Green's function) : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LC_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE, SIGMA_P
!           ** Green's function arrays ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_UU, PMULT_UD
!    ==> Other Input/Output arguments and other changes
!           ** for the multiple scatter source term linearization, use control flag DO_MSSTS
!           ** LInearizations of MSST functions (LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F), now output
!           ** BRDF Fourier inputs are defined locally for each Fourier component
!    -- (RTS 2/16/21). Use local bookkeeping variable LC_NSTK, as well as NSTOKES

              IF ( DO_UPWELLING ) THEN
                CALL UPUSER_COLUMNWF ( &
                  DO_USER_STREAMS,    DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,                   & ! Input flags (RT mode)
                  DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,                  & ! Input flags (sources)
                  DO_CLASSICAL_SOLUTION, DO_MSMODE_VLIDORT,    DO_INCLUDE_SURFACE,                    & ! Input flags (sources)
                  DO_LAMBERTIAN_SURFACE, DO_DBCORRECTION, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,   & ! Input flags (Surface)
                  DO_LAYER_SCATTERING, DO_MSSTS, FOURIER, IBEAM, NSTOKES, LC_NSTK, NSTREAMS, NLAYERS, & ! Input flags/numbers (basic)
                  N_PPSTREAMS, PPSTREAM_MASK, N_USER_LEVELS, N_TOTALCOLUMN_WFS,                     & ! Input numbers (basic)
                  LEVELMASK_UP, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,       & ! Input Level output control
                  TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,              & ! Input Taylor/Optical
                  FLUX_MULT, QUAD_WEIGHTS, QUAD_STRMWTS, USER_SECANTS, MUELLER_INDEX,               & ! Input Flux/quad/Mueller
                  BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_DISORDS,           & ! Input Beam parameterization
                  T_DELT_USERM, T_UTUP_USERM, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,          & ! Input UTrans/Surface.
                  RF_USER_DIRECT_BEAM, SL_USERTERM, LC_TRANS_ATMOS_FINAL, CUMSOURCE_UP,             & ! Input surf rads/Cumsource
                  K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, LCON, MCON,                & ! Input Homogeneous Solutions
                  ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_UU, PMULT_UD, SIGMA_P,            & ! Input Green's function vars
                  T_WLOWER, UHOM_UPDN, UHOM_UPUP, UPAR_UP_1, UPAR_UP_2, HMULT_1, HMULT_2,           & ! Input Thermal/Homog
                  EMULT_UP, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, UT_PMULT_UU, UT_PMULT_UD,        & ! Input multipliers  
                  LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, L_T_DELT_DISORDS,           & ! Lin BeamParm/DODTrans
                  L_T_DELT_USERM, L_T_UTUP_USERM, L_T_WLOWER, L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP,   & ! Lin Trans/Thermal
                  L_KEIGEN, L_T_DELT_EIGEN, L_SOLA_XPOS, L_SOLB_XNEG, L_UHOM_UPDN, L_UHOM_UPUP,     & ! Lin Homog solutions
                  L_HMULT_1, L_HMULT_2, LC_EMULT_UP, L_UT_HMULT_UU, L_UT_HMULT_UD, LC_UT_EMULT_UP,  & ! Lin multipliers
                  L_UPAR_UP_1, LC_UPAR_UP_2, L_ATERM_SAVE, L_BTERM_SAVE, L_WLOWER, NCON, PCON,      & ! Lin PI Solutions
                  L_BOA_THTONLY_SOURCE, COLUMNWF_F, LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F )               ! Output
              ENDIF

!  Post-processing for Column weighting functions (Downwelling)
!  ------------------------------------------------------------

!  Streamlined for  Version 2.8. 7/22/16

!  1/31/21. Version 2.8.3.
!    ==> New input arguments to the post processing (Green's function) : ==
!           ** Bookkeeping      inputs DO_CLASSICAL_SOLUTION, TAYLOR_ORDER
!           ** Linearization    arrays LC_AVERAGE_SECANT, L_ATERM_SAVE, L_BTERM_SAVE, SIGMA_M
!           ** Green's function arrays ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_DU, PMULT_DD
!    ==> Other Input/Output arguments and other changes
!           ** for the multiple scatter source term linearization, use control flag DO_MSSTS
!           ** LInearizations of MSST functions (LC_LAYER_MSSTS_F), now output
!           ** BRDF Fourier inputs are defined locally for each Fourier component
!    -- (RTS 2/16/21). Use local bookkeeping variable LC_NSTK, as well as NSTOKES
!mick fix 1/5/2021 - added NLAYERS to inputs

              IF ( DO_DNWELLING ) THEN
                CALL DNUSER_COLUMNWF ( &
                  DO_USER_STREAMS,   DO_OBSERVATION_GEOMETRY, DO_CLASSICAL_SOLUTION,                & ! Input flags (RT mode)
                  DO_SOLAR_SOURCES,  DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,                 & ! Input flags (sources)
                  DO_MSMODE_VLIDORT, DO_LAYER_SCATTERING,     DO_MSSTS,  FOURIER, IBEAM,            & ! Input flags (sources)
                  NSTOKES, LC_NSTK, NLAYERS, N_USER_LEVELS, N_TOTALCOLUMN_WFS,                      & ! Input flags/numbers (basic)
                  N_PPSTREAMS, PPSTREAM_MASK,                                                       & ! Input flags/numbers (basic)
                  LEVELMASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,       & ! Input Level output control
                  TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT, FLUX_MULT, USER_SECANTS,   & ! Input Taylor/Opt/Secants
                  BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,           & ! Input Beam param
                  T_DELT_USERM, T_UTDN_USERM, CUMSOURCE_DN,  K_REAL, K_COMPLEX, LCON, MCON,         & ! Input User/Cumsource/BVP
                  ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, PMULT_DU, PMULT_DD, SIGMA_M,            & ! Input Green's function vars
                  UHOM_DNDN, UHOM_DNUP, UPAR_DN_1, UPAR_DN_2, HMULT_1, HMULT_2,                     & ! Homog
                  EMULT_DN, UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, UT_PMULT_DU, UT_PMULT_DD,        & ! Input multipliers  
                  L_T_DELT_USERM, L_T_UTDN_USERM, LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Lin BeamParm/Trans
                  L_KEIGEN, L_UHOM_DNDN, L_UHOM_DNUP, L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN,               & ! Lin Homog/Thermal
                  L_HMULT_1, L_HMULT_2, LC_EMULT_DN, L_UT_HMULT_DU, L_UT_HMULT_DD, LC_UT_EMULT_DN,      & ! Lin multipliers
                  L_UPAR_DN_1, LC_UPAR_DN_2, L_ATERM_SAVE, L_BTERM_SAVE, NCON, PCON,                    & ! Lin PI Solutions
                  COLUMNWF_F, LC_LAYER_MSSTS_F )                                                          ! Output
              ENDIF

! @@@ Rob fix 1/31/11, - added FLUX_FACTOR argument
!    - Version 2.8, Cleaned-up I/O presentation. 7/22/16. Removed QUAD_OUTPUT flag.
!mick mod 9/19/2017 - replaced DO_INCLUDE_DIRECTBEAM with DO_LOCALBEAM(IBEAM)
!                   - added LC_LEVELS_SOLARTRANS & LC_PARTIALS_SOLARTRANS to input arguments

!  1/31/21. Version 2.8.3. Rewritten to incorporate new Linearized Green's function solution
!     -- New arguments include DO_CLASSICAL_SOLUTION, TAYLOR_ORDER, PARTAU_VERT, L_DELTAU_VERT
!     -- New Greens function variables ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, UT_GMULT_UP, UT_GMULT_DN
!     -- Logic changed to perform additional Partial-layer Green function multipliers initially
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

              IF ( DO_INCLUDE_MVOUTPUT ) THEN
                CALL MIFLUX_COLUMNWF ( &
                  DO_INCLUDE_MVOUTPUT, DO_LOCALBEAM(IBEAM), DO_CLASSICAL_SOLUTION,           & ! Input flags
                  DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,             & ! Input flags
                  IBEAM, N_TOTALCOLUMN_WFS, LC_NSTK, NSTREAMS, LC_NSTKNSTRM, NLAYERS,        & ! Input numbers
                  N_USER_LEVELS, N_DIRECTIONS, WHICH_DIRECTIONS, LEVELMASK_UP, LEVELMASK_DN, & ! Level/Dir output control
                  PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,              & ! Partial layer Output control
                  FLUX_MULT, FLUX_FACTOR, FLUXVEC, LOCAL_CSZA, TAYLOR_ORDER,                 & ! Input Flux/Angles/Taylor
                  QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, PARTAU_VERT, L_DELTAU_VERT,      & ! Input quadrature/Optical
                  BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! Beam parameterization
                  T_DELT_DISORDS, T_UTDN_DISORDS, T_UTUP_DISORDS, K_REAL, K_COMPLEX,         & ! DODTrans/Bookkeep
                  SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,            & ! Homog solution
                  ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, UT_GMULT_UP, UT_GMULT_DN,        & ! Green's function solution
                  WUPPER, LCON, MCON, T_WLOWER, T_WUPPER, BOA_THTONLY_SOURCE,                & ! Solar/BVP/Thermal
                  LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,                      & ! Lin Beam parameterization
                  LC_T_UTDN_MUBAR, LC_LEVELS_SOLARTRANS, LC_PARTIALS_SOLARTRANS,             & ! Lin Beam transmittances
                  L_T_DELT_DISORDS, L_T_UTDN_DISORDS, L_T_UTUP_DISORDS, L_KEIGEN,            & ! Lin DODTrans/Keigen
                  L_SOLA_XPOS, L_SOLB_XNEG, L_T_DELT_EIGEN, L_T_UTDN_EIGEN, L_T_UTUP_EIGEN,  & ! Lin Homog solution
                  L_WLOWER, L_WUPPER, NCON, PCON, L_ATERM_SAVE, L_BTERM_SAVE,                & ! Lin PI/Greens Solutons
                  L_UT_T_PARTIC, L_T_WLOWER, L_T_WUPPER, L_BOA_THTONLY_SOURCE,               & ! Lin Thermal solutions
                  MEANST_DIFFUSE_COLWF, DNMEANST_DIRECT_COLWF,                               & ! Lin Output Actinic Fluxes
                  FLUX_DIFFUSE_COLWF, DNFLUX_DIRECT_COLWF )                                    ! Lin Output Regular Fluxes
              ENDIF

!  End atmospheric column weighting functions

            ENDIF

!  Step 4. Surface Reflectance weighting functions
!  -----------------------------------------------

!mick fix 9/6/2012 - added (N_SURFACE_WFS > 0) condition to IF

!  5/24/21. Version 2.8.3. Must add SLEAVEWFS inclusion flag.

            IF ( DO_SURFACE_LINEARIZATION .AND. DO_INCLUDE_SURFACEWF .or. DO_INCLUDE_SLEAVEWFS ) THEN

              IF ( N_SURFACE_WFS > 0 ) THEN

!  4/9/19. Additional code to set the LS derivative of TRANS_ATMOS_FINAL
!     NOT ENABLED YET. Neglect this contribution, so set to zero       

                LS_TRANS_ATMOS_FINAL(IBEAM,1:N_SURFACE_WFS) = ZERO
          
!  Set local flag (same as in Intensity-only case)
!Rob Fix 5/27/19. DO_INCLUDE_DIRECTRF same as before (DO_INCLUDE_DIRECTBEAM)

!                DO_INCLUDE_DIRECTBEAM = &
!                  ( DO_UPWELLING .AND. DO_LOCALBEAM(IBEAM) ) .AND. .NOT. DO_MSMODE_VLIDORT
                 DO_INCLUDE_DIRECTRF = &
                ( DO_UPWELLING .AND. ( DO_LOCALBEAM(IBEAM).AND..NOT.DO_DBCORRECTION) ) .AND. .NOT. DO_MSMODE_VLIDORT

!  Surface WFs; Solve boundary value problem. Regular or telescoped problem
!mick mod 9/19/2017 - replaced DO_REFLECTED_DIRECTBEAM(IBEAM) with DO_LOCALBEAM(IBEAM) in
!                     both SURFACEWF_BVP_SOLUTION & SURFACEWF_BVPTEL_SOLUTION

!   4/9/19. Regular Solution. Introduce variability for adjusted water-leaving transmittance (handle only)

!    -- 1/31/21. Version 2.8.3. LS_BRDF_F, LS_BRDF_F_0 defined locally, drop MAXMOMENTS dimension
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

                IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN
                  CALL SURFACEWF_BVP_SOLUTION  ( FOURIER, &
                      DO_LOCALBEAM(IBEAM), DO_INCLUDE_SURFEMISS, DO_LAMBERTIAN_SURFACE,   & ! Input flags
                      DO_WATER_LEAVING, IBEAM, LC_NSTK, NSTREAMS, NLAYERS, N_SURFACE_WFS, & ! Input
                      LC_NSTKNSTRM, LC_NSTKNSTRM_2, LC_NTT, LC_NSUPD, LC_NSUBD,           & ! Input
                      MUELLER_INDEX, K_REAL, K_COMPLEX, QUAD_STRMWTS,                 & ! Input
                      SURFACE_FACTOR, ATMOS_ATTN, LS_TRANS_ATMOS_FINAL, SL_QUADTERM,  & ! input surface
                      SURFBB, LS_EMISSIVITY, LS_BRDF_F, LS_BRDF_F_0,                  & ! Input
                      T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, WLOWER,                     & ! Input
                      BANDMAT2, IPIVOT, SMAT2, SIPIVOT, LCON, MCON,                   & ! Input
                      NCON_SWF, PCON_SWF, STATUS_SUB, MESSAGE, TRACE_1 )                ! Output&
                ELSE  
                  CALL SURFACEWF_BVPTEL_SOLUTION &
                    ( DO_LOCALBEAM(IBEAM), DO_INCLUDE_SURFACE, FOURIER,               & ! Input Flags/indices
                      IBEAM, LC_NSTK, NSTREAMS, NLAYERS, N_SURFACE_WFS,               & ! Input Numbers
                      LC_NSTKNSTRM, LC_NSTKNSTRM_2,                                   & ! Input
                      N_BVTELMATRIX_SIZE, NLAYERS_TEL, ACTIVE_LAYERS,                 & ! BVPTel Control
                      N_BVTELMATRIX_SUPDIAG, N_BVTELMATRIX_SUBDIAG,                   & ! BVPTel Control
                      MUELLER_INDEX, K_REAL, K_COMPLEX, SURFACE_FACTOR, QUAD_STRMWTS, & ! Input
                      LS_BRDF_F, LS_BRDF_F_0, T_DELT_DISORDS, ATMOS_ATTN,             & ! Inputs
                      T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, WLOWER, LCON, MCON,         & ! Input
                      BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT,                         & ! Input BVPTel
                      NCON_SWF, PCON_SWF, STATUS_SUB, MESSAGE, TRACE_1 )                ! Output
                ENDIF

!  Exception handling

                IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN
                     TRACE_2 = 'Error return from SURFACEWF_BVP_SOLUTION, '// &
                               'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
                  ELSE
                     TRACE_2 = 'Error return from SURFACEWF_BVPTEL_SOLUTION, '// &
                               'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
                  ENDIF
                  STATUS = VLIDORT_SERIOUS ; RETURN
                ENDIF

!  Surface WFs: Postprocessing Master.  [DO_INCLUDE_DIRECTBEAM = .false. (automatic here)]
!  4/9/19. Separate terms from surface (reflected-DB and Sleave). Revised I/O list.

!  1/31/21. Version 2.8.3. Several minor changes to this call.
!    -- Drop LOCAL_UM_START; all UM loops start with 1
!    -- Use the post processing mask system, don't need separate OBSGEOM calculation
!    -- Introduce DO_MSSTS flag, calculate linearized MSSTS output LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F
!    -- BRDF arrays USER_BRDF_F, LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0 all defined locally each Fourier
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others
!mick fix 1/5/2021 - added NSTOKES back to inputs (in addition to LC_NSTK)

                CALL SURFACEWF_POSTPROCESS_MASTER &
                  ( DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_MSSTS,                       & ! Input flags (general)
                    DO_DBCORRECTION, DO_INCLUDE_MVOUTPUT, DO_OBSERVATION_GEOMETRY,               & ! Input flags (general)
                    DO_INCLUDE_SURFEMISS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,               & ! Input flags (thermal)
                    DO_LAMBERTIAN_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,             & ! Input flags (surface)
                    FOURIER, IBEAM, NSTOKES, LC_NSTK, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS, & ! Input Control numbers
                    N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_UP, LEVELMASK_DN,                      & ! Input Level output control
                    PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,                & ! input partlayer control
                    N_DIRECTIONS, WHICH_DIRECTIONS, MUELLER_INDEX,                               & ! Input bookkeeping
                    FLUX_MULT, FLUXVEC, SURFACE_FACTOR, QUAD_WEIGHTS, QUAD_STRMWTS,              & ! input flux/quads
                    SL_USERTERM, ALBEDO, USER_BRDF_F, SURFBB, LS_USER_EMISSIVITY, LS_EMISSIVITY, & ! input Surface stuff
                    LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_TRANS_ATMOS_FINAL,           & ! input Surface stuff
                    T_DELT_DISORDS, T_UTUP_DISORDS, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,    & ! Input trans (User/Disords)
                    K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG, T_DELT_EIGEN,                       & ! Input Homog. solutions
                    T_UTUP_EIGEN, T_UTDN_EIGEN, UHOM_UPDN, UHOM_UPUP, UHOM_DNDN, UHOM_DNUP,      & ! Input Homog. Solutions
                    HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,  UT_HMULT_DU, UT_HMULT_DD,       & ! Input Homog.Multipliers
                    ATMOS_ATTN, STOKES_DOWNSURF, NCON_SWF, PCON_SWF,                             & ! BOA downwelling, Linearized BVP
                    SURFACEWF_F, LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F, MEANST_DIFFUSE_SURFWF, FLUX_DIFFUSE_SURFWF)     ! Output

!  End clause for N_SURFACE_WFS > 0

              ENDIF

! Addition of SLEAVE weighting function,  R. Spurr, 22 August 2012

!    -- 5/24/21. Only works for external Water-leaving inputs.
!    -- 5/24/21. Add INCLUDE_SLEAVEWF flag (similar to DO_INCLUDE_SLEAVING in bvproblem)

              IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS .and. DO_INCLUDE_SLEAVEWFS ) then

!  Direct-beam Flag (repeated here, just so we know!)

                DO_INCLUDE_DIRECTSL = &
                ( DO_UPWELLING .AND. DO_LOCALBEAM(IBEAM).AND..NOT.DO_DBCORRECTION)  .AND. .NOT.DO_MSMODE_VLIDORT

!  4/9/19. Additional code to set the LSSL derivative of TRANS_ATMOS_FINAL
!     NOT ENABLED YET. Neglect this contribution, so set to zero       

                LSSL_TRANS_ATMOS_FINAL(IBEAM,1:N_SLEAVE_WFS) = ZERO
  
!  Rob Fix @@@ 11 Sep 12, Add New Line to the Argument list
!   Cleaned up for Version 2.8 7/22/16
!   4/9/19. Routine completely rewritten

!  1/31/21. Version 2.8.3. Some changes.
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier
!    -- Introduce DO_MSSTS flag (input) and linearized MSST output (LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F)
!    -- Drop LOCAL_UM_START; all UM loops start with 1
!    -- Use the post processing mask system, don't need separate OBSGEOM calculation
!    -- Extension to all Fourier components possible with water-leaving.
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others
!mick fix 1/5/2021 - added NSTOKES back to inputs (in addition to LC_NSTK)

                CALL VLIDORT_LSSL_WFS ( &
                    DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, DO_MSSTS,           & ! input flags (general)
                    DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_LAMBERTIAN_SURFACE,             & ! input flags (general)
                    DO_WATER_LEAVING, DO_SL_ISOTROPIC, DO_INCLUDE_DIRECTSL, NSTOKES,         & ! input flags (sleave)
                    LC_NSTK, NSTREAMS, NLAYERS, N_USER_LEVELS, N_PPSTREAMS, PPSTREAM_MASK,   & ! Input Numbers
                    N_SLEAVE_WFS, N_SURFACE_WFS, N_DIRECTIONS, WHICH_DIRECTIONS,             & ! Input Numbers
                    LC_NTT, LC_NSUBD, LC_NSUPD, LC_NSTKNSTRM, LC_NSTKNSTRM_2,                & ! Input Numbers
                    LEVELMASK_UP, LEVELMASK_DN, FLUX_MULT,                                   & ! Input level out
                    PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,            & ! Input Partial control
                    FOURIER, IBEAM, FLUX_FACTOR, DELTA_FACTOR, SURFACE_FACTOR,               & ! Inputs bookkeeping
                    ALBEDO, USER_BRDF_F, TRANS_ATMOS_FINAL, SL_QUADTERM, SL_USERTERM,        & ! Inputs surface
                    LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,            & ! Inputs Lin Sleave
                    QUAD_WEIGHTS, QUAD_STRMWTS, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Inputs Quads/BVP
                    LSSL_TRANS_ATMOS_FINAL, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,        & ! Inputs Transmittances
                    K_REAL, K_COMPLEX, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,             & ! Input RT solutions
                    SOLA_XPOS, SOLB_XNEG, UHOM_UPDN, UHOM_UPUP, UHOM_DNDN, UHOM_DNUP,        & ! Input  RT solutions
                    HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,    & ! Input Multipliers
                    SURFACEWF_F, LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F,                          & ! Output (Main)
                    MEANST_DIFFUSE_SURFWF, FLUX_DIFFUSE_SURFWF, STATUS_SUB, MESSAGE, TRACE_1 ) ! Output (Flux) and exceptions

! write(*,*)FOURIER,IBEAM,SURFACEWF_F(1:2,1,1,IBEAM,1,1),STOKES_F(1,1,IBEAM,1,1)

!  Exception handling

                IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  TRACE_2 = 'Error return from VLIDORT_LSSL_WFS, '// &
                     'Beam Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
                  STATUS = VLIDORT_SERIOUS ; RETURN
                ENDIF

!  End surface-leaving linearizations

              ENDIF

!  End of surface weighting functions

            ENDIF

!  New, 28 March 2014. Linearization for BLACKBODY
!   Enabled for Version 2.8, 11 July 2016
!  1/31/21. Version 2.8.3. BRDF arrays are defined locally for each Fourier
!   -- (RTS 2/16/21). Use local bookkeeping variables LC_NSTK, and others

            IF ( ( DO_ATMOS_LBBF .OR. DO_SURFACE_LBBF ) .AND. Fourier.EQ.0 ) THEN
              IF ( N_PARTLAYERS .GT. 0 ) THEN
                CALL VLIDORT_LBBF_JACOBIANS_WPARTIALS ( &
                  DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,      & ! Input
                  DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,              & ! Input
                  DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,   & ! input
                  DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,                 & ! input
                  LC_NSTK, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, & ! input
                  NMOMENTS, NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2,        & ! Input
                  LC_NTT, LC_NSUPD, LC_NSUBD,  MUELLER_INDEX,                & ! input
                  N_PARTLAYERS, PARTLAYERS_LAYERIDX,                         & ! Input
                  PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
                  LEVELMASK_UP, LEVELMASK_DN,                                & ! Input
                  USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input
                  QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                  & ! input
                  SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,               & ! input
                  EMISSIVITY, USER_EMISSIVITY, FLUX_MULT,                    & ! input
                  DELTAU_VERT, PARTAU_VERT, OMEGA_GREEK,                     & ! Input
                  T_DELT_DISORDS, T_UTDN_DISORDS, T_UTUP_DISORDS,            & ! input
                  PI_XQP, PI_XQP_PRE, PI_XUP, PI_XUM, SAB, DAB,              & ! input
                  K_REAL, K_COMPLEX, SOLA_XPOS, SOLB_XNEG,                   & ! input
                  T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,                  & ! input
                  T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                  & ! Input
                  BANDMAT2, IPIVOT, SMAT2, SIPIVOT, HMULT_1, HMULT_2,        & ! Input
                  UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,        & ! Input
                  UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                & ! input
                  ABBWFS_JACOBIANS, ABBWFS_FLUXES,                           & ! Output
                  SBBWFS_JACOBIANS, SBBWFS_FLUXES,                           & ! Output
                  STATUS_SUB, MESSAGE, TRACE_1 )                               ! Output
              ELSE
                CALL VLIDORT_LBBF_JACOBIANS_WHOLE (&
                  DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,      & ! Input
                  DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,              & ! Input
                  DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,   & ! input
                  DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,                 & ! input
                  LC_NSTK, NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, & ! input
                  NMOMENTS, NSTREAMS_2, LC_NSTKNSTRM, LC_NSTKNSTRM_2,        & ! Input
                  LC_NTT, LC_NSUPD, LC_NSUBD,  MUELLER_INDEX,                & ! input
                  LEVELMASK_UP, LEVELMASK_DN,                                & ! Input
                  USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,      & ! Input
                  QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                  & ! input
                  SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,               & ! input
                  EMISSIVITY, USER_EMISSIVITY,                               & ! input
                  FLUX_MULT, DELTAU_VERT, OMEGA_GREEK,                       & ! Input
                  T_DELT_DISORDS, T_DELT_USERM,                              & ! input
                  PI_XQP, PI_XQP_PRE, PI_XUP, PI_XUM, SAB, DAB,              & ! input
                  K_REAL, K_COMPLEX, T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG,     & ! input
                  BANDMAT2, IPIVOT, SMAT2, SIPIVOT, HMULT_1, HMULT_2,        & ! Input
                  UHOM_DNDN, UHOM_DNUP, UHOM_UPDN, UHOM_UPUP,                & ! input
                  ABBWFS_JACOBIANS, ABBWFS_FLUXES,                           & ! Output
                  SBBWFS_JACOBIANS, SBBWFS_FLUXES,                           & ! Output
                  STATUS_SUB, MESSAGE, TRACE_1 )                               ! Output
              ENDIF

!  Error handling

              IF ( STATUS_SUB .NE. VLIDORT_SUCCESS ) THEN
                write(CF,'(I2)')FOURIER
                 TRACE_2 = 'Error return from vlidort_lbbf_jacobians, thermal included in PI'// &
                           'Called in VLIDORT_LCS_FOURIER, Fourier # '//CF
                STATUS = VLIDORT_SERIOUS ; RETURN
              ENDIF

!  End of LBBF Jacobians

            ENDIF

!  Continuation point for avoiding weighting functions
!4000   CONTINUE, Removed Version 2.8

!  End Linearization clause

          ENDIF

!  End loop over beam solutions

        ENDIF
      ENDDO

!  debug write

      IF ( DO_DEBUG_WRITE ) THEN
        write(54,'(a,I3,a,100I3)') 'Fourier ',FOURIER, &
                      ' : # complex evalues by layer: ', (K_COMPLEX(LAYER),LAYER=1,NLAYERS)
      ENDIF

!  ######
!  finish
!  ######

      RETURN
      END SUBROUTINE VLIDORT_LCS_FOURIER

!  End Module

      END MODULE vlidort_lcs_fourier1_m

