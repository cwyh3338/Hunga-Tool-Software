module TONGA_RTInitialize_V2_m

!  Initialization for two radiance data sets now,
!  Use modules

   USE VLIDORT_PARS_m
   USE VLIDORT_IO_DEFS_m
   USE VLIDORT_LIN_IO_DEFS_m

   USE VLIDORT_INPUTS_m, Only   : VLIDORT_BRDF_Sup_Init, VLIDORT_SLEAVE_Sup_Init, VLIDORT_SS_Sup_Init
   USE VLIDORT_L_INPUTS_m, Only : VLIDORT_BRDF_LinSup_Init, VLIDORT_SLEAVE_LinSup_Init, VLIDORT_SS_LinSup_Init

private
public  :: TONGA_RTInitialize_V2

contains

subroutine TONGA_RTInitialize_V2 &
      ( SZA, VZA, SAZ, VAZ, nstokes, nstreams, nfinelayers, do_ssonly, eradius, albedo, npars, & ! Inputs
        VFixIn_Global, VModIn_Global, VSup_Global, VLinFixIn_Global, VLinModIn_Global, VLinSup_Global, fail, message ) 

!  double precision

   implicit none
   integer, parameter :: mpk = SELECTED_REAL_KIND(15)

!  top level inputs

   real(mpk), intent(in) :: SZA(2), VZA(2), SAZ(2), VAZ(2), eradius, albedo
   logical  , intent(in) :: do_ssonly
   integer  , intent(in) :: nstokes, nstreams, nfinelayers, npars

!  VLIDORT Global input structures, to be initialized

   TYPE(VLIDORT_Fixed_Inputs)   , intent(inout)   :: VFixIn_Global
   TYPE(VLIDORT_Modified_Inputs), intent(inout)   :: VModIn_Global
   TYPE(VLIDORT_Sup_InOut)      , intent(inout)   :: VSup_Global

!  VLIDORT Global linearized input structures

   TYPE(VLIDORT_Fixed_LinInputs)   , intent(inout)  :: VLinFixIn_Global
   TYPE(VLIDORT_Modified_LinInputs), intent(inout)  :: VLinModIn_Global
   TYPE(VLIDORT_LinSup_InOut)      , intent(inout)  :: VLinSup_Global

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local 

   real(mpk) :: RAZ
   integer   :: k

!  ##################################################################
!  ##################################################################
!      V L I D O R T    I N I T I A L I Z A T I O N   S E C T I O N
!  ##################################################################
!  ##################################################################

!  Intialize exception handling

   fail = .FALSE.
   message = ' '

!  Initialize type structures

!  Structure 1: Fixed Boolean inputs
!  ----------------------------------

!  Set MS-mode operations in VLIDORT, depenidng on first-order choices
!mick fix 8/15/2022 - initialized DO_DNWELLING

   if ( do_ssonly ) then
      VFixIn_Global%Bool%TS_DO_FULLRAD_MODE     = .FALSE.
   else
      VFixIn_Global%Bool%TS_DO_FULLRAD_MODE     = .TRUE.
   endif

   VFixIn_Global%Bool%TS_DO_THERMAL_EMISSION    = .FALSE.
   VFixIn_Global%Bool%TS_DO_SURFACE_EMISSION    = .FALSE.

   VFixIn_Global%Bool%TS_DO_PLANE_PARALLEL      = .FALSE.

   VFixIn_Global%Bool%TS_DO_UPWELLING           = .TRUE.
   VFixIn_Global%Bool%TS_DO_DNWELLING           = .FALSE.

   VFixIn_Global%Bool%TS_DO_LAMBERTIAN_SURFACE  = .TRUE.
   VFixIn_Global%Optical%TS_LAMBERTIAN_ALBEDO   = albedo

!  These are specialist operations, do not need them

   VFixIn_Global%Bool%TS_DO_TOA_CONTRIBS        = .FALSE.
   VFixIn_Global%Bool%TS_DO_SPECIALIST_OPTION_1 = .FALSE.
   VFixIn_Global%Bool%TS_DO_SPECIALIST_OPTION_2 = .FALSE.
   VFixIn_Global%Bool%TS_DO_SPECIALIST_OPTION_3 = .FALSE.

!  Not needed

   VFixIn_Global%Bool%TS_DO_SURFACE_LEAVING     = .FALSE.
   VFixIn_Global%Bool%TS_DO_SL_ISOTROPIC        = .FALSE.
   VFixIn_Global%Bool%TS_DO_FLUORESCENCE        = .FALSE.
   VFixIn_Global%Bool%TS_DO_WATER_LEAVING       = .FALSE.

   VFixIn_Global%Bool%TS_DO_WLADJUSTED_OUTPUT = .FALSE.
   VFixIn_Global%Bool%TS_DO_TOA_ILLUMINATION  = .FALSE.
   VFixIn_Global%Bool%TS_DO_BOA_ILLUMINATION  = .FALSE.
   VFixIn_Global%Bool%TS_DO_TF_ITERATION      = .FALSE.
   VFixIn_Global%Bool%TS_DO_ALBTRN_MEDIA      = .FALSE.
   VFixIn_Global%Bool%TS_DO_PLANETARY_PROBLEM = .FALSE.

!  1/31/21. Version 2.8.3, New variables must be set by hand

   VFixIn_Global%Cont%TS_ASYMTX_TOLERANCE     = 1.0d-12
   VFixIn_Global%Bool%TS_DO_MSSTS             = .FALSE.
   VFixIn_Global%Bool%TS_DO_FOURIER0_NSTOKES2 = .TRUE.

!  Structure 2: Modified Boolean inputs
!  ------------------------------------
!mick fix 8/15/2022 - initialized DO_CLASSICAL_SOLUTION & DO_EXTERNAL_WLEAVE

!  Computational mode

   VModIn_Global%MBool%TS_DO_CLASSICAL_SOLUTION  = .TRUE.

!  First-order control

   VModIn_Global%MBool%TS_DO_FOCORR              = .TRUE.
   VModIn_Global%MBool%TS_DO_FOCORR_NADIR        = .FALSE.
   VModIn_Global%MBool%TS_DO_FOCORR_OUTGOING     = .TRUE.

   VModIn_Global%MBool%TS_DO_FOCORR_EXTERNAL     = .FALSE. 
   VModIn_Global%MBool%TS_DO_SSCORR_TRUNCATION   = .FALSE.
   VModIn_Global%MBool%TS_DO_SSCORR_USEFMAT      = .FALSE.

!  Geometry

   VModIn_Global%MBool%TS_DO_SOLAR_SOURCES       = .TRUE.
   VModIn_Global%MBool%TS_DO_REFRACTIVE_GEOMETRY = .FALSE.
   VModIn_Global%MBool%TS_DO_CHAPMAN_FUNCTION    = .TRUE.

!  Need scaling and aerosols
!   7/28/22. Settings for Rayleigh-only and Dm Scaling will change, so just initialize here.

   VModIn_Global%MBool%TS_DO_RAYLEIGH_ONLY       = .FALSE.
   VModIn_Global%MBool%TS_DO_DOUBLE_CONVTEST     = .FALSE.
   VModIn_Global%MBool%TS_DO_DELTAM_SCALING      = .FALSE.
   VModIn_Global%MBool%TS_DO_SOLUTION_SAVING     = .FALSE.
   VModIn_Global%MBool%TS_DO_BVP_TELESCOPING     = .FALSE.

!  Always need this flag

   VModIn_Global%MBool%TS_DO_USER_VZANGLES       = .TRUE.

!  No thermal and flux

   VModIn_Global%MBool%TS_DO_ADDITIONAL_MVOUT    = .FALSE.
   VModIn_Global%MBool%TS_DO_MVOUT_ONLY          = .FALSE.
   VModIn_Global%MBool%TS_DO_THERMAL_TRANSONLY   = .FALSE.

!  No surface-leaving

   VModIn_Global%MBool%TS_DO_EXTERNAL_WLEAVE     = .FALSE.

!  Structure 3: Fixed control inputs
!  ---------------------------------

!  Set

   VFixIn_Global%Cont%TS_TAYLOR_ORDER = 3
   VFixIn_Global%Cont%TS_NSTOKES      = nstokes
   VFixIn_Global%Cont%TS_NSTREAMS     = nstreams

!  Check this input

   if ( VFixIn_Global%Cont%TS_NSTREAMS .gt. MAXSTREAMS ) then
       message = 'Number of discrete ordinates too large, exceeds VLIDORT dimension'
       fail = .TRUE. ; return
   endif

!  Initialize

   VFixIn_Global%Cont%TS_NLAYERS          = 0
   VFixIn_Global%Cont%TS_N_THERMAL_COEFFS = 0

!  Set 

   VFixIn_Global%Cont%TS_NFINELAYERS      = nfinelayers
   VFixIn_Global%Cont%TS_VLIDORT_ACCURACY = 0.00001d0

!  Not used

   VFixIn_Global%Cont%TS_NLAYERS_NOMS     = 0
   VFixIn_Global%Cont%TS_NLAYERS_CUTOFF   = 0
   VFixIn_Global%Cont%TS_TF_MAXITER       = 0
   VFixIn_Global%Cont%TS_TF_CRITERION     = ZERO
   VFixIn_Global%Cont%TS_TOA_ILLUMINATION = ZERO
   VFixIn_Global%Cont%TS_BOA_ILLUMINATION = ZERO

!  Structure 4: Modified control inputs
!  ------------------------------------

   VModIn_Global%MCont%TS_NGREEK_MOMENTS_INPUT = 0

!  Structure 5: Beam inputs
!  ------------------------

!  Flux Factor. Set to 1.0 if you want sun-normalized only.
!  Set later with solar spectrum (if used).

   VFixIn_Global%Sunrays%TS_FLUX_FACTOR = ONE

!  Structures 5/6: SZA and User Value inputs
!  -----------------------------------------

!  Observational geometry input. 7/28/22. Now 2 geometries

   VModIn_Global%MBool%TS_DO_OBSERVATION_GEOMETRY = .TRUE.
   VModIn_Global%MUserVal%TS_USER_RELAZMS         = ZERO
   VModIn_Global%MSunrays%TS_SZANGLES             = ZERO
   VModIn_Global%MUserVal%TS_USER_OBSGEOMS_INPUT  = ZERO
   VModIn_Global%MUserVal%TS_USER_VZANGLES_INPUT  = ZERO

   VModIn_Global%MUserVal%TS_N_USER_OBSGEOMS      = 2
   VModIn_Global%MSunrays%TS_N_SZANGLES           = 2
   VModIn_Global%MUserVal%TS_N_USER_RELAZMS       = 2
   VModIn_Global%MUserVal%TS_N_USER_VZANGLES      = 2

   do K = 1, 2
     RAZ = SAZ(K) + 180.0_mpk - VAZ(K)
     if(RAZ .lt. -180.0_mpk)RAZ = RAZ + 360.0_mpk
     if(RAZ .gt. +180.0_mpk)RAZ = RAZ - 360.0_mpk
     RAZ = ABS(RAZ)
     VModIn_Global%MUserVal%TS_USER_OBSGEOMS_INPUT(K,1)  = SZA(K)
     VModIn_Global%MUserVal%TS_USER_OBSGEOMS_INPUT(K,2)  = VZA(K)
     VModIn_Global%MUserVal%TS_USER_OBSGEOMS_INPUT(K,3)  = RAZ
     VModIn_Global%MSunrays%TS_SZANGLES(K)               = SZA(K)
     VModIn_Global%MUserVal%TS_USER_VZANGLES_INPUT(K)    = VZA(K)
     VModIn_Global%MUserVal%TS_USER_RELAZMS(K)           = RAZ
   enddo

!stop 'at USER_OBSGEOMS_INPUT check stop'

!  Levels

   VFixIn_Global%UserVal%TS_N_USER_LEVELS = 1           ! One level
   VModIn_Global%MUserVal%TS_USER_LEVELS  = ZERO        ! Zero for now

!  Geometry specification height = Bottom of height grid (set later, zeroed here)

   VModIn_Global%MUserval%TS_GEOMETRY_SPECHEIGHT = ZERO

!  Not needed

   VModIn_Global%MBool%TS_DO_DOUBLET_GEOMETRY = .FALSE.
   VModIn_Global%MUserVal%TS_N_USER_DOUBLETS  = 0
   VModIn_Global%MUserVal%TS_USER_DOUBLETS    = ZERO

!  Structure 7: Fixed Chapman function inputs
!  ------------------------------------------

!  Set later after PTH profiles have been formed, Zeroed here

   VFixIn_Global%Chapman%TS_HEIGHT_GRID       = ZERO
   
!  These are not required
   
   VFixIn_Global%Chapman%TS_PRESSURE_GRID     = ZERO
   VFixIn_Global%Chapman%TS_TEMPERATURE_GRID  = ZERO
   VFixIn_Global%Chapman%TS_FINEGRID          = 0

!  Earth radius

   VModIn_Global%MChapman%TS_EARTH_RADIUS     = eradius
   VFixIn_Global%Chapman%TS_RFINDEX_PARAMETER = ZERO

!  Structure 8: Fixed optical inputs
!  ---------------------------------

!  These are set later, zeroed here

   VFixIn_Global%Optical%TS_DELTAU_VERT_INPUT     = ZERO
   VModIn_Global%MOptical%TS_OMEGA_TOTAL_INPUT    = ZERO
   VFixIn_Global%Optical%TS_GREEKMAT_TOTAL_INPUT  = ZERO

!  Initialized here, but not used

   VFixIn_Global%Optical%TS_FMATRIX_UP            = ZERO
   VFixIn_Global%Optical%TS_FMATRIX_DN            = ZERO

!  No thermal values here

   VFixIn_Global%Optical%TS_THERMAL_BB_INPUT      = ZERO
   VFixIn_Global%Optical%TS_SURFACE_BB_INPUT      = ZERO

!   Atmospheric wavelength & index, set later on, zeroed here.
!mick mod 9/22/2022 - added input variable "ATMOS_INDEX"

   VFixIn_Global%Optical%TS_ATMOS_WAVELENGTH     = ZERO
   VFixIn_Global%Optical%TS_ATMOS_INDEX          = 0

!  Structure 9: Fixed write inputs
!  -------------------------------

!  No debug or write-to-file options

   VFixIn_Global%Write%TS_DO_DEBUG_WRITE          = .FALSE.

   VFixIn_Global%Write%TS_DO_WRITE_INPUT          = .FALSE.
   VFixIn_Global%Write%TS_INPUT_WRITE_FILENAME    = ' '

   VFixIn_Global%Write%TS_DO_WRITE_SCENARIO       = .FALSE.
   VFixIn_Global%Write%TS_SCENARIO_WRITE_FILENAME = ' '

   VFixIn_Global%Write%TS_DO_WRITE_FOURIER        = .FALSE.
   VFixIn_Global%Write%TS_FOURIER_WRITE_FILENAME  = ' '

   VFixIn_Global%Write%TS_DO_WRITE_RESULTS        = .FALSE.
   VFixIn_Global%Write%TS_RESULTS_WRITE_FILENAME  = ' '

!  Initialize linearized inputs. 2 parameters only.
!  ------------------------------------------------
!mick fix 8/15/2022 - initialized DO_ATMOS_LBBF, DO_SURFACE_LBBF,
!                                 DO_SLEAVE_WFS, N_SLEAVE_WFS, &
!                                 the 5 Lin Optical Inputs

   VLinModIn_Global%MCont%TS_DO_COLUMN_LINEARIZATION  = .TRUE.
   VLinModIn_Global%MCont%TS_DO_PROFILE_LINEARIZATION = .FALSE.
   VLinModIn_Global%MCont%TS_DO_ATMOS_LINEARIZATION   = .TRUE.
   VLinModIn_Global%MCont%TS_DO_SURFACE_LINEARIZATION = .FALSE.
   VLinModIn_Global%MCont%TS_DO_SIMULATION_ONLY       = .FALSE.
   VLinModIn_Global%MCont%TS_DO_LINEARIZATION         = .TRUE.

!  Not needed

   VLinModIn_Global%MCont%TS_DO_ATMOS_LBBF     = .FALSE.
   VLinModIn_Global%MCont%TS_DO_SURFACE_LBBF   = .FALSE.
   VLinModIn_Global%MCont%TS_DO_SLEAVE_WFS     = .FALSE.

!  Jacobian names

   VLinFixIn_Global%Cont%TS_COLUMNWF_NAMES(1)  = 'Jacobian w.r.t. Aerosol TotalOD'
   VLinFixIn_Global%Cont%TS_COLUMNWF_NAMES(2)  = 'Jacobian w.r.t. Aerosol Peakhgt'
   VLinFixIn_Global%Cont%TS_COLUMNWF_NAMES(3)  = ' '
   VLinFixIn_Global%Cont%TS_PROFILEWF_NAMES    = ' '

!  Layer & Jacobian/weighting function control

   VLinFixIn_Global%Cont%TS_LAYER_VARY_FLAG    = .TRUE.
   VLinFixIn_Global%Cont%TS_LAYER_VARY_NUMBER  = npars
   VLinFixIn_Global%Cont%TS_N_TOTALCOLUMN_WFS  = npars
   VLinFixIn_Global%Cont%TS_N_TOTALPROFILE_WFS = 0
   VLinFixIn_Global%Cont%TS_N_SURFACE_WFS      = 0
   VLinFixIn_Global%Cont%TS_N_SLEAVE_WFS       = 0

!  Optical inputs - set later, zeroed here

   VLinFixIn_Global%Optical%TS_L_DELTAU_VERT_INPUT    = ZERO
   VLinFixIn_Global%Optical%TS_L_OMEGA_TOTAL_INPUT    = ZERO
   VLinFixIn_Global%Optical%TS_L_GREEKMAT_TOTAL_INPUT = ZERO

!  Initialized here, but not used

   VLinFixIn_Global%Optical%TS_L_FMATRIX_UP = ZERO
   VLinFixIn_Global%Optical%TS_L_FMATRIX_DN = ZERO

!  Zero the supplemental variables
!  -------------------------------

!  Initialization now done with routines

   CALL VLIDORT_BRDF_Sup_Init   ( VSup_Global )
   CALL VLIDORT_SLEAVE_Sup_Init ( VSup_Global )
   CALL VLIDORT_SS_Sup_Init     ( VSup_Global )

   CALL VLIDORT_BRDF_LinSup_Init   ( VLinSup_Global )
   CALL VLIDORT_SLEAVE_LinSup_Init ( VLinSup_Global )
   CALL VLIDORT_SS_LinSup_Init     ( VLinSup_Global )

!  Done

   return
end subroutine TONGA_RTInitialize_V2

!  End module

end module TONGA_RTInitialize_V2_m

