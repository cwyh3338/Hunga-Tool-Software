
Module TONGA_Forward_Model_V4_m

!  7/28/22. Version 2 uses two measurement sets
!  ============================================
!mick mod 12/9/2022 - To accomodate a Bimodal PSD, added "RTSMie_Master_bimodal_m" & "RTSMie_Master_bimodal_plus_m"

!  Mie modules. Not going to do any linearized Mie yet

   USE RTSMie_parameters_m
   USE RTSMie_Inputs_Def_m
   USE RTSMie_Outputs_Def_m
   USE RTSMie_Lin_Outputs_Def_m
   USE RTSMie_IO_Readwrite_m

   !USE RTSMie_Master_m
   !USE RTSMie_Master_plus_m
   USE RTSMie_Master_bimodal_m
   USE RTSMie_Master_bimodal_plus_m

!  VLIDORT modules
!  7/28/22. Need MASTERS for the Rayleigh calculation

   USE VLIDORT_PARS_m
   USE VLIDORT_IO_DEFS_m
   USE VLIDORT_LIN_IO_DEFS_m

   USE VLIDORT_AUX_m, Only : VLIDORT_WRITE_STATUS
   USE VLIDORT_MASTERS_m
   USE VLIDORT_LCS_MASTERS_m

!  TONGA environment modules
!mick mod 9/22/2022 - added additional TONGA modules "Only" control

   USE TONGA_PARS_m, Only : MPK, E_MAXLAYERS
   USE TONGA_HPTORA_V5_m, Only : Create_HPTORA

private 
public ::  TONGA_Forward_Model_V4

contains

SUBROUTINE  TONGA_Forward_Model_V4 &
          ( TID, NUMIT, MAXA, MA, APARS, LISTA, I, WAV,             & ! Ret Info
            do_ssonly, do_Dump, HPTORA_Global,                      & ! Fwd Model Control, HPTORA inputs
            Do_Initial_Mie, Mie_Inputs_Global,                      & ! Mie I/O
            Bimodal_fraction, Mie_Outputs_Global,                   & ! Mie I/O
            VFixIn_Global, VModIn_Global, VSup_Global,              & ! VLIDORT inputs 
            VLinFixIn_Global, VLinModIn_Global, VLinSup_Global,     & ! VLIDORT inputs 
            SIMRAD, SIMJAC, SO2TAU, RadChk, fail, messages, trace_3 ) ! Output and exception handling

!  Double precision
!mick mod 9/12/2022  - changed HPTORA from "intent(inout)" to "intent(in)"
!mick mod 9/22/2022  - added inputs TID, NUMIT
!                      changed HPTORA to HPTORA_Global
!mick mod 12/9/2022  - To accomodate a Bimodal PSD, made "Mie_Inputs_Global" an array
!                      and added "Bimodal_fraction" and "trace_3"
!mick mod 1/24/2023  - converted "HPTORA_Global" from a scalar to vector of size 2
!mick mod 8/17/2023  - added LISTA
!mick mod 10/28/2024 - added SO2TAU

   implicit none
   !INTEGER, parameter :: mpk = SELECTED_REAL_KIND(15)

!  Retrieval info

   INTEGER  , intent(in) :: TID, NUMIT, MAXA, MA, I
   REAL(mpk), intent(in) :: APARS(MAXA), WAV
   INTEGER  , intent(in) :: LISTA(MAXA)

!  Local flags

   LOGICAL,   intent(in) :: do_ssonly
   LOGICAL,   intent(in) :: do_Dump

!  HPTORA Type structure

   TYPE(Create_HPTORA)             , intent(in) :: HPTORA_Global(2)

!  Mie type structures

   LOGICAL                         , intent(in) :: Do_Initial_Mie
   TYPE(RTSMie_Inputs)             , intent(in) :: Mie_Inputs_Global(2)
   REAL(dp)                        , intent(in) :: Bimodal_fraction
   TYPE(RTSMie_Outputs)            , intent(in) :: Mie_Outputs_Global

!  VLIDORT input structures

   TYPE(VLIDORT_Fixed_Inputs)      , intent(in) :: VFixIn_Global
   TYPE(VLIDORT_Modified_Inputs)   , intent(in) :: VModIn_Global
   TYPE(VLIDORT_Sup_InOut)         , intent(in) :: VSup_Global

!  VLIDORT linearized input structures

   TYPE(VLIDORT_Fixed_LinInputs)   , intent(in) :: VLinFixIn_Global
   TYPE(VLIDORT_Modified_LinInputs), intent(in) :: VLinModIn_Global
   TYPE(VLIDORT_LinSup_InOut)      , intent(in) :: VLinSup_Global

!  Output (simulated Intensity and Jacobians)
!mick mod 9/12/2022 - changed from "intent(inout)" to "intent(out)"

   REAL(mpk), intent(out) :: SIMRAD, SIMJAC(MAXA), SO2TAU, RadChk(2)
   
!  exception handling (already initialized)
!mick mod 9/12/2022 - changed from "intent(inout)" to "intent(out)"

   LOGICAL      , intent(out) :: fail
   CHARACTER*(*), intent(out) :: messages(3), trace_3

!  Local
!  -----

!  HPTORA Type structure
!mick mod 1/24/2023 - converted "HPTORA" from a scalar to vector of size 2

   TYPE(Create_HPTORA)      :: HPTORA(2)

!  Local Mie Input and Output structures
!mick mod 12/9/2022 - To accomodate a Bimodal PSD, made "Mie_Inputs" & "Mie_Lin_Outputs" an array

   TYPE(RTSMie_Inputs)      :: Mie_Inputs(2)
   TYPE(RTSMie_Outputs)     :: Mie_Outputs
   TYPE(RTSMie_Lin_Outputs) :: Mie_Lin_Outputs(2)

!  Local VLIDORT I/O Type structures

   TYPE(VLIDORT_Fixed_Inputs)       :: VFixIn
   TYPE(VLIDORT_Modified_Inputs)    :: VModIn
   TYPE(VLIDORT_Sup_InOut)          :: VSup

   TYPE(VLIDORT_Fixed_LinInputs)    :: VLinFixIn
   TYPE(VLIDORT_Modified_LinInputs) :: VLinModIn
   TYPE(VLIDORT_LinSup_InOut)       :: VLinSup

   TYPE(VLIDORT_Outputs)            :: VLIDORT_Out
   TYPE(VLIDORT_LinOutputs)         :: VLIDORT_LinOut

!  Help

   LOGICAL ::          do_debug_input, OPENFILEFLAG, Mie_Debug
   INTEGER ::          SCENE, NSCENES, NGREEKMAT_ENTRIES, L, N, Q, NLAYERS, NMOMS, NSTOKES
   INTEGER ::          K, CK, GK, RK, CMASK(8), GMASK(8), SMASK(8), RMASK(8)
   DOUBLE PRECISION :: WAER, RAYWT, AERWT
   DOUBLE PRECISION :: AERSCA, AEREXT, MOLEXT, MOLSCA, MOLABS, TOTSCA, TOTEXT, OMEGA, MOMRAT, OMGRAT, GASRATIO
   DOUBLE PRECISION :: PROBLEM_RAY(6,0:2), DEPOL, BETA2, SK, MOM, AERMOM, RAD2, RAD1, JAC2(MAXA)
   CHARACTER*2      :: C2

!  Additional variables for the Mie REF linearization (npars.eq.3)
!mick mod 8/17/2023 - added vector Q

   INTEGER          :: istatus, QA, Q2, QV(MAXA)
   DOUBLE PRECISION :: ext_R, ext_0, ext_w, L_ext_R, L_ext_0, L_ext_w, L_Waer, cutoff
   DOUBLE PRECISION :: L_Aerext, L_totext, L_aersca, L_totsca, L_Raywt, L_Aerwt, L_omega, L_aermom, L_mom

!  8/31/22. Additional variables to deal with SO2 absorption

   DOUBLE PRECISION :: T1, T2, dV, dW, SO2_Xsecs(3), SO2Abs(E_MAXLAYERS), So2Xcs, temp1, temp1sq, SO2
   DOUBLE PRECISION, PARAMETER :: DU_TO_CM2 = 2.6867e+16_mpk
   DOUBLE PRECISION, PARAMETER :: TZERO     = 273.15_mpk

!  Start of code
!  =============

   if ( Mod(I,20).eq.0 ) write(*,*) '  - Forward model progress run, #',I

!  Define some arrays
!  @@@ Rob Fix 6/13/13. CMASK is for the AEROSOLS, RMASK for RAYLEIGH.

   GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)
   RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)    !  This is for Rayleigh
   CMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)    !  This is for aerosol
   !CMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)
   SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)

!  Initialize error flag file

   OPENFILEFLAG = .false.

!  Initialize output

   SIMRAD = 0.0_mpk
   SIMJAC = 0.0_mpk

   fail = .false. ; messages = ' '

!  Copy HPTORA Global type scene data to local type

   HPTORA(1:2) = HPTORA_Global(1:2)

!  Define Jacobian indexing (based on state vector indexing)
!  mick mod 8/17/2023 - added this section

!     LISTA(1) - Aerosol column optical depth
!     LISTA(2) - Aerosol Pk Hgt
!     LISTA(3) - SO2 
!     LISTA(4) - PlumeHWHM
!     LISTA(5) - Fine mode PSD NReal or NImag (only one allowed in a given retrieval)

   qv = LISTA

!  Start scene loop - do RT for two scenes:
!    * 1st for clear scene        (Rayleigh & O3 only             / Rads only)
!    * 2nd for scene with aerosol (Rayleigh, O3, Aero & maybe SO2 / Rads & Jacobians)

   nscenes = 2
   do scene=1,nscenes

!  Rayleigh scattering law

      depol = HPTORA(scene)%RAYLEIGH_DEPOL(I)
      beta2 = ( 1.0_mpk - depol ) / ( 2.0_mpk + depol )
      PROBLEM_RAY = 0.0_mpk
      PROBLEM_RAY(1,0) =  1.0_mpk
      PROBLEM_RAY(4,1) =  3.0_mpk * ( 1.0_mpk - 2.0_mpk*depol ) / (2.0_mpk + depol )
      PROBLEM_RAY(1,2) =  beta2
      PROBLEM_RAY(5,2) =  - DSQRT(6.0_mpk) * beta2
      PROBLEM_RAY(2,2) =  6.0_mpk * beta2

!  Set/reset input data for each scene: copy Global types to local types
!mick mod 9/12/2022 - added VSup, VLinSup

      VFixIn    = VFixIn_Global
      VModIn    = VModIn_Global
      VSup      = VSup_Global

      VLinFixIn = VLinFixIn_Global
      VLinModIn = VLinModIn_Global
      VLinSup   = VLinSup_Global

!  Set number of geometries to 1. Zero the LOCAL ObsGeoms array.

      VModIn%MUserVal%TS_USER_OBSGEOMS_INPUT = ZERO
      VModIn%MUserVal%TS_N_USER_OBSGEOMS = 1

!  Set masking limits

      nstokes = VFixIn%Cont%TS_NSTOKES
      if ( nstokes .eq. 1 ) ngreekmat_entries = 1
      if ( nstokes .eq. 3 ) ngreekmat_entries = 5
      if ( nstokes .eq. 4 ) ngreekmat_entries = 8

!  Add wavelength Diagnostic (just a dummy here). Note: wvl in microns!
!mick mod 9/22/2022  - added input variable "ATMOS_INDEX"

      VFixIn%Optical%TS_ATMOS_WAVELENGTH = WAV * 0.001_mpk
      VFixIn%Optical%TS_ATMOS_INDEX = I

!  Number of layers and heights

      nlayers                                  = HPTORA(scene)%nlayers
      VFixIn%Cont%TS_NLAYERS                   = nlayers
      VFixIn%Chapman%TS_height_grid(0:nlayers) = HPTORA(scene)%Heights(0:nlayers)
      VModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT   = HPTORA(scene)%Heights(nlayers)

!  Initialize optical properties

      VFixIn%Optical%TS_deltau_vert_input    = zero
      VModIn%MOptical%TS_omega_total_input   = zero
      VFixIn%Optical%TS_greekmat_total_input = zero

!  Start scene IF block

      if (scene .eq. 1) then

!  First Call, without aerosols or SO2 (i.e. Rayleigh & O3 only)
!  =============================================================

!  RT control: Rayleigh-only flags, no delta-m

         VModIn%MBool%TS_DO_RAYLEIGH_ONLY   = .true.
         VModIn%MBool%TS_DO_DOUBLE_CONVTEST = .false.
         VModIn%MBool%TS_DO_DELTAM_SCALING  = .false.
         VModIn%MBool%TS_DO_SOLUTION_SAVING = .false.
         VModIn%MBool%TS_DO_BVP_TELESCOPING = .false.

!  Geometry. Locally, use first of Global geometry inputs for first (clear) scene

         VModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3) = VModIn_Global%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3)

!  Initialize # of Greek moments for Rayleigh

         VModIn%MCont%TS_ngreek_moments_input = 2

!  Define optical properties for each layer

         do n = 1, nlayers
            molabs = HPTORA(scene)%GasodO3(I,n) ; molsca = HPTORA(scene)%RayOD(I,n)
            molext = molabs + molsca
            VFixIn%Optical%TS_deltau_vert_input(n)  = molext
            VModIn%MOptical%TS_omega_total_input(n) = molsca / molext
            do K = 1, ngreekmat_entries
               gk = gmask(k) ; rk = rmask(k) ! ; ck = cmask(k)
               do L = 0, 2
                  VFixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(L,n,gk) = PROBLEM_RAY(rk,L)
               enddo
               VFixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0_mpk
            enddo
         enddo

!  DUMP 1. The Rayleigh/Ozone atmosphere.

         if ( do_Dump ) then
            if ( I.eq.1 ) then
               write(1,52)nlayers,VFixIn%Chapman%TS_height_grid(0)
               do n = 1, nlayers
                  write(1,53)n,VFixIn%Chapman%TS_height_grid(n),HPTORA(scene)%LayerAircolumns(n),HPTORA(scene)%LayerTemps(n)
               enddo
            endif
               write(1,53)i,Wav,HPTORA(scene)%RAYLEIGH_XSECS(I), HPTORA(scene)%RAYLEIGH_DEPOL(I)
               do n = 1, nlayers
                  write(1,54)n,VFixIn%Optical%TS_deltau_vert_input(n), VModIn%MOptical%TS_omega_total_input(n), &
                               VFixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:2,n,1) 
            enddo
52          format(I4,1x,1pe20.10)
53          format(I4,1x,1p3e20.10)
54          format(I4,1x,1p5e20.10)
         endif

!  VLIDORT Call (Regular, no Jacobians) --> for "denominator radiance"
!mick mod 9/12/2022 - changed VSup_Global to VSup

         do_debug_input = .false.
         call VLIDORT_MASTER ( do_debug_input, &
              VFixIn, VModIn, VSup, VLIDORT_Out )

!  Exception handling

         call VLIDORT_WRITE_STATUS ( 'TONGA_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )
         if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  .ne. VLIDORT_SUCCESS .or. &
              VLIDORT_Out%Status%TS_STATUS_CALCULATION .ne. VLIDORT_SUCCESS ) then
            messages(1) = 'VLIDORT aborted, first call: look in file TONGA_VLIDORT_Execution.log for diagnostics !!!'
            fail = .true. ; return
         endif

!  Retain first result (the "denominator radiance")

         Rad1 = VLIDORT_Out%Main%TS_STOKES(1,1,1,1)

      else

!  Second Call, with aerosols and maybe SO2
!  ========================================

!  RT control: using aerosols now, set the delta-m flags

         VModIn%MBool%TS_DO_RAYLEIGH_ONLY = .false.
         if ( .not.do_ssonly ) THEN
            VModIn%MBool%TS_DO_DOUBLE_CONVTEST = .true.
            VModIn%MBool%TS_DO_DELTAM_SCALING  = .true.
            VModIn%MBool%TS_DO_SOLUTION_SAVING = .true.
            VModIn%MBool%TS_DO_BVP_TELESCOPING = .true.
         endif

!  Geometry. Locally, use second of Global geometry inputs for second (aerosol) scene

         VModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3) = VModIn_Global%MUserVal%TS_USER_OBSGEOMS_INPUT(2,1:3)

!  Define aerosol scattering properties:

!  (1) Copy Global Mie input, set wavelength. Now we require expansion coefficients for this call.
!mick note 9/22/2022 - Lambda in microns

         Mie_Inputs              = Mie_Inputs_Global
         Mie_Inputs%Lambda       = WAV * 0.001_mpk
         Mie_Inputs%do_Expcoeffs = .true.

!  (2) Set Mie inputs

         if ( Do_Initial_Mie ) then

!      Pre-prepared Mie properties scenario: data just gets copied
!mick mod 8/17/2023 - changed IF condition of Mie write from "MA.eq.2" to ".not.(Retrieve_nReal.or.Retrieve_nImag)"

            Mie_Outputs = Mie_Outputs_Global
            Mie_Debug = .false. ! ; Mie_Debug = .true.
            if ( Mie_Debug ) then
               write(C2,'(I2.2)')I
               if ( .not.(HPTORA(scene)%Retrieve_nReal .or. HPTORA(scene)%Retrieve_nImag) ) &
                 call RTSMie_write_standard ( 'Mie_Debug_FileW'//C2//'.out', Mie_Inputs(1), Mie_Outputs )
               if ( I.eq.1) stop 'Mie Debug'
            endif

         else

!      Freshly-computed Mie properties scenario: compute Mie properties
!        for current wvl "from scratch"

!      Adjust n_real or n_imag, set derivatives flag
!mick mod 8/17/2023 - replaced "MA .eq. 3" IF conditions with "Retrieve_nReal .or. Retrieve_nImag" IF conditions
!                   - replaced APARS(3) with new refrac index location APARS(LISTA(5))

            if ( HPTORA(scene)%Retrieve_nReal .or. HPTORA(scene)%Retrieve_nImag ) then
              if ( HPTORA(scene)%Retrieve_nReal ) then
                 Mie_Inputs%n_real = APARS(LISTA(5)) ; qa = 1
              else if ( HPTORA(scene)%Retrieve_nImag ) then
                 Mie_Inputs%n_imag = APARS(LISTA(5)) ; qa = 2
              endif
            endif

!      Call Mie code

            if ( HPTORA(scene)%Retrieve_nReal .or. HPTORA(scene)%Retrieve_nImag ) then
               !call RTSMie_Master_plus & 
               !  ( Mie_Inputs, Mie_Outputs, Mie_Lin_Outputs,             & ! I/O type structures
               !    fail, istatus, messages(1), messages(2), messages(3) )  ! Exception handling
               call RTSMie_master_bimodal_plus &
                 ( Mie_Inputs(1), Mie_Inputs(2), Bimodal_fraction,      & ! Inputs
                   Mie_Outputs, Mie_Lin_Outputs(1), Mie_Lin_Outputs(2), & ! Outputs
                   fail, istatus, messages, trace_3 )
            else
               !call RTSMie_Master & 
               !  ( Mie_Inputs, Mie_Outputs,                              & ! I/O type structures
               !    fail, istatus, messages(1), messages(2), messages(3) )  ! Exception handling
               call RTSMie_master_bimodal &
                 ( Mie_Inputs(1), Mie_Inputs(2), Bimodal_fraction, & ! Inputs
                   Mie_Outputs, fail, istatus, messages, trace_3 )   ! Outputs
            endif
            if ( fail ) return

!      Use moment cutoff to set limit on number of Mie expansion coefficients

            cutoff = 1.0e-8_mpk
            nmoms = Mie_Outputs%Mie_NCoeffs ; L = 0 ; mom = 1.0e-1_mpk
            do while ( L.lt.nmoms .and. abs(mom).gt.cutoff ) 
               L =  L + 1 ; mom = Mie_Outputs%Mie_Expcoeffs(1,L)
            enddo
            Mie_Outputs%Mie_NCoeffs = L-1
            Mie_Outputs%Mie_Expcoeffs(:,L:max_Mie_angles) = d_zero

!      Mie debug
!mick mod 8/17/2023 - replaced "MA.eq.3" IF condition with "Retrieve_nReal .or. Retrieve_nImag" IF conditions

            Mie_Debug = .false. !; Mie_Debug = .true.
            if ( Mie_Debug ) then
               write(C2,'(I2.2)')I
               if ( HPTORA(scene)%Retrieve_nReal .or. HPTORA(scene)%Retrieve_nImag ) then
                 Call RTSMie_write_extended ( 'Mie_Debug_W'//C2//'.out', .TRUE., Mie_Inputs(1), &
                                               Mie_Outputs, Mie_Lin_Outputs(1) )
               else
                 Call RTSMie_write_standard ( 'Mie_Debug_W'//C2//'.out', Mie_Inputs(1), Mie_Outputs )
               endif
               if ( I.eq.1 ) stop 'Mie Debug'
            endif   

!  End of Mie setup

         endif

!  Initial setup up for VLIDORT RT
!  -------------------------------

!  Jacobian Counters

         q2 = 2

!  Proxies for the Mie calculations

         ext_0 = HPTORA(scene)%Aer_Refw_Extinction
         ext_w = Mie_Outputs%Mie_Bulk(1)
         waer  = Mie_Outputs%Mie_Bulk(3)
         nmoms = Mie_Outputs%Mie_NCoeffs
         ext_R = ext_w / ext_0

!  Linearized proxies for the Mie calculations -
!    Mie output derivatives must be SINGLE normalized
!mick mod 8/17/2023 - replaced "MA .eq. 3" IF conditions with "Retrieve_nReal .or. Retrieve_nImag" IF conditions
!                   - replaced APARS(3) with new refrac index location APARS(LISTA(5))

         if ( HPTORA(scene)%Retrieve_nReal .or. HPTORA(scene)%Retrieve_nImag ) then
           L_ext_0 = HPTORA(scene)%Aer_Refw_L_Extinction
           L_Ext_w = APARS(LISTA(5)) * Mie_Lin_Outputs(1)%LRFE_Mie_Bulk(1,qa)
           L_waer  = APARS(LISTA(5)) * Mie_Lin_Outputs(1)%LRFE_Mie_Bulk(3,qa)
           L_Ext_R = ( L_ext_w - ext_R * L_ext_0 ) / ext_0
         endif

!  Define SO2 Xsecs

         SO2_Xsecs(1) = HPTORA(scene)%SO2c0_xsecs(I) * du_to_cm2
         SO2_Xsecs(2) = HPTORA(scene)%SO2c1_xsecs(I)
         SO2_Xsecs(3) = HPTORA(scene)%SO2c2_xsecs(I)

!  Main optical property setups
!  ----------------------------

!  Initialize # of Greek moments for Rayleigh + aerosol
!mick mod 10/28/2024 - initialize all So2Abs layer values for new SO2TAU calc

         VModIn%MCont%TS_ngreek_moments_input = nmoms

         So2Abs = 0.0_mpk

!  Define optical properties for each layer

         do n = 1, nlayers

            if ( .not.HPTORA(scene)%aerlayerflags(n) ) then

!  Rayleigh layers with Ozone, No SO2:

               molabs = HPTORA(scene)%GasodO3(I,n)
               molsca = HPTORA(scene)%RayOD(I,n)
               molext = molabs + molsca
               VFixIn%Optical%TS_deltau_vert_input(n)  = molext
               VModIn%MOptical%TS_omega_total_input(n) = molsca / molext
               do K = 1, ngreekmat_entries
                  gk = gmask(k) ; rk = rmask(k) ! ; ck = cmask(k)
                  do L = 0, 2
                     VFixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(L,n,gk) = PROBLEM_RAY(rk,L)
                  enddo
                  VFixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0_mpk
               enddo

            else

!  Aerosol layers:

!  Molecular absorption may have the extra SO2 term

               molabs = HPTORA(scene)%GasODO3(I,n) ; So2Abs(n) = 0.0_mpk
               if ( HPTORA(scene)%Include_SO2 ) then
                  temp1     = 0.5_mpk * ( HPTORA(scene)%Temps(n-1) +  HPTORA(scene)%Temps(n) ) - TZERO
                  temp1sq   = temp1 * temp1
                  So2Xcs    = SO2_Xsecs(1) * ( 1.0_mpk + temp1 * SO2_Xsecs(2) + temp1sq * SO2_Xsecs(3) ) 
                  So2Abs(n) = HPTORA(scene)%SO2Du(n) * So2Xcs
                  molabs    = molabs + So2Abs(n)
               endif

!  Deltau and omega, bulk properties

               molsca = HPTORA(scene)%RayOD(I,n)
               molext = molabs + molsca
               aerext = HPTORA(scene)%Loading(n) * Ext_R ; aersca = aerext * waer
               totext = molext + aerext
               totsca = molsca + aersca
               raywt  = molsca / totsca
               aerwt  = 1.0_mpk - Raywt
               omega  = totsca / totext
               VFixIn%Optical%TS_deltau_vert_input(n)  = totext
               VModIn%MOptical%TS_omega_total_input(n) = omega

!  Derivatives of deltau & omega wrt loading (aerosol and, if included, SO2)
!   8/31/22. There was a bug here (trivial effect)
!      formerly, HPTORA(scene)%dLoading_da(n,q) was not multiplied by ext_R at (**)
!mick mod  8/17/2023 - added IF condition for T2
!mick note 8/17/2023 - using "q", not "q-1" for HPTORA(scene)%dSO2Du_da here

!  8/25/23. Here you just need the Include SO2, don't link it with the Retrieve SO2
!         -- otherwise you won't get the correct Jacobians.

               omgrat =  ( (waer/omega) - 1.0_mpk )
!               if ( HPTORA(scene)%Include_SO2 .and. HPTORA(scene)%Retrieve_SO2 ) then
               if ( HPTORA(scene)%Include_SO2 ) then
                  do q = 1, q2
                    T1 = HPTORA(scene)%dLoading_da(n,q) * ext_R   ! (**)
                    T2 = 0.0_mpk ; if ( q.gt.1 ) T2 = HPTORA(scene)%dSO2Du_da(n,q) * So2Xcs
                    dV = T1 + T2
                    dW = omgrat * T1 - T2
                    VLinFixIn%Optical%TS_L_deltau_vert_input(q,n) = dV / totext
                    VLinFixIn%Optical%TS_L_omega_total_input(q,n) = dW / totext
                  enddo
               else
                  do q = 1, q2
                    dV = HPTORA(scene)%dLoading_da(n,q) * ext_R  ! (**)
                    dW = omgrat * dV
                    VLinFixIn%Optical%TS_L_deltau_vert_input(q,n) = dV / totext
                    VLinFixIn%Optical%TS_L_omega_total_input(q,n) = dW / totext
                  enddo
                endif

!  REF derivative of deltau & omega (if flagged)
!mick mod 8/17/2023 - replaced "MA .eq. 3" IF condition with "Retrieve_nReal .or. Retrieve_nImag" IF conditions

               if ( HPTORA(scene)%Retrieve_nReal .or. HPTORA(scene)%Retrieve_nImag ) then
                  L_Aerext = HPTORA(scene)%Loading(n) * L_Ext_R
                  L_totext = L_Aerext
                  L_aersca = L_aerext * waer + aerext * L_waer
                  L_totsca = L_aersca
                  L_Raywt = - Raywt * L_totsca / totsca
                  L_Aerwt = - L_Raywt
                  L_omega = ( L_totsca - omega * L_totext ) / totext
                  VLinFixIn%Optical%TS_L_deltau_vert_input(qv(5),n) = L_totext / totext
                  VLinFixIn%Optical%TS_L_omega_total_input(qv(5),n) = L_omega  / omega
               endif

!  Scattering coefficients

               do K = 1, ngreekmat_entries
                  gk = gmask(k); sk = dble(smask(k)); ck = cmask(k) ; rk = rmask(k)

!  First moments including Rayleigh scattering
!  8/31/22. Small bug. Failed to multiply momrat by the factor ext_R
!mick mod 8/17/2023 - replaced "MA .eq. 3" IF condition with "Retrieve_nReal .or. Retrieve_nImag" IF conditions
!                   - replaced APARS(3) with new refrac index location APARS(LISTA(5))

                  do L = 0, 2
                     Aermom = SK * Mie_Outputs%Mie_Expcoeffs(ck,L)
                     MOM = RAYWT*PROBLEM_RAY(rk,L) + AERWT*Aermom
                     VFixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(L,N,GK) = MOM
                     if ( mom.ne.zero ) then
                        !momrat = ( (Aermom/Mom) - one )          ! Wrong
                        momrat = ( (Aermom/Mom) - one )  * ext_R  ! Right
                        VLinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(1:q2,L,N,GK) = &
                           momrat * HPTORA(scene)%dLoading_da(n,1:q2) * waer / totsca

                        if ( HPTORA(scene)%Retrieve_nReal .or. HPTORA(scene)%Retrieve_nImag ) then
                           L_Aermom = APARS(LISTA(5)) * SK * Mie_Lin_Outputs(1)%LRFE_Mie_Expcoeffs(ck,L,qa)
                           L_mom = L_RAYWT*PROBLEM_RAY(rk,L) + Aermom * L_Aerwt + Aerwt * L_Aermom
                           VLinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(qv(5),L,N,GK) = L_mom / Mom
                        endif
                     endif
                  enddo

!  Remaining moments aerosols only
!   8/31/22. Small bug. Failed to multiply momrat by the factor ext_R. Recast expression for safety.
!mick mod 8/17/2023 - replaced "MA .eq. 3" IF condition with "Retrieve_nReal .or. Retrieve_nImag" IF conditions
!                   - replaced APARS(3) with new refrac index location APARS(LISTA(5))

                  do L = 3, nmoms
                     Aermom = SK * Mie_Outputs%Mie_Expcoeffs(ck,L) ; MOM = AERWT * Aermom
                     VFixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(L,n,gk) = MOM
                     if ( mom.ne.zero ) then
                        momrat = ( (Aermom/Mom) - one )  * ext_R
                        !VLinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(1:q2,L,N,GK) = &
                        !   RAYWT * HPTOR_GlobalA%dLoading_da(n,1:q2) * waer / totsca
                        VLinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(1:q2,L,N,GK) = &
                            momrat * HPTORA(scene)%dLoading_da(n,1:q2) * waer / totsca

                        if ( HPTORA(scene)%Retrieve_nReal .or. HPTORA(scene)%Retrieve_nImag ) then
                           L_Aermom = APARS(LISTA(5)) * SK * Mie_Lin_Outputs(1)%LRFE_Mie_Expcoeffs(ck,L,qa)
                           L_mom =  Aermom * L_Aerwt + Aerwt * L_Aermom
                           VLinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(qv(5),L,N,GK) = L_mom / Mom
                        endif
                     endif
                  enddo

               enddo

!  Normalize

               VFixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0,n,1)           = 1.0_mpk
               VLinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(1:MA,0,n,1) = 0.0_mpk

!  Derivatives of deltau & omega wrt total column of SO2
!   8/25/23. (a) Use both Include and Retrieve SO2 flags
!   8/25/23. (b) Use HPTORA%dSO2Du_dA(n,1) * So2Xcs in place of So2Abs(n)
!           -- we have activated HPTORA%dSO2Du_dA(:,1) wrt total SO2 (see SO2 Loading)
!           -- They should be the same, i.e. 

!               if ( HPTORA(scene)%Retrieve_SO2 ) then
               if ( HPTORA(scene)%Include_SO2 .and. HPTORA(scene)%Retrieve_SO2 ) then
!                 gasratio = So2Abs(n) / totext
                 gasratio = HPTORA(scene)%dSO2Du_dA(n,1) * So2Xcs / totext
                 VLinFixIn%Optical%TS_L_deltau_vert_input(qv(3),n) = + gasratio
                 VLinFixIn%Optical%TS_L_omega_total_input(qv(3),n) = - gasratio
               endif

!  End aerosol layers

            endif

!  End layer loop

         enddo

!write(*,*)'SO2 pkht deriv',scene,HPTORA(scene)%SO2Du(67),HPTORA(scene)%dSO2Du_da(67,2)

!write(100,*)i,sum(SO2ABS(1:nlayers)),SUM(HPTORA(scene)%GasODO3(I,1:nlayers))
!if ( i.eq.139 ) stop


!  Save SO2 total tau
!mick mod 10/28/2024 - added this section

         if ( HPTORA(scene)%Include_SO2 ) SO2TAU = sum(SO2ABS(1:nlayers))

!  DUMP 2. The Aerosol atmosphere.

         if ( do_Dump ) then
            write(2,57)i,nmoms,wav
            do n = 1, nlayers
               write(2,'(L2)')HPTORA(scene)%aerlayerflags(n)
               if ( HPTORA(scene)%aerlayerflags(n) ) then
                 write(2,58)n, VFixIn%Optical%TS_deltau_vert_input(n), VModIn%MOptical%TS_omega_total_input(n),&
                               VFixIn%Optical%TS_greekmat_total_input(0:nmoms,n,1)
               endif
            enddo
57          format(I4,I5,1x,1pe20.10)
58          format(I4,1x,1p150e20.10)
         endif

!if ( i.eq.1) then
!            do n = 1, nlayers
!!               write(252,'(L2)')HPTORA(scene)%aerlayerflags(n)
!               if ( HPTORA(scene)%aerlayerflags(n) ) then
!                 write(252,258)n, VFixIn%Optical%TS_deltau_vert_input(n), VLinFixIn%Optical%TS_L_deltau_vert_input(1:3,n),&
!                                 VModIn%MOptical%TS_omega_total_input(n),VLinFixIn%Optical%TS_L_omega_total_input(1:3,n)
!     write(253,258)n, VFixIn%Optical%TS_greekmat_total_input(1,n,1),VLinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(1:3,1,N,1),&
!                      VFixIn%Optical%TS_greekmat_total_input(2,n,1),VLinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(1:3,2,N,1)
!               endif
!            enddo
!258         format(I4,1x,2(1p4e20.10))
!            stop '252/253'
!endif

!  VLIDORT Call (with Jacobians)
!mick mod 9/12/2022 - changed VSup_Global to VSup and VLinSup_Global to VLinSup

         do_debug_input = .false. !; do_debug_input = .true.
         call VLIDORT_LCS_MASTER ( do_debug_input,          &
              VFixIn,    VModIn,    VSup,    VLIDORT_Out,   &
              VLinFixIn, VLinModIn, VLinSup, VLIDORT_LinOut )

!  Exception handling

         call VLIDORT_WRITE_STATUS ( 'TONGA_VLIDORT_Execution.log', VLIDORT_ERRUNIT, OPENFILEFLAG, VLIDORT_Out%Status )
         if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  .ne. VLIDORT_SUCCESS .or. &
              VLIDORT_Out%Status%TS_STATUS_CALCULATION .ne. VLIDORT_SUCCESS ) then
            messages(1) = 'VLIDORT aborted, second call: look in file TONGA_VLIDORT_Execution.log for diagnostics !!!'
            fail = .true. ; return
         endif

!  Save second results (weighting functions should be unnormalized)
!mick note 2/2/2023 - recall VLIDORT outputs partially normalized WFs/Jacs (i.e. x*dI/dx);
!                     thus, here they are divided by "x" - in this case, the parameter being
!                     retrieved - to obtain the corresponding unnormalized WFs/Jacs (dI/dx).
!mick mod 8/17/2023 - ensure indexing of Jacs of different retrieved pars are the same as
!                     they are in the state vector [Q --> qv(Q)]

         Rad2 = VLIDORT_Out%Main%TS_STOKES(1,1,1,1)
         do Q = 1, MA
           Jac2(qv(Q)) = VLIDORT_LinOut%Col%TS_COLUMNWF(qv(Q),1,1,1,1) / APARS(qv(Q))

!mick check
!if (tid == 3) then
!  write(*,*)
!  write(*,*) 'at Jac check:'
!  write(*,*) 'I = ',I,' Q = ',Q
!  write(*,*) 'VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,1,1,1) = ',VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,1,1,1)
!  write(*,*) 'APARS(Q) = ',APARS(Q)
!  write(*,*) 'Jac2(Q)  = ',Jac2(Q)
!  write(*,*)
!endif

         enddo

!if (i == 139) then
!  write(*,*) 'at Jac check stop'
! stop
!endif

!  End scene IF block

      endif

!  End scene loop

   enddo

!  write(101,*) I, wav, rad1, rad2, rad2/rad1
!  write(*,*)   I, wav, rad1, rad2, rad2/rad1

!  Interpret output (SIMJAC weighting functions should be unnormalized)
!mick mod 8/17/2023 - ensure indexing of Jacs of different retrieved pars are the same as
!                     they are in the state vector [Q --> qv(Q)]

   SIMRAD = Rad2 / Rad1
   do Q = 1, MA
     SIMJAC(qv(Q)) = Jac2(qv(Q)) / Rad1
   enddo

   !For retrieval test purposes
   RadChk(1) = Rad1 ; RadChk(2) = Rad2

!  Done

   return
END SUBROUTINE TONGA_Forward_Model_V4

!  End module
 
End Module TONGA_Forward_Model_V4_m
