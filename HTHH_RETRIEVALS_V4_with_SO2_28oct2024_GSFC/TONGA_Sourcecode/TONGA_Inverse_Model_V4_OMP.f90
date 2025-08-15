
Module TONGA_Levenberg_marquardt_V4_m

!  7/28/22. Version 2 uses ratios of radiances
!  ===========================================
!mick mod 12/9/2022 - To accomodate a Bimodal PSD, added "RTSMie_Master_bimodal_m" & "RTSMie_Master_bimodal_plus_m"

!  Mie modules

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

   !USE VLIDORT_PARS_m
   USE VLIDORT_IO_DEFS_m
   USE VLIDORT_LIN_IO_DEFS_m

!  TONGA environment modules
!mick mod 9/22/2022 - added additional TONGA modules "Only" control
!                   - added TONGA_AUX_m
!                   - added module TONGA_Create_HPTORA_V2a_m

   USE TONGA_PARS_m, Only : MPK, ZERO, HALF, ONE
   USE TONGA_HPTORA_V5_m, Only : Create_HPTORA
   USE TONGA_Create_HPTORA_V5_m, Only : TONGA_Write_HPTORA_V5
   USE TONGA_AUX_m, Only : ASYMTX

!  Aerosol Loading module (changes the parameters)

   USE TONGA_Aerosol_Loading_m, Only : Tonga_Aerosol_Loading

!  SO2 Loading module

   USE TONGA_SO2_Loading_m, Only : Tonga_SO2_Loading

!  Forward model routine (7/28/22. Use ratios)

   USE TONGA_Forward_Model_V4_m, Only : TONGA_Forward_Model_V4

private :: GAUSSJ, TONGA_MRQCOF, COVSRT
public  :: TONGA_MRQMIN

contains

SUBROUTINE TONGA_MRQMIN ( DO_TIMING, NTHREADS,                                        & ! Forward model OMP inputs
                          NUMIT, MAXDATA, MAXA, NDATA, MA, WAVS, RADS, SIGS,          & ! Forward model inputs
                          do_ss, do_Dump, HPTORA_Global,                              & ! Forward model inputs + HPTORA
                          do_Initial_Mie, Initial_Mie_File,                           & ! Fwd model Mie inputs
                          Mie_Inputs_Global, Bimodal_fraction,                        & ! Fwd model Mie inputs
                          VFixIn_Global, VModIn_Global, VSup_Global,                  & ! Forward model VLIDORT inputs
                          VLinFixIn_Global, VLinModIn_Global, VLinSup_Global,         & ! Forward model VLIDORT inputs 
                          APARS, LISTA, COVAR, ALPHA, CHISQ, ALAMDA, SIMRADS, SIMJAC, & ! Inverse model input/output
                          OCHISQ, ATRY, BETA, DA, MFIT, SO2TAU, InvMod_Time, fail, messages, trace_3 ) ! Inv model output + status

!mick mod 9/12/2022  - added inputs DO_TIMING, NTHREADS, NUMIT
!mick mod 1/10/2023  - added output SIMJAC
!mick mod 10/28/2024 - added output SO2TAU

   implicit none
   !INTEGER, parameter :: mpk = SELECTED_REAL_KIND(15)

!  OpenMP Control & Timing

   LOGICAL, intent(in)      :: DO_TIMING
   INTEGER, intent(in)      :: NTHREADS

!  Arguments

   LOGICAL  , intent(inout) :: do_ss, do_Dump
   INTEGER  , intent(in)    :: MAXDATA, MAXA, NDATA, MA
   REAL(mpk), intent(in)    :: WAVS(MAXDATA), RADS(MAXDATA), SIGS(MAXDATA)

   INTEGER  , intent(inout) :: NUMIT, LISTA(MAXA), MFIT
   REAL(mpk), intent(inout) :: ALPHA(MAXA,MAXA), COVAR(MAXA,MAXA), CHISQ, ALAMDA, APARS(MAXA)
   REAL(mpk), intent(inout) :: BETA(MAXA), DA(MAXA), ATRY(MAXA), OCHISQ, SIMRADS(MAXDATA), SIMJAC(MAXA,MAXDATA), &
                               SO2TAU(MAXDATA)
   REAL     , intent(inout) :: InvMod_Time

!  Input HPTORA Global Type structure
!mick mod 1/24/2023 - converted "HPTORA_Global" from a scalar to vector of size 2

   TYPE(Create_HPTORA)          , intent(inout) :: HPTORA_Global(2)

!  Mie Global input structure, Mie file/flag
!mick mod 12/9/2022 - To accomodate a Bimodal PSD, made "Mie_Inputs_Global" an array
!                     and added "Bimodal_fraction" and "trace_3"

   LOGICAL                      , intent(in) :: do_Initial_Mie
   CHARACTER*(*)                , intent(in) :: Initial_Mie_File
   TYPE(RTSMie_Inputs)          , intent(in) :: Mie_Inputs_Global(2)
   REAL(dp)                     , intent(in) :: Bimodal_fraction

!  VLIDORT Global input structures

   TYPE(VLIDORT_Fixed_Inputs)   , intent(in)       :: VFixIn_Global
   TYPE(VLIDORT_Modified_Inputs), intent(in)       :: VModIn_Global
   TYPE(VLIDORT_Sup_InOut)      , intent(inout)    :: VSup_Global

!  VLIDORT Global linearized input structures

   TYPE(VLIDORT_Fixed_LinInputs)   , intent(in)     :: VLinFixIn_Global
   TYPE(VLIDORT_Modified_LinInputs), intent(in)     :: VLinModIn_Global
   TYPE(VLIDORT_LinSup_InOut)      , intent(inout)  :: VLinSup_Global

!  Exception handling (initialized in MRQCOF)

   LOGICAL      , intent(inout) :: fail
   CHARACTER*(*), intent(inout) :: messages(3), trace_3

!  Local variables
!  ---------------
!mick mod 9/22/2022 - added vars for ASYMTX & checking COVAR (i.e. the approx. to the Hessian)

   INTEGER   :: J, K
   REAL(mpk) :: DAG(MAXA,1)

!  For ASYMTX

   INTEGER   :: IER
   LOGICAL   :: ASYMTX_FAILURE
   DOUBLE PRECISION :: ASYMTX_TOLERANCE
   DOUBLE PRECISION :: EIGENMAT(MAXA,MAXA), EVEC(MAXA,MAXA), EVAL(MAXA), WK(4*MAXA)

!  For checking COVAR

   INTEGER   :: I, L, LIMIT, NLI
   REAL(mpk) :: ALAMDA_ORIG, STEPSIZE, MIN_STEPSIZE
   LOGICAL   :: BADVAL, CHECK_COVAR

!  Misc

   !LOGICAL   :: VERBOSE = .FALSE.
   LOGICAL   :: VERBOSE = .TRUE.

! Initial Guess
!mick mod 9/12/2022  - added inputs DO_TIMING, NTHREADS, & NUMIT to both calls to sub TONGA_MRQCOF below
!mick mod 12/9/2022  - To accomodate a Bimodal PSD, added "Bimodal_fraction" and "trace_3"
!mick mod 1/10/2023  - added output SIMJAC
!mick mod 10/28/2024 - added output SO2TAU

   IF (ALAMDA .LT. ZERO) THEN
      mfit = 0
      DO J=1,MA
         IF(LISTA(J).ne.0)mfit = mfit + 1
      ENDDO
!      ALAMDA=0.001_mpk
      ALAMDA=0.1_mpk
      CALL TONGA_MRQCOF ( DO_TIMING, NTHREADS,                                                & ! Forward model OMP inputs
                          NUMIT, MAXDATA, MAXA, NDATA, MA, WAVS, RADS, SIGS, do_ss, do_Dump,  & ! Forward model inputs
                          HPTORA_Global, do_Initial_Mie, Initial_Mie_File, Mie_Inputs_Global, & ! Forward model HPTORA + Mie inputs
                          Bimodal_fraction, VFixIn_Global, VModIn_Global, VSup_Global,        & ! Forward model Mie + VLIDORT inputs
                          VLinFixIn_Global, VLinModIn_Global, VLinSup_Global,                 & ! Forward model VLIDORT inputs 
                          APARS, LISTA, MFIT, ALPHA, BETA, CHISQ, SIMRADS, SIMJAC, SO2TAU,    & ! L-M I/O
                          InvMod_Time, fail, messages, trace_3 )                                ! L-M I/O

      if ( fail ) return
      do_Dump = .false.
      write(*,'(1x,a)')'Finished initial simulation call with ALAMBDA set now to 0.1 ----'

      OCHISQ=CHISQ
      DO J=1,MA
         ATRY(J) = APARS(J)
      ENDDO
   ENDIF

!  Loop to define the covariance matrix COVAR (i.e. the approx. to the Hessian) and, if desired,
!  (1) check it to ensure it's positive definite
!  (2) modifiy it if it isn't
!  mick mod 9/22/2022 - defined L-M par mod loop below and placed original COVAR/DA loops inside it

   ALAMDA_ORIG = ALAMDA

   LIMIT = 7
   DO L=1,LIMIT+1

      IF (VERBOSE) write(*,'(/1x,a,i1)') 'In TONGA_MRQMIN: doing covariance loop for L = ',L

!  Define covariance COVAR

      IF (L .EQ. LIMIT+1) THEN
         write(*,'(/a)') 'Unexpected condition in L-M performance COVAR code segment.  Proceeding with original code segment...'
         ALAMDA = ALAMDA_ORIG
      ENDIF

      IF (VERBOSE) write(*,'(1x,a,es12.5)') 'In TONGA_MRQMIN: current value of L-M parameter: ',ALAMDA
      DO J=1,MFIT
         DO K=1,MFIT
            COVAR(J,K)=ALPHA(J,K)
         ENDDO
         COVAR(J,J)=ALPHA(J,J)*(ONE+ALAMDA)
         DA(J)=BETA(J)
      ENDDO

      !Exit COVAR loop (performance code exception exit)
      IF (L .EQ. LIMIT+1) EXIT

!  If desired, check the eigenvalues of COVAR to ensure it's positive definite.
!  If it's not, adjust the L-M par and define COVAR again.

      CHECK_COVAR = .TRUE.
      !CHECK_COVAR = .FALSE.

      IF (CHECK_COVAR) THEN

!  Find the eigenvalues of COVAR
!  (use ASYMTX for now - replace with a pure symmetrix solver later if needed)

         ASYMTX_TOLERANCE = 1.0d-12
         EIGENMAT(1:MFIT,1:MFIT) = COVAR(1:MFIT,1:MFIT)
         CALL ASYMTX ( ASYMTX_TOLERANCE, EIGENMAT, MFIT, MAXA, MAXA, &
                       EVEC, EVAL, IER, WK, MESSAGES(2), FAIL )
         IF ( FAIL  ) THEN
            MESSAGES(1) = 'ASYMTX error in TONGA_MRQMIN'
            RETURN
         ENDIF

!  Check for existence of any negative eigenvalues and possibly the most negative one

         NLI = 0  !index of most neg eigenvalue
         DO I=1,MFIT
            IF (VERBOSE) write(*,'(1x,a,i2,a,es14.7)') 'In TONGA_MRQMIN: I = ',I,' EVAL(I) = ',EVAL(I)
            IF ( EVAL(I) < ZERO ) THEN
               IF ( NLI .EQ. 0 ) THEN
                 !First negative one
                 NLI = I
               ELSE
                 !More than one negative eigenvalue: define the most negative one
                 IF ( EVAL(I) < EVAL(NLI) ) NLI = I
               ENDIF
            ENDIF
         ENDDO
               
!  Modify L-M par if necessary
      
         IF ( NLI > 0 ) THEN
            !Found neg. eigenvalue: increase L-M par now, redefine COVAR, and recheck
            IF (VERBOSE) write(*,'(1x,a)') 'Found neg. eigenvalue: increasing L-M par and redefining COVAR'
            ALAMDA=10.0_mpk*ALAMDA
         ELSE
            !Exit COVAR loop (performance code normal exit)
            EXIT
         ENDIF

      ELSE
         !Exit COVAR loop (original code normal exit)
         EXIT
      ENDIF

   ENDDO

!  Linear algebra
!mick fix 8/15/2022 - limited 1st dim of DA & DAG

   !DAG(:,1) = DA(:)
   DAG(1:MFIT,1) = DA(1:MFIT)
   CALL GAUSSJ ( MAXA, 1, MA, 1, COVAR, DAG, fail, messages(1))
   if ( fail ) return
   !DA(:) = DAG(:,1)
   DA(1:MFIT) = DAG(1:MFIT,1)

   IF ( ALAMDA .EQ. ZERO )THEN
      CALL COVSRT ( MAXA, MA, MFIT, LISTA, COVAR )
      RETURN
   ENDIF

!  Next iteration
!mick mod 9/22/2022 - introduce step size and reduce if necessary to keep elements
!                     of vector "ATRY" in physical space or specified limits
!mick mod 8/17/2023 - added checks for additional pars / use LISTA for those checks

   !DO J=1,MFIT
   !   ATRY(LISTA(J))=APARS(LISTA(J))+DA(J)
   !ENDDO

   STEPSIZE = ONE
   MIN_STEPSIZE = 1.0e-8_mpk
   DO
      BADVAL = .FALSE.
      DO J=1,MFIT
         ATRY(LISTA(J)) = APARS(LISTA(J)) + STEPSIZE*DA(J)
         !General check
         !IF (ATRY(LISTA(J)) <= ZERO) BADVAL = .TRUE.

         IF (LISTA(J) .EQ. 1) THEN
           !AOD check: against specified lower and upper limits
           IF ( ATRY(LISTA(J)) < ZERO .OR. ATRY(LISTA(J)) > 1.0e2_mpk ) THEN
             BADVAL = .TRUE.
             write(*,*) 'AOD < 0.0 or AOD > 100.0: AOD = ',ATRY(LISTA(J))
           ENDIF
         ENDIF

         IF (LISTA(J) .EQ. 2) THEN
           !PeakHgt check: against specified lower and upper limits
           IF ( ATRY(LISTA(J)) < HPTORA_Global(2)%AerLowerLimit .OR. &
                ATRY(LISTA(J)) > HPTORA_Global(2)%AerUpperLimit ) THEN
             BADVAL = .TRUE.
             write(*,*) 'AeroPeakHgt < AerLowerLimit or AeroPeakHgt > AerUpperLimit: AeroPeakHgt = ',ATRY(LISTA(J))
           ENDIF
         ENDIF

         IF (LISTA(J) .EQ. 3) THEN
           !SO2 check: against specified lower and upper limits (in DU)
           IF (ATRY(LISTA(J)) < ZERO .OR. ATRY(LISTA(J)) > 1200.0_mpk) THEN
             BADVAL = .TRUE.
             write(*,*) 'SO2 < 0.0 or SO2 > 1200.0: SO2 = ',ATRY(LISTA(J))
           ENDIF
         ENDIF
      ENDDO

      IF (BADVAL) THEN
         IF (VERBOSE) write(*,'(1x,a,es14.7)') 'In TONGA_MRQMIN: Reducing step size, STEPSIZE = ',STEPSIZE
         STEPSIZE = HALF*STEPSIZE
         IF (STEPSIZE < MIN_STEPSIZE) THEN
            write(*,*) 'Problem: step size in subroutine TONGA_MRQMIN ~ 0.0.  Aborting ....'
            FAIL = .TRUE.
            RETURN
         ENDIF
      ELSE
         EXIT
      ENDIF
   ENDDO

!mick mod 12/9/2022  - To accomodate a Bimodal PSD, added "Bimodal_fraction" and "trace_3"
!mick mod 1/10/2023  - added output SIMJAC
!mick mod 10/28/2024 - added output SO2TAU

   CALL TONGA_MRQCOF ( DO_TIMING, NTHREADS,                                                & ! Forward model OMP inputs
                       NUMIT, MAXDATA, MAXA, NDATA, MA, WAVS, RADS, SIGS, do_ss, do_Dump,  & ! Forward model inputs
                       HPTORA_Global, do_Initial_Mie, Initial_Mie_File, Mie_Inputs_Global, & ! Forward model HPTORA + Mie inputs
                       Bimodal_fraction, VFixIn_Global, VModIn_Global, VSup_Global,        & ! Forward model Mie + VLIDORT inputs
                       VLinFixIn_Global, VLinModIn_Global, VLinSup_Global,                 & ! Forward model VLIDORT inputs 
                       ATRY, LISTA, MFIT, COVAR, DA, CHISQ, SIMRADS, SIMJAC, SO2TAU,       & ! L-M I/O
                       InvMod_Time, fail, messages, trace_3 )                                ! L-M I/O

   if ( fail ) return

!  Test

!mick debug
!write(*,*)
!write(*,*) 'inside TONGA_MRQMIN: before IF condition'
!write(*,*) 'chisq = ',chisq,' ochisq = ',ochisq

   IF ( CHISQ.LT.OCHISQ ) THEN
      ALAMDA = 0.10_mpk*ALAMDA
      OCHISQ = CHISQ
      DO J=1,MFIT
         DO K=1,MFIT
            ALPHA(J,K)=COVAR(J,K)
         ENDDO
         BETA(J)=DA(J)
         APARS(LISTA(J))=ATRY(LISTA(J))
      ENDDO
   ELSE
      ALAMDA=10.0_mpk*ALAMDA
      CHISQ=OCHISQ
   ENDIF


!mick debug
!write(*,*)
!write(*,*) 'inside TONGA_MRQMIN: after IF condition'
!write(*,*) 'chisq = ',chisq,' ochisq = ',ochisq

   RETURN
END SUBROUTINE TONGA_MRQMIN

!

SUBROUTINE TONGA_MRQCOF ( DO_TIMING, NTHREADS,                                                & ! Forward model OMP inputs
                          NUMIT, MAXDATA, MAXA, NDATA, MA, WAVS, RADS, SIGS, do_ss, do_Dump,  & ! Forward model inputs
                          HPTORA_Global, do_Initial_Mie, Initial_Mie_File, Mie_Inputs_Global, & ! Forward model HPTORA + Mie inputs
                          Bimodal_fraction, VFixIn_Global, VModIn_Global, VSup_Global,        & ! Forward model Mie + VLIDORT inputs
                          VLinFixIn_Global, VLinModIn_Global, VLinSup_Global,                 & ! Forward model VLIDORT inputs 
                          A, LISTA, MFIT, ALPHA, BETA, CHISQ, SIMRADS, SIMJAC, SO2TAU,        & ! L-M I/O
                          InvMod_Time, fail, messages, trace_3 )                                ! L-M I/O

!  Double precision

!mick mod 9/12/2022  - added inputs DO_TIMING, NTHREADS, NUMIT
!mick mod 12/9/2022  - To accomodate a Bimodal PSD, made "Mie_Inputs_Global" an array
!                      and added "Bimodal_fraction" and "trace_3"
!mick mod 1/10/2023  - added output SIMJAC
!mick mod 10/28/2024 - added output SO2TAU

   implicit none
   !INTEGER, parameter :: mpk = SELECTED_REAL_KIND(15)

!  OpenMP Control & Timing

   LOGICAL, intent(in)      :: DO_TIMING
   INTEGER, intent(in)      :: NTHREADS

!  Arguments

   LOGICAL  , intent(in)    :: do_ss, do_Dump
   INTEGER  , intent(in)    :: NUMIT, MAXDATA, MAXA, NDATA, MA, MFIT, LISTA(MAXA)
   REAL(mpk), intent(in)    :: WAVS(MAXDATA), RADS(MAXDATA), SIGS(MAXDATA), A(MAXA) 
   REAL(mpk), intent(inout) :: ALPHA(MAXA,MAXA), BETA(MAXA), CHISQ, SIMRADS(MAXDATA), SIMJAC(MAXA,MAXDATA), &
                               SO2TAU(MAXDATA)
   REAL     , intent(inout) :: InvMod_Time

!  Input HPTORA Global Type structure
!mick mod 1/24/2023 - converted "HPTORA_Global" from a scalar to vector of size 2

   Type(Create_HPTORA)          , intent(inout) :: HPTORA_Global(2)

!  Mie Global input structure, Mie file/flag

   LOGICAL                      , intent(in) :: do_Initial_Mie
   CHARACTER*(*)                , intent(in) :: Initial_Mie_File
   TYPE(RTSMie_Inputs)          , intent(in) :: Mie_Inputs_Global(2)
   REAL(dp)                     , intent(in) :: Bimodal_fraction

!  VLIDORT Global input structures

   TYPE(VLIDORT_Fixed_Inputs)   , intent(in)       :: VFixIn_Global
   TYPE(VLIDORT_Modified_Inputs), intent(in)       :: VModIn_Global
   TYPE(VLIDORT_Sup_InOut)      , intent(inout)    :: VSup_Global

!  VLIDORT Global linearized input structures

   TYPE(VLIDORT_Fixed_LinInputs)   , intent(in)    :: VLinFixIn_Global
   TYPE(VLIDORT_Modified_LinInputs), intent(in)    :: VLinModIn_Global
   TYPE(VLIDORT_LinSup_InOut)      , intent(inout) :: VLinSup_Global

!  Exception handling

   LOGICAL      , intent(out) :: fail
   CHARACTER*(*), intent(out) :: messages(3), trace_3

!  Local
!  -----

!  HPTORA type structure
!mick mod 9/22/2022 - replaced by HPTORA_Global

   !Type(Create_HPTORA)       :: HPTORA

!  Mie Input and Output structures
!mick mod 9/22/2022 - added "Mie_Outputs_Global" & "Mie_Outputs_Temp"
!mick mod 12/9/2022 - To accomodate a Bimodal PSD, made "Mie_Inputs" & "Mie_Lin_Outputs" an array

   Type(RTSMie_Inputs)       :: Mie_Inputs(2)
   Type(RTSMie_Outputs)      :: Mie_Outputs, Mie_Outputs_Global, Mie_Outputs_Temp(NDATA)
   Type(RTSMie_Lin_Outputs)  :: Mie_Lin_Outputs(2)

!  VLIDORT input structures

   TYPE(VLIDORT_Fixed_Inputs)    :: VFixIn
   TYPE(VLIDORT_Modified_Inputs) :: VModIn
   TYPE(VLIDORT_Sup_InOut)       :: VSup

!  VLIDORT linearized input structures

   TYPE(VLIDORT_Fixed_LinInputs)    :: VLinFixIn
   TYPE(VLIDORT_Modified_LinInputs) :: VLinModIn
   TYPE(VLIDORT_LinSup_InOut)       :: VLinSup

!  Help variables
!mick mod 9/12/2022 - removed old variable SIMRAD (code uses SIMRADS now)
!                   - changed SIMJAC from 1 to 2 dimensions
!mick mod 2/13/2023 - added RadChk for additional retrieval testing

   LOGICAL   :: Mie_Debug
   INTEGER   :: I, J, K, L, LDUM, N, nmoms, istatus
   REAL(mpk) :: WT, SIG2I, DY
   REAL      :: e1, e2, e3, FM_Time

   REAL(mpk) :: RadChk(2,MAXDATA)

!  Debug

   INTEGER          :: numloadings
   REAL(mpk)        :: loading, maxloading
   CHARACTER(LEN=2) :: char_numit

   INTEGER          :: I_in
   LOGICAL          :: Do_RT
   LOGICAL, SAVE    :: First_Call = .true.

!  OMP Variables
!  -------------
!mick mod 9/12/2022 - added this section

!  OpenMP (general)

   INTEGER   :: OMP_MAXTHREADS
   INTEGER   :: TID, OMP_NTHREADS

!  OpenMP functions

   INTEGER   :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

!  Initialize exception handling

   fail     = .false.
   messages = ' '

!  Initialize timing
!mick mod 9/12/2022 - added additional timing

   if ( DO_TIMING ) Call CPU_time(e1)

!  Dump Calls

   if ( do_Dump ) then
      Open(1,file='RayleighOnly_Dump.dat', status = 'replace' )
      Open(2,file='WithAerosols_Dump.dat', status = 'replace' )
   endif

!  Initialize

   DO J=1,MFIT
      DO K=1,J
         ALPHA(J,K)=ZERO
      ENDDO
      BETA(J)=ZERO
   ENDDO
   CHISQ=ZERO

!  Get the Loading (given A, set aerosol profile)
!mick mod 9/22/2022 - modified message passing to allow all 3 loading fail
!                     messages to get passed
!mick mod 1/24/2023 - apply to HPTORA_Global element 2 only (plume scene)
!mick mod 8/17/2023 - added LISTA

   CALL Tonga_Aerosol_Loading ( MAXA, MA, A, LISTA, HPTORA_Global(2), fail,  messages )
   if ( fail ) return

!  Initial debug aerosol
!   do i = 1, HPTORA_Global(2)%nlayers
!      write(76,*) HPTORA_Global(2)%aerlayerflags(i), HPTORA_Global(2)%heights(i), &
!                  HPTORA_Global(2)%Loading(i), HPTORA_Global(2)%dLoading_dA(i,1:3)
!   enddo
!write(*,*)A,MA,mFIT,MAXA,fail
!stop 'initial ok'

!  8/31/22 Get the SO2 Loading
!mick mod 1/24/2023 - apply to HPTORA_Global element 2 only (plume scene)
!mick mod 8/17/2023 - added LISTA

   if ( HPTORA_Global(2)%Include_SO2 ) then
      CALL Tonga_SO2_Loading ( MAXA, MA, A, LISTA, HPTORA_Global(2), fail,  messages )
      if ( fail ) return
   endif

!  Initial debug SO2

!   if ( HPTORA_Global(2)%Include_SO2 ) then
!      do i = 1, HPTORA_Global(2)%nlayers
!         write(77,*) HPTORA_Global(2)%SO2layerflags(i), HPTORA_Global(2)%heights(i), &
!                     HPTORA_Global(2)%SO2Du(i), HPTORA_Global(2)%dSO2Du_dA(i,1:3)
!      enddo
!      stop 'initial ok'
!   endif

!  Precomputed Mie properties scenario: for the reference wavelength,
!  get the extinction coefficient
!mick mod 9/12/2022 - For OpenMP, also moved reading of Mie data from wvl loop below to here

   if ( Do_Initial_Mie ) then

!  Open Mie file & read preamble (only get aerosol extinction at reference wavelength)
!mick mod 8/15/2022 - to protect Mie file, made it "read only"
!mick mod 1/24/2023 - apply to HPTORA_Global element 2 only

      open(45,file=Trim(Initial_Mie_File), status='old', action='read')
      read(45,*) ; read(45,*) ; read(45,*)
      read(45,'(f12.5,1pe24.12)')HPTORA_Global(2)%Aer_Refw,HPTORA_Global(2)%Aer_Refw_Extinction
      read(45,*) ; read(45,*)
      !write(*,*)' Aer_Refw_Extinction from File ', HPTORA_Global(2)%Aer_Refw, HPTORA_Global(2)%Aer_Refw_Extinction

!  Read pre-computed Mie properties FOR ALL WVLS here (array "Mie_Outputs_Temp" is fully filled)

      do I=1,NDATA
         read(45,*) ; read(45,*) ; read(45,*)
         read(45,'(1p4e24.12)') Mie_Outputs_Temp(I)%Mie_Bulk(1:3), Mie_Outputs_Temp(I)%Mie_Asymm
         read(45,'(I5)') Mie_Outputs_Temp(I)%Mie_NCoeffs ; nmoms = Mie_Outputs_Temp(I)%Mie_NCoeffs
         do L = 0, nmoms
            read(45,'(I5,1p6e24.12)') LDUM, Mie_Outputs_Temp(I)%Mie_Expcoeffs(1:6,L)
         enddo
         Mie_Outputs_Temp(I)%Mie_Expcoeffs(1:6,nmoms+1:max_Mie_angles) = ZERO
         !if ( I.eq.1 ) then
             Mie_Outputs_Temp(I)%Mie_Fmatrix  = ZERO
             Mie_Outputs_Temp(I)%Mie_dist1    = ZERO
             Mie_Outputs_Temp(I)%Mie_dist2    = ZERO
             Mie_Outputs_Temp(I)%Mie_distplot = ZERO
         !endif
      enddo

!  Close file

      close(45)

   endif

!  Freshly-computed Mie properties scenario: for REFERENCE wavelength,
!  ONLY get the extinction coefficient and possibly its derivative wrt n_real or n_imag.

   if ( .not. Do_Initial_Mie ) then

!  Copy global input, set wavelength. No expansion coefficients for this call.
!mick note 9/22/2022 - Lambda in microns

      Mie_Inputs        = Mie_Inputs_Global
      Mie_Inputs%Lambda = HPTORA_Global(2)%Aer_Refw * 0.001_mpk
      Mie_Inputs%do_Expcoeffs = .false.

!  Adjust n_real or n_imag, set derivatives flag.
!    Call to Mie code - either with or without linearization.
!mick mod 12/9/2022 - replace std and lin Mie calls with Bimodal ones
!mick mod 8/17/2023 - replaced "MA .eq. 3" IF condition with the Retrieve_nReal & Retrieve_nImag IF conditions
!                   - replaced A(3) with new refrac index location A(LISTA(5))

      if ( HPTORA_Global(2)%Retrieve_nReal .or. HPTORA_Global(2)%Retrieve_nImag ) then
         Mie_Inputs%do_LinearREF = .true.
         if ( HPTORA_Global(2)%Retrieve_nReal ) Mie_Inputs%n_real = A(LISTA(5))
         if ( HPTORA_Global(2)%Retrieve_nImag ) Mie_Inputs%n_imag = A(LISTA(5))
         !call RTSMie_Master_plus & 
         !  ( Mie_Inputs, Mie_Outputs, Mie_Lin_Outputs,            & ! I/O type structures
         !    fail, istatus, messages(1), messages(2), messages(3) ) ! Exception handling
         call RTSMie_master_bimodal_plus &
           ( Mie_Inputs(1), Mie_Inputs(2), Bimodal_fraction,      & ! Inputs
             Mie_Outputs, Mie_Lin_Outputs(1), Mie_Lin_Outputs(2), & ! Outputs
             fail, istatus, messages, trace_3 )
      else
         !call RTSMie_Master & 
         !  ( Mie_Inputs, Mie_Outputs,                             & ! I/O type structures
         !    fail, istatus, messages(1), messages(2), messages(3) ) ! Exception handling
         call RTSMie_master_bimodal &
            ( Mie_Inputs(1), Mie_Inputs(2), Bimodal_fraction, & ! Inputs
              Mie_Outputs, fail, istatus, messages, trace_3 )   ! Outputs
      endif
      if ( fail ) return

!  Mie debug
!mick mod 8/17/2023 - replaced "MA .eq. 3" IF condition with the Retrieve_nReal & Retrieve_nImag IF conditions

!      Mie_Debug = .true.
      Mie_Debug = .false.
      if ( Mie_Debug ) then
         if ( HPTORA_Global(2)%Retrieve_nReal .or. HPTORA_Global(2)%Retrieve_nImag ) then
            Call RTSMie_write_extended ( 'Mie_Debug_0.out', .TRUE., Mie_Inputs(1), Mie_Outputs, Mie_Lin_Outputs(1) )
         else
            Call RTSMie_write_standard ( 'Mie_Debug_0.out', Mie_Inputs(1), Mie_Outputs )
         endif
         stop 'Mie Debug'
      endif

!  Set the HPTORA Global aerosol extinction variables from Mie output. Derivatives should be single normalized.
!mick note 12/9/2022 - Mie_Bulk(1): extinction cross-section in um^2
!mick mod 8/17/2023 - replaced "MA .eq. 3" IF condition with the Retrieve_nReal & Retrieve_nImag IF conditions

      HPTORA_Global(2)%Aer_Refw_Extinction = Mie_Outputs%Mie_Bulk(1)
      if ( HPTORA_Global(2)%Retrieve_nReal .or. HPTORA_Global(2)%Retrieve_nImag ) then
         if ( HPTORA_Global(2)%Retrieve_NReal ) HPTORA_Global(2)%Aer_Refw_L_Extinction = A(3) &
                                                * Mie_Lin_Outputs(1)%LRFE_Mie_Bulk(1,1)
         if ( HPTORA_Global(2)%Retrieve_NImag ) HPTORA_Global(2)%Aer_Refw_L_Extinction = A(3) &
                                                * Mie_Lin_Outputs(1)%LRFE_Mie_Bulk(1,2)
      else
         HPTORA_Global(2)%Aer_Refw_L_Extinction = ZERO
      endif

!  This linearization checks out, 27 June 2022
!  write(*,*)Mie_Inputs%n_real, HPTORA_Global(2)%Aer_Refw_Extinction, HPTORA_Global(2)%Aer_Refw_L_Extinction
!   stop

!      write(*,*)' Aer_Refw_Extinction, from Scratch', HPTORA_Global(2)%Aer_Refw, HPTORA_Global(2)%Aer_Refw_Extinction

   endif

!  Check HPTORA inputs before passing to Forward model
!mick mod 9/22/2022 - added this call (on standby)
!mick mod 1/24/2023 - added loop

   !do i=1,2
   !   write(*,'(/1x,a,i1,a)') 'For HPTORA_Global(',i,')'
   !   CALL TONGA_Write_HPTORA_V5 ( 0, MA, HPTORA_Global(i) )
   !enddo

!  Get inverse model prep timing
!mick mod 9/12/2022 - added this section

   if ( DO_TIMING ) then
      Call CPU_time(e2)
      FM_Time = 0.0
      InvMod_Time = InvMod_Time + e2-e1
   endif

!  Main forward model calculation, over wavelength loop. 
!  -----------------------------------------------------
!mick mod 9/12/2022 - added OpenMP parallel region & OMP call below

   !Do_RT = .true.
   !Do_RT = .false.

   !if (First_Call) then
     !open (unit=10,file='radjac_test_file.dat')
     !open (unit=10,file='radjac_test_file_OMP_nt1.dat_2nd')
     !open (unit=10,file='radjac_test_file_OMP_nt2.dat_2nd')
     !open (unit=10,file='radjac_test_file_serial.dat_2nd')
     !open (unit=10,file='radjac_test_file_OMP_nt2.dat_3rd')
     !open (unit=10,file='radjac_test_file_dud.dat')
     !First_Call = .false.
   !endif

   !if (Do_RT) then

!  Set total number of threads (from NTHREADS input)
!     write(*,*) 'Setting total number of threads'

      OMP_MAXTHREADS = NTHREADS
      CALL OMP_SET_NUM_THREADS(OMP_MAXTHREADS)

!  Begin parallel region
!mick mod 12/9/2022  - To accomodate a Bimodal PSD, added "Bimodal_fraction"
!mick mod 8/17/2023  - added LISTA
!mick mod 10/28/2024 - added SO2TAU

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, &
!$OMP     NUMIT, MAXA, MA, A, LISTA, NDATA, WAVS,              & ! Ret Info
!$OMP     do_ss, do_Dump, HPTORA_Global,                       & ! Fwd model ctrl, HPTORA inputs
!$OMP     Do_Initial_Mie, Mie_Inputs_Global, Bimodal_fraction, & ! Mie I/O
!$OMP     Mie_Outputs_Temp,                                    & ! Mie I/O
!$OMP     VFixIn_Global, VModIn_Global, VSup_Global,           & ! VLIDORT inputs 
!$OMP     VLinFixIn_Global, VLinModIn_Global, VLinSup_Global,  & ! VLIDORT inputs
!$OMP     SIMRADS, SIMJAC, SO2TAU, RadChk)                       ! VLIDORT outputs

!  Obtain thread number

      tid = OMP_GET_THREAD_NUM()
      !tid = 0

!  Obtain and display total number of threads and local thread

      if (tid == 0) then
         omp_nthreads = OMP_GET_NUM_THREADS()
         !omp_nthreads = 1
         !write(*,*)
         !write(*,'(1x,a,i1)') 'total number of threads (inside parallel region) = ', omp_nthreads
         write(*,'(1x,a,i1,a)') 'Running forward model with ',omp_nthreads,' thread(s)'
         write(*,'(1x,a,i1/)')  ' *** nstokes = ',VFixIn_Global%Cont%TS_NSTOKES
      end if

      write(*,'(1x,a,i1,a)') 'Thread ',tid,' beginning its wavenumber loop'

!$OMP DO

!  Begin "wavenumber" loop

      DO I=1,NDATA

         !if (i == 1) write(*,*)
         !write(*,'(1x,a,i1,a,i5)') 'FM loop: wvl = ',i,' Thread  = ',tid

!  Get Mie data from Mie temporary output array

         if ( Do_Initial_Mie ) then

!  Pass pre-computed Mie properties for CURRENT WVL ONLY (array "Mie_Outputs_Global" is fully filled)

            Mie_Outputs_Global%Mie_Bulk(1:3) = Mie_Outputs_Temp(I)%Mie_Bulk(1:3)
            Mie_Outputs_Global%Mie_Asymm     = Mie_Outputs_Temp(I)%Mie_Asymm
            Mie_Outputs_Global%Mie_NCoeffs   = Mie_Outputs_Temp(I)%Mie_NCoeffs
            Mie_Outputs_Global%Mie_Expcoeffs(1:6,0:max_Mie_angles) = Mie_Outputs_Temp(I)%Mie_Expcoeffs(1:6,0:max_Mie_angles)
            !if ( I.eq.1 ) then
                Mie_Outputs_Global%Mie_Fmatrix  = Mie_Outputs_Temp(I)%Mie_Fmatrix
                Mie_Outputs_Global%Mie_dist1    = Mie_Outputs_Temp(I)%Mie_dist1
                Mie_Outputs_Global%Mie_dist2    = Mie_Outputs_Temp(I)%Mie_dist2
                Mie_Outputs_Global%Mie_distplot = Mie_Outputs_Temp(I)%Mie_distplot
            !endif
         endif

!  Call forward model results, generates radiance and Jacobians. 7/28/22. Version 2.
!mick mod 9/12/2022 - added inputs TID & NUMIT
!                   - modified error handling for OpenMP case
!mick mod 12/9/2022 - To accomodate a Bimodal PSD, added "Bimodal_fraction"
!mick mod 8/17/2023 - added LISTA

         CALL TONGA_Forward_Model_V4 &
            ( TID, NUMIT, MAXA, MA, A, LISTA, I, WAVS(I),         & ! Ret Info
              do_ss, do_Dump, HPTORA_Global,                      & ! Fwd Model Control, HPTORA inputs
              Do_Initial_Mie, Mie_Inputs_Global,                  & ! Mie I/O
              Bimodal_fraction, Mie_Outputs_Global,               & ! Mie I/O
              VFixIn_Global, VModIn_Global, VSup_Global,          & ! VLIDORT inputs 
              VLinFixIn_Global, VLinModIn_Global, VLinSup_Global, & ! VLIDORT inputs 
              SIMRADS(I), SIMJAC(:,I), SO2TAU(I), RadChk(:,I), fail, messages, trace_3 ) ! Outputs and exception handling

         if ( fail ) then
            write(*,'(/a)') ' Error: problem in TONGA_Forward_Model_V2a.  The error message(s) are:'
            do n = 1,3
              write(*,'(8x,a,i1,a)') 'Message # ',n,' : ' // trim(adjustl(messages(n)))
            enddo
            stop
         endif

         !write(*,*) 'Done point', I, WAVS(I), RADS(I), SIMRADS(I), SIGS(I), SIMJAC(1:MA,I)

!  End wavelength loop

      ENDDO

!  End OMP loop

!$OMP END DO

!  End parallel region

!$OMP END PARALLEL

   !DO I=1,NDATA
   !  write(10,*) I, SIMRADS(I), SIMJAC(1:MA,I)
   !ENDDO

!  Calculate retrieval Alpha, Beta and CHISQ
!mick mod 9/12/2022 - moved this calculation from wvl loop above to here for OpenMP

   DO I=1,NDATA
      SIG2I = ONE/(SIGS(I)*SIGS(I))
      DY = RADS(I) - SIMRADS(I)
      DO J=1,MFIT
         WT = SIMJAC(LISTA(J),I)*SIG2I
         DO K=1,J
            ALPHA(J,K) = ALPHA(J,K) + WT*SIMJAC(LISTA(K),I)
         ENDDO
         BETA(J) = BETA(J) + DY*WT
      ENDDO
      CHISQ = CHISQ + DY*DY*SIG2I
   ENDDO

!  Inverse model timing
!mick mod 9/12/2022 - added additional timing

   if ( DO_TIMING ) then
      Call CPU_time(e3)
      FM_Time = (e3-e2)/real(NTHREADS)
      write(*,'(/1x,a,es12.5)')'FM: One cycle, timing = ',FM_Time
      InvMod_Time = InvMod_Time + FM_Time
   endif

!  Covariance matrix symmetry

   DO J=2,MFIT
      DO K=1,J-1
         ALPHA(K,J)=ALPHA(J,K)
      ENDDO
   ENDDO

!  End Dump

   if ( do_Dump ) then
      close(1) ; close(2)
   endif

!  Done

   RETURN
END SUBROUTINE TONGA_MRQCOF
 
!

SUBROUTINE COVSRT ( MAXA, MA, MFIT, LISTA, COVAR )

!  double precision

   implicit none
   integer, parameter :: mpk = SELECTED_REAL_KIND(15)

!  arguments

   INTEGER  , intent(in)    :: MAXA, MA, MFIT, LISTA(MAXA)
   REAL(mpk), intent(inout) :: COVAR(MAXA,MAXA)

!  Local

   INTEGER   :: I, J
   REAL(mpk) :: SWAP

   DO J=1,MA-1
      DO I=J+1,MA
         COVAR(I,J)=ZERO
      ENDDO
   ENDDO
   DO I=1,MFIT-1
      DO J=I+1,MFIT
         IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
         ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
         ENDIF
      ENDDO
   ENDDO
   SWAP=COVAR(1,1)
   DO J=1,MA
      COVAR(1,J)=COVAR(J,J)
      COVAR(J,J)=ZERO
   ENDDO
   COVAR(LISTA(1),LISTA(1))=SWAP
   DO J=2,MFIT
      COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
   ENDDO
   DO J=2,MA
      DO I=1,J-1
         COVAR(I,J)=COVAR(J,I)
      ENDDO
   ENDDO
   RETURN
END SUBROUTINE COVSRT

!

SUBROUTINE GAUSSJ ( NP, MP, N, M, A, B, fail, message)

!  double precision

   implicit none
   integer, parameter :: mpk = SELECTED_REAL_KIND(15)

!  arguments

   INTEGER  , intent(in)    :: NP, MP, N, M
   REAL(mpk), intent(inout) :: A(NP,NP), B(NP,MP)
   logical      , intent(inout) :: fail
   character*(*), intent(inout) :: message

!  Local

   INTEGER, parameter :: NMAX = 50
   INTEGER   :: I, J, K, IROW, ICOL, L, LL
   INTEGER   :: IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
   REAL(mpk) :: BIG, DUM, PIVINV

   DO J=1,N
      IPIV(J)=0
   ENDDO
   DO I=1,N
      BIG=ZERO
      DO J=1,N
         IF(IPIV(J).NE.1)THEN
            DO K=1,N
               IF (IPIV(K).EQ.0) THEN
                  IF (ABS(A(J,K)).GE.BIG)THEN
                     BIG=ABS(A(J,K))
                     IROW=J
                     ICOL=K
                  ENDIF
               ELSE IF (IPIV(K).GT.1) THEN
                  message = 'GaussJ, IPIv(K) > 1 Singular matrix' ; fail = .true. ; return
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      IPIV(ICOL)=IPIV(ICOL)+1
      IF (IROW.NE.ICOL) THEN
         DO L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
         ENDDO
         DO L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
         ENDDO
      ENDIF
      INDXR(I)=IROW
      INDXC(I)=ICOL
      IF (A(ICOL,ICOL) .EQ. ZERO) then
         message = 'GaussJ, A(ICOL,ICOL) = 0, Singular matrix' ; fail = .true. ; return
      endif
      PIVINV=ONE/A(ICOL,ICOL)
      A(ICOL,ICOL)=ONE
      DO L=1,N
         A(ICOL,L)=A(ICOL,L)*PIVINV
      ENDDO
      DO L=1,M
         B(ICOL,L)=B(ICOL,L)*PIVINV
      ENDDO
      DO LL=1,N
         IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=ZERO
            DO L=1,N
               A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
            ENDDO
            DO L=1,M
               B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
            ENDDO
         ENDIF
      ENDDO
   ENDDO

   DO L=N,1,-1
      IF(INDXR(L).NE.INDXC(L))THEN
         DO K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
         ENDDO
      ENDIF
   ENDDO

   RETURN
END SUBROUTINE GAUSSJ

!

End Module TONGA_Levenberg_marquardt_V4_m

!subroutine Function_Simulator ( MAXA, MA, A, X, YMOD, DYDA )
!   implicit none
!   integer, parameter :: mpk = SELECTED_REAL_KIND(15)
!   INTEGER  , intent(in)  :: MAXA, MA
!   REAL(mpk), intent(in)  :: X, A(MAXA)
!   REAL(mpk), intent(inout) :: YMOD, DYDA(MAXA)
!   return
!end subroutine Function_Simulator


