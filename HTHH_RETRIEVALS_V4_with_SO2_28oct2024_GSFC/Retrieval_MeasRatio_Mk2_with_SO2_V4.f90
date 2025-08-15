program TONGA_Retrieval_MeasRatio_V4

!  5/26/22. Version extending the validity area for the plume from 20-45 km
!  5/31/22. Version with different aerosol properties
!     5/31/22 -- 6/6/22. Groups 2, 3, 4
!  6/9/22. Version with aerosol properties depending on wavelength (Group 5)

!  6/22/22. Version with Mie calculations from scratch and retrieval of refractive indices

!  7/28/22. Retrievals with ratioed measurements.
!  8/31/22. Introduce optional SO2 (total DU from TROPOMI, same plume location as aerosol)

!@@@@@@@@@@@@@@@@@@@@@@@@@@

!  VLIDORT modules

      USE VLIDORT_PARS_m, Only: PI4
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_LIN_IO_DEFS_m

!  TONGA environment modules
!mick mod 12/9/2022 - added TONGA_Mie_Initialize_m & RTSMie_master_bimodal_m

      USE TONGA_pars_m
      USE TONGA_HPTORA_V5_m

      USE TONGA_Mie_Initialize_m
      USE TONGA_RTInitialize_V2_m
      USE TONGA_Create_HPTORA_V5_m
      USE TONGA_Levenberg_marquardt_V4_m

!  Mie modules

      USE RTSMie_parameters_m
      USE RTSMie_Inputs_Def_m
      !USE RTSMie_Master_m
      USE RTSMie_master_bimodal_m

      implicit none

!  Type structures
!  ---------------

!  HPTORA Global scene input structure
!mick mod 1/24/2023 - converted "HPTORA_Global" from a scalar to vector of size 2

      Type(Create_HPTORA) :: HPTORA_Global(2)

!  VLIDORT Global input structures

      TYPE(VLIDORT_Fixed_Inputs)        :: VFixIn_Global
      TYPE(VLIDORT_Modified_Inputs)     :: VModIn_Global
      TYPE(VLIDORT_Sup_InOut)           :: VSup_Global

      TYPE(VLIDORT_Fixed_LinInputs)     :: VLinFixIn_Global
      TYPE(VLIDORT_Modified_LinInputs)  :: VLinModIn_Global
      TYPE(VLIDORT_LinSup_InOut)        :: VLinSup_Global

!  Mie Global input structure
!mick mod 12/9/2022 - To accomodate a Bimodal PSD, made "Mie_Inputs_Global" an array

      Type(RTSMie_Inputs)               :: Mie_Inputs_Global(2)

!  Mie Local Output structure

      Type(RTSMie_Outputs)              :: Mie_Outputs

!  Top-level configuration inputs
!  ------------------------------
!mick mod 12/9/2022 - To accomodate a Bimodel PSD, expanded "firstguess" inputs
!                     and added "Bimodal_fraction"
!mick mod 1/24/2023 - converted "HPTO3_FILENAME" from a scalar to vector of size 2
!mick mod 8/17/2023 - added "retrieve_plumehwhm"

      logical       :: do_ssonly
      integer       :: nstokes, nstreams, nssfine
      character*100 :: SO2XSec_Source
      character*100 :: atmos_datapath, rad_datapath, mie_datapath, physics_datapath
      character*100 :: HPTO3_FILENAME(2), radfile
      real(mpk)     :: lambda_1, lambda_2
      real(mpk)     :: plumetop, plumebot, plumehwhm_firstguess
      logical       :: retrieve_plumehwhm
      real(mpk)     :: AOD_refw, AOD_firstguess, Pkht_firstguess
      real(mpk)     :: nreal_firstguess(2), nimag_firstguess(2)
      real(mpk)     :: PSDpar1_firstguess(2), PSDpar2_firstguess(2)
      real(dp)      :: Bimodal_fraction
      logical       :: retrieve_nreal, retrieve_nimag

!  8/31/22. Add SO2 control
!mick mod 8/17/2023 - changed "SO2_Totalcolumn_DU" to "SO2_firstguess"
!                   - added "retrieve_SO2"

      real(mpk)     :: SO2_firstguess
      logical       :: include_SO2, retrieve_SO2

!  other settings
!mick mod 2/22/2024 - replaced "albedo" with "Surface_height" and "Surface_albedo"

      integer       :: n1, n2, nfine
      real(mpk)     :: eradius
      real(mpk)     :: Surface_height, Surface_albedo

!  Data quantities
!  ---------------

!  2 radiance spectra at the same wavelength grid, no irradiance
!mick mod 10/28/2024 - added "data_snrs"

      integer       :: data_nwavs
      real(mpk)     :: data_wavs   (E_MAXDATAWAVS)
      real(mpk)     :: data_rads   (E_MAXDATAWAVS,2)
      real(mpk)     :: data_snrs   (E_MAXDATAWAVS,2)
!      real(mpk)     :: data_irwavs (E_MAXDATAWAVS)
!      real(mpk)     :: data_irrads (E_MAXDATAWAVS)

!  Related to moving to 2 bands
!mick mod 7/31/2023  - added this 2-band section
!mick mod 10/28/2024 - added "data_snrs_in"

!mick test
      integer, parameter :: maxbands=2

      integer       :: ib, w_off, w_in, w0, w1
      integer       :: numbands,NumRadAvg
      integer       :: data_nwavs_in(maxbands),nbanddata(maxbands)
      integer       :: data_band_flag(E_MAXDATAWAVS),band_flag(E_MAXDATAWAVS) 

      real(mpk)     :: band_signal_to_noise
      real(mpk)     :: lat_in(2,maxbands), lon_in(2,maxbands), &
                       SZA_in(2,maxbands), VZA_in(2,maxbands), &
                       SAZ_in(2,maxbands), VAZ_in(2,maxbands)

      real(mpk)     :: data_wavs_in(E_MAXDATAWAVS,maxbands)
      real(mpk)     :: data_rads_in(E_MAXDATAWAVS,2,maxbands)
      real(mpk)     :: data_snrs_in(E_MAXDATAWAVS,2,maxbands)

      character*2   :: CharRadAvg(maxbands)
      character*3   :: CharRadBand(maxbands)
      character*100 :: which_radfile(2,maxbands)

!  2 of everything here

      real(mpk) :: lat(2), lon(2), SZA(2), VZA(2), SAZ(2), VAZ(2)

!  Buffered data. 1 wavelength grid, 2 radiance spectra, no irradiance
!mick mod 10/28/2024 - added "snrs"

      integer   :: ndata
      real(mpk) :: wavs (E_MAXWAVS)
      real(mpk) :: rads (E_MAXWAVS,2)
      real(mpk) :: snrs (E_MAXWAVS,2)
!      real(mpk) :: irrads (E_MAXWAVS)

!  Measured radiance-ratios, and simulated versions
!mick mod 1/10/2023 - added sigs_orig

      real(mpk) :: measrads (E_MAXWAVS)
      real(mpk) :: simrads (E_MAXWAVS), simrads_saved (50,E_MAXWAVS)
      real(mpk) :: sigs_orig (E_MAXWAVS), sigs (E_MAXWAVS), snr (E_MAXWAVS), snr_ratio (E_MAXWAVS)

!  Inverse Model
!  -------------

!  Parameter data input to L-M model
!mick fix 8/30/2023 - define MAXA by E_MAXPARS
!mick mod 9/1/2023  - changed name from CONV to PAR_CONV

      !real(mpk), parameter  :: CONV = 1.0d-04, CHISQ_CONV = 1.0d-04
      real(mpk), parameter  :: PAR_CONV = 1.0d-02, CHISQ_CONV = 1.0d-02
      integer  , parameter  :: MAXA = E_MAXPARS, MAXITER = 50

!  Forward model Jacobians
!mick mod 1/10/2023 - added simjac

      real(mpk) :: simjac(MAXA,E_MAXWAVS)

!  L-M method I/O and retrieval output
!mick mod 10/28/2024 - added SO2TAU & SO2AMF

      integer   :: LISTA(MAXA), MFIT, MA, IPARS, NPARS, CONVCOUNT
      real(mpk) :: ALPHA(MAXA,MAXA), COVAR(MAXA,MAXA), CHISQ, ALAMDA, APARS(MAXA), APARS_OLD(MAXA), &
                   Mom, Cutoff
      real(mpk) :: BETA(MAXA), DA(MAXA), ATRY(MAXA), OCHISQ, chisq_change, old_chisq, SDPARS(MAXA), &
                   CROSSCOR(MAXA,MAXA), VZERO(MAXA)
      real(mpk) :: SO2TAU(E_MAXWAVS), SO2AMF(E_MAXWAVS)

!  3/31/22.  Add signal to noise

      real(mpk) :: signal_to_noise, error
 
!  Mie flags

      logical       :: Do_Initial_Mie_Calculations
      logical       :: Use_Initial_Mie_Calculations 
      character*156 :: Initial_Mie_File

!  Local
!  -----
!mick mod 12/9/2022 - To accomodate a Bimodel PSD, added "Mie_trace_3"
!mick mod 8/17/2023 - added "CRetPar" & "CRetPar2" / modified "CHead"

      !integer       :: iscan, iscan_start, iscan_end, nscans
      integer       :: test, TestCfg
      integer       :: i, j, w, wc, q, qc, jsav, U33, U88, U99, U100, k, L, nmoms, &
                       Mie_istatus, NoiseOption, Test_option
      real(mpk)     :: avg_snr, avg_snr_ratio
      logical       :: trawl, fail, Region_extension, Mie_Fail

      character*1   :: CharTestCfg
      character*3   :: CWLo, CWHi
      character*4   :: Scene
      character*5   :: CSg, CNr, AeroOrbit, BackGndOrbit
      character*6   :: CBf, CRg, CNi
      character*7   :: SceneExt

      !character*10  :: HPTO3_file_prefix
      character*14  :: HPTO3_file_prefix

      character*13  :: FileScene, &
                       RaySceneAtmosIndex, AeroSceneAtmosIndex, &
                       RaySceneRadIndex,   AeroSceneRadIndex
      character*17  :: CRetPar(MAXA), CRetPar2(MAXA)
      character*20  :: RadFilePrefix
      character*40  :: SceneRadsLabel, FilePrefix, FileSceneExt
      character*120 :: messages(3), ExtraHeader, ExtraHeader_Res, Mie_messages(3), Mie_trace_3
      character*150 :: CHead

!  Retrieval orbit types

      character*2   :: Orbit22086_RetType, Orbit22087_RetType, RetType

!  Retrieval pre-screening
!mick mod 12/2/2022 - added this section

      integer       :: idx_lo, idx_hi, num_idx
      real(mpk)     :: CSI, CSI_threshold
      logical       :: no_aerosol

!  Including forward model error; adjusting chi-square

      integer       :: chi_loop
      real(mpk)     :: rndata, deg_freedom, std_dev, chisq_lower_lim, chisq_upper_lim, &
                       chisq_mu, varF, sigF(E_MAXWAVS), epsF(E_MAXWAVS)

!  Timing tests

      real :: e1, e2

!  Artificial Dump for the VLRRS runs

      logical :: do_Dump
      logical, parameter :: parameter_do_Dump = .false.
!      logical, parameter :: parameter_do_Dump = .true.

!  Start program

!  SETUPS
!  ======

!  Set up timing file
!mick mod 5/9/2023 - added this section

   U100 = 100
   Open(unit=U100,file='Tonga_ret_timing.log',status='replace')

!  Top level configuration read
!  ----------------------------
!  8/31/22. Add SO2 control; rename file (V2a)
!mick mod 12/2/2022 - added pre-screening section / changed config file name to V3
!mick mod 12/9/2022 - to accomodate a Bimodel PSD, expanded "firstguess" inputs
!                     and added "Bimodal_fraction"
!mick mod 8/17/2023 - modified config file & inputs to accomodate retrieving SO2
!mick fix 8/17/2023 - added "retrieve_plumehwhm" to be able to retrieve plume HWHM

      OPEN(1,file='Toplevel_Inputs_V4.cfg', STATUS = 'old')
      read(1,*) ! RETRIEVAL SCENE PRE-SCREENING
      read(1,*)idx_lo
      read(1,*)idx_hi
      read(1,*)CSI_threshold
      read(1,*) !       RT CONTROL
      read(1,*)nstokes               ! NSTOKES in VLIDORT
      read(1,*)nstreams              ! NSTREAMS in VLIDORT
      read(1,*)nssfine               ! fine layer control for FOCORR_OUTGOING in VLIDORT
      read(1,*)do_ssonly             ! No multiple scatter in VLIDORT ==> very fast !
      read(1,*) !       WAVELENGTH LIMITS
      read(1,*)lambda_1              ! Lower wavelength buffer limit [nm]
      read(1,*)lambda_2              ! Upper wavelength buffer limit [nm]
      read(1,*) !       PLUME CONTROL
      read(1,*)plumetop              ! Upper limit   of plume [km]
      read(1,*)plumebot              ! Lower limit   of plume [km]
      read(1,*)plumehwhm_firstguess  ! halfwidth half max of plume [km]
      read(1,*)retrieve_plumehwhm    ! Flag for retrieving plume hwhm
      read(1,*)AOD_refw              ! Reference wavelength [nm]
      read(1,*)AOD_firstguess        ! First-guess value of AOD
      read(1,*)Pkht_firstguess       ! First-guess value of Peak Height [km]
      read(1,*) !       SO2 CONTROL
      read(1,*)Include_SO2           ! Flag for including SO2
      read(1,*)SO2XSec_Source        ! SO2 XSec source file
      read(1,*)SO2_firstguess        ! First-guess value of total column SO2 (DU)
      read(1,*)retrieve_SO2          ! Flag for retrieving SO2
      read(1,*) !       MIE CONTROL
      read(1,*)retrieve_nreal        ! Flag for retrieving nreal
      read(1,*)retrieve_nimag        ! Flag for retrieving nimag
      read(1,*)Bimodal_fraction      ! Bimodal fraction (PSD #1)
      read(1,*)nreal_firstguess(1)   ! Mie program input: PSD #1 (fine mode)
      read(1,*)nimag_firstguess(1)   ! Mie program input: PSD #1 (fine mode)
      read(1,*)psdpar1_firstguess(1) ! Mie program input: PSD #1 (fine mode)
      read(1,*)psdpar2_firstguess(1) ! Mie program input: PSD #1 (fine mode)
      read(1,*)nreal_firstguess(2)   ! Mie program input: PSD #2 (coarse mode)
      read(1,*)nimag_firstguess(2)   ! Mie program input: PSD #2 (coarse mode)
      read(1,*)psdpar1_firstguess(2) ! Mie program input: PSD #2 (coarse mode)
      read(1,*)psdpar2_firstguess(2) ! Mie program input: PSD #2 (coarse mode)
      read(1,*)Do_Initial_Mie_Calculations  ! Performs desired Mie calculations and dumps Mie results to file ONLY
      read(1,*)Use_Initial_Mie_Calculations ! Uses dumped Mie data directly if TRUE; otherwise, Mie is computed on the fly!
      read(1,*) !       SURFACE CONTROL
      read(1,*)Surface_height        !Surface or cloudtop height [km]
      read(1,*)Surface_albedo        !Surface or cloudtop albedo
      close(1)

!  Begin test

!mick test
      !Single retrieval
      test=1 !22086
      !test=2 !22087

      write(*,'(/1x,a,i2)') 'Doing orbit = ',test

!  Begin Chi-Square Check Loop (this loop ends @ ~line 1630)
!mick mod 1/10/2023 - added this loop

      chi_loop = 0
      sigF = zero; epsF = zero

      DO

      chi_loop = chi_loop + 1
      write(*,'(/1x,a,i1)') 'Doing chi_loop = ',chi_loop

!  Start timing

      call cpu_time(e1)

!  Cannot have both Mie calculation flags

      if ( Do_Initial_Mie_Calculations .and. Use_Initial_Mie_Calculations ) then
         write(*,*)
         stop 'Cannot have both Mie calculation inputs TRUE - aborting program!!'
      endif

!  Cannot retrieve nreal or nimag using a precomputed Mie file (Aero ext at the ref wvl is not updated using the most
!  up-to-date nreal or nimag during the retrieval iteration if using a precomputed Mie file)
!mick fix 9/22/2022 - added this section

      if ( (retrieve_nreal .or. retrieve_nimag) .and. Use_Initial_Mie_Calculations ) then
         write(*,*)
         write(*,*) 'Cannot retrieve nreal or nimag using a precomputed Mie input file.'
         write(*,*) 'Must compute Mie inputs during retrieval - aborting program!!'
         stop
      endif

!  Define # of parameters in retrieval
!  -----------------------------------

      npars = 2

      if ( retrieve_SO2 ) npars = npars + 1

      if ( retrieve_plumehwhm .or. retrieve_nreal .or. retrieve_nimag ) &
           stop 'Not totally set up to retrieve Plume HWHM, Nreal, or NImag at this time ==> STOP PROGRAM!!!'

!      if ( retrieve_plumehwhm ) npars = npars + 1
!      if ( retrieve_nreal .or. retrieve_nimag ) npars = npars + 1

!  Define parameters in retrieval / set initial values
!  ---------------------------------------------------
!mick mod 8/17/2023 - modified to retrieve up to five pars
!                   - moved this section from near start of ret iteration to here
!                     to provide earlier information to user

!  The elements of vector "LISTA" tell you WHERE in state vector "APARS"
!  to find the following retrieval parameters during a given active
!  retrieval (if applicable):
!     LISTA(1) - Aerosol column optical depth
!     LISTA(2) - Aerosol Pk Hgt
!     LISTA(3) - SO2 
!     LISTA(4) - PlumeHWHM
!     LISTA(5) - Fine mode PSD NReal or NImag (only one allowed in a given retrieval)

      APARS = 0.0_mpk ; LISTA = 0

      LISTA(1)    = 1
      APARS(1)    = AOD_firstguess
      CRetPar(1)  = '     Ret AOD     '
      CRetPar2(1) = 'AOD'

      LISTA(2)    = 2 
      APARS(2)    = Pkht_firstguess
      CRetPar(2)  = ' Ret Aero Pk Hgt '
      CRetPar2(2) = 'Aero Peak Hgt'

      ipars = 2

      if ( Retrieve_SO2 ) then
        ipars = ipars + 1
        LISTA(3)        = ipars
        APARS(ipars)    = SO2_firstguess
        CRetPar(ipars)  = '     Ret SO2     '
        CRetPar2(ipars) = 'SO2'
      endif

      if ( Retrieve_PlumeHWHM ) then
        ipars = ipars + 1
        LISTA(4)        = ipars
        APARS(ipars)    = plumehwhm_firstguess
        CRetPar(ipars)  = '  Ret Aero HWHM  '
        CRetPar2(ipars) = 'Aero HWHM'
      endif

      if ( Retrieve_NReal ) then
        ipars = ipars + 1
        LISTA(5)        = ipars
        APARS(ipars)    = nreal_firstguess(1)
        !APARS(ipars)   = nreal_firstguess(1) * 0.99999_mpk  !Perturbation test
        CRetPar(ipars)  = '    Ret Nreal    '
        CRetPar2(ipars) = 'Aero PSD Nreal'
      endif

      if ( Retrieve_NImag ) then
        ipars = ipars + 1
        LISTA(5) = ipars
        APARS(ipars)   = nimag_firstguess(1)
        !APARS(ipars)   = nimag_firstguess(1) * 1.001_mpk  !Perturbation test
        CRetPar(ipars)  = '    Ret Nimag    '
        CRetPar2(ipars) = 'Aero PSD Nimag'
      endif

      if (ipars /= npars) then
        stop 'Error: mismatch in building state vector.  Aborting ....'
      endif

      MA = npars

!  Display parameters being retrieved

      write(*,'(/1x,a)') 'Parameters being retrieved'
      do ipars=1,npars
        write(*,'(1x,a,i1,a)') 'Par ',ipars,' : '//TRIM(ADJUSTL(CRetPar2(ipars)))
      enddo
      write(*,*)

!  Plume settings
!  --------------

!  5/26/22. Region extension, allows to define plume peak-height anywhere between:
!           Ht(n1) - Window, and Ht(n2) + Window, Window = 5.0
!     Use a fine-layer vertical resolution between these two levels

      Region_extension = .false.
      if ( Region_extension ) then
! Plume anywhere 20-45 km
        n1 = 26     ! 26 = Level height 50 km
        n2 = 61     ! 61 = Level height 15 km
        U33 = 331 ; U88 = 881 ; U99 = 991  ! Output file units
      else
!  Plume generally 24-34 km
        !n1 = 39     ! 39 = Level height 37 km
         n1 = 42     ! 42 = Level height 34 km

        !n2 = 48     ! 48 = Level height 28 km
        !n2 = 50     ! 50 = Level height 26 km
         n2 = 52     ! 52 = Level height 24 km
        !n2 = 56     ! 56 = Level height 20 km
        !n2 = 61     ! 61 = Level height 15 km (below tropical tropopause)
        U33 = 330 ; U88 = 880 ; U99 = 991  ! Output file units
      endif

!  Set fine atmos resolution for plume

!mick test
      nfine = 4  ! 4 subdivisions = 0.25 km resolution
!     nfine = 2  ! 2 subdivisions = 0.5 km resolution
!     nfine = 1  ! 1 subdivision  = 1.0 km resolution (used for 151-layer atmos)

!  Check plume lies within this fine-resolution regime

      if ( plumetop.gt.real(76-n1,mpk) )  then
         stop ' plume top is outside appropriate fine-resolution height-grid'
      else if ( plumebot.lt.real(76-n2,mpk) )  then
         stop ' plume bottom is outside appropriate fine-resolution height-grid'
      endif

!  Old Code. Choice of O3 climo
!      month  = 1  ! January
!      o3zone = 8  ! 10-20 S latidude zone in the Labow climo
!      o3zone = 7  ! 20-30 S latidude zone in the Labow climo

!  Other fixed settings: Earth Radius

      eradius = 6371.0_mpk

!  Old Code. Start Scan loop
!   do iscan = iscan_start, iscan_end
!      which_radfile = 'rad.xxx.dat'
!      write(which_radfile(5:7),'(I3.3)')iscan
!  Screen shot and progress file
!      write(*,       '(/A,I3/)')'       @@@@@@@@@@@@@@@@@@@@@@@ DOING SCAN #', iscan
!      if ( nscans.gt.1 ) then
!        write(F33_unit,'(/A,I3/)')'       @@@@@@@@@@@@@@@@@@@@@@@ DOING SCAN #', iscan
!      endif

!  Headers
!mick mod 1/10/2023 - added bimodal fraction (wrt PSD #1)

      write(CBf,'(F6.4)')Bimodal_fraction
      write(CRg,'(F6.4)')psdpar1_firstguess(1)
      write(CSg,'(F5.3)')psdpar2_firstguess(1)
      write(CNr,'(F5.3)')NReal_firstguess(1)
      write(CNi,'(F6.4)')NImag_firstguess(1)
      write(CWLo,'(I3)') Int(Lambda_1)
      write(CWHi,'(I3)') Int(Lambda_2)
      Extraheader = '_Bf'//CBf//'_Rg'//CRg//'_Sg'//CSg//'_Nr'//CNr//'_Ni'//CNi//'_WLo'//CWLo//'_WHi'//CWHi
      IF ( Include_SO2 ) THEN
         Extraheader_Res = Trim(ExtraHeader)//'_IncSO2_'
      ELSE
         Extraheader_Res = Trim(ExtraHeader)//'_O3Only_'
      ENDIF

!  Data paths
!mick mod 10/31/2022 -  defined new data paths for flexibility in defining
!                       the input file names below for additional testing

      atmos_datapath   = 'orbit_data/atmos/'
      mie_datapath     = 'orbit_data/mie/'
      rad_datapath     = 'orbit_data/rad/'

      physics_datapath = 'physics_data/'

!  File names
!mick mod 7/31/2023 - modified this section for 2-band retrievals

      Test_option = 4 !Eun Su orbit set #2 (orbit overlap): TROPOMI 2-band studies

      if (Test_option == 4) then

        !New setup
        !=========

        !Normal sets
        !Orbit22086_RetType = '4' ; Orbit22087_RetType = '-1'
        !Orbit22086_RetType = '-1'; Orbit22087_RetType = '4'

!mick test
        !TROPOMI 2-band cases
        if ( (test>=1) .and. (test<=1) ) then
          !Orbit 22086
          Orbit22086_RetType = '5' ; Orbit22087_RetType = '-1'
        elseif ( (test>=2) .and. (test<=2) ) then
          !Orbit 22087
          Orbit22086_RetType = '-1'; Orbit22087_RetType = '5'
        endif

        !Tonga Ret FM sim cases
        !Orbit22086_RetType = '6' ; Orbit22087_RetType = '-1'
        !Orbit22086_RetType = '-1'; Orbit22087_RetType = '6'


        !Echo orbit test input & check
        if (chi_loop == 1) then
           write(*,*) 'trim(adjustl(Orbit22086_RetType)) = |'//trim(adjustl(Orbit22086_RetType))//'|'
           write(*,*) 'trim(adjustl(Orbit22087_RetType)) = |'//trim(adjustl(Orbit22087_RetType))//'|'
           !read(*,*)
        endif

        !Perform additional check
        if ( .not. ( trim(adjustl(Orbit22086_RetType)) /= '-1' .or. &
                     trim(adjustl(Orbit22087_RetType)) /= '-1' ) ) then
           write(*,'(/1x,a)') 'Error: either Orbit22086_RetType or Orbit22087_RetType should be -1.'
           write(*,'(1x,a/)') 'Stopping ....'
           stop
        endif

        !Save retrieval type for output file name
        if ( trim(adjustl(Orbit22086_RetType)) /= '-1' ) then
           RetType = Orbit22086_RetType
        elseif ( trim(adjustl(Orbit22087_RetType)) /= '-1' ) then
           RetType = Orbit22087_RetType
        endif

        !Orbit 22086
        !-----------

        if ( trim(adjustl(Orbit22086_RetType)) == '5' ) then
           !Ret Type #5  - using Haffner scene HPTZ file only
!mick test
           !RaySceneAtmosIndex  = '22086_22087_1' !76 layers
           !AeroSceneAtmosIndex = '22086_22087_1'

           !RaySceneAtmosIndex  = 'MLS' !76 layers
           !AeroSceneAtmosIndex = 'MLS'

           RaySceneAtmosIndex  = '22085_2071_03'
           AeroSceneAtmosIndex = '22086_2071_03'

           RadFilePrefix       = 'TROPOMI_' !TROPOMI test

           RaySceneRadIndex    = '22085_2071_03'
           AeroSceneRadIndex   = '22086_2071_03'

        endif

        !Orbit 22087
        !-----------
        if ( trim(adjustl(Orbit22087_RetType)) == '5' ) then
           !Ret Type #5  - using Haffner scene HPTZ file only
!mick test
           !RaySceneAtmosIndex  = '22086_22087_1' !76 layers
           !AeroSceneAtmosIndex = '22086_22087_1'

           !RaySceneAtmosIndex  = 'MLS' !76 layers
           !AeroSceneAtmosIndex = 'MLS'

           RaySceneAtmosIndex  = '22085_1963_73'
           AeroSceneAtmosIndex = '22087_1963_73'

           RadFilePrefix       = 'TROPOMI_' !TROPOMI test

           !RaySceneRadIndex    = '22085_1962_73'
           !AeroSceneRadIndex   = '22087_1962_73'

           RaySceneRadIndex    = '22085_1963_73'
           AeroSceneRadIndex   = '22087_1963_73'

        endif

        !Define retrieval input files
        !----------------------------
!mick test - set to use 2 bands
        !Number of spectral bands used
        numbands = maxbands

        !Atmos file prefix
        !HPTO3_file_prefix = 'make_HPTZ_' !default
        HPTO3_file_prefix = 'make_063_HPTZ_' !using M2SCREEM data

        !(1) Rayleigh scene

!mick test - set to use spatial smoothing
!          - set to use spectral smoothing x

        HPTO3_filename(1) = Trim(atmos_datapath)//Trim(HPTO3_file_prefix)//Trim(RaySceneAtmosIndex)//'.dat'
        do ib=1,numbands
          !CharRadAvg(ib)  = '01'
          if (ib == 1) then
            !CharRadAvg(ib)  = '09'
            CharRadAvg(ib)  = '03'
            CharRadBand(ib) = 'BD1'
          elseif (ib == 2) then
            !CharRadAvg(ib)  = '27'
            CharRadAvg(ib)  = '21'
            CharRadBand(ib) = 'BD2'
          endif
          which_radfile(1,ib) = Trim(rad_datapath)//Trim(RadFilePrefix)//CharRadAvg(ib)//'_'// &
                                CharRadBand(ib)//'_'//RaySceneRadIndex
        enddo

        !(2) Plume (Rayleigh + aerosol) scene

!mick test

        !For 2-band operation: general
        Initial_Mie_file  = Trim(mie_datapath)//'MIE_'//AeroSceneRadIndex//'_TROPOMI_2Band'

        !For 2-band operation: mode radius focus
        !Initial_Mie_file  = Trim(mie_datapath)//'MIE_'//AeroSceneRadIndex//'_TROPOMI_2Band'//'_'//CRg

        !For band 1 operation: mode radius focus
        !Initial_Mie_file  = Trim(mie_datapath)//'MIE_'//AeroSceneRadIndex//'_TROPOMI_Band1'//'_'//CRg

        !For band 1 operation: bimodal fraction focus
        !Initial_Mie_file  = Trim(mie_datapath)//'MIE_'//AeroSceneRadIndex//'_TROPOMI_Band1'//'_'//CBf

        !For band 1 operation: refractive index (NReal) focus
        !Initial_Mie_file  = Trim(mie_datapath)//'MIE_'//AeroSceneRadIndex//'_TROPOMI_Band1'//'_'//CNr

        !For band 1 operation: Eun-Su Yang set
        !Initial_Mie_file  = Trim(mie_datapath)//'MIE_'//AeroSceneRadIndex//'_Band1_Yang'

!mick test - set to use spatial smoothing
!          - set to use spectral smoothing x

        HPTO3_filename(2) = Trim(atmos_datapath)//Trim(HPTO3_file_prefix)//Trim(AeroSceneAtmosIndex)//'.dat'
        do ib=1,numbands
          !CharRadAvg(ib)  = '01'
          if (ib == 1) then
            !CharRadAvg(ib)  = '09'
            CharRadAvg(ib)  = '03'
            CharRadBand(ib) = 'BD1'
          elseif (ib == 2) then
            !CharRadAvg(ib)  = '27'
            CharRadAvg(ib)  = '21'
            CharRadBand(ib) = 'BD2'
          endif
          which_radfile(2,ib) = Trim(rad_datapath)//Trim(RadFilePrefix)//CharRadAvg(ib)//'_'// &
                                CharRadBand(ib)//'_'//AeroSceneRadIndex
        enddo

        !Define output file name extension
        FileScene = AeroSceneRadIndex

      endif

!mick mod 7/31/2023 - modified this section for 2-band retrievals

      !Index guide
      !  Rayleigh scene k=1
      !  Aerosol scene  k=2
      do k = 1, 2
        if (k == 2) write(*,'(/1x,a)') 'Mie file      = |'//Trim(Adjustl(Initial_Mie_file))//'|'
        write(*,'(/1x,a,i1,a)')        'HPTO3 file(',k,') = |'//Trim(Adjustl(HPTO3_filename(k)))//'|'
        do ib=1,numbands
          write(*,'(1x,a,2(i1,a))')    'Rad file(',k,',',ib,') = |'//Trim(Adjustl(which_radfile(k,ib)))//'|'
        enddo
      enddo

!mick test
!write(*,*) 'at rad file check'
!read(*,*)

!  Label to describe scene radiances

      !SceneRadsLabel = 'TROPOMI__22085_2071_03__22086_2071_03'
      !SceneRadsLabel = 'TROPOMI__22085_1962_73__22087_1962_73'
      !SceneRadsLabel = 'TROPOMI__' // RaySceneRadIndex // '__' // AeroSceneRadIndex
      SceneRadsLabel = Trim(RadFilePrefix) // '_' // RaySceneRadIndex // '__' // AeroSceneRadIndex

!  Set surface albedo
!mick mod 2/22/2024 - retrieval config file input now

      !Tonga aerosol scenes
      !albedo = 0.06_mpk

!  Data Reads
!  -----------

!  Old Code. data path for OMPS DATA
!      data_nwavs  = 151
!      rad_datapath = '../DHAFFNER_data_07Mar22/2022m0117_N20_NP_o21577/'

!  Radiance data from TROPOMI & others
!    -- wavelength grids the same
!mick mod 7/31/2023  - modified this section for 2-band retrievals
!mick mod 10/28/2024 - modified to also read and process SNR values of input radiances

      !TROPOMI files
      !data_nwavs = 497
      !OMI Simulator Case files
      !data_nwavs = 100
      !Tonga Retrieval FM Case files
      !data_nwavs = 30

      do k = 1, 2
        !Read band files for current scene
        do ib=1,numbands
          radfile = which_radfile(k,ib)
          open(ib,file=Trim(radfile),err=98,status='old')
          read(ib,*) data_nwavs_in(ib)
          read(ib,*) lat_in(k,ib), lon_in(k,ib), SZA_in(k,ib), VZA_in(k,ib), SAZ_in(k,ib), VAZ_in(k,ib)
          !write(*,*)
          !write(*,*) 'for scene ',k,' band ',ib
          !write(*,*) 'data_nwavs_in(ib) = ',data_nwavs_in(ib)
          !write(*,*) 'geos: '
          !write(*,*) lat_in(k,ib), lon_in(k,ib), SZA_in(k,ib), VZA_in(k,ib), SAZ_in(k,ib), VAZ_in(k,ib)
          !write(*,*) 'wvls & rads:'
          do w = 1, data_nwavs_in(ib)
             read(ib,*) data_wavs_in(w,ib), data_rads_in(w,k,ib), data_snrs_in(w,k,ib)
             !write(*,*) w, data_wavs_in(w,ib), data_rads_in(w,k,ib), data_snrs_in(w,k,ib)
          enddo
          close(ib)
        enddo
        !if (k==2) stop 'at spot input data check'

        !Contruct average data from files:
        !(1) location & geometry
        lat(k) = sum(lat_in(k,1:numbands))/real(numbands,mpk)
        lon(k) = sum(lon_in(k,1:numbands))/real(numbands,mpk)
        SZA(k) = sum(SZA_in(k,1:numbands))/real(numbands,mpk)
        VZA(k) = sum(VZA_in(k,1:numbands))/real(numbands,mpk)
        SAZ(k) = sum(SAZ_in(k,1:numbands))/real(numbands,mpk)
        VAZ(k) = sum(VAZ_in(k,1:numbands))/real(numbands,mpk)

        !(2) wvls & rads
        if (k==1) data_nwavs = sum( data_nwavs_in(1:numbands) )
        do ib=1,numbands
          if (ib==1) then
            w_off = 0
          else
            w_off = data_nwavs_in(ib-1)
          endif
          do w_in = 1, data_nwavs_in(ib)
            w = w_off + w_in
            if (k==1) then
              data_wavs(w) = data_wavs_in(w_in,ib)
              data_band_flag(w) = ib
            endif
            data_rads(w,k) = data_rads_in(w_in,k,ib)
            data_snrs(w,k) = data_snrs_in(w_in,k,ib)
          enddo
        enddo
      enddo

      !Check input data if desired
      !write(*,'(/1x,a,i4)') 'data_nwavs = ',data_nwavs
      !write(*,'(/5(1x,a))') ' w ','   data_wavs  ',' data_rads(1) ',' data_rads(2) ','data_band_flag'
      !write(*,'( 5(1x,a))') '---','--------------','--------------','--------------','--------------'
      !do w = 1, data_nwavs
      !   write(*,'(1x,i3.3,3(1x,es14.7),8x,i1)') w,data_wavs(w),data_rads(w,1:2),data_band_flag(w)
      !enddo
      !stop 'at input data check'

!  Step 4. Buffer the data
!mick note 3/6/2023  - define actual data & number of data used from the files
!mick mod  3/6/2023  - display that information
!                    - changed relational operator for Lambda_1 from ".gt." to ".ge."
!mick mod 7/31/2023  - define "band_flag"
!mick mod 10/28/2024 - added "snrs"

      do k = 1, 2
        wc = 0 
        do w = 1, data_nwavs
          if ( data_wavs(w).ge.lambda_1 .and. data_wavs(w).le.lambda_2 ) then
            wc = wc + 1 ; if ( wc.gt.E_Maxwavs ) stop 'Not enough points in E_MAXWAVS!'
            wavs(wc)      = data_wavs(w)
            band_flag(wc) = data_band_flag(w)
            rads(wc,k)    = data_rads(w,k)
            snrs(wc,k)    = data_snrs(w,k)
          endif
        enddo
        ndata = wc
      enddo

      if (chi_loop == 1) then
        write(*,'(/1x,a,2(f7.3,a))') 'spectral range of rad files: [',data_wavs(1),',',data_wavs(data_nwavs),']'
        write(*,'(1x,a,i4)')         '# of wvls in rad files     :  ',data_nwavs
        write(*,'(1x,a,2(f7.3,a))')  'spectral range used        : [',lambda_1,',',lambda_2,']'
        write(*,'(1x,a,i4)')         '# of wvls used             :  ',ndata 
!        read(*,*)
      endif

!  Define the ratioed measurements, their SNR, and estimated noise
!mick mod 9/22/2022  - modified code to more formally create the two noise options:
!                      (1) one to use fixed measurement noise (used for early testing)
!                      (2) one to use fixed SNR (used more now)
!mick mod 7/31/2023  - modified this section for 2-band retrievals: added 3rd noise option
!mick mod 10/28/2024 - added noise option #4

      do w = 1, ndata
        measrads(w) = rads(w,2) / rads(w,1)
      enddo

      !NoiseOption = 1 !use fixed measurement noise
      !NoiseOption = 2 !use fixed SNR
      !NoiseOption = 3 !use fixed SNRs (2 bands)
      NoiseOption = 4 !use variables SNRs (2 bands)

      !Note: "snr"       - SNR of individual rads (info only)
      !      "snr_ratio" - SNR of the rad ratios  (info only)
      !      "sigs_orig" - std dev of rad ratio noise (used in retrieval)
      if (NoiseOption .eq. 1) then
        !"error" - assumed error of rad ratios (delQ)
        error              = 0.1_mpk

        snr(1:ndata)       = (sqrt(2.0_mpk) * abs(measrads(1:ndata))) / error
        snr_ratio(1:ndata) = abs(measrads(1:ndata)) / error
        sigs_orig(1:ndata) = error

      elseif (NoiseOption .eq. 2) then
        !"signal_to_noise" - assumed SNR of individual rads
        signal_to_noise    = 150.0_mpk !Initial suggestion by Dave Haffner
        !signal_to_noise    = 53.3_mpk  !Used for helping estimate Fwd Model Error

        snr(1:ndata)       = signal_to_noise
        snr_ratio(1:ndata) = signal_to_noise / sqrt(2.0_mpk)
        sigs_orig(1:ndata) = (sqrt(2.0_mpk) * abs(measrads(1:ndata))) / signal_to_noise

      elseif (NoiseOption .eq. 3) then
        !"signal_to_noise" - assumed SNR of individual rads
        signal_to_noise    = 150.0_mpk !Initial suggestion by Dave Haffner
        !signal_to_noise    = 53.3_mpk  !Used for helping estimate Fwd Model Error

        do ib=1,numbands
          !boost the original SNR by spatial averaging done
          read(CharRadAvg(ib),'(i2.2)') NumRadAvg
          band_signal_to_noise = sqrt(real(NumRadAvg,mpk))*signal_to_noise

          !find number of wvls in current wvl subset being used in the retrieval
          !to whom this band SNR applies
          nbanddata(ib) = 0
          do w =1, ndata 
            if (band_flag(w) == ib) nbanddata(ib) = nbanddata(ib) + 1
          enddo

          !define the bounds of this wvl subset within the "(1:ndata)" wvl set
          if (ib==1) then
            w0 = 1; w1 = nbanddata(ib)
          else
            w0 = nbanddata(ib-1) + 1; w1 = nbanddata(ib-1) + nbanddata(ib)
          endif

          !define the SNR and error quantities for this wvl subset
          snr(w0:w1)       = band_signal_to_noise
          snr_ratio(w0:w1) = band_signal_to_noise / sqrt(2.0_mpk)
          sigs_orig(w0:w1) = (sqrt(2.0_mpk) * abs(measrads(w0:w1))) / band_signal_to_noise
        enddo

      elseif (NoiseOption .eq. 4) then
        do ib=1,numbands
          !find number of wvls in current wvl subset being used in the retrieval
          !to whom this band SNR applies
          nbanddata(ib) = 0
          do w =1, ndata 
            if (band_flag(w) == ib) nbanddata(ib) = nbanddata(ib) + 1
          enddo

          !define the bounds of this wvl subset within the "(1:ndata)" wvl set
          if (ib==1) then
            w0 = 1; w1 = nbanddata(ib)
          else
            w0 = nbanddata(ib-1) + 1; w1 = nbanddata(ib-1) + nbanddata(ib)
          endif

          !boost the original SNR by spectral averaging done
          !  --> ASSUME, for the moment, this has already been done in the calibrated-rad file pre-processing
          !read(CharRadAvg(ib),'(i2.2)') NumRadAvg
          !snrs(w0:w1,:) = sqrt(real(NumRadAvg,mpk))*snrs(w0:w1,:)

          !define the SNR and error quantities for this wvl subset
          snr_ratio(w0:w1) = sqrt( (snrs(w0:w1,1)**2 * snrs(w0:w1,2)**2) / (snrs(w0:w1,1)**2 + snrs(w0:w1,2)**2 ) )
          sigs_orig(w0:w1) = abs(measrads(w0:w1)) / snr_ratio(w0:w1)
        enddo

      endif

!  Check effective SNR if desired
!mick mod 10/28/2024 - added noise option #4

      if (NoiseOption.ge.1 .and. NoiseOption.le.3) then
        !write(*,'(/5(1x,a))') ' w ','   measrads   ','      snr     ','  snr_ratio   ','     sigs     '
        !write(*,'( 5(1x,a))') '---','--------------','--------------','--------------','--------------'
        !do w = 1, ndata
        !   write(*,'(1x,i3.3,4(1x,es14.7))') w,measrads(w),snr(w),snr_ratio(w),sigs_orig(w)
        !enddo

        if (NoiseOption .eq. 1) then
          avg_snr = sum(snr(1:ndata))/real(ndata,mpk)
          write(*,'(/1x,a,es14.7)') 'avg snr of individual rads = ',avg_snr
          avg_snr_ratio = sum(snr_ratio(1:ndata))/real(ndata,mpk)
          write(*,'(/1x,a,es14.7)') 'avg snr of individual rads = ',avg_snr_ratio
        endif

      elseif (NoiseOption .eq. 4) then
        write(*,'(/7(1x,a))') ' w ','      wvl     ','   measrads   ','    snr(1)    ','    snr(2)    ',&
                                    '  snr_ratio   ','     sigs     '
        write(*,'( 7(1x,a))') '---','--------------','--------------','--------------','--------------',&
                                    '--------------','--------------'
        do w = 1, ndata
           write(*,'(1x,i3.3,6(1x,es14.7))') w,wavs(w),measrads(w),snrs(w,1),snrs(w,2),snr_ratio(w),sigs_orig(w)

        enddo

        !special rad/snr/sig file output
        open(45,file='rad_snr_sig_table.dat',action='write',status='replace')
        write(45,'(/7(1x,a))') ' w ','      wvl     ','   measrads   ','    snr(1)    ','    snr(2)    ',&
                                     '  snr_ratio   ','     sigs     '
        write(45,'( 7(1x,a))') '---','--------------','--------------','--------------','--------------',&
                                     '--------------','--------------'
        do w = 1, ndata
           write(45,'(1x,i3.3,6(1x,es14.7))') w,wavs(w),measrads(w),snrs(w,1),snrs(w,2),snr_ratio(w),sigs_orig(w)
        enddo
        close(45)

      endif
      !stop 'at measurement error/SNR check stop'

!  Redefine the std devs of Seps (i.e. the "sigs"), based initially on the std devs of
!  measurement noise alone, to include std devs related to estimated forward model error
!  if required
!mick mod 1/10/2023 - added this section

      sigs(1:ndata) = sigs_orig(1:ndata) + sigF(1:ndata)

      !Check effect of adding std devs of estimated forward model error if desired
      !write(*,'(/5(1x,a))') ' w ','   sigs_orig  ','     sigF     ','     sigs     ','  Rel % Diff  '
      !write(*,'( 5(1x,a))') '---','--------------','--------------','--------------','--------------'
      !do w = 1, ndata
      !   write(*,'(1x,i3.3,4(1x,es14.7))') w,sigs_orig(w),sigF(w),sigs(w),(sigs(w)/sigs_orig(w) - one)*100.0_mpk
      !enddo
      !stop 'at sigF check stop'

!  Pre-screen scene for aerosol using the radiance ratios.
!  If aerosol present, proceed; if not, exit retrieval.
!mick mod 12/2/2022 - added this section
!mick mod 3/6/2023  - added method #2

      !method #1: avg rad ratio of band
      !num_idx = idx_hi - idx_lo + 1
      !CSI = sum(measrads(idx_lo:idx_hi))/real(num_idx,mpk)

      !method #2: rad ratio at 1st wvl at low end of band
      CSI = measrads(idx_lo)

      no_aerosol = .false.
      if (CSI < CSI_threshold) then
        !CSI too close to 1.0: aerosol not present --> exit
        no_aerosol = .true.
        alamda = 0.0_mpk
        go to 455
      endif

!  Tonga retrieval input check stop

   !write(*,*)
   !write(*,*) 'At Tonga retrieval input check stop: hit enter to continue'
   !read(*,*)

!  Initialize atmospheric properties input structure "HPTORA_Global"
!  -----------------------------------------------------------------
!mick mod 9/22/2022 - added this subroutine
!mick mod 1/24/2023 - added loop

      do i=1,2
         CALL TONGA_Initialize_HPTORA_V5 ( HPTORA_Global(i) )
      enddo

!  Initialize Mie inputs
!  ---------------------
!mick mod 12/9/2022 - added this subroutine to set up Bimodel PSD
!mick mod 8/17/2023 - removed input "npars"
!                   - added inputs "retrieve_nreal" and "retrieve_nimag"

      CALL TONGA_Mie_Initialize &
         ( retrieve_nreal, retrieve_nimag, AOD_refw, & !Inputs
           PSDpar1_firstguess, PSDpar2_firstguess,   & !Inputs
           nreal_firstguess, nimag_firstguess,       & !Inputs
           Mie_Inputs_Global )                         !Outputs

!  Initial Mie calls to store data
!  -------------------------------

      if ( Do_Initial_Mie_Calculations .and. .not.(retrieve_nreal .or. retrieve_nimag )) then

!  Open file

         open(44,file=Trim(Initial_Mie_file),status = 'replace')

!  First call for extinction coefficient at the Reference wavelength
!mick mod 12/9/2022 - replace std Mie call with Bimodal one

         write(*,'(/1x,a)') 'Doing Mie calculation for ref wavelength'

         !call RTSMie_Master &
         !   ( Mie_Inputs_Global, Mie_Outputs, & ! I/O type structures
         !     Mie_fail, Mie_istatus, Mie_messages(1), Mie_messages(2), Mie_messages(3) )  ! Exception handling
         !if ( Mie_fail ) go to 677

         call RTSMie_master_bimodal &
            ( Mie_Inputs_Global(1), Mie_Inputs_Global(2), Bimodal_fraction, & ! Inputs
              Mie_Outputs, Mie_fail, Mie_istatus, Mie_messages, Mie_trace_3 ) ! Outputs
         if ( Mie_fail ) go to 677

!  Write to file
!mick note 12/9/2022 - Mie_Bulk(1): extinction cross-section in um^2
!mick mod 1/24/2023  - apply to HPTORA_Global 2nd element

         HPTORA_Global(2)%Aer_Refw_Extinction = Mie_Outputs%Mie_Bulk(1)
         write(44,'(/A/)')'  Aerosol extinction at reference wavelength'
         write(44,'(f12.5,1pe24.12)')AOD_refw,HPTORA_Global(2)%Aer_Refw_Extinction
         write(44,'(/A,I4,A)')'  Full aerosol properties for ',ndata,' wavelengths'

!  Main Mie data loop

         do w = 1, ndata

!  Need coefficients as well

            Mie_Inputs_Global(1:2)%lambda = wavs(w) * 0.001_mpk  ! Microns. Set later on, according to UV wavelengths
            Mie_Inputs_Global(1:2)%do_Expcoeffs = .true.         ! Will be set true for the main Mie calculations

!  Main call at each wavelength
!mick mod 12/9/2022 - replace std Mie call with Bimodal one

            write(*,'(/1x,a,es24.12)') 'Doing Mie calculation for wavelength ',wavs(w)

            !call RTSMie_Master & 
            !   ( Mie_Inputs_Global, Mie_Outputs, & ! I/O type structures
            !     Mie_fail, Mie_istatus, Mie_messages(1), Mie_messages(2), Mie_messages(3) )  ! Exception handling
            !if ( Mie_fail ) go to 678

            call RTSMie_master_bimodal &
               ( Mie_Inputs_Global(1), Mie_Inputs_Global(2), Bimodal_fraction, & ! Inputs
                 Mie_Outputs, Mie_fail, Mie_istatus, Mie_messages, Mie_trace_3 ) ! Outputs
            if ( Mie_fail ) go to 678

!  Use moment cutoff

            cutoff = 1.0d-8
            nmoms = Mie_Outputs%Mie_NCoeffs ; L = 0 ; mom = 1.0d-1
            do while ( L.lt.nmoms .and. abs(mom).gt.cutoff ) 
               L =  L + 1 ; mom = Mie_Outputs%Mie_Expcoeffs(1,L)
            enddo
            Mie_Outputs%Mie_NCoeffs = L-1
            Mie_Outputs%Mie_Expcoeffs(:,L:max_Mie_angles) = d_zero

!  Write to file

            write(44,'(/A,F12.5/)')' Mie data at wavelength ', wavs(w)
            write(44,'(1p4e24.12)')Mie_Outputs%Mie_Bulk(1:3),Mie_Outputs%Mie_Asymm
            write(44,'(I5)')Mie_Outputs%Mie_NCoeffs
            do L = 0, Mie_Outputs%Mie_NCoeffs
               write(44,'(I5,1p6e24.12)')L,Mie_Outputs%Mie_Expcoeffs(1:6,L)
            enddo
 
!  End Mie data loop

         enddo

!  Finish Mie calculation (close file and stop program)

         close(44)
         stop ' -- Finished Mie calculations '
      endif

!  Main program
!  ============

!  Initialize VLIDORT
!  ------------------

!  Initialize with 2 geometries in the global input file

      CALL TONGA_RTInitialize_V2 &
         ( SZA, VZA, SAZ, VAZ, nstokes, nstreams, nssfine,     & ! Inputs
           do_ssonly, eradius, Surface_albedo, npars,          & ! Inputs
           VFixIn_Global, VModIn_Global, VSup_Global,          & ! Outputs
           VLinFixIn_Global, VLinModIn_Global, VLinSup_Global, & ! Outputs
           fail, messages(1) )                                   ! Exception handling

      if ( fail ) then
         write(*,*)'RT Initialization V2failed - here is the error message ---'
         write(*,*)'  ---- ',Trim(messages(1))
         write(*,*)'stop program !!'; stop
      endif

!  Open output files
!mick mod 1/10/2023 - modified ret result file names to increase info content
!                     while trying to still limit length of string
!                   - simplified file name setup some
!                   - added capability to add "RetType" to file name

      if ( npars.eq.2 ) then
        FilePrefix = 'RetRes_2par_AOD_PKHT'
      else if ( npars.eq.3 ) then
        FilePrefix = 'RetRes_3par_AOD_PKHT_SO2'
      endif

      if (Test_option == 4) then
         !Add on "RetType" for option #3
         FileSceneExt = FileScene // '_RType' // trim(adjustl(RetType))

         !For atmos background file tests
         !write(CharTestCfg,'(i1)') TestCfg
         !FileSceneExt = FileScene // '_RType' // trim(adjustl(RetType)) // '_Cfg' // CharTestCfg
      else
         !Just copy "FileScene"
         FileSceneExt = FileScene
      endif

      !Retrieval iteration log
      Open(U33,file=trim(FilePrefix)//trim(ExtraHeader_Res)//trim(adjustl(FileSceneExt))// &
                       '.LOG',status = 'unknown')
      !Retrieval results
      Open(U88,file=trim(FilePrefix)//trim(ExtraHeader_Res)//trim(adjustl(FileSceneExt))// &
                       '.DAT',status = 'unknown')
      !Retrieval spectra for plotting
      Open(U99,file=trim(FilePrefix)//trim(ExtraHeader_Res)//trim(adjustl(FileSceneExt))// &
                       '.PLT',status = 'unknown')

!  Screen shot and progress file

      write(*,  '(/A,I3/)')'       @@@@@@@@@@@@@@@@@@@@@@@ DOING '//trim(adjustl(SceneRadsLabel))
      write(U33,'(/A,I3/)')'       @@@@@@@@@@@@@@@@@@@@@@@ DOING '//trim(adjustl(SceneRadsLabel))

!  Create the initial properties (HPTORA_Global) using buffering
!    5/31/22. Use dedicated aerosol file.
!    8/31/22. Introduce SO2 control
!mick mod 1/24/2023  - added loop
!mick mod 8/17/2023  - replaced "SO2_totalcolumn_DU" with "SO2_firstguess"
!                    - replaced "plumehwm" with "plumehwhm_firstguess"
!                    - added "retrieve_SO2"
!mick mod 2/22/2024  - moved to V5; added inputs "Surface_height" & "Surface_albedo"
!mick fix 10/28/2024 - added scene index "i" to inputs

      do i=1,2
         write(*,'(/1x,a,i1)') 'For scene ',i !1=Rayleigh scene; 2=Plume scene
         CALL TONGA_Create_HPTORA_V5 &
            ( i, ndata, wavs, n1, n2, nfine, Surface_height, Surface_Albedo, &
              Include_SO2, SO2_firstguess, SO2XSec_Source, retrieve_SO2, &
              HPTO3_filename(i), physics_datapath, &
              plumetop, plumebot, plumehwhm_firstguess, HPTORA_Global(i), fail, messages(1) )

         if ( fail ) then
            write(*,*)'Create the initial properties (HPTORA_Global) using buffered data, failed - here is the error message ---'
            write(*,*)'  ---- ',Trim(messages(1))
            write(*,*)'stop program !!'; stop
         endif
      enddo

!  Remaining settings for the aerosols
!mick mod 1/24/2023 - apply to HPTORA_Global 2nd element
!mick mod 8/17/2023 - added "retrieve_plumehwhm"

      HPTORA_Global(2)%Aer_Refw           = AOD_Refw
      HPTORA_Global(2)%Retrieve_PlumeHWHM = retrieve_plumehwhm
      HPTORA_Global(2)%Retrieve_NReal     = retrieve_nreal
      HPTORA_Global(2)%Retrieve_NImag     = retrieve_nimag

!  RETRIEVAL
!  =========

!  Inverse model Levenberg-Marquardt iteration
!  iteration starting values

      J = 0 ;  ALAMDA = -0.1_mpk ; trawl = .true. ; ochisq = 1.0d+30 ; convcount = 0

!  File Headers and initial guess write-up
!mick mod 8/17/2023 - modified to display up to five pars

      CHead = '  Iteration  CHANGE_CHISQ         CHISQ       ALAMDA        '
      do ipars=1,npars
        CHead = Trim(CHead)//CRetPar(ipars)
      enddo
      !write(*,'(1x,a)') '|'//Trim(CHead)//'|'

      write(*,  '(/A//4x,i2,36x,e10.2,5(2x,F13.8,2x))')Trim(CHead), j, alamda, APARS(1:npars)
      write(U33,'(/A//4x,i2,36x,e10.2,5(2x,F13.8,2x))')Trim(CHead), j, alamda, APARS(1:npars)

!  Inverse model Levenberg-Marquardt iteration
!  -------------------------------------------

!  Start inverse loop

      do while ( trawl .and. j .lt. MAXITER )

!  Increase iteration count by 1 and call the inverse model
!mick mod 9/22/2022  - enhanced error handling
!mick mod 12/9/2022  - To accomodate a Bimodal PSD, added "Bimodal_fraction" and "Mie_trace_3"
!mick mod 1/10/2023  - added output SIMJAC
!mick mod 10/28/2024 - added output SO2TAU

         j = j + 1 ; do_Dump = .false. ; if ( j.eq.1 ) do_Dump = parameter_do_Dump
         APARS_OLD = APARS ; old_chisq = ochisq
         CALL TONGA_MRQMIN ( J, E_MAXWAVS, MAXA, NDATA, MA, WAVS, MEASRADS, SIGS,               & !Fwd model Inputs
                             do_ssonly, do_Dump, HPTORA_Global,                                 & !Fwd model Inputs + HPTORA
                             Use_Initial_Mie_Calculations, Initial_Mie_File,                    & !Fwd model Mie inputs
                             Mie_Inputs_Global, Bimodal_fraction,                               & !Fwd model Mie inputs
                             VFixIn_Global, VModIn_Global, VSup_Global,                         & !Fwd model VLIDORT input
                             VLinFixIn_Global, VLinModIn_Global, VLinSup_Global,                & !Fwd model VLIDORT input  
                             APARS, LISTA, COVAR, ALPHA, CHISQ, ALAMDA, SIMRADS, SIMJAC,        & !Inv model Input/Output
                             OCHISQ, ATRY, BETA, DA, MFIT, SO2TAU, fail, messages, Mie_trace_3 )  !Inv model Output+Status

!  Retrieval error: skip to 455

         IF ( FAIL ) THEN
            write(*,'(/A,I3)')' TROPOMI double-measurement scans failed'
            write(*,'(A,I2,A)')'   TONGA_MRQMIN failed for iteration # ',J,'; - here is the error message(s) ---'
            do i=1,3
              write(*,*)'  ---- ','    **** '//Trim(messages(i))
            enddo
            GO TO 455
!            write(*,*)'stop program !!'; stop
         ENDIF

!  Save simulations

         SIMRADS_SAVED(J,1:NDATA) = SIMRADS(1:NDATA)

!  Determine convergence
!mick note 9/22/2022 - define the magnitude of the relative reduction in the objective function
!                      (i.e. cost function) to be minimized
!mick mod 8/17/2023  - modified to display up to five pars
!mick note 9/1/2023  - convergence test: compare the relative reduction with the user-defined
!                      convergence threshold
!mick mod 9/1/2023   - modified convergence criteria

         !Compute convergence parameters

         !part of old criteria
         chisq_change = abs ((chisq/old_chisq) - 1.0_mpk )
         !part of new criteria (hold)
         !chisq_change = (chisq/old_chisq) - 1.0_mpk
         qc = 0
         do q = 1, MA
            if ( abs ((APARS(q)/APARS_OLD(q)) - 1.0_mpk ) .lt. par_conv ) qc = qc + 1
         enddo

         !Write-up of information
         write(*  ,'(2x,2i2,1x,2f18.8,1x,e10.2,5(2x,F13.8,2x))') j, qc, chisq_change, chisq, alamda, APARS(1:npars)
         write(U33,'(2x,2i2,1x,2f18.8,1x,e10.2,5(2x,F13.8,2x))') j, qc, chisq_change, chisq, alamda, APARS(1:npars)

         !Chisq and parameter convergence

         !part of old criteria
         if ( chisq_change .lt. chisq_conv ) then
           !if ( npars.eq.2 .and. qc.eq.2 ) convcount = convcount + 1
           !if ( npars.eq.3 .and. qc.eq.3 ) convcount = convcount + 1
           if ( npars.eq.qc ) convcount = convcount + 1

           if ( convcount .eq. 2 ) then
             JSAV = J
             trawl = .false. ; write(*,*)'   ---- Converged'; write(U33,*)'   ---- Converged'
           endif
         endif

         !part of new criteria (hold)
         !if (chisq .lt. old_chisq) then
         !  if ( qc .eq. npars ) then
         !    JSAV = J
         !    trawl = .false. ; write(*,*)'   ---- Converged'; write(U33,*)'   ---- Converged'
         !  endif
         !endif

!  do not use Chisq as a convergence criterion
!         if ( CHISQ .lt. CONV ) trawl = .false.

!  End inverse-model iteration

      enddo

!  One more call to get the covariance matrix (set ALAMDA = 0.0_mpk )
!mick mod 9/22/2022  - enhanced error handling
!mick mod 12/9/2022  - To accomodate a Bimodal PSD, added "Bimodal_fraction" and "Mie_trace_3"
!mick mod 1/10/2023  - added output SIMJAC
!mick mod 10/28/2024 - added output SO2TAU

      ALAMDA = 0.0_mpk
      CALL TONGA_MRQMIN ( J, E_MAXWAVS, MAXA, NDATA, MA, WAVS, MEASRADS, SIGS,               & !Fwd model Inputs
                          do_ssonly, do_Dump, HPTORA_Global,                                 & !Fwd model Inputs + HPTORA
                          Use_Initial_Mie_Calculations, Initial_Mie_File,                    & !Fwd model Mie inputs
                          Mie_Inputs_Global, Bimodal_fraction,                               & !Fwd model Mie inputs
                          VFixIn_Global, VModIn_Global, VSup_Global,                         & !Fwd model VLIDORT inputs
                          VLinFixIn_Global, VLinModIn_Global, VLinSup_Global,                & !Fwd model VLIDORT inputs  
                          APARS, LISTA, COVAR, ALPHA, CHISQ, ALAMDA, SIMRADS, SIMJAC,        & !Inv model Input/Output
                          OCHISQ, ATRY, BETA, DA, MFIT, SO2TAU, fail, messages, Mie_trace_3 )  !Inv model Output+Status

!  Retrieval error: skip to 455

      if ( fail ) then
         write(*,'(/A,I3)')' TROPOMI double-measurement scans failed'
         write(*,'(A,I2,A)')'   TONGA_MRQMIN failed on last call to TONGA_MRQMIN - here is the message(s) ---'
         do i=1,3
           write(*,*)'  ---- ','    **** '//Trim(messages(i))
         enddo
         GO TO 455
!         write(*,*)'stop program !!'; stop
      endif

!  Set the error and cross-correlations
!mick mod 8/17/2023 - modified for up to "npars" pars

      SDPARS = 0.0_mpk ; CROSSCOR = 0.0_mpk
      do i=1,npars
        SDPARS(i) = SQRT(COVAR(i,i))
      enddo
      do i=1,npars
        do j=1,npars
          if (i .eq. j) then
            CROSSCOR(i,j) = 1.0_mpk
          else
            CROSSCOR(i,j) = COVAR(i,j) / ( SDPARS(i)*SDPARS(j) )
          endif
        enddo
      enddo

!  Continuation point for skipped output

455   continue

!  Write results
!  =============
!mick mod 8/15/2022 - modified writes & formats some to bring writes for both success & fail cases into agreement
!mick mod 12/2/2022 - added "no aerosol" section

!  Main results

      VZERO = 0.0_mpk
      if ( no_aerosol ) then
         write(U88,222)trim(adjustl(SceneRadsLabel)),0,.false.,0.000,0.000,alamda, &
                       (VZERO(i),VZERO(i),i=1,npars)
      elseif ( fail ) then
         write(U88,222)trim(adjustl(SceneRadsLabel)),j,.not.trawl,0.000,0.000,alamda, &
                       (VZERO(i),VZERO(i),i=1,npars)
      else
         write(U88,222)trim(adjustl(SceneRadsLabel)),jsav,.not.trawl,chisq_conv,chisq,alamda, &
                       (APARS(i),SDPARS(i),i=1,npars)
         write(U88,*)
         do i=1,npars
            write(U88,223) (CROSSCOR(i,j),j=1,npars)
         enddo
      endif
222   format(A,2x,i3,2x,L2,2x,2F16.8,e10.2,5(3x,F13.8,1x,F13.8))
223   format(5(1x,F13.8))

!  Simulations (to fort.400)
!mick mod 12/2/2022 - added "no aerosol" condition
!mick mod 1/10/2023 - turned off write to unit=400 for now (see chi-square section
!                     below where similar results are output (plus more)

      !if ( .not.(no_aerosol .or. fail) ) then
      !   do w = 1, NDATA
      !      write(400,*)w,wavs(w),measrads(w),simrads_saved(1:JSav,w)
      !   enddo
      !endif

!  Close files

      Close(U33)
      Close(U88)

!  Display timing for this retrieval pass
!mick mod 5/9/2023  - added passing of ret timing to timing output file

      call cpu_time(e2)
      write(*,*)
      write(*,'(1x,a,f6.2,a)')    'retrieval time = ',e2-e1,' sec'
      write(U100,'(1x,a,f6.2,a)') 'retrieval time = ',e2-e1,' sec'

!  Check chi-square of retrieval.  If out of range, estimate forward model error (epsF)
!  and std dev (sigF), go back to top of driver, and do retrieval again - this time
!  including non-zero std dev (sigF)
!mick mod 1/10/2023  - added this section & end of loop
!mick note 1/10/2023 - numbers to the right of program lines for a UV range of [289.0,298.0]
!mick fix 8/30/2023  - close U99 also for the ELSE case below
!mick mod 10/28/2024 - added output of simple (Column-Jac-based) SO2 AMFs

      rndata = real(ndata,mpk) !139
      deg_freedom = rndata - one !138
      std_dev = sqrt(two*deg_freedom) !~16.6
      chisq_lower_lim = deg_freedom - three*std_dev !~88.2
      chisq_upper_lim = deg_freedom + three*std_dev !~187.8

      write(*,*)
      write(*,*) 'chi_loop        = ',chi_loop
      write(*,*) 'rndata          = ',rndata
      write(*,*) 'deg_freedom     = ',deg_freedom
      write(*,*) 'std_dev         = ',std_dev
      write(*,*) 'chisq           = ',chisq
      write(*,*) 'chisq_lower_lim = ',chisq_lower_lim
      write(*,*) 'chisq_upper_lim = ',chisq_upper_lim

      !if ( no_aerosol .or. fail .or. &
      !    (chisq_lower_lim < chisq) .and. (chisq < chisq_upper_lim) .or. &
      !     chi_loop.eq.2 ) then
      !if ( no_aerosol .or. fail .or. &
      !    (chisq_lower_lim < chisq) .and. (chisq < chisq_upper_lim) ) then
      if ( no_aerosol .or. fail .or. (chisq < chisq_upper_lim) ) then

         write(*,'(/1x,a,i1)') 'Retrieval done at chi_loop = ',chi_loop

         !Only write this output if retrieval was performed and successful
         if ( .not.(no_aerosol .or. fail) ) then
            if ( Retrieve_SO2 ) then
              !define simple SO2 AMFs

              !write(*,'(/1x,a)') 'SO2TAU check:'
              !do w = 1, ndata
              !  write(*,*) 'w = ',w,' SO2TAU(w) = ',SO2TAU(w)
              !enddo
              SO2AMF(1:ndata) = -(APARS(3)*simjac(3,1:ndata))/(SO2TAU(1:ndata)*simrads(1:ndata))

              !write output w/SO2 AMFs
              do w = 1, ndata
                 write(U99,'(1x,i3,11es14.5e3)')w,wavs(w),measrads(w),simrads(w),simjac(1:MA,w),&
                                                sigs_orig(w),sigF(w),sigs(w),epsF(w),SO2AMF(w)
              enddo
            else
              !write standard output
              do w = 1, ndata
                 write(U99,'(1x,i3,10es14.5e3)')w,wavs(w),measrads(w),simrads(w),simjac(1:MA,w),&
                                                sigs_orig(w),sigF(w),sigs(w),epsF(w)
              enddo
            endif
         endif
         Close(U99)

         exit !exit chi-square check loop
      else
         chisq_mu = (rndata - one)/rndata !~0.9928
         do w = 1, ndata
           varF = (measrads(w) - simrads(w))**2/chisq_mu - sigs(w)**2

           !if (w == 1) write(*,*)
           !write(*,*) 'w = ',w
           !write(*,*) '(measrads(w) - simrads(w))**2/chisq_mu = ',(measrads(w) - simrads(w))**2/chisq_mu
           !write(*,*) 'sigs(w)**2                             = ',sigs(w)**2
           !write(*,*) 'varF                                   = ',varF

           !Case where forward model error is within envelope
           !of measurement noise, so zeroize it
           if (varF < zero) then
             varF = zero
             !write(*,*) 'varF modified                          = ',varF
           endif

           sigF(w) = sigF(w) + sqrt(varF)
           epsF(w) = sigF(w)/simrads(w)
         enddo
         Close(U99)
      endif

!  End Chi-Square Check Loop (this loop starts @ ~line 500)

      ENDDO

!  Close timing file

      Close(U100)

!  Finish
!mick fix 9/22/2022 - added error section
!mick mod 12/2/2022 - added "no aerosol" section

      if ( no_aerosol ) then
        write(*,*)' **** Retrieval message: pre-screening indicates no aerosol, exiting ...'
      elseif (fail) then
        write(*,*)' **** Retrieval error: see above message(s)'
      else
        write(*,*)' **** Retrieval successful finish'
      endif
      stop

!  Radfile failure stop

98    continue
      write(*,*)' **** Radfile not found, stop !!!!!!!!!!!!!!!' ; stop
      stop

!  Mie failure stops

677   Continue
      close(44)
      write(*,*)' Failure from initial Mie calculation at AOD_refw, here are the messages --'
      do j = 1, 3
         write(*,'(A,I2,A)')'  - Message # ', j,' : '//Trim(Mie_messages(j))
      enddo
      write(*,'(A)')'  - Message #  4 : '//Trim(Mie_trace_3)
      stop

678   Continue
      close(44)
      write(*,*)' Failure from Main Mie calculation, here are the messages --'
      do j = 1, 3
         write(*,'(A,I2,A)')'  - Message # ', j,' : '//Trim(Mie_messages(j))
      enddo
      write(*,'(A)')'  - Message #  4 : '//Trim(Mie_trace_3)
      stop

!  End program

end program TONGA_Retrieval_MeasRatio_V4

