Module TONGA_Create_HPTORA_V5_m

!  5/31/22. Add subroutine for Tests with different aerosol properties
!  7/28/22. Adapted for the V2 retrievals
!  8/31/22. Introduce optional SO2 (total DU from TROPOMI, same plume location as aerosol)
!mick mod 9/22/2022 - added additional TONGA_pars variable control
!                   - added subs "TONGA_Initialize_HPTORA" & "TONGA_Write_HPTORA"

!Robmod 1/25/2024. - Add surface height as input, and regrid original PTHO3 to this height
!                  - Regridding is done before the aerosol plume finelayer regridding
!                  - initial surface-height regridding is always below the aerosol plume
!                  - Add Surface height and albedo to HPTORA type structure (now V5)

    USE TONGA_pars_m, Only : MPK, E_MAXLAYERS, E_MAXPARS, E_MAXWAVS, &
                             DWFL, DWFL1, DWFI, DWFI1, DWFR, DWFR1, DWFR2, &
                             ZERO
    USE TONGA_HPTORA_V5_m

!private subs - Get_HPTO3_TROPOMI, GET_O3XSEC_3, GET_SO2XSEC_2, GET_SO2XSEC_3, &
!               RAYLEIGH_FUNCTION, SPLINE, SPLINT

PRIVATE
PUBLIC  :: TONGA_Initialize_HPTORA_V5, TONGA_Create_HPTORA_V5, TONGA_Write_HPTORA_V5

CONTAINS

SUBROUTINE TONGA_Initialize_HPTORA_V5 ( HPTORA )

      implicit none

!  Type structure to be filled

      Type(Create_HPTORA), intent(out) :: HPTORA

!Robmod 1/25/2024. - Add Surface height and albedo to HPTORA type structure (now V5

      HPTORA%Surface_Albedo = ZERO
      HPTORA%Surface_Height = ZERO

!  HPT, air and O3 data

      HPTORA%nlayers    = 0
      HPTORA%Heights    = ZERO
      HPTORA%Temps      = ZERO
      HPTORA%Airdensity = ZERO
      HPTORA%O3density  = ZERO
      HPTORA%O3Du       = ZERO

!  8/31/22. Include SO2 total and SO2 profile in DU
!mick mod 8/17/2023 - added "Retrieve_SO2"

      HPTORA%Include_SO2   = .FALSE.
      HPTORA%TotalSO2      =  ZERO
      HPTORA%SO2Du         =  ZERO
      HPTORA%dSO2Du_dA     =  ZERO
      HPTORA%SO2Layerflags = .FALSE.
      HPTORA%Retrieve_SO2  = .FALSE.

!  Additional for VLRRS work

      HPTORA%LayerAirColumns = ZERO
      HPTORA%LayerTemps      = ZERO

!  Cross sections and Rayleigh
!    --8/31/22. Include SO2 cross-sections

      HPTORA%nwavs          = 0
      HPTORA%SO2c0_xsecs    = ZERO
      HPTORA%SO2c1_xsecs    = ZERO
      HPTORA%SO2c2_xsecs    = ZERO
      HPTORA%o3c0_xsecs     = ZERO
      HPTORA%o3c1_xsecs     = ZERO
      HPTORA%o3c2_xsecs     = ZERO
      HPTORA%RAYLEIGH_XSECS = ZERO
      HPTORA%RAYLEIGH_DEPOL = ZERO

!  Molecular optical depths
!    --8/31/22. Separate O3 GasOD from Total

      HPTORA%Gasod   = ZERO
      HPTORA%GasodO3 = ZERO
      HPTORA%Rayod   = ZERO

!  Aerosol extinction coefficient + linearization, at reference wavelength

      HPTORA%Aer_Refw              =  ZERO
      HPTORA%Aer_Refw_Extinction   =  ZERO
      HPTORA%Aer_Refw_L_Extinction =  ZERO
      HPTORA%Retrieve_NReal        = .FALSE.
      HPTORA%Retrieve_NImag        = .FALSE.

!  Aerosol load propereties
!mick mod 8/17/2023 - added "Retrieve_PlumeHWHM"

      HPTORA%HWHMPar            =  ZERO
      HPTORA%AerUpperLimit      =  ZERO
      HPTORA%AerLowerLimit      =  ZERO
      HPTORA%aerlayerflags      = .FALSE.
      HPTORA%Loading            =  ZERO
      HPTORA%dLoading_dA        =  ZERO
      HPTORA%Retrieve_PlumeHWHM = .FALSE.

END SUBROUTINE TONGA_Initialize_HPTORA_V5

!

SUBROUTINE TONGA_Create_HPTORA_V5 &
                ( scene, nwavs, wavs, n1, n2, nfine, Surface_height, Surface_Albedo, &
                  Include_SO2, SO2_firstguess, SO2XSec_Source, retrieve_SO2,  &
                  HPTO3_filename, physics_datapath,                           &
                  plumetop, plumebot, plumehwhm_firstguess, HPTORA, fail, message )

!mick mod 8/17/2023 - replaced "SO2_totalcolumn_DU" with "SO2_firstguess"
!                   - replaced "plumehwm" with "plumehwhm_firstguess"
!                   - added "retrieve_SO2"

!Robmod 1/25/2024. - Add surface height as input, and regrid original PTHO3 to this height
!                  - Regridding is done before the aerosol plume finelayer regridding
!                  - initial surface-height regridding is always below the aerosol plume
!                  - Add Surface height and albedo to HPTORA type structure (now V5)

      implicit none

!  inputs from the initial setup

      integer  , intent(in) :: scene, nwavs, n1, n2, nfine
      real(mpk), intent(in) :: wavs (E_MAXWAVS)
      real(mpk), intent(in) :: plumetop, plumebot, plumehwhm_firstguess

!  1/25/24. Add surface height/albedo input
!mick mod 10/28/2024 - modified intent of "Surface_height" to "inout"

      real(mpk), intent(inout) :: Surface_height
      real(mpk), intent(in)    :: Surface_Albedo

!  8/31/22. Add SO2 control

      real(mpk), intent(in) :: SO2_firstguess
      logical  , intent(in) :: include_SO2, retrieve_SO2
      character*(*), intent(in) :: SO2XSec_Source, HPTO3_filename, physics_datapath

!  Type structure to be filled

      Type(Create_HPTORA),  intent(inout) :: HPTORA

!  exception handling

      logical      , intent(out) :: fail
      character*(*), intent(out) :: message

!  Local help
!mick fix 3/6/2023 - changed max dim of "sumgas1" & "sumgas2" from 220 to E_MAXWAVS

      integer   :: w, n, nc
      real(mpk) :: hdiff, sig1, sig2, temp1, temp2, temp1sq, temp2sq
      real(mpk) :: sumgas1(E_MAXWAVS), sumgas2(E_MAXWAVS)

!  parameters

      real(mpk), parameter :: CO2_PPMV_MIXRATIO = 410.0_mpk
      real(mpk), parameter :: TZERO             = 273.15_mpk
      real(mpk), parameter :: DU_TO_CM2         = 2.6867e+16_mpk

!  Filenames

      character*256 :: O3_xsecs_filename, SO2_xsecs_filename

!  START CODE
!  @@@@@@@@@@

!  File name(s)

!mick test - use Won TROPOMI-convolved O3 xsecs
      !O3_xsecs_filename  = Trim(physics_datapath)//'o3abs_brion_195_660_vacfinal.dat'
      O3_xsecs_filename  = Trim(physics_datapath)//'o3abs_brion_195_660_vacfinal.dat.conv'

!  Number of points

      HPTORA%nwavs = nwavs

!  1/25/24. Set surface height/albedo input in Type structure
!mick mod 10/28/2024 - moved defining of these HPTORA values until after call to 
!                      subroutine Get_HPTO3_TROPOMI_V5

      !HPTORA%Surface_height = Surface_height
      !HPTORA%Surface_Albedo = Surface_Albedo

!  Get the Rayleigh

      CALL RAYLEIGH_FUNCTION  &
          ( E_MAXWAVS, CO2_PPMV_MIXRATIO, NWAVS, WAVS,   &
            HPTORA%RAYLEIGH_XSECS, HPTORA%RAYLEIGH_DEPOL )

!  Get the ozone cross-sections

      CALL GET_O3XSEC_3 ( O3_xsecs_filename, E_MAXWAVS, NWAVS, WAVS, &
                          HPTORA%o3c0_xsecs, HPTORA%o3c1_xsecs, HPTORA%o3c2_xsecs, fail, message )
      if ( fail ) return

!  debug
!      do w = 1, nwavs
!         write(*,*)w,wavs(w),HPTORA%RAYLEIGH_XSECS(w), HPTORA%RAYLEIGH_DEPOL(w), &
!                          HPTORA%o3c0_xsecs(w), HPTORA%o3c1_xsecs(w), HPTORA%o3c2_xsecs(w)
!      enddo
!pause

!  8/31/22. Add SO2 cross-sections (if flagged)
!   - 2 different sources.

      if ( Include_SO2 ) then
         if ( Trim(SO2XSec_Source) .eq. 'Bogumil' ) then

!mick test - use Won TROPOMI-convolved SO2 xsecs
            !SO2_xsecs_filename = Trim(physics_datapath)//'Bogumil_so2_xsec_tmp_270_400nm.dat'
            SO2_xsecs_filename = Trim(physics_datapath)//'Bogumil_so2_xsec_tmp_270_400nm.dat.conv'

            CALL GET_SO2XSEC_3 ( SO2_xsecs_filename, E_MAXWAVS, NWAVS, WAVS, &
                                 HPTORA%SO2c0_xsecs, HPTORA%SO2c1_xsecs, HPTORA%SO2c2_xsecs, fail, message )
         else
            HPTORA%SO2c1_xsecs = 0.0d0 ;  HPTORA%SO2c1_xsecs = 0.0d0
            SO2_xsecs_filename = Trim(physics_datapath)//'SO2_298_BISA.nm'
            CALL GET_SO2XSEC_2 ( SO2_xsecs_filename, E_MAXWAVS, NWAVS, WAVS, &
                                 HPTORA%SO2c0_xsecs, fail, message )
         endif
         if ( fail ) return
      else 
         HPTORA%SO2c0_xsecs = 0.0d0 ; HPTORA%SO2c1_xsecs = 0.0d0 ; HPTORA%SO2c2_xsecs = 0.0d0
      endif

!  debug
!     do w = 1, nwavs
!         write(136,*)w,wavs(w),HPTORA%RAYLEIGH_XSECS(w), HPTORA%RAYLEIGH_DEPOL(w), HPTORA%so2c0_xsecs(w), HPTORA%o3c0_xsecs(w)
!      enddo
!      pause

!  get PTH/O3. Cannot get SO2, as this is plume-related

!Robmod 1/25/2024. - Add surface height as input, and regrid original PTHO3 to this height
!                  - Regridding is done before the aerosol plume finelayer regridding
!                  - initial surface-height regridding is always below the aerosol plume
!mick fix 10/28/2024 - added scene index "scene" to inputs

      CALL Get_HPTO3_TROPOMI_V5 &
         ( scene, HPTO3_filename, n1, n2, nfine, Surface_Height, &
           HPTORA%nlayers, HPTORA%Heights, HPTORA%Temps, HPTORA%Airdensity, HPTORA%O3density, HPTORA%O3du, fail, message )
      if ( fail ) return

!  1/25/24. Set surface height/albedo input in Type structure
!mick mod 10/28/2024 - moved defining of these HPTORA values from above

      HPTORA%Surface_height = Surface_height
      HPTORA%Surface_Albedo = Surface_Albedo

!  8/2/22. A couple of additional entries (necessary for the VLRRS code)

      do n = 1, HPTORA%nlayers
         nc = n - 1 ; hdiff = 1.0d+05 * (HPTORA%heights(nc) - HPTORA%heights(n)) * 0.5d0
         HPTORA%LayerTemps(n)      = 0.5d0 * ( HPTORA%Temps(nc) +  HPTORA%Temps(n) )
         HPTORA%LayerAirColumns(n) = hdiff * ( HPTORA%AirDensity(nc) + HPTORA%AirDensity(n) )
      enddo

!  Plume limits and half width (dummy if not retrieved)

      HPTORA%AerUpperLimit = Plumetop
      HPTORA%AerLowerLimit = PlumeBot
      HPTORA%HWHMPar       = Log(3.0_mpk + sqrt(8.0_mpk)) / plumehwhm_firstguess

!   8/31/22. Cannot include SO2 profile at this stage, as it is plume-related, but set parameter
!mick mod 8/17/2023 - added Retrieve_SO2

      HPTORA%TotalSO2    = 0.0d0 ; HPTORA%SO2Du = 0.0d0 ; HPTORA%dSO2Du_dA = 0.0d0
      HPTORA%Include_SO2 = Include_SO2
      IF ( Include_SO2 ) THEN
        HPTORA%TotalSO2      = SO2_firstguess
        HPTORA%Retrieve_SO2  = retrieve_SO2
      ENDIF

!  Apply O3 Cross-sections and determine molecular optical depths
!  (1) This section uses cruder hydrostatic equation

      do w = 1, nwavs
         do n = 1, HPTORA%nlayers
            nc = n - 1 ; hdiff = 1.0d+05 * (HPTORA%heights(nc) - HPTORA%heights(n)) * 0.5d0
            if ( n.eq.1) then
               temp1 = HPTORA%Temps(nc) - TZERO ; temp1sq = temp1 * temp1
               sig1  = HPTORA%o3c0_xsecs(w) + temp1 * HPTORA%o3c1_xsecs(w) + temp1sq * HPTORA%o3c2_xsecs(w)
               temp2 = HPTORA%Temps(n) - TZERO ; temp2sq = temp2 * temp2
               sig2  = HPTORA%o3c0_xsecs(w) + temp2 * HPTORA%o3c1_xsecs(w) + temp2sq * HPTORA%o3c2_xsecs(w)
            else
               sig1  = sig2
               temp2 = HPTORA%Temps(n) - TZERO ; temp2sq = temp2 * temp2
               sig2  = HPTORA%o3c0_xsecs(w) + temp2 * HPTORA%o3c1_xsecs(w) + temp2sq * HPTORA%o3c2_xsecs(w)
            endif
            HPTORA%GasodO3(w,n)= hdiff * ( sig1 * HPTORA%O3Density(nc) + sig2 * HPTORA%O3Density(n) ) 
            HPTORA%Rayod(w,n)  = hdiff * ( HPTORA%AirDensity(nc) + HPTORA%AirDensity(n) ) * HPTORA%Rayleigh_xsecs(w)
         enddo
         sumgas1(w) = sum(HPTORA%GasodO3(w,1:HPTORA%nlayers))
      enddo

!  (2) This section uses column o3 amounts directly with average temperature. Better.
!mick fix 8/15/2022 - in "temp1" calc, replaced "nc-1" with "n"
!mick mod 2/9/2023  - moved calc of "sumgas2" outside layer loop as in section above

      do w = 1, nwavs
         do n = 1, HPTORA%nlayers
            nc = n - 1 ; hdiff = 1.0d+05 * (HPTORA%heights(nc) - HPTORA%heights(n)) * 0.5d0
            !temp1 = 0.5d0 * ( HPTORA%Temps(nc) +  HPTORA%Temps(nc-1) ) - TZERO ; temp1sq = temp1 * temp1
            temp1 = 0.5d0 * ( HPTORA%Temps(nc) +  HPTORA%Temps(n) ) - TZERO ; temp1sq = temp1 * temp1
            sig1  = HPTORA%o3c0_xsecs(w) + temp1 * HPTORA%o3c1_xsecs(w) + temp1sq * HPTORA%o3c2_xsecs(w)
            HPTORA%GasodO3(w,n)  = du_to_cm2 * sig1 * HPTORA%O3Du(n)
            !sumgas2(w) = sum(HPTORA%GasodO3(w,1:HPTORA%nlayers))
            HPTORA%Rayod(w,n)  = hdiff * ( HPTORA%AirDensity(nc) + HPTORA%AirDensity(n) ) * HPTORA%Rayleigh_xsecs(w)
         enddo
         sumgas2(w) = sum(HPTORA%GasodO3(w,1:HPTORA%nlayers))
      enddo

!  debug optical
!mick note 2/9/20203 - note that "sumgas1" & "sumgas2" only used for diagnostics (here!)
!  do w = 1, nwavs
!     write(146,'(i3,f11.5,3f10.5)')w,wavs(w),sumgas1(w),sumgas2(w)*253.01/235.40,sum(HPTORA%Rayod(w,1:HPTORA%nlayers))
!  enddo

!mick check
!write(*,'(/4(1x,a))') 'Idx','    Wvl    ','O3 Col Tau #1','O3 Col Tau #2'
!write(*,'( 4(1x,a))') '---','-----------','-------------','-------------'
!do w = 1, nwavs
!   !write(146,'(1x,i3,2x,f9.5,1x,2(3x,f8.5,3x))') w,wavs(w),sumgas1(w),sumgas2(w)
!    write(*,  '(1x,i3,2x,f9.5,1x,2(3x,f8.5,3x))') w,wavs(w),sumgas1(w),sumgas2(w)
!enddo
!stop

!  8/31/22. Set GasOD = GasodO3, if no SO2

      if ( .not. Include_SO2 ) then
        do w = 1, nwavs
          do n = 1, HPTORA%nlayers
            HPTORA%Gasod(w,n) = HPTORA%GasodO3(w,n)
          enddo
        enddo
      endif

!  normal return

      return

END SUBROUTINE TONGA_Create_HPTORA_V5

!

SUBROUTINE TONGA_Write_HPTORA_V5 ( TID, NPARS, HPTORA )

!mick mod 8/17/2023 - added "Retrieve_SO2" and "Retrieve_PlumeHWHM"
!Rob mod 1/25/2024. - Add surface height and albedo to HPTORA, write them out here

      implicit none

!  Thread index

      INTEGER, intent(in) :: TID, NPARS

!  Type structure to be filled

      Type(Create_HPTORA), intent(in) :: HPTORA

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: LAY, WAV, PAR, NLAYERS, NWAVS
      CHARACTER (LEN=2) :: TI

!  Open output file

      OUTUNIT = 101
      WRITE(TI,'(i2.2)') TID
      OPEN (OUTUNIT,file = 'TONGA_HPTORA_'//TI//'.dbg',status = 'replace')

!Robmod 1/25/2024. surface height and albedo in HPTORA, write them out here

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'HPTORA%Surface_height = ',HPTORA%Surface_height
      WRITE(OUTUNIT,DWFR) 'HPTORA%Surface_albedo = ',HPTORA%Surface_albedo

!  HPT, air and O3 data

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI) 'HPTORA%nlayers = ',HPTORA%nlayers

      NLAYERS = HPTORA%nlayers

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%Heights(LAY)    = ',HPTORA%Heights(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%Temps(LAY)      = ',HPTORA%Temps(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%Airdensity(LAY) = ',HPTORA%Airdensity(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%O3density(LAY)  = ',HPTORA%O3density(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%O3Du(LAY)       = ',HPTORA%O3Du(LAY)
      END DO

!  SO2 data

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL) 'HPTORA%Include_SO2 = ',HPTORA%Include_SO2

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'HPTORA%TotalSO2    = ',HPTORA%TotalSO2

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%SO2Du(LAY) = ',HPTORA%SO2Du(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO PAR=1,NPARS
          WRITE(OUTUNIT,DWFR2)  'LAY = ',LAY,' PAR = ',PAR,&
            ' HPTORA%dSO2Du_dA(LAY,PAR) = ',HPTORA%dSO2Du_dA(LAY,PAR)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFL1)  'LAY = ',LAY,' HPTORA%SO2Layerflags(LAY) = ',HPTORA%SO2Layerflags(LAY)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL) 'HPTORA%Retrieve_SO2 = ',HPTORA%Retrieve_SO2

!  Additional for VLRRS work

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%LayerAirColumns(LAY) = ',HPTORA%LayerAirColumns(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%LayerTemps(LAY)      = ',HPTORA%LayerTemps(LAY)
      END DO

!  Cross sections and Rayleigh

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI) 'HPTORA%nwavs = ',HPTORA%nwavs

      NWAVS = HPTORA%nwavs

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%SO2c0_xsecs(LAY) = ',HPTORA%SO2c0_xsecs(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%SO2c1_xsecs(LAY) = ',HPTORA%SO2c1_xsecs(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%SO2c2_xsecs(LAY) = ',HPTORA%SO2c2_xsecs(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%o3c0_xsecs(LAY) = ',HPTORA%o3c0_xsecs(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%o3c1_xsecs(LAY) = ',HPTORA%o3c1_xsecs(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%o3c2_xsecs(LAY) = ',HPTORA%o3c2_xsecs(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%RAYLEIGH_XSECS(LAY) = ',HPTORA%RAYLEIGH_XSECS(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%RAYLEIGH_DEPOL(LAY) = ',HPTORA%RAYLEIGH_DEPOL(LAY)
      END DO

!  Molecular optical depths

      WRITE(OUTUNIT,*)
      DO WAV=1,NWAVS
        DO LAY=1,NLAYERS
          WRITE(OUTUNIT,DWFR2)  ' WAV = ',WAV,' LAY = ',LAY,&
            ' HPTORA%Gasod(WAV,LAY) = ',HPTORA%Gasod(WAV,LAY)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO WAV=1,NWAVS
        DO LAY=1,NLAYERS
          WRITE(OUTUNIT,DWFR2)  ' WAV = ',WAV,' LAY = ',LAY,&
            ' HPTORA%GasodO3(WAV,LAY) = ',HPTORA%GasodO3(WAV,LAY)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO WAV=1,NWAVS
        DO LAY=1,NLAYERS
          WRITE(OUTUNIT,DWFR2)  ' WAV = ',WAV,' LAY = ',LAY,&
            ' HPTORA%Rayod(WAV,LAY) = ',HPTORA%Rayod(WAV,LAY)
        END DO
      END DO

!  Aerosol extinction coefficient + linearization, at reference wavelength

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'HPTORA%Aer_Refw              = ',HPTORA%Aer_Refw

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'HPTORA%Aer_Refw_Extinction   = ',HPTORA%Aer_Refw_Extinction

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'HPTORA%Aer_Refw_L_Extinction = ',HPTORA%Aer_Refw_L_Extinction

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL) 'HPTORA%Retrieve_NReal        = ',HPTORA%Retrieve_NReal

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL) 'HPTORA%Retrieve_NImag        = ',HPTORA%Retrieve_NImag

!  Aerosol load propereties

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'HPTORA%HWHMPar       = ',HPTORA%HWHMPar

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'HPTORA%AerUpperLimit = ',HPTORA%AerUpperLimit

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'HPTORA%AerLowerLimit = ',HPTORA%AerLowerLimit

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFL1) 'LAY = ',LAY,' HPTORA%aerlayerflags(LAY) = ',HPTORA%aerlayerflags(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) 'LAY = ',LAY,' HPTORA%Loading(LAY) = ',HPTORA%Loading(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO PAR=1,NPARS
          WRITE(OUTUNIT,DWFR2)  'LAY = ',LAY,' PAR = ',PAR,&
            ' HPTORA%dLoading_dA(LAY,PAR) = ',HPTORA%dLoading_dA(LAY,PAR)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL) 'HPTORA%Retrieve_PlumeHWHM = ',HPTORA%Retrieve_PlumeHWHM

!  Close output file

      CLOSE (OUTUNIT)

END SUBROUTINE TONGA_Write_HPTORA_V5

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!                                      PRIVATE ROUTINES
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

SUBROUTINE Get_HPTO3_TROPOMI_V5 &
         ( scene, HPTO3_TROPOMI_filename, n1, n2, nfine, Surface_height, &
           nlayers, Heights, Temps, Airdensity, O3density, O3Du, fail, message )

!Rob mod 1/25/2024 - Add surface height as input, and regrid original PTHO3 to this height
!                  - Regridding is done before the aerosol plume finelayer regridding
!                  - initial surface-height regridding is always below the aerosol plume

      implicit none

!  Scene index

      integer  , intent(in) :: scene

!  Filename

      character*(*), intent(in) :: HPTO3_TROPOMI_filename

!  regridding control

      integer  , intent(in) :: n1, n2, nfine

!  1/25/24. Add surface height input
!mick mod 10/28/2024 - modified intent to "inout"

      real(mpk), intent(inout) :: Surface_height

!  Heights and temps, air and O3 density

      INTEGER  , intent(inout) :: nlayers
      REAL(MPK), intent(inout) :: Heights     (0:E_MAXLAYERS)
      REAL(MPK), intent(inout) :: Temps       (0:E_MAXLAYERS)
      REAL(MPK), intent(inout) :: Airdensity  (0:E_MAXLAYERS)
      REAL(MPK), intent(inout) :: O3density   (0:E_MAXLAYERS)
      REAL(MPK), intent(inout) :: O3Du          (E_MAXLAYERS)

!  exception handling

      logical      , intent(out) :: fail
      character*(*), intent(out) :: message

!  Local arrays

      REAL(MPK) :: Press        (0:E_MAXLAYERS)
      REAL(MPK) :: O3vmrs       (0:E_MAXLAYERS)
      REAL(MPK) :: O3Du_Calc      (E_MAXLAYERS)

      REAL(MPK) :: Hold_Press   (0:E_MAXLAYERS)
      REAL(MPK) :: Hold_heights (0:E_MAXLAYERS)
      REAL(MPK) :: Hold_Temps   (0:E_MAXLAYERS)
      REAL(MPK) :: Hold_O3vmrs  (0:E_MAXLAYERS)
      REAL(MPK) :: Hold_O3du      (E_MAXLAYERS)

!  Loschmidt's number (particles/cm3), STP values

      real(mpk)   , PARAMETER :: RHO_STANDARD = 2.68675D+19
      real(mpk)   , PARAMETER :: DU_TO_CM2    = 2.6867D+16
      real(mpk)   , PARAMETER :: PZERO = 1013.25D0
      real(mpk)   , PARAMETER :: TZERO = 273.15D0
      real(mpk)   , PARAMETER :: CONSTANT = RHO_STANDARD * TZERO / PZERO

!  help variables

      logical   :: trawl
      integer   :: n, nc, nm, ndum, j, stat
      real(mpk) :: qgrad, vgrad, hdiff, regrid_res, tgrad, const, delp

!  Initialize exception handling

      fail = .false.
      message = ' '

!  Initialize

      heights = 0.0d0 ; press = 0.0d0 ; temps = 0.0d0 ; o3vmrs = 0.0d0 ; o3du = 0.0d0

!  read the TROPOMI file
!mick fix 5/2/2023 - replaced original read: obtain "nlayers" from file

      Open(1,file=Trim(HPTO3_TROPOMI_filename), err=97, status='old' )

      !nlayers = 76
      !do n = 1, nlayers
      !   read(1,*) ndum,heights(n),press(n),temps(n),O3du(n),o3vmrs(n) ; o3vmrs(n) = 1.0d-06 * o3vmrs(n)
      !enddo

      n = 0
      do
         n = n + 1
         read(1,*,iostat=stat) ndum,heights(n),press(n),temps(n),O3du(n),o3vmrs(n) ; o3vmrs(n) = 1.0d-06 * o3vmrs(n)
         if (stat .eq. -1) then
           !At EOF: define "nlayers"
           nlayers = n - 1
           exit
         endif
      enddo
      close(1)

!  Set TOA at 80 km.

      heights(0) = 80.0d0
      qgrad      = log(press(1)/press(2)) / ( heights(1) - heights(2) )
      press(0)   = exp ( log(press(1)) + qgrad * ( heights(0) - heights(1) ) )
      temps(0)   = temps(1)
      o3vmrs(0)  = 0.0d0
!      vgrad      = ( o3vmrs(1) - o3vmrs(2) ) / ( heights(1)-heights(2) )
!      o3vmrs(0)  = o3vmrs(1) + vgrad * ( heights(0) - heights(1) )

!  debug this section
!     do n = 0, nlayers
!       write(145,'(i3,1p4e14.5)')n,o3vmrs(n),temps(n),press(n),heights(n)
!     enddo

!  First re-gridding of original PTHO3 to surface height
!  -----------------------------------------------------
!  1/25/24. New Section. Regrid done before second plume finelayer regridding (below)

!  1/25/24. Check initial surface-height regridding is always below the aerosol plume

      if ( Surface_height .ge. heights(n2) ) then
        message = 'Surface height is inside the aerosol plume - not allowed!!!'
        fail = .true. ; return
      endif

!  1/25/24. Check surface height is below the lowest height level in the data set
!mick fix 10/28/2024 - modified surface height check to focus on surface height of
!                      the plume scene

      !if ( Surface_height .lt. heights(nlayers) ) then
      !  message = 'Surface height is below lowest height level in data set - not allowed!!!'
      !  fail = .true. ; return
      !endif

      if ( (scene .eq. 2) .and. (Surface_height .lt. heights(nlayers)) ) then
        write(*,*)
        write(*,*) 'Alert: Retrieval cfg file input surface height was below lowest height level'
        write(*,*) '       in atmosphere file describing the plume scene - adjusted surface'
        write(*,*) '       height to that level.'
        Surface_height = heights(nlayers)
      endif

!  Regridding not necessary if surface is at the lowest height level in the data set
!mick mod 10/28/2024 - added ELSE section and writes to give feedback to user

      if ( heights(nlayers).eq.Surface_height ) then
        write(*,'(/1x,a)') 'Albedo of atmospheric lower boundary is at the surface.'
        go to 556
      else
        write(*,'(/1x,a)') 'Alert: albedo of atmospheric lower boundary is above the surface'
        write(*,'( 1x,a)') '       --> lower level cloud and/or aerosol present.'
      endif

!  Find the layer containing the surface height

      trawl = .true. ; nc = 0
      do while (trawl .and.nc.lt.nlayers)
         nc = nc + 1
         If (heights(nc) .le. Surface_height ) trawl = .false.
      enddo
      nlayers = nc

!  Interpolate to the surface height

      n = nlayers
      nm = n - 1 ; hdiff = 1.0_mpk / (heights(nm) - heights(n))
      tgrad = hdiff * (temps(nm)  - temps(n))
      qgrad = hdiff * Log(press(nm)/press(n))
      vgrad = hdiff * (o3vmrs(nm) - o3vmrs(n))
      delp  = press(nm) - press(n)
      temps(n)   = temps(n)           + tgrad * ( Surface_height - heights(n) ) 
      press(n)   = exp( Log(press(n)) + qgrad * ( Surface_height - heights(n) ) )
      o3vmrs(n)  = o3vmrs(n)          + vgrad * ( Surface_height - heights(n) )
      o3du(n)    = o3du(n) * ( press(nc-1) - press(nc) ) / delp
      heights(n) = Surface_height

!  debug this section

!     do n = 0, nlayers
!       write(*,'(i3,1p4e14.5)')n,o3vmrs(n),temps(n),press(n),heights(n)
!     enddo
!     stop ' first regrid debug'

!  continuation point for avoidng first regridding

556  continue

!   Second re-gridding due to presence of plume
!   -------------------------------------------

!mick check
write(*,'(/1x,a)')     'Checking Total Column O3 (in DU)'
write(*,'(1x,a,i3)')   'nlayers before Aerosol Plume regridding  = ',nlayers
write(*,'(1x,a,f6.2)') 'Total Column O3 before plume re-gridding = ',sum(o3du(1:nlayers))

!  Copy to holding arrays

      hold_heights(1:nlayers) = heights(1:nlayers)
      hold_temps  (1:nlayers) = temps  (1:nlayers)
      hold_press  (1:nlayers) = press  (1:nlayers)
      hold_o3vmrs (1:nlayers) = o3vmrs (1:nlayers)
      hold_o3du   (1:nlayers) = o3du   (1:nlayers)

!  fill down to n1

      heights(1:n1) = hold_heights(1:n1)
      temps(1:n1)   = hold_temps(1:n1)
      press(1:n1)   = hold_press(1:n1)
      o3vmrs(1:n1)  = hold_o3vmrs(1:n1)

!  fine-grid buffer zone

      nc = n1 ; regrid_res = 1.0_mpk / real(nfine,mpk)
      do n = n1 + 1, n2
         nm = n - 1 ; hdiff = 1.0_mpk / (hold_heights(nm) - hold_heights(n))
         tgrad = hdiff * (hold_temps(nm)  - hold_temps(n))
         qgrad = hdiff * Log(hold_press(nm)/hold_press(n))
         vgrad = hdiff * (hold_o3vmrs(nm) - hold_o3vmrs(n))
         delp  = hold_press(nm) - hold_press(n)
         do j = 1, nfine - 1
            nc = nc + 1
            heights(nc) = hold_heights(nm) - real(j,mpk) * regrid_res
            temps(nc)   = hold_temps(n)           + tgrad * ( heights(nc) - hold_heights(n) ) 
            press(nc)   = exp( Log(hold_press(n)) + qgrad * ( heights(nc) - hold_heights(n) ) )
            o3vmrs(nc)  = hold_o3vmrs(n)          + vgrad * ( heights(nc) - hold_heights(n) )
            o3du(nc)    = hold_o3du(n) * ( press(nc-1) - press(nc) ) / delp
         enddo
         nc = nc + 1
         heights(nc) = hold_heights(n)
         temps(nc)   = hold_temps(n) 
         press(nc)   = hold_press(n)
         o3vmrs(nc)  = hold_o3vmrs(n)
         o3du  (nc)  = hold_o3du(n) * ( press(nc-1) - press(nc) ) / delp
      enddo

!  fill up remainder

      do n = n2 + 1, nlayers
         nc = nc + 1
         heights(nc) = hold_heights(n)
         temps(nc)   = hold_temps(n) 
         press(nc)   = hold_press(n)
         o3vmrs(nc)  = hold_o3vmrs(n)
         o3du(nc)    = hold_o3du(n)
      enddo

!  Final nlayers

      nlayers = nc

!mick check
write(*,'(/1x,a)')     'Checking Total Column O3 (in DU) #2'
write(*,'(1x,a,i3)')   'nlayers after  = ',nlayers
write(*,'(1x,a,f6.2)') 'Total Column O3 after re-gridding  = ',sum(o3du(1:nlayers))

!stop

!  Check HTP

      write(*,*)
      write(*,'(1x,a,3(1x,a))') 'level','   height  ','    temp   ','   press   '
      write(*,'(1x,a,3(1x,a))') '-----','-----------','-----------','-----------'
      do nc = 0, nlayers
        write(*,'(3x,i3.3,1x,3(1x,es11.4))') nc,heights(nc),temps(nc),press(nc)
      enddo
      write(*,*)
      !stop 'at Get_HPTO3_TROPOMI_V5 stop'

!  Air and O3 density

      const = 0.5d0 * 1.0d+5 / du_to_cm2
      do n = 0, nlayers
        airdensity(n) = constant * press(n)/temps(n)
        O3density (n) = O3vmrs(n) * airdensity(n)
        if ( n.gt.0) then
           o3du_calc(n) = const * ( heights(n-1) - heights(n) ) * ( O3density (n) + O3density (n-1) ) 
!write(*,'(i4,f10.5,1p4e15.7)')n,heights(n),airdensity(n),O3density (n), o3du_calc(n) , o3du(n) 
        endif
      enddo

!mick check
!write(*,*)
!write(*,*) 'In Get_HPTO3_TROPOMI_V5'
!write(*,*) 'Checking Total Column O3 (in DU)'
!write(*,'(/2(1x,a))') ' Grid Using O3 VMR  ','   Grid Using DU    '
!write(*,'( 2(1x,a))') '--------------------','--------------------'
!write(*,'(2(8x,f6.2,6x))') sum(o3du_calc(1:nlayers)),sum(o3du(1:nlayers))
!write(*,*)
!stop

!  Normal return

      return

!  error returns

97    continue
      message = 'Create HPTORA TROPOMI_V5 - HPTO3_TROPOMI_filename file not found' ; fail = .true.
      return


END SUBROUTINE Get_HPTO3_TROPOMI_V5


SUBROUTINE GET_O3XSEC_3 ( FILENAME, MAXLAMBDAS, NLAMBDAS, LAMBDAS, &
                          o3c0_xsecs, o3c1_xsecs, o3c2_xsecs, fail, message )

      implicit none

!  Dimensioning

      INTEGER       :: MAXLAMBDAS

!  Filename

      character*(*) :: FILENAME

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(mpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Output arguments
!  ----------------

      REAL(mpk),    DIMENSION( MAXLAMBDAS ) :: o3c0_xsecs 
      REAL(mpk),    DIMENSION( MAXLAMBDAS ) :: o3c1_xsecs 
      REAL(mpk),    DIMENSION( MAXLAMBDAS ) :: o3c2_xsecs 

!  Exceptiopn handling

      logical             :: FAIL
      character*(*)       :: MESSAGE

!  Local stuff
!  -----------

   integer, parameter :: maxspline = 7500  ! This covers the range 250-315 nm

   logical       :: reading
   integer       :: nbuff, ndata, n, maxdata
   real(mpk)     :: lamb1, lamb2, wav, val, c1, c2, conversion, xsec, wav1, wav2
   real(mpk)     :: x(maxspline), y(maxspline), y2(maxspline)
   real(mpk)     :: yc1(maxspline), yc2(maxspline)
   character*5   :: c5

!  initialize exception handling

   fail    = .false.
   message = ' '

!  Initialize output

   o3c0_xsecs = 0.0d0;  o3c1_xsecs = 0.0d0;  o3c2_xsecs = 0.0d0

!  Buffer control, initialization

   nbuff = 0
   ndata = 1
   lamb1 = lambdas(1)        - 0.5d0
   lamb2 = lambdas(nlambdas) + 0.5d0
   maxdata    = 150000
   conversion = 1.0d+20

!  Read the data, check that window is covered

   open(1,file=adjustl(trim(filename)),err=90,status='old')

!  .. First line is a dummy
         read(1,*)
!  .. Read second line, set lower limit
         read(1,*)wav, val, c1, c2 ; wav1 = wav
!  .. Read more lines, saving to spline buffer, check it does not get too big!
         reading = .true.
         do while (reading .and. ndata .lt. maxdata )
            ndata = ndata + 1
            read(1,*) wav, val, c1, c2
            if ( wav .ge. lamb1 ) then
               if ( wav .le. lamb2 ) then
                  nbuff = nbuff + 1
                  if (nbuff .le. maxspline) then
                     x(nbuff)   = wav
                     y(nbuff)   = val
                     yc1(nbuff) = c1
                     yc2(nbuff) = c2
                  endif
               else if ( wav .gt. lamb2 ) then
                  reading = .false.
                  wav2 = wav
               endif
            endif
         enddo
         if (ndata .eq. maxdata) wav2 = wav
         close(1)

!  Buffer too small

        if (nbuff .gt. maxspline) then
           write(c5,'(I5)')nbuff-1
           message = 'Error: You need to increase maxspline in GET_O3XSEC_3 to at least '//c5
           fail = .true. ; return
        endif

! Spline data if they exist

         if ( nbuff.eq.0) then

            o3c0_xsecs(:) = 0.0d0
            o3c1_xsecs(:) = 0.0d0
            o3c2_xsecs(:) = 0.0d0

         else

            call spline(maxspline,x,y,nbuff,0.0d0,0.0d0,y2)
            do n = 1, nlambdas
               if ( lambdas(n).ge.wav1 .and.lambdas(n).le.wav2 ) then
                  call splint(maxspline,x,y,y2,nbuff,lambdas(n),xsec)
                  o3c0_xsecs(n) = xsec / conversion
               endif
            enddo

            call spline(maxspline,x,yc1,nbuff,0.0d0,0.0d0,y2)
            do n = 1, nlambdas
               if ( lambdas(n).ge.wav1 .and. lambdas(n).le.wav2 ) then
                  call splint(maxspline,x,yc1,y2,nbuff,lambdas(n),xsec)
                  o3c1_xsecs(n) = xsec / conversion
               endif
            enddo

            call spline(maxspline,x,yc2,nbuff,0.0d0,0.0d0,y2)
            do n = 1, nlambdas
               if ( lambdas(n).ge.wav1 .and. lambdas(n).le.wav2) then
                  call splint(maxspline,x,yc2,y2,nbuff,lambdas(n),xsec)
                  o3c2_xsecs(n) = xsec / conversion
               endif
            enddo

         endif

!  normal return

   return

!  error return

90 continue
   fail    = .true.
   message = 'Open failure for Xsec file = '//filename(1:LEN(filename))
   return

!  Finish

  return
END SUBROUTINE GET_O3XSEC_3

!

SUBROUTINE GET_SO2XSEC_2 ( FILENAME, MAXLAMBDAS, NLAMBDAS, LAMBDAS, &
                           so2_crosssecs, fail, message )

      implicit none

!  Dimensioning

      INTEGER       :: MAXLAMBDAS

!  Filename

      character*(*) :: FILENAME

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(mpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(mpk),    DIMENSION( MAXLAMBDAS ) :: SO2_CROSSSECS 

!  Exceptiopn handling

   logical             :: FAIL
   character*(*)       :: MESSAGE

!  Local stuff
!  -----------

   integer, parameter :: maxspline = 7500
   logical       :: reading
   integer       :: nbuff, ndata, n, maxdata
   real(mpk)     :: lamb1, lamb2, wav, val, conversion, xsec
   real(mpk)     :: x(maxspline), y(maxspline), y2(maxspline)

!  initialize exception handling

   fail    = .false.
   message = ' '

!  Buffer control, initialization

   nbuff = 0
   ndata = 1
   lamb1 = lambdas(1)        - 0.05d0
   lamb2 = lambdas(nlambdas) + 0.05d0
   maxdata    = 23450
   conversion = 1.0d+20

!  Read the data, check that window is covered

    open(1,file=adjustl(trim(filename)),err=90,status='old')

!  .. Read First line, and check window is covered
         read(1,*)wav, val
         if (lamb1 .lt. wav ) then
            fail = .true.
            message = 'SO2 Xsec data does not cover input window at lower end'
            return
         endif
!  .. Read more lines, saving to spline buffer
         reading = .true.
         do while (reading .and. ndata .lt. maxdata )
            ndata = ndata + 1
            read(1,*) wav, val
            if ( wav .gt. lamb1 ) then
               if ( wav .lt. lamb2 ) then
                  nbuff = nbuff + 1
                  x(nbuff)   = wav
                  y(nbuff)   = val
               endif
               if ( wav .ge. lamb2 ) reading = .false.
            endif
         enddo
!  .. Check if last data line is present, then window not covered
         if (ndata .eq. maxdata) then
            fail = .true.
            message = 'SO2 Xsec data does not cover input window at upper end'
            return
         endif
         close(1)

         if (nbuff .gt. maxspline) then
            fail = .true.
            message = 'need to increase maxspline in GET_SO2XSEC_2'
            return
         endif
         
!        write(*,*)nbuff, lamb1, lamb2
!         do n = 1, nbuff
!            write(*,*)x(n),y(n),yc1(n)
!         enddo
!         pause
         
!  Spline the data

   call spline(maxspline,x,y,nbuff,0.0d0,0.0d0,y2)
   do n = 1, nlambdas
      call splint(maxspline,x,y,y2,nbuff,lambdas(n),xsec)
      so2_crosssecs(n) = xsec
   enddo

!  normal return

   return

!  error return

90 continue
   fail    = .true.
   message = 'Open failure for Xsec file = '//filename(1:LEN(filename))
   return

!  Finish

   return
END SUBROUTINE GET_SO2XSEC_2

!

SUBROUTINE GET_SO2XSEC_3 ( FILENAME, MAXLAMBDAS, NLAMBDAS, LAMBDAS, &
                          so2c0_xsecs, so2c1_xsecs, so2c2_xsecs, fail, message )

      implicit none

!  Dimensioning

      INTEGER       :: MAXLAMBDAS

!  Filename

      character*(*) :: FILENAME

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(mpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Output arguments
!  ----------------

      REAL(mpk),    DIMENSION( MAXLAMBDAS ) :: so2c0_xsecs 
      REAL(mpk),    DIMENSION( MAXLAMBDAS ) :: so2c1_xsecs 
      REAL(mpk),    DIMENSION( MAXLAMBDAS ) :: so2c2_xsecs 

!  Exceptiopn handling

      logical             :: FAIL
      character*(*)       :: MESSAGE

!  Local stuff
!  -----------

   integer, parameter :: maxspline = 7500  ! This covers the range 270-315 nm

   logical       :: reading
   integer       :: nbuff, ndata, n, maxdata, k
   real(mpk)     :: lamb1, lamb2, wav, val, c1, c2, conversion, xsec, wav1, wav2
   real(mpk)     :: x(maxspline), y(maxspline), y2(maxspline)
   real(mpk)     :: yc1(maxspline), yc2(maxspline)
   character*5   :: c5

!  initialize exception handling

   fail    = .false.
   message = ' '

!  Initialize output

   so2c0_xsecs = 0.0d0;  so2c1_xsecs = 0.0d0;  so2c2_xsecs = 0.0d0

!  Buffer control, initialization

   nbuff = 0
   ndata = 1
   lamb1 = lambdas(1)        - 0.5d0
   lamb2 = lambdas(nlambdas) + 0.5d0
   maxdata    = 13001
   conversion = 1.0d+20

!  Read the data, check that window is covered

   open(1,file=adjustl(trim(filename)),err=90,status='old')

!  .. First 17 lines header
         do k = 1, 17
            read(1,*)
         enddo
!  .. Read first good line, set lower limit
         read(1,*)wav, val, c1, c2 ; wav1 = wav
!  .. Read more lines, saving to spline buffer, check it does not get too big!
         reading = .true.
         do while (reading .and. ndata .lt. maxdata )
            ndata = ndata + 1
            read(1,*) wav, val, c1, c2
            if ( wav .ge. lamb1 ) then
               if ( wav .le. lamb2 ) then
                  nbuff = nbuff + 1
                  if (nbuff .le. maxspline) then
                     x(nbuff)   = wav
                     y(nbuff)   = val
                     yc1(nbuff) = c1
                     yc2(nbuff) = c2
                  endif
               else if ( wav .gt. lamb2 ) then
                  reading = .false.
                  wav2 = wav
               endif
            endif
         enddo
         if (ndata .eq. maxdata) wav2 = wav
         close(1)

!  Buffer too small

        if (nbuff .gt. maxspline) then
           write(c5,'(I5)')nbuff-1
           message = 'Error: You need to increase maxspline in GET_SO2XSEC_3 to at least '//c5
           fail = .true. ; return
        endif

! Spline data if they exist

         if ( nbuff.gt.0) then

            call spline(maxspline,x,y,nbuff,0.0d0,0.0d0,y2)
            do n = 1, nlambdas
               if ( lambdas(n).ge.wav1 .and.lambdas(n).le.wav2 ) then
                  call splint(maxspline,x,y,y2,nbuff,lambdas(n),xsec)
                  so2c0_xsecs(n) = xsec / conversion
               endif
            enddo

            call spline(maxspline,x,yc1,nbuff,0.0d0,0.0d0,y2)
            do n = 1, nlambdas
               if ( lambdas(n).ge.wav1 .and. lambdas(n).le.wav2 ) then
                  call splint(maxspline,x,yc1,y2,nbuff,lambdas(n),xsec)
                  so2c1_xsecs(n) = xsec
               endif
            enddo

            call spline(maxspline,x,yc2,nbuff,0.0d0,0.0d0,y2)
            do n = 1, nlambdas
               if ( lambdas(n).ge.wav1 .and. lambdas(n).le.wav2) then
                  call splint(maxspline,x,yc2,y2,nbuff,lambdas(n),xsec)
                  so2c2_xsecs(n) = xsec
               endif
            enddo

         endif

!  normal return

   return

!  error return

90 continue
   fail    = .true.
   message = 'Open failure for Xsec file = '//filename(1:LEN(filename))
   return

!  Finish

  return
END SUBROUTINE GET_SO2XSEC_3

!

SUBROUTINE RAYLEIGH_FUNCTION                       &
          ( FORWARD_MAXLAMBDAS, CO2_PPMV_MIXRATIO, &
            FORWARD_NLAMBDAS,   FORWARD_LAMBDAS,   &
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Rayleigh cross sections and depolarization ratios
!     Bodhaine et. al. (1999) formulae
!     Module is stand-alone.
!     Wavelengths in nm (formerly Angstroms)

!  Input arguments
!  ---------------

!  wavelength
 
      INTEGER     :: FORWARD_MAXLAMBDAS, FORWARD_NLAMBDAS
      real(fpk),    dimension ( FORWARD_MAXLAMBDAS ) :: FORWARD_LAMBDAS 

!  CO2 mixing ratio

      real(fpk)    :: CO2_PPMV_MIXRATIO

!  Output arguments
!  ----------------

!  cross-sections and depolarization output

      real(fpk),    dimension ( FORWARD_MAXLAMBDAS ) :: RAYLEIGH_XSEC 
      real(fpk),    dimension ( FORWARD_MAXLAMBDAS ) :: RAYLEIGH_DEPOL

!  Local variables
!  ---------------

      INTEGER      :: W
      real(fpk)    :: MASS_DRYAIR
      real(fpk)    :: NMOL, PI, CONS
      real(fpk)    :: MO2,MN2,MARG,MCO2,MAIR
      real(fpk)    :: FO2,FN2,FARG,FCO2,FAIR
      real(fpk)    :: LAMBDA_C,LAMBDA_M,LPM2,LP2
      real(fpk)    :: N300M1,NCO2M1,NCO2
      real(fpk)    :: NCO2SQ, NSQM1,NSQP2,TERM

!  data statements and parameters
!  ------------------------------

      DATA MO2  / 20.946D0 /
      DATA MN2  / 78.084D0 /
      DATA MARG / 0.934D0 /

      real(fpk),    PARAMETER ::        S0_A = 15.0556D0
      real(fpk),    PARAMETER ::        S0_B = 28.9595D0

      real(fpk),    PARAMETER ::        S1_A = 8060.51D0
      real(fpk),    PARAMETER ::        S1_B = 2.48099D+06
      real(fpk),    PARAMETER ::        S1_C = 132.274D0
      real(fpk),    PARAMETER ::        S1_D = 1.74557D+04
      real(fpk),    PARAMETER ::        S1_E = 39.32957D0

      real(fpk),    PARAMETER ::        S2_A = 0.54D0

      real(fpk),    PARAMETER ::        S3_A = 1.034D0
      real(fpk),    PARAMETER ::        S3_B = 3.17D-04
      real(fpk),    PARAMETER ::        S3_C = 1.096D0
      real(fpk),    PARAMETER ::        S3_D = 1.385D-03
      real(fpk),    PARAMETER ::        S3_E = 1.448D-04

!  Start of code
!  -------------

!  constants

      NMOL = 2.546899D19
      PI   = DATAN(1.0D0)*4.0D0
      CONS = 24.0D0 * PI * PI * PI

!  convert co2.  mass of dry air: Eq.(17) of BWDS

!  Older code
      MCO2 = 1.0D-06 * CO2_PPMV_MIXRATIO
      MASS_DRYAIR = S0_A * MCO2 + S0_B

! VMG correction
!      xco2 = 1.0d-06 * co2_ppmv_mixratio
!      mco2 = 100.0d0 * xco2
!      MASS_DRYAIR = S0_A * MCO2 + S0_B

!  start loop

      DO W = 1, FORWARD_NLAMBDAS

!  wavelength in micrometers

      LAMBDA_M = 1.0D-03 * FORWARD_LAMBDAS(W)
      LAMBDA_C = 1.0D-07 * FORWARD_LAMBDAS(W)
!      LAMBDA_M = 1.0D-04 * FORWARD_LAMBDAS(W)             ! Angstrom input
!      LAMBDA_C = 1.0D-08 * FORWARD_LAMBDAS(W)             ! Angstrom input
      LPM2     = 1.0D0 / LAMBDA_M / LAMBDA_M

!  step 1: Eq.(18) of BWDS

      N300M1 = S1_A + ( S1_B / ( S1_C - LPM2 ) ) + ( S1_D / ( S1_E - LPM2 ) )
      N300M1 = N300M1 * 1.0D-08

!  step 2: Eq.(19) of BWDS

! Older code
      nco2m1 = n300m1 * ( 1.0d0 + s2_a * ( mco2  - 0.0003d0 ) )

!  VMG correction
!      NCO2M1 = N300M1 * ( 1.0D0 + S2_A * ( MCO2  - 0.0003D0 ) )

      NCO2   = NCO2M1 + 1.0d0
      NCO2SQ = NCO2 * NCO2

!  step 3: Eqs. (5&6) of BWDS (Bates' results)

      FN2  = S3_A + S3_B * LPM2
      FO2  = S3_C + S3_D * LPM2 + S3_E * LPM2 * LPM2

!  step 4: Eq.(23) of BWDS
!     ---> King factor and depolarization ratio

      FARG = 1.0D0
      FCO2 = 1.15D0
      MAIR = MN2 + MO2 + MARG + MCO2
      FAIR = MN2*FN2 + MO2*FO2 + MARG*FARG + MCO2*FCO2
      FAIR = FAIR / MAIR
      RAYLEIGH_DEPOL(W) = 6.0D0*(FAIR-1.0D0)/(3.0D0+7.0D0*FAIR)

!  step 5: Eq.(22) of BWDS
!     ---> Cross section

      LP2  = LAMBDA_C * LAMBDA_C
      NSQM1 = NCO2SQ - 1.0D0
      NSQP2 = NCO2SQ + 2.0D0
      TERM = NSQM1 / LP2 / NMOL / NSQP2
      RAYLEIGH_XSEC(W) =  CONS * TERM * TERM * FAIR

!  end loop

      ENDDO

!  finish
!  ------

      RETURN
END SUBROUTINE RAYLEIGH_FUNCTION

!

SUBROUTINE spline(nmax,x,y,n,yp1,ypn,y2)
      INTEGER            ::  n,nmax
      REAL(kind=mpk)     :: yp1,ypn,x(nmax),y(nmax),y2(nmax)
      INTEGER            :: i,k
      REAL(kind=mpk)     :: p,qn,sig,un,u(NMAX)

      if (yp1.gt..99e30) then
        y2(1)=0.0d0
        u(1)=0.0d0
      else
        y2(1)=-0.5d0
        u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.0d0*( (y(i+1)-y(i)) / (x(i+1)-x(i)) - &
                     (y(i)-y(i-1)) / (x(i)-x(i-1))   &
                   ) / (x(i+1)-x(i-1)) - sig*u(i-1)  &
            )/p
      enddo

      if (ypn.gt..99d30) then
        qn=0.0d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo

      return
 END SUBROUTINE SPLINE


 SUBROUTINE splint(nmax,xa,ya,y2a,n,x,y)
      INTEGER            ::  n, nmax
      REAL(kind=mpk)     :: x, y, xa(nmax),ya(nmax),y2a(nmax)
      INTEGER            ::  k,khi,klo
      REAL(kind=mpk)     :: a,b,h

      klo=1
      khi=n

 1    if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
            khi=k
          else
            klo=k
          endif
       goto 1
       endif
       h=xa(khi)-xa(klo)
!       if (h.eq.0.0d0) pause
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+ &
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0

      return
 END SUBROUTINE SPLINT

!  End module

END MODULE TONGA_Create_HPTORA_V5_m

