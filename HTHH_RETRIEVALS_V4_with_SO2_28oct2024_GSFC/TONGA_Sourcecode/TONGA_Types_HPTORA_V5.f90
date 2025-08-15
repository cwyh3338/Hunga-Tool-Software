Module TONGA_HPTORA_V5_m

!   Basic property setup for the TONGA retrieval
!mick mod 9/22/2022 - added additional TONGA_pars variable control

!Robmod 1/25/2024. - Add surface height as input, and regrid original PTHO3 to this height
!                  - Regridding is done before the aerosol plume finelayer regridding
!                  - initial surface-height regridding is always below the aerosol plume
!                  - Add Surface height and albedo to HPTORA type structure (now V5)

      USE TONGA_pars_m, Only : MPK, E_MAXLAYERS, E_MAXPARS, E_MAXWAVS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE Create_HPTORA

!  1/25/24. Add Surface height and albedo to HPTORA type structure (now V5)

      REAL(MPK) :: Surface_height
      REAL(MPK) :: Surface_albedo

!  HPT, air and O3 data.

      INTEGER   :: nlayers
      REAL(MPK) :: Heights     (0:E_MAXLAYERS)
      REAL(MPK) :: Temps       (0:E_MAXLAYERS)
      REAL(MPK) :: Airdensity  (0:E_MAXLAYERS)
      REAL(MPK) :: O3density   (0:E_MAXLAYERS)
      REAL(MPK) :: O3Du          (E_MAXLAYERS)

!  8/31/22. Include SO2 total and SO2 profile in DU
!mick mod 8/17/2023 - added Retrieve_SO2

      LOGICAL   :: Include_SO2
      REAL(MPK) :: TotalSO2
      REAL(MPK) :: SO2Du         (E_MAXLAYERS)
      REAL(MPK) :: dSO2Du_dA     (E_MAXLAYERS, E_MAXPARS)
      LOGICAL   :: SO2Layerflags (E_MAXLAYERS)
      LOGICAL   :: Retrieve_SO2

!  Additional for VLRRS work

      REAL(MPK) :: LayerAirColumns (E_MAXLAYERS)
      REAL(MPK) :: LayerTemps      (E_MAXLAYERS)

!  Cross sections and Rayleigh
!    --8/31/22. Include SO2 cross-sections

      INTEGER   :: nwavs
      REAL(MPK) :: SO2c0_xsecs    (E_MAXWAVS)
      REAL(MPK) :: SO2c1_xsecs    (E_MAXWAVS)
      REAL(MPK) :: SO2c2_xsecs    (E_MAXWAVS)
      REAL(MPK) :: o3c0_xsecs     (E_MAXWAVS)
      REAL(MPK) :: o3c1_xsecs     (E_MAXWAVS)
      REAL(MPK) :: o3c2_xsecs     (E_MAXWAVS)
      REAL(MPK) :: RAYLEIGH_XSECS (E_MAXWAVS)
      REAL(MPK) :: RAYLEIGH_DEPOL (E_MAXWAVS)

!  Molecular optical depths
!    --8/31/22. Separate O3 GasOD from Total

      REAL(MPK) :: Gasod   (E_MAXWAVS,E_MAXLAYERS)
      REAL(MPK) :: GasodO3 (E_MAXWAVS,E_MAXLAYERS)
      REAL(MPK) :: Rayod   (E_MAXWAVS,E_MAXLAYERS)

!  Aerosol extinction coefficient + linearization, at reference wavelength

      REAL(MPK) :: Aer_Refw
      REAL(MPK) :: Aer_Refw_Extinction
      REAL(MPK) :: Aer_Refw_L_Extinction
      LOGICAL   :: Retrieve_NReal
      LOGICAL   :: Retrieve_NImag

!  Aerosol load propereties
!mick mod 8/17/2023 - added Retrieve_PlumeHWHM

      REAL(MPK) :: HWHMPar
      REAL(MPK) :: AerUpperLimit
      REAL(MPK) :: AerLowerLimit
      LOGICAL   :: aerlayerflags(E_MAXLAYERS)
      REAL(MPK) :: Loading      (E_MAXLAYERS)
      REAL(MPK) :: dLoading_dA  (E_MAXLAYERS, E_MAXPARS)
      LOGICAL   :: Retrieve_PlumeHWHM

      END TYPE Create_HPTORA

! #####################################################################
! #####################################################################

! Public

      PRIVATE
      PUBLIC :: Create_HPTORA

End Module TONGA_HPTORA_V5_m
