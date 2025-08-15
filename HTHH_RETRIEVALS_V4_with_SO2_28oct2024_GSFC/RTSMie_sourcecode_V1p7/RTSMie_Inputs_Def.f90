Module RTSMie_Inputs_Def_m

   use RTSMie_parameters_m, only : dp, max_Mie_angles

   implicit none

! #####################################################################
! #####################################################################

   TYPE RTSMie_Inputs

!  List of Inputs
!  ==============

!  Flag inputs
!  -----------

!  Do_Expcoeffs      - Boolean flag for computing Expansion Coefficients
!  Do_Fmatrix        - Boolean flag for computing F-matrix at equal-angles

   LOGICAL :: Do_Expcoeffs
   LOGICAL :: Do_Fmatrix

!  Do_Monodisperse   - Boolean flag for Doing a Monodisperse calculation.
!                      If set, the PSD stuff will be turned off internally

   LOGICAL :: Do_Monodisperse

!  Linearization control

!  Do_LinearRef      - Boolean Flag for doing Refractive Index linearization
!  Do_LinearPSD      - Boolean Flag for doing PSD linearization.
!                        This is checked and turned off for Monodisperse

   LOGICAL :: Do_LinearRef
   LOGICAL :: Do_LinearPSD

!  PSD inputs
!  ----------

!  FixR1R2 : If set, Use "R1R2_cutoff" for smallest particle size
!                    then Internal routine to calculate R1 and R2 (outputs)
!            If Not set, Use Input R1 and R2 for PSD limits.

   LOGICAL :: FixR1R2

!  R1, R2  : Minimum and Maximum radii (Microns)

   REAL    (KIND=dp)  :: R1, R2

!  Limiting particle size value. Set to 10000.0 default.
!    If you exceed this, program will tell you to increase dimensioning.

   REAL    (KIND=dp) :: xparticle_limit

!      Monoradius     - Monodisperse radius size (Microns)

   REAL    (KIND=dp)  :: Monoradius

!      psd_Index      - Index for particle size distribution of spheres
!      psd_pars       - Parameters characterizing PSD (up to 3 allowed)

!  Mie inputs (distribution index, PSD parameters)
!    PSD_index = 1 : TWO-PARAMETER GAMMA with alpha and b given
!    PSD_index = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given
!    PSD_index = 3 : BIMODAL GAMMA with equal mode weights
!    PSD_index = 4 : LOG-NORMAL with rg and sigma given
!    PSD_index = 5 : LOG-NORMAL with reff and veff given
!    PSD_index = 6 : POWER LAW
!    PSD_index = 7 : MODIFIED GAMMA with alpha, rc and gamma given
!    PSD_index = 8 : MODIFIED GAMMA with alpha, b and gamma given

   INTEGER            :: PSD_Index
   REAL    (KIND=dp)  :: PSD_pars (3)

!  PSD quadrature control
!  ----------------------

!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.
!    R1R2_cutoff particle size for setting R1 and R2 internally

   INTEGER           :: nblocks
   INTEGER           :: nweights
   REAL    (KIND=dp) :: R1R2_cutoff

!  Optical: Wavelength, refractive index
!  -------------------------------------

!      LAMBDA         - wavelength of light (microns)
!      N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)

   REAL    (KIND=dp)  :: lambda, n_real, n_imag

!  F-matrix Angular control input
!  ------------------------------

!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER           :: n_Fmatrix_angles
   REAL    (KIND=dp) :: Fmatrix_angles(max_Mie_angles)

   END TYPE RTSMie_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: RTSMie_Inputs

end module RTSMie_Inputs_def_m
