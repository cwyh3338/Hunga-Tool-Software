
module TONGA_Mie_Initialize_m

!  Mie modules

      USE RTSMie_parameters_m
      USE RTSMie_Inputs_Def_m

private
public  :: TONGA_Mie_Initialize

contains

subroutine TONGA_Mie_Initialize &
      ( retrieve_nreal, retrieve_nimag, AOD_refw, & !Input
        PSDpar1_firstguess, PSDpar2_firstguess,   & !Input
        nreal_firstguess, nimag_firstguess,       & !Input
        Mie_Inputs_Global )                         !Output

!mick note 12/9/2022:
!  This subroutine used to set up Mie calculations for a Bimodal PSD; however,
!  it is still set up to only RETRIEVE nreal & nimag for one PSD (in this case
!  PSD #1).
!mick mod 8/17/2023 - removed input "npars"
!                   - added inputs "retrieve_nreal" and "retrieve_nimag"

   implicit none

   integer, parameter :: mpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   logical,   intent(in) :: retrieve_nreal, retrieve_nimag
   real(mpk), intent(in) :: AOD_refw
   real(mpk), intent(in) :: PSDpar1_firstguess(2), PSDpar2_firstguess(2)
   real(mpk), intent(in) :: nreal_firstguess(2), nimag_firstguess(2)

!  Outputs
!  -------

!  Mie Global input structure

   Type(RTSMie_Inputs), intent(out) :: Mie_Inputs_Global(2)

!  Initialize Mie code
!  ===================

!  Inputs for PSD #1 (Fine Mode)
!  -----------------------------

!  Mie General flags

   Mie_Inputs_Global(1)%do_Expcoeffs    = .false.      ! Will be set true for the main Mie calculations
   Mie_Inputs_Global(1)%do_Fmatrix      = .false.      ! Not needed
   Mie_Inputs_Global(1)%do_monodisperse = .false.      ! Not needed
   Mie_Inputs_Global(1)%MonoRadius      =  0.0d0       ! Not needed
   Mie_Inputs_Global(1)%do_LinearREF    = (retrieve_nreal .or. retrieve_nimag) 
   Mie_Inputs_Global(1)%do_LinearPSD    = .false.      ! Not needed (So Far?!?)

!  Mie PSD is Lognormal with reff and veff given

   Mie_Inputs_Global(1)%PSD_index   = 4                     ! 2-Par Lognormal
   Mie_Inputs_Global(1)%PSD_pars(1) = psdpar1_firstguess(1) ! R_g
   Mie_Inputs_Global(1)%PSD_pars(2) = psdpar2_firstguess(1) ! S_g
   Mie_Inputs_Global(1)%PSD_pars(3) = 0.0d0                 ! Not needed

!  Wavelength and refractive index, initial global settings. These will change locally.

   Mie_Inputs_Global(1)%lambda      = AOD_refw * 0.001d0  ! Microns. Set later on, according to UV wavelengths
   Mie_Inputs_Global(1)%n_real      = nreal_firstguess(1)
   Mie_Inputs_Global(1)%n_imag      = nimag_firstguess(1)

!  PSD integration. Fairly standard set of conditions

   Mie_Inputs_Global(1)%nblocks         = 40
   Mie_Inputs_Global(1)%nweights        = 32
   Mie_Inputs_Global(1)%xparticle_limit = 10000.0d0
   Mie_Inputs_Global(1)%R1R2_cutoff     = 1.0d-06
   Mie_Inputs_Global(1)%FixR1R2         = .false.    ! False uses the next 2 inputs
   Mie_Inputs_Global(1)%R1              = 0.001d0    ! starting  value PSD radius (Microns)
   Mie_Inputs_Global(1)%R2              = 10.0d0     ! finishing value PSD radius (Microns)

!  Angular outputs (Not needed)

   Mie_Inputs_Global(1)%n_Fmatrix_angles = 0
   Mie_Inputs_Global(1)%Fmatrix_angles   = 0.0d0

!  Inputs for PSD #2 (Coarse Mode)
!  -----------------------------

!mick note:

!  Mie General flags

   Mie_Inputs_Global(2)%do_Expcoeffs    = .false.   ! Will be set true for the main Mie calculations
   Mie_Inputs_Global(2)%do_Fmatrix      = .false.   ! Not needed
   Mie_Inputs_Global(2)%do_monodisperse = .false.   ! Not needed
   Mie_Inputs_Global(2)%MonoRadius      =  0.0d0    ! Not needed
   Mie_Inputs_Global(2)%do_LinearREF    = .false.   ! Not retrieving ref index of PSD #2
   Mie_Inputs_Global(2)%do_LinearPSD    = .false.   ! Not needed (So Far?!?)

!  Mie PSD is Lognormal with reff and veff given

   Mie_Inputs_Global(2)%PSD_index   = 4                     ! 2-Par Lognormal
   Mie_Inputs_Global(2)%PSD_pars(1) = psdpar1_firstguess(2) ! R_g
   Mie_Inputs_Global(2)%PSD_pars(2) = psdpar2_firstguess(2) ! S_g
   Mie_Inputs_Global(2)%PSD_pars(3) = 0.0d0                 ! Not needed

!  Wavelength and refractive index, initial global settings. These will change locally.

   Mie_Inputs_Global(2)%lambda      = AOD_refw * 0.001d0  ! Microns. Set later on, according to UV wavelengths
   Mie_Inputs_Global(2)%n_real      = nreal_firstguess(2)
   Mie_Inputs_Global(2)%n_imag      = nimag_firstguess(2)

!  PSD integration. Fairly standard set of conditions
!mick note 12/9/20222 - some of these coarse mode inputs may need adjusted

   Mie_Inputs_Global(2)%nblocks         = 40
   Mie_Inputs_Global(2)%nweights        = 32
   Mie_Inputs_Global(2)%xparticle_limit = 10000.0d0
   Mie_Inputs_Global(2)%R1R2_cutoff     = 1.0d-06
   Mie_Inputs_Global(2)%FixR1R2         = .false.    ! False uses the next 2 inputs
   Mie_Inputs_Global(2)%R1              = 0.001d0    ! starting  value PSD radius (Microns)
   Mie_Inputs_Global(2)%R2              = 10.0d0     ! finishing value PSD radius (Microns)

!  Angular outputs (Not needed)

   Mie_Inputs_Global(2)%n_Fmatrix_angles = 0
   Mie_Inputs_Global(2)%Fmatrix_angles   = 0.0d0

end subroutine TONGA_Mie_Initialize

!  End module

end module TONGA_Mie_Initialize_m