Module RTSMie_Lin_Outputs_Def_m

   use RTSMie_parameters_m, only : dp, max_Mie_angles

   implicit none

! #####################################################################
! #####################################################################

   TYPE RTSMie_Lin_Outputs

!  Output arguments
!  ================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

!  linearizations w.r.t. PSD parameters

   real(KIND=dp) :: LPSD_Mie_bulk (3,3)

!  linearizations w.r.t. Refractive index parameters

   real(KIND=dp) :: LRFE_Mie_bulk (3,2)

!  with respect to Fraction (Bimodal only)

   real(KIND=dp) :: LFRC_Mie_bulk (3)

!  Expansion coefficients and Asymmetry parameter, optional output
!  ---------------------------------------------------------------

!  linearizations w.r.t. PSD parameters

   real(KIND=dp) :: LPSD_Mie_expcoeffs (6,0:max_Mie_angles,3)
   real(KIND=dp) :: LPSD_Mie_asymm(3)

!  linearizations w.r.t. Refractive index parameters

   real(KIND=dp) :: LRFE_Mie_expcoeffs (6,0:max_Mie_angles,2)
   real(KIND=dp) :: LRFE_Mie_asymm(2)

!  with respect to Fraction (Bimodal only)

   real(KIND=dp) :: LFRC_Mie_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp) :: LFRC_Mie_asymm

!  F-matrix, optional output
!  -------------------------

!  Linearizations of F-matrix

   real(KIND=dp) :: LPSD_Mie_Fmatrix(4,max_Mie_angles,3)
   real(KIND=dp) :: LRFE_Mie_Fmatrix(4,max_Mie_angles,2)
   real(KIND=dp) :: LFRC_Mie_Fmatrix(4,max_Mie_angles)

!  Distribution parameters.
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(KIND=dp) :: LPSD_Mie_dist (5,3)

   END TYPE RTSMie_Lin_Outputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: RTSMie_Lin_Outputs

end module RTSMie_Lin_Outputs_def_m
