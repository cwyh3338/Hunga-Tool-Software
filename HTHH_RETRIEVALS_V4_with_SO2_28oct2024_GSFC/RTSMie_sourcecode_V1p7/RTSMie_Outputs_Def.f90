Module RTSMie_Outputs_Def_m

   use RTSMie_parameters_m, only : dp, max_Mie_angles

   implicit none

! #####################################################################
! #####################################################################

   TYPE RTSMie_Outputs

!  Output arguments
!  ================

!  Bulk distribution parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(KIND=dp) :: Mie_bulk (3)

!  Expansion coefficients and Asymmetry parameter, optional output

   integer       :: Mie_ncoeffs
   real(KIND=dp) :: Mie_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp) :: Mie_asymm

!  F-matrix, optional output

   real(KIND=dp) :: Mie_Fmatrix(4,max_Mie_angles)

!  Distribution parameters. dist1 is the main output. Dist2 is required for the Bimodal.

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(KIND=dp) :: Mie_dist1 (5)
   real(KIND=dp) :: Mie_dist2 (5)

!  Debug

   real(KIND=dp) :: Mie_DistPlot (2,max_Mie_angles)

   END TYPE RTSMie_Outputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: RTSMie_Outputs

end module RTSMie_Outputs_def_m
