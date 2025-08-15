Module Tonga_Aerosol_Loading_m

  USE TONGA_pars_m
  USE TONGA_HPTORA_V5_m
  USE Loading_General_Routines_m

!  Preliminary Version  7/1/21 ( based on GEMSTOOL and VLIDORT Version 2.8.1. )
!  First Version with VLIDORT 2.8.3, 7/30/21.

private :: Aer_Loading_Assign
public  :: Tonga_Aerosol_Loading

contains

SUBROUTINE Tonga_Aerosol_Loading ( MAXA, MA, APARS, LISTA, HPTORA, fail, message )

!  double precision

   implicit none
   integer, parameter :: mpk = SELECTED_REAL_KIND(15)

!  Arguments
!mick mod 8/17/2023 - added LISTA

   INTEGER  , intent(in) :: MAXA, MA
   REAL(mpk), intent(in) :: APARS(MAXA)
   INTEGER  , intent(in) :: LISTA(MAXA)

   Type(Create_HPTORA) , intent(inout) :: HPTORA

!  Exception handling (already initialized)
!mick mod 9/22/2022 - modified message passing to allow all 3 loading fail
!                     messages to get passed

   logical      , intent(inout) :: fail
   character*(*), intent(inout) :: message(3)

!  Local
!  -----

!  Local Assignments, aerosol loading control
!  TOTAL AOT = aerosol optical thickness (value at reference wavelength )
!    Case 1 = uniform aerosol 
!    Case 2 = Exponential Loading
!    Case 3 = GDF loading

   logical  , parameter :: do_Jacobians = .true.
   integer  , parameter :: acase        = 3
   real(mpk), parameter :: dumzero = 0.0_mpk

   integer    :: i
   real(mpk)  :: aertau_input_w0
   real(mpk)  :: gdfloading_peakheight, gdfloading_halfwidth
   real(mpk)  :: Const
   real(mpk)  :: loading_upperboundary
   real(mpk)  :: loading_lowerboundary

!  Loading output and derivatives, flags

    real(mpk),    dimension ( E_maxlayers )    :: Loading
    logical,      dimension ( E_maxlayers )    :: AERLAYERFLAGS 
    real(mpk),    dimension ( E_maxlayers )    :: DLoading_Dtau
    real(mpk),    dimension ( E_maxlayers, 2 ) :: Dloading_Dpars

!  Exception handling

   logical           :: fail1
   character*100     :: message_Loading(3)

!  debug
!   integer           :: n

!  Start Code
!  ==========

!  Constant for conversion of input HWHM to tgdfloading_halfwidth 

   CONST = log(3.0_mpk + sqrt(8.0_mpk))

!  The retrieval parameters are:
!mick note 9/22/2022 - retrieval state vector elements applied here!
!mick mod 8/17/2023  - adjusted GDF par index for HWHM for new location in state vector APARS

   aertau_input_w0       = APARS(LISTA(1))
   gdfloading_peakheight = APARS(LISTA(2))

   if ( HPTORA%Retrieve_PlumeHWHM ) then
      gdfloading_halfwidth = CONST / APARS(LISTA(4))
   else
      gdfloading_halfwidth = HPTORA%HWHMPar  ! Already converted, assumed known in this case
   endif

!  Boundaries from HPTORA

   loading_upperboundary = HPTORA%AerUpperLimit
   loading_lowerboundary = HPTORA%AerLowerLimit

!  Alternative idea

!   Window = 10.0d0 * CONST / HPTORA%HWHMPar ; Window = 5.0d0
!   loading_upperboundary = gdfloading_peakheight + Window
!   loading_lowerboundary = gdfloading_peakheight - Window 

!  Make the loading assignment call. Subroutine lifted from NGST

   CALL Aer_Loading_Assign &
      ( E_maxlayers, HPTORA%nlayers, HPTORA%heights, acase, do_Jacobians, & ! INPUT, Grid and case
        loading_upperboundary, loading_lowerboundary, aertau_input_w0,    & ! INPUT, Loading control
        dumzero, gdfloading_peakheight, gdfloading_halfwidth,             & ! INPUT, Loading control
        aerlayerflags, Loading, Dloading_Dtau, Dloading_Dpars,            & ! OUTPUT, aerosol Loading
        fail1, Message_Loading )                                            ! Exception handling

!  Error handling
!mick mod 9/22/2022 - modified message passing to allow all 3 loading fail
!                     messages to get passed

   if ( FAIL1 ) then
      fail = .true.
      do i=1,3
         message(i) = Trim(message_Loading(i))
      enddo
   endif

!  Interpret output. Derivatives from above subroutine are unnormalized, need to be normalized for VLIDORT inputs
!   Loading routine uses variable x = C / h, where H = APARS(LISTA(4)). Thus h.dL/dh = h.dL/dx.dx/dh = - (C/h).dL/dx
!   *****--> C = Log(3.0_mpk + sqrt(8.0_mpk))
!mick mod 8/17/2023 - replaced "MA .eq. 3" IF condition with "Retrieve_PlumeHWHM" condition
!                   - replaced "APARS(3)" with "APARS(LISTA(4))" to reflect new location in state vector

   HPTORA%aerlayerflags    = aerlayerflags
   HPTORA%Loading          = Loading
   HPTORA%dLoading_dA(:,1) = Dloading_Dtau(:)    * APARS(LISTA(1)) ! Aero profile depends on total AOD
   HPTORA%dLoading_dA(:,2) = Dloading_DPars(:,1) * APARS(LISTA(2)) ! Aero profile depends on Peak height

!  If you are retrieving the half width, then we have additional dependency

   if ( HPTORA%Retrieve_PlumeHWHM ) then
      HPTORA%dLoading_dA(:,3) = - Dloading_DPars(:,2) * CONST / APARS(LISTA(4)) ! Aero profile depends on HWHM
   else
      HPTORA%dLoading_dA(:,3) = 0.0d0 
   endif

!  Debug and plot to 76
!   write(*,*)MA,apars(1:3)
!   do n = 1, HPTORA%nlayers
!      write(76,'(i3, f8.3, 1p4e18.8)')n, HPTORA%heights(n), HPTORA%Loading(n), HPTORA%dLoading_dA(n,1:3)
!   enddo
!   stop ' to check Jacobians'

!  done

   RETURN
END SUBROUTINE Tonga_Aerosol_Loading

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        P R I V A T E    R O U T I N E 
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine Aer_Loading_Assign &
      ( maxlayers, nlayers, height_grid, Loading_case, do_Aer_Jacobians,    & ! INPUT, Grid and case
        loading_upperboundary, loading_lowerboundary, aertau_input_w0,      & ! INPUT, Loading control
        exploading_relaxation, gdfloading_peakheight, gdfloading_halfwidth, & ! INPUT, Loading control
        aerlayerflags, Loading, Dloading_Dtau, Dloading_Dpars,              & ! OUTPUT, aerosol Loading
        fail1, Message_Loading )                                              ! Exception handling

      implicit none

!  Precision

      integer, parameter :: mpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ======

!  Grid Numbers

   integer  , INTENT (IN) :: maxlayers, nlayers
   real(mpk), INTENT (IN) :: height_grid(0:maxlayers)

!  Aerosol loading control
!  TOTAL AOT = aerosol optical thickness (value at reference wavelength )
!    Case 1 = uniform aerosol 
!    Case 2 = Exponential Loading
!    Case 3 = GDF loading

   logical  , INTENT (IN) :: do_Aer_Jacobians
   integer  , INTENT (IN) :: loading_case
   real(mpk), INTENT (IN) :: aertau_input_w0

   real(mpk), INTENT (IN) :: exploading_relaxation
   real(mpk), INTENT (IN) :: gdfloading_peakheight
   real(mpk), INTENT (IN) :: gdfloading_halfwidth
   real(mpk), INTENT (IN) :: loading_upperboundary
   real(mpk), INTENT (IN) :: loading_lowerboundary

!  output
!  ======

!  Loading output

    real(mpk),    dimension ( maxlayers ),  INTENT (OUT)      :: Loading

!  aerosol layer flags

    LOGICAL,      DIMENSION ( maxlayers ), intent(out)  :: AERLAYERFLAGS 

!  Loading derivatives

    real(mpk),    dimension ( maxlayers )   , intent(out)    :: DLoading_Dtau
    real(mpk),    dimension ( maxlayers, 2 ), intent(out)    :: Dloading_Dpars

!  Exception handling

   logical,        INTENT (OUT)           :: fail1
   character*(*),  INTENT (OUT)           :: message_Loading(3)

!  Local

   integer :: n

!  Code
!  ====

!  Initialize exception handling

   fail1 = .false.
   Message_Loading   = ' '

!  initialize Aerosol Loading

   Loading        = 0.0d0
   Dloading_Dtau  = 0.0d0
   Dloading_Dpars = 0.0d0
   aerlayerflags = .false.

!  Now form the aerosol Loading
!  ----------------------------

!  @@@ Notes: 18 February 2013
!        profiles_lidar routine added (loading case 4)
!         - Read LIDAR Extinction [km-1] and height [km] profiles from FILE
!         - Parcel the entire LIDAR profile into the output array
!         - Ignores z_upperlimit, z_lowerlimit
!         - Only one offline test so far (19 february, 2013)

!    Loading = optical depth profile

!    Derivatives : 
!       Cases 1-3 : dloading_Dtau      = derivative of profile w.r.t the total aerosol optical depth at wavelength w0
!       Case 2    : dloading_Dpars(1)  = derivative of profile w.r.t the relaxation parameter (Exponential)
!       Case 3    : dloading_Dpars(1)  = derivative of profile w.r.t the GDF Peak height
!       Case 3    : dloading_Dpars(2)  = derivative of profile w.r.t the GDF Half Width 
!       Case 4    : dloading_Dpars     = 0 (No Derivatives allowed, LIDAR profile)

!  Case 1: Uniform layer of aerosol
!  ********************************

   if ( loading_case .eq. 1 ) then

      CALL profiles_uniform                                                 &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,              & ! Inputs
            loading_upperboundary, loading_lowerboundary, aertau_input_w0,  & ! Inputs
            Loading, Dloading_Dtau,                                         & ! Output
            fail1, message_Loading(1), message_Loading(2) )                   ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'Uniform aerosol Loading failed'

!  Case 2: Exponential decay profile
!  *********************************

   else if ( loading_case .eq. 2 ) then

      CALL profiles_expone                                                  &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,              & ! Inputs
            loading_upperboundary, loading_lowerboundary,                   & ! Inputs
            exploading_relaxation, aertau_input_w0,                         & ! Inputs
            Loading, Dloading_Dtau, Dloading_Dpars(:,1),                    & ! Output
            fail1, message_Loading(1), message_Loading(2) )                   ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'Exponential aerosol Loading failed'

!  Case 3: GDF (quasi-Gaussian) profile
!  ************************************

   else if ( loading_case .eq. 3 ) then

      CALL profiles_gdfone                                                        &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,                    & ! Inputs
            loading_upperboundary, gdfloading_peakheight, loading_lowerboundary,  & ! Inputs
            gdfloading_halfwidth, aertau_input_w0,                                & ! Inputs
            Loading, Dloading_Dtau, Dloading_Dpars(:,1), Dloading_Dpars(:,2),     & ! Output
            fail1, message_Loading(1), message_Loading(2) )                         ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'GDF aerosol Loading failed'

   endif

!  Flags

   do n = 1, nlayers
      aerlayerflags(n) = (Loading(n).gt.1.0d-15)
   enddo

!  Check loading

   !do n = 1, nlayers
   !   write(*,*) 'n = ',n,' Aero Loading(n) = ', Loading(n)
   !enddo
   !stop 'at aerosol loading stop'

!  Return if failure at this stage

   Return
end subroutine Aer_Loading_Assign

!  End module
 
End Module Tonga_Aerosol_Loading_m

