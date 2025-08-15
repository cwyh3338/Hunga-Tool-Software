Module Tonga_SO2_Loading_m

  USE TONGA_pars_m
  USE TONGA_HPTORA_V5_m
  USE Loading_General_Routines_m

!  8/31/22. New

private :: SO2_Loading_Assign
public  :: Tonga_SO2_Loading

contains

SUBROUTINE Tonga_SO2_Loading ( MAXA, MA, APARS, LISTA, HPTORA, fail, message )

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
!mick mod 8/17/2023 - modified message passing to allow all 3 loading fail
!                     messages to get passed

   logical      , intent(inout) :: fail
   character*(*), intent(inout) :: message(3)

!  Local
!  -----

!  Local Assignments, SO2 loading control
!    Case 1 = uniform SO2 
!    Case 2 = Exponential Loading
!    Case 3 = GDF loading

   logical  , parameter :: do_Jacobians = .true.
   integer  , parameter :: acase        = 3
   real(mpk), parameter :: dumzero = 0.0_mpk

   integer    :: i
   real(mpk)  :: TotalSO2
   real(mpk)  :: gdfloading_peakheight, gdfloading_halfwidth
   real(mpk)  :: Const
   real(mpk)  :: loading_upperboundary
   real(mpk)  :: loading_lowerboundary

!  Loading output and derivatives, flags

    real(mpk),    dimension ( E_maxlayers )       :: Loading
    LOGICAL,      DIMENSION ( E_maxlayers )       :: Layerflags 
    real(mpk),    dimension ( E_maxlayers )       :: Dloading_Dtot
    real(mpk),    dimension ( E_maxlayers, 2 )    :: Dloading_Dpars

!  Exception handling

   logical           :: fail1
   character*100     :: message_Loading(3)

!  debug
!   integer           :: n

!  Start Code
!  ==========

!  Constant for conversion of input HWHM to tgdfloading_halfwidth 

   CONST = Log(3.0_mpk + sqrt(8.0_mpk))

!  The retrieval parameters are
!mick note 8/17/2023 - retrieval state vector elements applied here!
!mick mod 8/17/2023  - modified to retrieve SO2
!                    - adjusted GDF par index for HWHM for new location in state vector APARS

   gdfloading_peakheight = APARS(LISTA(2))

   if ( HPTORA%Retrieve_SO2 ) then
     TotalSO2 = APARS(LISTA(3))
   else
     TotalSO2 = HPTORA%TotalSO2
   endif

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

   CALL SO2_Loading_Assign &
      ( E_maxlayers, HPTORA%nlayers, HPTORA%heights, acase, do_Jacobians, & ! INPUT, Grid and case
        loading_upperboundary, loading_lowerboundary, TotalSO2,           & ! INPUT, Loading control
        dumzero, gdfloading_peakheight, gdfloading_halfwidth,             & ! INPUT, Loading control
        Layerflags, Loading, Dloading_Dtot, Dloading_Dpars,               & ! OUTPUT, SO2 Loading
        fail1, Message_Loading )                                            ! Exception handling

!  Error handling
!mick mod 8/17/2023 - modified message passing to allow all 3 loading fail
!                     messages to get passed

   if ( FAIL1 ) then
      fail = .true.
      do i=1,3
         message(i) = Trim(message_Loading(i))
      enddo
   endif

!  Interpret output. Derivatives from above subroutine are unnormalized, need to be normalized for VLIDORT inputs
!   Loading routine uses variable x = C / h, where h = APARS(LISTA(4)). Thus h.dL/dh = h.dL/dx.dx/dh = - (C/h).dL/dx
!   *****--> C = Log(3.0_mpk + sqrt(8.0_mpk))
!mick mod 8/17/20203 - activated HPTORA%dSO2Du_dA(:,1) wrt total SO2

   HPTORA%SO2layerflags  = Layerflags
   HPTORA%SO2Du          = Loading
   HPTORA%dSO2Du_dA(:,1) = Dloading_Dtot(:)    * APARS(LISTA(3)) ! SO2 profile depends on total SO2 (in DU)
   HPTORA%dSO2Du_dA(:,2) = Dloading_DPars(:,1) * APARS(LISTA(2)) ! SO2 profile depends on Peak height

!  If you are retrieving the half width, then we have additional dependency
!mick mod 8/17/20203 - activated this section

   if ( HPTORA%Retrieve_PlumeHWHM ) then
      HPTORA%dSO2Du_dA(:,3) = - Dloading_DPars(:,2) * CONST / APARS(LISTA(4)) ! SO2 profile depends on HWHM
   else
      HPTORA%dSO2Du_dA(:,3) = 0.0d0
   endif

!  done

   RETURN
END SUBROUTINE Tonga_SO2_Loading

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        P R I V A T E    R O U T I N E 
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine SO2_Loading_Assign &
      ( maxlayers, nlayers, height_grid, Loading_case, do_SO2_Jacobians,    & ! INPUT, Grid and case
        loading_upperboundary, loading_lowerboundary, TotalSO2,             & ! INPUT, Loading control
        exploading_relaxation, gdfloading_peakheight, gdfloading_halfwidth, & ! INPUT, Loading control
        Layerflags, Loading, Dloading_Dtot, Dloading_Dpars,                 & ! OUTPUT, SO2 Loading
        fail1, Message_Loading )                                              ! Exception handling

      implicit none

!  Precision

      integer, parameter :: mpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ======

!  Grid Numbers

   integer  , INTENT (IN) :: maxlayers, nlayers
   REAL(mpk), INTENT (IN) :: height_grid(0:maxlayers)

!  SO2 loading control
!  TOTAL AOT = SO2 optical thickness (value at reference wavelength )
!    Case 1 = uniform SO2 
!    Case 2 = Exponential Loading
!    Case 3 = GDF loading

   Logical  , INTENT(IN)  :: do_SO2_Jacobians
   integer  , INTENT(IN)  :: loading_case
   real(mpk), INTENT (IN) :: TotalSO2

   real(mpk), INTENT (IN) :: exploading_relaxation
   real(mpk), INTENT (IN) :: gdfloading_peakheight
   real(mpk), INTENT (IN) :: gdfloading_halfwidth
   real(mpk), INTENT (IN) :: loading_upperboundary
   real(mpk), INTENT (IN) :: loading_lowerboundary

!  output
!  ======

!  Loading output

    real(mpk),    dimension ( maxlayers ), INTENT (OUT) :: Loading

!  SO2 layer flags

    LOGICAL,      DIMENSION ( maxlayers ), intent(out)  :: LAYERFLAGS 

!  Loading derivatives

    real(mpk),    dimension ( maxlayers )   , intent(out)    :: dloading_Dtot
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

!  initialize SO2 Loading

   Loading        = 0.0d0
   dloading_Dtot  = 0.0d0
   Dloading_Dpars = 0.0d0
   Layerflags     = .false.

!  Now form the SO2 Loading
!  ------------------------

!    Loading = Dobson Unit profile
!    Derivatives : 
!       Cases 1-3 : dloading_Dtot      = derivative of profile w.r.t the total SO2 Amount
!       Case 2    : dloading_Dpars(1)  = derivative of profile w.r.t the relaxation parameter (Exponential)
!       Case 3    : dloading_Dpars(1)  = derivative of profile w.r.t the GDF Peak height
!       Case 3    : dloading_Dpars(2)  = derivative of profile w.r.t the GDF Half Width 

!  Case 1: Uniform layer of SO2
!  ********************************

   if ( loading_case .eq. 1 ) then

      CALL profiles_uniform                                                 &
          ( maxlayers, nlayers, height_grid, do_SO2_Jacobians,              & ! Inputs
            loading_upperboundary, loading_lowerboundary, TotalSO2,         & ! Inputs
            Loading, Dloading_Dtot,                                         & ! output
            fail1, message_Loading(1), message_Loading(2) )                   ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'Uniform SO2 Loading failed'

!  Case 2: Exponential decay profile
!  *********************************

   else if ( loading_case .eq. 2 ) then

      CALL profiles_expone                                                  &
          ( maxlayers, nlayers, height_grid, do_SO2_Jacobians,              & ! Inputs
            loading_upperboundary, loading_lowerboundary,                   & ! Inputs
            exploading_relaxation, TotalSO2,                                & ! Inputs
            Loading, Dloading_Dtot, Dloading_Dpars(:,1),                    & ! output
            fail1, message_Loading(1), message_Loading(2) )                   ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'Exponential SO2 Loading failed'

!  Case 3: GDF (quasi-Gaussian) profile
!  ************************************

   else if ( loading_case .eq. 3 ) then

      CALL profiles_gdfone &
          ( maxlayers, nlayers, height_grid, do_SO2_Jacobians,                    & ! Inputs
            loading_upperboundary, gdfloading_peakheight, loading_lowerboundary,  & ! Inputs
            gdfloading_halfwidth, TotalSO2,                                       & ! Inputs
            Loading, Dloading_Dtot, Dloading_Dpars(:,1), Dloading_Dpars(:,2),     & ! output
            fail1, message_Loading(1), message_Loading(2) )                         ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'GDF SO2 Loading failed'

   endif

!  Flags

   do n = 1, nlayers
      Layerflags(n) = (Loading(n).gt.1.0d-15)
   enddo
      
!  Check loading

   !do n = 1, nlayers
   !   write(*,*) 'n = ',n,' Aero Loading(n) = ', Loading(n)
   !enddo
   !stop 'at SO2 loading stop'

!  Return if failure at this stage

   Return
end subroutine SO2_Loading_Assign

!  End module
 
End Module Tonga_SO2_Loading_m

