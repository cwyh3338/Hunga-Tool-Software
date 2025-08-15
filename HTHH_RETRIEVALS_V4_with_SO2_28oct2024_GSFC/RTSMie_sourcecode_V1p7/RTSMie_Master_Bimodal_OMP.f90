module RTSMie_master_bimodal_OMP_m

!  This is the Bimodal master #1 for RTS MIE code
!    ** RT Solutions, Version 1.4, 30 June     2011 (Bimodal control Tmatrix)
!    ** RT Solutions, Version 1.5, 17 August   2011 (Bimodal control Mie)
!  mick mod 12/9/2022 - modified RTSMie_master_bimodal_OMP to use
!                       (1) array version of "Mie_Inputs_Mode"
!                       (2) OpenMP

  use RTSMie_parameters_m
  use RTSMie_Inputs_Def_m
  use RTSMie_Outputs_Def_m
  use RTSMie_Master_m

!  Everything PUBLIC here
!  ----------------------

public

contains

subroutine RTSMie_master_bimodal_OMP &
       ( Mie_Inputs_Mode, Bimodal_fraction, & ! I
         Mie_Outputs_BiModal,               & ! O
         fail, istatus, Bmessages, trace_3 )  ! O

   USE RTSMie_parameters_m

!  implicit none statement

   IMPLICIT NONE

!  I/O Type structures

   TYPE(RTSMie_Inputs) , INTENT (inout)    :: Mie_Inputs_Mode(2)
   TYPE(RTSMie_Outputs), INTENT (inout)    :: Mie_Outputs_BiModal

!  Bimodal Fraction (by number weight, First PSD mode)

   REAL(KIND=DP), INTENT(IN)  :: Bimodal_fraction

!  Exception handling

   LOGICAL          , INTENT (OUT)   :: fail
   INTEGER          , INTENT (OUT)   :: istatus
   CHARACTER(LEN=*) , INTENT (OUT)   :: Bmessages(3), trace_3

!  Local Variables
!  ---------------

!  Individual mode Mie outputs

   TYPE(RTSMie_Outputs)    :: Mie_Outputs_Mode(2)

!  Proxies

   integer       :: n_Fmatrix_Angles
   integer       :: Mie1_ncoeffs, Mie2_ncoeffs
   Logical       :: do_Expcoeffs
   Logical       :: do_FMatrix

!  Help

   integer       :: L, n, npsd
   real(KIND=dp) :: FF1, FF2, WW1, WW2
   real(KIND=dp) :: Csca1_FF1, Csca2_FF2, Csca_total
   character(len=1) :: chpsd

!  OMP Variables
!  -------------
!mick mod 12/9/2022 - added this section

!  OpenMP (general)

   INTEGER   :: OMP_MAXTHREADS
   INTEGER   :: TID, OMP_NTHREADS

!  OpenMP functions

   INTEGER   :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

!  Zero the output
!  ---------------

!  Main output

   Mie_Outputs_BiModal%Mie_dist1  = d_zero
   Mie_Outputs_BiModal%Mie_dist2  = d_zero
   Mie_Outputs_BiModal%Mie_bulk   = d_zero
   Mie_Outputs_BiModal%Mie_asymm  = d_zero

   Mie_Outputs_BiModal%Mie_Fmatrix   = d_zero
   Mie_Outputs_BiModal%Mie_ncoeffs   = 0
   Mie_Outputs_BiModal%Mie_expcoeffs = d_zero

!  Exception handling

   fail           = .false.
   istatus        = 0
   Bmessages(1:3) = ' '
   trace_3        = ' '

!  Checks and proxy settings
!  =========================

!  fraction proxy

   FF1 = Bimodal_Fraction
   FF2 = d_one - FF1

!  Check: No Monodisperse here !

   if ( Mie_Inputs_Mode(1)%do_monodisperse .or. Mie_Inputs_Mode(2)%do_monodisperse ) then
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error - MONODISPERSE FLAG must be Off!'
      return
   endif

!  Lambda check

   if ( Mie_Inputs_Mode(1)%lambda .ne. Mie_Inputs_Mode(2)%lambda ) then
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error - Wavelengths must be equal!'
      return
   endif

!  Coefficients check
!mick fix 12/9/2022 - added ELSEIF section

   if ( Mie_Inputs_Mode(1)%do_Expcoeffs .and. Mie_Inputs_Mode(2)%do_Expcoeffs ) then
      do_ExpCoeffs = .true.
   elseif ( .not.Mie_Inputs_Mode(1)%do_Expcoeffs .and. .not.Mie_Inputs_Mode(2)%do_Expcoeffs ) then
      do_ExpCoeffs = .false.
   else
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error - do_ExpCoeffs flags must be the same !'
      return
   endif

!  F-Matrix check
!mick fix 12/9/2022 - added ELSEIF section

   if ( Mie_Inputs_Mode(1)%do_Fmatrix .and. Mie_Inputs_Mode(2)%do_Fmatrix ) then
      do_Fmatrix = .true.
      if ( Mie_Inputs_Mode(1)%n_Fmatrix_Angles .ne. Mie_Inputs_Mode(2)%n_Fmatrix_Angles ) then
         fail = .true.; istatus = 1
         trace_3 = 'RTSMie_master_bimodal module: Input error - F_Matrix angles must be the same !'
         return
      else
         n_Fmatrix_Angles = Mie_Inputs_Mode(1)%n_Fmatrix_Angles
      endif
   elseif ( .not.Mie_Inputs_Mode(1)%do_Fmatrix .and. .not.Mie_Inputs_Mode(2)%do_Fmatrix ) then
      do_Fmatrix = .false.
   else
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error - do_Fmatrix flags must be the same !'
      return
   endif

!  RTS Mie calls
!  =============
!mick mod 12/9/2022 - added OpenMP parallel region

   OMP_MAXTHREADS = 2
   CALL OMP_SET_NUM_THREADS(OMP_MAXTHREADS)

!  Begin parallel region

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (Mie_Inputs_Mode, Mie_Outputs_Mode)  ! Mie I/O

!  Obtain thread number

      tid = OMP_GET_THREAD_NUM()
      !tid = 0

!  Obtain and display total number of threads and local thread

      if (tid == 0) then
         omp_nthreads = OMP_GET_NUM_THREADS()
         !omp_nthreads = 1
         !write(*,*)
         !write(*,'(1x,a,i1)') 'total number of threads (inside parallel region) = ', omp_nthreads
         write(*,'(1x,a,i1,a)') 'Running RTSMie_master_bimodal_OMP with ',omp_nthreads,' thread(s)'
      end if

!$OMP DO

      do npsd=1,2
         write(chpsd,'(i1)') npsd
         write(*,'(1x,a,i1,a)') '** Thread ',tid,' doing Mie for PSD mode ' // chpsd // ' ----------------------'
         CALL RTSMie_Master  & 
            ( Mie_Inputs_Mode(npsd), Mie_Outputs_Mode(npsd),              & ! I/O type
              fail, istatus, Bmessages(1), Bmessages(2), Bmessages(3) )   ! Exception handling

!  Exception handling

         if ( fail ) then
            write(*,'(/1x,a)') 'RTSMie_master_bimodal_OMP: error during RTSMie_Master call for PSD ' // chpsd
            write(*,'( 1x,a)') 'Error message(s) are:'
            do n = 1,3
              write(*,'(8x,a,i1,a)') 'Message # ',n,' : ' // trim(adjustl(Bmessages(n)))
            enddo
            if ( Istatus .eq. 1 ) stop
         endif
      enddo

!$OMP END DO

!  End parallel region

!$OMP END PARALLEL


!  Bimodal determination
!  =====================

!  Bulk properties
!  ---------------

!  Revision 20 September 2011
!    Correct definition for Expcoeffs/Fmatrix: WW1/WW2 in place of FF1/FF2
!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong

   Csca1_FF1  = FF1 * Mie_Outputs_Mode(1)%Mie_bulk(2)   
   Csca2_FF2  = FF2 * Mie_Outputs_Mode(2)%Mie_bulk(2)   
   Csca_total =  Csca1_FF1 + Csca2_FF2

   WW1 = Csca1_FF1 / Csca_total
   WW2 = Csca2_FF2 / Csca_total

   Mie_Outputs_Bimodal%Mie_bulk(1:2) = FF1 * Mie_Outputs_Mode(1)%Mie_bulk(1:2) + &
                                       FF2 * Mie_Outputs_Mode(2)%Mie_bulk(1:2)
   Mie_Outputs_Bimodal%Mie_bulk(3)   = Mie_Outputs_Bimodal%Mie_bulk(2) / Mie_Outputs_Bimodal%Mie_bulk(1)

!  Coefficients
!  ------------

   if ( Do_Expcoeffs ) then
      Mie_Outputs_Bimodal%Mie_asymm = WW1 * Mie_Outputs_Mode(1)%Mie_asymm + &
                                      WW2 * Mie_Outputs_Mode(2)%Mie_asymm
      Mie1_ncoeffs = Mie_Outputs_Mode(1)%Mie_ncoeffs
      Mie2_ncoeffs = Mie_Outputs_Mode(2)%Mie_ncoeffs
      Mie_Outputs_Bimodal%Mie_ncoeffs = max(Mie1_ncoeffs,Mie2_ncoeffs)

!mick check
!write(*,*)
!write(*,*) 'Mie_Outputs_Bimodal%Mie_ncoeffs = ',Mie_Outputs_Bimodal%Mie_ncoeffs
!stop 'at Bimodal ncoeffs check'

      do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
         Mie_Outputs_Bimodal%Mie_expcoeffs(1:6,L) = WW1 * Mie_Outputs_Mode(1)%Mie_expcoeffs(1:6,L) + & 
                                                    WW2 * Mie_Outputs_Mode(2)%Mie_expcoeffs(1:6,L)
      enddo
      if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
          do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
              Mie_Outputs_Bimodal%Mie_expcoeffs(1:6,L) = WW2 * Mie_Outputs_Mode(2)%Mie_expcoeffs(1:6,L)
          enddo
      else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
          do L = Mie2_ncoeffs + 1, Mie1_ncoeffs
              Mie_Outputs_Bimodal%Mie_expcoeffs(1:6,L) = WW1 * Mie_Outputs_Mode(1)%Mie_expcoeffs(1:6,L)
          enddo
      endif
   endif

!  Fmatrix
!  -------

   if ( Do_Fmatrix ) then
      do L = 1, n_Fmatrix_angles
         Mie_Outputs_Bimodal%Mie_Fmatrix(1:4,L) = WW1 * Mie_Outputs_Mode(1)%Mie_Fmatrix(1:4,L) + & 
                                                  WW2 * Mie_Outputs_Mode(2)%Mie_Fmatrix(1:4,L)
      enddo
   endif   

!  Distributions
!  -------------

    Mie_Outputs_Bimodal%Mie_Dist1 = Mie_Outputs_Mode(1)%Mie_Dist1
    Mie_Outputs_Bimodal%Mie_Dist2 = Mie_Outputs_Mode(2)%Mie_Dist1

!  Finish

   return
end subroutine RTSMie_Master_Bimodal_OMP

!  End module

end module RTSMie_Master_Bimodal_OMP_m

