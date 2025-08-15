module RTSMie_Master_bimodal_plus_m

!  This is the Bimodal master #1 for RTS MIE code
!    ** RT Solutions, Version 1.4, 30 June     2011 (Bimodal control Tmatrix)
!    ** RT Solutions, Version 1.5, 17 August   2011 (Bimodal control Mie)

  use RTSMie_parameters_m
  use RTSMie_Inputs_Def_m
  use RTSMie_Outputs_Def_m
  use RTSMie_Lin_Outputs_Def_m
  use RTSMie_Master_plus_m

!  Everything PUBLIC here
!  ----------------------

public

contains

subroutine RTSMie_master_bimodal_plus &
       ( Mie_Inputs_Mode1, Mie_Inputs_Mode2, Bimodal_fraction,                  & ! I
         Mie_Outputs_BiModal, Mie_Lin_Outputs_BiModal_Mode1, Mie_Lin_Outputs_BiModal_Mode2, &
         fail, istatus, Bmessages, trace_3 )         ! O

   USE RTSMie_parameters_m

!  implicit none statement

   IMPLICIT NONE

!  I/O Type structures

   TYPE(RTSMie_Inputs) , INTENT (inout)    :: Mie_Inputs_Mode1
   TYPE(RTSMie_Inputs) , INTENT (inout)    :: Mie_Inputs_Mode2
   TYPE(RTSMie_Outputs), INTENT (inout)    :: Mie_Outputs_BiModal
   TYPE(RTSMie_Lin_Outputs), INTENT (inout)    :: Mie_Lin_Outputs_BiModal_Mode1
   TYPE(RTSMie_Lin_Outputs), INTENT (inout)    :: Mie_Lin_Outputs_BiModal_Mode2

!  Bimodal Fraction (by number weight, First PSD mode)

   real    (KIND=dp), intent(in)  :: Bimodal_fraction

!  Exception handling

   LOGICAL          , INTENT (OUT)   :: fail
   INTEGER          , INTENT (OUT)   :: istatus
   CHARACTER*(*)    , INTENT (OUT)   :: Bmessages(3), trace_3

!  Local Variables
!  ---------------

!  Individual mode Mie outputs

   TYPE(RTSMie_Outputs)     :: Mie_Outputs_Mode1
   TYPE(RTSMie_Outputs)     :: Mie_Outputs_Mode2
   TYPE(RTSMie_Lin_Outputs) :: Mie_Lin_Outputs_Mode1
   TYPE(RTSMie_Lin_Outputs) :: Mie_Lin_Outputs_Mode2

!  Proxies

   integer       :: n_Fmatrix_Angles
   integer       :: Mie1_ncoeffs, Mie2_ncoeffs
   Logical       :: do_Expcoeffs
   Logical       :: do_FMatrix
   logical       :: Do_LinearRef
   logical       :: Do_LinearPSD



!  Local Arrays
!  ============

   real(KIND=dp) :: Dasy, DCoeff(6), DFmat(4), L_Omega
   real(KIND=dp) :: Bulk(3), LBulk1(3,3), LBulk2(3,3), LBulkf(3)

!  Other local variables

   integer       :: k, L, q, nlin
   real(KIND=dp) :: FF1, FF2, WW1, WW2, TERM1, TERM2
   real(KIND=dp) :: Csca1_FF1, Csca2_FF2, Csca_total, D1_Csca1(3), D2_Csca2(3)
   real(KIND=dp) :: D1_WW1(3), D2_WW1(3), D1_WW2(3), D2_WW2(3), DF_WW1

!  Zero the output
!  ---------------

   Mie_Outputs_BiModal%Mie_dist1  = d_zero
   Mie_Outputs_BiModal%Mie_dist2  = d_zero
   Mie_Outputs_BiModal%Mie_bulk   = d_zero
   Mie_Outputs_BiModal%Mie_asymm  = d_zero

   Mie_Outputs_BiModal%Mie_Fmatrix    = d_zero
   Mie_Outputs_BiModal%Mie_ncoeffs    = 0
   Mie_Outputs_BiModal%Mie_expcoeffs  = d_zero

   Mie_Lin_Outputs_BiModal_Mode1%LPSD_Mie_bulk      = d_zero
   Mie_Lin_Outputs_BiModal_Mode1%LPSD_Mie_Fmatrix   = d_zero
   Mie_Lin_Outputs_BiModal_Mode1%LPSD_Mie_expcoeffs = d_zero
   Mie_Lin_Outputs_BiModal_Mode1%LPSD_Mie_asymm     = d_zero

   Mie_Lin_Outputs_BiModal_Mode1%LRFE_Mie_bulk      = d_zero
   Mie_Lin_Outputs_BiModal_Mode1%LRFE_Mie_Fmatrix   = d_zero
   Mie_Lin_Outputs_BiModal_Mode1%LRFE_Mie_expcoeffs = d_zero
   Mie_Lin_Outputs_BiModal_Mode1%LRFE_Mie_asymm     = d_zero

   Mie_Lin_Outputs_BiModal_Mode1%LFRC_Mie_bulk      = d_zero
   Mie_Lin_Outputs_BiModal_Mode1%LFRC_Mie_Fmatrix   = d_zero
   Mie_Lin_Outputs_BiModal_Mode1%LFRC_Mie_expcoeffs = d_zero
   Mie_Lin_Outputs_BiModal_Mode1%LFRC_Mie_asymm     = d_zero

   Mie_Lin_Outputs_BiModal_Mode2%LPSD_Mie_bulk      = d_zero
   Mie_Lin_Outputs_BiModal_Mode2%LPSD_Mie_Fmatrix   = d_zero
   Mie_Lin_Outputs_BiModal_Mode2%LPSD_Mie_expcoeffs = d_zero
   Mie_Lin_Outputs_BiModal_Mode2%LPSD_Mie_asymm     = d_zero

   Mie_Lin_Outputs_BiModal_Mode2%LRFE_Mie_bulk      = d_zero
   Mie_Lin_Outputs_BiModal_Mode2%LRFE_Mie_Fmatrix   = d_zero
   Mie_Lin_Outputs_BiModal_Mode2%LRFE_Mie_expcoeffs = d_zero
   Mie_Lin_Outputs_BiModal_Mode2%LRFE_Mie_asymm     = d_zero

   Mie_Lin_Outputs_BiModal_Mode2%LFRC_Mie_bulk      = d_zero
   Mie_Lin_Outputs_BiModal_Mode2%LFRC_Mie_Fmatrix   = d_zero
   Mie_Lin_Outputs_BiModal_Mode2%LFRC_Mie_expcoeffs = d_zero
   Mie_Lin_Outputs_BiModal_Mode2%LFRC_Mie_asymm     = d_zero

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

   if ( Mie_Inputs_Mode1%do_monodisperse .or. Mie_Inputs_Mode2%do_monodisperse ) then
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error: MONODISPERSE FLAG must be Off!'
      return
   endif

!  Lambda check

   if ( Mie_Inputs_Mode1%lambda .ne. Mie_Inputs_Mode2%lambda ) then
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error: Wavelengths must be equal!'
      return
   endif

!  Coefficients check
!mick fix 12/9/2022 - added ELSEIF section

   if ( Mie_Inputs_Mode1%do_Expcoeffs .and. Mie_Inputs_Mode2%do_Expcoeffs ) then
      do_ExpCoeffs = .true.
   elseif ( .not.Mie_Inputs_Mode1%do_Expcoeffs .and. .not.Mie_Inputs_Mode2%do_Expcoeffs ) then
      do_ExpCoeffs = .false.
   else
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error: do_ExpCoeffs flags must be the same !'
      return
   endif

!  F-Matrix check
!mick fix 12/9/2022 - added ELSEIF section

   if ( Mie_Inputs_Mode1%do_Fmatrix .and. Mie_Inputs_Mode2%do_Fmatrix ) then
      do_Fmatrix = .true.
      if ( Mie_Inputs_Mode1%n_Fmatrix_Angles .ne. Mie_Inputs_Mode2%n_Fmatrix_Angles ) then
         fail = .true.; istatus = 1
         trace_3 = 'RTSMie_master_bimodal module: Input error: F_Matrix angles must be the same !'
         return
      else
         n_Fmatrix_Angles = Mie_Inputs_Mode1%n_Fmatrix_Angles
      endif
   elseif ( .not.Mie_Inputs_Mode1%do_Fmatrix .and. .not.Mie_Inputs_Mode2%do_Fmatrix ) then
      do_Fmatrix = .false.
   else
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error: do_Fmatrix flags must be the same !'
      return
   endif

!  Linearization control
!mick fix 12/9/2022 - added ELSEIF section to both IF blocks

   if ( Mie_Inputs_Mode1%do_LinearRef .and. Mie_Inputs_Mode2%do_LinearRef ) then
      do_LinearRef = .true.
   elseif ( .not.Mie_Inputs_Mode1%do_LinearRef .and. .not.Mie_Inputs_Mode2%do_LinearRef ) then
      do_LinearRef = .false.
   else
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error: REF linearization must be the same !'
      return
   endif

   if ( Mie_Inputs_Mode1%do_LinearPSD .and. Mie_Inputs_Mode2%do_LinearPSD ) then
      do_LinearPSD = .true.
   elseif ( .not.Mie_Inputs_Mode1%do_LinearPSD .and. .not.Mie_Inputs_Mode2%do_LinearPSD ) then
      do_LinearPSD = .false.
   else
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error: PSD linearization must be the same !'
      return
   endif

!  Calls 
!  =====

!  First Call

   write(*,*)' ** Doing Mie Linearized for PSD # 1 ----------------------'

   CALL RTSMie_Master_plus & 
       ( Mie_Inputs_Mode1, Mie_Outputs_Mode1, Mie_Lin_Outputs_Mode1, & ! I/O type structures
         fail, istatus, Bmessages(1), Bmessages(2), Bmessages(3) ) ! O

!  Exception handling

   if ( fail ) then
      trace_3 = 'RTSMie_master_bimodal_PLUS module: First PSD call, Error'
      if ( Istatus .eq. 1 ) return
   endif

!  Second call
!  -----------

   write(*,*)' ** Doing Mie Linearized for PSD # 2 ----------------------'
   CALL RTSMie_Master_plus & 
       ( Mie_Inputs_Mode2, Mie_Outputs_Mode2, Mie_Lin_Outputs_Mode2, & ! I/O type structures
         fail, istatus, Bmessages(1), Bmessages(2), Bmessages(3) ) ! O

!  Exception handling

   if ( fail ) then
      trace_3 = 'RTSMie_master_bimodal_PLUS module: second PSD call, Error'
      if ( Istatus .eq. 1 ) return
   endif

!  Bimodal determination
!  =====================


!  Bulk properties
!  ---------------

!  Revision 20 September 2011
!    Correct definition for Expcoeffs/Fmatrix: WW1/WW2 in place of FF1/FF2
!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong

   Csca1_FF1  = FF1 * Mie_Outputs_Mode1%Mie_bulk(2)   
   Csca2_FF2  = FF2 * Mie_Outputs_Mode2%Mie_bulk(2)   
   Csca_total =  Csca1_FF1 + Csca2_FF2
   WW1   = Csca1_FF1 / Csca_total
   WW2   = Csca2_FF2 / Csca_total

   Mie_Outputs_Bimodal%Mie_bulk(1:2) = FF1 * Mie_Outputs_Mode1%Mie_bulk(1:2) + FF2 * Mie_Outputs_Mode2%Mie_bulk(1:2)
   Mie_Outputs_Bimodal%Mie_bulk(3)   = Mie_Outputs_Bimodal%Mie_bulk(2) / Mie_Outputs_Bimodal%Mie_bulk(1)

!  Coefficients
!  ------------

   if ( Do_Expcoeffs ) then
      Mie_Outputs_Bimodal%Mie_asymm = WW1 * Mie_Outputs_Mode1%Mie_asymm + WW2 * Mie_Outputs_Mode2%Mie_asymm
      Mie1_ncoeffs = Mie_Outputs_Mode1%Mie_ncoeffs
      Mie2_ncoeffs = Mie_Outputs_Mode2%Mie_ncoeffs
      Mie_Outputs_Bimodal%Mie_ncoeffs = max(Mie1_ncoeffs,Mie2_ncoeffs)
      do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
         Mie_Outputs_Bimodal%Mie_expcoeffs(1:6,L) = WW1 * Mie_Outputs_Mode1%Mie_expcoeffs(1:6,L) + & 
                                                    WW2 * Mie_Outputs_Mode2%Mie_expcoeffs(1:6,L)
      enddo
      if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
          do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
              Mie_Outputs_Bimodal%Mie_expcoeffs(1:6,L) = WW2 * Mie_Outputs_Mode2%Mie_expcoeffs(1:6,L)
          enddo
      else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
          do L = Mie2_ncoeffs + 1, Mie1_ncoeffs
              Mie_Outputs_Bimodal%Mie_expcoeffs(1:6,L) = WW1 * Mie_Outputs_Mode1%Mie_expcoeffs(1:6,L)
          enddo
      endif
   endif

!  Fmatrix
!  -------

   if ( Do_Fmatrix ) then
      do L = 1, n_Fmatrix_angles
         Mie_Outputs_Bimodal%Mie_Fmatrix(1:4,L) = WW1 * Mie_Outputs_Mode1%Mie_Fmatrix(1:4,L) + & 
                                                  WW2 * Mie_Outputs_Mode2%Mie_Fmatrix(1:4,L)
      enddo
   endif   

!  Distributions
!  -------------

    Mie_Outputs_Bimodal%Mie_Dist1 = Mie_Outputs_Mode1%Mie_Dist1
    Mie_Outputs_Bimodal%Mie_Dist2 = Mie_Outputs_Mode2%Mie_Dist1

!  LRFE Linearizations
!  ===================

   if ( do_LinearRef ) then
      nlin = 2

!  help variables
!  --------------

      D1_Csca1(1:nlin) = Mie_Lin_Outputs_Mode1%LRFE_Mie_bulk(2,1:nlin)
      D2_Csca2(1:nlin) = Mie_Lin_Outputs_Mode2%LRFE_Mie_bulk(2,1:nlin)
      D1_WW1(1:nlin) = FF1 * D1_Csca1(1:nlin) * WW2 / Csca_total
      D2_WW2(1:nlin) = FF2 * D2_Csca2(1:nlin) * WW1 / Csca_total
      D1_WW2(1:nlin) = - D1_WW1(1:nlin)
      D2_WW1(1:nlin) = - D2_WW2(1:nlin)

!  Bulk variables
!  --------------

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong

      Mie_Lin_Outputs_Bimodal_Mode1%LRFE_Mie_bulk(1:2,1:nlin) = FF1 * Mie_Lin_Outputs_Mode1%LRFE_Mie_bulk(1:2,1:nlin)
      Mie_Lin_Outputs_Bimodal_Mode2%LRFE_Mie_bulk(1:2,1:nlin) = FF2 * Mie_Lin_Outputs_Mode2%LRFE_Mie_bulk(1:2,1:nlin)
      Bulk               = Mie_Outputs_Bimodal%Mie_bulk
      LBulk1(1:2,1:nlin) = Mie_Lin_Outputs_Mode1%LRFE_Mie_bulk(1:2,1:nlin)
      LBulk2(1:2,1:nlin) = Mie_Lin_Outputs_Mode2%LRFE_Mie_bulk(1:2,1:nlin)
      Mie_Lin_Outputs_Bimodal_Mode1%LRFE_Mie_bulk(3,1:nlin) = ( LBulk1(2,1:nlin) - bulk(3) * LBulk1(1,1:nlin) ) / bulk(1)
      Mie_Lin_Outputs_Bimodal_Mode2%LRFE_Mie_bulk(3,1:nlin) = ( LBulk2(2,1:nlin) - bulk(3) * LBulk2(1,1:nlin) ) / bulk(1)

!  original code
!      LRFE_BMie_bulk(1:3,1:nlin,1) = FF1 * LRFE_Mie1_bulk(1:3,1:nlin)
!      LRFE_BMie_bulk(1:3,1:nlin,2) = FF2 * LRFE_Mie2_bulk(1:3,1:nlin)

!  Asymmetry parameter
!  -------------------

      if ( Do_Expcoeffs ) then
         do q = 1, nlin
            TERM1 = WW1 * Mie_Lin_Outputs_Mode1%LRFE_Mie_asymm(q)
            TERM2 = D1_WW1(q) * Mie_Outputs_Mode1%Mie_asymm + D1_WW2(q) * Mie_Outputs_Mode2%Mie_asymm
            Mie_Lin_Outputs_Bimodal_Mode1%LRFE_Mie_asymm(q) = TERM1 + TERM2
            TERM1 = WW2 * Mie_Lin_Outputs_Mode2%LRFE_Mie_asymm(q)
            TERM2 = D2_WW1(q) * Mie_Outputs_Mode1%Mie_asymm + D2_WW2(q) * Mie_Outputs_Mode2%Mie_asymm
            Mie_Lin_Outputs_Bimodal_Mode2%LRFE_Mie_asymm(q) = TERM1 + TERM2
         enddo
      endif

!  Expansion coefficients
!  ----------------------

      if ( Do_Expcoeffs ) then

!  w.r.t  microphysical parameters for Particles in Mode 1

         do q = 1, nlin
            do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW1 * Mie_Lin_Outputs_Mode1%LRFE_Mie_expcoeffs(k,L,q)
                  TERM2 = D1_WW1(q) * Mie_Outputs_Mode1%Mie_expcoeffs(k,L) + &
                          D1_WW2(q) * Mie_Outputs_Mode2%Mie_expcoeffs(k,L)
                  Mie_Lin_Outputs_Bimodal_Mode1%LRFE_Mie_expcoeffs(k,L,q) = TERM1 + TERM2
               enddo
            enddo
            if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
               do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
                  do k = 1, 6
                     TERM2 = D1_WW2(q) * Mie_Outputs_Mode2%Mie_expcoeffs(k,L)
                     Mie_Lin_Outputs_Bimodal_Mode1%LRFE_Mie_expcoeffs(k,L,q) = TERM2
                  enddo
               enddo  
            else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
               do L = Mie2_ncoeffs + 1,Mie1_ncoeffs
                  do k = 1, 6
                     TERM1 = WW1       * Mie_Lin_Outputs_Mode1%LRFE_Mie_expcoeffs(k,L,q)
                     TERM2 = D1_WW1(q) * Mie_Outputs_Mode1%Mie_expcoeffs(k,L)
                     Mie_Lin_Outputs_Bimodal_Mode1%LRFE_Mie_expcoeffs(k,L,q) = TERM1 + TERM2
                  enddo
               enddo  
            endif
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW2 * Mie_Lin_Outputs_Mode2%LRFE_Mie_expcoeffs(k,L,q)
                  TERM2 = D2_WW1(q) * Mie_Outputs_Mode1%Mie_expcoeffs(k,L) + &
                          D2_WW2(q) * Mie_Outputs_Mode2%Mie_expcoeffs(k,L)
                  Mie_Lin_Outputs_Bimodal_Mode2%LRFE_Mie_expcoeffs(k,L,q) = TERM1 + TERM2
               enddo
            enddo
            if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
               do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
                  do k = 1, 6
                     TERM1 = WW2       * Mie_Lin_Outputs_Mode2%LRFE_Mie_expcoeffs(k,L,q)
                     TERM2 = D2_WW2(q) * Mie_Outputs_Mode2%Mie_expcoeffs(k,L)
                     Mie_Lin_Outputs_Bimodal_Mode2%LRFE_Mie_expcoeffs(k,L,q) = TERM1 + TERM2
                  enddo
               enddo  
            else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
               do L = Mie2_ncoeffs + 1,Mie1_ncoeffs
                  do k = 1, 6
                     TERM1 = D2_WW1(q) * Mie_Outputs_Mode1%Mie_expcoeffs(k,L)
                     Mie_Lin_Outputs_Bimodal_Mode2%LRFE_Mie_expcoeffs(k,L,q) = TERM1
                  enddo
               enddo  
            endif
         enddo

      endif

!  Fmatrix
!  -------

      if ( Do_Fmatrix ) then

!  w.r.t  microphysical parameters for Particles in Mode 1

         do q = 1, nlin
            do L = 1, n_Fmatrix_angles
               do k = 1, 4
                  TERM1 = WW1 * Mie_Lin_Outputs_Mode1%LRFE_Mie_Fmatrix(k,L,q)
                  TERM2 = D1_WW1(q) * Mie_Outputs_Mode1%Mie_Fmatrix(k,L) + &
                          D1_WW2(q) * Mie_Outputs_Mode2%Mie_Fmatrix(k,L)
                  Mie_Lin_Outputs_Bimodal_Mode1%LRFE_Mie_Fmatrix(k,L,q) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 1, n_Fmatrix_angles
               do k = 1, 4
                  TERM1 = WW2 * Mie_Lin_Outputs_Mode2%LRFE_Mie_Fmatrix(k,L,q)
                  TERM2 = D2_WW1(q) * Mie_Outputs_Mode1%Mie_Fmatrix(k,L) + &
                          D2_WW2(q) * Mie_Outputs_Mode2%Mie_Fmatrix(k,L)
                  Mie_Lin_Outputs_Bimodal_Mode2%LRFE_Mie_Fmatrix(k,L,q) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  End clauses

      endif   
   endif

!  LPSD Linearizations
!  ===================

   if ( do_LinearPSD ) then

!  help variables
!  --------------

!  Bug corrected 5/25/20. Thanks to Lukas Tirpitz.
!  Same Bug discovered here by R. Spurr, June 2014, but not implemented generally
!   -- Ordering of these two lines swapped.
!      D2_WW1(1:nlin) = - D2_WW2(1:nlin)
!      D2_WW2(1:nlin) = FF2 * D2_Csca2(1:nlin) * WW1 / Csca_total

      nlin = 3
      D1_Csca1(1:nlin) = Mie_Lin_Outputs_Mode1%LPSD_Mie_bulk(2,1:nlin)
      D2_Csca2(1:nlin) = Mie_Lin_Outputs_Mode2%LPSD_Mie_bulk(2,1:nlin)
      D1_WW1(1:nlin) = FF1 * D1_Csca1(1:nlin) * WW2 / Csca_total
      D1_WW2(1:nlin) = - D1_WW1(1:nlin)
      D2_WW2(1:nlin) = FF2 * D2_Csca2(1:nlin) * WW1 / Csca_total
      D2_WW1(1:nlin) = - D2_WW2(1:nlin)

!  Bulk variables
!  --------------

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong

      do q = 1, nlin
        Mie_Lin_Outputs_Bimodal_Mode1%LPSD_Mie_bulk(1:2,q) = FF1 * Mie_Lin_Outputs_Mode1%LPSD_Mie_bulk(1:2,q)
        Mie_Lin_Outputs_Bimodal_Mode2%LPSD_Mie_bulk(1:2,q) = FF2 * Mie_Lin_Outputs_Mode2%LPSD_Mie_bulk(1:2,q)
        Bulk          = Mie_Outputs_Bimodal%Mie_bulk
        LBulk1(1:2,q) = Mie_Lin_Outputs_Mode1%LPSD_Mie_bulk(1:2,q)
        LBulk2(1:2,q) = Mie_Lin_Outputs_Mode2%LPSD_Mie_bulk(1:2,q)
        Mie_Lin_Outputs_Bimodal_Mode1%LPSD_Mie_bulk(3,q) = ( LBulk1(2,q) - bulk(3) * LBulk1(1,q) ) / bulk(1)
        Mie_Lin_Outputs_Bimodal_Mode2%LPSD_Mie_bulk(3,q) = ( LBulk2(2,q) - bulk(3) * LBulk2(1,q) ) / bulk(1)
      enddo

!  original code
!      LPSD_BMie_bulk(1:3,1:nlin,1) = FF1 * LPSD_Mie1_bulk(1:3,1:nlin)
!      LPSD_BMie_bulk(1:3,1:nlin,2) = FF2 * LPSD_Mie2_bulk(1:3,1:nlin)

!  Asymmetry parameter
!  -------------------

      if ( Do_Expcoeffs ) then
         do q = 1, nlin
            TERM1 = WW1 * Mie_Lin_Outputs_Mode1%LPSD_Mie_asymm(q)
            TERM2 = D1_WW1(q) * Mie_Outputs_Mode1%Mie_asymm + D1_WW2(q) * Mie_Outputs_Mode2%Mie_asymm
            Mie_Lin_Outputs_Bimodal_Mode1%LPSD_Mie_asymm(q) = TERM1 + TERM2
            TERM1 = WW2 * Mie_Lin_Outputs_Mode2%LPSD_Mie_asymm(q)
            TERM2 = D2_WW1(q) * Mie_Outputs_Mode1%Mie_asymm + D2_WW2(q) * Mie_Outputs_Mode2%Mie_asymm
            Mie_Lin_Outputs_Bimodal_Mode2%LPSD_Mie_asymm(q) = TERM1 + TERM2
         enddo
      endif

!  Expansion coefficients
!  ----------------------

      if ( Do_Expcoeffs ) then

!  w.r.t  microphysical parameters for Particles in Mode 1

         do q = 1, nlin
            do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW1 * Mie_Lin_Outputs_Mode1%LPSD_Mie_expcoeffs(k,L,q)
                  TERM2 = D1_WW1(q) * Mie_Outputs_Mode1%Mie_expcoeffs(k,L) + &
                          D1_WW2(q) * Mie_Outputs_Mode2%Mie_expcoeffs(k,L)
                  Mie_Lin_Outputs_Bimodal_Mode1%LPSD_Mie_expcoeffs(k,L,q) = TERM1 + TERM2
               enddo
            enddo
            if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
               do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
                  do k = 1, 6
                     TERM2 = D1_WW2(q) * Mie_Outputs_Mode2%Mie_expcoeffs(k,L)
                     Mie_Lin_Outputs_Bimodal_Mode1%LPSD_Mie_expcoeffs(k,L,q) = TERM2
                  enddo
               enddo  
            else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
               do L = Mie2_ncoeffs + 1,Mie1_ncoeffs
                  do k = 1, 6
                     TERM1 = WW1       * Mie_Lin_Outputs_Mode1%LPSD_Mie_expcoeffs(k,L,q)
                     TERM2 = D1_WW1(q) * Mie_Outputs_Mode1%Mie_expcoeffs(k,L)
                     Mie_Lin_Outputs_Bimodal_Mode1%LPSD_Mie_expcoeffs(k,L,q) = TERM1 + TERM2
                  enddo
               enddo  
            endif
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW2 * Mie_Lin_Outputs_Mode2%LPSD_Mie_expcoeffs(k,L,q)
                  TERM2 = D2_WW1(q) * Mie_Outputs_Mode1%Mie_expcoeffs(k,L) + &
                          D2_WW2(q) * Mie_Outputs_Mode2%Mie_expcoeffs(k,L)
                  Mie_Lin_Outputs_Bimodal_Mode2%LPSD_Mie_expcoeffs(k,L,q) = TERM1 + TERM2
               enddo
            enddo
            if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
               do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
                  do k = 1, 6
                     TERM1 = WW2       * Mie_Lin_Outputs_Mode2%LPSD_Mie_expcoeffs(k,L,q)
                     TERM2 = D2_WW2(q) * Mie_Outputs_Mode2%Mie_expcoeffs(k,L)
                     Mie_Lin_Outputs_Bimodal_Mode2%LPSD_Mie_expcoeffs(k,L,q) = TERM1 + TERM2
                  enddo
               enddo  
            else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
               do L = Mie2_ncoeffs + 1,Mie1_ncoeffs
                  do k = 1, 6
                     TERM1 = D2_WW1(q) * Mie_Outputs_Mode1%Mie_expcoeffs(k,L)
                     Mie_Lin_Outputs_Bimodal_Mode2%LPSD_Mie_expcoeffs(k,L,q) = TERM1
                  enddo
               enddo  
            endif
         enddo

      endif

!  Fmatrix
!  -------

      if ( Do_Fmatrix ) then

!  w.r.t  microphysical parameters for Particles in Mode 1

         do q = 1, nlin
            do L = 1, n_Fmatrix_angles
               do k = 1, 4
                  TERM1 = WW1 * Mie_Lin_Outputs_Mode1%LPSD_Mie_Fmatrix(k,L,q)
                  TERM2 = D1_WW1(q) * Mie_Outputs_Mode1%Mie_Fmatrix(k,L) + &
                          D1_WW2(q) * Mie_Outputs_Mode2%Mie_Fmatrix(k,L)
                  Mie_Lin_Outputs_Bimodal_Mode1%LPSD_Mie_Fmatrix(k,L,q) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 1, n_Fmatrix_angles
               do k = 1, 4
                  TERM1 = WW2 * Mie_Lin_Outputs_Mode2%LPSD_Mie_Fmatrix(k,L,q)
                  TERM2 = D2_WW1(q) * Mie_Outputs_Mode1%Mie_Fmatrix(k,L) + &
                          D2_WW2(q) * Mie_Outputs_Mode2%Mie_Fmatrix(k,L)
                  Mie_Lin_Outputs_Bimodal_Mode2%LPSD_Mie_Fmatrix(k,L,q) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  End clauses

      endif   
   endif

!  Fractional linearization. NOT NORMALIZED
!  ========================================

!  Help variables

   DF_WW1 = ( WW1 * Mie_Outputs_Mode2%Mie_Bulk(2) + WW2 * Mie_Outputs_Mode1%Mie_Bulk(2) ) / Csca_total 

!  Bulk

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong

   Bulk        = Mie_Outputs_Bimodal%Mie_Bulk
   LBulkf(1:2) = Mie_Outputs_Mode1%Mie_Bulk(1:2) - Mie_Outputs_Mode2%Mie_Bulk(1:2)
   Mie_Lin_Outputs_Bimodal_Mode1%LFRC_Mie_Bulk(1:2) = LBulkf(1:2)
   Mie_Lin_Outputs_Bimodal_Mode2%LFRC_Mie_Bulk(1:2) = LBulkf(1:2)

   L_omega = ( LBulkf(2) - Bulk(3) * LBulkf(1) ) / Bulk(1)
   Mie_Lin_Outputs_Bimodal_Mode1%LFRC_Mie_Bulk(3) = L_omega
   Mie_Lin_Outputs_Bimodal_Mode2%LFRC_Mie_Bulk(3) = L_omega

!  Coefficients, asymmetry parameter

   if ( Do_Expcoeffs ) then
      Dasy = Mie_Outputs_Mode1%Mie_Asymm - Mie_Outputs_Mode2%Mie_Asymm
      Mie_Lin_Outputs_Bimodal_Mode1%LFRC_Mie_asymm = DF_WW1 * Dasy
      Mie_Lin_Outputs_Bimodal_Mode2%LFRC_Mie_asymm = DF_WW1 * Dasy
      do L = 0, Mie_Outputs_Bimodal%Mie_ncoeffs
         Dcoeff(1:6) = Mie_Outputs_Mode1%Mie_expcoeffs(1:6,L) - Mie_Outputs_Mode2%Mie_expcoeffs(1:6,L)
         Mie_Lin_Outputs_Bimodal_Mode1%LFRC_Mie_expcoeffs(1:6,L) = DF_WW1 * Dcoeff(1:6)
         Mie_Lin_Outputs_Bimodal_Mode2%LFRC_Mie_expcoeffs(1:6,L) = DF_WW1 * Dcoeff(1:6)
      enddo
   endif

!  Fmatrix

   if ( Do_Fmatrix ) then
      do L = 1, n_Fmatrix_angles
         DFmat(1:4) = Mie_Outputs_Mode1%Mie_Fmatrix(1:4,L) - Mie_Outputs_Mode2%Mie_Fmatrix(1:4,L)
         Mie_Lin_Outputs_Bimodal_Mode1%LFRC_Mie_Fmatrix(1:4,L) = DFmat(1:4)
         Mie_Lin_Outputs_Bimodal_Mode2%LFRC_Mie_Fmatrix(1:4,L) = DFmat(1:4)
      enddo
   endif   

!  Finish

   return
end subroutine RTSMie_master_bimodal_plus

!  End module

end module RTSMie_master_bimodal_plus_m



