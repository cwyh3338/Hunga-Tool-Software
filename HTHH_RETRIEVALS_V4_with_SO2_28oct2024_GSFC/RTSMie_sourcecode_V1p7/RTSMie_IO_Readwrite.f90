module RTSMie_IO_Readwrite_m

!  Mie modules

  use RTSMie_parameters_m
  use RTSMie_Inputs_Def_m
  use RTSMie_Outputs_Def_m

!  Auxiliary stand-alone routines for Reading and writing

!    RTSMie_read_configfile.          Reads the Regular configuration file
!    RTSMie_read_configbimodal.       Reads the Bimodal configuration file

!    RTSMie_write_standard.           Standard output
!    RTSMie_write_extended.           Standard + Linearized  output

!    RTSMie_write_bimodal_standard.   Standard output, bimodal
!    RTSMie_write_bimodal_extended.   Standard + Linearized output, bimodal

public

contains

subroutine RTSMie_read_configfile ( filename, Mie_Inputs )

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles
   use RTSMie_Inputs_Def_m

   implicit none

!  Filename

   character*(*), intent(in)  :: filename

!  Input Type structure to be filled

   TYPE(RTSMie_Inputs) , INTENT (inout)    :: Mie_Inputs

!  Local
!  -----

   integer       :: j
   real(kind=dp) :: div

!  Open file and read
!  ------------------

   open(1,file=TRIM(ADJUSTL(filename)),status = 'old')
   read(1,*)
   read(1,*)
   read(1,*) ! First group (Boolean flags)
   read(1,*)
   read(1,*) Mie_Inputs%Do_expcoeffs
   read(1,*) Mie_Inputs%Do_Fmatrix
   read(1,*) Mie_Inputs%Do_monodisperse
   read(1,*) Mie_Inputs%Do_LinearRef
   read(1,*) Mie_Inputs%Do_LinearPSD
   read(1,*)
   read(1,*) ! Second group (PSD-related variables)
   read(1,*)
   read(1,*) Mie_Inputs%psd_Index
   read(1,*) Mie_Inputs%psd_pars(1)
   read(1,*) Mie_Inputs%psd_pars(2)
   read(1,*) Mie_Inputs%psd_pars(3)
   read(1,*) Mie_Inputs%FixR1R2
   read(1,*) Mie_Inputs%R1
   read(1,*) Mie_Inputs%R2
   read(1,*) Mie_Inputs%R1R2_cutoff
   read(1,*) Mie_Inputs%nblocks
   read(1,*) Mie_Inputs%nweights
   read(1,*)
   read(1,*) ! Third group (Monodisperse and optical inputs)
   read(1,*)
   read(1,*) Mie_Inputs%xparticle_limit
   read(1,*) Mie_Inputs%Monoradius
   read(1,*) Mie_Inputs%lambda
   read(1,*) Mie_Inputs%n_real
   read(1,*) Mie_Inputs%n_imag
   read(1,*)
   read(1,*) ! Fourth group (Fmatrix inputs)
   read(1,*)
   read(1,*) Mie_Inputs%n_Fmatrix_angles
   close(1)

!  Fmatrix angles are regular here between 0 and 180.

   if ( Mie_Inputs%do_Fmatrix ) then
      div = 180.0d0 / dble( Mie_Inputs%n_Fmatrix_angles - 1 )
      do j = 1, Mie_Inputs%n_Fmatrix_angles
         Mie_Inputs%Fmatrix_angles(j) = dble(j-1) * div
      enddo
   else
      Mie_Inputs%Fmatrix_angles = 0.0d0
   endif

!  return

   return
end subroutine RTSMie_read_configfile

subroutine RTSMie_read_configbimodal ( filename, Mie_Inputs_Mode1, Mie_Inputs_Mode2, Bimodal_fraction )

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles
   use RTSMie_Inputs_Def_m

   implicit none

!  Filename

   character*(*), intent(in)  :: filename

!  Input Type structure to be filled

   TYPE(RTSMie_Inputs) , INTENT (inout)    :: Mie_Inputs_Mode1
   TYPE(RTSMie_Inputs) , INTENT (inout)    :: Mie_Inputs_Mode2

!  Bimodal fractional weight for PSD 1 

   REAL(KIND=dp), intent(out)  :: Bimodal_fraction

!  Local
!  -----

   integer       :: j
   real(kind=dp) :: div

!  Open file and read
!  ------------------

   open(1,file=TRIM(ADJUSTL(filename)),status = 'old')
   read(1,*)
   read(1,*)
   read(1,*) ! First group (Boolean flags)
   read(1,*)
   read(1,*) Mie_Inputs_Mode1%Do_expcoeffs     ; Mie_Inputs_Mode2%Do_expcoeffs     = Mie_Inputs_Mode1%Do_expcoeffs
   read(1,*) Mie_Inputs_Mode1%Do_Fmatrix       ; Mie_Inputs_Mode2%Do_Fmatrix       = Mie_Inputs_Mode1%Do_Fmatrix
   read(1,*) Mie_Inputs_Mode1%Do_monodisperse  ; Mie_Inputs_Mode2%Do_monodisperse  = Mie_Inputs_Mode1%Do_monodisperse
   read(1,*) Mie_Inputs_Mode1%Do_LinearRef     ; Mie_Inputs_Mode2%Do_LinearRef     = Mie_Inputs_Mode1%Do_LinearRef
   read(1,*) Mie_Inputs_Mode1%Do_LinearPSD     ; Mie_Inputs_Mode2%Do_LinearPSD     = Mie_Inputs_Mode1%Do_LinearPSD
   read(1,*)
   read(1,*) ! Second group (PSD-related variables)
   read(1,*)
   read(1,*) Mie_Inputs_Mode1%psd_Index,   Mie_Inputs_Mode2%psd_Index
   read(1,*) Mie_Inputs_Mode1%psd_pars(1), Mie_Inputs_Mode2%psd_pars(1)
   read(1,*) Mie_Inputs_Mode1%psd_pars(2), Mie_Inputs_Mode2%psd_pars(2)
   read(1,*) Mie_Inputs_Mode1%psd_pars(3), Mie_Inputs_Mode2%psd_pars(3)
   read(1,*) Mie_Inputs_Mode1%FixR1R2,     Mie_Inputs_Mode2%FixR1R2
   read(1,*) Mie_Inputs_Mode1%R1,          Mie_Inputs_Mode2%R1
   read(1,*) Mie_Inputs_Mode1%R2,          Mie_Inputs_Mode2%R2
   read(1,*) Mie_Inputs_Mode1%R1R2_cutoff, Mie_Inputs_Mode2%R1R2_cutoff
   read(1,*) Mie_Inputs_Mode1%nblocks,     Mie_Inputs_Mode2%nblocks
   read(1,*) Mie_Inputs_Mode1%nweights,    Mie_Inputs_Mode2%nweights
   read(1,*)
   read(1,*) ! Third group (Monodisperse and optical inputs)
   read(1,*)
   read(1,*) Mie_Inputs_Mode1%xparticle_limit ; Mie_Inputs_Mode2%xparticle_limit = Mie_Inputs_Mode1%xparticle_limit
   read(1,*) Mie_Inputs_Mode1%Monoradius      ; Mie_Inputs_Mode2%Monoradius      = Mie_Inputs_Mode1%Monoradius
   read(1,*) Mie_Inputs_Mode1%lambda          ; Mie_Inputs_Mode2%lambda          = Mie_Inputs_Mode1%lambda
   read(1,*) Mie_Inputs_Mode1%n_real, Mie_Inputs_Mode2%n_real
   read(1,*) Mie_Inputs_Mode1%n_imag, Mie_Inputs_Mode2%n_imag
   read(1,*)
   read(1,*) ! Fourth group (Fmatrix inputs)
   read(1,*)
   read(1,*) Mie_Inputs_Mode1%n_Fmatrix_angles ; Mie_Inputs_Mode2%n_Fmatrix_angles = Mie_Inputs_Mode1%n_Fmatrix_angles
   read(1,*)
   read(1,*) ! Fifth group (fractional input)
   read(1,*)
   read(1,*) Bimodal_fraction
   close(1)

!  Some settings fixed

   Mie_Inputs_Mode1%Do_monodisperse = .false.
   Mie_Inputs_Mode1%Monoradius      = 0.0d0
   Mie_Inputs_Mode2%Do_monodisperse = .false.
   Mie_Inputs_Mode2%Monoradius      = 0.0d0

!  Fmatrix angles are regular here between 0 and 180.

   if (  Mie_Inputs_Mode1%do_Fmatrix ) then
      div = 180.0d0 / dble(  Mie_Inputs_Mode1%n_Fmatrix_angles - 1 )
      do j = 1,  Mie_Inputs_Mode1%n_Fmatrix_angles
         Mie_Inputs_Mode1%Fmatrix_angles(j) = dble(j-1) * div
         Mie_Inputs_Mode2%Fmatrix_angles(j) = dble(j-1) * div
      enddo
   else
      Mie_Inputs_Mode1%Fmatrix_angles = 0.0d0
      Mie_Inputs_Mode2%Fmatrix_angles = 0.0d0
   endif

!  return

   return
end subroutine RTSMie_read_configbimodal


subroutine RTSMie_write_standard ( filename, Mie_Inputs, Mie_Outputs )

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles
   use RTSMie_Inputs_Def_m
   use RTSMie_Outputs_Def_m

   implicit none

!  routine inputs
!  ==============

!  Filename

   character*(*), intent(in)  :: filename

!  I/O Type structures

   TYPE(RTSMie_Inputs) , INTENT (in)    :: Mie_Inputs
   TYPE(RTSMie_Outputs), INTENT (in)    :: Mie_Outputs

!  Local

   integer       :: L, k, cmask(6),fmask(4)

!  Old indexing

!   data cmask / 1, 3, 4, 6, 2, 5 /
!   data fmask / 1, 3, 2, 4 /

!  New indexing

   data cmask / 1, 2, 3, 4, 5, 6 /
   data fmask / 1, 2, 3, 4 /

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   if ( Mie_Inputs%Do_monodisperse ) then
      write(35,'(/a)')' *** MONODISPERSE OUTPUT ONLY '
   else
      write(35,'(/a/)')' PSD output ------------- '
      write(35,'(A,T25,1pe15.6)')'Number density     = ',Mie_Outputs%Mie_dist1 (1)
      write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Mie_Outputs%Mie_dist1 (2)
      write(35,'(A,T25,1pe15.6)')'Volume             = ',Mie_Outputs%Mie_dist1 (3)
      write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Mie_Outputs%Mie_dist1 (4)
      write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Mie_Outputs%Mie_dist1 (5)
   endif

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction Coefficient = ',Mie_Outputs%Mie_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering Coefficient = ',Mie_Outputs%Mie_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Mie_Outputs%Mie_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Mie_Outputs%Mie_asymm

!  Expansion coefficients
  
   if ( Mie_Inputs%Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',Mie_Outputs%Mie_ncoeffs
      do L = 0, Mie_Outputs%Mie_ncoeffs
        write(35,'(i5,1p6e20.11)')L,(Mie_Outputs%Mie_expcoeffs (cmask(k),L),k=1,6)
      enddo
   endif

!  F-matrix

   if ( Mie_Inputs%Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',Mie_Inputs%n_Fmatrix_angles
      do L = 1, Mie_Inputs%n_Fmatrix_angles
        write(35,'(f6.2,4F20.11)')Mie_Inputs%Fmatrix_angles(L),(Mie_Outputs%Mie_Fmatrix(fmask(k),L),k=1,4)
      enddo
   endif

!  Close file

   close(35)

!  Finish

   return
end subroutine RTSMie_write_standard

subroutine RTSMie_write_extended ( filename, donorm, Mie_Inputs, Mie_Outputs, Mie_Lin_Outputs )

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles
   use RTSMie_Inputs_Def_m
   use RTSMie_Outputs_Def_m
   use RTSMie_Lin_Outputs_Def_m

   implicit none

!  routine inputs
!  ==============

!  Filename

   character*(*), intent(in)  :: filename

!  Normalization flag

   logical, intent(in) :: donorm

!  I/O Type structures

   TYPE(RTSMie_Inputs)     , INTENT (in)    :: Mie_Inputs
   TYPE(RTSMie_Outputs)    , INTENT (in)    :: Mie_Outputs
   TYPE(RTSMie_Lin_Outputs), INTENT (in)    :: Mie_Lin_Outputs

!  Local

   integer       :: L, k, M, MS, cmask(6),fmask(4)
   character*1   :: c1
   real(kind=dp) :: pf

!  Old indexing

!   data cmask / 1, 3, 4, 6, 2, 5 /
!   data fmask / 1, 3, 2, 4 /

!  New indexing

   data cmask / 1, 2, 3, 4, 5, 6 /
   data fmask / 1, 2, 3, 4 /

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   if ( Mie_Inputs%Do_monodisperse ) then
      write(35,'(/a)')' *** MONODISPERSE OUTPUT ONLY '
   else
      write(35,'(/a/)')' PSD output ------------- '
      write(35,'(A,T25,1pe15.6)')'Number density     = ',Mie_Outputs%Mie_dist1 (1)
      write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Mie_Outputs%Mie_dist1 (2)
      write(35,'(A,T25,1pe15.6)')'Volume             = ',Mie_Outputs%Mie_dist1 (3)
      write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Mie_Outputs%Mie_dist1 (4)
      write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Mie_Outputs%Mie_dist1 (5)
   endif

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction Coefficient = ',Mie_Outputs%Mie_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering Coefficient = ',Mie_Outputs%Mie_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Mie_Outputs%Mie_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Mie_Outputs%Mie_asymm

!  Expansion coefficients
  
   if ( Mie_Inputs%Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',Mie_Outputs%Mie_ncoeffs
      do L = 0, Mie_Outputs%Mie_ncoeffs
        write(35,'(i5,1p6e20.11)')L,(Mie_Outputs%Mie_expcoeffs (cmask(k),L),k=1,6)
      enddo
   endif

!  F-matrix

   if ( Mie_Inputs%Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',Mie_Inputs%n_Fmatrix_angles
      do L = 1, Mie_Inputs%n_Fmatrix_angles
        write(35,'(f6.2,4F20.11)')Mie_Inputs%Fmatrix_angles(L),(Mie_Outputs%Mie_Fmatrix(fmask(k),L),k=1,4)
      enddo
   endif

!  Close file

   close(35)

!  write PSD-linearized results
!    Filenames derived from Input name, with '_LPSD_*' added (where * = 1, 2, or 3)
 
   if ( .not. Mie_Inputs%do_monodisperse ) then
      if ( Mie_Inputs%Do_LinearPSD ) then
         MS = 0 ; if (Mie_Inputs%psd_index.eq.4) MS = 2 ; if (Mie_Inputs%psd_index.ne.4) MS = 3
         do M = 1, MS
            pf = 1.0d0 ; if (donorm) pf = Mie_Inputs%psd_pars(m)
            write(C1,'(I1)')M
            open(36,file=TRIM(ADJUSTL(filename))//'_LPSD_'//C1,status = 'unknown' )

            write(36,'(/a/)')' LPSD Distribution output ------------- '
            write(36,'(A,T25,1pe15.6)')'Number density     = ',pf*Mie_Lin_Outputs%LPSD_Mie_dist (1,m) 
            write(36,'(A,T25,1pe15.6)')'Cross-section      = ',pf*Mie_Lin_Outputs%LPSD_Mie_dist (2,m)
            write(36,'(A,T25,1pe15.6)')'Volume             = ',pf*Mie_Lin_Outputs%LPSD_Mie_dist (3,m)
            write(36,'(A,T25,1pe15.6)')'Effective Radius   = ',pf*Mie_Lin_Outputs%LPSD_Mie_dist (4,m)
            write(36,'(A,T25,1pe15.6)')'Effective Variance = ',pf*Mie_Lin_Outputs%LPSD_Mie_dist (5,m)

            write(36,'(/a/)')' LPSD Bulk property output ------------- '
            write(36,'(A,T40,1pe20.11)')'Extinction Coefficient = ',pf*Mie_Lin_Outputs%LPSD_Mie_bulk (1,m)
            write(36,'(A,T40,1pe20.11)')'Scattering Coefficient = ',pf*Mie_Lin_Outputs%LPSD_Mie_bulk (2,m)
            write(36,'(A,T40,1pe20.11)')'Singlescattering albedo= ',pf*Mie_Lin_Outputs%LPSD_Mie_bulk (3,m)
            write(36,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',pf*Mie_Lin_Outputs%LPSD_Mie_asymm(m)

            if ( Mie_Inputs%Do_expcoeffs ) then
               write(36,'(/a)')'LPSD Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
               write(36,'(a,I5/)')' - Number of coefficients =  ',Mie_Outputs%Mie_ncoeffs
               do L = 0, Mie_Outputs%Mie_ncoeffs
                  write(36,'(i5,1p6e20.11)')L,(pf*Mie_Lin_Outputs%LPSD_Mie_expcoeffs (cmask(k),L,M),k=1,6)
               enddo
            endif

            if ( Mie_Inputs%Do_Fmatrix ) then
               write(36,'(/a)')' LPSD F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
               write(36,'(a,I5/)')' - Number of angles =  ',Mie_Inputs%n_Fmatrix_angles
               do L = 1, Mie_Inputs%n_Fmatrix_angles
                  write(36,'(f6.2,4F20.11)')Mie_Inputs%Fmatrix_angles(L),&
                                      (pf*Mie_Lin_Outputs%LPSD_Mie_Fmatrix (fmask(k),L,m),k=1,4)
               enddo
            endif

            close(36)

         enddo
      endif
   endif

!  write RFE-linearized results
!    Filenames derived from Input name, with '_LRFE_*' added (where * = 1, 2)

   if ( Mie_Inputs%Do_LinearRef ) then
      do M = 1, 2
         pf = 1.0d0 ; if (donorm.and.m.eq.1) pf = Mie_Inputs%n_real ; if (donorm.and.m.eq.2) pf = Mie_Inputs%n_imag

         write(C1,'(I1)')M
         open(37,file=TRIM(ADJUSTL(filename))//'_LRFE_'//C1,status = 'unknown' )

         write(37,'(/a/)')' LRFE Bulk property output ------------- '
         write(37,'(A,T40,1pe20.11)')'Extinction Coefficient = ',pf*Mie_Lin_Outputs%LRFE_Mie_bulk (1,m)
         write(37,'(A,T40,1pe20.11)')'Scattering Coefficient = ',pf*Mie_Lin_Outputs%LRFE_Mie_bulk (2,m)
         write(37,'(A,T40,1pe20.11)')'Singlescattering albedo= ',pf*Mie_Lin_Outputs%LRFE_Mie_bulk (3,m)
         write(37,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',pf*Mie_Lin_Outputs%LRFE_Mie_asymm(m)

         if ( Mie_Inputs%Do_expcoeffs ) then
            write(37,'(/a)')'LRFE Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
            write(37,'(a,I5/)')' - Number of coefficients =  ',Mie_Outputs%Mie_ncoeffs
            do L = 0, Mie_Outputs%Mie_ncoeffs
               write(37,'(i5,1p6e20.11)')L,(pf*Mie_Lin_Outputs%LRFE_Mie_expcoeffs (cmask(k),L,M),k=1,6)
            enddo
         endif

         if ( Mie_Inputs%Do_Fmatrix ) then
            write(37,'(/a)')' LRFE F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
            write(37,'(a,I5/)')' - Number of angles =  ',Mie_Inputs%n_Fmatrix_angles
            do L = 1, Mie_Inputs%n_Fmatrix_angles
               write(37,'(f6.2,4F20.11)')Mie_Inputs%Fmatrix_angles(L),(pf*Mie_Lin_Outputs%LRFE_Mie_Fmatrix (fmask(k),L,m),k=1,4)
            enddo
         endif

         close(37)

      enddo
   endif

!  Finish

   return
end subroutine RTSMie_write_extended

subroutine RTSMie_write_bimodal_standard &
         ( filename, Bimodal_Fraction, do_Expcoeffs, do_Fmatrix, n_Fmatrix_angles, FMatrix_Angles, Mie_Outputs )

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles
   use RTSMie_Inputs_Def_m
   use RTSMie_Outputs_Def_m

   implicit none

!  routine inputs
!  ==============

!  Filename

   character*(*), intent(in)  :: filename

!  Bimodal Fraction

   REAL(KIND=dp), intent(in) :: Bimodal_Fraction

!  Input controls

   logical, intent(in)     :: Do_Expcoeffs
   logical, intent(in)     :: Do_Fmatrix

!  Input Fmatrix angles

   INTEGER          , intent(in) :: n_Fmatrix_angles
   REAL    (KIND=dp), intent(in) :: Fmatrix_angles(max_Mie_angles)

!  I/O Type structures

   TYPE(RTSMie_Outputs), INTENT (in)    :: Mie_Outputs


!   WRITES OUT BOTH DISTRIBUTIONS

!  Local

   integer       :: L, k, cmask(6),fmask(4)

!  Old indexing

!   data cmask / 1, 3, 4, 6, 2, 5 /
!   data fmask / 1, 3, 2, 4 /

!  New indexing

   data cmask / 1, 2, 3, 4, 5, 6 /
   data fmask / 1, 2, 3, 4 /

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   write(35,'(/a/)')' PSD output,  first mode  ------------- '
      write(35,'(A,T25,1pe15.6)')'Number density     = ',Mie_Outputs%Mie_dist1 (1)
      write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Mie_Outputs%Mie_dist1 (2)
      write(35,'(A,T25,1pe15.6)')'Volume             = ',Mie_Outputs%Mie_dist1 (3)
      write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Mie_Outputs%Mie_dist1 (4)
      write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Mie_Outputs%Mie_dist1 (5)

   write(35,'(/a/)')' PSD output, second mode  ------------- '
      write(35,'(A,T25,1pe15.6)')'Number density     = ',Mie_Outputs%Mie_dist2 (1)
      write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Mie_Outputs%Mie_dist2 (2)
      write(35,'(A,T25,1pe15.6)')'Volume             = ',Mie_Outputs%Mie_dist2 (3)
      write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Mie_Outputs%Mie_dist2 (4)
      write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Mie_Outputs%Mie_dist2 (5)

   write(35,'(/a/)')' PSD output, Bimodal fraction  ------------- '
   write(35,'(A,T25,f10.5)')  'Bi-modal fraction  = ',Bimodal_fraction

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' TOTAL Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction Coefficient = ',Mie_Outputs%Mie_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering Coefficient = ',Mie_Outputs%Mie_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Mie_Outputs%Mie_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Mie_Outputs%Mie_asymm

!  Expansion coefficients
  
   if ( Do_expcoeffs ) then
      write(35,'(/a)')' TOTAL Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',Mie_Outputs%Mie_ncoeffs
      do L = 0, Mie_Outputs%Mie_ncoeffs
        write(35,'(i5,1p6e20.11)')L,(Mie_Outputs%Mie_expcoeffs (cmask(k),L),k=1,6)
      enddo
   endif

!  F-matrix

   if ( Do_Fmatrix ) then
      write(35,'(/a)')' TOTAL F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
      do L = 1, n_Fmatrix_angles
        write(35,'(f6.2,4F20.11)')Fmatrix_angles(L),(Mie_Outputs%Mie_Fmatrix(fmask(k),L),k=1,4)
      enddo
   endif

!  Close file

   close(35)

!  Finish

   return
end subroutine RTSMie_write_bimodal_standard

subroutine RTSMie_write_bimodal_extended &
         ( filename, Bimodal_Fraction, do_LinearREF, do_LinearPSD,      &
           do_Expcoeffs, do_Fmatrix, n_Fmatrix_angles, FMatrix_Angles,  &
           Mie_Outputs, Mie_Lin_Outputs_M1, Mie_Lin_Outputs_M2 )

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles
   use RTSMie_Inputs_Def_m
   use RTSMie_Outputs_Def_m
   use RTSMie_Lin_Outputs_Def_m

   implicit none

!  routine inputs
!  ==============

!  Filename

   character*(*), intent(in)  :: filename

!  Bimodal Fraction

   REAL(KIND=dp), intent(in) :: Bimodal_Fraction

!  Input controls

   logical, intent(in)     :: Do_Expcoeffs
   logical, intent(in)     :: Do_Fmatrix
   logical, intent(in)     :: do_LinearREF
   logical, intent(in)     :: do_LinearPSD

!  Input Fmatrix angles

   INTEGER          , intent(in) :: n_Fmatrix_angles
   REAL    (KIND=dp), intent(in) :: Fmatrix_angles(max_Mie_angles)

!  Type structures

   TYPE(RTSMie_Outputs)    , INTENT (in)    :: Mie_Outputs
   TYPE(RTSMie_Lin_Outputs), INTENT (in)    :: Mie_Lin_Outputs_M1
   TYPE(RTSMie_Lin_Outputs), INTENT (in)    :: Mie_Lin_Outputs_M2

!  Local

   integer       :: L, k, M, J, cmask(6),fmask(4)
   character*1   :: c1, CJ
   character*22  :: CPSD

!  Local structure

   TYPE(RTSMie_Lin_Outputs)  :: Mie_Lin_Outputs

!  Old indexing
!   data cmask / 1, 3, 4, 6, 2, 5 /
!   data fmask / 1, 3, 2, 4 /

!  New indexing

   data cmask / 1, 2, 3, 4, 5, 6 /
   data fmask / 1, 2, 3, 4 /

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

   write(35,'(/a/)')' PSD output,  first mode  ------------- '
      write(35,'(A,T25,1pe15.6)')'Number density     = ',Mie_Outputs%Mie_dist1 (1)
      write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Mie_Outputs%Mie_dist1 (2)
      write(35,'(A,T25,1pe15.6)')'Volume             = ',Mie_Outputs%Mie_dist1 (3)
      write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Mie_Outputs%Mie_dist1 (4)
      write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Mie_Outputs%Mie_dist1 (5)

   write(35,'(/a/)')' PSD output, second mode  ------------- '
      write(35,'(A,T25,1pe15.6)')'Number density     = ',Mie_Outputs%Mie_dist2 (1)
      write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Mie_Outputs%Mie_dist2 (2)
      write(35,'(A,T25,1pe15.6)')'Volume             = ',Mie_Outputs%Mie_dist2 (3)
      write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Mie_Outputs%Mie_dist2 (4)
      write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Mie_Outputs%Mie_dist2 (5)

   write(35,'(/a/)')' PSD output, Bimodal fraction  ------------- '
   write(35,'(A,T25,f10.5)')  'Bi-modal fraction  = ',Bimodal_fraction

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction Coefficient = ',Mie_Outputs%Mie_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering Coefficient = ',Mie_Outputs%Mie_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Mie_Outputs%Mie_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Mie_Outputs%Mie_asymm

!  Expansion coefficients
  
   if ( Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',Mie_Outputs%Mie_ncoeffs
      do L = 0, Mie_Outputs%Mie_ncoeffs
        write(35,'(i5,1p6e20.11)')L,(Mie_Outputs%Mie_expcoeffs (cmask(k),L),k=1,6)
      enddo
   endif

!  F-matrix

   if ( Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
      do L = 1, n_Fmatrix_angles
        write(35,'(f6.2,4F20.11)')Fmatrix_angles(L),(Mie_Outputs%Mie_Fmatrix(fmask(k),L),k=1,4)
      enddo
   endif

!  write PSD-linearized results
!    Filenames derived from Input name, with '_LPSD_*' added (where * = 1, 2, or 3)
 
   DO J = 1, 2
      write(CJ,'(I1)')J
      if ( J.eq.1 ) Mie_Lin_Outputs = Mie_Lin_Outputs_M1
      if ( J.eq.2 ) Mie_Lin_Outputs = Mie_Lin_Outputs_M2
      if ( Do_LinearPSD ) then
         do M = 1, 3 
            write(C1,'(I1)')M
            CPSD = 'PSD # '//CJ//', parameter # '//C1
            open(36,file=TRIM(ADJUSTL(filename))//'_LPSD_output_PSD#'//CJ//'_PAR#'//C1,status='unknown')
          
            write(36,'(/a/)')' LPSD Distribution output, PSD # '//CJ//' ------------- '
            write(36,'(A,T25,1pe15.6)')'Number density     = ',Mie_Lin_Outputs%LPSD_Mie_dist (1,m) 
            write(36,'(A,T25,1pe15.6)')'Cross-section      = ',Mie_Lin_Outputs%LPSD_Mie_dist (2,m)
            write(36,'(A,T25,1pe15.6)')'Volume             = ',Mie_Lin_Outputs%LPSD_Mie_dist (3,m)
            write(36,'(A,T25,1pe15.6)')'Effective Radius   = ',Mie_Lin_Outputs%LPSD_Mie_dist (4,m)
            write(36,'(A,T25,1pe15.6)')'Effective Variance = ',Mie_Lin_Outputs%LPSD_Mie_dist (5,m)

            write(36,'(/a/)')' LPSD Bulk property output, PSD # '//CJ//' ------------- '
            write(36,'(A,T40,1pe20.11)')'Extinction Coefficient = ',Mie_Lin_Outputs%LPSD_Mie_bulk (1,m)
            write(36,'(A,T40,1pe20.11)')'Scattering Coefficient = ',Mie_Lin_Outputs%LPSD_Mie_bulk (2,m)
            write(36,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Mie_Lin_Outputs%LPSD_Mie_bulk (3,m)
            write(36,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Mie_Lin_Outputs%LPSD_Mie_asymm(m)

            if ( Do_expcoeffs ) then
               write(36,'(/a)')'LPSD Expansion coefficient output (a1,a2,a3,a4,b1,b2), PSD # '//CJ
               write(36,'(a,I5/)')' - Number of coefficients =  ',Mie_Outputs%Mie_ncoeffs
               do L = 0, Mie_Outputs%Mie_ncoeffs
                  write(36,'(i5,1p6e20.11)')L,(Mie_Lin_Outputs%LPSD_Mie_expcoeffs (cmask(k),L,M),k=1,6)
               enddo
            endif

            if ( Do_Fmatrix ) then
               write(36,'(/a)')' LPSD F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs), PSD # '//CJ
               write(36,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
               do L = 1, n_Fmatrix_angles
                  write(36,'(f6.2,4F20.11)')Fmatrix_angles(L),(Mie_Lin_Outputs%LPSD_Mie_Fmatrix (fmask(k),L,m),k=1,4)
               enddo
            endif
            close(36)
         enddo
      endif

!  write RFE-linearized results
!    Filenames derived from Input name, with '_LRFE_*' added (where * = 1, 2)

      if ( Do_LinearRef ) then
         do M = 1, 2
            write(C1,'(I1)')M
            CPSD = 'PSD # '//CJ//', parameter # '//C1
            open(37,file=TRIM(ADJUSTL(filename))//'_LRFE_output_PSD#'//CJ//'_PAR#'//C1,status = 'unknown' )
          
            write(37,'(/a/)')' LRFE Bulk property output, RFE # '//CJ//' ------------- '
            write(37,'(A,T40,1pe20.11)')'Extinction Coefficient = ',Mie_Lin_Outputs%LRFE_Mie_bulk (1,m)
            write(37,'(A,T40,1pe20.11)')'Scattering Coefficient = ',Mie_Lin_Outputs%LRFE_Mie_bulk (2,m)
            write(37,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Mie_Lin_Outputs%LRFE_Mie_bulk (3,m)
            write(37,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Mie_Lin_Outputs%LRFE_Mie_asymm(m)

            if ( Do_expcoeffs ) then
               write(37,'(/a)')'LRFE Expansion coefficient output (a1,a2,a3,a4,b1,b2), RFE # '//CJ
               write(37,'(a,I5/)')' - Number of coefficients =  ',Mie_Outputs%Mie_ncoeffs
               do L = 0, Mie_Outputs%Mie_ncoeffs
                  write(37,'(i5,1p6e20.11)')L,(Mie_Lin_Outputs%LRFE_Mie_expcoeffs (cmask(k),L,M),k=1,6)
               enddo
            endif

            if ( Do_Fmatrix ) then
               write(37,'(/a)')' LRFE F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs), RFE # '//CJ
               write(37,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
               do L = 1, n_Fmatrix_angles
                  write(37,'(f6.2,4F20.11)')Fmatrix_angles(L),(Mie_Lin_Outputs%LRFE_Mie_Fmatrix (fmask(k),L,m),k=1,4)
               enddo
            endif

            close(37)

         enddo
      endif

!  End double distribution loop

   enddo

!  write FRC-linearized results (w.r.t. fraction)

   open(38,file=TRIM(ADJUSTL(filename))//'_LFRC',status = 'unknown' )

   write(38,'(/a/)')' LFRC Bulk property output ------------- '
   write(38,'(A,T40,1pe20.11)')'Extinction Coefficient = ',Mie_Lin_Outputs_M1%LFRC_Mie_bulk(1)
   write(38,'(A,T40,1pe20.11)')'Scattering Coefficient = ',Mie_Lin_Outputs_M1%LFRC_Mie_bulk(2)
   write(38,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Mie_Lin_Outputs_M1%LFRC_Mie_bulk(3)
   write(38,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Mie_Lin_Outputs_M1%LFRC_Mie_asymm

   if ( Do_expcoeffs ) then
      write(38,'(/a)')'LFRC Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(38,'(a,I5/)')' - Number of coefficients =  ',Mie_Outputs%Mie_ncoeffs
      do L = 0, Mie_Outputs%Mie_ncoeffs
         write(38,'(i5,1p6e20.11)')L,(Mie_Lin_Outputs_M1%LFRC_Mie_expcoeffs (cmask(k),L),k=1,6)
      enddo
   endif

   if ( Do_Fmatrix ) then
      write(38,'(/a)')' LFRC F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
      write(38,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
      do L = 1, n_Fmatrix_angles
         write(38,'(f6.2,4F20.11)')Fmatrix_angles(L),(Mie_Lin_Outputs_M1%LFRC_Mie_Fmatrix (fmask(k),L),k=1,4)
      enddo
   endif

   close(38)

!  Finish

   return
end subroutine RTSMie_write_bimodal_extended

!  Finish Module

end module RTSMie_IO_Readwrite_m

