program Mie_2pars_Sensitivity

!  Mie modules

      use RTSMie_parameters_m
      use RTSMie_Inputs_Def_m
      use RTSMie_Outputs_Def_m
      use RTSMie_Master_m

      implicit none

      integer, parameter :: mpk = SELECTED_REAL_KIND(15)

!  Type structures
!  ---------------

!  Mie Global input structure

      Type(RTSMie_Inputs)   :: Mie_Inputs
      Type(RTSMie_Outputs)  :: Mie_Outputs

!  top level configuration inputs
!  ------------------------------

      real(mpk)     :: AOD_refw
      real(mpk)     :: nreal_firstguess, nimag_firstguess
      real(mpk)     :: Rg_firstguess, Sg_firstguess

!  Local
!  -----

      integer :: istatus, i, j
      logical :: fail
      character*100 :: message, trace, action
      real(mpk)     :: Mie_Distplot(17,2500), Rgvals(16)
      data Rgvals / 0.20, 0.23, 0.25, 0.27, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.37, 0.40, 0.43, 0.45, 0.47, 0.50 /

!  Top level configuration read
!  ----------------------------

      OPEN(1,file='MieDists_2pars_Sensitivity.inp', STATUS = 'old')
      read(1,*)AOD_refw          ! Reference wavelength [nm]. Should be 300.0
      read(1,*)nreal_firstguess  ! Obvious
      read(1,*)nimag_firstguess  ! Obvious
      read(1,*)Rg_firstguess     ! Obvious
      read(1,*)Sg_firstguess     ! Obvious
      close(1)

!  Loop

   do j = 1, 16

!  Vary

!      Sg_firstguess = 1.2d0 + 0.05d0 * dble(j-1)
      Rg_firstguess = Rgvals(j)


!  Mie General flags

      Mie_Inputs%do_Expcoeffs    = .false.   ! 
      Mie_Inputs%do_Fmatrix      = .false.   ! Not needed
      Mie_Inputs%do_monodisperse = .false.   ! Not needed
      Mie_Inputs%MonoRadius      = 0.0d0     ! Not needed
      Mie_Inputs%do_LinearREF    = .false.   ! T
      Mie_Inputs%do_LinearPSD    = .false.   ! 

!  Mie PSD is Lonormal with Rg and Sg given

      Mie_Inputs%PSD_index    = 4               ! Lognormal
      Mie_Inputs%PSD_pars(1)  = Rg_firstguess ! Rg
      Mie_Inputs%PSD_pars(2)  = Sg_Firstguess ! Sg
      Mie_Inputs%PSD_pars(3)  = 0.0d0           ! Not needed

!  wavelength and refractive index, initial global settings. These will change locally.

      Mie_Inputs%lambda     = AOD_refw * 0.001d0  ! Microns. Set later on, according to UV wavelengths
      Mie_Inputs%n_real     = nreal_firstguess
      Mie_Inputs%n_imag     = nimag_firstguess

!  PSD integration. Fairly standard set of conditions

      Mie_Inputs%nblocks  = 40
      Mie_Inputs%nweights = 32
      Mie_Inputs%xparticle_limit = 10000.0d0
      Mie_Inputs%R1R2_cutoff  = 1.0d-06
      Mie_Inputs%FixR1R2      = .false.
      Mie_Inputs%R1 = 0.001d0
      Mie_Inputs%R2 = 10.0d0

!  Angular outputs (Not needed)

      Mie_Inputs%n_Fmatrix_angles = 0
      Mie_Inputs%Fmatrix_angles   = 0.0d0

!  Call

     CALL RTSMie_Master  & 
       ( Mie_Inputs, Mie_Outputs,                 & ! I/O type structures
         fail, istatus, message, trace, action )    ! Exception handling

!  Exception handling

      if ( fail ) then
         write(*,*)'Mie failure - here are the error messages ---'
         write(*,*)'  ---- ',Trim(message)
         write(*,*)'  ---- ',Trim(trace)
         write(*,*)'  ---- ',Trim(action)
         write(*,*)'stop program !!'; stop
      endif

!      write(*,'(A,f6.3,A,1pe17.8)')'Done Sg = ',Sg_firstguess,' ; DIST NORM = ',Mie_Outputs%Mie_Dist1(1)
      write(*,'(A,f6.3,A,1pe17.8)')'Doing Rg = ',Rg_firstguess,' ; DIST NORM = ',Mie_Outputs%Mie_Dist1(1)

!  results

      do i = 1, Mie_Inputs%nblocks * Mie_Inputs%nweights
         Mie_DistPlot(1,i)   = Mie_Outputs%Mie_DistPlot(1,i)
         Mie_DistPlot(1+j,i) = Mie_Outputs%Mie_DistPlot(2,i)
!         write(98,'(i5,2x,1p2e15.6)')i,Mie_Outputs%Mie_DistPlot(1,i),Mie_Outputs%Mie_DistPlot(2,i)
      enddo

   enddo

!  write results

!   open(98,file='Mie_Sg_Distributions.dat',status='replace')
   open(98,file='Mie_Rg_Distributions.dat',status='replace')
      do i = 1, Mie_Inputs%nblocks * Mie_Inputs%nweights
         write(98,'(i4,x,1pe15.6,3x,1p16e15.6)')i,Mie_DistPlot(1,i),Mie_DistPlot(2:17,i)
      enddo
   close(98)

!  Normal finish

   write(*,*)' **** Successful finish '
   stop

!  done

end program Mie_2pars_Sensitivity

