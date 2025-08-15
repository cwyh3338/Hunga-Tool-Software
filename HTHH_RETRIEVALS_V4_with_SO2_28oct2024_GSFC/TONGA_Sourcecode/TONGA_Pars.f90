
!  File of constants for TONGA retrieval.
!mick mod 9/12/2022 - added some number parameters to have options independent of VLIDORT_pars

      MODULE TONGA_pars_m

      IMPLICIT NONE

!  Real number type definitions

      INTEGER, PARAMETER :: MPK = SELECTED_REAL_KIND(15)

!  Maximum number of layers

      !INTEGER, PARAMETER :: E_MAXLAYERS   = 120
      INTEGER, PARAMETER :: E_MAXLAYERS   = 160

!  Maximum number of wavelengths

      !INTEGER, PARAMETER :: E_MAXDATAWAVS = 500
      INTEGER, PARAMETER :: E_MAXDATAWAVS = 1000

      !INTEGER, PARAMETER :: E_MAXWAVS     = 300
      INTEGER, PARAMETER :: E_MAXWAVS     = 600

!  Maximum number of retrieval state vector elements

      INTEGER, PARAMETER :: E_MAXPARS     = 5

!  Aerosol data
!    -- 7/28/22. Increase the number of coefficients

      INTEGER, PARAMETER :: E_MAXAERWAVS  = 37
      INTEGER, PARAMETER :: E_MAXAERCOEFS = 1000

!  Debug write format (DWF) constants

      CHARACTER (LEN=*), PARAMETER :: DWFL  = '(A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL1 = '(A,I3,A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL2 = '(2(A,I3),A,L1)'

      CHARACTER (LEN=*), PARAMETER :: DWFI  = '(A,I5)'
      CHARACTER (LEN=*), PARAMETER :: DWFI1 = '(A,I3,A,I5)'
      CHARACTER (LEN=*), PARAMETER :: DWFI2 = '(2(A,I3),A,I5)'

      CHARACTER (LEN=*), PARAMETER :: DWFR  = '(A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR1 = '(A,I3,A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR2 = '(2(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR3 = '(3(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR4 = '(4(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR5 = '(5(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR6 = '(6(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR7 = '(7(A,I3),A,ES13.6E2)'

      CHARACTER (LEN=*), PARAMETER :: DWFR1_3 = '(A,I3,3(A,ES13.6E2))'

      !"Enhanced # of digits" parameter set for checking reals
      !CHARACTER (LEN=*), PARAMETER :: DWFR  = '(A,ES21.14E2)'
      !CHARACTER (LEN=*), PARAMETER :: DWFR1 = '(A,I3,A,ES21.14E2)'
      !CHARACTER (LEN=*), PARAMETER :: DWFR2 = '(2(A,I3),A,ES21.14E2)'
      !CHARACTER (LEN=*), PARAMETER :: DWFR3 = '(3(A,I3),A,ES21.14E2)'
      !CHARACTER (LEN=*), PARAMETER :: DWFR4 = '(4(A,I3),A,ES21.14E2)'
      !CHARACTER (LEN=*), PARAMETER :: DWFR5 = '(5(A,I3),A,ES21.14E2)'
      !CHARACTER (LEN=*), PARAMETER :: DWFR6 = '(6(A,I3),A,ES21.14E2)'
      !CHARACTER (LEN=*), PARAMETER :: DWFR7 = '(7(A,I3),A,ES21.14E2)'

      !CHARACTER (LEN=*), PARAMETER :: DWFR1_3 = '(A,I3,3(A,ES21.14E2))'

      CHARACTER (LEN=*), PARAMETER :: DWFC  = '(2A)'
      CHARACTER (LEN=*), PARAMETER :: DWFC1 = '(A,I3,2A)'
      CHARACTER (LEN=*), PARAMETER :: DWFC2 = '(2(A,I3),2A)'

!  Numbers

      REAL(mpk), PARAMETER :: &
        ZERO = 0.0_mpk, QUARTER = 0.25_mpk, &
        HALF = 0.5_mpk, ONE     = 1.0_mpk,  &
        TWO  = 2.0_mpk, THREE   = 3.0_mpk,  &
        FOUR = 4.0_mpk
      REAL(mpk), PARAMETER :: &
        MINUS_TWO = -TWO, MINUS_ONE = -ONE

!  End of file

      END MODULE TONGA_pars_m
