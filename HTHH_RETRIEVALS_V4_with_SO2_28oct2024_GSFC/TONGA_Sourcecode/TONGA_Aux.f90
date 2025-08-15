
!  File of auxiliary subroutines for TONGA retrieval.
!mick mod 9/12/2022 - created this module
!                   - added sub ASYMTX from VLIDORT (note: uses TONGA parameter E_MAXPARS)

      MODULE TONGA_aux_m

!  Parameter types

      USE TONGA_PARS_m, Only : E_MAXPARS, ZERO, ONE

      PRIVATE
      PUBLIC :: ASYMTX

      CONTAINS

      SUBROUTINE ASYMTX ( TOL, AAD, M, IA, IEVEC, &
                          EVECD, EVALD, IER, WKD, &
                          MESSAGE, BAD_STATUS )

!        AAD  :  input asymmetric matrix, destroyed after solved
!        M    :  order of  A
!       IA    :  first dimension of  A
!    IEVEC    :  first dimension of  EVECD

!   O U T P U T    V A R I A B L E S:

!       EVECD :  (unnormalized) eigenvectors of  A
!                   ( column J corresponds to EVALD(J) )

!       EVALD :  (unordered) eigenvalues of  A ( dimension at least M )

!       IER   :  if .NE. 0, signals that EVALD(IER) failed to converge;
!                   in that case eigenvalues IER+1,IER+2,...,M  are
!                   correct but eigenvalues 1,...,IER are set to zero.



!  1/31/21. Version 2.8.3. Add the tolerance variable

!    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

!       Solves eigenfunction problem for real asymmetric matrix
!       for which it is known a priori that the eigenvalues are real.

!       This is an adaptation of a subroutine EIGRF in the IMSL
!       library to use real instead of complex arithmetic, accounting
!       for the known fact that the eigenvalues and eigenvectors in
!       the discrete ordinate solution are real.  Other changes include
!       putting all the called subroutines in-line, deleting the
!       performance index calculation, updating many DO-loops
!       to Fortran77, and in calculating the machine precision
!       TOL instead of specifying it in a data statement.

!       EIGRF is based primarily on EISPACK routines.  The matrix is
!       first balanced using the parlett-reinsch algorithm.  Then
!       the Martin-Wilkinson algorithm is applied.

!       References:
!          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
!             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
!             Sources and Development of Mathematical Software,
!             Prentice-Hall, Englewood Cliffs, NJ
!         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
!             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
!         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
!             Clarendon Press, Oxford

!   I N P U T    V A R I A B L E S:

!        AAD  :  input asymmetric matrix, destroyed after solved
!        M    :  order of  A
!       IA    :  first dimension of  A
!    IEVEC    :  first dimension of  EVECD

!   O U T P U T    V A R I A B L E S:

!       EVECD :  (unnormalized) eigenvectors of  A
!                   ( column J corresponds to EVALD(J) )

!       EVALD :  (unordered) eigenvalues of  A ( dimension at least M )

!       IER   :  if .NE. 0, signals that EVALD(IER) failed to converge;
!                   in that case eigenvalues IER+1,IER+2,...,M  are
!                   correct but eigenvalues 1,...,IER are set to zero.

!  Include file for dimensioning

      implicit none

      LOGICAL            BAD_STATUS
      CHARACTER*(*)      MESSAGE

!   S C R A T C H   V A R I A B L E S:

!       WKD    :  WORK AREA ( DIMENSION AT LEAST 2*M )
!+---------------------------------------------------------------------+

!  input/output arguments
!    - 1/31/21. Version 2.8.3. Add the tolerance variable
!mick mod 9/12/2022 - replaced VLIDORT max variable "MAXSTRMSTKS"
!                     with Tonga max variable "E_MAXPARS"

      INTEGER              M, IA, IEVEC, IER
      DOUBLE PRECISION TOL, &
               AAD(E_MAXPARS,E_MAXPARS), WKD(4*E_MAXPARS), &
               EVALD(E_MAXPARS), EVECD(E_MAXPARS,E_MAXPARS)

!  local variables (explicit declaration

      LOGICAL           NOCONV, NOTLAS
      INTEGER              I, J, L, K, KKK, LLL
      INTEGER              N, N1, N2, IN, LB, KA, II
      DOUBLE PRECISION  C1, C2, C3, C4, C5, C6
      DATA              C1 / 0.4375D0 /
      DATA              C2 / 0.5D0 /
      DATA              C3 / 0.75D0 /
      DATA              C4 / 0.95D0 /
      DATA              C5 / 16.0D0 /
      DATA              C6 / 256.0D0 /
      DOUBLE PRECISION  DISCRI, SGN, RNORM, W, F, G, H, P, Q, R
      DOUBLE PRECISION  REPL, COL, ROW, SCALE, T, X, Z, S, Y, UU, VV

      IER = 0
      BAD_STATUS = .FALSE.
      MESSAGE = ' '

!  1/31/21. Version 2.8.3. Removed this section
!---4/18/08. Use 1.0d-12 (1.0d-24 has convergence issues in Rayleigh lay
!      TOL = 0.0000001
!      TOL = 1.0D-12
!      TOL = 1.0D-24
!      TOL = 1.0D-20

!       Here change to bypass D1MACH:
!        TOL = D1MACH(4)

      IF ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M ) THEN
        MESSAGE = 'ASYMTX--bad input variable(s)'
        BAD_STATUS = .TRUE.
        RETURN
      ENDIF
!                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES
      IF ( M.EQ.1 )  THEN
         EVALD(1) = AAD(1,1)
         EVECD(1,1) = 1.0D0
         RETURN
      ELSE IF ( M.EQ.2 )  THEN
         DISCRI = ( AAD(1,1) - AAD(2,2) )**2 + 4.0D0*AAD(1,2)*AAD(2,1)
         IF ( DISCRI.LT.ZERO ) THEN
           MESSAGE = 'ASYMTX--COMPLEX EVALS IN 2X2 CASE'
           BAD_STATUS = .TRUE.
           RETURN
         ENDIF
         SGN = ONE
         IF ( AAD(1,1).LT.AAD(2,2) )  SGN = - ONE
         EVALD(1) = 0.5D0*( AAD(1,1) + AAD(2,2) + SGN*DSQRT(DISCRI) )
         EVALD(2) = 0.5D0*( AAD(1,1) + AAD(2,2) - SGN*DSQRT(DISCRI) )
         EVECD(1,1) = ONE
         EVECD(2,2) = ONE
         IF ( AAD(1,1).EQ.AAD(2,2) .AND. &
               (AAD(2,1).EQ.ZERO.OR.AAD(1,2).EQ.ZERO) ) THEN
            RNORM = DABS(AAD(1,1))+DABS(AAD(1,2))+ &
                      DABS(AAD(2,1))+DABS(AAD(2,2))
            W = TOL * RNORM
            EVECD(2,1) = AAD(2,1) / W
            EVECD(1,2) = - AAD(1,2) / W
         ELSE
            EVECD(2,1) = AAD(2,1) / ( EVALD(1) - AAD(2,2) )
            EVECD(1,2) = AAD(1,2) / ( EVALD(2) - AAD(1,1) )
         ENDIF
         RETURN
      END IF
!                                        ** INITIALIZE OUTPUT VARIABLES
      DO 20 I = 1, M
         EVALD(I) = ZERO
         DO 10 J = 1, M
            EVECD(I,J) = ZERO
10       CONTINUE
         EVECD(I,I) = ONE
20    CONTINUE
!                  ** BALANCE THE INPUT MATRIX AND REDUCE ITS NORM BY
!                  ** DIAGONAL SIMILARITY TRANSFORMATION STORED IN WK;
!                  ** THEN SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!                  ** AND PUSH THEM DOWN
      RNORM = ZERO
      L  = 1
      K  = M

!mick check
!write(*,*)
!write(*,*) 'K (=M=MFIT) = ',K
!do J=1,K
!  do I=1,K
!    write(*,*) 'J = ',J,' I = ',I,' AAD(J,I) = ',AAD(J,I)
!  enddo
!enddo
!write(*,*) 'balancing the matrix ...'

30    KKK = K
         DO 70  J = KKK, 1, -1
            ROW = ZERO
            DO 40 I = 1, K

               !IF ( I.NE.J ) ROW = ROW + DABS( AAD(J,I) )
               IF ( I.NE.J ) THEN
                 !write(*,*) 'J = ',J,' I = ',I
                 !write(*,*) 'AAD(J,I) = ',AAD(J,I)
                 ROW = ROW + DABS( AAD(J,I) )
               ENDIF

40          CONTINUE
            IF ( ROW.EQ.ZERO ) THEN
               WKD(K) = J
               IF ( J.NE.K ) THEN
                  DO 50 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,K)
                     AAD(I,K) = REPL
50                CONTINUE
                  DO 60 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(K,I)
                     AAD(K,I) = REPL
60                CONTINUE
               END IF
               K = K - 1
               GO TO 30
            END IF
70       CONTINUE
!                                     ** SEARCH FOR COLUMNS ISOLATING AN
!                                       ** EIGENVALUE AND PUSH THEM LEFT
80    LLL = L
         DO 120 J = LLL, K
            COL = ZERO
            DO 90 I = L, K
               IF ( I.NE.J ) COL = COL + DABS( AAD(I,J) )
90          CONTINUE
            IF ( COL.EQ.ZERO ) THEN
               WKD(L) = J
               IF ( J.NE.L ) THEN
                  DO 100 I = 1, K
                     REPL   = AAD(I,J)
                     AAD(I,J) = AAD(I,L)
                     AAD(I,L) = REPL
100               CONTINUE
                  DO 110 I = L, M
                     REPL   = AAD(J,I)
                     AAD(J,I) = AAD(L,I)
                     AAD(L,I) = REPL
110               CONTINUE
               END IF
               L = L + 1
               GO TO 80
            END IF
120      CONTINUE
!                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH K
      DO 130 I = L, K
         WKD(I) = ONE
130   CONTINUE

140   NOCONV = .FALSE.
         DO 200 I = L, K
            COL = ZERO
            ROW = ZERO
            DO 150 J = L, K
               IF ( J.NE.I ) THEN
                  COL = COL + DABS( AAD(J,I) )
                  ROW = ROW + DABS( AAD(I,J) )
               END IF
150         CONTINUE
            F = ONE
            G = ROW / C5
            H = COL + ROW
160         IF ( COL.LT.G ) THEN
               F   = F * C5
               COL = COL * C6
               GO TO 160
            END IF
            G = ROW * C5
170         IF ( COL.GE.G ) THEN
               F   = F / C5
               COL = COL / C6
               GO TO 170
            END IF
!                                                         ** NOW BALANCE
            IF ( (COL+ROW)/F .LT. C4*H ) THEN
               WKD(I)  = WKD(I) * F
               NOCONV = .TRUE.
               DO 180 J = L, M
                  AAD(I,J) = AAD(I,J) / F
180            CONTINUE
               DO 190 J = 1, K
                  AAD(J,I) = AAD(J,I) * F
190            CONTINUE
            END IF
200      CONTINUE

      IF ( NOCONV ) GO TO 140
!                                  ** IS -A- ALREADY IN HESSENBERG FORM?
      IF ( K-1 .LT. L+1 ) GO TO 350
!                                   ** TRANSFER -A- TO A HESSENBERG FORM
      DO 290 N = L+1, K-1
         H        = ZERO
         WKD(N+M) = ZERO
         SCALE    = ZERO
!                                                        ** SCALE COLUMN
         DO 210 I = N, K
            SCALE = SCALE + DABS(AAD(I,N-1))
210      CONTINUE
         IF ( SCALE.NE.ZERO ) THEN
            DO 220 I = K, N, -1
               WKD(I+M) = AAD(I,N-1) / SCALE
               H = H + WKD(I+M)**2
220         CONTINUE
            G = - SIGN( DSQRT(H), WKD(N+M) )
            H = H - WKD(N+M) * G
            WKD(N+M) = WKD(N+M) - G
!                                                 ** FORM (I-(U*UT)/H)*A
            DO 250 J = N, M
               F = ZERO
               DO 230  I = K, N, -1
                  F = F + WKD(I+M) * AAD(I,J)
230            CONTINUE
               DO 240 I = N, K
                  AAD(I,J) = AAD(I,J) - WKD(I+M) * F / H
240            CONTINUE
250         CONTINUE
!                                    ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 280 I = 1, K
               F = ZERO
               DO 260  J = K, N, -1
                  F = F + WKD(J+M) * AAD(I,J)
260            CONTINUE
               DO 270 J = N, K
                  AAD(I,J) = AAD(I,J) - WKD(J+M) * F / H
270            CONTINUE
280         CONTINUE
            WKD(N+M)  = SCALE * WKD(N+M)
            AAD(N,N-1) = SCALE * G
         END IF
290   CONTINUE

      DO 340  N = K-2, L, -1
         N1 = N + 1
         N2 = N + 2
         F  = AAD(N1,N)
         IF ( F.NE.ZERO ) THEN
            F  = F * WKD(N1+M)
            DO 300 I = N2, K
               WKD(I+M) = AAD(I,N)
300         CONTINUE
            IF ( N1.LE.K ) THEN
               DO 330 J = 1, M
                  G = ZERO
                  DO 310 I = N1, K
                     G = G + WKD(I+M) * EVECD(I,J)
310               CONTINUE
                  G = G / F
                  DO 320 I = N1, K
                     EVECD(I,J) = EVECD(I,J) + G * WKD(I+M)
320               CONTINUE
330            CONTINUE
            END IF
         END IF
340   CONTINUE

350   CONTINUE
      N = 1
      DO 370 I = 1, M
         DO 360 J = N, M
            RNORM = RNORM + DABS(AAD(I,J))
360      CONTINUE
         N = I
         IF ( I.LT.L .OR. I.GT.K ) EVALD(I) = AAD(I,I)
370   CONTINUE
      N = K
      T = ZERO
!                                         ** SEARCH FOR NEXT EIGENVALUES
380   IF ( N.LT.L ) GO TO 530
      IN = 0
      N1 = N - 1
      N2 = N - 2
!                          ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
390   CONTINUE
      DO 400 I = L, N
         LB = N+L - I
         IF ( LB.EQ.L ) GO TO 410
         S = DABS( AAD(LB-1,LB-1) ) + DABS( AAD(LB,LB) )
         IF ( S.EQ.ZERO ) S = RNORM
         IF ( DABS(AAD(LB,LB-1)) .LE. TOL*S ) GO TO 410
400   CONTINUE

410   X = AAD(N,N)
      IF ( LB.EQ.N ) THEN
!                                        ** ONE EIGENVALUE FOUND
         AAD(N,N)  = X + T
         EVALD(N) = AAD(N,N)
         N = N1
         GO TO 380
      END IF

      Y = AAD(N1,N1)
      W = AAD(N,N1) * AAD(N1,N)
      IF ( LB.EQ.N1 ) THEN
!                                        ** TWO EIGENVALUES FOUND
         P = (Y-X) * C2
         Q = P**2 + W
         Z = DSQRT( DABS(Q) )
         AAD(N,N) = X + T
         X = AAD(N,N)
         AAD(N1,N1) = Y + T
!                                        ** REAL PAIR
         Z = P + SIGN(Z,P)
         EVALD(N1) = X + Z
         EVALD(N)  = EVALD(N1)
         IF ( Z.NE.ZERO ) EVALD(N) = X - W / Z
         X = AAD(N,N1)
!                                  ** EMPLOY SCALE FACTOR IN CASE
!                                  ** X AND Z ARE VERY SMALL
         R = SQRT( X*X + Z*Z )
         P = X / R
         Q = Z / R
!                                             ** ROW MODIFICATION
         DO 420 J = N1, M
            Z = AAD(N1,J)
            AAD(N1,J) = Q * Z + P * AAD(N,J)
            AAD(N,J)  = Q * AAD(N,J) - P * Z
420      CONTINUE
!                                             ** COLUMN MODIFICATION
         DO 430 I = 1, N
            Z = AAD(I,N1)
            AAD(I,N1) = Q * Z + P * AAD(I,N)
            AAD(I,N)  = Q * AAD(I,N) - P * Z
430      CONTINUE
!                                          ** ACCUMULATE TRANSFORMATIONS
         DO 440 I = L, K
            Z = EVECD(I,N1)
            EVECD(I,N1) = Q * Z + P * EVECD(I,N)
            EVECD(I,N)  = Q * EVECD(I,N) - P * Z
440      CONTINUE

         N = N2
         GO TO 380
      END IF

      IF ( IN.EQ.30 ) THEN
!                    ** NO CONVERGENCE AFTER 30 ITERATIONS; SET ERROR
!                    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUE
         IER = N
         GO TO 670
      END IF
!                                                          ** FORM SHIFT
      IF ( IN.EQ.10 .OR. IN.EQ.20 ) THEN
         T = T + X
         DO 450 I = L, N
            AAD(I,I) = AAD(I,I) - X
450      CONTINUE
         S = DABS(AAD(N,N1)) + DABS(AAD(N1,N2))
         X = C3 * S
         Y = X
         W = - C1 * S**2
      END IF

      IN = IN + 1
!                ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS

      DO 460 J = LB, N2
         I = N2+LB - J
         Z = AAD(I,I)
         R = X - Z
         S = Y - Z
         P = ( R * S - W ) / AAD(I+1,I) + AAD(I,I+1)
         Q = AAD(I+1,I+1) - Z - R - S
         R = AAD(I+2,I+1)
         S = DABS(P) + DABS(Q) + DABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF ( I.EQ.LB ) GO TO 470
         UU = DABS( AAD(I,I-1) ) * ( DABS(Q) + DABS(R) )
         VV = DABS(P)*(DABS(AAD(I-1,I-1))+DABS(Z)+DABS(AAD(I+1,I+1)))
         IF ( UU .LE. TOL*VV ) GO TO 470
460   CONTINUE

470   CONTINUE
      AAD(I+2,I) = ZERO
      DO 480 J = I+3, N
         AAD(J,J-2) = ZERO
         AAD(J,J-3) = ZERO
480   CONTINUE

!             ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO N

      DO 520 KA = I, N1
         NOTLAS = KA.NE.N1
         IF ( KA.EQ.I ) THEN
            S = SIGN( DSQRT( P*P + Q*Q + R*R ), P )
            IF ( LB.NE.I ) AAD(KA,KA-1) = - AAD(KA,KA-1)
         ELSE
            P = AAD(KA,KA-1)
            Q = AAD(KA+1,KA-1)
            R = ZERO
            IF ( NOTLAS ) R = AAD(KA+2,KA-1)
            X = DABS(P) + DABS(Q) + DABS(R)
            IF ( X.EQ.ZERO ) GO TO 520
            P = P / X
            Q = Q / X
            R = R / X
            S = SIGN( DSQRT( P*P + Q*Q + R*R ), P )
            AAD(KA,KA-1) = - S * X
         END IF
         P = P + S
         X = P / S
         Y = Q / S
         Z = R / S
         Q = Q / P
         R = R / P
!                                                    ** ROW MODIFICATION
         DO 490 J = KA, M
            P = AAD(KA,J) + Q * AAD(KA+1,J)
            IF ( NOTLAS ) THEN
               P = P + R * AAD(KA+2,J)
               AAD(KA+2,J) = AAD(KA+2,J) - P * Z
            END IF
            AAD(KA+1,J) = AAD(KA+1,J) - P * Y
            AAD(KA,J)   = AAD(KA,J)   - P * X
490      CONTINUE
!                                                 ** COLUMN MODIFICATION
         DO 500 II = 1, MIN0(N,KA+3)
            P = X * AAD(II,KA) + Y * AAD(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * AAD(II,KA+2)
               AAD(II,KA+2) = AAD(II,KA+2) - P * R
            END IF
            AAD(II,KA+1) = AAD(II,KA+1) - P * Q
            AAD(II,KA)   = AAD(II,KA) - P
500      CONTINUE
!                                          ** ACCUMULATE TRANSFORMATIONS
         DO 510 II = L, K
            P = X * EVECD(II,KA) + Y * EVECD(II,KA+1)
            IF ( NOTLAS ) THEN
               P = P + Z * EVECD(II,KA+2)
               EVECD(II,KA+2) = EVECD(II,KA+2) - P * R
            END IF
            EVECD(II,KA+1) = EVECD(II,KA+1) - P * Q
            EVECD(II,KA)   = EVECD(II,KA) - P
510      CONTINUE

520   CONTINUE
      GO TO 390
!                     ** ALL EVALS FOUND, NOW BACKSUBSTITUTE REAL VECTOR
530   CONTINUE
      IF ( RNORM.NE.ZERO ) THEN
         DO 560  N = M, 1, -1
            N2 = N
            AAD(N,N) = ONE
            DO 550  I = N-1, 1, -1
               W = AAD(I,I) - EVALD(N)
               IF ( W.EQ.ZERO ) W = TOL * RNORM
               R = AAD(I,N)
               DO 540 J = N2, N-1
                  R = R + AAD(I,J) * AAD(J,N)
540            CONTINUE
               AAD(I,N) = - R / W
               N2 = I
550         CONTINUE
560      CONTINUE
!                      ** END BACKSUBSTITUTION VECTORS OF ISOLATED EVALS

         DO 580 I = 1, M
            IF ( I.LT.L .OR. I.GT.K ) THEN
               DO 570 J = I, M
                  EVECD(I,J) = AAD(I,J)
570            CONTINUE
            END IF
580      CONTINUE
!                                   ** MULTIPLY BY TRANSFORMATION MATRIX
         IF ( K.NE.0 ) THEN
            DO 600  J = M, L, -1
               DO 600 I = L, K
                  Z = ZERO
                  DO 590 N = L, MIN0(J,K)
                     Z = Z + EVECD(I,N) * AAD(N,J)
590               CONTINUE
                  EVECD(I,J) = Z
600         CONTINUE
         END IF

      END IF

      DO 620 I = L, K
         DO 620 J = 1, M
            EVECD(I,J) = EVECD(I,J) * WKD(I)
620   CONTINUE
!                           ** INTERCHANGE ROWS IF PERMUTATIONS OCCURRED
      DO 640  I = L-1, 1, -1
         J = INT(WKD(I))
         IF ( I.NE.J ) THEN
            DO 630 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
630         CONTINUE
         END IF
640   CONTINUE

      DO 660 I = K+1, M
         J = INT(WKD(I))
         IF ( I.NE.J ) THEN
            DO 650 N = 1, M
               REPL       = EVECD(I,N)
               EVECD(I,N) = EVECD(J,N)
               EVECD(J,N) = REPL
650         CONTINUE
         END IF
660   CONTINUE
!
  670 CONTINUE

      RETURN
      END SUBROUTINE  ASYMTX

!  End module

      END MODULE TONGA_aux_m
