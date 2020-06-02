module m_tbessy
!************************************************************************
!*                                                                      *
!*    Program to calculate the second kind Bessel function of integer   *
!*    order N, for any REAL X, using the function BESSY(N,X).           *
!*                                                                      *
!* -------------------------------------------------------------------- *
!*                                                                      *
!*    SAMPLE RUN:                                                       *
!*                                                                      *
!*    (Calculate Bessel function for N=2, X=0.75).                      *
!*                                                                      *
!*    Second kind Bessel Function of order  2 for X =  0.7500:          *
!*                                                                      *
!*         Y = -0.26297459E+01                                          *
!*                                                                      *
!* -------------------------------------------------------------------- *
!*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
!*                                                                      *
!*                               F90 Release 1.0 By J-P Moreau, Paris.  *
!*                                     all variables declared           *
!*                                        (www.jpmoreau.fr)             *
!************************************************************************
private
public :: BESSY

contains

real(kind=8) FUNCTION BESSY(N,X)
! ------------------------------------------------------------------
!     This subroutine calculates the second kind Bessel Function of
!     integer order N, for any real X. We use here the classical
!     recursive formula.
! ------------------------------------------------------------------
IMPLICIT NONE
INTEGER N, J
real(kind=8) :: X,TOX,BY,BYM,BYP

IF (N.EQ.0) THEN
   BESSY = BESSY0(X)
   RETURN
ENDIF

IF (N.EQ.1) THEN
   BESSY = BESSY1(X)
   RETURN
ENDIF

IF (X.EQ.0.) THEN
   BESSY = -1.E30
   RETURN
ENDIF

TOX = 2./X
BY  = BESSY1(X)
BYM = BESSY0(X)
DO J = 1,N-1
   BYP = J*TOX*BY-BYM
   BYM = BY
   BY  = BYP
enddo
BESSY = BY

end function

! ---------------------------------------------------------------------------
real(kind=8) FUNCTION BESSYP (N,X)
IMPLICIT NONE
INTEGER N
real(kind=8) ::  X,BESSY

BESSYP=0.0d0
IF (N.EQ.0) THEN
   BESSYP=-BESSY(1,X)
ELSE IF(X.EQ.0.D0) THEN
   X=1.D-30
ELSE
   BESSYP=BESSY(N-1,X)-(DFLOAT(N)/X)*BESSY(N,X)
ENDIF

END function

! ---------------------------------------------------------------------------
real(kind=8) FUNCTION BESSY0 (X)
IMPLICIT NONE
real(kind=8) :: X,FS,FR,Z,FP,FQ,XX
! ---------------------------------------------------------------------
!     This subroutine calculates the Second Kind Bessel Function of
!     order 0, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
! ---------------------------------------------------------------------
real(kind=8) :: Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3,  &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
DATA R1,R2,R3,R4,R5,R6 /-2957821389.D0,7062834065.D0, &
      -512359803.6D0,10879881.29D0,-86327.92757D0,228.4622733D0 /
DATA S1,S2,S3,S4,S5,S6 /40076544269.D0,745249964.8D0, &
      7189466.438D0,47447.26470D0,226.1030244D0,1.D0 /

BESSY0 = -1.E30
IF (X.EQ.0.D0) return

IF (X.LT.8.D0) THEN
   Y = X*X
   FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
   FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
   BESSY0 = FR/FS+.636619772D0*BESSJ0(X)*LOG(X)
ELSE
   Z = 8.D0/X
   Y = Z*Z
   XX = X-.785398164D0
   FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
   FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
   BESSY0 = SQRT(.636619772D0/X)*(FP*SIN(XX)+Z*FQ*COS(XX))
ENDIF

END function

! ---------------------------------------------------------------------------
real(kind=8) FUNCTION BESSY1 (X)
IMPLICIT NONE
real(kind=8) :: X,FR,FS,Z,FP,FQ,XX
! ----------------------------------------------------------------------
!     This subroutine calculates the Second Kind Bessel Function of
!     order 1, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
! ----------------------------------------------------------------------
real(kind=8) :: Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6,S7
DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4, &
      .2457520174D-5,-.240337019D-6 /
DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,  &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
DATA R1,R2,R3,R4,R5,R6 /-.4900604943D13,.1275274390D13,   &
      -.5153438139D11,.7349264551D9,-.4237922726D7,.8511937935D4 /
DATA S1,S2,S3,S4,S5,S6,S7 /.2499580570D14,.4244419664D12, &
      .3733650367D10,.2245904002D8,.1020426050D6,.3549632885D3,1.D0 /

BESSY1 = -1.E30
IF (X.EQ.0.) return

IF (X.LT.8.) THEN
   Y = X*X
   FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
   FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*(S6+Y*S7)))))
   BESSY1 = X*(FR/FS)+.636619772*(BESSJ1(X)*LOG(X)-1./X)
ELSE
   Z = 8./X
   Y = Z*Z
   XX = X-2.356194491
   FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
   FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
   BESSY1 = SQRT(.636619772/X)*(SIN(XX)*FP+Z*COS(XX)*FQ)
ENDIF

END function

! ----------------------------------------------------------------------
real(kind=8) FUNCTION BESSJ0(X)
IMPLICIT NONE
real(kind=8) :: X,AX,FR,FS,Z,FP,FQ,XX
! --------------------------------------------------------------------
!     This subroutine calculates the First Kind Bessel Function of
!     order 0, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
! ---------------------------------------------------------------------
real(kind=8) :: Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /

BESSJ0 = 1.D0
IF(X.EQ.0.D0) return

AX = ABS (X)
IF (AX.LT.8.) THEN
   Y = X*X
   FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
   FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
   BESSJ0 = FR/FS
ELSE
   Z = 8./AX
   Y = Z*Z
   XX = AX-.785398164
   FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
   FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
   BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
ENDIF

END function


! ---------------------------------------------------------------------------
real(kind=8) FUNCTION BESSJ1(X)
IMPLICIT NONE
real(kind=8) ::  X,AX,FR,FS,Z,FP,FQ,XX
! ---------------------------------------------------------------------
!     This subroutine calculates the First Kind Bessel Function of
!     order 1, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
! ---------------------------------------------------------------------
real(kind=8) :: Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, &
      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

BESSJ1=0.0d0
AX = ABS(X)
IF (AX.LT.8.) THEN
   Y = X*X
   FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
   FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
   BESSJ1 = X*(FR/FS)
ELSE
   Z = 8./AX
   Y = Z*Z
   XX = AX-2.35619491
   FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
   FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
   BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
ENDIF

END function


end module m_tbessy
