module m_invert
use m_determinant, only : determinant
use m_constants, only : identity

interface invert
  module procedure :: invert_real,invert_complex
end interface invert

private
public :: invert

contains

subroutine invert_real(in,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
implicit none
integer, intent(in) :: n
real(8), intent(in) :: in(:,:)
real(8), intent(out) :: c(:,:)
! internal
real(kind=8) :: L(n,n), U(n,n), b(n), d(n), x(n), a(n,n)
real(kind=8) :: coeff
integer :: i, j, k

! step -1 try the transpose as the inverse
c=transpose(in)
a=matmul(c,in)
if (all(abs(a-identity(1.0d0,n)).lt.1.0d-8)) return

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0
a=in
if (determinant(1.0d-8,n,in).lt.1.0d-8) STOP 'DET=0 matrix found and can not be inverted'
! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do

end subroutine invert_real

!***********************************************
!* SOLVING A COMPLEX LINEAR MATRIX SYSTEM AX=B *
!* with Gauss-Jordan method using full pivoting*
!* at each step. During the process, original  *
!* A and B matrices are destroyed to spare     *
!* storage location.                           *
!* ------------------------------------------- *
!* INPUTS:    A    COMPLEX MATRIX N*N          *
!*            B    COMPLEX MATRIX N*M          *
!* ------------------------------------------- *
!* OUTPUTS:   A    INVERSE OF A N*N            *
!*            DET  COMPLEX DETERMINANT OF A    *
!*            B    SOLUTION MATRIX N*M         *
!* ------------------------------------------- *
!* NOTA - If M=0 inversion of A matrix only.   *
!***********************************************
SUBROUTINE invert_complex(N,M,AA,BB,DET)
implicit none
integer, intent(in) :: N,M
COMPLEX(kind=8), intent(inout) :: AA(:,:),BB(:,:)    !AA(N,N),BB(N,M)
COMPLEX(kind=8), intent(out) :: DET
! internal variables
REAL(kind=8), PARAMETER :: EPSMACH=2.E-12
INTEGER, POINTER :: PC(:), PL(:)
COMPLEX(kind=8), POINTER :: CS(:)
COMPLEX(kind=8) :: PV, TT
REAL(kind=8) :: PAV
integer :: I, IK, J, JK, K, ialloc

!Initializations :
allocate(PC(1:N),stat=ialloc)
allocate(PL(1:N),stat=ialloc)
allocate(CS(1:N),stat=ialloc)
DET=CMPLX(1.0,0.0)
DO I=1,N
   PC(I)=0
   PL(I)=0
   CS(I)=CMPLX(0.0,0.0)
END DO
!main loop
DO K=1,N
!Searching greatest pivot:
   PV=AA(K,K)
   IK=K
   JK=K
   PAV=ABS(PV)
   DO I=K,N
      DO J=K,N
         IF (ABS(AA(I,J)).GT.PAV) THEN
            PV=AA(I,J)
            PAV=ABS(PV)
            IK=I
            JK=J
         ENDIF
      ENDDO
   ENDDO

!Search terminated, the pivot is in location I=IK, J=JK.
!Memorizing pivot location: :
   PC(K)=JK
   PL(K)=IK

!Determinant DET is actualised
!If DET=0, ERROR MESSAGE and STOP
!Machine dependant EPSMACH equals here 2.E-12

   IF (IK.NE.K) DET=-DET
   IF (JK.NE.K) DET=-DET
   DET=DET*PV
   !Error message and Stop
   IF (ABS(DET).LT.EPSMACH) then
     write(6,'(f12.8)') ABS(DET)
     stop '  DETERMINANT EQUALS ZERO, NO SOLUTION!'
   endif


!POSITIONNING PIVOT IN K,K:
   IF(IK.NE.K) THEN
      DO I=1,N
!EXCHANGE LINES IK and K:
         TT=AA(IK,I)
         AA(IK,I)=AA(K,I)
         AA(K,I)=TT
      ENDDO
   ENDIF
   IF (M.NE.0) THEN
      DO I=1,M
         TT=BB(IK,I)
         BB(IK,I)=BB(K,I)
         BB(K,I)=TT
      ENDDO
   ENDIF
!Pivot is at correct line
   IF(JK.NE.K) THEN
      DO I=1,N
!Exchange columns JK and K of matrix AA
         TT=AA(I,JK)
         AA(I,JK)=AA(I,K)
         AA(I,K)=TT
      ENDDO
   ENDIF
!Pivot is at correct column and located in K,K

!Store column K in vector CS then set column K to zero
   DO I=1,N
      CS(I)=AA(I,K)
      AA(I,K)=CMPLX(0.0,0.0)
   ENDDO
!
   CS(K)=CMPLX(0.0,0.0)
   AA(K,K)=CMPLX(1.0,0.0)
!Modify line K :
   IF(ABS(PV).LT.EPSMACH) STOP '  PIVOT TOO SMALL - STOP'
   DO I=1,N
      AA(K,I)=AA(K,I)/PV
   ENDDO
   IF (M.NE.0) THEN
      DO I=1,M
         BB(K,I)=BB(K,I)/PV
      ENDDO
   ENDIF
!Modify other lines of matrix AA:
   DO J=1,N
      IF (J.EQ.K) cycle
         DO I=1,N
!Modify line J of matrix AA :
            AA(J,I)=AA(J,I)-CS(J)*AA(K,I)
         ENDDO
      IF (M.NE.0) THEN
         DO I=1,M
            BB(J,I)=BB(J,I)-CS(J)*BB(K,I)
         ENDDO
      ENDIF
   ENDDO
!Line K is ready.
ENDDO
!End of K loop

!The matrix AA is inverted - Rearrange AA

!Exchange lines
   DO I=N,1,-1
      IK=PC(I)
      IF (IK.EQ.I) cycle
!EXCHANGE LINES I AND PC(I) OF AA:
      DO J=1,N
         TT=AA(I,J)
         AA(I,J)=AA(IK,J)
         AA(IK,J)=TT
      ENDDO
      IF (M.NE.0) THEN
         DO J=1,M
           TT=BB(I,J)
           BB(I,J)=BB(IK,J)
           BB(IK,J)=TT
         ENDDO
      ENDIF
!NO MORE EXCHANGE NEEDED
!GO TO NEXT LINE
   ENDDO

!EXCHANGE COLUMNS

   DO J=N,1,-1
     JK=PL(J)
     IF (JK.EQ.J) cycle
!EXCHANGE COLUMNS J AND PL(J) OF AA :
     DO I=1,N
        TT=AA(I,J)
        AA(I,J)=AA(I,JK)
        AA(I,JK)=TT
     ENDDO
!NO MORE EXCHANGE NEEDED
!GO TO NEXT COLUMN
   ENDDO


END subroutine invert_complex


end module m_invert
