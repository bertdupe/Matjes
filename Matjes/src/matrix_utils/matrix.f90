module m_matrix

!
! contains utilities for basic matrix operations
!
!
!

interface invert
  module procedure :: invert_real,invert_complex
end interface invert

interface reduce
  module procedure :: reduce_real,reduce_op,reduce_Nreal
end interface reduce

interface Jacobi
  module procedure :: eigenvalvec_2D
end interface

interface determinant
  module procedure :: det_real,det_complex
end interface determinant

interface TSGT
  module procedure :: TSRGT,TSCGT
end interface TSGT

private
public :: invert,reduce,Jacobi,determinant

contains


!
! take a matric of rank 2 and size (n,n**2) and a matrix of size (n) and return
! a matric a vector of size n. Used typically for the Hamiltonian of rank N > 2
!
subroutine reduce_real(matrix_in,mode_out,mode,n)
implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: matrix_in(:,:),mode(:)
real(kind=8), intent(inout) :: mode_out(:)
! internal functions
integer :: i,shape_in(2),j,k,test,ncolumn,nline
real(kind=8) :: matrix_int(n,n**3),S_int(n)

matrix_int=0.0d0
shape_in=shape(matrix_in)
matrix_int(:,1:shape_in(2))=matrix_in

ncolumn=shape_in(1)
nline=shape_in(2)
S_int=mode(1:n)

test=nline/ncolumn
do while (test.ne.0)

  do i=1,nline
    j=mod(i-1,ncolumn)+1
    k=(i-1)/ncolumn+1
    matrix_int(j,k)=dot_product( matrix_int(:,i) , S_int )
  enddo

  nline=test
  test=nline/ncolumn
enddo

mode_out=matrix_int(:,1)

end subroutine

!
! take a matrix of rank 2 and size (n,n**2) and a matrix of size (n,N) and return
! a matric a vector of size n. Used typically for the Hamiltonian of rank N > 2
!
subroutine reduce_Nreal(matrix_in,mode_out,mode,n)
implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: matrix_in(:,:),mode(:,:)
real(kind=8), intent(out) :: mode_out(:)
! internal functions
integer :: i,shape_in(2),j,k,test,ncolumn,nline,i_loop
real(kind=8) :: matrix_int(n,n**3),S_int(n)

matrix_int=0.0d0
shape_in=shape(matrix_in)
matrix_int(:,1:shape_in(2))=matrix_in

ncolumn=shape_in(1)
nline=shape_in(2)

test=nline/ncolumn
i_loop=0

do while (test.ne.0)

  i_loop=i_loop+1
  S_int=mode(:,i_loop)

  do i=1,nline
    j=mod(i-1,ncolumn)+1
    k=(i-1)/ncolumn+1
    matrix_int(j,k)=dot_product( matrix_int(:,i) , S_int )
  enddo

  nline=test
  test=nline/ncolumn

enddo

mode_out=matrix_int(:,1)

end subroutine

!
! take a table of operater matrices of rank P and size (n,n**P) and a matrix of size (n) and return
! a matric a vector of size n. Used typically for the Hamiltonian of rank N > 2
!
! energy%value(i,iomp),energy%value(i,iomp)%num,S_int,spin(iomp)%w,spin(j)%w
subroutine reduce_op(matrix_in,n_order,mode_out,mode_1,mode_2,dim_ham)
use m_basic_types, only : Op_real_order_N
implicit none
integer, intent(in) :: dim_ham,n_order
real(kind=8), intent(in) :: mode_1(:),mode_2(:)
type(Op_real_order_N), intent(in) :: matrix_in
real(kind=8), intent(inout) :: mode_out(:)
! internal functions
integer :: i
real(kind=8) :: big_mode(dim_ham,n_order)
real(kind=8) :: mode_int(dim_ham)

mode_out=0.0d0
big_mode(:,1)=mode_2

call reduce_real(matrix_in%order_op(1)%Op_loc,mode_out,big_mode(:,1),dim_ham)

big_mode(:,1)=mode_1
do i=2,n_order
   big_mode(:,2)=mode_2
   call reduce_Nreal(matrix_in%order_op(2)%Op_loc,mode_int,big_mode,dim_ham)
   mode_out=mode_out+mode_int
enddo

end subroutine









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
real(kind=8), intent(in) :: in(:,:)
real(kind=8), intent(out) :: c(:,:)
! internal
real(kind=8) :: L(n,n), U(n,n), b(n), d(n), x(n), a(n,n)
real(kind=8) :: coeff
integer :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0
a=in
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
COMPLEX(kind=16), intent(inout) :: AA(:,:),BB(:,:)    !AA(N,N),BB(N,M)
COMPLEX(kind=16), intent(out) :: DET
! internal variables
REAL(kind=8), PARAMETER :: EPSMACH=2.E-12
INTEGER, POINTER :: PC(:), PL(:)
COMPLEX(kind=16), POINTER :: CS(:)
COMPLEX(kind=16) :: PV, TT
REAL(kind=16) :: PAV
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







subroutine eigenvalvec_2D(a,x,abserr,n)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x
! method: Jacoby method for symmetric matrices
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: abserr
real(kind=8), intent(inout) :: a(:,:),x(:,:)
! internal
integer :: i, j, k
real(kind=8) :: b2, bar
real(kind=8) :: beta, coeff, c, s, cs, sc

! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
x = 0.0
do i=1,n
  x(i,i) = 1.0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.0
do i=1,n
  do j=1,n
    if (i.ne.j) b2 = b2 + a(i,j)**2
  end do
end do

if (b2 <= abserr) return

! average for off-diagonal elements /2
bar = 0.5*b2/float(n*n)

do while (b2.gt.abserr)
  do i=1,n-1
    do j=i+1,n
      if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
      b2 = b2 - 2.0*a(j,i)**2
      bar = 0.5*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
      coeff = 0.5*beta/sqrt(1.0+beta**2)
      s = sqrt(max(0.5+coeff,0.0))
      c = sqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
      do k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      end do
! new matrix a_{k+1} from a_{k}, and eigenvectors
      do k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      end do
    end do
  end do
end do

end subroutine






!The function DMGT returns the determinant of a real square matrix
!A(n,n) by Gauss method with full pivoting.
!----------------------------------------------------------------------------
!  Input parameters:
!  eps        precision (real*8)
!  n          size of A matrix (integer)
!  A          pointer to input real square matrix
!  Output parameters:
!  None
!-----------------------------------------------------------------------------
!The procedure TSRGT is used to reduce A matrix to an upper triangular matrix.
!Output variables are it(integer), C(n,n), Kp(n) and Lp(n).
!If it=0, matrix A is singular, if it=1, matrix A is regular. Table C contains
!at location i,j (j>=i) the corresponding element of the upper triangular matrix.
!Tables Lp and Kp contain informations relative to exchanges of line or column
!that occured during the process. For instance, the element number k of Lp is
!an integer <> k if an exchange of line has been made at step k (k=1,2,...,n).
!The number of exchanges of lines and columns is stored in integer L. the
!determinant of A matrix is stored in d0 (real*8).
!-----------------------------------------------------------------------------
Function det_real(eps, n, A) result(DMGT)
implicit none
integer :: n
real(kind=8) :: eps, A(:,:),DMGT,d0

real(kind=8), pointer :: C(:,:)
integer,pointer :: Kp(:), Lp(:)
integer :: it,k,l,ialloc

!allocate local matrix C and vectors Kp, Lp
  allocate(C(n,n),STAT=ialloc)
  allocate(Kp(n),STAT=ialloc)
  allocate(Lp(n),STAT=ialloc)

  call TSGT(eps,n,A,it,C,Kp,Lp)  !call triangularization subroutine

  if (it==0) write(6,'(a)') 'WARNING: Values are very small in the matrix'

  d0=1.d0
  do k=1, n
     d0=d0*C(k,k)
  end do
  l=0
  do k=1, n-1
     if (Lp(k).ne.k)  l=l+1
     if (Kp(k).ne.k)  l=l+1
  end do
  if (MOD(l,2).ne.0) d0=-d0  !l is odd

  DMGT=d0   !return determinant

End



!****************************************************************
 !* The DCGT procedure calculates the complex determinant of a   *
 !* complex square matrix by the Gauss method with full pivoting *
 !* ------------------------------------------------------------ *
 !* INPUTS:                                                      *
 !*         eps:  required precision                             *
 !*          N :  size of matrix A                               *
 !*          A :  complex matrix of size N x N                   *
 !* OUTPUT:                                                      *
 !*         det:  complex determinant.                           *
 !****************************************************************
  function det_complex(eps, N, A) result(det)
    implicit none
    real(kind=8) :: eps
    integer :: N
    Complex(kind=16) :: A(:,:), det
    Complex(kind=16) :: C0,Z1
    Complex(kind=16) :: C(N,N)
    integer :: KP(N), LP(N),it,k,l

    Z1 = 1.

    call TSGT(eps,N,A,it,C,KP,LP)

    if (it==0) write(6,'(a)') 'WARNING: Values are very small in the matrix'

    det=Z1
    do k=1, N
        C0=det
        det = C0 * C(k,k)
    end do
    l=0
    do K=1, N-1
       if (LP(k).ne.k) l=l+1
       if (KP(k).ne.k) l=l+1
    end do
    if (Mod(l,2).ne.0)  det = -det
  end  !DCGT

!The subroutine TSRGT applies to input real square matrix A(n,n) the upper
!triangularization algorithm of Gauss method with full pivoting and keeps
!trace of successive transformations done in integer vectors KP and LP.
!-----------------------------------------------------------------------------
!  Input parameters:
!  eps        precision (real*8)
!  n          size of A matrix (integer)
!  A          pointer to input real square matrix (real*8)
!  Output parameters:
!  it         flag=1 if A matrix ok, =0 if A matrix is singular (integer)
!  C          pointer to table storing main diagonal elements and supra-
!             diagonal elements of upper triangular matrix and the multi-
!             plying coefficients used during triangularization process
!  KP         table storing informations concerning the column exchanges
!             during process (integer)
!  LP         table storing informations concerning the line exchanges
!             during process (integer)
!-----------------------------------------------------------------------------
!The table C is first initialized to A matrix, then receives at each step k
!of the triangularization process, usefull elements of A matrix at step k for
!k=1,2,...n.
!The variables po(real*8), lo and ko(integer) store respectively pivot at step k,
!its line number and its column number.
!------------------------------------------------------------------------------
Subroutine TSRGT(eps, n, A, it, C, Kp, Lp)
  implicit none
  real(kind=8), intent(in) :: eps
  integer, intent(in) :: n
  real(kind=8), intent(in) :: A(:,:)
  real(kind=8), intent(out) :: C(:,:)
  integer, intent(out) :: Kp(:),Lp(:),it
  ! internal
  real(kind=8) ::  po,t0
  integer :: k,i,j,lo,ko

  C=A; it=1; k=1
  do while (it==1.and.k<n)
    po=C(k,k); lo=k; ko=k
    do i=k, n
      do j=k, n
        if (dabs(C(i,j))>dabs(po)) then
          po=C(i,j); lo=i; ko=j
        end if
      end do
    end do
    Lp(k)=lo; Kp(k)=ko
    if (dabs(po)<eps) then
      it=0
    else
      if (lo.ne.k) then
        do j=k, n
          t0=C(k,j); C(k,j)=C(lo,j); C(lo,j)=t0
        end do
      end if
      if (ko.ne.k) then
        do i=1, n
          t0=C(i,k); C(i,k)=C(i,ko); C(i,ko)=t0
        end do
      end if
      do i=k+1, n
        C(i,k)=C(i,k)/po
        do j=k+1, n
          C(i,j)=C(i,j)-C(i,k)*C(k,j)
        end do
      end do
      k=k+1
    end if
  end do
  if (it==1.and.dabs(C(n,n))<eps)  it=0

End !TSRGT

!*****************************************************************
 !* TSCGT procedure implements the triangularization algorithm of *
 !* Gauss with full pivoting at each step for a complex matrix, A *
 !* and saves the made transformations in KP and LP.              *
 !* ------------------------------------------------------------- *
 !* INPUTS:                                                       *
 !*          N:   size of complex matrix A                        *
 !*          A:   complex matrix of size N x N                    *
 !* OUTPUTS;                                                      *
 !*          it:  =0 if A is singular, else =1.                   *
 !*           C:  contains the upper triangular matrix and the    *
 !*               multipliers used during the process.            *
 !*          KP:  contains the column exchanges.                  *
 !*          LP:  contains the line exchanges.                    *
 !*****************************************************************
  Subroutine TSCGT(eps, N, A, it, C, KP, LP)
    implicit none
    real(kind=8), intent(in) :: eps
    integer, intent(in) :: N
    Complex(kind=16), intent(in) ::  A(:,:)
    Complex(kind=16), intent(out) :: C(:,:)
    Integer, intent(out) :: KP(:), LP(:), it
    ! internal
    Complex(kind=16) :: C0,C1,P0,T0
    integer :: k,i,j,k0,l0


    C=A
    it=1; K=1
    do while (it==1.and.k<N)
      P0=C(k,k); l0=k; k0=k
      do i=k, N
        do j=1, N
          if (ABS(C(i,j)) > ABS(P0)) then
            P0=C(i,j)
            l0=i; k0=j
          end if
        end do
      end do
      LP(k)=l0; KP(k)=k0
      if (ABS(P0) < eps) then
        it=0
      else
        if (l0.ne.k) then
          do j=k, N
            T0=C(k,j)
            C(k,j)=C(l0,j)
            C(l0,j)=T0
          end do
        end if
        if (k0.ne.k) then
          do i=1, N
            T0=C(i,k)
            C(i,k)=C(i,k0)
            C(i,k0)=T0
          end do
        end if
        do i=k+1, N
          C0=C(i,k)
          C(i,k) = C0/P0
          do j=k+1, N
            C0=C(i,j)
            C1 = C(i,k)*C(k,j)
            C(i,j) = C0 - C1
          end do
        end do
        k=k+1
      end if
    end do
    if (it==1.and.ABS(C(N,N)) < eps)  it=0

  End  !TSCGT

end module m_matrix
