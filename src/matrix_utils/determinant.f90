module m_determinant
implicit none

private
public :: determinant

interface determinant
  module procedure :: det_real,det_complex
end interface determinant

interface TSGT
  module procedure :: TSRGT,TSCGT
end interface TSGT

contains

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
    complex(kind=8) :: A(:,:), det
    complex(kind=8) :: C0,Z1
    complex(kind=8) :: C(N,N)
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
    complex(kind=8), intent(in) ::  A(:,:)
    complex(kind=8), intent(out) :: C(:,:)
    Integer, intent(out) :: KP(:), LP(:), it
    ! internal
    complex(kind=8) :: C0,C1,P0,T0
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

end module m_determinant
