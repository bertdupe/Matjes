module m_matrix

!
! contains utilities for basic matrix operations
!
!
!

interface invert
  module procedure :: invert_real
end interface invert

interface reduce
  module procedure :: reduce_real,reduce_op,reduce_Nreal
end interface reduce

interface Jacobi
  module procedure :: eigenvalvec_2D
end interface

private
public :: invert,reduce,Jacobi

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
    if (sum( abs( matrix_int(:,i) ) ).lt.1.0d-8) cycle
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
    if (sum( abs( matrix_int(:,i) ) ).lt.1.0d-8) cycle
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



end module m_matrix
