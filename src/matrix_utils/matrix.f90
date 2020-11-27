module m_matrix

!
! contains utilities for basic matrix operations
!
!
!

interface reduce
  module procedure :: reduce_real,reduce_Nreal
end interface reduce


private
public :: reduce

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

end module m_matrix
