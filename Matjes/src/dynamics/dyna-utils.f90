module m_dyna_utils
use m_derived_types

interface init_lattice
 module procedure init_spin_lattice_5d
end interface init_lattice

interface copy_lattice
 module procedure copy_Mat_Point
 module procedure copy_Mat_Mat
 module procedure copy_1D_1D
 module procedure copy_1Dp_1Dp
 module procedure copy_1Dp_1Dreal
end interface copy_lattice

interface test_size_mat
 module procedure Test_1D
end interface test_size_mat

private
public :: init_lattice,copy_lattice
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the spin lattice to 0.0 to prepare the matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_spin_lattice_5d(matrix)
implicit none
real(kind=8), intent(inout) :: matrix(:,:,:,:,:)
! internal variables
integer :: i,j,k,l,N(5)

N=shape(matrix)

do l=1,N(5)
  do k=1,N(4)
    do j=1,N(3)
      do i=1,N(2)
      matrix(:,i,j,k,l)=0.0d0
      enddo
    enddo
  enddo
enddo

end subroutine init_spin_lattice_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copy the content of a 1D vector into another one
! This routine must be fast because it is used all the time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine copy_1D_1D(in,out)
implicit none
type(vec),intent(out) :: out(:)
type(vec),intent(in) :: in(:)
! internal variables
integer :: i,N

N=size(out)

do i=1,N
  out(i)%w=in(i)%w
enddo

end subroutine copy_1D_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copy the content of a 1D vector pointer into another one
! This routine must be fast because it is used all the time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine copy_1Dp_1Dp(in,out)
implicit none
type(vec_point),intent(out) :: out(:)
type(vec_point),intent(in) :: in(:)
! internal variables
integer :: i,N

N=size(out)

do i=1,N
  out(i)%w=in(i)%w
enddo

end subroutine copy_1Dp_1Dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copy the content of a 1D vector into a 1D vector pointer
! This routine must be fast because it is used all the time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine copy_1Dp_1Dreal(in,out)
implicit none
type(vec_point),intent(out) :: out(:)
real(kind=8),intent(in) :: in(:,:)
! internal variables
integer :: i,N

N=size(out)
if (size(out(1)%w).ne.size(in(:,1))) stop 'size of matrices do not match - copy_1Dp_1D'

do i=1,N
  out(i)%w=in(:,i)
enddo

end subroutine copy_1Dp_1Dreal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copy the content of the pointer into the spin matrix
! or whatever matrix that comes in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine copy_Mat_Point(matrix,my_lattice)
implicit none
real(kind=8),intent(out) :: matrix(:,:,:,:,:)
type(lattice),intent(in) :: my_lattice
! internal variables
integer :: i,j,k,l,N_point(4),N_mat(5)

N_point=shape(my_lattice%l_modes)
N_mat=shape(matrix)

!!!!$omp do ordered private(i,j,k,l) schedule(auto) collapse(4)
do l=1,N_point(4)
  do k=1,N_point(3)
    do j=1,N_point(2)
      do i=1,N_point(1)
      matrix(:,i,j,k,l)=my_lattice%l_modes(i,j,k,l)%w
      enddo
    enddo
  enddo
enddo
!!!!$omp end do

end subroutine copy_Mat_Point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copy the content of the pointer into the spin matrix
! or whatever matrix that comes in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine copy_Mat_Mat(matrix1,matrix2)
implicit none
real(kind=8),intent(out) :: matrix1(:,:,:,:,:)
real(kind=8),intent(in) :: matrix2(:,:,:,:,:)
! internal variables
integer :: i,j,k,l,N_mat1(5),N_mat2(5)
logical :: test

N_mat1=shape(matrix1)
N_mat2=shape(matrix2)

test=test_size_mat(N_mat1,N_mat2)
if (.not.test) stop 'error'

!!!$omp do ordered private(i,j,k,l) schedule(auto) collapse(4)
do l=1,N_mat1(5)
  do k=1,N_mat1(4)
    do j=1,N_mat1(3)
      do i=1,N_mat1(2)
      matrix1(:,i,j,k,l)=matrix2(:,i,j,k,l)
      enddo
    enddo
  enddo
enddo
!!!!$omp end do

end subroutine copy_Mat_Mat

logical function Test_1D(N1,N2)
implicit none
integer, intent(in) :: N1(:),N2(:)
! internal
integer :: i,N(size(N1)),Nup,N_test

Nup=size(N1)
Test_1D=.True.

do i=1,Nup
  N(i)=N1(i)-N2(i)
enddo

N_test=sum(N)

if (N_test.ne.0) then
   write(6,'(a)') 'error in copy_Mat_Mat'
   write(6,'(a)') 'The matrix and the matrix of pointer to not have the same dimension'
   Test_1D=.False.
endif

end function Test_1D

end module m_dyna_utils
