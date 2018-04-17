module m_dyna_utils
use m_derived_types

interface init_lattice
 module procedure init_spin_lattice_5d
end interface init_lattice

interface copy_lattice
 module procedure copy_Mat_Point
 module procedure copy_Mat_Mat
end interface copy_lattice

interface test_size_mat
 module procedure Test_1D
end interface test_size_mat

interface associate_pointer
 module procedure asso_pointer_to_matrix
end interface associate_pointer

private
public :: init_lattice,copy_lattice,associate_pointer
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

do l=1,N_point(4)
  do k=1,N_point(3)
    do j=1,N_point(2)
      do i=1,N_point(1)
      matrix(:,i,j,k,l)=my_lattice%l_modes(i,j,k,l)%w
      enddo
    enddo
  enddo
enddo

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

do l=1,N_mat1(5)
  do k=1,N_mat1(4)
    do j=1,N_mat1(3)
      do i=1,N_mat1(2)
      matrix1(:,i,j,k,l)=matrix2(:,i,j,k,l)
      enddo
    enddo
  enddo
enddo

end subroutine copy_Mat_Mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! associate the pointer to the matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine asso_pointer_to_matrix(my_lattice,matrix2)
implicit none
type(lattice),intent(inout) :: my_lattice
real(kind=8),target,intent(in) :: matrix2(:,:,:,:,:)
! internal variables
integer :: i,j,k,l,N_mat2(5),N_mat1(4)
logical :: test

! check that the dimensions are equal
N_mat1=shape(my_lattice%l_modes)
N_mat2=shape(matrix2)

test=test_size_mat(N_mat1,N_mat2(2:5))
if (.not.test) stop 'error'

do l=1,N_mat1(4)
  do k=1,N_mat1(3)
    do j=1,N_mat1(2)
      do i=1,N_mat1(1)
      my_lattice%l_modes(i,j,k,l)%w=>matrix2(:,i,j,k,l)
      enddo
    enddo
  enddo
enddo

end subroutine asso_pointer_to_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! associate the pointer to the matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
