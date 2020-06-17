module m_eval_Beff
use m_dipolar_field, only : i_dip,get_dipole_B
use m_internal_fields_commons, only : B_total
! the version norm of the calculation of Beff is 90% of the running time. Really too slow



interface calculate_Beff
    module procedure normal,order_2only!,order_N
end interface calculate_Beff

! all vectors on one line (must be updated at each line)
real(kind=8), allocatable, dimension(:) :: all_vectors

! all B matrice on one line (must be done once since the table is always done in the same order)
type all_order_all_B
  real(kind=8), allocatable, dimension(:,:) :: B_all_shell
  integer :: nline,n_column
end type all_order_all_B
type(all_order_all_B), allocatable, dimension(:) :: all_B

private
public :: calculate_Beff,get_B_matrix,kill_B_matrix
contains
! subroutine that calculates the field
! dE/DM
!
!
!
!--------------------------------------------------------------
! for normal (very slow)
!
subroutine normal(B,iomp,spin,dim_mode)
use m_basic_types, only : vec_point
use m_modes_variables, only : point_shell_mode
use m_matrix, only : reduce
implicit none
! input variable
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: iomp,dim_mode
! output of the function
real(kind=8), intent(out) :: B(:)
! internals
integer :: N,i,j
real(kind=8) :: S_int(dim_mode)

!N=B_total%ncolumn
N=size(B_total%line(:,iomp))
B=0.0d0
S_int=0.0d0

do i=1,N

   j=B_total%line(i,iomp)

   call reduce(B_total%value(i,iomp),size(B_total%value(i,iomp)%order_op),S_int,spin(iomp)%w,spin(j)%w,dim_mode)

   B=B+S_int

enddo

if (i_dip) call get_dipole_B(B,iomp)

end subroutine normal


!
!--------------------------------------------------------------
! store all vectors and all Matrices in a contiguous matrix before doing the calculation
!
subroutine order_2only(B,iomp,spin)
use m_basic_types, only : vec_point
implicit none
! input variable
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: iomp
! output of the function
real(kind=8), intent(out) :: B(:)
! internals
integer :: N,i,j,k
integer :: dim_mode

B=0.0d0
N=size(B_total%line(:,iomp))
dim_mode=size(B)

do i=1,N
  j=B_total%line(i,iomp)
  all_vectors((i-1)*dim_mode+1:i*dim_mode)=spin(j)%w
enddo

B=matmul(all_B(1)%B_all_shell,all_vectors)

if (i_dip) call get_dipole_B(B,iomp)

end subroutine order_2only

!--------------------------------------------------------------
! store all vectors and all Matrices in a contiguous matrix before doing the calculation
!
subroutine order_N(B,iomp,spin)
use m_basic_types, only : vec_point
implicit none
! input variable
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: iomp
! output of the function
real(kind=8), intent(out) :: B(:)
! internals
integer :: N,i,j,k,counter
integer :: dim_mode,size_B
real(kind=8),allocatable :: B_transfer(:)

B=0.0d0
N=size(B_total%line(:,iomp))
dim_mode=size(B)
size_B=size(all_B)

do k=size_B,1,-1
  allocate(B_transfer(dim_mode**2))
  counter=0
  do i=1,N
    if (size(B_total%value(i,1)%order_op).ge.k) then
      counter=counter+1
      j=iomp
      if (k.eq.1)j=B_total%line(i,iomp)
      all_vectors((counter-1)*dim_mode+1:counter*dim_mode)=spin(j)%w
    endif
  enddo

    B_transfer=matmul(all_B(k)%B_all_shell,all_vectors(1:counter*dim_mode))
    write(*,*) B_transfer
   pause
!  B=matmul(all_B(k)%B_all_shell,all_vectors(1:counter*dim_mode))
enddo

if (i_dip) call get_dipole_B(B,iomp)

end subroutine order_N

subroutine get_B_matrix(dim_mode)
implicit none
integer, intent(in) :: dim_mode
! internal
integer :: N,i,j,counter
integer :: max_order,N_site
integer, allocatable :: order(:)
logical :: high_order_ham

high_order_ham=.false.
max_order=1
N=size(B_total%line(:,1))

!
! check if the variables were already allocated
!

if (allocated(all_B).and.allocated(all_vectors)) return

allocate(order(N))

do i=1,N
 order(i)=size(B_total%value(i,1)%order_op)
 if (order(i).gt.1) high_order_ham=.true.
enddo

if (high_order_ham) max_order=maxval(order)

! allocate the value correponding to the max order
allocate(all_B(max_order))
do i=1,max_order
  if (i.eq.1) then
    allocate(all_B(i)%B_all_shell(dim_mode,dim_mode*N))
    all_B(i)%B_all_shell=0.0d0
  else
    N_site=count(order.eq.i)
    allocate(all_B(i)%B_all_shell(dim_mode**i,dim_mode*N_site))
    all_B(i)%B_all_shell=0.0d0
  endif
enddo

! the B_total is always read in the same direction so we can fill it only once
! one has to do a transpose here

do j=1,max_order
  counter=0
  do i=1,N
    if (order(i).ge.j) then
      counter=counter+1
      all_B(j)%B_all_shell(:,(counter-1)*dim_mode+1:counter*dim_mode)=transpose(B_total%value(i,1)%order_op(j)%Op_loc)
    endif
  enddo
enddo

! we always use the same vector but only part of it
allocate(all_vectors(dim_mode*N))
all_vectors=0.0d0

end subroutine get_B_matrix

subroutine kill_B_matrix()
implicit none

deallocate(all_vectors,all_B)
write(6,'(a)') 'Effective field matrix deallocated'

end subroutine kill_B_matrix

end module m_eval_Beff
