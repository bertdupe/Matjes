module m_eval_Beff
use m_dipolar_field, only : i_dip,get_dipole_B
use m_internal_fields_commons, only : B_total
! the version norm of the calculation of Beff is 90% of the running time. Really too slow



interface calculate_Beff
    module procedure normal,optimized
end interface calculate_Beff

! all vectors on one line (must be updated at each line)
real(kind=8), allocatable, dimension(:) :: all_vectors

! all B matrice on one line (must be done once since the table is always done in the same order)
real(kind=8), allocatable, dimension(:,:) :: all_B



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
subroutine optimized(B,iomp,spin)
use m_basic_types, only : vec_point
implicit none
! input variable
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: iomp
! output of the function
real(kind=8), intent(out) :: B(:)
! internals
integer :: N,i,j,dim_mode

B=0.0d0
N=size(B_total%line(:,iomp))
dim_mode=size(B)

do i=1,N
   j=B_total%line(i,iomp)
   all_vectors((i-1)*dim_mode+1:i*dim_mode)=spin(j)%w
enddo

B=matmul(all_B,all_vectors)

if (i_dip) call get_dipole_B(B,iomp)

end subroutine optimized


subroutine get_B_matrix(dim_mode)
implicit none
integer, intent(in) :: dim_mode
! internal
integer :: N,i,j

N=size(B_total%line(:,1))

allocate(all_vectors(dim_mode*N),all_B(dim_mode,dim_mode*N))
all_vectors=0.0d0
all_B=0.0d0

! the B_total is always read in the same direction so we can fill it only once
! one has to do a transpose here

do i=1,N
   all_B(:,(i-1)*dim_mode+1:i*dim_mode)=transpose(B_total%value(i,1)%order_op(1)%Op_loc)
enddo

end subroutine get_B_matrix

subroutine kill_B_matrix()
implicit none

deallocate(all_vectors,all_B)
write(6,'(a)') 'Effective field matrix deallocated'

end subroutine kill_B_matrix

end module m_eval_Beff
