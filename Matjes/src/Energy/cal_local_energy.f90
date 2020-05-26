module m_local_energy
use m_basic_types, only : vec_point
use m_derived_types, only : point_shell_Operator
use m_modes_variables, only : point_shell_mode

interface local_energy
  module procedure  local_energy_pointer,local_energy_optimized
end interface

! all vectors on one line (must be updated at each line)
real(kind=8), allocatable, dimension(:) :: all_vectors

! all E matrice on one line (must be done once since the table is always done in the same order)
real(kind=8), allocatable, dimension(:,:) :: all_E

private
public :: local_energy,get_E_matrix,kill_E_matrix

contains

!
! really too slow
!
subroutine local_energy_pointer(E_int,iomp,spin,dim_mode)
use m_energy_commons, only : energy
use m_dipole_energy
use m_dipolar_field, only : i_dip
use m_matrix, only : reduce
implicit none
! input
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: iomp,dim_mode
! ouput
real(kind=8), intent(out) :: E_int
! internal
integer :: i,N,j
real(kind=8) :: S_int(dim_mode)

N=size(energy%line(:,iomp))
E_int=0.0d0

do i=1,N

   j=energy%line(i,iomp)

   call reduce(energy%value(i,iomp),size(energy%value(i,iomp)%order_op),S_int,spin(iomp)%w,spin(j)%w,dim_mode)

   E_int=E_int+dot_product( spin(iomp)%w , S_int )

enddo

if (i_dip) E_int=E_int+get_dipole_E(iomp)

end subroutine local_energy_pointer


!
! much much faster
!

subroutine local_energy_optimized(E_int,iomp,spin)
use m_energy_commons, only : energy
use m_dipole_energy
use m_dipolar_field, only : i_dip
use m_matrix, only : reduce
implicit none
! input
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: iomp
! ouput
real(kind=8), intent(out) :: E_int
! internal
integer :: i,N,j,dim_mode

N=size(energy%line(:,iomp))
E_int=0.0d0
dim_mode=size(all_vectors)/N

do i=1,N
   j=energy%line(i,iomp)
   all_vectors((i-1)*dim_mode+1:i*dim_mode)=spin(j)%w
enddo

E_int=dot_product( spin(iomp)%w , matmul( all_E , all_vectors ) )

if (i_dip) E_int=E_int+get_dipole_E(iomp)

end subroutine local_energy_optimized





subroutine get_E_matrix(dim_mode)
use m_energy_commons, only : energy
implicit none
integer, intent(in) :: dim_mode
! internal
integer :: N,i,j

N=size(energy%line(:,1))

allocate(all_vectors(dim_mode*N),all_E(dim_mode,dim_mode*N))
all_vectors=0.0d0
all_E=0.0d0

! the B_total is always read in the same direction so we can fill it only once
! one has to do a transpose here

do i=1,N
   all_E(:,(i-1)*dim_mode+1:i*dim_mode)=transpose(energy%value(i,1)%order_op(1)%Op_loc)
enddo

end subroutine get_E_matrix

subroutine kill_E_matrix()
implicit none

deallocate(all_vectors,all_E)
write(6,'(a)') 'Energy matrix deallocated'

end subroutine kill_E_matrix

end module m_local_energy
