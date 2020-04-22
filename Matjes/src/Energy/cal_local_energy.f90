module m_local_energy
use m_basic_types, only : vec_point
use m_derived_types, only : point_shell_Operator
use m_modes_variables, only : point_shell_mode

interface local_energy
  module procedure  local_energy_pointer
end interface
private
public :: local_energy

contains

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
integer :: i,N,j,dim_ham
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

end module m_local_energy
