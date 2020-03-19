module m_local_energy
use m_derived_types

interface local_energy
   module procedure local_energy_pointer
end interface local_energy

private
public :: local_energy

contains

subroutine local_energy_pointer(E_int,iomp,spin,E_line)
use m_energy_commons
use m_dipole_energy
use m_dipolar_field, only : i_dip
implicit none
! input
type(point_shell_mode), intent(in) :: spin
type(point_shell_Operator), intent(in) :: E_line
integer, intent(in) :: iomp
! ouput
real(kind=8), intent(out) :: E_int
! internal
integer :: i,N

N=size(spin%shell)
E_int=0.0d0

do i=1,N

   E_int=E_int+dot_product( spin%shell(1)%w , matmul(E_line%shell(i)%Op_loc,spin%shell(i)%w) )

enddo

if (i_dip) E_int=E_int+get_dipole_E(iomp)

end subroutine local_energy_pointer

end module m_local_energy
