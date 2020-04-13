module m_total_energy
use m_derived_types, only : operator_real,lattice,point_shell_Operator
use m_Hamiltonian_variables, only : Coeff_Ham
use m_modes_variables, only : point_shell_mode

!
! this module contains the routines that calculate the total energy for all the energy terms
!

private
public :: total_energy
contains

real(kind=8) function total_energy(N,E_column,E_line)
use m_local_energy, only : local_energy
integer, intent(in) :: N
type(point_shell_Operator), intent(in) :: E_line(:)
type(point_shell_mode), intent(in) :: E_column(:)
! internal
integer :: i

total_energy=0.0d0

do i=1,N
   call local_energy(total_energy,i,E_column(i),E_line(i))
enddo

end function total_energy

end module
