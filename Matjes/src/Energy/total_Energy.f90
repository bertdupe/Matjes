module m_total_energy

!
! this module contains the routines that calculate the total energy for all the energy terms
!

private
public :: total_energy
contains

real(kind=8) function total_energy(N,lat)
use m_type_lattice, only : lattice
use m_local_energy, only : sum_energy
implicit none
integer, intent(in) :: N
type(lattice), intent(in) :: lat

Call sum_energy(total_energy,lat)

end function total_energy

end module
