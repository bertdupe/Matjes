module m_total_energy

!
! this module contains the routines that calculate the total energy for all the energy terms
!

private
public :: total_energy
contains

real(kind=8) function total_energy(N,spin)
use m_basic_types, only : vec_point
use m_local_energy, only : local_energy
implicit none
integer, intent(in) :: N
type(vec_point), intent(in) :: spin(:)
! internal
integer :: i

total_energy=0.0d0

do i=1,N
   call local_energy(total_energy,i,spin)
enddo

end function total_energy

end module
