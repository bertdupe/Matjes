module m_dipole_energy
use m_dipolar_field
!
! contains all the routines to calculate the dipole dipole energy interaction
!
!

private
public :: get_dipole_E

contains

!
! calculate the double sum of the dipole dipole energy
!
real(kind=8) function get_dipole_E(iomp)
use m_constants, only : mu_B
implicit none
! input
integer, intent(in) :: iomp
! internal variable

! internal variable
real(kind=8) :: B(3)

get_dipole_E=0.0d0
B=0.0d0

call get_dipole_B(B,iomp)

!
! the dipolar is 2 times larger to take into account the sommation on the atoms
! so the energy should be divided by 2
!

get_dipole_E=-dot_product(mode_dipole(iomp)%w,B)/2.0d0

end function get_dipole_E

end module m_dipole_energy
